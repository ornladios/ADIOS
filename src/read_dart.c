/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */


/**************************************************/
/* Read method for DART memory-to-memory coupling */
/**************************************************/

#include <stdlib.h>
#include <string.h>
#include "adios.h"
#include "bp_utils.h"
#include "bp_types.h"
#include "adios_types.h"
#include "adios_read.h"
#include "adios_read_hooks.h"
#include "adios_error.h"
#include "futils.h"

#include "dc_gspace.h"
#include "ss_data.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

// DART needs an application ID. Suppose 1 is the writer, so we elect to be 2.
static const int DART_APPLICATION_ID = 2; 
static struct dcg_space *dcg = 0;
static enum storage_type st = column_major;

struct adios_read_dart_data_struct {
    int access_version;        // counting the access
    int disconnect_at_fclose;  // disconnect from DART in fclose()
    int mpi_rank;              // just for debug prints
};

// Declarations
static int adios_read_dart_get (const char * varname, enum ADIOS_DATATYPES vartype, 
                                struct adios_read_dart_data_struct * ds, 
                                int * offset, int * readsize, void * data);

/* If init is used, we connect to DART here, otherwise we connect in fopen.
   If multiple fopen..fclose cycles are used, init/finalize must be used too to
   avoid multiple connection/disconnection in fopen/fclose.
*/
int adios_read_dart_init (MPI_Comm comm) 
{ 
    int  nproc, dart_rank, dart_peers;
    int  rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    /* Fortran call dart_init (nproc, 2): */
    printf("-- %s, rank %d: connect to dart with nproc=%d and appid=%d\n", __func__, rank, nproc, DART_APPLICATION_ID);
    dcg = dcg_alloc(nproc, DART_APPLICATION_ID);
    if (!dcg) {
        error(err_connection_failed, "Failed to connect with DART\n");
        return -err_connection_failed;
    }

    /* Fortran call dart_rank(dartrank): */
    dart_rank = dcg_get_rank(dcg);
    /* Fortran call dart_peers(dartpeers) */
    dart_peers = dcg_get_num_peers(dcg);

    return 0; 
}

int adios_read_dart_finalize () 
{ 
    if (dcg) {
        /* Fortran call dart_finalize */
        dcg_free(dcg);
        dcg = NULL;
    }
    printf("-- %s: disconnect from dart\n", __func__);
    return 0; 
}

ADIOS_FILE * adios_read_dart_fopen (const char * fname, MPI_Comm comm)
{
    ADIOS_FILE * fp;
    struct adios_read_dart_data_struct * ds;
    int i;    

    adios_errno = 0;

    ds = (struct adios_read_dart_data_struct *) malloc (sizeof(struct adios_read_dart_data_struct));
    if (!ds) {
        error( err_no_memory, "Cannot allocate memory for file info.");
        return NULL;
    }

    /* fill out dart method specific struct */
    ds->access_version = 0;
    MPI_Comm_rank(comm, &ds->mpi_rank);

    /* if not connected to DART, connect now (and disconnect in fclose) */
    if (!dcg) {
        printf("-- %s, rank %d: call init first\n", __func__, ds->mpi_rank);
        if (!adios_read_dart_init(comm)) {
            return NULL;
        }
        ds->disconnect_at_fclose = 1;
    } else {
        ds->disconnect_at_fclose = 0;
    }

    /* Try to get variable with fname. If it does not exists, we get an error, which means
       the data does not exist. So we return an error just like with real files */
    int offset[] = {0,0,0}, readsize[3] = {1,1,1};
    int time_index;
    int err;
    enum ADIOS_DATATYPES time_index_type = adios_integer;
    printf("-- %s, rank %d: Get variable %s\n", __func__, ds->mpi_rank, fname);
    printf("   rank %d: call dcg_lock_on_read()\n", ds->mpi_rank);
    dcg_lock_on_read();
    printf("   rank %d: dart_get %s\n", ds->mpi_rank, fname);
    err = adios_read_dart_get(fname, time_index_type, ds, offset, readsize, &time_index);
    printf("   rank %d: call dcg_unlock_on_read()\n", ds->mpi_rank);
    dcg_unlock_on_read();
    if (err) {
        error(err_file_not_found_error, "Data of '%s' does not exist in DataSpaces\n", fname);
        return NULL;
    } else {
        printf("-- %s, rank %d: data of '%s' exists, time index = %d\n", __func__, ds->mpi_rank, fname, time_index);
    }


    fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    if (!fp) {
        error( err_no_memory, "Cannot allocate memory for file info.");
        return NULL;
    }

    /* fill out ADIOS_FILE struct */
    fp->fh = (uint64_t) ds;
    fp->groups_count = 1;
    fp->vars_count = 0;
    fp->attrs_count = 0;
    fp->tidx_start = 0;
    fp->ntimesteps = 1;
    fp->file_size = 0;
    fp->version = 1;
    fp->endianness = 0; /* FIXME: not always Little Endian. Does it matter? */
    alloc_namelist (&fp->group_namelist,fp->groups_count); 
    for (i=0;i<fp->groups_count;i++) {
        if (!fp->group_namelist[i]) {
            error(err_no_memory, "Could not allocate buffer for %d strings in adios_fopen()", fp->groups_count);
            adios_read_dart_fclose(fp);
            return NULL;
        }
        else  {
            strcpy(fp->group_namelist[i],"dart");
        }
    }
    printf("-- %s, rank %d: done fp=%x, fp->fh=%x\n", __func__, ds->mpi_rank, fp, fp->fh);
    return fp;
}

int adios_read_dart_fclose (ADIOS_FILE *fp) 
{
    struct adios_read_dart_data_struct * ds = (struct adios_read_dart_data_struct *) fp->fh;
    int i,j;

    adios_errno = 0;

    printf("-- %s, rank %d: fp=%x\n", __func__, ds->mpi_rank, fp);
    /* Disconnect from DART if we connected in fopen() */
    if (ds && ds->disconnect_at_fclose) {
        adios_read_dart_finalize();
    }

    free_namelist ((fp->group_namelist),fp->groups_count);
    free(ds);
    free(fp);
    return 0;
}


ADIOS_GROUP * adios_read_dart_gopen (ADIOS_FILE *fp, const char * grpname)
{
    /* DART has no groups, so any grpname is accepted and the same empty stuff is returned */
    return adios_read_dart_gopen_byid(fp, 0);
}

ADIOS_GROUP * adios_read_dart_gopen_byid (ADIOS_FILE *fp, int grpid)
{
    struct adios_read_dart_data_struct * ds = (struct adios_read_dart_data_struct *) fp->fh;
    ADIOS_GROUP * gp;

    /* DART has no groups, so any grpid is accepted and the same empty stuff is returned */
    /* But we perform dart_lock_on_read here and release it in gclose */

    adios_errno = 0;
    gp = (ADIOS_GROUP *) malloc(sizeof(ADIOS_GROUP));
    if (!gp) {
        error( err_no_memory, "Could not allocate memory for group info");
        return NULL;
    }

    /* fill out ADIOS_GROUP struct */
    gp->grpid = grpid;
    gp->gh = (uint64_t) 0;
    gp->fp = fp;
    gp->vars_count = 0;
    gp->attrs_count = 0;
    gp->var_namelist = 0;
    gp->attr_namelist = 0;
    
    printf("-- %s, rank %d: call dcg_lock_on_read() gp=%x fp=%x\n", __func__, ds->mpi_rank, gp, gp->fp);
    dcg_lock_on_read();
    printf("-- %s, rank %d: returned from dcg_lock_on_read()\n", __func__);

    return gp;
}
                   
int adios_read_dart_gclose (ADIOS_GROUP *gp)
{
    struct adios_read_dart_data_struct * ds = (struct adios_read_dart_data_struct *) gp->fp->fh;

    adios_errno = 0;
    /* Release read lock locked in gopen */
    dcg_unlock_on_read();

    free_namelist ((gp->var_namelist),gp->vars_count);
    free_namelist ((gp->attr_namelist),gp->attrs_count);
    free(gp);
    return 0;
}



int adios_read_dart_get_attr (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type,
                    int * size, void ** data)
{
    /* DART does not support attributes */
    error(err_invalid_attrname, "DART read method does not support attributes!");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}

int adios_read_dart_get_attr_byid (ADIOS_GROUP * gp, int attrid, 
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    /* DART does not support attributes */
    error(err_invalid_attrid, "DART read method does not support attributes!");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}


ADIOS_VARINFO * adios_read_dart_inq_var (ADIOS_GROUP *gp, const char * varname) 
{
    /* DART has no inquiry capability, report somthing dummy */
    return adios_read_dart_inq_var_byid(gp, 0);
}

ADIOS_VARINFO * adios_read_dart_inq_var_byid (ADIOS_GROUP *gp, int varid)
{
    struct adios_read_dart_data_struct * ds = (struct adios_read_dart_data_struct *) gp->fp->fh;
    ADIOS_VARINFO * vi;
    int i,k;

    adios_errno = 0;
    vi = (ADIOS_VARINFO *) malloc(sizeof(ADIOS_VARINFO));
    if (!vi) {
        error( err_no_memory, "Could not allocate memory for variable info.");
        return NULL;
    }

    /* DART has no inquiry capability, report somthing dummy */
    vi->varid = varid;
    vi->type = adios_unknown;
    vi->ndim = 0;
    vi->dims = NULL;
    vi->timedim = -1;
    vi->value = NULL;
    vi->gmin = NULL;
    vi->gmax = NULL;
    vi->mins = NULL;
    vi->maxs = NULL;
    
    return vi;
}

void adios_read_dart_free_varinfo (ADIOS_VARINFO *vp)
{
    if (vp) {
        if (vp->dims)   free(vp->dims);
        if (vp->value)  free(vp->value);
        if (vp->gmin && vp->gmin != vp->value)   free(vp->gmin);
        if (vp->gmax && vp->gmax != vp->value)   free(vp->gmax);
        if (vp->mins)   free(vp->mins);
        if (vp->maxs)   free(vp->maxs);
        free(vp);
    }
}

static int adios_read_dart_get (const char * varname, enum ADIOS_DATATYPES vartype, 
                                struct adios_read_dart_data_struct * ds, 
                                int * offset, int * readsize, void * data)
{
    struct obj_data *od;
    int elemsize = common_read_type_size(vartype, NULL);
    int err;
    /* Fortran call dart_get( name, version, {bounding box}, data) */
    struct obj_descriptor odsc = {
        .version = ds->access_version, 
        .owner = -1,
        .size = elemsize,
        .st = st, // column_major,
        .bb = {
            .num_dims = 3,
            .lb.c = {offset[1], offset[0], offset[2]}, // !! flip 1st and 2nd dim
            .ub.c = {offset[1]+readsize[1]-1, 
                     offset[0]+readsize[0]-1,
                     offset[2]+readsize[2]-1,}
        }
    };
    strncpy(odsc.name, varname, sizeof(odsc.name));

    od = obj_data_alloc_no_data(&odsc, data);
    if (!od) {
        error(err_no_memory, "DART obj_data_alloc_no_data() failed");  
        return -err_no_memory;
    }

    printf("-- %s, rank %d: get data: varname=%s odsc={.version=%d, lb=(%d,%d,%d) ub=(%d,%d,%d)}\n",
        __func__, ds->mpi_rank, varname, odsc.version, odsc.bb.lb.c[0], odsc.bb.lb.c[1], odsc.bb.lb.c[2],
        odsc.bb.ub.c[0], odsc.bb.ub.c[1], odsc.bb.ub.c[2]);
    err = dcg_obj_get(od);
    free(od);

    if (err) {
        error(err_corrupted_variable, "DART dcg_obj_get() failed to read variable %s.", varname);  
        return -err_corrupted_variable;
    }

    return 0;
}

int64_t adios_read_dart_read_var (ADIOS_GROUP * gp, const char * varname,
                        const uint64_t * start, const uint64_t * count,
                        void * data)
{
    int64_t total_size;
    int offset[3], readsize[3], Toffset[3], Treadsize[3];
    struct adios_read_dart_data_struct * ds = (struct adios_read_dart_data_struct *) gp->fp->fh;
    char dart_type_var_name[128]; 
    enum ADIOS_DATATYPES vartype;
    int elemsize;
    int err;
    int i;

    /* DART uses integers for boundaries */
    total_size = 1;
    for (i=0; i<3; i++) {
        offset[i]    = (int) start[i];
        Toffset[i]   = 0;
        readsize[i]  = (int) count[i];
        Treadsize[i] = 1;
        total_size   = total_size * count[i];
    }

    /* Get type information for the variable from DataSpaces:
       type variable name = T<varname>
    */
    snprintf(dart_type_var_name, 128, "T%s", varname);
    err = adios_read_dart_get(dart_type_var_name, adios_integer, ds, Toffset, Treadsize, &vartype);
    if (err)
        return err;

    elemsize = common_read_type_size(vartype, NULL);
    printf("-- %s, rank %d: get data: varname=%s type=%d (%s) elemsize=%d}\n",
        __func__, ds->mpi_rank, varname, vartype, common_read_type_to_string(vartype), elemsize);

    total_size *= elemsize; 

    printf("-- %s, rank %d: get data: varname=%s start=(%lld,%lld,%lld) count=(%lld,%lld,%lld)}\n",
        __func__, ds->mpi_rank, varname, start[0], start[1], start[2], count[0], count[1], count[2]);
    printf("-- %s, rank %d: get data: varname=%s offset=(%d,%d,%d) readsize=(%d,%d,%d)}\n",
        __func__, ds->mpi_rank, varname, offset[0], offset[1], offset[2], readsize[0], readsize[1], readsize[2]);

    err = adios_read_dart_get(varname, vartype, ds, offset, readsize, data);
    if (err)
        return err;

    return total_size;
}

int64_t adios_read_dart_read_var_byid (ADIOS_GROUP    * gp,
                             int              varid,
                             const uint64_t  * start,
                             const uint64_t  * count,
                             void           * data)
{
    error( err_invalid_varid, "DART does not know variable indicies, only variable names can be used.");
    return -err_invalid_varid;
}


