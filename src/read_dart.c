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
//#include <errno.h>  /* ENOMEM */
#include "adios.h"
#include "bp_utils.h"
#include "bp_types.h"
#include "adios_types.h"
#include "adios_read.h"
#include "adios_read_hooks.h"
#include "adios_error.h"
#include "futils.h"
#include "globals.h"

#include "dart.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#if 1
#   define DBG_PRINTF printf
#else
#   define DBG_PRINTF(a,...) 
#endif

#define MAXDARTNAMELEN 128

/*#define DART_DO_VERSIONING   define it at configure as -DDART_DO_VERSIONING in CFLAGS */
static int number_of_fopens = 0;  /* for versioning, works only if one file is fopened (in a loop) in the application */

struct adios_read_dart_data_struct {
    char *fname;               // path of file 
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
    int  nproc, drank, dpeers;
    int  rank, err;
    int  appid, was_set;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    /* Connect to DART, but only if we are not yet connected (from Write API) */
    if (!globals_adios_is_dart_connected()) {
        appid = globals_adios_get_application_id (&was_set);
        if (!was_set) 
            appid = 2;
        DBG_PRINTF("-- %s, rank %d: connect to dart with nproc=%d and appid=%d\n", __func__, rank, nproc, appid);
        err = dart_init(nproc, appid);
        if (err < 0) {
            error(err_connection_failed, "Failed to connect with DART\n");
            return -err_connection_failed;
        }

        //drank = dart_rank(dcg);
        //dpeers = dart_peers(dcg);
    }
    globals_adios_set_dart_connected_from_reader();
    return 0; 
}

int adios_read_dart_finalize () 
{ 
    // disconnect from DART only if we the reader is connected (the writer not anymore)
    if (globals_adios_is_dart_connected_from_reader() && 
        !globals_adios_is_dart_connected_from_both()) 
    {
        dart_finalize();
        DBG_PRINTF("-- %s: disconnected from dart\n", __func__);
    }
    globals_adios_set_dart_disconnected_from_reader();
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
    ds->fname = strdup(fname);
#ifdef DART_DO_VERSIONING
    ds->access_version = number_of_fopens;  /* Read data of separate versions from DataSpaces */
#else
    ds->access_version = 0;    /* Data in DataSpaces is always overwritten (read same version) */
#endif
    MPI_Comm_rank(comm, &ds->mpi_rank);

    /* if not connected to DART, connect now (and disconnect in fclose) */
    if (!globals_adios_is_dart_connected_from_reader()) {
        DBG_PRINTF("-- %s, rank %d: call init first\n", __func__, ds->mpi_rank);
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
    char dart_fname[MAXDARTNAMELEN];
    snprintf(dart_fname, MAXDARTNAMELEN, "FILE@%s",fname);
    DBG_PRINTF("-- %s, rank %d: Get variable %s\n", __func__, ds->mpi_rank, dart_fname);
    /* We perform dart_lock_on_read here and release it in fclose */
    DBG_PRINTF("   rank %d: call dcg_lock_on_read(%s)\n", ds->mpi_rank, fname);
    dart_lock_on_read(fname);
    DBG_PRINTF("   rank %d: dart_get %s\n", ds->mpi_rank, dart_fname);
    err = adios_read_dart_get(dart_fname, time_index_type, ds, offset, readsize, &time_index);
    if (err) {
        error(err_file_not_found_error, "Data of '%s' does not exist in DataSpaces\n", dart_fname);
        DBG_PRINTF("   rank %d: call dcg_unlock_on_read(%s)\n", ds->mpi_rank, fname);
        dart_unlock_on_read(fname);
        free(ds);
        return NULL;
    } else {
        DBG_PRINTF("-- %s, rank %d: data of '%s' exists, time index = %d\n", __func__, ds->mpi_rank, dart_fname, time_index);
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
    DBG_PRINTF("-- %s, rank %d: done fp=%x, fp->fh=%x\n", __func__, ds->mpi_rank, fp, fp->fh);
#ifdef DART_DO_VERSIONING
    number_of_fopens++;
#endif
    return fp;
}

int adios_read_dart_fclose (ADIOS_FILE *fp) 
{
    struct adios_read_dart_data_struct * ds = (struct adios_read_dart_data_struct *) fp->fh;
    int i,j;

    adios_errno = 0;

    DBG_PRINTF("-- %s, rank %d: fp=%x\n", __func__, ds->mpi_rank, fp);

    /* Release read lock locked in fopen */
    dart_unlock_on_read(ds->fname);
    
    /* Disconnect from DART if we connected in fopen() */
    if (ds && ds->disconnect_at_fclose) {
        adios_read_dart_finalize();
    }

    free_namelist ((fp->group_namelist),fp->groups_count);
    if (ds->fname) { free(ds->fname); ds->fname = 0; }
    free(ds);
    free(fp);
    return 0;
}

/* This function can be called if user places 
   the wrong sequences of dims for a var 
*/
void adios_read_dart_reset_dimension_order (ADIOS_FILE *fp, int is_fortran)
{
    /* unimplemented */
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
    
    return gp;
}
                   
int adios_read_dart_gclose (ADIOS_GROUP *gp)
{
    struct adios_read_dart_data_struct * ds = (struct adios_read_dart_data_struct *) gp->fp->fh;

    adios_errno = 0;

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

    DBG_PRINTF("-- %s, rank %d: get data: varname=%s version=%d, lb=(%d,%d,%d) ub=(%d,%d,%d)}\n",
        __func__, ds->mpi_rank, varname, ds->access_version, offset[1], offset[0], offset[2],
        offset[1]+readsize[1]-1, offset[0]+readsize[0]-1, offset[2]+readsize[2]-1);

    err =  dart_get (varname, ds->access_version, elemsize, 
                     offset[1], offset[0], offset[2],
                     offset[1]+readsize[1]-1,
                     offset[0]+readsize[0]-1,
                     offset[2]+readsize[2]-1,
                     data
                    );
    /*if (err == -ENOMEM) {
        error(err_no_memory, "Not enough memory for DART to perform dart_get()");  
        return -err_no_memory;
    } 
    else*/ if (err) {
        error(err_corrupted_variable, "DART failed to read variable %s.", varname);  
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
    enum ADIOS_DATATYPES vartype;
    int elemsize;
    int err;
    int i;
    char dart_name[MAXDARTNAMELEN];

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
       type variable name = TYPE@<filename>/<varname>
    */
    snprintf(dart_name, MAXDARTNAMELEN, "TYPE@%s/%s", ds->fname, varname);
    err = adios_read_dart_get(dart_name, adios_integer, ds, Toffset, Treadsize, &vartype);
    if (err)
        return err;
    DBG_PRINTF("-- %s, rank %d: get type: varname=%s type=%d (%s)}\n",
        __func__, ds->mpi_rank, dart_name, vartype, common_read_type_to_string(vartype)) ;

    elemsize = common_read_type_size(vartype, NULL);
    DBG_PRINTF("-- %s, rank %d: get data: varname=%s type=%d (%s) elemsize=%d}\n",
        __func__, ds->mpi_rank, dart_name, vartype, common_read_type_to_string(vartype), elemsize);

    total_size *= elemsize; 

    DBG_PRINTF("-- %s, rank %d: get data: varname=%s start=(%lld,%lld,%lld) count=(%lld,%lld,%lld)}\n",
        __func__, ds->mpi_rank, dart_name, start[0], start[1], start[2], count[0], count[1], count[2]);
    DBG_PRINTF("-- %s, rank %d: get data: varname=%s offset=(%d,%d,%d) readsize=(%d,%d,%d)}\n",
        __func__, ds->mpi_rank, dart_name, offset[0], offset[1], offset[2], readsize[0], readsize[1], readsize[2]);

    snprintf(dart_name, MAXDARTNAMELEN, "%s/%s", ds->fname, varname);
    err = adios_read_dart_get(dart_name, vartype, ds, offset, readsize, data);
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


