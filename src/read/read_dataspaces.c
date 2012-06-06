/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */


/**************************************************/
/* Read method for DATASPACES memory-to-memory coupling */
/**************************************************/

#include <stdlib.h>
#include <string.h>
//#include <errno.h>  /* ENOMEM */
#include "public/adios_types.h"
#include "public/adios_read.h"
#include "public/adios_error.h"
#include "core/globals.h"
#include "core/util.h"
#include "core/adios_logger.h"
//#include "core/bp_types.h"
#include "core/adios_read_hooks.h"
#include "core/futils.h"
#include "core/ds_metadata.h"

#include "dart_interface.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#include <sys/time.h>
struct timeval adios_timer_tp;
static inline double time_get() 
{ 
    gettimeofday(&adios_timer_tp, NULL); \
    return  ((double)adios_timer_tp.tv_sec + ((double)adios_timer_tp.tv_usec)/1000000.0); 
}

#define time_current(timer) { gettimeofday(&adios_timer_tp, NULL); \
    timer =  ((double)adios_timer_tp.tv_sec + ((double)adios_timer_tp.tv_usec)/1000000.0); }

#define timer_start(timer) { gettimeofday(&adios_timer_tp, NULL); \
    timer -= ((double)adios_timer_tp.tv_sec + ((double)adios_timer_tp.tv_usec)/1000000.0); }

#define timer_end(timer) { gettimeofday(&adios_timer_tp, NULL); \
    timer += ((double)adios_timer_tp.tv_sec + ((double)adios_timer_tp.tv_usec)/1000000.0); }


#define MAX_DS_NAMELEN 128
/* Maximum number of different filenames allowed per process during the whole run */
#define MAXNFILE 10 
/*#define DATASPACES_NO_VERSIONING   define it at configure as -DDATASPACES_NO_VERSIONING in CFLAGS */

struct dataspaces_fileversions_struct { // describes one variable (of one group)
    char      * filename[MAXNFILE];
    int         version[MAXNFILE];  /* for versioning of one given filename */
};
static struct dataspaces_fileversions_struct file_versions;
static int n_filenames; /* number of filenames appeared during the run */


struct dataspaces_var_struct { // describes one variable (of one group)
    char                 * name;
    enum ADIOS_DATATYPES   type;
    int                    hastime; // 0: no, 1:yes (time dimension is not stored in dataspaces)
    int                    ndims;
    uint64_t               dims[4]; // we have max 3+time dims in DataSpaces
    void                 * value;
};

struct dataspaces_attr_struct { // describes one attribute (of one group)
    char                 * name;
    enum ADIOS_DATATYPES   type;
    void                 * value;
};

struct dataspaces_group_struct { // accessible as fp->fh->groups[grpid]
    int group_index_len;         // length of group index in GROUP@fn/gn variable 
    int nvars;                   // number of vars in this group
    int nattrs;                  // number of attrs in this group
};

struct dataspaces_data_struct { // accessible as fp->fh
    char *fname;                // path of file 
    int access_version;         // counting the access
    int disconnect_at_close;    // disconnect from DATASPACES in fclose()
    int mpi_rank;               // just for debug prints
    struct dataspaces_var_struct  * vars;  // number of vars is ADIOS_FILE->nvars
    struct dataspaces_attr_struct * attrs; // number of attrs is ADIOS_FILE->nattrs
    /* Group info */
    int ngroups;                // number of groups (=1 right now)
    char **group_namelist;      // names of the groups
    struct dataspaces_group_struct * groups;
};

// Declarations
static int adios_read_dataspaces_get (const char * varname, enum ADIOS_DATATYPES vartype, 
                                struct dataspaces_data_struct * ds, 
                                int ndims, int is_fortran_ordering, 
                                int * offset, int * readsize, void * data);

/* If init is used, we connect to DATASPACES here, otherwise we connect in fopen.
   If multiple fopen..fclose cycles are used, init/finalize must be used too to
   avoid multiple connection/disconnection in fopen/fclose.
*/
int adios_read_dataspaces_init_method (MPI_Comm comm, const char * parameters) 
{ 
    int  nproc, drank, dpeers;
    int  rank, err;
    int  appid, was_set;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    PairStruct *params = text_to_name_value_pairs (parameters);
    PairStruct *p = params;
    while (p) {
        if (!strcasecmp (p->name, "appid")) {
            appid = strtol(p->value, NULL, 10);
            if (appid > 0) {
                log_debug ("Appid parameter set to %d for DATASPACES read method\n", appid);
                globals_adios_set_application_id (appid);
            } else {
                log_error ("Invalid appid parameter given to the DATASPACES read method: '%s'\n", p->value);
            }
        } else if (!strcasecmp (p->name, "max_memory")) {
            log_debug ("max_memory parameter set for DATASPACES read method\n");

        } else if (!strcasecmp (p->name, "verbose")) {
            log_debug ("Parameter name %s is not recognized by the DATASPACES read method\n", p->name);

        } else {
            log_error ("Parameter name %s is not recognized by the DATASPACES read method\n", p->name);
        }
    }
    free_name_value_pairs (params);


    /* Connect to DATASPACES, but only if we are not yet connected (from Write API) */
    if (!globals_adios_is_dart_connected()) {
        appid = globals_adios_get_application_id (&was_set);
        if (!was_set) 
            appid = 2;
        log_debug("-- %s, rank %d: connect to dataspaces with nproc=%d and appid=%d\n", __func__, rank, nproc, appid);
        err = dart_init(nproc, appid);
        if (err < 0) {
            adios_error (err_connection_failed, "Failed to connect with DATASPACES\n");
            return err_connection_failed;
        }

        //drank = dart_rank(dcg);
        //dpeers = dart_peers(dcg);
    }
    globals_adios_set_dart_connected_from_reader();
    n_filenames = 0;
    return 0; 
}

int adios_read_dataspaces_finalize_method () 
{ 
    // disconnect from DATASPACES only if we the reader is connected (the writer not anymore)
    if (globals_adios_is_dart_connected_from_reader() && 
        !globals_adios_is_dart_connected_from_both()) 
    {
        dart_finalize();
        log_debug("-- %s: disconnected from dataspaces\n", __func__);
    }
    globals_adios_set_dart_disconnected_from_reader();
}

ADIOS_FILE * ds_unpack_file_info (char * buf, int buf_len, 
                                  /* OUT */ struct dataspaces_data_struct * ds)
{
    ADIOS_FILE *fp;
    struct dataspaces_group_struct * groups; // will be added to ds->groups
    char * b = buf;
    int blen, glen, i;

    if (!buf || buf_len < 21)
        return NULL;

    fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    if (!fp) {
        adios_error (err_no_memory, "Cannot allocate memory for file info.");
        return NULL;
    }

    // add the single group to ds
    groups = (struct dataspaces_group_struct *) malloc (sizeof(struct dataspaces_data_struct));
    if (!groups) {
        adios_error (err_no_memory, "Cannot allocate memory for group struct in file info.");
        return NULL;
    }
    ds->ngroups = 1;  // FIXME: more groups per file in the future?
    ds->groups = groups;

    blen = *(int*)b;  // buf len again from buffer itself
    if (blen != buf_len) {
        log_debug("WARNING: %s(): expected file info  buffer length is %d but buffer head says %d\n", 
                    __func__, buf_len, blen);
    }
    b += sizeof(int); // skip buf len
    fp->current_step = *(int*)b; // time index
    fp->last_step = fp->current_step;
    b += sizeof(int);
    fp->nvars = *(int*)b; // number of variables
    b += sizeof(int);
    fp->nattrs = *(int*)b; // number of attributes
    b += sizeof(int);
    groups[0].group_index_len = *(int*)b; // length of group index
    b += sizeof(int);
    glen = *(int*)b; // length of (only) group name
    b += sizeof(int);

    fp->fh = (uint64_t)ds; // note: repeated in the caller fopen(), where is the right place?
    fp->nsubfiles = 1; /* FIXME: get the number of DataSpaces servers here */
    fp->nwriters = 1; /* FIXME: get the number of writers into the packed into */
    fp->file_size = 0;
    fp->version = 1;
    fp->endianness = 0; // FIXME: not always Little Endian. Does it matter? 
    alloc_namelist (&ds->group_namelist, ds->ngroups);
    for (i=0;i<ds->ngroups;i++) {
        if (!ds->group_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate buffer for %d strings in adios_fopen()", 
                        ds->ngroups);
            if (i>0) free_namelist ((ds->group_namelist),i);
            free(fp);
            return NULL;
        }
        else  {
            strcpy(ds->group_namelist[i],b);
        }
    }

    return fp;
}


ADIOS_FILE * adios_read_dataspaces_open_file (const char * fname, MPI_Comm comm)
{
    adios_error (err_operation_not_supported, 
                 "DATASPACES staging method does not support file mode for reading. "
                 "Use adios_read_open_stream() to open a dataset.\n");
    return NULL;
}


ADIOS_FILE * adios_read_dataspaces_open_stream (const char * fname, 
                                                MPI_Comm comm, 
                                                enum ADIOS_LOCKMODE lock_mode, 
                                                int timeout_msec)
{
    ADIOS_FILE * fp;
    struct dataspaces_data_struct * ds;
    int i;    

    adios_errno = 0;

    ds = (struct dataspaces_data_struct *) malloc (sizeof(struct dataspaces_data_struct));
    if (!ds) {
        adios_error (err_no_memory, "Cannot allocate memory for file info.");
        return NULL;
    }

    double t1 = time_get();
    double t2;

    /* fill out dataspaces method specific struct */
#ifdef DATASPACES_NO_VERSIONING
    ds->access_version = 0;    /* Data in DataSpaces is always overwritten (read same version) */
#else
   // check this filename's version number
    int fidx;
    for (fidx=0; fidx<n_filenames; fidx++) {
        if (!strcmp(fname, file_versions.filename[fidx]))
            break;
    }
    if (fidx == n_filenames) {
        if (n_filenames < MAXNFILE) {
            file_versions.filename[ n_filenames ] = strdup(fname);
            file_versions.version [ n_filenames ] = -1;
            n_filenames++;
        } else {
            log_error("rank %d: Too many different filenames has been used for adios_fopen().\n\tDATASPACES method allows max %d files\n", ds->mpi_rank, MAXNFILE);
            adios_error (err_too_many_files, "Too many different filenames has been used for adios_fopen().\n\tDATASPACES method allows max %d files", MAXNFILE);
            return NULL;
        }
    }
    ds->access_version = file_versions.version[fidx]+1;  /* Read data of separate versions from DataSpaces */
    log_debug("open version filename=%s version=%d fidx=%d\n", fname, ds->access_version, fidx);
#endif
    MPI_Comm_rank(comm, &ds->mpi_rank);

    /* if not connected to DATASPACES, connect now (and disconnect in adios_read_close) */
    if (!globals_adios_is_dart_connected_from_reader()) {
        log_debug("-- %s, rank %d: call init first\n", __func__, ds->mpi_rank);
        if (!adios_read_dataspaces_init(comm)) {
            return NULL;
        }
        ds->disconnect_at_close = 1;
    } else {
        ds->disconnect_at_close = 0;
    }

    /* Try to get variable with fname. If it does not exists, we get an error, which means
       the data does not exist. So we return an error just like with real files */
    int offset[] = {0,0,0}, readsize[3] = {128,1,1};
    char file_info_buf[128];
    int err;
    char ds_fname[MAX_DS_NAMELEN];
    snprintf(ds_fname, MAX_DS_NAMELEN, "FILE@%s",fname);
    log_debug("-- %s, rank %d: Get variable %s\n", __func__, ds->mpi_rank, ds_fname);
    /* We perform dart_lock_on_read here and release it in fclose (using fname as id) */
    log_debug("   rank %d: call dcg_lock_on_read(%s)\n", ds->mpi_rank, fname);
    dart_lock_on_read(fname);

    /* While loop for handling timeout
       timeout > 0: wait up to this long to open the stream
       timeout = 0: wait forever
       timeout < 0: return immediately
    */
    int stay_in_loop = 1;
    int found_stream = 0;
    while (stay_in_loop) {
        log_debug("   rank %d: dart_get %s\n", ds->mpi_rank, ds_fname);
        err = adios_read_dataspaces_get(ds_fname, adios_byte, ds, 1, 0, offset, readsize, file_info_buf);
        if (!err) {
            stay_in_loop = 0;
            log_debug("-- %s, rank %d: data of '%s' exists\n", __func__, ds->mpi_rank, ds_fname);

            fp = ds_unpack_file_info (file_info_buf, 128, ds);
            if (fp) {
                found_stream = 1;
                fp->fh = (uint64_t) ds;
                fp->path = strdup(fname);
                log_debug("-- %s, rank %d: done fp=%x, fp->fh=%x\n", __func__, 
                        ds->mpi_rank, fp, fp->fh);
#ifndef DATASPACES_NO_VERSIONING
                file_versions.version[fidx]++; // next open will use new version
#endif

            } else {
                adios_error (err_no_memory, "Cannot allocate memory for file info.");
            }

        } else {

            /* Check if this file is finalized by the writer */
            char ds_cname[MAX_DS_NAMELEN];
            int coffset[] = {0,0,0}, creadsize[3] = {1,1,1};
            int cvalue; 
            int save_version = ds->access_version;
            snprintf(ds_cname, MAX_DS_NAMELEN, "CLOSED@%s",fname);
            ds->access_version = 0; // this variable is stored with version 0 only
            log_debug("   rank %d: dart_get %s\n", ds->mpi_rank, ds_cname);
            err = adios_read_dataspaces_get(ds_cname, adios_integer, ds, 1, 0, coffset, creadsize, &cvalue);
            ds->access_version = save_version; 
            if (!err) {
                /* This file was finalized */
                adios_error (err_end_of_stream, "End of file '%s'. Writer finalized this file.\n", fname);
                stay_in_loop = 0;
            } else {
                /* Not finalized, just no new timestep available */

                /* FIXME: implement code to distinguish if new versions are available.
                   i.e. the first version is gone and we need to open the first available.

                   fp->current_step should be the actual version opened
                   fp->last_step should be the last available step already in DataSpaces

                */


                log_debug ("Data of '%s' does not exist in DataSpaces (yet)\n", ds_fname);
                //adios_error (err_file_not_found_error, "Data of '%s' does not exist in DataSpaces\n", ds_fname);
            }
            log_debug("   rank %d: call dcg_unlock_on_read(%s)\n", ds->mpi_rank, fname);
            free(ds);
        }

        // check if we need to stay in loop 
        if (stay_in_loop) {
            if (timeout_msec < 0) 
                stay_in_loop = 0;
            else if (timeout_msec > 0 && (time_get()-t1 > timeout_msec*1000.0))
                stay_in_loop = 0;
        }

    } // while (stay_in_loop)

    if (!found_stream) {
        // release lock here if we have not opened it (otherwise release at close)
        dart_unlock_on_read(fname);
    }

    return fp;
}

int adios_read_dataspaces_close (ADIOS_FILE *fp) 
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    int i,j;

    adios_errno = 0;

    log_debug("-- %s, rank %d: fp=%x\n", __func__, ds->mpi_rank, fp);

    /* Release read lock locked in fopen */
    dart_unlock_on_read(fp->path);

    /* Disconnect from DATASPACES if we connected in fopen() */
    if (ds && ds->disconnect_at_close) {
        adios_read_dataspaces_finalize_method();
    }

    free_namelist ((ds->group_namelist),ds->ngroups);
    /* FIXME: free ds-> vars, attrs, groups and their stuff inside */
    free(ds);
    if (fp->path) { free(fp->path); fp->path = 0; }
    free(fp);
    return 0;
}

/* This function can be called if user places 
   the wrong sequences of dims for a var 
 */
void adios_read_dataspaces_reset_dimension_order (ADIOS_FILE *fp, int is_fortran)
{
    /* unimplemented */
}

ADIOS_GROUP * ds_unpack_group_info (char * buf, struct dataspaces_group_struct * group)
        
{
    ADIOS_GROUP *gp;
    char * b = buf;
    int i, j, k, blen, namelen, extrabyte;
    int datasize;
    struct dataspaces_var_struct * vars;
    struct dataspaces_attr_struct * attrs;
    int buf_len = group->group_index_len;
    uint64_t dims[3]; // all variables has 3 dimension values in the index 
    int didx[3]; // dimension reordering 
    

    if (!buf || buf_len < 12)
        return NULL;

    gp = (ADIOS_GROUP *) malloc(sizeof(ADIOS_GROUP));
    if (!gp) {
        adios_error (err_no_memory, "Could not allocate memory for group info");
        return NULL;
    }

    
    log_debug("   %s: buffer length = %d, content:\n", __func__, buf_len);
    for (i=0; i<buf_len; i+=16) {
        for (j=0; j<4; j++) {
            log_debug("%3.3d %3.3d %3.3d %3.3d    ", b[i+4*j], b[i+4*j+1], b[i+4*j+2], b[i+4*j+3]);
        }
        log_debug("\n");
    }
    

    blen = *(int*)b;  // buf len again from buffer itself
    if (blen != buf_len) {
        log_debug("WARNING: %s(): expected group info buffer length is %d but buffer head says %d\n", 
                    __func__, buf_len, blen);
    }
    b += sizeof(int); 
    gp->vars_count = *(int*)b; // number of variables
    b += sizeof(int);
    gp->attrs_count = *(int*)b; // number of attributes
    b += sizeof(int);

    log_debug("   %s(): vars count = %d, attrs count = %d\n", __func__, gp->vars_count, gp->attrs_count);
    gp->var_namelist = (char **) calloc (sizeof(char*), gp->vars_count);
    gp->attr_namelist = (char **) calloc (sizeof(char*),  gp->attrs_count);
    if (!gp->var_namelist || !gp->attr_namelist) {
        adios_error (err_no_memory, "Could not allocate space for variable/attribute names when opening a group");
        free(gp);
        if (gp->var_namelist) free(gp->var_namelist);
        return NULL;
    }

    vars = (struct dataspaces_var_struct *) 
              malloc (gp->vars_count * sizeof(struct dataspaces_var_struct));
    attrs = (struct dataspaces_attr_struct *) 
              malloc (gp->attrs_count * sizeof(struct dataspaces_attr_struct));

    if (!vars || !attrs) {
        adios_error (err_no_memory, "Could not allocate space for variable/attribute metadata");
        free(gp);
        free(gp->var_namelist);
        free(gp->attr_namelist);
        if (vars) free(vars);
        return NULL;
    }

    group->vars = vars;
    group->attrs = attrs;

    // extract each variable
    log_debug("    Extract variables\n");
    for (i=0;i<gp->vars_count;i++) {
        log_debug("      var %d, b = %d\n", i, b);
        namelen = *(int*)b; // lenght of name
        b += sizeof(int);
        log_debug("        namelen = %d, b = %d\n", namelen, b);
        gp->var_namelist[i] = (char *) malloc (namelen+1);
        if (!gp->var_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate space for variable name when opening a group"); 
            if (i>0) free_namelist ((gp->var_namelist),i);
            free(gp);
            return NULL;
        }
        else  {
            memcpy(gp->var_namelist[i], b, namelen);
            gp->var_namelist[i][namelen] = '\0';
        }
        b += namelen;
        vars[i].name = strdup(gp->var_namelist[i]);  
        log_debug("        name = %s, b = %d\n", vars[i].name, b);
        // type
        vars[i].type = *(enum ADIOS_DATATYPES*)b; 
        b += sizeof(int);
        log_debug("        type = %d, b = %d\n", (int)vars[i].type, b);
        // hastime
        vars[i].hastime = *(int*)b; 
        b += sizeof(int);
        log_debug("        hastime = %d, b = %d\n", vars[i].hastime, b);
        // dimensions
        vars[i].ndims = *(int*)b; 
        b += sizeof(int);
        log_debug("        ndims w/o time = %d, b = %d\n", vars[i].ndims, b);
        for (j=0; j < 3; j++) {
            dims[j] = *(uint64_t*)b; 
            b += 8;
            log_debug("          unordered dim[%d] = %lld, b = %d\n", j, dims[j], b);
        }
        // reorder DS dimensions to Fortran/C dimensions
        ds_dimension_ordering (vars[i].ndims, futils_is_called_from_fortran(), 
                               1 /*unpack*/, didx);
        // handle time too in this reordering
        k=0;
        if (vars[i].hastime && !futils_is_called_from_fortran()) {
            // C reader: time dimension is first. 
            k=1;
            vars[i].dims[0] = 1;
            log_debug("          dim[0] = %lld (time)\n", vars[i].dims[0], b);
        }
        for (j=0; j < vars[i].ndims; j++) {
            vars[i].dims[j+k] = dims[didx[j]];
            log_debug("          dim[%d] = %lld, b = %d\n", j+k, vars[i].dims[j+k], b);
        }
        if (vars[i].hastime && futils_is_called_from_fortran()) {
            // Fortran reader: time dimension is last
            vars[i].dims[ vars[i].ndims ] = 1;
            log_debug("          dim[%d] = %lld (time)", 
                    vars[i].ndims, vars[i].dims[vars[i].ndims]);
        }

        if (vars[i].hastime) {
            vars[i].ndims++;
            log_debug("        ndims = %d (with time)\n", vars[i].ndims);
        }
        if (vars[i].ndims == 0) { // && !vars[i].hastime) {
            // need to get scalar value too
            if (vars[i].type != adios_string) {
                datasize = common_read_type_size(vars[i].type, NULL);
                extrabyte=0;
            } else {
                memcpy (&datasize, b, sizeof(int));
                b += sizeof(int);
                extrabyte=1;
            }
            vars[i].value = (void *) malloc (datasize+extrabyte);
            if (vars[i].value) {
                memcpy (vars[i].value, b, datasize);
                if (vars[i].type == adios_string) 
                    ((char *)vars[i].value)[datasize] = '\0';
            } else {
                log_debug("ERROR: cannot allocate %d bytes to store the value of variable %s\n",
                        datasize, vars[i].name);
            }
            
            b += datasize;
            log_debug("        value read, b = %d\n", b);
        } else {
            vars[i].value = NULL;
        }
    }

    // extract each attribute
    for (i=0;i<gp->attrs_count;i++) {
        namelen = *(int*)b; // lenght of name
        b += sizeof(int);
        gp->attr_namelist[i] = (char *) malloc (namelen+1);
        if (!gp->attr_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate space for attribute name when opening a group"); 
            if (i>0) free_namelist ((gp->attr_namelist),i);
            free(gp);
            return NULL;
        }
        else  {
            memcpy(gp->attr_namelist[i], b, namelen);
            gp->attr_namelist[i][namelen] = '\0';
        }
        b += namelen;
        attrs[i].name = strdup(gp->var_namelist[i]);  
        // type
        attrs[i].type = *(enum ADIOS_DATATYPES*)b; 
        b += sizeof(int);
        // get attribute value 
        if (attrs[i].type != adios_string) {
            datasize = common_read_type_size(attrs[i].type, NULL);
            extrabyte=0;
        } else {
            memcpy (&datasize, b, sizeof(int));
            b += sizeof(int);
            extrabyte=1;
        }
        
        attrs[i].value = (void *) malloc (datasize+extrabyte);
        if (attrs[i].value) {
            memcpy (attrs[i].value, b, datasize);
            if (attrs[i].type == adios_string) 
                ((char*)attrs[i].value)[datasize] = '\0';
            log_debug("        value read, b = %d\n", b);
        } else {
            log_debug("ERROR: cannot allocate %d bytes to store the value of attribute %s\n",
                    datasize, attrs[i].name);
        }
        
        b += datasize;
    }
    return gp;
}

ADIOS_GROUP * adios_read_dataspaces_gopen (ADIOS_FILE *fp, const char * grpname)
{
    ADIOS_GROUP * gp;
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    int grpid;

    adios_errno = 0;
    /* check if grpname is found in the group name list */
    for (grpid=0; grpid < fp->groups_count; grpid++) {
        if (!strcmp (grpname, fp->group_namelist[grpid]))
            break;
    }
    if (grpid >= fp->groups_count) {
        adios_error (err_invalid_group, "Invalid group name %s for file %s", grpname, fp->path);
        return NULL;
    }

    /* Try to get group metadata from DataSpaces. If it does not exists, we get an error, which means
       the data does not exist. */
    int offset[] = {0,0,0}, readsize[3] = {ds->groups[0].group_index_len,1,1};
    char * group_info_buf = malloc (ds->groups[0].group_index_len);
    if (!group_info_buf) {
            adios_error (err_no_memory, "%s: Could not allocate buffer for group info buffer of %d bytes", 
                        __func__, ds->groups[0].group_index_len);
        return NULL;
    }
    int err;
    char ds_name[MAX_DS_NAMELEN];
    snprintf (ds_name, MAX_DS_NAMELEN, "GROUP@%s/%s",fp->path, grpname);
    log_debug("-- %s, rank %d: Get variable %s with size %d\n", __func__, ds->mpi_rank, ds_name, readsize[0]);
    err = adios_read_dataspaces_get(ds_name, adios_byte, ds, 1, 0, offset, readsize, group_info_buf);
    if (err) {
        adios_error (err_invalid_group, "Invalid group name %s for file %s. "
                     "Entity %s could not be retrieved from DataSpaces\n",
                     grpname, fp->path, ds_name);
        return NULL;
    } else {
        log_debug("-- %s, rank %d: data of '%s' exists\n", __func__, ds->mpi_rank, ds_name);
    }


    gp = ds_unpack_group_info (group_info_buf, &(ds->groups[grpid]));
    if (!gp) {
        log_debug("-- %s, rank %d: unpacking group index failed\n", __func__, ds->mpi_rank);
        return NULL;
    }

    /* fill out ADIOS_GROUP struct */
    gp->grpid = grpid;
    gp->gh = (uint64_t) 0;
    gp->fp = fp;
    gp->timestep = ds->access_version;
    gp->lasttimestep = ds->access_version;
    
    free (group_info_buf);
    return gp;
}

ADIOS_GROUP * adios_read_dataspaces_gopen_byid (ADIOS_FILE *fp, int grpid)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    if ( grpid < 0 || grpid > fp->groups_count-1) {
        adios_error (err_invalid_group, "Invalid group id %s. There are %d groups in %s", 
                     grpid, fp->groups_count, fp->path);
        return NULL;
    }
    return adios_read_dataspaces_gopen (fp, fp->group_namelist[grpid]);
}
                   
int adios_read_dataspaces_gclose (ADIOS_GROUP *gp)
{
    struct dataspaces_data_struct * ds = 
                (struct dataspaces_data_struct *) gp->fp->fh;
    struct dataspaces_var_struct * vars = 
                (struct dataspaces_var_struct *) ds->groups[gp->grpid].vars;
    struct dataspaces_attr_struct * attrs = 
                (struct dataspaces_attr_struct *) ds->groups[gp->grpid].attrs;
    int i;

    adios_errno = 0;

    free_namelist ((gp->var_namelist),gp->vars_count);
    free_namelist ((gp->attr_namelist),gp->attrs_count);
    for (i=0; i<gp->vars_count; i++) {
        free (vars[i].name);
        free (vars[i].value);
    }
    free(ds->groups[gp->grpid].vars);
    for (i=0; i<gp->attrs_count; i++) {
        free (attrs[i].name);
        free (attrs[i].value);
    }
    free(ds->groups[gp->grpid].attrs);
    free(gp);
    return 0;
}


/*
int adios_read_dataspaces_get_attr (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type,
                    int * size, void ** data)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) gp->fp->fh;
    struct dataspaces_attr_struct * attrs = ds->groups[gp->grpid].attrs;

    int attrid, vstartpos=0, fstartpos;

    // find names with or without beginning /
    if (attrname[0] == '/') vstartpos = 1;

    for (attrid=0; attrid < gp->attrs_count; attrid++) {
        if (attrs[attrid].name[0] == '/') vstartpos = 1;
        else fstartpos = 0;
        if (!strcmp (attrs[attrid].name+fstartpos, attrname+vstartpos))
            break; // found this name
    }

    if (attrid == gp->attrs_count) {
        adios_error (err_invalid_attrname, 
                     "attribute name '%s' not found among %d attributes of group '%s'!",
                     attrname, gp->attrs_count, gp->fp->group_namelist[gp->grpid]
                    );
        return adios_errno;
    }

    return adios_read_dataspaces_get_attr_byid(gp, attrid, type, size, data);
}
*/

int adios_read_dataspaces_get_attr_byid (ADIOS_FILE * fp, int attrid, 
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    struct dataspaces_attr_struct * attrs = ds->groups[0].attrs;

    if (attrid < 0 || attrid > fp->nattrs) {
        adios_error (err_invalid_attrid, 
                     "File %s has %d attributes. Invalid attrid %d",
                     fp->path, fp->nattrs, attrid);
        return adios_errno;
    }

    adios_errno = 0;
    *type = attrs[attrid].type;
    *size = common_read_type_size(*type, attrs[attrid].value);
    *data = (void *) malloc (*size);
    if (*data) {
        memcpy (*data, attrs[attrid].value, *size);
    } else {
        adios_error (err_no_memory, "Could not allocate memory for variable info.");
        return adios_errno;
    }
    return 0; 
}

/** Search name in list with or without leading /, return position or -1 if not found */
static int adios_read_dataspaces_find_varname (int nvars, struct dataspaces_var_struct * vars, const char * name) 
{
    int id, vstartpos=0, fstartpos;

    // find names with or without beginning /
    if (name[0] == '/') vstartpos = 1;

    log_debug("   find var %s, startpos=%d\n", name, vstartpos);

    for (id=0; id < nvars; id++) {
        if (vars[id].name[0] == '/') fstartpos = 1;
        else fstartpos = 0;
        log_debug("     check %s, startpos=%d\n", vars[id].name, fstartpos);
        if (!strcmp (vars[id].name+fstartpos, name+vstartpos))
            break; // found this name
    }

    if (id == nvars) {
        return -1;
    }
    return id;
}

ADIOS_VARINFO * adios_read_dataspaces_inq_var (ADIOS_GROUP *gp, const char * varname) 
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) gp->fp->fh;
    struct dataspaces_var_struct * vars = ds->groups[gp->grpid].vars;
    int varid; 

    varid = adios_read_dataspaces_find_varname (gp->vars_count, vars, varname);
    if (varid == -1) {
        adios_error (err_invalid_varname, 
                     "Variable name '%s' not found among %d variables of group '%s'!",
                     varname, gp->vars_count, gp->fp->group_namelist[gp->grpid]
                    );
        return NULL;
    }

    return adios_read_dataspaces_inq_var_byid(gp, varid);
}

ADIOS_VARINFO * adios_read_dataspaces_inq_var_byid (ADIOS_GROUP *gp, int varid)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) gp->fp->fh;
    struct dataspaces_var_struct * vars = ds->groups[gp->grpid].vars;
    ADIOS_VARINFO * vi;
    int i,k;

    if (varid < 0 || varid > gp->vars_count) {
        adios_error (err_invalid_varid, "Group %s has %d variables. Invalid varid %d",
                    gp->fp->group_namelist[gp->grpid], gp->vars_count, varid);
        return NULL;
    }

    adios_errno = 0;
    vi = (ADIOS_VARINFO *) malloc(sizeof(ADIOS_VARINFO));
    if (!vi) {
        adios_error (err_no_memory, "Could not allocate memory for variable info.");
        return NULL;
    }

    vi->varid = varid;
    vi->type = vars[varid].type;
    vi->ndim = vars[varid].ndims;
    vi->dims = vars[varid].dims;
    if (vars[varid].hastime)
        vi->timedim = 0;
    else
        vi->timedim = -1;
    vi->value = vars[varid].value;
    vi->gmin = NULL;
    vi->gmax = NULL;
    vi->mins = NULL;
    vi->maxs = NULL;
    
    return vi;
}

void adios_read_dataspaces_free_varinfo (ADIOS_VARINFO *vp)
{
    if (vp) {
        //if (vp->dims)   free(vp->dims);  // link to vars[varid].dims
        //if (vp->value)  free(vp->value); // link to vars[varid].value
        if (vp->gmin && vp->gmin != vp->value)   free(vp->gmin);
        if (vp->gmax && vp->gmax != vp->value)   free(vp->gmax);
        if (vp->mins)   free(vp->mins);
        if (vp->maxs)   free(vp->maxs);
        free(vp);
    }
}

static int adios_read_dataspaces_get (const char * varname, enum ADIOS_DATATYPES vartype, 
                                struct dataspaces_data_struct * ds, 
                                int ndims, int is_fortran_ordering, 
                                int * offset, int * readsize, void * data)
{

    struct obj_data *od;
    int elemsize = common_read_type_size(vartype, NULL);
    int i, err;
    int didx[3];
    int lb[3] = {0,0,0};
    int ub[3] = {0,0,0};

    // reorder DS dimensions to Fortran/C dimensions
    ds_dimension_ordering (ndims, is_fortran_ordering, 0 /*pack*/, didx);
    for (i=0; i<3; i++) {
        lb[i] = offset[didx[i]];
        ub[i] = offset[didx[i]]+readsize[didx[i]]-1;
    }

    log_debug("-- %s, rank %d: get data: varname=%s version=%d, lb=(%d,%d,%d) ub=(%d,%d,%d)}\n",
        __func__, ds->mpi_rank, varname, ds->access_version, lb[0], lb[1], lb[2], 
        ub[0], ub[1], ub[2]);

    err =  dart_get (varname, ds->access_version, elemsize, 
                     lb[0], lb[1], lb[2],
                     ub[0], ub[1], ub[2],
                     data
                    );
    /*if (err == -ENOMEM) {
        adios_error (err_no_memory, "Not enough memory for DATASPACES to perform dart_get()");  
        return -err_no_memory;
    } 
    else*/ if (err) {
        adios_error (err_corrupted_variable, "DATASPACES failed to read variable %s.", varname);  
        return -err_corrupted_variable;
    }

    return 0;
}

int64_t adios_read_dataspaces_read_var (ADIOS_GROUP * gp, const char * varname,
                        const uint64_t * start, const uint64_t * count,
                        void * data)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) gp->fp->fh;
    struct dataspaces_var_struct * vars = ds->groups[gp->grpid].vars;
    int varid; 

    varid = adios_read_dataspaces_find_varname (gp->vars_count, vars, varname);
    if (varid == -1) {
        adios_error (err_invalid_varname,
                "Variable name '%s' not found among %d variables of group '%s'!",
                varname, gp->vars_count, gp->fp->group_namelist[gp->grpid]
                );
        return -adios_errno;
    }

    return adios_read_dataspaces_read_var_byid(gp, varid, start, count, data);
}


int64_t adios_read_dataspaces_read_var_byid (ADIOS_GROUP    * gp,
                             int              varid,
                             const uint64_t  * start,
                             const uint64_t  * count,
                             void           * data)
{
    int64_t total_size;
    int offset[3] = {0,0,0};
    int readsize[3] = {1,1,1};
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) gp->fp->fh;
    struct dataspaces_var_struct * vars = ds->groups[gp->grpid].vars;
    int elemsize;
    int err;
    int i,k,ndims,tidx;
    char ds_name[MAX_DS_NAMELEN];

    if (varid < 0 || varid > gp->vars_count) {
        adios_error (err_invalid_varid, "Group %s has %d variables. Invalid varid %d",
                    gp->fp->group_namelist[gp->grpid], gp->vars_count, varid);
        return -adios_errno;
    }

    elemsize = common_read_type_size(vars[varid].type, vars[varid].value);

    // handle scalars first (no need to read from space again)
    if (vars[varid].ndims == 0) { // && !vars[varid].hastime) {
        if (data) {
            memcpy (data, vars[varid].value, elemsize);
        } else {
            adios_error (err_no_memory, "adios_read_var() expects an allocated array to store data.");
            return -adios_errno;
        }
        return elemsize;
    }
        
    adios_errno = 0;
    k=0;
    ndims = vars[varid].ndims;
    /* DATASPACES uses integers for boundaries */
    if (vars[varid].hastime) {
        /* DATASPACES has no time dimensions stored in space. Remove from dims */
        ndims = vars[varid].ndims - 1;
        if (futils_is_called_from_fortran()) {
            // Fortran: time dim is last
            k = 0;
            tidx = ndims;
        } else {
            // C: time dim is first
            k = 1;
            tidx = 0;
        }
        if (start[tidx] != gp->timestep || count[tidx] != 1) {
            adios_error (err_out_of_bound, 
                         "offset/readsize is out of bound in time dimension %d for variable %s\n"
                         "you can read 1 element from offset (=timestep) %lld\n" 
                         "you provided start=%lld, count=%lld\n",
                         tidx, vars[varid].name, gp->timestep, start[tidx], count[tidx]
                         );
            return -adios_errno;
        }
    } else {
        k = 0;
    }
    total_size = elemsize;
    for (i=0; i<ndims; i++) {
        if ( start[i+k] < 0 || start[i+k]+count[i+k] > vars[varid].dims[i+k]) {
            adios_error (err_out_of_bound, 
                         "offset/readsize is out of bound in dimension %d for variable %s\n"
                         "size of dimension %d = %lld; you provided start=%lld, count=%lld",
                         i+k, vars[varid].name, i+k, vars[varid].dims[i+k], start[i+k], count[i+k]
                         );
            return -adios_errno;
        }
        offset[i]    = (int) start[i+k];
        readsize[i]  = (int) count[i+k];
        total_size   = total_size * count[i+k];
    }

    snprintf(ds_name, MAX_DS_NAMELEN, "%s/%s/%s", 
             fp->path, gp->fp->group_namelist[gp->grpid], vars[varid].name);
    log_debug("-- %s, rank %d: get data: varname=%s start=(%lld,%lld,%lld,%lld) count=(%lld,%lld,%lld,%lld)}\n",
        __func__, ds->mpi_rank, ds_name, start[0], start[1], start[2], start[3], count[0], count[1], count[2], count[3]);
    log_debug("-- %s, rank %d: get data: varname=%s offset=(%d,%d,%d) readsize=(%d,%d,%d)}\n",
        __func__, ds->mpi_rank, ds_name, offset[0], offset[1], offset[2], readsize[0], readsize[1], readsize[2]);

    err = adios_read_dataspaces_get(ds_name, vars[varid].type, ds, 
                              ndims, futils_is_called_from_fortran(),
                              offset, readsize, data);
    if (err)
        return err;

    return total_size;
}


/* Tell the DataSpaces order of dimensions for a 1-3 dim array written from Fortran or C.
   unpack=1: the reverse of packing (to retrieve the original order).
   didx should be an int [3] array in any case.
*/
void ds_dimension_ordering(int ndims, int is_app_fortran, int unpack, int *didx)
{
    /* Order of dimensions: in DataSpaces: slow, fast, slowest
       Fortran: i,j,k --> j, i, k  = lb[1], lb[0], lb[2]
                i,j   --> j, i     = lb[1], lb[0], lb[2]=1
                i     --> 1, i     = lb[1]=1, lb[0], lb[2]=1
       C:       i,j,k --> j, k, i  = lb[1], lb[2], lb[0]
                i,j   --> i, j     = lb[0], lb[1], lb[2]=1
                i     --> 1, i     = lb[1]=1, lb[0], lb[2]=1 (same as Fortran)

       unpack: C(i,j,k) ordering applied twice does not result in the original order
               so we need to have a reverse mapping for this case.
               For all the other cases, applying twice results in the same order
               even for packing from Fortran and unpacking to C, and vice versa.
               F(i,j,k) -(pack)-> DS(j,i,k) -(unpack)-> C(k,j,i) or F(i,j,k)
               C(i,j,k) -(pack)-> DS(j,k,i) -(unpack)-> C(i,j,k) or F(k,j,i)
               F(i,j)   -(pack)-> DS(j,i)   -(unpack)-> C(j,i)   or F(i,j)
               C(i,j)   -(pack)-> DS(i,j)   -(unpack)-> C(i,j)   or F(j,i)
               F(i)     -(pack)-> DS(1,i)   -(unpack)-> C(i,(1)) or F(i,(1))
               C(i)     -(pack)-> DS(1,i)   -(unpack)-> C(i,(1)) or F(i,(1))
    */

    if (ndims == 0) {
        didx[0] = 0; didx[1] = 1; didx[2] = 2;
    } else if (is_app_fortran || ndims == 1) {
        /* Flip 1st and 2nd dimension for DataSpaces representation
           for any Fortran writings and for any 1D array :
           Fortran: i,j,k --> j, i, k  = lb[1], lb[0], lb[2]
           C:       i     --> 1, i     = lb[1]=1, lb[0], lb[2]=1 
        */
        didx[0] = 1; didx[1] = 0; didx[2] = 2;
    } else if (ndims == 2) {
        /* C: i,j   --> i, j     = lb[0], lb[1], lb[2]=1 */
        didx[0] = 0; didx[1] = 1; didx[2] = 2;
    } else { // (ndims == 3) 
        if (!unpack) {
            /* C: i,j,k --> j, k, i  = lb[1], lb[2], lb[0] */
            didx[0] = 1; didx[1] = 2; didx[2] = 0;
        } else {
            /* DataSpaces x,y,z --> z,x,y  (ijk->kij, or jki->ijk) */
            didx[0] = 2; didx[1] = 0; didx[2] = 1;
        }
    }
}


int64_t adios_read_dataspaces_read_local_var (ADIOS_GROUP * gp, const char * varname,
                                      int vidx, const uint64_t * start,
                                      const uint64_t * count, void * data)
{  
    adios_error (err_operation_not_supported, "adios_read_local_var() is not supported with DATASPACES method.");
    return -adios_errno;
}
