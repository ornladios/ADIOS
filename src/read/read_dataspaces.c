/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */


/**************************************************/
/* Read method for DATASPACES memory-to-memory coupling */
/**************************************************/

#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <errno.h>  /* errno */
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
#include "core/common_read.h" // common_read_selection_* functions

#include "core/transforms/adios_transforms_common.h" // NCSU ALACRITY-ADIOS

#include "dataspaces.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define MYFREE(p) {free(p); p=NULL;}

void adios_read_dataspaces_release_step (ADIOS_FILE *fp);

/*
#include <time.h> // nanosleep
static void adios_nanosleep (int sec, int nanosec) 
{
#if HAVE_NANOSLEEP
    struct timespec treq = {.tv_sec=sec, .tv_nsec=nanosec}; 
    struct timespec trem;
    int r;
    r = nanosleep(&treq, &trem);
    while (r == -1 && errno == EINTR) {
        treq.tv_sec = trem.tv_sec;
        treq.tv_nsec = trem.tv_nsec;
        r = nanosleep (&trem, &treq);
    }
#else
    if (sec>0)
        sleep(sec);
    else 
        sleep(1);

#endif
}

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
*/

#define MAX_DS_NAMELEN 128
/* Maximum number of different filenames allowed per process during the whole run */
#define MAXNFILE 20 
/*#define DATASPACES_NO_VERSIONING   define it at configure as -DDATASPACES_NO_VERSIONING in CFLAGS */
/* Length of the 1D array representing the file in DataSpaces (contains metadata) */
#define FILEINFO_BUFLEN 128

static int chunk_buffer_size = 1024*1024*16; // 16MB default size for reading in data in chunking mode
static char *chunk_buffer = 0;

static int poll_interval_msec = 10; // how much to wait between polls when timeout is used
static int enable_read_meta_collective = 1; // when enabled, meta data reading becomes collective. One reader process would fetch meta data from DataSpaces and broadcast to other procseses using MPI_Bcast

struct dataspaces_fileversions_struct { // current opened version of each stream/file
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
    uint64_t               dims[MAX_DS_NDIM]; // we have max MAX_DS_NDIM dims in DataSpaces
    void                 * value;
};

struct dataspaces_attr_struct { // describes one attribute (of one group)
    char                 * name;
    enum ADIOS_DATATYPES   type;
    void                 * value;
};

/*struct dataspaces_group_struct { // accessible as fp->fh->groups[grpid]
    int nvars;                   // number of vars in this group
    int nattrs;                  // number of attrs in this group
};*/

struct dataspaces_data_struct { // accessible as fp->fh
    int current_step;           // counting the access
    int disconnect_at_close;    // disconnect from DATASPACES in fclose()
    MPI_Comm comm;              // communicator saved for lock/unlock operations
    int mpi_rank;               // rank of this process
    int nproc;                  // number of processes opening this stream
    int locked_fname;           // 1: locked 'fname' in DATASPACES, 0: unlocked
    int freed_mem;              // 1: freed all memory used for internal data
    int file_index;             // index to file_versions[] array
    enum ADIOS_LOCKMODE lock_mode; // locking mode requested by user
    struct dataspaces_var_struct  * vars;  // number of vars is ADIOS_FILE->nvars
    struct dataspaces_attr_struct * attrs; // number of attrs is ADIOS_FILE->nattrs
    /* Group info */
    int group_index_len;         // length of group index in GROUP@fn/gn variable 
    char *group_name;            // name of the group
    /* Read requests */
    read_request * req_list;     // list of scheduled requests
    read_request * req_list_tail; // tail of list of scheduled requests (to speed up insert)
    int nreq;                    // number of scheduled requests
    /* single chunk */
    ADIOS_VARCHUNK *chunk;       // the single chunk to store the last read and serve the user
};

// Declarations
static int adios_read_dataspaces_get (const char * varname, enum ADIOS_DATATYPES vartype, 
                                int version, int rank, 
                                int ndims, int is_fortran_ordering, 
                                uint64_t * offset, uint64_t * readsize, uint64_t * dims, void * data);
static int adios_read_dataspaces_get_meta(const char * varname, enum ADIOS_DATATYPES vartype,
                                int version, int rank,
                                int ndims, int is_fortran_ordering,
                                uint64_t * offset, uint64_t * readsize, void * data);
static int adios_read_dataspaces_get_meta_collective(const char * varname,
                                enum ADIOS_DATATYPES vartype,
                                int version, int rank,
                                int ndims, int is_fortran_ordering,
                                uint64_t * offset, uint64_t * readsize, void * data, MPI_Comm comm);

static char* get_chunk_buffer()
{
    if (!chunk_buffer) {
        chunk_buffer = (char *) malloc (chunk_buffer_size);
        if (!chunk_buffer) {
            adios_error (err_no_memory, 
                    "Could not allocate chunk buffer of size %dMB\n", chunk_buffer_size/1024/1024);
        }
    }
    return chunk_buffer;
}

static int get_meta(const char * varname,
                    enum ADIOS_DATATYPES vartype,
                    int version, int rank,
                    int ndims, int is_fortran_ordering,
                    uint64_t * offset, uint64_t * readsize, void * data, MPI_Comm comm)
{
    if (enable_read_meta_collective) {
        return adios_read_dataspaces_get_meta_collective(varname, vartype, version,
            rank, ndims, is_fortran_ordering, offset, readsize, data, comm);
    } else {
        return adios_read_dataspaces_get_meta(varname, vartype, version, rank, ndims,
            is_fortran_ordering, offset, readsize, data);
    }
}

/* If init is used, we connect to DATASPACES here, otherwise we connect in fopen.
   If multiple fopen..fclose cycles are used, init/finalize must be used too to
   avoid multiple connection/disconnection in fopen/fclose.
*/
int adios_read_dataspaces_init_method (MPI_Comm comm, PairStruct * params) 
{ 
    int  nproc, drank, dpeers;
    int  rank, err;
    int  appid, max_chunk_size, pollinterval, was_set;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nproc);

    PairStruct *p = params;
    while (p) {
        if (!strcasecmp (p->name, "app_id")) {
            errno = 0;
            appid = strtol(p->value, NULL, 10);
            if (appid > 0 && !errno) {
                log_debug ("App ID parameter set to %d for DATASPACES read method\n", 
                            appid);
                globals_adios_set_application_id (appid);
            } else {
                log_error ("Invalid 'app_id' parameter given to the DATASPACES read "
                           "method: '%s'\n", p->value);
            }
        } else if (!strcasecmp (p->name, "max_chunk_size")) {
            errno = 0;
            max_chunk_size = strtol(p->value, NULL, 10);
            if (max_chunk_size > 0 && !errno) {
                log_debug ("max_chunk_size set to %dMB for DATASPACES read method\n", 
                            max_chunk_size);
                chunk_buffer_size = max_chunk_size * 1024 * 1024;
            } else {
                log_error ("Invalid 'max_chunk_size' parameter given to the DATASPACES "
                            "read method: '%s'\n", p->value);
            }
        } else if (!strcasecmp (p->name, "poll_interval")) {
            errno = 0;
            pollinterval = strtol(p->value, NULL, 10);
            if (pollinterval > 0 && !errno) {
                log_debug ("poll_interval set to %d millisecs for DATASPACES read method\n", 
                            pollinterval);
                poll_interval_msec = pollinterval;
            } else {
                log_error ("Invalid 'poll_interval' parameter given to the DATASPACES "
                            "read method: '%s'\n", p->value);
            }
        } else if (!strcasecmp (p->name, "disable_collective_read_meta")) {
            errno = 0;
            enable_read_meta_collective = 0;
            log_debug("Set 'disable_collective_read_meta' for DATASPACES read method\n");
        } else { 
            log_error ("Parameter name %s is not recognized by the DATASPACES read "
                        "method\n", p->name);
        }
        p = p->next;
    }


    /* Connect to DATASPACES, but only if we are not yet connected (from Write API) */
    if (!globals_adios_is_dataspaces_connected()) {
        appid = globals_adios_get_application_id (&was_set);
        if (!was_set) 
            appid = 2;
        log_debug("-- %s, rank %d: connect to dataspaces with nproc=%d and appid=%d\n", 
                    __func__, rank, nproc, appid);
        err = dspaces_init(nproc, appid, &comm, NULL);
        if (err < 0) {
            adios_error (err_connection_failed, "Failed to connect with DATASPACES\n");
            return err_connection_failed;
        }

        //drank = dspaces_rank(dcg);
        //dpeers = dspaces_peers(dcg);
    }
    globals_adios_set_dataspaces_connected_from_reader();
    n_filenames = 0;
    log_info("Connected to DATASPACES\n");
    return 0; 
}

int adios_read_dataspaces_finalize_method () 
{ 
    // disconnect from DATASPACES only if we the reader is connected (the writer not anymore)
    if (globals_adios_is_dataspaces_connected_from_reader() && 
        !globals_adios_is_dataspaces_connected_from_both()) 
    {
        dspaces_finalize();
        log_info("Disconnected from DATASPACES\n");
    }
    globals_adios_set_dataspaces_disconnected_from_reader();
    if (chunk_buffer) {
        free (chunk_buffer); 
        chunk_buffer = NULL;
    }
}

static int ds_unpack_file_info (ADIOS_FILE *fp, char * buf, int buf_len)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    char * b = buf;
    int blen, glen, i;

    if (!buf || buf_len < 21)
        return 1;

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
    ds->group_index_len = *(int*)b; // length of group index
    b += sizeof(int);
    glen = *(int*)b; // length of (only) group name
    b += sizeof(int);

    fp->file_size = 0;
    fp->version = 1;
    fp->endianness = 0; // FIXME: not always Little Endian. Does it matter? 
    ds->group_name = (char *) malloc (strlen(b)+1);
    if (!ds->group_name) {
        adios_error (err_no_memory, 
                "Could not allocate buffer for group name in adios_read_open()\n");
        free(fp);
        return 1;
    }
    else  {
        strcpy(ds->group_name,b);
    }

    return 0;
}


static int ds_unpack_group_info (ADIOS_FILE *fp, char * buf)
        
{
    char * b = buf;
    int i, j, k, blen, namelen, extrabyte;
    int datasize;
    struct dataspaces_var_struct * vars;
    struct dataspaces_attr_struct * attrs;
    uint64_t dims[MAX_DS_NDIM]; // all variables has MAX_DS_NDIM dimension values in the index 
    int didx[MAX_DS_NDIM]; // dimension reordering 
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    int buf_len = ds->group_index_len;
    

    if (!buf || buf_len < 12)
        return 1;

    
    log_debug("   %s: buffer length = %d, content:\n", __func__, buf_len);
    for (i=0; i<buf_len; i+=16) {
        for (j=0; j<4; j++) {
            log_debug_cont("%3.3hhu %3.3hhu %3.3hhu %3.3hhu    ", 
                    b[i+4*j], b[i+4*j+1], b[i+4*j+2], b[i+4*j+3]);
        }
        log_debug_cont("\n");
    }
    

    blen = *(int*)b;  // buf len again from buffer itself
    if (blen != buf_len) {
        log_debug("WARNING: %s(): expected group info buffer length is %d but buffer head says %d\n", 
                    __func__, buf_len, blen);
    }
    b += sizeof(int); 
    int vars_count = *(int*)b; // number of variables (should be = fp->nvars)
    b += sizeof(int);
    if (vars_count != fp->nvars) {
        log_debug("WARNING: %s(): expected number of variables in group %s is %d "
                  "but group metadata says %d\n", 
                  __func__, ds->group_name, fp->nvars, vars_count);
    }
    int attrs_count = *(int*)b; // number of attributes (= fp->nattrs)
    b += sizeof(int);
    if (attrs_count != fp->nattrs) {
        log_debug("WARNING: %s(): expected number of attributes in group %s is %d "
                  "but group metadata says %d\n", 
                  __func__, ds->group_name, fp->nattrs, attrs_count);
    }

    log_debug("   %s(): vars count = %d, attrs count = %d\n", __func__, vars_count, attrs_count);
    fp->var_namelist = (char **) calloc (sizeof(char*), vars_count);
    fp->attr_namelist = (char **) calloc (sizeof(char*),  attrs_count);
    if (!fp->var_namelist || !fp->attr_namelist) {
        adios_error (err_no_memory, "Could not allocate space for variable/attribute names when opening a group\n");
        if (fp->var_namelist) free(fp->var_namelist);
        return 1;
    }

    vars = (struct dataspaces_var_struct *) 
              malloc (vars_count * sizeof(struct dataspaces_var_struct));
    attrs = (struct dataspaces_attr_struct *) 
              malloc (attrs_count * sizeof(struct dataspaces_attr_struct));

    if (!vars || !attrs) {
        adios_error (err_no_memory, "Could not allocate space for variable/attribute metadata\n");
        free(fp->var_namelist);
        free(fp->attr_namelist);
        if (vars) free(vars);
        return 1;
    }

    ds->vars = vars;
    ds->attrs = attrs;

    // extract each variable
    log_debug("    Extract variables\n");
    for (i=0;i<vars_count;i++) {
        log_debug("      var %d, b = %d\n", i, b);
        namelen = *(int*)b; // lenght of name
        b += sizeof(int);
        log_debug("        namelen = %d, b = %d\n", namelen, b);
        fp->var_namelist[i] = (char *) malloc (namelen+1);
        if (!fp->var_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate space for variable name when opening a group\n"); 
            if (i>0) a2s_free_namelist ((fp->var_namelist),i);
            return 1;
        }
        else  {
            memcpy(fp->var_namelist[i], b, namelen);
            fp->var_namelist[i][namelen] = '\0';
        }
        b += namelen;
        vars[i].name = strdup(fp->var_namelist[i]);  
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
        for (j=0; j < MAX_DS_NDIM; j++) {
            dims[j] = *(uint64_t*)b; 
            b += 8;
            log_debug("          unordered dim[%d] = %lld, b = %d\n", j, dims[j], b);
        }
        // reorder DS dimensions to Fortran/C dimensions
        ds_dimension_ordering (vars[i].ndims, futils_is_called_from_fortran(), 
                               1 /*unpack*/, didx);
        for (j=0; j < vars[i].ndims; j++) {
            vars[i].dims[j] = dims[didx[j]];
            log_debug("          dim[%d] = %lld, b = %d\n", j, vars[i].dims[j], b);
        }

        if (vars[i].ndims == 0 && !vars[i].hastime) {
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
                log_error("Cannot allocate %d bytes to store the value of variable %s\n",
                        datasize, vars[i].name);
                return 1;
            }
            
            b += datasize;
            log_debug("        value read, b = %d\n", b);
        } else {
            vars[i].value = NULL;
        }
    }

    // extract each attribute
    for (i=0;i<attrs_count;i++) {
        namelen = *(int*)b; // lenght of name
        b += sizeof(int);
        fp->attr_namelist[i] = (char *) malloc (namelen+1);
        if (!fp->attr_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate space for attribute name when opening a group\n"); 
            if (i>0) a2s_free_namelist ((fp->attr_namelist),i);
            return 1;
        }
        else  {
            memcpy(fp->attr_namelist[i], b, namelen);
            fp->attr_namelist[i][namelen] = '\0';
        }
        b += namelen;
        attrs[i].name = strdup(fp->attr_namelist[i]);  
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
            log_error("Cannot allocate %d bytes to store the value of attribute %s\n",
                    datasize, attrs[i].name);
            return 1;
        }
        
        b += datasize;
    }
    return 0;
}

static int get_groupdata (ADIOS_FILE *fp)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;

    /* Try to get group metadata from DataSpaces. If it does not exists, we get an error, which means
       the data does not exist. */
    uint64_t offset[MAX_DS_NDIM], readsize[MAX_DS_NDIM];
    offset[0] = 0; readsize[0] = ds->group_index_len;
    int ndim = 1;
    char * group_info_buf = malloc (ds->group_index_len);
    if (!group_info_buf) {
            adios_error (err_no_memory, 
                    "%s: Could not allocate buffer for group info buffer of %d bytes\n", 
                    __func__, ds->group_index_len);
        return 1;
    }
    int err;
    char ds_name[MAX_DS_NAMELEN];
    snprintf (ds_name, MAX_DS_NAMELEN, "GROUP@%s/%s",fp->path, ds->group_name);
    log_debug("-- %s, rank %d: Get variable %s with size %lld\n", __func__, 
              ds->mpi_rank, ds_name, readsize[0]);
    err = get_meta (ds_name, adios_byte, ds->current_step, ds->mpi_rank, 
                     ndim, 0, offset, readsize, group_info_buf, ds->comm);
    if (err) {
        adios_error (err_invalid_group, "Invalid group name %s for file %s. "
                     "Entity %s could not be retrieved from DataSpaces\n",
                     ds->group_name, fp->path, ds_name);
        return 1;
    } else {
        log_debug("-- %s, rank %d: data of '%s' exists\n", __func__, ds->mpi_rank, ds_name);
    }


    err = ds_unpack_group_info (fp, group_info_buf);
    if (err) {
        log_debug("-- %s, rank %d: unpacking group index failed\n", __func__, ds->mpi_rank);
        return 1;
    }

    free (group_info_buf);
    return 0;
}

static void lock_file (ADIOS_FILE *fp, struct dataspaces_data_struct *ds)
{
    if (!ds->locked_fname) {
        log_debug("   rank %d: call dspaces_lock_on_read(%s)\n", ds->mpi_rank, fp->path);
        dspaces_lock_on_read(fp->path, &ds->comm);
        ds->locked_fname = 1;
    //} else {
    //    log_error("   rank %d: lock_file called with the lock already in place\n", ds->mpi_rank);
    }
}

static void unlock_file (ADIOS_FILE *fp, struct dataspaces_data_struct *ds)
{
    if (ds->locked_fname) {
        log_debug("   rank %d: call dspaces_unlock_on_read(%s)\n", ds->mpi_rank, fp->path);
        dspaces_unlock_on_read(fp->path, &ds->comm);
        ds->locked_fname = 0;
    //} else {
    //    log_error("   rank %d: unlock_file called with the lock already released\n", ds->mpi_rank);
    }
}


enum WHICH_VERSION {NEXT_VERSION, NEXT_AVAILABLE_VERSION, LAST_VERSION};
static const char * which_version_str[3] = {"current", "next available", "last"};
enum STEP_STATUS {STEP_OK, STEP_STREAMNOTFOUND, STEP_STEPNOTREADY, STEP_STEPDISAPPEARED, STEP_STREAMTERMINATED, STEP_OTHERERROR};

static int get_step (ADIOS_FILE *fp, int step, enum WHICH_VERSION which_version, float timeout_sec)
{
    /* Try to get variable with fname. If it does not exists, we get an error, which means
       the data does not exist. So we return an error just like with real files */
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    uint64_t offset[MAX_DS_NDIM], readsize[MAX_DS_NDIM];
    int ndim;
    char file_info_buf[FILEINFO_BUFLEN];
    int version_info_buf[2]; // 0: last version, 1: terminated?
    int err, i;
    char ds_vname[MAX_DS_NAMELEN];
    char ds_fname[MAX_DS_NAMELEN];
    double t1 = adios_gettime();
    enum STEP_STATUS step_status = STEP_OK;

    snprintf(ds_vname, MAX_DS_NAMELEN, "VERSION@%s",fp->path);
    snprintf(ds_fname, MAX_DS_NAMELEN, "FILE@%s",fp->path);
    //log_debug("-- %s, rank %d: Get variable %s\n", __func__, ds->mpi_rank, ds_fname);
    ds->freed_mem = 0;

    /* While loop for handling timeout
       timeout >= 0: wait up to this long to open the stream
       timeout <  0: wait forever
    */
    int stay_in_poll_loop = 1;
    int found_stream = 0;
    int nversions, *versions;
    while (stay_in_poll_loop) {
        lock_file (fp, ds);
        step_status = STEP_OK;

        log_debug("   rank %d: dspaces_get %s\n", ds->mpi_rank, ds_vname);
        ndim = 1;
        offset[0] = 0; readsize[0] = 2; //*sizeof(int); // VERSION%name is 2 integers only
        err = get_meta (ds_vname, adios_integer, 0, ds->mpi_rank, ndim, 0, 
                         offset, readsize, version_info_buf, ds->comm);

        if (!err) {
            int last_version = version_info_buf[0];
            int terminated = version_info_buf[1];
            log_debug("   rank %d: version info: last=%d, terminated=%d\n", 
                      ds->mpi_rank, last_version, terminated);

            if (last_version < step) {
                // we have no more new steps
                if (terminated) {
                    // stream is gone, we read everything 
                    step_status = STEP_STREAMTERMINATED;
                    stay_in_poll_loop = 0;
                } else {
                    // a next step may come 
                    step_status = STEP_STEPNOTREADY;
                    // we may stay in poll loop
                }
            } else {
                // Try to get the version the user wants
                if (which_version == LAST_VERSION)
                    step = last_version;
                readsize[0] = FILEINFO_BUFLEN; // FILE%name is FILEINFO_BUFLEN bytes long

                int max_check_version = last_version;
                if (which_version == NEXT_VERSION) 
                    max_check_version = step;

                // Loop until we find what we need or go past the last version
                do {
                    log_debug("   rank %d: dspaces_get %s\n", ds->mpi_rank, ds_fname);
                    err = get_meta(ds_fname, adios_byte, step, ds->mpi_rank, 
                                 ndim, 0, offset, readsize, file_info_buf, ds->comm);
                    step++; // value will go over the target with 1
                } while (err && step <= max_check_version);

                if (!err) {
                    /* Found object with this access version */
                    step--; // undo the last increment above
                    ds->current_step = step;
                    stay_in_poll_loop = 0;
                    log_debug("   rank %d: step %d of '%s' exists\n", 
                            ds->mpi_rank, ds->current_step, ds_fname);

                    err = ds_unpack_file_info (fp, file_info_buf, FILEINFO_BUFLEN);
                    if (!err) {
                        found_stream = 1;
                        fp->current_step = ds->current_step;
                        fp->last_step = last_version;

                        /* Get the variables and attributes the (only) group separately */
                        err = get_groupdata (fp);
                        if (err) {
                            // something went wrong with the group(s)
                            step_status = STEP_OTHERERROR;
                        }
                    } else {
                        // something went wrong with the file metadata
                        step_status = STEP_OTHERERROR;
                    }

                } else {
                    if (which_version == NEXT_VERSION) 
                    {
                        if (step < last_version) {
                            step_status = STEP_STEPDISAPPEARED;
                            stay_in_poll_loop = 0;
                        } else {
                            step_status = STEP_STEPNOTREADY;
                            // we may stay in poll loop
                        }
                    } 
                    else if (which_version == LAST_VERSION || 
                             which_version == NEXT_AVAILABLE_VERSION) 
                    {
                        step_status = STEP_OTHERERROR;
                        stay_in_poll_loop = 0;
                        log_warn ("DATASPACES method: Unexpected state: found last version %d"
                                "of dataset but then could not read it.\n", step);
                    }
                }
            }

        } else {
            // This stream does not exist yet
            log_info ("Data of '%s' does not exist (yet) in DataSpaces\n", fp->path);
            step_status = STEP_STREAMNOTFOUND;
        }

        if (step_status != STEP_OK)
            unlock_file (fp, ds);
        

        // check if we need to stay in loop 
        if (stay_in_poll_loop) {
            if (timeout_sec >= 0.0 && (adios_gettime()-t1 > timeout_sec))
                stay_in_poll_loop = 0;
            else
                adios_nanosleep (poll_interval_msec/1000, 
                     (int)(((uint64_t)poll_interval_msec * 1000000L)%1000000000L)); 
        }

    } // while (stay_in_poll_loop)

    // generate the appropriate error if needed
    switch (step_status) {
    case STEP_STREAMNOTFOUND:
            adios_error (err_file_not_found, 
                    "Data of '%s' does not exist in DataSpaces\n", fp->path);
            break;
    case STEP_STREAMTERMINATED:
            adios_errno = err_end_of_stream; // not an error that has to be printed
            log_debug ("Stream '%s' has been terminated. No more steps available\n", fp->path);
            break;
    case STEP_STEPNOTREADY:
            adios_errno = err_step_notready; 
            log_debug ("Step %d in stream '%s' is not yet available\n", step, fp->path);
            break;
    case STEP_STEPDISAPPEARED:
            adios_errno = err_step_disappeared; 
            log_debug ("Step %d in stream '%s' is not available anymore\n", step, fp->path);
            break;
    default:
            adios_errno = err_no_error; // clear temporary error during polling
            break;

    }

    return (step_status != STEP_OK); // 0 on success
}


ADIOS_FILE * adios_read_dataspaces_open_file (const char * fname, MPI_Comm comm)
{
    adios_error (err_operation_not_supported, 
                 "DATASPACES staging method does not support file mode for reading. "
                 "Use adios_read_open() to open a dataset.\n");
    return NULL;
}


ADIOS_FILE * adios_read_dataspaces_open (const char * fname, 
                                         MPI_Comm comm, 
                                         enum ADIOS_LOCKMODE lock_mode, 
                                         float timeout_sec)
{
    ADIOS_FILE * fp;
    struct dataspaces_data_struct * ds;
    int i;    

    ds = (struct dataspaces_data_struct *) malloc (sizeof(struct dataspaces_data_struct));
    if (!ds) {
        adios_error (err_no_memory, "Cannot allocate memory for file info.\n");
        return NULL;
    }

    fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    if (!fp) {
        adios_error (err_no_memory, "Cannot allocate memory for file info.\n");
        free(ds);
        return NULL;
    }

    fp->fh = (uint64_t) ds;
    fp->path = strdup(fname);
    ds->comm = comm;
    MPI_Comm_rank(comm, &ds->mpi_rank);
    MPI_Comm_size(comm, &ds->nproc);
    ds->current_step = -1;
    ds->lock_mode = lock_mode;
    ds->locked_fname = 0;
    ds->req_list = NULL;
    ds->req_list_tail = NULL;
    ds->nreq = 0;
    ds->chunk = NULL;

    /* if not connected to DATASPACES, connect now (and disconnect in adios_read_close) */
    if (!globals_adios_is_dataspaces_connected_from_reader()) {
        log_debug("-- %s, rank %d: call init first\n", __func__, ds->mpi_rank);
        if (!adios_read_dataspaces_init_method(comm, NULL)) {
            free(ds);
            free(fp->path);
            free(fp);
            return NULL;
        }
        ds->disconnect_at_close = 1;
    } else {
        ds->disconnect_at_close = 0;
    }

    /* fill out dataspaces method specific struct */
#ifndef DATASPACES_NO_VERSIONING
   // check this filename's version number
    int fidx;
    for (fidx=0; fidx<n_filenames; fidx++) {
        if (!strcmp(fname, file_versions.filename[fidx]))
            break;
    }
    if (fidx == n_filenames) {
        if (n_filenames < MAXNFILE) {
            file_versions.filename[ n_filenames ] = strdup(fname);
            file_versions.version [ n_filenames ] = 0;
            n_filenames++;
        } else {
            adios_error (err_too_many_files, "Too many different filenames has been used for adios_fopen().\n\tDATASPACES method allows max %d files\n", MAXNFILE);
            if (ds->disconnect_at_close) 
                adios_read_dataspaces_finalize_method();
            free(ds);
            free(fp->path);
            free(fp);
            return NULL;
        }
    }
    ds->file_index = fidx;
#endif
    log_debug("open stream filename=%s fidx=%d\n", fname, fidx);

    int err = get_step (fp, 0, NEXT_AVAILABLE_VERSION, timeout_sec);
    if (err) {
        free(ds);
        free(fp);
        fp = NULL;
    } else {
        log_debug("opened version %d of filename=%s\n", ds->current_step, fname);
    }

    return fp;
}

static void free_step_data (ADIOS_FILE *fp) 
{
    struct dataspaces_data_struct * ds = 
                (struct dataspaces_data_struct *) fp->fh;
    struct dataspaces_var_struct * vars = 
                (struct dataspaces_var_struct *) ds->vars;
    struct dataspaces_attr_struct * attrs = 
                (struct dataspaces_attr_struct *) ds->attrs;

    int i;

    if (fp->nvars) {
        for (i=0; i<fp->nvars; i++) {
            free (vars[i].name);
            free (vars[i].value);
        }
        MYFREE(ds->vars);
        a2s_free_namelist ((fp->var_namelist),fp->nvars);
        fp->nvars = 0;
    }
    if (fp->nattrs) {
        for (i=0; i<fp->nattrs; i++) {
            free (attrs[i].name);
            free (attrs[i].value);
        }
        MYFREE (ds->attrs);
        a2s_free_namelist ((fp->attr_namelist),fp->nattrs);
        fp->nattrs = 0;
    }
    MYFREE (ds->group_name);
}


int adios_read_dataspaces_close (ADIOS_FILE *fp) 
{
    struct dataspaces_data_struct * ds = 
                (struct dataspaces_data_struct *) fp->fh;

    log_debug("-- %s, rank %d: fp=%x\n", __func__, ds->mpi_rank, fp);

    /* Release read lock locked in fopen */
    unlock_file (fp, ds);

    /* Disconnect from DATASPACES if we connected at open() */
    if (ds && ds->disconnect_at_close) {
        adios_read_dataspaces_finalize_method();
    }

    free_step_data (fp);

    if (ds->chunk) MYFREE(ds->chunk);
    free (ds);
    if (fp->path) MYFREE(fp->path);
    free (fp);
    return 0;
}

int adios_read_dataspaces_peek_ahead (ADIOS_FILE *fp)
{
    struct dataspaces_data_struct * ds = 
                (struct dataspaces_data_struct *) fp->fh;

    log_debug("peek ahead from version %d of filename=%s, last available step %d\n", 
              ds->current_step, fp->path, fp->last_step);

    uint64_t offset[MAX_DS_NDIM], readsize[MAX_DS_NDIM];
    int ndim;
    int version_info_buf[2]; // 0: last version, 1: terminated?
    int err;
    char ds_vname[MAX_DS_NAMELEN];
    enum STEP_STATUS step_status = STEP_OK;
    snprintf(ds_vname, MAX_DS_NAMELEN, "VERSION@%s",fp->path);

    log_debug("   rank %d: dspaces_get %s\n", ds->mpi_rank, ds_vname);
    ndim = 1;
    offset[0] = 0; readsize[0] = 2; //*sizeof(int); // VERSION%name is 2 integers only
    err = get_meta(ds_vname, adios_integer, 0, ds->mpi_rank, ndim, 0, 
            offset, readsize, version_info_buf, ds->comm);

    if (!err) {
        int last_version = version_info_buf[0];
        int terminated = version_info_buf[1];
        log_debug("   rank %d: version info: last=%d, terminated=%d\n", 
                ds->mpi_rank, last_version, terminated);

        fp->last_step = last_version; // update to new last_step
        if (last_version <= fp->current_step) {
            // we have no more new steps
            if (terminated) {
                // stream is gone, we read everything 
                adios_errno = err_end_of_stream;
                log_debug ("Stream '%s' has been terminated. No more steps available\n", fp->path);
            } else {
                // a next step may come 
                adios_errno = err_step_notready;
                log_debug ("No new step in stream '%s' is not yet available\n", fp->path);
            }
        }
    }

    return adios_errno;
}

int adios_read_dataspaces_advance_step (ADIOS_FILE *fp, int last, float timeout_sec)
{
    struct dataspaces_data_struct * ds = 
                (struct dataspaces_data_struct *) fp->fh;

    enum WHICH_VERSION which_version;

#ifdef DATASPACES_NO_VERSIONING
    ds->current_step = -1;    /* Data in DataSpaces is always overwritten (read same version) */
    which_version = NEXT_VERSION;
#else
    if (last)
        which_version = LAST_VERSION; 
    else if (ds->lock_mode == ADIOS_LOCKMODE_ALL)
        which_version = NEXT_VERSION; 
    else
        which_version = NEXT_AVAILABLE_VERSION; 
    log_debug("advance from version %d of filename=%s which_version=%s\n", 
              ds->current_step, fp->path,
              which_version_str[which_version]);
#endif

    if (ds->locked_fname) {
        /* Release previous step (app did not call release_step() */
        adios_read_dataspaces_release_step (fp);
    }

    int err = get_step (fp, ds->current_step+1, which_version, timeout_sec); // content of fp and ds changes!

    if (!err) {
        file_versions.version [ ds->file_index ] = ds->current_step; // why do we store this?
        log_debug("advanced to version %d of filename=%s\n", ds->current_step, fp->path);
    } else {
        if (!adios_errno) {
            adios_error (err_unspecified, "Unspecified error during adios_advance_step()\n");
        }
    }

    return adios_errno;
}


void adios_read_dataspaces_release_step (ADIOS_FILE *fp)
{
    struct dataspaces_data_struct * ds = 
                (struct dataspaces_data_struct *) fp->fh;
    struct dataspaces_var_struct * vars = 
                (struct dataspaces_var_struct *) ds->vars;
    struct dataspaces_attr_struct * attrs = 
                (struct dataspaces_attr_struct *) ds->attrs;

    /* Release read lock locked in fopen */
    unlock_file (fp, ds);

    free_step_data (fp);
}


ADIOS_VARINFO * adios_read_dataspaces_inq_var_byid (const ADIOS_FILE *fp, int varid)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    struct dataspaces_var_struct * vars = ds->vars;
    ADIOS_VARINFO * vi;
    int i;
    int datasize;

    if (varid < 0 || varid > fp->nvars) {
        adios_error (err_invalid_varid, "Stream %s has %d variables. Invalid variable id %d\n",
                    fp->path, fp->nvars, varid);
        return NULL;
    }

    vi = (ADIOS_VARINFO *) malloc(sizeof(ADIOS_VARINFO));
    if (!vi) {
        adios_error (err_no_memory, "Could not allocate memory for variable info.\n");
        return NULL;
    }

    vi->varid = varid;
    vi->type = vars[varid].type;
    vi->nsteps = 1;

    /* Copy the dimensions (adios_free_varinfo() will free the copy */
    vi->ndim = vars[varid].ndims;
    if (vi->ndim) {
        vi->dims = (uint64_t *) malloc (vi->ndim*sizeof(uint64_t));
        memcpy (vi->dims, vars[varid].dims, vi->ndim*sizeof(uint64_t));
    } else {
        vi->dims = NULL;
    }

    /* Copy the value */
    if (vars[varid].value) {
        datasize = common_read_type_size(vi->type, vars[varid].value);
        vi->value = (void *) malloc (datasize);
        memcpy (vi->value, vars[varid].value, datasize);
    } else {
        vi->value = NULL;
    }

    vi->global = 1;
    vi->nblocks = (int *) malloc (sizeof(int));
    vi->nblocks[0] = 1;
    vi->sum_nblocks = vi->nblocks[0];
    vi->statistics = NULL;
    vi->blockinfo = NULL;
    
    return vi;
}

int adios_read_dataspaces_inq_var_stat (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo, int per_step_stat, int per_block_stat)
{
    /* FIXME: store and retrieve statistics from DataSpaces */
    varinfo->statistics = NULL;
    /*
    varinfo->statistics->min = NULL;
    varinfo->statistics->max = NULL;
    varinfo->statistics->avg = NULL;
    varinfo->statistics->std_dev = NULL;
    varinfo->statistics->steps = NULL;
    varinfo->statistics->blocks = NULL;
    varinfo->statistics->histogram = NULL;
    */
    return 0;
}

ADIOS_TRANSINFO* adios_read_dataspaces_inq_var_transinfo(const ADIOS_FILE *gp, const ADIOS_VARINFO *vi)
{
    ADIOS_TRANSINFO *trans = malloc(sizeof(ADIOS_TRANSINFO));
    memset(trans, 0, sizeof(ADIOS_TRANSINFO));
    trans->transform_type = adios_transform_none;
    return trans;
}



int adios_read_dataspaces_inq_var_blockinfo (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    /* FIXME: return the actual block decomposition by the writers
              but we need to store that in DataSpaces in the first place 
    */
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    int i;
    varinfo->blockinfo = (ADIOS_VARBLOCK *) malloc (sizeof(ADIOS_VARBLOCK)); // just one block
    varinfo->blockinfo->start = (uint64_t *) malloc (ds->vars[varinfo->varid].ndims * sizeof(uint64_t));
    varinfo->blockinfo->count = (uint64_t *) malloc (ds->vars[varinfo->varid].ndims * sizeof(uint64_t));
    for (i = 0; i<ds->vars[varinfo->varid].ndims; i++) {
        varinfo->blockinfo->start[i] = 0;
        varinfo->blockinfo->count[i] = ds->vars[varinfo->varid].dims[i];
    }
    varinfo->blockinfo->process_id = 0;
    varinfo->blockinfo->time_index = fp->current_step;
    return 0;
}

int adios_read_dataspaces_inq_var_trans_blockinfo(const ADIOS_FILE *gp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti)
{
    adios_error(err_operation_not_supported, "DataSpaces does not yet support transforms: trans_blockinfo.\n");
    return 0;
}


static int adios_read_dataspaces_get (const char * varname, enum ADIOS_DATATYPES vartype, 
                                int version, int rank,
                                int ndims, int is_fortran_ordering, 
                                uint64_t * offset, uint64_t * readsize, uint64_t * dims, void * data)
{

    struct obj_data *od;
    int elemsize = common_read_type_size(vartype, NULL);
    int i, err;
    int didx[MAX_DS_NDIM];
    uint64_t gdims[MAX_DS_NDIM], lb[MAX_DS_NDIM], ub[MAX_DS_NDIM];
    char gdims_str[256], lb_str[256], ub_str[256];

    // reorder DS dimensions to Fortran/C dimensions
    ds_dimension_ordering (ndims, is_fortran_ordering, 0 /*pack*/, didx);
    for (i=0; i<ndims; i++) {
        lb[i] = offset[didx[i]];
        ub[i] = offset[didx[i]]+readsize[didx[i]]-1;
        gdims[i] = dims[didx[i]];
    }

    ds_int64s_to_str(ndims, gdims, gdims_str);
    ds_int64s_to_str(ndims, lb, lb_str);
    ds_int64s_to_str(ndims, ub, ub_str);
    log_debug("-- %s, rank %d: get data: varname=%s version=%d, lb=(%s) ub=(%s) gdims=(%s)\n",
        __func__, rank, varname, version, lb_str, ub_str, gdims_str);

    dspaces_define_gdim (varname, ndims, gdims);
    err =  dspaces_get (varname, version, elemsize,
                ndims, lb, ub, data); 

    /*if (err == -ENOMEM) {
        adios_error (err_no_memory, "Not enough memory for DATASPACES to perform dspaces_get()");  
        return err_no_memory;
    } 
    else*/ if (err) {
        adios_error (err_corrupted_variable, "DATASPACES failed to read variable %s.\n", varname);  
        return err_corrupted_variable;
    }

    return 0;
}

static int adios_read_dataspaces_get_meta(const char * varname, enum ADIOS_DATATYPES vartype,
                                int version, int rank,
                                int ndims, int is_fortran_ordering,
                                uint64_t * offset, uint64_t * readsize, void * data)
{

    int elemsize = common_read_type_size(vartype, NULL);
    int i, err;
    int didx[MAX_DS_NDIM];
    uint64_t gdims[MAX_DS_NDIM], lb[MAX_DS_NDIM], ub[MAX_DS_NDIM];
    char gdims_str[256], lb_str[256], ub_str[256];

    // reorder DS dimensions to Fortran/C dimensions
    ds_dimension_ordering (ndims, is_fortran_ordering, 0 /*pack*/, didx);
    for (i=0; i<ndims; i++) {
        lb[i] = offset[didx[i]];
        ub[i] = offset[didx[i]]+readsize[didx[i]]-1;
        gdims[i] = (ub[i]-lb[i]+1) * dspaces_get_num_space_server();
    }

    ds_int64s_to_str(ndims, lb, lb_str);
    ds_int64s_to_str(ndims, ub, ub_str); 
    log_debug("-- %s, rank %d: get data: varname=%s version=%d, lb=(%s) ub=(%s)\n",
        __func__, rank, varname, version, lb_str, ub_str);

    dspaces_define_gdim(varname, ndims, gdims);
    err =  dspaces_get (varname, version, elemsize, ndims, lb, ub, data);
    /*if (err == -ENOMEM) {
        adios_error (err_no_memory, "Not enough memory for DATASPACES to perform dspaces_get()");  
        return err_no_memory;
    } 
    else*/ if (err) {
        log_debug ("DATASPACES failed to read metadata variable %s.\n", varname);
    }

    return err;
}

static int adios_read_dataspaces_get_meta_collective(const char * varname,
                                enum ADIOS_DATATYPES vartype,
                                int version, int rank,
                                int ndims, int is_fortran_ordering,
                                uint64_t * offset, uint64_t * readsize, void * data, MPI_Comm comm)
{
    int elemsize = common_read_type_size(vartype, NULL);
    int i, err;
    int didx[MAX_DS_NDIM];
    uint64_t gdims[MAX_DS_NDIM], lb[MAX_DS_NDIM], ub[MAX_DS_NDIM];
    char gdims_str[256], lb_str[256], ub_str[256];
    size_t datasize;
    size_t padsize = common_read_type_size(adios_integer, NULL);
    void *buf;
    int *pad;
    int root = 0;

    // calculate size of the meta data 
    for (i = 0, datasize = 1; i < ndims; i++) {
        datasize *= (readsize[i]*elemsize);
    }
    buf = malloc(padsize + datasize);

    // reorder DS dimensions to Fortran/C dimensions
    ds_dimension_ordering (ndims, is_fortran_ordering, 0 /*pack*/, didx);
    for (i=0; i<ndims; i++) {
        lb[i] = offset[didx[i]];
        ub[i] = offset[didx[i]]+readsize[didx[i]]-1;
        gdims[i] = (ub[i]-lb[i]+1) * dspaces_get_num_space_server();
    }

    ds_int64s_to_str(ndims, lb, lb_str);
    ds_int64s_to_str(ndims, ub, ub_str);
    log_debug("-- %s, rank %d: get data: varname=%s version=%d, lb=(%s) ub=(%s)\n",
        __func__, rank, varname, version, lb_str, ub_str);

    if (rank == root) {
        dspaces_define_gdim(varname, ndims, gdims);
        err =  dspaces_get (varname, version, elemsize, ndims, lb, ub, data);
        if (err) log_debug ("DATASPACES failed to read metadata variable %s.\n", varname);

        // set pad to indicate if root rank successfully fetch the meta data
        pad = (int*)buf;
        pad[0] = err;
        // copy meta data
        memcpy(buf+padsize, data, datasize);
    }

    // broadcast the meta data        
    if (MPI_Bcast(buf, datasize+padsize, MPI_BYTE, root, comm) == MPI_SUCCESS) {
        if (rank != root) {
            pad = (int*)buf;
            err = pad[0];
            memcpy(data, buf+padsize, datasize);
        }
    } else {
        err = -1;
    }

    free(buf);
    return err;
}

int adios_read_dataspaces_schedule_read_byid (const ADIOS_FILE * fp, 
                                              const ADIOS_SELECTION * sel, 
                                              int varid, 
                                              int from_steps, 
                                              int nsteps, 
                                              void * data) 
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    struct dataspaces_var_struct * var = &ds->vars[varid];
    read_request * r;
    uint64_t *s, *c, ld0, off0;
    uint64_t reqsize;
    int i;

    if (nsteps != 1) {
        adios_error (err_invalid_timestep, 
                     "Only one step can be read from a stream at a time. "
                     "You requested %d steps in adios_schedule_read()\n", nsteps);
        return err_invalid_timestep;
    }

    r = (read_request *) malloc (sizeof (read_request));
    if (!r) {
        adios_error (err_no_memory, "Could not allocate memory when scheduling a read request.\n");
        return err_no_memory;
    }
    r->sel = NULL;

    reqsize = common_read_type_size (var->type, NULL);
    // process, check selection and create target selection for read
    if (sel) {
        switch (sel->type) {

            case ADIOS_SELECTION_BOUNDINGBOX:

                if (var->ndims != sel->u.bb.ndim) {
                    adios_error (err_out_of_bound, 
                        "Number of dimensions in ADIOS_SELECTION = %d should be equal "
                        "to the number of dimensions of the variable %s = %d\n",
                        var->ndims, fp->var_namelist[varid], sel->u.bb.ndim);
                    return err_out_of_bound;
                }
                for (i=0; i<sel->u.bb.ndim; i++) {
                    if (sel->u.bb.start[i] + sel->u.bb.count[i] > var->dims[i]) {
                        adios_error (err_out_of_bound, 
                                "offset/readsize is out of bound in dimension %d for variable %s\n"
                                "size of dimension = %lld; you provided start=%lld, count=%lld\n",
                                i, fp->var_namelist[varid], var->dims[i], 
                                sel->u.bb.start[i], sel->u.bb.count[i]
                                );
                    }
                    reqsize *= sel->u.bb.count[i];
                }
                r->sel = copy_selection (sel);
                break;

            case ADIOS_SELECTION_POINTS:

                if (var->ndims != sel->u.bb.ndim) {
                    adios_error (err_out_of_bound, 
                        "Number of dimensions in ADIOS_SELECTION = %d should be equal "
                        "to the number of dimensions of the variable %s = %d\n",
                        var->ndims, fp->var_namelist[varid], sel->u.bb.ndim);
                    return err_out_of_bound;
                }
                reqsize *= sel->u.points.npoints;
                r->sel = copy_selection (sel);
                break;

            case ADIOS_SELECTION_WRITEBLOCK:

                /* We cannot do this with DataSpaces yet (fp->nwriter == 1) */
                /* Read the whole variable */
                s = (uint64_t *) calloc (var->ndims, sizeof(uint64_t));
                r->sel = common_read_selection_boundingbox(var->ndims, s, var->dims);
                for (i=0; i<var->ndims; i++) 
                    reqsize *= var->dims[i];
                break;

            case ADIOS_SELECTION_AUTO:

                /* We determine here what to read for this process.
                   Let's do a simple 1D domain decomposition
                   FIXME: should be smarter and do multi-dim decomp if needed
                */
                s = (uint64_t *) calloc (var->ndims, sizeof(uint64_t));
                c = (uint64_t *) malloc (var->ndims * sizeof(uint64_t));
                memcpy (c, var->dims, var->ndims*sizeof(uint64_t));
                if (var->ndims) {
                    ld0 = var->dims[0]/ds->nproc;
                    if (ld0 != 0) {
                        off0 = ds->mpi_rank*ld0;
                        if (ds->mpi_rank == ds->nproc-1) {
                            /* last reader reads the rest */
                            ld0 = var->dims[0] - (int)(ld0*(ds->nproc-1));
                        }
                    } else {
                        /* Only first dims[0] processes read one piece, the rest does nothing */
                        if (ds->mpi_rank < var->dims[0]) {
                            ld0 = 1;
                            off0 = ds->mpi_rank;
                        }
                    }
                    if (ld0 > 0) {
                        s[0] = off0;
                        c[0] = ld0;
                        r->sel = common_read_selection_boundingbox(
                                var->ndims, s, c);
                    }
                    for (i=0; i<var->ndims; i++) 
                        reqsize *= c[i];
                } else {
                    /* Scalar: just read it for each process */
                    r->sel = common_read_selection_boundingbox(0, 0, 0);
                }

                break;
        } // switch
    } else {
        // NULL selection means the whole variable
        s = (uint64_t *) calloc (var->ndims, sizeof(uint64_t));
        r->sel = common_read_selection_boundingbox(var->ndims, s, var->dims);
        for (i=0; i<var->ndims; i++) 
            reqsize *= var->dims[i];
    }

    if (r->sel) {
        r->varid = varid;
        r->from_steps = 0; // we read the current step anyway
        r->nsteps = 1;
        r->data = data;
        r->datasize = reqsize;
        r->priv = 0;
        r->next = 0;
        if (ds->req_list == NULL) {
            list_append_read_request_list (&ds->req_list, r);
            ds->req_list_tail = ds->req_list;
        } else {
            // just speed up insert directly after the tail
            list_append_read_request_list (&ds->req_list_tail, r);
        }
    } else {
        free(r);
    }

    return 0;
}

static ADIOS_VARCHUNK * read_var (const ADIOS_FILE *fp, read_request * r)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    struct dataspaces_var_struct * var = &ds->vars[r->varid];
    //int64_t total_size;
    uint64_t offset[MAX_DS_NDIM], readsize[MAX_DS_NDIM], dims[MAX_DS_NDIM];
    char offset_str[256], readsize_str[256];
    int elemsize;
    int err;
    int i,k,tidx;
    char ds_name[MAX_DS_NAMELEN];

    // set global dimension
    for (i = 0; i < var->ndims; i++) {
        dims[i] = var->dims[i];
    }

    if (!ds->chunk)
        ds->chunk = (ADIOS_VARCHUNK *) malloc (sizeof (ADIOS_VARCHUNK));
    if (!ds->chunk) {
        adios_error (err_no_memory, 
                "Could not allocate buffer for group name in adios_read_open()\n");
        return NULL;
    }
    ds->chunk->varid = r->varid;
    ds->chunk->type = var->type;
    ds->chunk->sel = NULL;
    elemsize = common_read_type_size (var->type, var->value);

    // handle scalars first (no need to read from space again)
    if (var->ndims == 0 && !var->hastime) { 
        if (r->data) {
            memcpy (r->data, var->value, elemsize);
            ds->chunk->data = r->data;
        } else {
            ds->chunk->data = var->value;
        }
        log_debug("-- %s, rank %d: scalar: varname=%s value=%d\n",
                __func__, ds->mpi_rank, var->name, *(int*)var->value); 
        return ds->chunk;
    }
        
    if (r->data) {
        // read into user allocated memory
        ds->chunk->data = r->data;
    } else {
        // read into method-allocated memory
        ds->chunk->data = get_chunk_buffer();
        if (!ds->chunk->data) {
            return NULL;
        }
    }

    snprintf(ds_name, MAX_DS_NAMELEN, "%s/%s/%s", fp->path, ds->group_name, var->name);

    if (r->sel->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        for (i=0; i<var->ndims; i++) {
            offset[i]    = r->sel->u.bb.start[i];
            readsize[i]  = r->sel->u.bb.count[i];
        }

        ds_int64s_to_str(var->ndims, offset, offset_str);
        ds_int64s_to_str(var->ndims, readsize, readsize_str);
        log_debug("-- %s, rank %d: get data: varname=%s offset=(%s) readsize=(%s)\n",
                __func__, ds->mpi_rank, ds_name, offset_str, readsize_str);

        err = adios_read_dataspaces_get (ds_name, var->type, ds->current_step, ds->mpi_rank, 
                var->ndims, futils_is_called_from_fortran(),
                offset, readsize, dims, ds->chunk->data);
    }
    else if (r->sel->type == ADIOS_SELECTION_POINTS)
    {
        err = 0;
        k = 0;
        while (!err && k < r->sel->u.points.npoints) {
            // pick and read k-th point individually
            for (i=0; i<var->ndims; i++) {
                offset[i]    = r->sel->u.points.points[k*var->ndims+i];
                readsize[i]  = 1;
            }

            err = adios_read_dataspaces_get (ds_name, var->type, ds->current_step, ds->mpi_rank, 
                    var->ndims, futils_is_called_from_fortran(),
                    offset, readsize, dims, ds->chunk->data+k*elemsize);
            k++;
        }
    }
    else
    {
        log_error ("Programming error: in Dataspaces method's read_var(), there "
                   "should be no selection type other than bounding box and points. "
                   "type = %d, var = %s\n", r->sel->type, var->name);
        err = 1;

    }

    if (err) {
        return NULL;
    }

    return ds->chunk;
}

int adios_read_dataspaces_perform_reads (const ADIOS_FILE *fp, int blocking)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    read_request * r;

    /* 1. prepare all reads */
    // check if all user memory is provided for blocking read
    if (blocking) {
        r = ds->req_list;
        while (r != NULL) {
            if (r->data == NULL) {
                adios_error (err_operation_not_supported, 
                    "Blocking mode at adios_perform_reads() requires that user "
                    "provides the memory for each read request. Request for "
                    "variable %s was scheduled without user-allocated memory\n",
                    ds->vars[r->varid].name);
                return err_operation_not_supported;
            }
            r = r->next;
        }
    }

    /* 2. if blocking, do all reads here, otherwise do it one-by-one in check_reads */
    if (!blocking) 
        return 0;

    while (ds->req_list != NULL && adios_errno == err_no_error) {
        read_var (fp, ds->req_list);
    
        // remove head from list
        r = ds->req_list;
        ds->req_list = ds->req_list->next;
        // FIXME: if we allocated start/count arrays in schedule read for r->sel,
        // we need to manually free them here
        common_read_selection_delete(r->sel);
        free(r);
        ds->nreq--;
    }
    ds->req_list_tail = NULL;

    return adios_errno;
}

int adios_read_dataspaces_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    read_request * r;
    int retval;

    if (ds->req_list != NULL) 
    {
        if (!ds->req_list->data && ds->req_list->datasize > chunk_buffer_size) {
            /* Request size does not fit into the chunk buffer.
               Chop up the request into multiple smaller chunks here.
            */
            /* FIXME: do this chunking */
            log_error ("DATASPACES method cannot do chunking at this moment. "
                        "Choose a max_chunk_size=N parameter so that each variable "
                        "request fits into the buffer\n");

        }

        *chunk = read_var (fp, ds->req_list);

        if (*chunk)
            retval = 1;
        else
            retval = adios_errno;
    
        // remove head from list
        r = ds->req_list;
        ds->req_list = ds->req_list->next;
        free(r);
        ds->nreq--;
    } 
    else 
    {
        // no more chunks (variables) to be read
        retval = 0; 
        ds->req_list_tail = NULL;
    }
    return retval;
}


/* Tell the DataSpaces order of dimensions for a 1-MAX_DS_NDIM dim array written from Fortran or C.
   unpack=1: the reverse of packing (to retrieve the original order).
   didx should be an int [MAX_DS_NDIM] array in any case.
*/
void ds_dimension_ordering(int ndims, int is_app_fortran, int unpack, int *didx)
{
    /* Order of dimensions: in DataSpaces: fast --> slow --> slowest
       For example: 
       Fortran: i,j,k --> i, j, k  = lb[0], lb[1], lb[2]
                i,j   --> i, j     = lb[0], lb[1]
                i     --> i        = lb[0]
       C:       i,j,k --> k, j, i  = lb[2], lb[1], lb[0]
                i,j   --> j, i     = lb[1], lb[0]
                i     --> i        = lb[0] 
    */
    int i;
    // initialize didx[]
    for (i = 0; i < MAX_DS_NDIM; i++) {
        didx[i] = i;
    }

    if (ndims == 0) return; 
    if (is_app_fortran) {
        for (i = 0; i < ndims; i++) {
            didx[i] = i;
        }
    } else {
        for (i = 0; i < ndims; i++) {
            didx[i] = ndims-1-i;
        }
    }

    return;
}

void ds_ints_to_str (int ndim, int *values, char *s)
{
    int i;
    char v[32];
    if (!ndim) {
        s[0] = '\0';
        return;
    }
    sprintf(s,"%d", values[0]);
    for (i=1; i<ndim; i++)
    {
        sprintf (v,",%d", values[i]);
        strcat (s,v);
    }
}

void ds_int64s_to_str(int ndim, uint64_t *values, char *s)
{
    int i;
    char v[32];
    if (!ndim) {
        s[0] = '\0';
        return;
    }
    sprintf(s,"%llu", values[0]);
    for (i=1; i<ndim; i++)
    {
        sprintf (v,",%llu", values[i]);
        strcat (s,v);
    }
}

int adios_read_dataspaces_get_attr_byid (const ADIOS_FILE * fp, int attrid, 
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    struct dataspaces_data_struct * ds = (struct dataspaces_data_struct *) fp->fh;
    struct dataspaces_attr_struct * attrs = ds->attrs;

    if (attrid < 0 || attrid > fp->nattrs) {
        adios_error (err_invalid_attrid, 
                     "File %s has %d attributes. Invalid attribute id %d\n",
                     fp->path, fp->nattrs, attrid);
        return adios_errno;
    }

    *type = attrs[attrid].type;
    *size = common_read_type_size(*type, attrs[attrid].value);
    *data = (void *) malloc (*size);
    if (*data) {
        memcpy (*data, attrs[attrid].value, *size);
    } else {
        adios_error (err_no_memory, "Could not allocate memory for attribute info.\n");
        return adios_errno;
    }
    return 0; 
}


int adios_read_dataspaces_get_dimension_order (const ADIOS_FILE *fp)
{
    return 0;
}


void adios_read_dataspaces_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{
    /* not implemented */
}

void adios_read_dataspaces_get_groupinfo (const ADIOS_FILE *fp, int *ngroups, 
            char ***group_namelist, uint32_t **nvars_per_group, uint32_t **nattrs_per_group) 
{
    struct dataspaces_data_struct * ds;
    if (fp) {
        ds = (struct dataspaces_data_struct *) fp->fh;
        *ngroups = 1;
        *group_namelist = (char **) malloc (sizeof (char*));
        *group_namelist[0] = strdup (ds->group_name);
    }
}

int adios_read_dataspaces_is_var_timed (const ADIOS_FILE *fp, int varid)
{
    struct dataspaces_data_struct * ds;
    int retval = 0;
    if (fp) {
        ds = (struct dataspaces_data_struct *) fp->fh;
        if (varid > 0 && varid < fp->nvars)
            retval = ds->vars[varid].hastime;
    }
    return retval;
}
