#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// xml parser
#include <mxml.h>

// for dataspaces ???
#include <sys/uio.h>

// see if we have MPI or other tools
#include "config.h"

#include "public/adios.h"
#include "public/adios_types.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_internals.h"
#include "core/adios_internals_mxml.h"
#include "core/util.h"
#include "core/ds_metadata.h"
#include "core/adios_logger.h"

#include "dimes_interface.h"
#include "dataspaces.h"

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#  define START_TIMER(t) adios_timing_go (fd->group->timing_obj, (t) ) 
#else
#  define START_TIMER(t) ; 
#endif

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#  define STOP_TIMER(t) adios_timing_stop (fd->group->timing_obj, (t) )
#else
#  define STOP_TIMER(t) ;
#endif

// Indices for the timer object
#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
static int T_DIMES_PUT   = ADIOS_TIMING_MAX_USER_TIMERS + 0;
static int T_GETLOCK     = ADIOS_TIMING_MAX_USER_TIMERS + 1;
static int T_MPI_BARRIER = ADIOS_TIMING_MAX_USER_TIMERS + 2;
static int T_MD          = ADIOS_TIMING_MAX_USER_TIMERS + 3;
static int T_AD_OPEN     = ADIOS_TIMING_MAX_USER_TIMERS + 4;
static int T_AD_WRITE    = ADIOS_TIMING_MAX_USER_TIMERS + 5;
static int T_AD_CLOSE    = ADIOS_TIMING_MAX_USER_TIMERS + 6;

static int timer_count = 7;
static char * timer_names[] = {
        "dimes_put",
        "dimes_lock",
        "mpi_barrier",
        "metadata",
        "ad_open",
        "ad_write",
        "ad_close",
};
#endif

/*#define DATASPACES_NO_VERSIONING  define it at configure as -DDATASPACES_NO_VERSIONING in CFLAGS */

static int adios_dimes_initialized = 0;
#define MAX_DS_NAMELEN 128
static char ds_var_name[MAX_DS_NAMELEN];
static unsigned int adios_dimes_verbose = 3;
static int check_read_status = 2; // 0: disable, 1: at every step (not supported yet), 2: at finalize (default value)
static double check_read_status_timeout_sec = 30;
static int check_read_status_poll_interval_ms = 100;
// count the number of inits/finalizes (one per adios group using this method
static unsigned int number_of_inits = 0;


struct adios_dimes_data_struct
{
    int rank;   // dataspaces rank or MPI rank if MPI is available
    int peers;  // from xml parameter or group communicator
    int appid;  // from xml parameter or 1
    int n_writes; // how many times adios_write has been called
#if HAVE_MPI
    MPI_Comm mpi_comm; // for use in open..close
    MPI_Comm mpi_comm_init; // for use in init/finalize
#endif
};


/**********************************************************************************
* Functions to manage the set of "files" or streams opened for all ADIOS groups
* We store all names (with version info, and responsible "rank 0" master process id
* to be used in adios_dataspaces_finalize().
**********************************************************************************/
#define MAX_NUM_OF_STREAMS 20
struct adios_dimes_stream_info
{
    char *name;         // file name passed in adios_open()
    int  time_index;    // versioning, start from 0
    int  iam_rank0;     // 1: current process has been rank 0 for this stream
                        //    rank 0 of communicator for this file does extra work in finalize
};
static struct adios_dimes_stream_info stream_info[MAX_NUM_OF_STREAMS];
static int  num_of_streams = 0; // how many files do we have with this method (in total for entire run)
                                // i.e. this variable never decreases
//char *fnames[MAX_NUM_OF_STREAMS];  // names of files (needed at finalize)
//int  fversions[MAX_NUM_OF_STREAMS];   // last steps of files (needed at finalize)
//int  mpi_ranks[MAX_NUM_OF_STREAMS];   // mpi rank of current process for each written file (needed at finalize)

static void free_dimes_stream_info()
{
    int i;
    struct adios_dimes_stream_info *info;
    for (i = 0; i < num_of_streams; i++) {
        info = &stream_info[i];
        if (info->name) {
            free(info->name);
        }
        info->name = NULL;
        info->time_index = -1; // time_index (dataspaces versioning) starts from 0
        info->iam_rank0 = 0;
    }
    return;
}

static struct adios_dimes_stream_info* lookup_dimes_stream_info(const char* fname)
{
    int i;
    // search from last to first
    for (i = num_of_streams-1; i >= 0; i--)
    {
        if (stream_info[i].name != NULL &&
            strcmp(stream_info[i].name, fname) == 0)
        {
            log_debug ("Stream %s is going to be continued... num_of_streams=%d\n",
                fname, num_of_streams);
            return &stream_info[i];
        }
    }
    // not found, add new opened stream to list
    if (num_of_streams < MAX_NUM_OF_STREAMS)
    {
        log_debug ("New stream %s added.  num_of_streams=%d\n",
                fname, num_of_streams);
        i = num_of_streams;
        num_of_streams++;
        stream_info[i].name = strdup(fname);
        stream_info[i].time_index = -1;
        return &stream_info[i];
    }
    else
    {
        // we cannot add more
        adios_error (err_too_many_files,
                     "ERROR: Max %d different files can be written by one application "
                     "using the same ADIOS group when using the DATASPACES method.\n",
                     MAX_NUM_OF_STREAMS);
    }

    return NULL;
}



static int check_read_status_var(const char* fname, int last_version)
{
    int stay_in_poll_loop = 1;
    double t1 = adios_gettime();

    uint64_t lb[MAX_DS_NDIM], ub[MAX_DS_NDIM], gdims[MAX_DS_NDIM];
    int elemsize, ndim;
    int read_status_buf[1] = {-1};
    int read_status_buf_len = 1;

    while (stay_in_poll_loop) {
        snprintf(ds_var_name, MAX_DS_NAMELEN, "READ_STATUS@%s", fname);
        elemsize = sizeof(int); ndim = 1;
        lb[0] = 0; ub[0] = read_status_buf_len-1;
        gdims[0] = (ub[0]-lb[0]+1) * dspaces_get_num_space_server();
        dspaces_define_gdim(ds_var_name, ndim, gdims);
        int err = dspaces_get(ds_var_name, 0, elemsize, ndim, lb, ub, read_status_buf);
        if (!err) {
            int version = read_status_buf[0];
            log_debug("%s: ds_var_name %s read_status_buf = {%d}\n",
                __func__, ds_var_name, version);
            if (version == last_version) {
                stay_in_poll_loop = 0;
            }            
        } else {
            log_error("%s: failed to read ds_var_name %s from space\n",
                __func__, ds_var_name); 
        }

        // check if we need to stay in loop
        if (stay_in_poll_loop) {
            double elapsed_time = adios_gettime() - t1;
            if (check_read_status_timeout_sec >= 0.0 &&
                elapsed_time > check_read_status_timeout_sec) {
                stay_in_poll_loop = 0;
            } else {
                adios_nanosleep(check_read_status_poll_interval_ms/1000,
                    (int)(((uint64_t)check_read_status_poll_interval_ms * 1000000L)%1000000000L));      
            }
        }
    }

    return 0;
}


static int connect_to_dimes (struct adios_dimes_data_struct *md, MPI_Comm comm)
{
    int ret = 0;
    int num_peers;

    if (!globals_adios_is_dimes_connected()) {

        MPI_Comm_rank (comm, &(md->rank));
        MPI_Comm_size (comm, &num_peers);

        // Application ID should be set by the application calling adios_set_application_id()
        int was_set;
        md->appid = globals_adios_get_application_id (&was_set);
        if (!was_set)
            md->appid = 1;

        log_debug ("adios_dimes: rank=%d connect to DATASPACES, peers=%d, appid=%d \n",
                md->rank, num_peers, md->appid);

        //Init the dart client
        ret = dspaces_init (num_peers, md->appid, &md->mpi_comm_init, NULL);
        if (ret) {
            log_error ("adios_dimes: rank=%d Failed to connect to DATASPACES: err=%d,  rank=%d\n", md->rank, ret);        
            return ret;
        }

#if ! HAVE_MPI
        dspaces_rank (&(md->rank));
        dspaces_peers (&(md->peers));
#else
        md->peers = num_peers;
#endif

        log_debug ("adios_dimes: rank=%d connected to DATASPACES: peers=%d\n", md->rank, md->peers);        
    }

    globals_adios_set_dimes_connected_from_writer();
    return ret;
}


void adios_dimes_init (const PairStruct * parameters,
                     struct adios_method_struct * method
                     )
{
    struct adios_dimes_data_struct *md = 0;
    if (!adios_dimes_initialized)
    {
        adios_dimes_initialized = 1;
    }
   
    method->method_data = calloc (1, sizeof (struct adios_dimes_data_struct));
    md = (struct adios_dimes_data_struct*)method->method_data;
   
    int check_read; 
    int index, i;
    char temp[64];

    //Init the static data structure
    md->peers = 1;
    md->appid = -1;
    md->n_writes = 0;
#if HAVE_MPI
    md->mpi_comm = MPI_COMM_NULL;
    md->mpi_comm_init = method->init_comm;
#endif

    // process user parameters
    const PairStruct *p = parameters;
    while (p) {
        if (!strcasecmp (p->name, "app_id")) {
            errno = 0;
            md->appid = strtol(p->value, NULL, 10);
            if (md->appid > 0 && !errno) {
                log_debug ("App ID parameter set to %d for DIMES write method\n",
                            md->appid);
                globals_adios_set_application_id (md->appid);
            } else {
                log_error ("Invalid 'app_id' parameter given to the DIMES write "
                           "method: '%s'\n", p->value);
            }
        } else if (!strcasecmp(p->name, "check_read_status")) {
            errno = 0;
            check_read = strtol(p->value, NULL, 10);
            if (!errno && (check_read == 0 || check_read == 2)) {
                check_read_status = check_read;
                log_debug("check_read_status set to %d for DIMES write method\n", 
                    check_read_status);
            } else {
                log_error("Invalid 'check_read_status' parameter given to the DIMES "
                            "write method: '%s'\n", p->value);
                log_error("check_read_status=<value>, 0: disable, 1: at every step "
                            " (not supported yet), 2: at finalize (default value).\n");
            }
        } else if (!strcasecmp(p->name, "check_read_status_timeout_sec")) {
            errno = 0;
            double timeout = strtof(p->value, NULL);
            if (timeout > 0.0 && !errno) {
                log_debug("check_read_status_timeout_sec set to %f seconds for DIMES write method\n", timeout);
                check_read_status_timeout_sec = timeout;
            } else {
                log_error("Invalid 'check_read_status_timeout_sec' parameter given to the DIMES "
                        "write method: '%s'\n", p->value);
            }   
        } else if (!strcasecmp(p->name, "check_read_status_poll_interval")) {
            errno = 0;
            int pollinterval = strtol(p->value, NULL, 10);
            if (pollinterval > 0 && !errno) {
                log_debug("check_read_status_poll_interval set to %d milliseconds for DIMES write method\n", pollinterval);
                check_read_status_poll_interval_ms = pollinterval;
            } else {
                log_error("Invalid 'check_read_status_poll_interval' parameter given to DIMES "
                    "write method: '%s'\n", p->value);
            }
        } else {
            log_error ("Parameter name %s is not recognized by the DIMES write "
                        "method\n", p->name);
        }
        p = p->next;
    }
    connect_to_dimes (md, method->init_comm);
    number_of_inits++;

    log_info ("adios_dimes_init: called the %d. time\n", number_of_inits);
   
}



int adios_dimes_open (struct adios_file_struct * fd,
                    struct adios_method_struct * method,
                    MPI_Comm comm
                    )
{
    int ret = 0;
    struct adios_dimes_data_struct *md = (struct adios_dimes_data_struct *)
                                                method->method_data;
    if (fd->mode == adios_mode_read)
    {
        adios_error (err_operation_not_supported,
                "DIMES transport method does not support old adios_read() calls. "
                "Use the ADIOS read API and it's DATASPACES method.\n");
        return adios_errno;
    }

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
    // Ensure both timing objects exist
    // timing_obj should get created at every open for the same file
    // prev_timing_obj should only be created at the first open 
    if (fd->group)
    {
        if (!fd->group->timing_obj)
            fd->group->timing_obj = adios_timing_create (timer_count, timer_names);

        if (!fd->group->prev_timing_obj)
            fd->group->prev_timing_obj = adios_timing_create (timer_count, timer_names);
    }
#endif
    START_TIMER (T_AD_OPEN);

    struct adios_dimes_stream_info *info = lookup_dimes_stream_info(fd->name);
    if (!info) {
        return adios_errno;
    }

    /* Increment the time index, start with 0. Not good to increment in close, because 
     * finalize needs the same time_index as the last writing cycle */
    info->time_index++;
    log_info ("adios_dimes_open: open %s, mode=%d, time_index=%d \n",
                        fd->name, fd->mode, info->time_index);


#if HAVE_MPI
    // if we have MPI and a communicator, we can get the exact size of this application
    // that we need to tell DATASPACES
    md->mpi_comm = comm;
    MPI_Comm_rank (md->mpi_comm, &(md->rank));
    MPI_Comm_size (md->mpi_comm, &(md->peers));
#endif

    info->iam_rank0 = (md->rank == 0);

    log_debug ("adios_dimes_open: rank=%d call write lock...\n", md->rank);       
    START_TIMER (T_GETLOCK);
    dspaces_lock_on_write (fd->name, &md->mpi_comm);  
    STOP_TIMER (T_GETLOCK);
    log_debug ("adios_dimes_open: rank=%d got write lock\n", md->rank);        
    // Free data objects written in the previous steps
    dimes_put_sync_group(fd->name, info->time_index);
    dimes_put_set_group(fd->name, info->time_index);    

    STOP_TIMER (T_AD_OPEN);
    return ret;
}

enum BUFFERING_STRATEGY adios_dimes_should_buffer (struct adios_file_struct * fd
                                                  ,struct adios_method_struct * method
                                                  )
{
    return no_buffering;  
}


void adios_dimes_write (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,const void * data
                      ,struct adios_method_struct * method
                      )
{
    if (fd->mode == adios_mode_read) {
        return;
    }
    START_TIMER (T_AD_WRITE);

    struct adios_dimes_data_struct *md = (struct adios_dimes_data_struct *)
                                                            method->method_data;
    struct adios_group_struct *group = fd->group;
    struct adios_dimes_stream_info *info = lookup_dimes_stream_info(fd->name);
    //Get var size
    //  FIXME: type size of a string >2GB does not fit to int. 
    //  adios_get_type_size returns uint64_t but dspaces_put handles only int
    //  as element size
    int var_type_size = (int) adios_get_type_size(v->type, v->data);
    //Get var name
    char * var_name = v->name;
    int err;

    char lb_str[256], ub_str[256], gdims_str[256], dims_str[256], didx_str[256];
    //Get two offset coordinate values
    unsigned int version;
    uint64_t dims[MAX_DS_NDIM], gdims[MAX_DS_NDIM], lb[MAX_DS_NDIM], ub[MAX_DS_NDIM]; /* lower and upper bounds for DataSpaces */
    int didx[MAX_DS_NDIM]; // for reordering the dimensions
    int ndims = 0;
    int hastime = 0;
    uint64_t nelems = 1;
    gdims[0] = 0;
    struct adios_dimension_struct* var_dimensions = v->dimensions;
    // Calculate lower and upper bounds for each available dimension (up to 3 dims)
    while( var_dimensions && ndims < MAX_DS_NDIM)
    {
        dims[ndims] = adios_get_dim_value (&(var_dimensions->dimension));
        gdims[ndims] = adios_get_dim_value (&(var_dimensions->global_dimension));
        lb[ndims] = adios_get_dim_value (&(var_dimensions->local_offset));
        if (gdims[ndims] > 0 && dims[ndims] > 0)  {
            ub[ndims] = lb[ndims] + dims[ndims] - 1;
            nelems *= dims[ndims];
            ndims++;
        } else if (gdims[ndims] > 0 && dims[ndims] == 0) {
            // piece of array with 0 elements, skip it
            nelems *= dims[ndims];
        }   else {
            // time dimension (ldim=0 indicates this). Leave out from the dimensions.
            //ub[ndims] = lb[ndims]; 
            hastime = 1;
        }
        var_dimensions = var_dimensions->next;
    }

#ifdef DATASPACES_NO_VERSIONING
    version = 0;              /* Update/overwrite data in DataSpaces  (we write time_index as a variable at close)*/
#else
    version = info->time_index;  /* Add new data as separate to DataSpaces */
#endif
    
    if (v->path != NULL && v->path[0] != '\0' && strcmp(v->path,"/")) 
        snprintf(ds_var_name, MAX_DS_NAMELEN, "%s/%s/%s/%s", fd->name, fd->group->name, v->path, v->name);
    else if (!strcmp(v->path,"/"))
        snprintf(ds_var_name, MAX_DS_NAMELEN, "%s/%s//%s", fd->name, fd->group->name, v->name);
    else
        snprintf(ds_var_name, MAX_DS_NAMELEN, "%s/%s/%s", fd->name, fd->group->name, v->name);

    //snprintf(dspaces_type_var_name, MAX_DS_NAMELEN, "TYPE@%s", ds_var_name);
    
    /* The next line is just to fix adios_build_index_v1() call not to abort on 
     * offset=0 variables (in case of files, written variables have offset > 0)
     */
    v->write_offset = 1; 

    /* non-global variables are put in space ONLY by rank = 0 process */
    if (gdims[0] == 0 && md->rank != 0) {
        //fprintf(stderr, "rank=%d var_name=%s is not global. Skip\n", md->rank, ds_var_name);
        return;
    }

    
    //if (fd->shared_buffer == adios_flag_no)
    //{
        // var payload sent for sizing information
        //adios_write_var_header_v1 (fd, v);
    //}
    
     
    /* This is not needed here, this is already called in common_adios_write() 
    adios_generate_var_characteristics_v1 (fd, v); // characteristics will be included in build index
    adios_write_var_characteristics_v1 (fd, v);
    */

    dimes_int64s_to_str(ndims, lb, lb_str);
    dimes_int64s_to_str(ndims, ub, ub_str);
    dimes_int64s_to_str(ndims, dims, dims_str);
    dimes_int64s_to_str(ndims, gdims, gdims_str);
    log_debug ("var_name=%s, type=%s(%d) elemsize=%d, version=%d, ndims=%d, size=(%s), gdim=(%s), lb=(%s), ub=(%s)\n",
            ds_var_name, adios_type_to_string_int(v->type), v->type, var_type_size, version, ndims,
            dims_str, gdims_str, lb_str, ub_str);    


    /* If variable is empty, do not push to space */
    if (nelems == 0) {
        log_debug ("rank=%d var_name=%s is empy variable piece. Skip\n", md->rank, ds_var_name);
        return;
    }

    /* non-timed scalars are written in the metadata at close(), not here */
    if (ndims == 0 && !hastime)
        return;

    /* Put type info as T<varname>, integer in 0,0,0,0,0,0 position */
    //err = dspaces_put(dspaces_type_var_name, version, 4, 0,0,0,0,0,0, &(v->type)); 

    dimes_dimension_ordering(ndims,
            group->adios_host_language_fortran == adios_flag_yes, 
            0 /*pack*/, didx);

    uint64_t lb_in[MAX_DS_NDIM], ub_in[MAX_DS_NDIM], gdims_in[MAX_DS_NDIM];
    int i;
    for (i = 0; i < ndims; i++) {
        lb_in[i] = lb[didx[i]];
        ub_in[i] = ub[didx[i]];
        gdims_in[i] = gdims[didx[i]];
    }
    dimes_define_gdim(ds_var_name, ndims, gdims_in);
    START_TIMER (T_DIMES_PUT);
    dimes_put(ds_var_name, version, var_type_size, ndims, lb_in, ub_in, data);
    STOP_TIMER (T_DIMES_PUT);

    dimes_ints_to_str(ndims, didx, didx_str);
    dimes_int64s_to_str(ndims, gdims_in, gdims_str);
    dimes_int64s_to_str(ndims, lb_in, lb_str);
    dimes_int64s_to_str(ndims, ub_in, ub_str);
    log_debug ("var_name=%s, dimension ordering=(%s), gdims=(%s), lb=(%s), ub=(%s)\n",
            ds_var_name, didx_str, gdims_str, lb_str, ub_str);
    STOP_TIMER (T_AD_WRITE);
}

void adios_dimes_get_write_buffer (struct adios_file_struct * fd
                                 ,struct adios_var_struct * v
                                 ,uint64_t * size
                                 ,void ** buffer
                                 ,struct adios_method_struct * method
                                 )
{
    uint64_t mem_allowed;

    if (*size == 0)
    {
        *buffer = 0;

        return;
    }

    if (v->adata && v->free_data == adios_flag_yes)
    {
        adios_method_buffer_free (v->data_size);
        free (v->adata);
        v->data = v->adata = NULL;
    }

    mem_allowed = adios_method_buffer_alloc (*size);
    if (mem_allowed == *size)
    {
        *buffer = malloc (*size);
        if (!*buffer)
        {
            adios_method_buffer_free (mem_allowed);
            log_error ("ERROR: Out of memory allocating %llu bytes for %s in %s:%s()\n"
                    ,*size, v->name, __FILE__, __func__
                    );
            v->got_buffer = adios_flag_no;
            v->free_data = adios_flag_no;
            v->data_size = 0;
            v->data = 0;
            *size = 0;
            *buffer = 0;
        }
        else
        {
            v->got_buffer = adios_flag_yes;
            v->free_data = adios_flag_yes;
            v->data_size = mem_allowed;
            v->data = *buffer;
        }
    }
    else
    {
        adios_method_buffer_free (mem_allowed);
        log_error ("OVERFLOW: Cannot allocate requested buffer of %llu "
                         "bytes for %s in %s:%s()\n"
                ,*size
                ,v->name
                ,__FILE__, __func__
                );
        *size = 0;
        *buffer = 0;
    }
}

/* NOT IMPLEMENTED. Use the Read API to read variables */
void adios_dimes_read (struct adios_file_struct * fd
                     ,struct adios_var_struct * v, void * buffer
                     ,uint64_t buffer_size
                     ,struct adios_method_struct * method
                     )
{
}

/* Gather var/attr indices from all processes to rank 0 */
static void adios_dimes_gather_indices (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               ,struct adios_index_struct_v1 * index
                               )
{
    struct adios_dimes_data_struct *md = (struct adios_dimes_data_struct *)
                                                method->method_data;
    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;
    
    // build local index first appending to any existing index
    adios_build_index_v1 (fd, index);

    log_debug ("%s index after first build is pg=%x vars=%x attrs=%x\n", 
                __func__, index->pg_root, index->vars_root, index->attrs_root);
#if 0
#if HAVE_MPI
    // gather all on rank 0
    if (md->mpi_comm != MPI_COMM_NULL)
    {                                
        if (md->rank == 0)           
        {                            
            int * index_sizes = malloc (4 * md->peers);
            int * index_offsets = malloc (4 * md->peers);
            char * recv_buffer = 0;
            uint32_t size = 0;
            uint32_t total_size = 0;
            int i;
            struct adios_bp_buffer_struct_v1 b;

            MPI_Gather (&size, 1, MPI_INT
                    ,index_sizes, 1, MPI_INT
                    ,0, md->mpi_comm
                    );

            for (i = 0; i < md->peers; i++)
            {
                index_offsets [i] = total_size;
                total_size += index_sizes [i];
            }                    

            recv_buffer = malloc (total_size);

            MPI_Gatherv (&size, 0, MPI_BYTE
                    ,recv_buffer, index_sizes, index_offsets
                    ,MPI_BYTE, 0, md->mpi_comm
                    );

            for (i = 1; i < md->peers; i++)
            {
                b.buff = recv_buffer + index_offsets [i];
                b.length = index_sizes [i];
                b.offset = 0;

                adios_parse_process_group_index_v1 (&b
                        ,&new_pg_root
                        );
                adios_parse_vars_index_v1 (&b, &new_vars_root, NULL, NULL);
                adios_parse_attributes_index_v1 (&b
                        ,&new_attrs_root
                        );
                adios_merge_index_v1 (index, new_pg_root, 
                                      new_vars_root, new_attrs_root, 0);
                new_pg_root = 0;
                new_vars_root = 0;
                new_attrs_root = 0;
            }

            free (recv_buffer);
            free (index_sizes);
            free (index_offsets);
        }
        else
        {
            char * buffer = 0; 
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;

            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                    ,0, index);

            uint32_t tmp_buffer_size = (uint32_t) buffer_size;
            MPI_Gather (&tmp_buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                    ,0, md->mpi_comm
                    );
            MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                    ,0, 0, 0, MPI_BYTE
                    ,0, md->mpi_comm
                    );
            free (buffer);
        }
    }

#endif
#endif

    log_debug ("%s index after gathering is pg=%x vars=%x attrs=%x\n", 
                __func__, index->pg_root, index->vars_root, index->attrs_root);
}

static int dimes_get_full_name_len (char * path, char * name)
{
    int len;
    // make full name
    if (!path || !path[0]) { 
        // no path, just name + leading /
        len = strlen(name);
    } else if (!strcmp (path, "/")) {
        len = strlen(name)+1;
    } else {
        len = strlen(path) + strlen(name) + 1;
    }
    return len;
}

static int dimes_get_full_name (char * path, char * name, int maxlen,
                            /*OUT*/char * out)
{
    int len;
    // make full name
    if (!path || !path[0]) { 
        // no path, just name + leading /
        len = strlen(name);
        strncpy(out, name, maxlen);
    } else if (!strcmp (path, "/")) {
        len = strlen(name) + 1;
        out[0] = '/';
        strncpy(out+1, name, maxlen-1);
    } else {
        len = strlen(path);
        strncpy(out, path, maxlen-1);  // path +
        out[len] = '/';                         //   /  +
        strncpy(out+len+1, name, maxlen-len-1); // name
        len += strlen(name) + 1;
    }
    return len;
}

void dimes_pack_group_info (struct adios_file_struct *fd
                                  ,struct adios_method_struct * method
                                  ,struct adios_index_struct_v1 *index
                                  ,char ** buffer, int *buffer_size, int *nvars, int *nattrs
                                  )
{
    struct adios_dimes_data_struct *md = (struct adios_dimes_data_struct *)
                                                method->method_data;
    struct adios_index_var_struct_v1 * v = index->vars_root;
    struct adios_index_attribute_struct_v1 * a = index->attrs_root;
    int size;
    int ndims; // whatever the type of v->characteristics->dims.count is, we write an int to buffer
    int hastime; // true if variable has time dimension
    uint64_t ldims[MAX_DS_NDIM], gdims[MAX_DS_NDIM];
    *nvars = 0;
    *nattrs = 0;
    int didx[MAX_DS_NDIM]; // dimension ordering indices

    log_debug ("%s entered\n", __func__);

    /* First cycle: count the size of info to allocate index buffer */
    size = 3*sizeof(int); //header for buffer: length, nvars, nattrs
    while (v) {
        size += 4*sizeof(int) // name len, type, hastime, number of dims 
                + dimes_get_full_name_len (v->var_path, v->var_name) // full path
                + MAX_DS_NDIM * 8; // always write 3 dimensions in the index (even for scalars)
        if (v->characteristics->dims.count == 0) {
            // For scalars, we write the value into the index
            if (v->type != adios_string)
                size += adios_get_type_size(v->type, NULL);
            else
                size += adios_get_type_size(v->type, v->characteristics->value) + sizeof(int);
        }
        log_debug (" var %s/%s, size = %d\n", v->var_path, v->var_name, size);
        (*nvars)++;
        v = v->next;
    }

    while (a) {
        size += sizeof(int) 
                + dimes_get_full_name_len (a->attr_path, a->attr_name)
                + sizeof(int); // type
        if (a->type != adios_string)
            size += adios_get_type_size(a->type, NULL);
        else
            size += adios_get_type_size(a->type, a->characteristics->value) + sizeof(int);

        log_debug (" attr %s/%s, size = %d\n", a->attr_path, a->attr_name, size);
        (*nattrs)++;
        a = a->next;
    }

    // Required for Cray Gemini: align buffer to 8 bytes boundaries
    int align_bytes = 0; // number of extra bytes at the end
    if (size % 8) {
        align_bytes = 8 - (size % 8);
        size += align_bytes;
        log_debug (" after alignment, size = %d, align_bytes = %d\n", size, align_bytes);
    }

    *buffer = (char *) malloc (size);
    *buffer_size = size;

    /* Second cycle: fill up the buffer */
    v = index->vars_root;
    a = index->attrs_root;
    char * b = *buffer;
    int i, j, namelen;
    char name[256];

    //header for buffer: length, nvars, nattrs
    memcpy (b, buffer_size, sizeof(int));  
    b += sizeof(int); 
    memcpy (b, nvars, sizeof(int));  
    b += sizeof(int); 
    memcpy (b, nattrs, sizeof(int));  
    b += sizeof(int); 
    while (v) {
        namelen = dimes_get_full_name (v->var_path, v->var_name, sizeof(name), name);
        memcpy (b, &namelen, sizeof(int));  // length of full path
        b += sizeof(int); 
        memcpy (b, name, namelen);          // full path
        b += namelen;
        memcpy (b, &(v->type), sizeof(int)); // type 
        b += sizeof(int); 
        //ndims = MAX(v->characteristics->dims.count,3); // convert whatever type to int
        //memcpy (b, &(v->characteristics->dims.count), sizeof(int)); // number of dimensions
        log_debug("Variable %s, total dims = %d\n", name, v->characteristics->dims.count);
        j = 0; // will drop the time dimension
        hastime = 0;
        for (i = 0; i<v->characteristics->dims.count; i++) {
            ldims[j] = v->characteristics->dims.dims[j*3];  // ith dimension 
            gdims[j] = v->characteristics->dims.dims[j*3+1];  // ith dimension 
            log_debug("           , ldim = %lld gdim = %lld)\n", ldims[j], gdims[j]);
            if (gdims[j] == 0 && ldims[j] == 1) {
                // time dimension's global=0, local=1, skip
                // FIXME: This is true for a local array of length 1 (not defined as global)
                log_debug("               skip this dimension )\n");
                hastime = 1;
                continue;
            }
            j++;
        }
        for (i=j; i<MAX_DS_NDIM; i++) {
            // fill up dimensions up to MAX_DS_NDIM dim
            ldims[i] = 1;
            gdims[i] = 1;
        }
        ndims = (j < MAX_DS_NDIM ? j : MAX_DS_NDIM); // we can have max MAX_DS_NDIM dimensions in DataSpaces
        memcpy (b, &hastime, sizeof(int)); // has time dimension?
        log_debug("             has time = %d (%d)\n", hastime, *(int*)b);
        b += sizeof(int); 
        memcpy (b, &ndims, sizeof(int)); // number of dimensions
        log_debug("             ndims = %d (%d)\n", ndims, *(int*)b);
        b += sizeof(int); 
        dimes_dimension_ordering(ndims, 
                fd->group->adios_host_language_fortran == adios_flag_yes, 
                0 /*pack*/, didx);
        for (i = 0; i < MAX_DS_NDIM; i++) {
            if (gdims[didx[i]]) { 
                // global variable
                memcpy (b, &(gdims[didx[i]]), 8);  // ith dimension 
            } else { 
                // a local variable has no global dimensions
                // in space, its local dimensions become the global dimensions
                memcpy (b, &(ldims[didx[i]]), 8);  // ith dimension 
            }
            b += 8; 
        }
        if (v->characteristics->dims.count == 0) {
            // NOTE: ndims = 0 can mean a timed scalar, which has no characteristics->value!
            // store scalar value too
            if (v->type != adios_string) {
                size = adios_get_type_size(v->type, NULL);
                memcpy (b, v->characteristics->value, size); 
                b += size; 
            } else {
                size = adios_get_type_size(v->type, v->characteristics->value);
                memcpy (b, &size, sizeof(int)); 
                b += sizeof(int); 
                memcpy (b, v->characteristics->value, size); 
                b += size; 
            }
        }
        v = v->next;
    }

    while (a) {
        namelen = dimes_get_full_name (a->attr_path, a->attr_name, sizeof(name), name);
        memcpy (b, &namelen, sizeof(int));  // length of full path
        b += sizeof(int); 
        memcpy (b, name, namelen);          // full path
        b += namelen;
        memcpy (b, &(a->type), sizeof(int)); // type 
        b += sizeof(int); 
        // store scalar value too
        if (a->type != adios_string) {
            size = adios_get_type_size(a->type, NULL);
            memcpy (b, a->characteristics->value, size); 
            b += size; 
        } else {
            size = adios_get_type_size(a->type, a->characteristics->value);
            memcpy (b, &size, sizeof(int)); 
            b += sizeof(int); 
            memcpy (b, a->characteristics->value, size); 
            b += size; 
        }
        a = a->next;
    }

    // alignment
    if (align_bytes) {
        uint64_t zero = 0;
        memcpy (b, &zero, align_bytes); 
        b += align_bytes;
    }

    // sanity check
    if ( (int)(b-*buffer) > *buffer_size) {
        log_error ("ERROR in %s. Calculated group index buffer size as %d, but filled after that with %d bytes\n",
            __func__, *buffer_size, (int)(b-*buffer));
    }
    // written buffer might be shorter than calculated since we skip time dimensions.
    // set the correct size now
    *buffer_size = (int)(b-*buffer);
    memcpy (*buffer, buffer_size, sizeof(int));  

    
    log_debug("   %s: buffer length = %d, content:\n", __func__, *buffer_size);
    b = *buffer;
    for (i=0; i<*buffer_size; i+=16) {
        for (j=0; j<4; j++) {
            log_debug_cont ("%3.3hhu %3.3hhu %3.3hhu %3.3hhu    ", 
                            b[i+4*j], b[i+4*j+1], b[i+4*j+2], b[i+4*j+3]);
        }
        log_debug_cont("\n");
    }
    

    log_debug ("%s exit\n", __func__);
}

/* FIXME: put this function into ds_metadata.c */
/* buff is allocated and must be freed after use */
void dimes_pack_file_info (int time, int nvars, int nattrs, int group_index_len, char * groupname, 
                    /*OUT*/char **buf, /*OUT*/int *buf_len)
{
    *buf_len = 128;
    *buf = (char *) malloc (*buf_len);

    char *b = *buf;
    int namelen = strlen(groupname);
    memcpy (b, buf_len, sizeof(int));  /* 0-: length of this buffer */
    b += sizeof(int);
    memcpy (b, &time, sizeof(int));  /* 4-: time */
    b += sizeof(int); 
    memcpy (b, &nvars, sizeof(int));  /* 8-: number of variables */
    b += sizeof(int);
    memcpy (b, &nattrs, sizeof(int));  /* 12-: number of attributes */
    b += sizeof(int);
    memcpy (b, &group_index_len, sizeof(int));  /* 16-: length of group index*/
    b += sizeof(int);
    memcpy (b, &namelen, sizeof(int));  /* 20-: length of group name */
    b += sizeof(int);
    memcpy (b, groupname, namelen);  /* 24-: group name */
    b[namelen] = 0;
}

void adios_dimes_buffer_overflow (struct adios_file_struct * fd, 
                                  struct adios_method_struct * method)
{
    // this call never happens without shared buffering
}

void adios_dimes_close (struct adios_file_struct * fd
                      ,struct adios_method_struct * method
                      )
{
    if (fd->mode == adios_mode_read) {
        return;
    }
    START_TIMER (T_AD_CLOSE);

    struct adios_dimes_data_struct *md = (struct adios_dimes_data_struct *)
                                                method->method_data;
    struct adios_index_struct_v1 * index = adios_alloc_index_v1(1);
    struct adios_attribute_struct * a = fd->group->attributes;
    struct adios_dimes_stream_info *info = lookup_dimes_stream_info(fd->name);
    uint64_t gdims[MAX_DS_NDIM], lb[MAX_DS_NDIM], ub[MAX_DS_NDIM];
    int didx[MAX_DS_NDIM]; // for reordering DS dimensions
    int elemsize, ndim;
    unsigned int version;

    START_TIMER (T_MD);
    // finalize variable info in fd buffer, next we call build_index
    while (a) {
        a->write_offset = 1; // only attributes with !=0 offset will be included in build index
        a=a->next;
    }

    //adios_write_close_vars_v1 (fd);
    /* Gather var/attr indices from all processes to rank 0 */
    adios_dimes_gather_indices (fd, method, index);
    STOP_TIMER (T_MD);

    // make sure all processes have finished putting data to the space 
    // before we put metadata from rank 0
    START_TIMER (T_MPI_BARRIER);
    MPI_Barrier (md->mpi_comm); 
    STOP_TIMER (T_MPI_BARRIER);

    if (md->rank == 0) {

        START_TIMER (T_MD);
        /* Write two adios specific variables with the name of the file and name of the group into the space */
        /* ADIOS Read API fopen() checks these variables to see if writing already happened */
#ifdef DATASPACES_NO_VERSIONING
        version = 0;              /* Update/overwrite data in DataSpaces */
#else
        version = info->time_index;  /* Add new data as separate to DataSpaces */
#endif

        /* Make metadata from indices */
        char * indexbuf;
        int    indexlen;
        int    nvars, nattrs;
        dimes_pack_group_info (fd, method, index, 
                &indexbuf, &indexlen, &nvars, &nattrs);


        /* Put GROUP@fn/gn header into space */
        snprintf(ds_var_name, MAX_DS_NAMELEN, "GROUP@%s/%s", fd->name, fd->group->name);
        log_debug ("%s: put %s buflen=%d (bytes) into space\n", __func__, ds_var_name, indexlen);
        elemsize = 1; ndim = 1;
        lb[0] = 0; ub[0] = indexlen-1;
        gdims[0] = (ub[0]-lb[0]+1) * dspaces_get_num_space_server();
        dspaces_define_gdim(ds_var_name, ndim, gdims);
        dspaces_put(ds_var_name, version, elemsize, ndim, lb, ub, indexbuf);
        free (indexbuf);

        /* Create and put FILE@fn header into space */
        char * file_info_buf; /* store FILE@fn's group list */
        int    file_info_buf_len; /* = 128 currently */
        snprintf (ds_var_name, MAX_DS_NAMELEN, "FILE@%s", fd->name);
        dimes_pack_file_info (info->time_index, nvars, nattrs, indexlen,
                fd->group->name, &file_info_buf, &file_info_buf_len);
        log_debug ("%s: put %s buflen=%d (bytes) time=%d nvars=%d nattr=%d index=%d name=%d:%s into space\n",
                __func__, ds_var_name, 
                *(int*)file_info_buf, *(int*)(file_info_buf+4), 
                *(int*)(file_info_buf+8), *(int*)(file_info_buf+12),
                *(int*)(file_info_buf+16), *(int*)(file_info_buf+20),
                file_info_buf+24);
        dspaces_put_sync(); //wait on previous put to finish
        elemsize = 1; ndim = 1;
        lb[0] = 0; ub[0] = file_info_buf_len-1;
        gdims[0] = (ub[0]-lb[0]+1) * dspaces_get_num_space_server();
        dspaces_define_gdim(ds_var_name, ndim, gdims);
        dspaces_put(ds_var_name, version, elemsize, ndim, lb, ub, file_info_buf);

        /* Create and put VERSION@fn version info into space */
        int version_buf[2] = {version, 0}; /* last version put in space; not terminated */
        int version_buf_len = 2; 
        snprintf (ds_var_name, MAX_DS_NAMELEN, "VERSION@%s", fd->name);
        log_debug ("%s: put %s buf= [%d,%d] buflen=%d (integers) into space\n", 
                __func__, ds_var_name, version_buf[0], version_buf[1], version_buf_len);
        dspaces_put_sync(); //wait on previous put to finish
        elemsize = sizeof(int); ndim = 1;
        lb[0] = 0; ub[0] = version_buf_len-1;
        gdims[0] = (ub[0]-lb[0]+1) * dspaces_get_num_space_server();
        dspaces_define_gdim(ds_var_name, ndim, gdims);
        dspaces_put(ds_var_name, 0, elemsize, ndim, lb, ub, version_buf);
        dspaces_put_sync(); //wait on previous put to finish
        STOP_TIMER (T_MD);
    }


    // free allocated index lists
    adios_clear_index_v1 (index);
    adios_free_index_v1 (index);

    // rank=0 may be in put_sync when others call unlock, which is a global op
    START_TIMER (T_MPI_BARRIER);
    MPI_Barrier (md->mpi_comm); 
    STOP_TIMER (T_MPI_BARRIER);
    //log_debug("%s: call dspaces_put_sync()\n", __func__);
    //dspaces_put_sync();
    dimes_put_unset_group();
    log_debug("%s: call dspaces_unlock_on_write(%s)\n", __func__, fd->name);
    dspaces_unlock_on_write(fd->name, &md->mpi_comm);

    STOP_TIMER (T_AD_CLOSE);

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
    //Finished timing this cycle, swap the timing buffers
    adios_timing_destroy(fd->group->prev_timing_obj);
    fd->group->prev_timing_obj = fd->group->timing_obj;
    fd->group->timing_obj = 0;
    // prev_timing_obj points to unwritten timing info, timing_obj is
    // ready to allocate at the next open
#endif

    log_info ("%s: exit\n", __func__);
}

void adios_dimes_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_dimes_data_struct *md = (struct adios_dimes_data_struct *)
        method->method_data;
    struct adios_dimes_stream_info *info;
    int i;
    char ds_var_name[MAX_DS_NAMELEN];
    uint64_t gdims[MAX_DS_NDIM], lb[MAX_DS_NDIM], ub[MAX_DS_NDIM];
    int elemsize, ndim;
    int value[2] = {0, 1}; // integer to be written to space (terminated=1)

    log_debug("%s: called the %d. time, rank=%d\n", __func__, number_of_inits, mype);

    number_of_inits--;
    if (number_of_inits == 0)
    {
        // tell the readers which files are finalized
        for (i=0; i<num_of_streams; i++) {
            info = &stream_info[i];
            /* Put VERSION@fn into space. Indicates that this file will not be extended anymore.  */
            if (info->iam_rank0) {
                if (check_read_status == 2) {
                    check_read_status_var(info->name, info->time_index);
                }
                MPI_Comm mpi_comm = MPI_COMM_SELF;
                log_debug("%s: call dspaces_lock_on_write(%s), rank=%d\n", __func__, info->name, mype);
                dspaces_lock_on_write(info->name, &mpi_comm); // lock is global operation in DataSpaces

                value[0] = info->time_index;
                snprintf(ds_var_name, MAX_DS_NAMELEN, "VERSION@%s", info->name);
                log_debug ("%s: update %s in the space [%d, %d]\n", 
                        __func__, ds_var_name, value[0], value[1] );
                elemsize = sizeof(int); ndim = 1;
                lb[0] = 0; ub[0] = 1;
                gdims[0] = (ub[0]-lb[0]+1) * dspaces_get_num_space_server();
                dspaces_define_gdim(ds_var_name, ndim, gdims);
                dspaces_put(ds_var_name, 0, elemsize, ndim, lb, ub, &value);
                log_debug("%s: call dspaces_put_sync()\n", __func__);
                dspaces_put_sync();

                log_debug("%s: call dspaces_unlock_on_write(%s), rank=%d\n", __func__, info->name, mype);
                dspaces_unlock_on_write(info->name, &mpi_comm);
            }
        }

        free_dimes_stream_info();

        if (check_read_status == 2) {
            // Note: dspaces_lock_on_write() above is only called by single process (whose md->rank == 0). MPI_Barrier ensures all writer processes to wait until reader application fetches data of last version. 
            MPI_Barrier(md->mpi_comm_init);
        }
        // Free all previsouly allocated RDMA buffers
        dimes_put_sync_all();

        // disconnect from dataspaces if we are connected from writer but not anymore from reader
        if (globals_adios_is_dimes_connected_from_writer() && 
                !globals_adios_is_dimes_connected_from_both())
        {
            log_debug ("%s: call MPI Barrier on all connected processes(), rank=%d\n", __func__,mype);
            MPI_Barrier (md->mpi_comm_init); 
            log_debug ("%s: call dspaces_finalize(), rank=%d\n", __func__, mype);
            dspaces_finalize();
        }
        globals_adios_set_dimes_disconnected_from_writer();
        adios_dimes_initialized = 0;
    }

    log_debug("%s: done, remaining groups = %d, rank=%d\n", __func__, number_of_inits, mype);
}

void adios_dimes_end_iteration (struct adios_method_struct * method)
{
}

void adios_dimes_start_calculation (struct adios_method_struct * method)
{
}

void adios_dimes_stop_calculation (struct adios_method_struct * method)
{
}
