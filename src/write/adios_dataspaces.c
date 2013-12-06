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

#include "dataspaces.h"

/*#define DATASPACES_NO_VERSIONING  define it at configure as -DDATASPACES_NO_VERSIONING in CFLAGS */

static int adios_dataspaces_initialized = 0;
#define MAX_DS_NAMELEN 128
#define MAX_NUM_OF_FILES 20
//static char ds_type_var_name[MAX_DS_NAMELEN];
static char ds_var_name[MAX_DS_NAMELEN];
static unsigned int adios_dataspaces_verbose = 3;

struct adios_dspaces_file_info {
    char *name;
    int time_index;
};

struct adios_ds_data_struct
{
    int rank;   // dataspaces rank or MPI rank if MPI is available
    int peers;  // from xml parameter or group communicator
    int appid;  // from xml parameter or 1
    int time_index; // versioning in DataSpaces, start from 0
    int n_writes; // how many times adios_write has been called
    struct adios_dspaces_file_info file_info[MAX_NUM_OF_FILES];
#if HAVE_MPI
    MPI_Comm mpi_comm;
#endif
    int  num_of_files; // how many files do we have with this method
    char *fnames[MAX_NUM_OF_FILES];  // names of files (needed at finalize)
    int  fversions[MAX_NUM_OF_FILES];   // last steps of files (needed at finalize)
};

static int init_dspaces_file_info(struct adios_ds_data_struct *p)
{
    int i;
    for (i = 0; i < MAX_NUM_OF_FILES; i++) {
        p->file_info[i].name = NULL;
        p->file_info[i].time_index = 0;
    }
}

static void free_dspaces_file_info(struct adios_ds_data_struct *p)
{
    int i;
    for (i = 0; i < MAX_NUM_OF_FILES; i++) {
        if (p->file_info[i].name) {
            free(p->file_info[i].name);
        }
    }

    return;
}

static struct adios_dspaces_file_info* lookup_dspaces_file_info(struct adios_ds_data_struct *p, const char* fname)
{
    int i;    for (i = 0; i < MAX_NUM_OF_FILES; i++) {
        if (p->file_info[i].name != NULL &&
            strcmp(p->file_info[i].name, fname) == 0) {
            return &p->file_info[i];
        }
    }
    for (i = 0; i < MAX_NUM_OF_FILES; i++) {
        if (p->file_info[i].name == NULL) {
            p->file_info[i].name = malloc(strlen(fname)+1);
            strcpy(p->file_info[i].name, fname);
            return &p->file_info[i];
        }
    }

    return NULL;
}

static int connect_to_dspaces (struct adios_ds_data_struct * p, MPI_Comm comm)
{
    int ret = 0;
    int num_peers;

    if (!globals_adios_is_dataspaces_connected()) {

        MPI_Comm_rank (comm, &(p->rank));
        MPI_Comm_size (comm, &num_peers);

        // Application ID should be set by the application calling adios_set_application_id()
        int was_set;
        p->appid = globals_adios_get_application_id (&was_set);
        if (!was_set)
            p->appid = 1;

        log_debug ("adios_dataspaces: rank=%d connect to DATASPACES, peers=%d, appid=%d \n",
                p->rank, num_peers, p->appid);

        //Init the dart client
        ret = dspaces_init (num_peers, p->appid);
        if (ret) {
            log_error ("adios_dataspaces: rank=%d Failed to connect to DATASPACES: err=%d,  rank=%d\n", p->rank, ret);        
            return ret;
        }

#if ! HAVE_MPI
        dspaces_rank (&(p->rank));
        dspaces_peers (&(p->peers));
#endif

        log_debug ("adios_dataspaces: rank=%d connected to DATASPACES: peers=%d\n", p->rank, p->peers);        
    }

    globals_adios_set_dataspaces_connected_from_writer();
    return ret;
}


void adios_dataspaces_init (const PairStruct * parameters,
                     struct adios_method_struct * method
                     )
{
    struct adios_ds_data_struct *p = 0;
    if (!adios_dataspaces_initialized)
    {
        adios_dataspaces_initialized = 1;
    }
   
    method->method_data = calloc (1, sizeof (struct adios_ds_data_struct));
    p = (struct adios_ds_data_struct*)method->method_data;
    
    int index, i;
    char temp[64];

    //Init the static data structure
    p->peers = 1;
    p->appid = -1;
    p->time_index = 0;
    p->n_writes = 0;
#if HAVE_MPI
    p->mpi_comm = MPI_COMM_NULL;
#endif
    p->num_of_files = 0;

    init_dspaces_file_info(p);
    connect_to_dspaces (p, method->init_comm);

    log_info ("adios_dataspaces_init: done\n");
   
}



int adios_dataspaces_open (struct adios_file_struct * fd,
                    struct adios_method_struct * method,
                    MPI_Comm comm
                    )
{
    int ret = 0;
    struct adios_ds_data_struct *p = (struct adios_ds_data_struct *)
                                                method->method_data;
    struct adios_dspaces_file_info *info = lookup_dspaces_file_info(p,fd->name);
    log_info ("adios_dataspaces_open: open %s, mode=%d, time_index=%d \n",
                        fd->name, fd->mode, info->time_index);

#if HAVE_MPI
    // if we have MPI and a communicator, we can get the exact size of this application
    // that we need to tell DATASPACES
    p->mpi_comm = comm;
    MPI_Comm_rank (p->mpi_comm, &(p->rank));
    MPI_Comm_size (p->mpi_comm, &(p->peers));
#endif

    // connect to DATASPACES at the very first adios_open(), disconnect in adios_finalize()
    // connect only if the READ API has not connected yet
    /*
    ret = connect_to_dspaces (p, p->mpi_comm);
    if (ret)
        return ret;
    */

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        log_debug ("adios_dataspaces_open: rank=%d call write lock...\n", p->rank);        
        dspaces_lock_on_write (fd->name, &p->mpi_comm);  
        log_debug ("adios_dataspaces_open: rank=%d got write lock\n", p->rank);        
    }
    else if (fd->mode == adios_mode_read)
    {
        dspaces_lock_on_read (fd->name, &p->mpi_comm);
    } 
  
    return ret;
}

enum ADIOS_FLAG adios_dataspaces_should_buffer (struct adios_file_struct * fd
                                         ,struct adios_method_struct * method
                                         )
{
    
    //if (fd->shared_buffer == adios_flag_no && fd->mode != adios_mode_read)
    //{
        // write the process group header
        //adios_write_process_group_header_v1 (fd, fd->write_size_bytes);
        //adios_write_open_vars_v1 (fd);
    //} else {
    //    log_warn("WARNING: %s expects that fd->shared_buffer is false\n", __func__);
    //}
    

    return adios_flag_no;  // this will take care of it
}


void adios_dataspaces_write (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * data
                      ,struct adios_method_struct * method
                      )
{
    struct adios_ds_data_struct *p = (struct adios_ds_data_struct *)
                                                            method->method_data;
    struct adios_group_struct *group = fd->group;
    struct adios_dspaces_file_info *info = lookup_dspaces_file_info(p,fd->name);
    //Get var size
    //  FIXME: type size of a string >2GB does not fit to int. 
    //  adios_get_type_size returns uint64_t but dspaces_put handles only int
    //  as element size
    int var_type_size = (int) adios_get_type_size(v->type, v->data);
    //Get var name
    char * var_name = v->name;
    int err;

    //Get two offset coordinate values
    unsigned int version;

    int dims[3]={1,1,1}, gdims[3]={0,0,0}, lb[3]={0,0,0}, ub[3]={0,0,0}; /* lower and upper bounds for DataSpaces */
    int didx[3]; // for reordering the dimensions
    int ndims = 0;
    int hastime = 0;
    struct adios_dimension_struct* var_dimensions = v->dimensions;
    // Calculate lower and upper bounds for each available dimension (up to 3 dims)
    while( var_dimensions && ndims < 3)
    {
        dims[ndims] = adios_get_dim_value (&(var_dimensions->dimension));
        gdims[ndims] = adios_get_dim_value (&(var_dimensions->global_dimension));
        lb[ndims] = adios_get_dim_value (&(var_dimensions->local_offset));
        if (dims[ndims] > 0)  {
            ub[ndims] = lb[ndims] + dims[ndims] - 1;
            ndims++;
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
    
    /* non-global variables are put in space ONLY by rank = 0 process */
    if (gdims[0] == 0 && p->rank != 0) {
        //fprintf(stderr, "rank=%d var_name=%s is not global. Skip\n", p->rank, ds_var_name);
        return;
    }

    
    //if (fd->shared_buffer == adios_flag_no)
    //{
        // var payload sent for sizing information
        //adios_write_var_header_v1 (fd, v);
    //}
    
     
    v->write_offset = 1; // only !=0 offsets will be included in build index
    adios_generate_var_characteristics_v1 (fd, v); // characteristics will be included in build index
    adios_write_var_characteristics_v1 (fd, v);
    

    log_debug ("var_name=%s, type=%s(%d) elemsize=%d, version=%d, ndims=%d, size=(%d,%d,%d), gdim=(%d,%d,%d), lb=(%d,%d,%d), ub=(%d,%d,%d)\n",
            ds_var_name, adios_type_to_string_int(v->type), v->type, var_type_size, version, ndims,
            dims[0], dims[1], dims[2], gdims[0], gdims[1], gdims[2], lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);

    /* non-timed scalars are written in the metadata at close(), not here */
    if (ndims == 0 && !hastime)
        return;

    /* Put type info as T<varname>, integer in 0,0,0,0,0,0 position */
    //err = dspaces_put(dspaces_type_var_name, version, 4, 0,0,0,0,0,0, &(v->type)); 

    ds_dimension_ordering(ndims,
            group->adios_host_language_fortran == adios_flag_yes, 
            0 /*pack*/, didx);

    dspaces_put(ds_var_name, version, var_type_size, 
             lb[didx[0]], lb[didx[1]], lb[didx[2]], 
             ub[didx[0]], ub[didx[1]], ub[didx[2]], 
             data);
    
    log_debug ("var_name=%s, dimension ordering=(%d,%d,%d), gdims=(%d,%d,%d), lb=(%d,%d,%d), ub=(%d,%d,%d)\n",
            ds_var_name, 
            didx[0], didx[1], didx[2], 
            gdims[didx[0]], gdims[didx[1]], gdims[didx[2]], 
            lb[didx[0]], lb[didx[1]], lb[didx[2]], 
            ub[didx[0]], ub[didx[1]], ub[didx[2]]);
    dspaces_put_sync();
}

void adios_dataspaces_get_write_buffer (struct adios_file_struct * fd
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

    if (v->data && v->free_data == adios_flag_yes)
    {
        adios_method_buffer_free (v->data_size);
        free (v->data);
        v->data = NULL;
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
void adios_dataspaces_read (struct adios_file_struct * fd
                     ,struct adios_var_struct * v, void * buffer
                     ,uint64_t buffer_size
                     ,struct adios_method_struct * method
                     )
{
    struct adios_ds_data_struct *p = (struct adios_ds_data_struct *)
                                                            method->method_data;
    uint64_t var_type_size = adios_get_type_size(v->type, v->data);

    //Get var name
    char * var_name = v->name;

    //Get two offset coordinate values
    int version, offset1[3],offset2[3];
    int dim_size[3];
    memset(offset1, 0, 3*sizeof(int));
    memset(offset2, 0, 3*sizeof(int));
    memset(dim_size, 0, 3*sizeof(int));

    struct adios_dspaces_file_info *info = lookup_dspaces_file_info(p,fd->name);
    version = info->time_index;
    //dspaces_lock_on_read_();

    //dspaces_get

    //dspaces_unlock_on_read_();
}

/* Gather var/attr indices from all processes to rank 0 */
static void adios_dataspaces_gather_indices (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               ,struct adios_index_struct_v1 * index
                               )
{
    struct adios_ds_data_struct *p = (struct adios_ds_data_struct *)
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
    if (p->mpi_comm != MPI_COMM_NULL)
    {                                
        if (p->rank == 0)           
        {                            
            int * index_sizes = malloc (4 * p->peers);
            int * index_offsets = malloc (4 * p->peers);
            char * recv_buffer = 0;
            uint32_t size = 0;
            uint32_t total_size = 0;
            int i;
            struct adios_bp_buffer_struct_v1 b;

            MPI_Gather (&size, 1, MPI_INT
                    ,index_sizes, 1, MPI_INT
                    ,0, p->mpi_comm
                    );

            for (i = 0; i < p->peers; i++)
            {
                index_offsets [i] = total_size;
                total_size += index_sizes [i];
            }                    

            recv_buffer = malloc (total_size);

            MPI_Gatherv (&size, 0, MPI_BYTE
                    ,recv_buffer, index_sizes, index_offsets
                    ,MPI_BYTE, 0, p->mpi_comm
                    );

            for (i = 1; i < p->peers; i++)
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
                                      new_vars_root, new_attrs_root);
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
                    ,0, p->mpi_comm
                    );
            MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                    ,0, 0, 0, MPI_BYTE
                    ,0, p->mpi_comm
                    );
            free (buffer);
        }
    }

#endif
#endif

    log_debug ("%s index after gathering is pg=%x vars=%x attrs=%x\n", 
                __func__, index->pg_root, index->vars_root, index->attrs_root);
}

static int ds_get_full_name_len (char * path, char * name)
{
    int len;
    // make full name
    if (!path || !path[0]) { 
        // no path, just name
        len = strlen(name);
    } else if (!strcmp (path, "/")) {
        len = strlen(name)+1;
    } else {
        len = strlen(path) + strlen(name) + 1;
    }
    return len;
}

static int ds_get_full_name (char * path, char * name, int maxlen,
                            /*OUT*/char * out)
{
    int len;
    // make full name
    if (!path || !path[0] || !strcmp (path, "/")) { 
        // no path, just name 
        len = strlen(name);
        strncpy(out, name, maxlen);
    } else if (!strcmp (path, "/")) {
        len = strlen(name)+1;
        out[0]='/';
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

void ds_pack_group_info (struct adios_file_struct *fd
                        ,struct adios_method_struct * method
                        ,struct adios_index_struct_v1 *index
                        ,char ** buffer, int *buffer_size, int *nvars, int *nattrs
                        )
{
    struct adios_ds_data_struct *p = (struct adios_ds_data_struct *)
                                                method->method_data;
    struct adios_index_var_struct_v1 * v = index->vars_root;
    struct adios_index_attribute_struct_v1 * a = index->attrs_root;
    int size;
    int ndims; // whatever the type of v->characteristics->dims.count is, we write an int to buffer
    int hastime; // true if variable has time dimension
    uint64_t ldims[10], gdims[10]; // we can write only 3 dimensions, will drop time dim
    *nvars = 0;
    *nattrs = 0;
    int didx[3]; // dimension ordering indices

    log_debug ("%s entered\n", __func__);

    /* First cycle: count the size of info to allocate index buffer */
    size = 3*sizeof(int); //header for buffer: length, nvars, nattrs
    while (v) {
        size += 4*sizeof(int) // name len, type, hastime, number of dims 
                + ds_get_full_name_len (v->var_path, v->var_name) // full path
                + 3 * 8; // always write 3 dimensions in the index (even for scalars)
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
                + ds_get_full_name_len (a->attr_path, a->attr_name)
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
        namelen = ds_get_full_name (v->var_path, v->var_name, sizeof(name), name);
        memcpy (b, &namelen, sizeof(int));  // length of full path
        b += sizeof(int); 
        memcpy (b, name, namelen);          // full path
        b += namelen;
        memcpy (b, &(v->type), sizeof(int)); // type 
        b += sizeof(int); 
        //ndims = MAX(v->characteristics->dims.count,3); // convert whatever type to int
        //memcpy (b, &(v->characteristics->dims.count), sizeof(int)); // number of dimensions
        log_debug("Variable %s, total dims = %d\n", name, v->characteristics->dims.count);
        j = 0; // we can write only 3 dims, will drop the time dimension
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
        for (i=j; i<3; i++) {
            // fill up dimensions up to 3rd dim
            ldims[i] = 1;
            gdims[i] = 1;
        }
        ndims = (j < 3 ? j : 3); // we can have max 3 dimensions in DataSpaces
        memcpy (b, &hastime, sizeof(int)); // has time dimension?
        log_debug("             has time = %d (%d)\n", hastime, *(int*)b);
        b += sizeof(int); 
        memcpy (b, &ndims, sizeof(int)); // number of dimensions
        log_debug("             ndims = %d (%d)\n", ndims, *(int*)b);
        b += sizeof(int); 
        ds_dimension_ordering(ndims, 
                fd->group->adios_host_language_fortran == adios_flag_yes, 
                0 /*pack*/, didx);
        for (i = 0; i < 3; i++) {
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
        namelen = ds_get_full_name (a->attr_path, a->attr_name, sizeof(name), name);
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
void ds_pack_file_info (int time, int nvars, int nattrs, int group_index_len, char * groupname, 
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

void adios_dataspaces_close (struct adios_file_struct * fd
                      ,struct adios_method_struct * method
                      )
{
    struct adios_ds_data_struct *p = (struct adios_ds_data_struct *)
                                                method->method_data;
    struct adios_index_struct_v1 * index = adios_alloc_index_v1(1);
    struct adios_attribute_struct * a = fd->group->attributes;
    struct adios_dspaces_file_info *info = lookup_dspaces_file_info(p,fd->name);
    int lb[3], ub[3], didx[3]; // for reordering DS dimensions
    unsigned int version;

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        // finalize variable info in fd buffer, next we call build_index
        while (a) {
            a->write_offset = 1; // only attributes with !=0 offset will be included in build index
            a=a->next;
        }

        //adios_write_close_vars_v1 (fd);
        /* Gather var/attr indices from all processes to rank 0 */
        adios_dataspaces_gather_indices (fd, method, index);

        // make sure all processes have finished putting data to the space 
        // before we put metadata from rank 0
        MPI_Barrier (p->mpi_comm); 

        if (p->rank == 0) {

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
            ds_pack_group_info (fd, method, index, 
                                   &indexbuf, &indexlen, &nvars, &nattrs);

            
            /* Put GROUP@fn/gn header into space */
            snprintf(ds_var_name, MAX_DS_NAMELEN, "GROUP@%s/%s", fd->name, fd->group->name);
            log_debug ("%s: put %s with buf len %d into space\n", __func__, ds_var_name, indexlen);
            ub[0] = indexlen-1; ub[1] = 0; ub[2] = 0;
            ds_dimension_ordering(1, 0, 0, didx); // C ordering of 1D array into DS
            dspaces_put(ds_var_name, version, 1,    0, 0, 0, /* lb 0..2 */
                     ub[didx[0]], ub[didx[1]], ub[didx[2]],  indexbuf); 
            free (indexbuf);

            /* Create and put FILE@fn header into space */
            char * file_info_buf; /* store FILE@fn's group list */
            int    file_info_buf_len; /* = 128 currently */
            snprintf (ds_var_name, MAX_DS_NAMELEN, "FILE@%s", fd->name);
            ds_pack_file_info (info->time_index, nvars, nattrs, indexlen,
                          fd->group->name, &file_info_buf, &file_info_buf_len);
            log_debug ("%s: put %s = buflen=%d time=%d nvars=%d nattr=%d index=%d name=%d:%s into space\n",
                __func__, ds_var_name, 
                *(int*)file_info_buf, *(int*)(file_info_buf+4), 
                *(int*)(file_info_buf+8), *(int*)(file_info_buf+12),
                *(int*)(file_info_buf+16), *(int*)(file_info_buf+20),
                file_info_buf+24);
            /* Flip 1st and 2nd dimension for DataSpaces representation for a 1D array*/
            ub[0] = file_info_buf_len-1; ub[1] = 0; ub[2] = 0;
            ds_dimension_ordering(1, 0, 0, didx); // C ordering of 1D array into DS
            dspaces_put_sync(); //wait on previous put to finish
            dspaces_put(ds_var_name, version, 1,    0, 0, 0, /* lb 0..2 */
                     ub[didx[0]], ub[didx[1]], ub[didx[2]], file_info_buf); 

            /* Create and put VERSION@fn version info into space */
            int version_buf[2] = {version, 0}; /* last version put in space; not terminated */
            int version_buf_len = 2; 
            snprintf (ds_var_name, MAX_DS_NAMELEN, "VERSION@%s", fd->name);
            log_debug ("%s: put %s with buf = [%d,%d] (len=%d integers) into space\n", 
                       __func__, ds_var_name, version_buf[0], version_buf[1], version_buf_len);
            ub[0] = version_buf_len-1; ub[1] = 0; ub[2] = 0;
            ds_dimension_ordering(1, 0, 0, didx); // C ordering of 1D array into DS
            dspaces_put_sync(); //wait on previous put to finish
            dspaces_put(ds_var_name, 0, sizeof(int),    0, 0, 0, /* lb 0..2 */
                     ub[didx[0]], ub[didx[1]], ub[didx[2]],  version_buf); 
            dspaces_put_sync(); //wait on previous put to finish
            
        }

        // remember this filename and its version for finalize
        int i;
        for (i=0; i<p->num_of_files; i++) {
            if (!strcmp(fd->name, p->fnames[i]))
                break;
        }
        if (i == p->num_of_files) {
            if (p->num_of_files < MAX_NUM_OF_FILES) {
                p->fnames[ p->num_of_files ] = strdup(fd->name);
                p->num_of_files++;
            } else {
                log_error ("%s: Max %d files can be written by one application "
                        "using the DATASPACES method\n",
                        __func__, MAX_NUM_OF_FILES);
            }
        }
        if (i < p->num_of_files) {
            p->fversions[i] = version;
        }


        // free allocated index lists
        adios_clear_index_v1 (index);
        adios_free_index_v1 (index);

        // rank=0 may be in put_sync when others call unlock, which is a global op
        MPI_Barrier (p->mpi_comm); 
        //log_debug("%s: call dspaces_put_sync()\n", __func__);
        //dspaces_put_sync();
        log_debug("%s: call dspaces_unlock_on_write(%s)\n", __func__, fd->name);
        dspaces_unlock_on_write(fd->name, &p->mpi_comm);
    }
    else if( fd->mode == adios_mode_read )
    {
        dspaces_unlock_on_read(fd->name, &p->mpi_comm);
    } 

    /* Increment the time index */
    info->time_index++;


    log_info ("%s: exit\n", __func__);
}

void adios_dataspaces_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_ds_data_struct *p = (struct adios_ds_data_struct *)
        method->method_data;
    int i;
    char ds_var_name[MAX_DS_NAMELEN];
    int lb[3] = {0,0,0}; 
    int ub[3] = {1,0,0}; // we put 2 integers to space, 
    int didx[3]; // for reordering DS dimensions
    int value[2] = {0, 1}; // integer to be written to space (terminated=1)

    free_dspaces_file_info(p);

    // tell the readers which files are finalized
    ds_dimension_ordering(1, 0, 0, didx); // C ordering of 1D array into DS
    for (i=0; i<p->num_of_files; i++) {
        /* Put VERSION@fn into space. Indicates that this file will not be extended anymore.  */
        log_debug("%s: call dspaces_lock_on_write(%s), rank=%d\n", __func__, p->fnames[i], mype);
        dspaces_lock_on_write(p->fnames[i], &p->mpi_comm); // lock is global operation in DataSpaces
        if (p->rank == 0) {
            value[0] = p->fversions[i];
            snprintf(ds_var_name, MAX_DS_NAMELEN, "VERSION@%s", p->fnames[i]);
            log_debug ("%s: update %s in the space [%d, %d]\n", 
                    __func__, ds_var_name, value[0], value[1] );
            dspaces_put(ds_var_name, 0, sizeof(int),   
                    lb[didx[0]], lb[didx[1]], lb[didx[2]], 
                    ub[didx[0]], ub[didx[1]], ub[didx[2]],  
                    &value); 
            log_debug("%s: call dspaces_put_sync()\n", __func__);
            dspaces_put_sync();
        }
        log_debug("%s: call dspaces_unlock_on_write(%s), rank=%d\n", __func__, p->fnames[i], mype);
        dspaces_unlock_on_write(p->fnames[i], &p->mpi_comm);
        free (p->fnames[i]);
    }

    // disconnect from dataspaces if we are connected from writer but not anymore from reader
    if (globals_adios_is_dataspaces_connected_from_writer() && 
            !globals_adios_is_dataspaces_connected_from_both())
    {
        log_debug ("%s: call dspaces_barrier(), rank=%d\n", __func__,mype);
        dspaces_barrier();
        log_debug ("%s: call dspaces_finalize(), rank=%d\n", __func__,mype);
        dspaces_finalize();

    }
    globals_adios_set_dataspaces_disconnected_from_writer();

    adios_dataspaces_initialized = 0;

    log_info("%s: exit\n", __func__);
}

void adios_dataspaces_end_iteration (struct adios_method_struct * method)
{
}

void adios_dataspaces_start_calculation (struct adios_method_struct * method)
{
}

void adios_dataspaces_stop_calculation (struct adios_method_struct * method)
{
}

