#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// xml parser
#include <mxml.h>

// dart
#include <sys/uio.h>

// see if we have MPI or other tools
#include "config.h"

#include "adios.h"
#include "adios_types.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"
#include "adios_internals_mxml.h"
#include "ds_metadata.h"

#include "dart.h"

/*#define DART_DO_VERSIONING  define it at configure as -DDART_DO_VERSIONING in CFLAGS */

#define adios_logger(verbose_level, ...) if (adios_dart_verbose >= verbose_level) fprintf (stderr, __VA_ARGS__); 

#define log_error(...) adios_logger(0, __VA_ARGS__)
#define log_warn(...) adios_logger(1, __VA_ARGS__)
#define log_info(...) adios_logger(2, __VA_ARGS__)
#define log_debug(...) adios_logger(3, __VA_ARGS__)

static int adios_dart_initialized = 0;
#define MAXDARTNAMELEN 128
//static char dart_type_var_name[MAXDARTNAMELEN];
static char dart_var_name[MAXDARTNAMELEN];
static unsigned int adios_dart_verbose = 3;

struct adios_DART_data_struct
{
    int rank;   // dart rank or MPI rank if MPI is available
    int peers;  // from xml parameter or group communicator
    int appid;  // from xml parameter or 1
    int time_index; // versioning in DataSpaces, start from 0
    int n_writes; // how many times adios_write has been called
#if HAVE_MPI
    MPI_Comm mpi_comm;
#endif
    int  num_of_files; // how many files do we have with this method
    char *fnames[20];  // names of files (needed at finalize)
};



int get_dim_rank_value(struct adios_dimension_item_struct * dim_info, struct adios_group_struct *group)
{
    if(!dim_info)
       return 0;

    if(dim_info->id != 0)
    {
        struct adios_var_struct *dim_var = NULL;
        dim_var = adios_find_var_by_id(group->vars, dim_info->id);
        if(!dim_var || !dim_var->data)
           return 0;
        
        int rank = 0;
        switch ( dim_var->type )
        {
        case adios_unsigned_byte:
             rank = (*(uint8_t*) dim_var->data);
             break;
        case adios_byte:
             rank = (*(int8_t*) dim_var->data);
             break;
        case adios_unsigned_short:
             rank = (*(uint16_t*) dim_var->data);
             break;
        case adios_short:
             rank = (*(int16_t*) dim_var->data);
             break;
        case adios_unsigned_integer:
             rank = (*(uint32_t*) dim_var->data);
             break;
        case adios_integer:
             rank = (*(int32_t*) dim_var->data);
             break;
        case adios_unsigned_long:
             rank = (*(uint64_t*) dim_var->data);
             break;
        case adios_long:
             rank = (*(int64_t*) dim_var->data);
             break;

        default: break;
        }

        return rank;
    }
    else
    {
        return dim_info->rank;
    }
}

void adios_dart_init (const char * parameters,
                     struct adios_method_struct * method
                     )
{
    struct adios_DART_data_struct *p = 0;
    if (!adios_dart_initialized)
    {
        adios_dart_initialized = 1;
    }
   
    method->method_data = calloc (1, sizeof (struct adios_DART_data_struct));
    p = (struct adios_DART_data_struct*)method->method_data;
    
    int index, i;
    char temp[64];
    int num_peers;
    int appid;
    int was_set;
    
    num_peers = 1;
    // Application ID should be set by the application calling adios_set_application_id()
    appid = globals_adios_get_application_id (&was_set);
    if (!was_set)
        appid = 1;

    //Init the static data structure
    p->peers = num_peers;
    p->appid = appid;
    p->time_index = 0;
    p->n_writes = 0;
#if HAVE_MPI
    p->mpi_comm = MPI_COMM_NULL;
#endif
    p->num_of_files = 0;

    log_info ("adios_dart_init: appid=%d\n", p->appid);
   
}


static void adios_dart_var_to_comm  (const char * comm_name
                                    ,enum ADIOS_FLAG host_language_fortran
                                    ,void * data
                                    ,MPI_Comm * comm
                                    )
{
    if (data)
    {
        int t = *(int *) data;

        if (!comm_name || !strcmp (comm_name, ""))
        {
            if (!t)
            {
                log_error ("ERROR: communicator not provided and none "
                                 "listed in XML.  Defaulting to "
                                 "MPI_COMM_SELF\n"
                        );

                *comm = MPI_COMM_SELF;
            }
            else
            {
                if (host_language_fortran == adios_flag_yes)
                {
                    *comm = MPI_Comm_f2c (t);
                }
                else
                {
                    *comm = *(MPI_Comm *) data;
                }
            }
        }
        else
        {
            if (!t)
            {
                log_error ("ERROR: communicator not provided but one "
                                 "listed in XML.  Defaulting to "
                                 "MPI_COMM_WORLD\n"
                        );

                *comm = MPI_COMM_WORLD;
            }
            else
            {
                if (host_language_fortran == adios_flag_yes)
                {
                    *comm = MPI_Comm_f2c (t);
                }
                else
                {
                    *comm = *(MPI_Comm *) data;
                }
            }
        }
    }
    else
    {
        log_error ("ERROR: coordination-communication not provided. "
                         "Using MPI_COMM_WORLD instead\n"
                );

        *comm = MPI_COMM_WORLD;
    }
}

int adios_dart_open (struct adios_file_struct * fd,
                    struct adios_method_struct * method,
                    void *comm
                    )
{
    int ret = 0;
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
                                                method->method_data;
    int num_peers = p->peers;
  
    log_info ("adios_dart_open: open %s, mode=%d, time_index=%d \n",
                        fd->name, fd->mode, p->time_index);

    // connect to DART at the very first adios_open(), disconnect in adios_finalize()
    // connect only if the READ API has not connected yet
    if (!globals_adios_is_dart_connected()) {

#if HAVE_MPI
    // if we have MPI and a communicator, we can get the exact size of this application
    // that we need to tell DART
        MPI_Comm group_comm;
        if (comm) {
            adios_dart_var_to_comm (fd->group->group_comm, fd->group->adios_host_language_fortran,
                                    comm, &group_comm);
            MPI_Comm_rank ( group_comm, &(p->rank));
            MPI_Comm_size ( group_comm, &num_peers);
            p->peers = num_peers;
            p->mpi_comm = group_comm;
        }
#endif

        log_debug ("adios_dart_open: rank=%d connect to DART, peers=%d, appid=%d \n",
                        p->rank, num_peers, p->appid);

        //Init the dart client
        ret = dart_init (num_peers, p->appid);
        if (ret) {
            log_error ("adios_dart_open: rank=%d Failed to connect to DART: err=%d,  rank=%d\n", p->rank, ret);        
            return ret;
        }

#if ! HAVE_MPI
        dart_rank (&(p->rank));
        dart_peers (&(p->peers));
#endif

        log_debug ("adios_dart_open: rank=%d connected to DART: peers=%d\n", p->rank, p->peers);        
    }
    globals_adios_set_dart_connected_from_writer();
   
    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        log_debug ("adios_dart_open: rank=%d call write lock...\n", p->rank);        
        dart_lock_on_write (fd->name);  
        log_debug ("adios_dart_open: rank=%d got write lock\n", p->rank);        
    }
    else if (fd->mode == adios_mode_read)
    {
        dart_lock_on_read (fd->name);
    } 
  
    return ret;
}

enum ADIOS_FLAG adios_dart_should_buffer (struct adios_file_struct * fd
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


void adios_dart_write (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * data
                      ,struct adios_method_struct * method
                      )
{
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
                                                            method->method_data;
    struct adios_group_struct *group = fd->group;
    //Get var size
    //  FIXME: type size of a string >2GB does not fit to int. 
    //  adios_get_type_size returns uint64_t but dart_put handles only int
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
        dims[ndims] = get_dim_rank_value(&(var_dimensions->dimension), group);
        gdims[ndims] = get_dim_rank_value(&(var_dimensions->global_dimension), group);
        lb[ndims] = get_dim_rank_value(&(var_dimensions->local_offset), group);
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

#ifdef DART_DO_VERSIONING
    version = p->time_index;  /* Add new data as separate to DataSpaces */
#else
    version = 0;              /* Update/overwrite data in DataSpaces  (we write time_index as a variable at close)*/
#endif
    
    if (v->path != NULL && v->path[0] != '\0' && strcmp(v->path,"/")) 
        snprintf(dart_var_name, MAXDARTNAMELEN, "%s/%s/%s/%s", fd->name, fd->group->name, v->path, v->name);
    else 
        snprintf(dart_var_name, MAXDARTNAMELEN, "%s/%s//%s", fd->name, fd->group->name, v->name);

    //snprintf(dart_type_var_name, MAXDARTNAMELEN, "TYPE@%s", dart_var_name);
    
    /* non-global variables are put in space ONLY by rank = 0 process */
    if (gdims[0] == 0 && p->rank != 0) {
        //fprintf(stderr, "rank=%d var_name=%s is not global. Skip\n", p->rank, dart_var_name);
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
            dart_var_name, adios_type_to_string_int(v->type), v->type, var_type_size, version, ndims,
            dims[0], dims[1], dims[2], gdims[0], gdims[1], gdims[2], lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);

    /* non-timed scalars are written in the metadata at close(), not here */
    if (ndims == 0 && !hastime)
        return;

    /* Put type info as T<varname>, integer in 0,0,0,0,0,0 position */
    //err = dart_put(dart_type_var_name, version, 4, 0,0,0,0,0,0, &(v->type)); 

    ds_dimension_ordering(ndims,
            group->adios_host_language_fortran == adios_flag_yes, 
            0 /*pack*/, didx);

    dart_put(dart_var_name, version, var_type_size, 
             lb[didx[0]], lb[didx[1]], lb[didx[2]], 
             ub[didx[0]], ub[didx[1]], ub[didx[2]], 
             data);
    
    log_debug ("var_name=%s, dimension ordering=(%d,%d,%d), gdims=(%d,%d,%d), lb=(%d,%d,%d), ub=(%d,%d,%d)\n",
            dart_var_name, 
            didx[0], didx[1], didx[2], 
            gdims[didx[0]], gdims[didx[1]], gdims[didx[2]], 
            lb[didx[0]], lb[didx[1]], lb[didx[2]], 
            ub[didx[0]], ub[didx[1]], ub[didx[2]]);
}

void adios_dart_get_write_buffer (struct adios_file_struct * fd
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
void adios_dart_read (struct adios_file_struct * fd
                     ,struct adios_var_struct * v, void * buffer
                     ,uint64_t buffer_size
                     ,struct adios_method_struct * method
                     )
{
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
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

    version = p->time_index;

    //dart_lock_on_read_();

    //dart_get

    //dart_unlock_on_read_();
}

/* Gather var/attr indices from all processes to rank 0 */
static void adios_dart_gather_indices (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               ,struct adios_index_process_group_struct_v1 **pg_root 
                               ,struct adios_index_var_struct_v1 **vars_root
                               ,struct adios_index_attribute_struct_v1 ** attrs_root
                               )
{
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
                                                method->method_data;
    struct adios_index_process_group_struct_v1 * my_pg_root = 0;
    struct adios_index_var_struct_v1 * my_vars_root = 0;
    struct adios_index_attribute_struct_v1 * my_attrs_root = 0;
    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;
    
    // build local index first appending to any existing index
    adios_build_index_v1 (fd, &my_pg_root, &my_vars_root, &my_attrs_root);

    log_debug ("%s index after first build is pg=%x vars=%x attrs=%x\n", 
                __func__, my_pg_root, my_vars_root, my_attrs_root);
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
                adios_parse_vars_index_v1 (&b, &new_vars_root);
                adios_parse_attributes_index_v1 (&b
                        ,&new_attrs_root
                        );
                adios_merge_index_v1 (&my_pg_root
                        ,&my_vars_root
                        ,&my_attrs_root
                        ,new_pg_root, new_vars_root
                        ,new_attrs_root
                        );
                adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);
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
                    ,0, my_pg_root ,my_vars_root ,my_attrs_root);

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

    *pg_root = my_pg_root;
    *vars_root = my_vars_root;
    *attrs_root = my_attrs_root;
    log_debug ("%s index after gathering is pg=%x vars=%x attrs=%x\n", 
                __func__, my_pg_root, my_vars_root, my_attrs_root);
}

static int ds_get_full_name (char * path, char * name, int maxlen,
                            /*OUT*/char * out)
{
    int len;
    // make full name
    if (!path || !path[0] || !strcmp (path, "/")) { 
        // no path, just name + leading /
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

void ds_pack_group_info (struct adios_file_struct *fd
                                  ,struct adios_method_struct * method
                                  ,struct adios_index_var_struct_v1 *vars_root
                                  ,struct adios_index_attribute_struct_v1 * attrs_root
                                  ,char ** buffer, int *buffer_size, int *nvars, int *nattrs
                                  )
{
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
                                                method->method_data;
    struct adios_index_var_struct_v1 * v = vars_root;
    struct adios_index_attribute_struct_v1 * a = attrs_root;
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
                + strlen(v->var_name) + strlen(v->var_path) + 1  // full path
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
        size += sizeof(int) + strlen(a->attr_name) + strlen(a->attr_path) + 1 // name len + name
                + sizeof(int); // type
        if (a->type != adios_string)
            size += adios_get_type_size(a->type, NULL);
        else
            size += adios_get_type_size(a->type, a->characteristics->value) + sizeof(int);

        log_debug (" attr %s/%s, size = %d\n", a->attr_path, a->attr_name, size);
        (*nattrs)++;
        a = a->next;
    }

    *buffer = (char *) malloc (size);
    *buffer_size = size;

    /* Second cycle: fill up the buffer */
    v = vars_root;
    a = attrs_root;
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
            memcpy (b, &(a->characteristics->value), size); 
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
            log_debug("%3.3d %3.3d %3.3d %3.3d    ", b[i+4*j], b[i+4*j+1], b[i+4*j+2], b[i+4*j+3]);
        }
        log_debug("\n");
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

void adios_dart_close (struct adios_file_struct * fd
                      ,struct adios_method_struct * method
                      )
{
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
                                                method->method_data;
    struct adios_index_process_group_struct_v1 * pg_root;
    struct adios_index_var_struct_v1 * vars_root;
    struct adios_index_attribute_struct_v1 * attrs_root;
    struct adios_attribute_struct * a = fd->group->attributes;
    int lb[3], ub[3], didx[3]; // for reordering DS dimensions

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        // finalize variable info in fd buffer, next we call build_index
        while (a) {
            a->write_offset = 1; // only attributes with !=0 offset will be included in build index
            a=a->next;
        }

        //adios_write_close_vars_v1 (fd);
        /* Gather var/attr indices from all processes to rank 0 */
        adios_dart_gather_indices (fd, method, &pg_root, &vars_root ,&attrs_root);

        if (p->rank == 0) {

            /* Write two adios specific variables with the name of the file and name of the group into the space */
            /* ADIOS Read API fopen() checks these variables to see if writing already happened */
            unsigned int version;
#ifdef DART_DO_VERSIONING
            version = p->time_index;  /* Add new data as separate to DataSpaces */
#else
            version = 0;              /* Update/overwrite data in DataSpaces */
#endif

            /* Make metadata from indices */
            char * indexbuf;
            int    indexlen;
            int    nvars, nattrs;
            ds_pack_group_info (fd, method, vars_root, attrs_root, 
                                   &indexbuf, &indexlen, &nvars, &nattrs);

            
            /* Put GROUP@fn/gn header into space */
            snprintf(dart_var_name, MAXDARTNAMELEN, "GROUP@%s/%s", fd->name, fd->group->name);
            log_debug ("%s: put %s with buf len %d into space\n", __func__, dart_var_name, indexlen);
            ub[0] = indexlen-1; ub[1] = 0; ub[2] = 0;
            ds_dimension_ordering(1, 0, 0, didx); // C ordering of 1D array into DS
            dart_put(dart_var_name, version, 1,    0, 0, 0, /* lb 0..2 */
                     ub[didx[0]], ub[didx[1]], ub[didx[2]],  indexbuf); 

            /* Create and put FILE@fn header into space */
            char * file_info_buf; /* store FILE@fn's group list */
            int    file_info_buf_len; /* = 128 currently */
            snprintf (dart_var_name, MAXDARTNAMELEN, "FILE@%s", fd->name);
            ds_pack_file_info (p->time_index, nvars, nattrs, indexlen, fd->group->name, 
                               &file_info_buf, &file_info_buf_len);
            log_debug ("%s: put %s = buflen=%d time=%d nvars=%d nattr=%d index=%d name=%d:%s into space\n",
                __func__, dart_var_name, 
                *(int*)file_info_buf, *(int*)(file_info_buf+4), 
                *(int*)(file_info_buf+8), *(int*)(file_info_buf+12),
                *(int*)(file_info_buf+16), *(int*)(file_info_buf+20),
                file_info_buf+24);
            /* Flip 1st and 2nd dimension for DataSpaces representation for a 1D array*/
            ub[0] = file_info_buf_len-1; ub[1] = 0; ub[2] = 0;
            ds_dimension_ordering(1, 0, 0, didx); // C ordering of 1D array into DS
            dart_put(dart_var_name, version, 1,    0, 0, 0, /* lb 0..2 */
                     ub[didx[0]], ub[didx[1]], ub[didx[2]], file_info_buf); 

            free (indexbuf);

            // remember this filename for finalize
            int i;
            for (i=0; i<p->num_of_files; i++) {
                if (!strcmp(fd->name, p->fnames[i]))
                    break;
            }
            if (i == p->num_of_files) {
                if (p->num_of_files < 20) {
                    p->fnames[ p->num_of_files ] = strdup(fd->name);
                    p->num_of_files++;
                } else {
                    log_error ("%s: Max 20 files can be written by one application using the DART method\n",__func__);
                }
            }
        }

        // free allocated index lists
        adios_clear_index_v1 (pg_root, vars_root, attrs_root);

        log_debug("%s: call dart_put_sync()\n", __func__);
        dart_put_sync();
        log_debug("%s: call dart_unlock_on_write(%s)\n", __func__, fd->name);
        dart_unlock_on_write(fd->name);
    }
    else if( fd->mode == adios_mode_read )
    {
        dart_unlock_on_read(fd->name);
    } 

    /* Increment the time index */
    p->time_index++;


    log_info ("%s: exit\n", __func__);
}

void adios_dart_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
        method->method_data;
    int i;
    char dart_var_name[MAXDARTNAMELEN];

    if (p->rank == 0) {
        // tell the readers which files are finalized
        //dart_lock_on_write("CLOSE");
        for (i=0; i<p->num_of_files; i++) {
            /* Put CLOSED@fn into space. Indicates that this file will not be extended anymore. */
            snprintf(dart_var_name, MAXDARTNAMELEN, "CLOSED@%s", p->fnames[i]);
            log_debug ("%s: put %s with %d bytes into the space\n", __func__, dart_var_name, sizeof(int));
            dart_put(dart_var_name, 0, sizeof(int),    0, 0, 0,    0, 0, 0, &i); 
            free (p->fnames[i]);
        }
        //dart_put_sync();
        //dart_unlock_on_write("CLOSE");
    }

    // disconnect from dart if we are connected from writer but not anymore from reader
    if (globals_adios_is_dart_connected_from_writer() && 
            !globals_adios_is_dart_connected_from_both())
    {
        log_debug ("%s: call dart_barrier(), rank=%d\n", __func__,mype);
        dart_barrier();
        log_debug ("%s: call dart_finalize(), rank=%d\n", __func__,mype);
        dart_finalize();

    }
    globals_adios_set_dart_disconnected_from_writer();

    adios_dart_initialized = 0;

    log_info("%s: exit\n", __func__);
}

void adios_dart_end_iteration (struct adios_method_struct * method)
{
}

void adios_dart_start_calculation (struct adios_method_struct * method)
{
}

void adios_dart_stop_calculation (struct adios_method_struct * method)
{
}

