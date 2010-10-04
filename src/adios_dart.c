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

// mpi
#if HAVE_MPI
#include "mpi.h"
#endif

#include "adios.h"
#include "adios_types.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"
#include "adios_internals_mxml.h"

#include "dart.h"

/*#define DART_DO_VERSIONING  define it at configure as -DDART_DO_VERSIONING in CFLAGS */


static int adios_dart_initialized = 0;
#define MAXDARTNAMELEN 128
static char dart_type_var_name[MAXDARTNAMELEN];
static char dart_var_name[MAXDARTNAMELEN];

struct adios_DART_data_struct
{
    int rank;   // dart rank or MPI rank if MPI is available
    int peers;  // from xml parameter or group communicator
    int appid;  // from xml parameter or 1
    int time_index; // versioning in DataSpaces, start from 0
    int n_writes; // how many times adios_write has been called
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

    fprintf(stderr, "adios_dart_init: appid=%d\n", p->appid);
   
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
                fprintf (stderr, "communicator not provided and none "
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
                fprintf (stderr, "communicator not provided but one "
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
        fprintf (stderr, "coordination-communication not provided. "
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
  
    fprintf(stderr, "adios_dart_open: open %s, mode=%d, time_index=%d \n",
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
        }
#endif

        fprintf(stderr, "adios_dart_open: rank=%d connect to DART, peers=%d, appid=%d \n",
                        p->rank, num_peers, p->appid);

        //Init the dart client
        ret = dart_init (num_peers, p->appid);
        if (ret) {
            fprintf(stderr, "adios_dart_open: rank=%d Failed to connect to DART: err=%d,  rank=%d\n", p->rank, ret);        
            return ret;
        }

#if ! HAVE_MPI
        dart_rank (&(p->rank));
        dart_peers (&(p->peers));
#endif

        fprintf(stderr, "adios_dart_open: rank=%d connected to DART: peers=%d\n", p->rank, p->peers);        
    }
    globals_adios_set_dart_connected_from_writer();
   
    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        fprintf(stderr, "adios_dart_open: rank=%d call write lock...\n", p->rank);        
        dart_lock_on_write (fd->name);  
        fprintf(stderr, "adios_dart_open: rank=%d got write lock\n", p->rank);        
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
    int ndims = 0;
    struct adios_dimension_struct* var_dimensions = v->dimensions;
    // Calculate lower and upper bounds for each available dimension (up to 3 dims)
    while( var_dimensions && ndims < 3)
    {
        dims[ndims] = get_dim_rank_value(&(var_dimensions->dimension), group);
        gdims[ndims] = get_dim_rank_value(&(var_dimensions->global_dimension), group);
        lb[ndims] = get_dim_rank_value(&(var_dimensions->local_offset), group);
        if (dims[ndims] > 0) 
            ub[ndims] = lb[ndims] + dims[ndims] - 1;
        else 
            ub[ndims] = lb[ndims]; // a time dimension is 0, so we need to handle this
        var_dimensions = var_dimensions->next;
        ndims++;
    }

#ifdef DART_DO_VERSIONING
    version = p->time_index;  /* Add new data as separate to DataSpaces */
#else
    version = 0;              /* Update/overwrite data in DataSpaces  (we write time_index as a variable at close)*/
#endif
    
    if (v->path != NULL && v->path[0] != '\0' && strcmp(v->path,"/")) 
        snprintf(dart_var_name, MAXDARTNAMELEN, "%s/%s/%s", fd->name, v->path, v->name);
    else 
        snprintf(dart_var_name, MAXDARTNAMELEN, "%s/%s", fd->name, v->name);

    snprintf(dart_type_var_name, MAXDARTNAMELEN, "TYPE@%s", dart_var_name);
    
    /* non-glboal variables are put in space ONLY by rank = 0 process */
    if (gdims[0] == 0 && p->rank != 0) {
        //fprintf(stderr, "rank=%d var_name=%s is not global. Skip\n", p->rank, dart_var_name);
        return;
    }

//    if (var_dimensions > 0 || p->rank == 0) {
    
        fprintf(stderr, "var_name=%s, type=%s(%d) elemsize=%d, version=%d, size=(%d,%d,%d), gdim=(%d,%d,%d), lb=(%d,%d,%d), ub=(%d,%d,%d)\n",
                dart_var_name, adios_type_to_string_int(v->type), v->type, var_type_size, version,
                dims[0], dims[1], dims[2], gdims[0], gdims[1], gdims[2], lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);

        /* Put type info as T<varname>, integer in 0,0,0,0,0,0 position */
        err = dart_put(dart_type_var_name, version, 4, 0,0,0,0,0,0, &(v->type)); 

        if (group->adios_host_language_fortran == adios_flag_yes || ndims == 1) {
            /* Flip 1st and 2nd dimension for DataSpaces representation
               for any Fortran writings and for any 1D array*/
            dart_put(dart_var_name, version, var_type_size, lb[1], lb[0], lb[2], ub[1], ub[0], ub[2], data);
        } else {
            /* Keep dimension order in case of C writer of 2-3D arrays for DataSpaces representation */
            dart_put(dart_var_name, version, var_type_size, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], data);
        }
//    }
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
            fprintf (stderr, "Out of memory allocating %llu bytes for %s\n"
                    ,*size, v->name
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
        fprintf (stderr, "OVERFLOW: Cannot allocate requested buffer of %llu "
                         "bytes for %s\n"
                ,*size
                ,v->name
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

void adios_dart_close (struct adios_file_struct * fd
                      ,struct adios_method_struct * method
                      )
{
    //TODO
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
                                                method->method_data;

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        if (p->rank == 0) {
            /* Write two adios specific variables with the name of the file and name of the group into the space */
            /* ADIOS Read API fopen() checks these variables to see if writing already happened */
            unsigned int version;
#ifdef DART_DO_VERSIONING
            version = p->time_index;  /* Add new data as separate to DataSpaces */
#else
            version = 0;              /* Update/overwrite data in DataSpaces */
#endif
            snprintf(dart_var_name, MAXDARTNAMELEN, "FILE@%s", fd->name);
            printf("%s: put %s = %d into space\n", __func__, dart_var_name, p->time_index);
            dart_put(dart_var_name, version, 4, 0, 0, 0, 0, 0, 0, &(p->time_index)); 

            snprintf(dart_var_name, MAXDARTNAMELEN, "GROUP@%s", fd->group->name);
            printf("%s: put %s = %d into space\n", __func__, dart_var_name, p->time_index);
            dart_put(dart_var_name, version, 4, 0, 0, 0, 0, 0, 0, &(p->time_index)); 
        }

        printf("%s: call dart_put_sync()\n", __func__);
        dart_put_sync();
        printf("%s: call dart_unlock_on_write(%s)\n", __func__, fd->name);
        dart_unlock_on_write(fd->name);
    }
    else if( fd->mode == adios_mode_read )
    {
        dart_unlock_on_read(fd->name);
    } 

    /* Increment the time index */
    p->time_index++;


    printf("%s: exit\n", __func__);
}

void adios_dart_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
                                                method->method_data;

    // disconnect from dart if we are connected from writer but not anymore from reader
    if (globals_adios_is_dart_connected_from_writer() && 
        !globals_adios_is_dart_connected_from_both())
    {
        printf("%s: call dart_barrier(), rank=%d\n", __func__,mype);
        dart_barrier();
        printf("%s: call dart_finalize(), rank=%d\n", __func__,mype);
        dart_finalize();

    }
    globals_adios_set_dart_disconnected_from_writer();

    adios_dart_initialized = 0;

    printf("%s: exit\n", __func__);
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

