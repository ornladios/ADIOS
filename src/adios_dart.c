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



#if HAVE_DART
    void dart_init_(int *num_peers, int *appid);
    void dart_rank_(int *rank);
    void dart_peers_(int *peers);
    void dart_lock_on_write_();
    void dart_unlock_on_write_();
    void dart_lock_on_read_();
    void dart_unlock_on_read_();
    void dart_put_(char *var_name, unsigned int* ver, int *size,
                    int *xl, int *yl, int *zl,
                    int *xu, int *yu, int *zu,
                    void *data, int len);
    void dart_put_sync_();
    void dart_get_(char *var_name, unsigned int *ver, int *size,
                    int *xl, int *yl, int *zl,
                    int *xu, int *yu, int *zu,
                    void *data, int len);
    void dart_barrier_();
    void dart_finalize_();
#endif

static int adios_dart_initialized = 0;
static int adios_dart_connected_to_dart = 0;

struct adios_DART_data_struct
{
    int rank;
    int peers;  // from xml parameter or group communicator
    int appid;  // from xml parameter or 1
    int time_index; // versioning in DataSpaces, start from 1
    int n_writes; // how many times adios_write has been called
    int sync_at_close; // bool: whether to call a put_sync in adios_close
};


#if HAVE_DART
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
    
    index=0;
    for(i=0; *(parameters+index)!=','; index++,i++ )
    {
        temp[i] = *(parameters+index);
    }
    temp[i] = 0;
    //Get peers num info from parameters
    int num_peers = atoi(temp);

    for( index++, i=0; *(parameters+index)!=0; index++,i++ )
    {
      temp[i] = *(parameters+index);
    } 
    temp[i] = 0;
    //How to get app id?
    int appid = atoi(temp);

    //Init the static data structure
    p->peers = num_peers;
    p->appid = appid;
    p->time_index = 1;
    p->n_writes = 0;
    p->sync_at_close = 0;

    fprintf(stderr, "adios_dart_init:  param peers=%d, param appid=%d\n", p->peers, p->appid);
   
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

        if (!comm_name && !strcmp (comm_name, ""))
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
    if (!adios_dart_connected_to_dart) {

#if HAVE_MPI
    // if we have MPI and a communicator, we can get the exact size of this application
    // that we need to tell DART
    MPI_Comm group_comm;
    if (comm) {
        adios_dart_var_to_comm (fd->group->group_comm, fd->group->adios_host_language_fortran,
                                comm, &group_comm);
        MPI_Comm_size ( group_comm, &num_peers);
    }
#endif

        fprintf(stderr, "adios_dart_open: connect to DART, peers=%d, appid=%d \n",
                        num_peers, p->appid);

        //Init the dart client
        dart_init_(&num_peers, &(p->appid));
        dart_rank_(&(p->rank));
        dart_peers_(&(p->peers));

        fprintf(stderr, "adios_dart_open: connected: peers=%d,  rank=%d\n", p->peers, p->rank);        

        adios_dart_connected_to_dart = 1;
    }
   
    if( fd->mode == adios_mode_write )
    {
        fprintf(stderr, "adios_dart_open: rank=%d call write lock...\n", p->rank);        
        dart_lock_on_write_();        
        fprintf(stderr, "adios_dart_open: rank=%d got write lock\n", p->rank);        
    }
    else if( fd->mode == adios_mode_read )
    {
        dart_lock_on_read_();                    
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
    char dart_type_var_name[128];
    int  i0 = 0, i4 = 4;

    //Get two offset coordinate values
    unsigned int version;

    int dims[3]={1,1,1}, lb[3]={0,0,0}, ub[3]={0,0,0}; /* lower and upper bounds for DataSpaces */
    int index = 0;
    struct adios_dimension_struct* var_dimensions = v->dimensions;
    // Calculate lower and upper bounds for each available dimension (up to 3 dims)
    while( var_dimensions && index < 3)
    {
        dims[index] = get_dim_rank_value(&(var_dimensions->dimension), group);
        lb[index] = get_dim_rank_value(&(var_dimensions->local_offset), group);
        ub[index] = lb[index] + dims[index] - 1;
        var_dimensions = var_dimensions->next;
        index++;
    }

    /*version = p->time_index; */ /* Add new data as separate to DataSpaces */
    version = 0;                  /* Update/overwrite data in DataSpaces  (we write time_index as a variable at close)*/

    snprintf(dart_type_var_name, 128, "T%s", var_name);
    
    fprintf(stderr, "var_name=%s, type=%s(%d) elemsize=%d, version=%d, size_x=%d, size_y=%d, size_z=%d, (%d,%d,%d), (%d,%d,%d)\n",
                         var_name, adios_type_to_string_int(v->type), v->type, var_type_size, version,
                         dims[0], dims[1], dims[2], lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);

    /* Put type info as T<varname>, integer in 0,0,0,0,0,0 position */
    dart_put_(dart_type_var_name, &version, &i4, &i0, &i0, &i0, &i0, &i0, &i0, 
                &(v->type), strlen(dart_type_var_name)); 
    /* Flip 1st and 2nd dimension for DataSpaces representation */
    dart_put_(var_name, &version, &var_type_size, 
              lb+1, lb+0, lb+2,
              ub+1, ub+0, ub+2, 
              data, strlen(var_name));

#if 0
    p->n_writes++;
    /* have a sync after every 5th variable writing */
    printf("%s: This was write %d, mod5 = %d\n", __func__, p->n_writes, p->n_writes%5);
    if ( ! (p->n_writes % 5) ) {
        printf("%s: call dart_put_sync()\n", __func__);
        dart_put_sync_();
        p->sync_at_close = 0;
    } else {
        // no sync after this write, should sync at adios_close if this was the last write
        p->sync_at_close = 1;
    }
#else
    p->sync_at_close = 1;
#endif

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

    if( fd->mode == adios_mode_write )
    {
        /* Write two adios specific variables with the name of the file and name of the group into the space */
        /* ADIOS Read API fopen() checks these variables to see if writing already happened */
        unsigned int version = 0;
        int i4 = 4, i0 = 0;
        printf("%s: put %s = %d into space\n", __func__, fd->name, p->time_index);
        dart_put_(fd->name, &version, &i4, &i0, &i0, &i0, &i0, &i0, &i0, 
                  &(p->time_index), strlen(fd->name)); 
        printf("%s: put %s = %d into space\n", __func__, fd->group->name, p->time_index);
        dart_put_(fd->group->name, &version, &i4, &i0, &i0, &i0, &i0, &i0, &i0, 
                  &(p->time_index), strlen(fd->group->name)); 

        if (p->sync_at_close) {
            printf("%s: call dart_put_sync()\n", __func__);
            dart_put_sync_();
        }
        printf("%s: call dart_unlock_on_write()\n", __func__);
        dart_unlock_on_write_();        
    }
    else if( fd->mode == adios_mode_read )
    {
        dart_unlock_on_read_();                    
    } 

    /* Increment the time index */
    p->time_index++;


    printf("%s: exit\n", __func__);
}

void adios_dart_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_DART_data_struct *p = (struct adios_DART_data_struct *)
                                                method->method_data;
    if (adios_dart_connected_to_dart)
    {
        //TODO
        printf("%s: call dart_barrier()\n", __func__);
        dart_barrier_();
        printf("%s: call dart_finalize()\n", __func__);
        dart_finalize_();

        adios_dart_connected_to_dart = 0;
    }

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

#endif
