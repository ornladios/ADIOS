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

#include "dimes.h"
/*#define DIMES_DO_VERSIONING define it at configure as -DDART_DO_VERSIONING in CFLAGS*/

static int adios_dimes_initialized = 0;
#define MAX_DIMES_NAMELEN 128
static char dimes_type_var_name[MAX_DIMES_NAMELEN];
static char dimes_var_name[MAX_DIMES_NAMELEN];

static int dimes_sync_id = 0;

struct adios_DIMES_data_struct
{
	int rank; //dart rank or MPI rank if MPI is available
	int peers; //from xml parameter or group communicator
	int appid; //from xml parameter or 1
	int time_index; //versioning in DIMES, start from 0
	int n_writes; //how many times adioes_write has been called
};

int adios_dimes_get_dim_rank_value(struct adios_dimension_item_struct * dim_info, struct adios_group_struct *group)
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

void adios_dimes_init(const char *parameters,
					struct adios_method_struct *method)
{
	struct adios_DIMES_data_struct *p = NULL;
	if(!adios_dimes_initialized){
		adios_dimes_initialized = 1;
	}
	
    	method->method_data = calloc (1, sizeof (struct adios_DIMES_data_struct));
    	p = (struct adios_DIMES_data_struct*)method->method_data;
	
	int index, i;
	char temp[64];
	int num_peers;
	int appid;
	int was_set;
	
	num_peers = 1; //Init value
	//application ID should be set by the application calling adios_set_application_id()
	appid = globals_adios_get_application_id(&was_set);
	if(!was_set)
		appid = 1;
	
	//init the static data structure
	p->peers = num_peers;
	p->appid = appid;
	p->time_index = 0;
	p->n_writes = 0;
	
	fprintf(stderr, "adios_dimes_init: appid=%d, peers=%d\n", p->appid, p->peers);
   
}

static void adios_dimes_var_to_comm(const char *comm_name,
								enum ADIOS_FLAG host_language_fortran,
								void *data,
								MPI_Comm *comm)
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

int adios_dimes_open(struct adios_file_struct *fd,
					struct adios_method_struct *method,
					void *comm)
{
	int ret = 0;
	struct adios_DIMES_data_struct *p = (struct adios_DIMES_data_struct*)
									method->method_data;
	int num_peers = p->peers;
	
        fprintf(stderr, "adios_dimes_open: open %s, mode=%d, time_index=%d \n",
                        fd->name, fd->mode, p->time_index);	
	
#if HAVE_MPI
		//Get applicaitons size / process rank info from MPI communicator,
		//which is necessary for dimes
		MPI_Comm group_comm;
		if(comm) {
			adios_dimes_var_to_comm(fd->group->group_comm,fd->group->adios_host_language_fortran,
									comm, &group_comm);
			MPI_Comm_rank(group_comm, &(p->rank));
			MPI_Comm_size(group_comm, &num_peers);
			p->peers = num_peers;
		}
#endif
	
	//connect to DIMES index srv at the very first adios_open(),disconnect in adios_finalize()
	if(!globals_adios_is_dimes_connected()){
		//Init dimes client
		//int num_total_peers = 64+16+1;
		ret = dimes_init(num_peers,num_peers,p->appid);
		if(ret<0){
			fprintf(stderr, "adios_dimes_open: rank=%d Failed to connect to index srv: err=%d,  rank=%d\n", p->rank, ret);        
            		return ret;
		}
		
#if !HAVE_MPI
		p->rank = dimes_rank();
		p->peers = dimes_peers();
#endif

		fprintf(stderr,"adios_dimes_open:rank=%d connected to DIMES index srv: peers=%d, appid=%d\n",
			p->rank, p->peers, p->appid);
	}
	globals_adios_set_dimes_connected_from_writer();
	
	return ret;
}

enum ADIOS_FLAG adios_dimes_should_buffer(struct adios_file_struct *fd,
										struct adios_method_struct *method)
{
	return adios_flag_no;
}

void adios_dimes_write(struct adios_file_struct *fd,
			struct adios_var_struct *v,
			void *data,
			struct adios_method_struct *method)
{
	struct adios_DIMES_data_struct *p = (struct adios_DIMES_data_struct *)
							method->method_data;
	
	struct adios_group_struct *group = fd->group;
	
	//Get var size
	int var_type_size = (int)adios_get_type_size(v->type,v->data);
	//Get var name
	char *var_name = v->name;
	int err;
	
	//Get two offset coordinate values
	unsigned int version;
	
	int dims[3]={1,1,1},lb[3]={0,0,0},ub[3]={0,0,0};/*lower and upper bounds*/
	int ndims = 0;
	struct adios_dimension_struct *var_dimensions = v->dimensions;
	//Calculate lower and upper bounds for each var(up to 3 dims)
	while( var_dimensions && ndims<3 ){
		dims[ndims] = adios_dimes_get_dim_rank_value(&(var_dimensions->dimension), group);
		lb[ndims] = adios_dimes_get_dim_rank_value(&(var_dimensions->local_offset),group);
		ub[ndims] = lb[ndims] + dims[ndims] - 1;
		var_dimensions = var_dimensions->next;
		ndims++;
	}
	
#ifdef DIMES_DO_VERSIONING
	version = p->time_index;
#else
	version = 0;
#endif

	//Construct var name and var type name	
	snprintf(dimes_var_name, MAX_DIMES_NAMELEN, "%s/%s", fd->name, v->name);
	snprintf(dimes_type_var_name, MAX_DIMES_NAMELEN, "TYPE@%s", dimes_var_name);

	/*Scalar variables are put into space only by rank=0 process*/
	if(p->rank == 0){
		err = dimes_put_scalar(dimes_type_var_name,version,4,0,0,0,0,0,0,v->type);
	}

	/*	
	if(group->adios_host_language_fortran == adios_flag_yes || ndims == 1){
		//Flip 1st and 2nd dimension for any Fortran writings and any 1D array
		dimes_put(dimes_var_name, version, var_type_size,lb[1],lb[0],lb[2],ub[1],ub[0],ub[2],data,1);
	} else {
		//Keep dimension order in case of C writer of 2-3D arrays
		dimes_put(dimes_var_name, version, var_type_size,lb[0],lb[1],lb[2],ub[0],ub[1],ub[2],data,1);
	}
	*/
	if( lb[0]==0 && lb[1]==0 && lb[2]==0 && ub[0]==0 && ub[1]==0 && ub[2]==0){
		if(p->rank == 0){
			/*Scalar variables are put into space only by rank=0 process*/
			int scalar_data = *(int *)data;
			dimes_put_scalar(dimes_var_name, version, 4, 0,0,0,0,0,0, scalar_data);

			fprintf(stderr, "var_name=%s, type=%s(%d) elemsize=%d, version=%d, size_x=%d, size_y=%d, size_z=%d, (%d,%d,%d), (%d,%d,%d)\n",
				dimes_var_name, adios_type_to_string_int(v->type), v->type, var_type_size, version,
				dims[0], dims[1], dims[2], lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);
		}
	}
	else {
		/*Multi-dimension(1-D,2-D or 3-D) array is put into space*/
		if(group->adios_host_language_fortran == adios_flag_yes || ndims == 1){
			//Flip 1st and 2nd dimension for any Fortran writings and any 1D array
			dimes_put(dimes_var_name, version, var_type_size,lb[1],lb[0],lb[2],ub[1],ub[0],ub[2],data,1);
		} else {
			//Keep dimension order in case of C writer of 2-3D arrays
			dimes_put(dimes_var_name, version, var_type_size,lb[0],lb[1],lb[2],ub[0],ub[1],ub[2],data,1);
		}
	
		/*
		fprintf(stderr, "var_name=%s, type=%s(%d) elemsize=%d, version=%d, size_x=%d, size_y=%d, size_z=%d, (%d,%d,%d), (%d,%d,%d)\n",
			dimes_var_name, adios_type_to_string_int(v->type), v->type, var_type_size, version,
			dims[0], dims[1], dims[2], lb[0], lb[1], lb[2], ub[0], ub[1], ub[2]);
		*/
	}
}

void adios_dimes_get_write_buffer(struct adios_file_struct *fd,
				struct adios_var_struct *v,
				uint64_t *size,
				void **buffer,
				struct adios_method_struct *method)
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

/*NOT IMPLEMENTED. Use READ API to read variables*/
void adios_dimes_read(struct adios_file_struct *fd,
					struct adios_var_struct *v, void *buffer,
					uint64_t buffer_size,
					struct adios_method_struct *method)
{
	//NOT TO DO
}

void adios_dimes_close(struct adios_file_struct *fd,
					struct adios_method_struct *method)
{
	struct adios_DIMES_data_struct *p = (struct adios_DIMES_data_struct *)
											method->method_data;

	unsigned int version;	
#ifdef DIMES_DO_VERSIONING
	version = p->time_index;
	dimes_sync_id = version;
#else
	version = 0;
	dimes_sync_id++;
#endif


	if(fd->mode == adios_mode_write)
	{
	    //put file name & group name scalar info
       	    if (p->rank == 0) {
            /* Write two adios specific variables with the name of the file and name of the group into the space */
            /* ADIOS Read API fopen() checks these variables to see if writing already happened */
            snprintf(dimes_var_name, MAX_DIMES_NAMELEN, "FILE@%s", fd->name);
            printf("%s: put %s = %d into space\n", __func__, dimes_var_name, p->time_index);
            dimes_put_scalar(dimes_var_name, version, 4, 0, 0, 0, 0, 0, 0, p->time_index);

	    }
	}
	else if(fd->mode == adios_mode_read)
	{
		fprintf(stderr,"adios_dimes_open: should not happens when fd->mode==adios_mode_read\n");
	}

	//Synchronization
	if( fd->mode == adios_mode_write )
	{
		//printf("%s: call dimes_sync()\n",__func__);
		dimes_sync(dimes_sync_id);
		//dimes_wait();	
	}
	else if( fd->mode == adios_mode_read )
	{
		fprintf(stderr,"adios_dimes_close: should not happens when fd->mode==adios_mode_read\n");
	}
	
	/*Incrementa of time index*/
	p->time_index++;
	
	printf("%s: exit\n", __func__);
}

void adios_dimes_finalize(int mype, struct adios_method_struct *method)
{
	struct adios_DIMES_data_struct *p =(struct adios_DIMES_data_struct *)
											method->method_data;
	
	//disconnect from dimes if it is connected from writer but not anymore from reader??
	if(globals_adios_is_dimes_connected_from_writer() &&
		!globals_adios_is_dimes_connected_from_both())
	{
		printf("%s: call dimes_barrier(), rank=%d\n",__func__,mype);
		dimes_barrier();
		printf("%s: call dimes_finalize(), rank=%d\n",__func__,mype);
		dimes_finalize();
	}
	globals_adios_set_dimes_disconnected_from_writer();
	
	adios_dimes_initialized = 0;
	
	printf("%s: exit\n",__func__);
}

void adios_dimes_end_iteration(struct adios_method_struct *method)
{
}

void adios_dimes_start_calculation (struct adios_method_struct * method)
{
}

void adios_dimes_stop_calculation (struct adios_method_struct * method)
{
}
