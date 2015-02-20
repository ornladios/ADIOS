/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

// see if we have MPI or other tools
#include "config.h"

// xml parser
#include <mxml.h>

#ifdef HAVE_XPMEM	        

//xpmem headers
#include "public/adios_xpmem.h"


#include "public/adios_mpi.h" // MPI or dummy MPI
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"
#include "public/adios_error.h"
#include "core/adios_logger.h"

//initialization checks
static int adios_xpmem_initialized = 0;

//xpmem helper functions
typedef struct _xmeminfo
{
	xpmem_segid_t buffer_id;
	uint64_t size;
	shared_data *sp; //pointer to the shared data
	char *b; //buffer that contains a copy of the shared data
}adios_xpmem_data_struct;





 
void adios_xpmem_init (const PairStruct * parameters
                      ,struct adios_method_struct * method
                      )
{
	adios_xpmem_data_struct *p = NULL;
    if (!adios_xpmem_initialized)
    {
        adios_xpmem_initialized = 1;
    }

    adios_logger_open(NULL, 0);
    
    method->method_data = malloc (sizeof (adios_xpmem_data_struct));
    memset(method->method_data, 0, sizeof(adios_xpmem_data_struct));

    //store the data struct in the method 
    p = (adios_xpmem_data_struct*)method->method_data;
    
    //now create the shared memory space
    //fake the buffer size to 10 MB
    p->buffer_id = make_share(&p->b, share_size);
    memset(p->b, 0, share_size);
    p->size = share_size;

    log_debug("xpmem initialized\tbuffer_id = %llu \t writing to disk",
              p->buffer_id);

    p->sp = (shared_data*)p->b;

   
    memset(p->sp, 0, sizeof(shared_data));

    p->sp->offset = (uint64_t)p->sp->buffer - (uint64_t)p->b;
        

    log_debug("xpmem data offset = %d",
              p->sp->offset);
        
    write_segid(p->buffer_id, "xpmem.data");

}



int adios_xpmem_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method, MPI_Comm comm
                     )
{
    char * name;
    adios_xpmem_data_struct * p = (adios_xpmem_data_struct *)
                                                          method->method_data;

    if(fd->mode == adios_mode_read || fd->mode == adios_mode_append)
    {
	    adios_error(err_operation_not_supported,
	                "xpmem does not support old read api\n");
	    return 0;
    }

    //check if there is a reader attached

    //wait for readcount to be 1
    if(p->sp->version != 0)
	    while(p->sp->readcount != 1)
		    adios_nanosleep(0, 100000000);

    while(p->sp->version != 0)
	    adios_nanosleep(0, 100000000);
    
    p->sp->readcount = 0;
    p->sp->version = 0;

    log_debug("xpmem_open completed\n");

    return 1;
}

enum ADIOS_FLAG adios_xpmem_should_buffer (struct adios_file_struct * fd
                                          ,struct adios_method_struct * method
                                          )
{
     adios_xpmem_data_struct * p = ( adios_xpmem_data_struct *)
                                                          method->method_data;

     log_debug("xpmem should buffer\n");
     
    //we alwasy return yes so we don't have to deal with the buffering shit
    return adios_flag_yes;
    
}

void adios_xpmem_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,void * data
                       ,struct adios_method_struct * method
                       )
{
     adios_xpmem_data_struct * p = ( adios_xpmem_data_struct *)
                                                          method->method_data;


    if(fd->shared_buffer != adios_flag_yes)
    {
	    log_error("xpmem relys on shared buffer\n");	    
    }

    // if(v->got_buffer != adios_flag_yes)
    // {
	//     log_error("xpmem relys on shared buffer\n");
    // }
    
}

//unknown functionality, just let it be

void adios_xpmem_get_write_buffer (struct adios_file_struct * fd
                                  ,struct adios_var_struct * v
                                  ,uint64_t * size
                                  ,void ** buffer
                                  ,struct adios_method_struct * method
                                  )
{
	adios_error(err_operation_not_supported,
	            "xpmem does not support get_write_buffer\n");
}

void adios_xpmem_read (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * buffer
                      ,uint64_t buffer_size
                      ,struct adios_method_struct * method
                      )
{
	log_error("xpmem does not support old read api\n");
	adios_error(err_operation_not_supported,
	            "xpmem does not support old read api\n");
}


void adios_xpmem_close (struct adios_file_struct * fd
                        ,struct adios_method_struct *method )
{
	adios_xpmem_data_struct * p = ( adios_xpmem_data_struct *)
		method->method_data;
	struct adios_attribute_struct * a = fd->group->attributes;
	struct adios_index_struct_v1 * index;
    
	switch (fd->mode)
	{
	case adios_mode_write:
	{
		// buffering or not, write the index
		char * buffer = 0;
		uint64_t buffer_size = 0;
		uint64_t buffer_offset = 0;
		uint64_t index_start = fd->base_offset + fd->offset;
		char *data_buffer;
		char *index_buffer;

		// build index
		index = adios_alloc_index_v1(1); // with hashtables
		adios_build_index_v1 (fd, index);
		adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset,
		                      index_start, index);
		adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

		//now copy the data buffer into the shared memory area
		data_buffer = &p->sp->buffer[0];
		//we can do a sanity check here
		log_info("data_buffer = %p\tb+offset = %p\n",
		         data_buffer, &p->b[p->sp->offset]);

		//copy the raw data
		memcpy(data_buffer, fd->buffer, fd->bytes_written);
		//copy the built index
		memcpy(&data_buffer[fd->bytes_written], buffer,
		       buffer_offset);

		//confirm buffer locations
		log_debug("bytes_writter = %d, buffer_offset = %d, total = %d\n",
		          fd->bytes_written, buffer_offset, 
		          fd->bytes_written + buffer_offset);
		
		//set the sizes for the data
		p->sp->size = fd->bytes_written + buffer_offset;
            
		//now set the version to 1
		ATOMIC_INCREMENT(p->sp->version);
		
		adios_free_index_v1(index);
		free (buffer);

		//now update the version variable
            

		break;
	}

	default:
	{
		fprintf (stderr, "Unsupported file mode: %d\n", fd->mode);

		return;
	}
	}

	// adios_posix_close_internal (&p->b);
	// adios_clear_index_v1 (index);

}

void adios_xpmem_finalize (int mype, struct adios_method_struct * method)
{
	adios_xpmem_data_struct * p = (adios_xpmem_data_struct *)
		method->method_data;

	//set the finalize flag on the shared segment
	p->sp->finalized = 1;
	
	//now loop over the readcount
	while(p->sp->readcount != 1)
		adios_nanosleep(0, 100000000);

	//now loop on finalized until it is set to 2
	while(p->sp->finalized != 2)
		adios_nanosleep(0, 100000000);

	//unmake the segment
	unmake_share(p->buffer_id, p->b);

	//now we can delete shit and exit
	free(p);
	
	if (adios_xpmem_initialized)
		adios_xpmem_initialized = 0;
}

void adios_xpmem_end_iteration (struct adios_method_struct * method)
{
}

void adios_xpmem_start_calculation (struct adios_method_struct * method)
{
}

void adios_xpmem_stop_calculation (struct adios_method_struct * method)
{
}
#endif	        
