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
	shared_data *d;
	char *buffer;
}xmeminfo;


typedef struct _adios_xpmem_data_struct
{
	xmeminfo data;
	xmeminfo index;
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

    adios_logger_open("xpmem.log", 0);
    
    method->method_data = malloc (sizeof (adios_xpmem_data_struct));
    memset(method->method_data, 0, sizeof(adios_xpmem_data_struct));

    p = (adios_xpmem_data_struct*)method->method_data;
    
    //now create the shared memory space
    //fake the buffer size to 10 MB
    p->data.buffer_id = make_share(&p->data.buffer, share_size);
    p->index.buffer_id = make_share(&p->index.buffer, index_share_size);

    p->data.size = share_size;
    p->index.size = index_share_size;

    log_debug("xpmem initialized\tbuffer_id = %llu index_id = %llu\nwriting to disk",
              p->data.buffer_id, p->index.buffer_id);

    p->data.d = (shared_data*)p->data.buffer;
    p->index.d = (shared_data*)p->index.buffer;

    p->data.d->offset = (uint64_t)p->data.d->buffer - (uint64_t)p->data.buffer;
    p->index.d->offset = (uint64_t)p->index.d->buffer - (uint64_t)p->index.buffer;

    memset(p->data.d, 0, sizeof(shared_data));
    memset(p->index.d, 0, sizeof(shared_data));

    log_debug("xpmem data offset = %d\t index offset = %d\n",
              p->data.d->offset, p->index.d->offset);
    
    
    write_segid(p->data.buffer_id, "xpmem.data");
    write_segid(p->index.buffer_id, "xpmem.index");

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
    if(p->data.d->version != 0)
    while(p->data.d->readcount != 1)
	    adios_nanosleep(0, 100000000);

    if(p->index.d->version != 0)
    while(p->index.d->readcount != 1)
	    adios_nanosleep(0, 100000000);

    p->data.d->readcount = 0;
    p->index.d->readcount = 0;

    p->data.d->version = 0;
    p->index.d->version = 0;
    
    log_debug("xpmem_open completed\n");

    return 1;
}

enum ADIOS_FLAG adios_xpmem_should_buffer (struct adios_file_struct * fd
                                          ,struct adios_method_struct * method
                                          )
{
     adios_xpmem_data_struct * p = ( adios_xpmem_data_struct *)
                                                          method->method_data;

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

    if(v->got_buffer != adios_flag_yes)
    {
	    log_error("xpmem relys on shared buffer\n");
    }
    
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

// static void adios_xpmem_do_write (struct adios_file_struct * fd
//                                  ,struct adios_method_struct * method
//                                  ,char * buffer
//                                  ,uint64_t buffer_size
//                                  )
// {
//     struct adios_XPMEM_data_struct * p = (struct adios_XPMEM_data_struct *)
//                                                           method->method_data;
//     int32_t to_write;
//     uint64_t bytes_written = 0;

//     if (fd->shared_buffer == adios_flag_yes)
//     {
//         lseek (p->b.f, p->b.end_of_pgs, SEEK_SET);
//         if (p->b.end_of_pgs + fd->bytes_written > fd->pg_start_in_file + fd->write_size_bytes)
//             fprintf (stderr, "adios_xpmem_write exceeds pg bound. File is corrupted. "
//                              "Need to enlarge group size. \n");

//         if (fd->bytes_written > MAX_MPIWRITE_SIZE)
//         {
//             to_write = MAX_MPIWRITE_SIZE;
//         }
//         else
//         {
//             to_write = (int32_t) fd->bytes_written;
//         }

//         while (bytes_written < fd->bytes_written)
//         {
//             write (p->b.f, fd->buffer, to_write);
//             bytes_written += to_write;
//             if (fd->bytes_written > bytes_written)
//             {
//                 if (fd->bytes_written - bytes_written > MAX_MPIWRITE_SIZE)
//                 {
//                     to_write = MAX_MPIWRITE_SIZE;
//                 }
//                 else
//                 {
//                     to_write = fd->bytes_written - bytes_written;
//                 }
//             }
//         }
//     }

//     // index location calculation:
//     // for buffered, base_offset = 0, fd->offset = write loc
//     // for unbuffered, base_offset = write loc, fd->offset = 0
//     // for append buffered, base_offset = start, fd->offset = size
//     lseek (p->b.f, fd->base_offset + fd->offset, SEEK_SET);
//     write (p->b.f, buffer, buffer_size);
// }


void adios_xpmem_close (struct adios_file_struct * fd
                       ,struct adios_method_struct * method
                       )
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
            data_buffer = &p->data.buffer[p->data.d->offset];
            memcpy(data_buffer, fd->buffer, fd->bytes_written);

            //copy the index into the index area
            index_buffer = &p->index.buffer[p->index.d->offset];           
            memcpy(index_buffer, buffer, buffer_offset);

            log_debug("xpmem copied data into %p index into %p\n",
                      data_buffer, index_buffer);

            //set the sizes for the data
            p->data.d->size = fd->bytes_written;
            p->index.d->size = buffer_offset;

            log_debug("xpmem sizes = %d, %d\n",
                      fd->bytes_written, buffer_offset);
            
            //now set the version to 1
            p->data.d->version = 1;
            p->index.d->version = 1;
            
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

}

void adios_xpmem_finalize (int mype, struct adios_method_struct * method)
{
	adios_xpmem_data_struct * p = (adios_xpmem_data_struct *)
		method->method_data;
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
