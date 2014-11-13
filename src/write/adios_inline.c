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

#include "public/adios.h" // MPI or dummy MPI
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"

static int adios_insitu_initialized = 0;

typedef struct adios_insitu_buffer_
{

}adios_insitu_buffer;

typedef struct adios_insitu_data_struct
{
    adios_insitu_buffer b;

}insitu_data;

void adios_insitu_init (const PairStruct * parameters, struct adios_method_struct * method)
{
    insitu_data *mdata = NULL;
    if(!adios_insitu_initialized)
    {
	adios_insitu_initialized = 1;
    }
    else
	return;

    mdata = (insitu_data*) malloc(sizeof(insitu_data));
    method->method_data = (void*) mdata;

    //initialize the buffer we are going to write into
    adios_buffer_struct_init(&p->b);

    //we need to create the memory that we use here?

}

int adios_insitu_open (struct adios_file_struct * fd,
		       struct adios_method_struct * method, void * comm)
{

    /*
     * We don't need to open any file here.
     * But we need to create a block of memory that we will write into
     * And we need to create a structure to identify the block of memory
     * A hash table(?) to translate the file name into 
     * the buffer pointer
     */

    insitu_data *mdata = (insitu_data*) method->method_data;
    
    switch(fd->mode)
    {
    case adios_mode_read:
	break;
    case adios_mode_write:
	break;
    case adios_mode_append:
	break;

    default:
	break;
    }

    return -1;
}

enum ADIOS_FLAG adios_insitu_should_buffer (struct adios_file_struct * fd,
					    struct adios_method_struct * method)
{
    return adios_flag_no;
}

void adios_insitu_write (struct adios_file_struct * fd,
			 struct adios_var_struct * v,
			 void * data,
			 struct adios_method_struct * method)
{

}

void adios_insitu_get_write_buffer (struct adios_file_struct * fd,
				    struct adios_var_struct * v,
				    uint64_t * size,
				    void ** buffer,
				    struct adios_method_struct * method)
{
    *buffer = 0;
    *size = 0;
}

void adios_insitu_read (struct adios_file_struct * fd,
			struct adios_var_struct * v,
			void * buffer,
			uint64_t buffer_size,
			struct adios_method_struct * method)
{

}

static void adios_insitu_do_write (struct adios_file_struct * fd,
				   struct adios_method_struct * method,
				   char * buffer,
				   uint64_t buffer_size)
{

}

static void adios_insitu_do_read (struct adios_file_struct * fd,
				  struct adios_method_struct * method)
{

}

void adios_insitu_close (struct adios_file_struct * fd,
			 struct adios_method_struct * method)
{

}

void adios_insitu_finalize (int mype, struct adios_method_struct * method)
{

}

void adios_insitu_end_iteration (struct adios_method_struct * method)
{
}

void adios_insitu_start_calculation (struct adios_method_struct * method)
{
}

void adios_insitu_stop_calculation (struct adios_method_struct * method)
{
}
