/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// xml parser
#include <mxml.h>

// dart
#include "dc_gspace.h"

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"

static int adios_dart_initialized = 0;

struct adios_DART_data_struct
{
    struct dcg_space *dcg;
    int version;            /* Versioning of objects put in space */
};

void adios_dart_init (const char * parameters
                     ,struct adios_method_struct * method
                     )
{
    struct adios_DART_data_struct * p = 0;
    if (!adios_dart_initialized)
    {
        adios_dart_initialized = 1;
    }
    method->method_data = calloc (1, sizeof (struct adios_DART_data_struct));
    p = (struct adios_DART_data_struct *) method->method_data;
    p->dcg = 0;
    p->version = 0;
}

int adios_dart_open (struct adios_file_struct * fd
                    ,struct adios_method_struct * method, void * comm
                    )
{
    struct adios_DART_data_struct * p = (struct adios_DART_data_struct *)
                                                              method->method_data;
    switch (fd->mode) 
    {
        case adios_mode_read:
        {
            fprintf (stderr, "ADIOS DART: Read mode is not supported\n");
        }

        case adios_mode_write:
        case adios_mode_append:
        {
        }
    } // switch (fd->mode)


    return 1;
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
    int dst_tcp = 2; // magic number
    int count = 1; // number of copies of buffer to send(?)
    struct adios_DART_data_struct * d =
                  (struct adios_DART_data_struct *) method->method_data;

}

void adios_dart_get_write_buffer (struct adios_file_struct * fd
                                 ,struct adios_var_struct * v
                                 ,uint64_t * size
                                 ,void ** buffer
                                 ,struct adios_method_struct * method
                                 )
{
}

void adios_dart_read (struct adios_file_struct * fd
                     ,struct adios_var_struct * v, void * buffer
                     ,uint64_t buffer_size
                     ,struct adios_method_struct * method
                     )
{
}

void adios_dart_close (struct adios_file_struct * fd
                      ,struct adios_method_struct * method
                      )
{
}

void adios_dart_finalize (int mype, struct adios_method_struct * method)
{
    if (adios_dart_initialized)
    {
        finish ();
        adios_dart_initialized = 0;
    }
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
