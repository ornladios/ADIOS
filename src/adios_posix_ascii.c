#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// xml parser
#include <mxml.h>

#include "binpack-general.h"
#include "binpack-utils.h"
#include "br-utils.h"
#include "bw-utils.h"
#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

static int adios_posix_ascii_initialized = 0;

void adios_posix_ascii_init (const char * parameters
                            ,struct adios_method_struct * method
                            )
{
    if (!adios_posix_ascii_initialized)
    {
        adios_posix_ascii_initialized = 1;
    }
    method->method_data = 0;
}

int adios_posix_ascii_open (struct adios_file_struct * fd
                           ,struct adios_method_struct * method
                           )
{
    return 0;
}

enum ADIOS_FLAG adios_posix_ascii_should_buffer (struct adios_file_struct * fd
                                          ,struct adios_method_struct * method
                                          ,void * comm
                                          )
{
    return adios_flag_yes;   // as far as we care, buffer
}

void adios_posix_ascii_write (struct adios_file_struct * fd
                             ,struct adios_var_struct * v
                             ,void * data
                             ,struct adios_method_struct * method
                             )
{
}

void adios_posix_ascii_get_write_buffer (struct adios_file_struct * fd
                                        ,struct adios_var_struct * v
                                        ,uint64_t * size
                                        ,void ** buffer
                                        ,struct adios_method_struct * method
                                        )
{
    *buffer = 0;
}

void adios_posix_ascii_read (struct adios_file_struct * fd
                            ,struct adios_var_struct * v
                            ,void * buffer
                            ,struct adios_method_struct * method
                            )
{
    v->data = buffer;
}

void adios_posix_ascii_close (struct adios_file_struct * fd
                             ,struct adios_method_struct * method
                             )
{
}

void adios_posix_ascii_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
    if (adios_posix_ascii_initialized)
        adios_posix_ascii_initialized = 0;
}

void adios_posix_ascii_end_iteration (struct adios_method_struct * method)
{
}

void adios_posix_ascii_start_calculation (struct adios_method_struct * method)
{
}

void adios_posix_ascii_stop_calculation (struct adios_method_struct * method)
{
}
