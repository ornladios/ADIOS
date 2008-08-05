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

static int adios_vtk_initialized = 0;

void adios_vtk_init (const char * parameters
                    ,struct adios_method_struct * method
                    )
{
    if (!adios_vtk_initialized)
    {
        adios_vtk_initialized = 1;
    }
    method->method_data = 0;
}

int adios_vtk_open (struct adios_file_struct * fd
                   ,struct adios_method_struct * method
                   )
{
    return 1;
}

enum ADIOS_FLAG adios_vtk_should_buffer (struct adios_file_struct * fd
                                        ,struct adios_method_struct * method
                                        ,void * comm
                                        )
{
    return adios_flag_yes;   // as far as we care, buffer
}

void adios_vtk_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,void * data
                     ,struct adios_method_struct * method
                     )
{
}

void adios_vtk_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
    *buffer = 0;
}

void adios_vtk_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v
                    ,void * buffer
                    ,struct adios_method_struct * method
                    )
{
}

void adios_vtk_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    if (fd->mode == adios_mode_read)
    {
        //adios_vtk_do_read (fd, method);
        struct adios_var_struct * v = fd->group->vars;
        while (v)
        {
            v->data = 0;
            v = v->next;
        }
    }
}

void adios_vtk_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
    if (adios_vtk_initialized)
        adios_vtk_initialized = 0;
}

void adios_vtk_end_iteration (struct adios_method_struct * method)
{
}

void adios_vtk_start_calculation (struct adios_method_struct * method)
{
}

void adios_vtk_stop_calculation (struct adios_method_struct * method)
{
}
