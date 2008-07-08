#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

// dart
#include <sys/uio.h>
#if USE_PORTALS
#include "sender.h"
#endif

#include "binpack-general.h"
#include "binpack-utils.h"
#include "br-utils.h"
#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"

static int adios_dart_initialized = 0;

struct adios_DART_data_struct
{
    struct iovec ioarr;
};

#if USE_PORTALS
void adios_dart_init (const char * parameters
                     ,struct adios_method_struct * method
                     )
{
    if (!adios_dart_initialized)
    {
        adios_dart_initialized = 1;

        init (DEFAULT_CONF);
    }
    method->method_data = calloc (1, sizeof (struct adios_DART_data_struct));
}

void adios_dart_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
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

    d->ioarr.iov_base = buffer + start;
    d->ioarr.iov_len = end - start;

    ptl_cn_send_sync (&d->ioarr, count, dst_tcp);
//    ptl_cn_send (&d->ioarr, count, dst_tcp);
}

void adios_dart_get_write_buffer (struct adios_file_struct * fd
                                 ,struct adios_var_struct * v
                                 ,unsigned long long * size
                                 ,void ** buffer
                                 ,struct adios_method_struct * method
                                 )
{
}

void adios_dart_read (struct adios_file_struct * fd
                     ,struct adios_var_struct * v, void * buffer
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
#else
void adios_dart_init (const char * parameters
                     ,struct adios_method_struct * method
                     )
{
    fprintf (stderr, "adios_dart_init: DART disabled, no portals support\n");
    if (!adios_dart_initialized)
    {
        adios_dart_initialized = 1;
    }
    method->method_data = 0;
}

void adios_dart_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    fprintf (stderr, "adios_dart_open: DART disabled, no portals support\n");
}

void adios_dart_write (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * data
                      ,struct adios_method_struct * method
                      )
{
    fprintf (stderr, "adios_dart_write: DART disabled, no portals support\n");
}

void adios_dart_get_write_buffer (struct adios_file_struct * fd
                                 ,struct adios_var_struct * v
                                 ,unsigned long long * size
                                 ,void ** buffer
                                 ,struct adios_method_struct * method
                                 )
{
    fprintf (stderr, "adios_dart_get_write_buffer: DART disabled, "
                     "no portals support\n"
            );
}

void adios_dart_read (struct adios_file_struct * fd
                     ,struct adios_var_struct * v, void * buffer
                     ,struct adios_method_struct * method
                     )
{
    fprintf (stderr, "adios_dart_read: DART disabled, no portals support\n");
}

void adios_dart_close (struct adios_file_struct * fd
                      ,struct adios_method_struct * method
                      )
{
    fprintf (stderr, "adios_dart_close: DART disabled, no portals support\n");
}

void adios_dart_finalize (int mype, struct adios_method_struct * method)
{
    if (adios_dart_initialized)
    {
        adios_dart_initialized = 0;
    }
    fprintf (stderr, "adios_dart_finalize: "
                     "DART disabled, no portals support\n"
            );
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
