/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS output functions for write.c */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include "output.h"
#include "timing.h"

#include "mpi.h"

static int rank;
static int fh;
static int nbytes = 0;
static char filename[256];
static int last_step;

int output_init(MPI_Comm comm, int bufsizeMB)
{
    MPI_Comm_rank (comm, &rank);
    snprintf (filename, sizeof(filename), "data%6.6d.posix",rank);
    Tio_open[0] = MPI_Wtime();
    fh = creat (filename, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH); 
    if (fh == -1) {
        fprintf (stderr, "Error creating file '%s' in output_posix, errno=%d\n",
                filename, errno);
    }
    Tio_open[0] = MPI_Wtime() - Tio_open[0];
    return errno;
}

int output_define(int nx, int ny, int gnx, int gny, int offsx, int offsy)
{
    nbytes = sizeof(double) * nx * (ssize_t) ny;
    return 0;
}

int output_dump(char *filename, int step, void *data)
{
    ssize_t nbytes_written;

    if (step>0)
        Tio_open[step] = 0.0;
    Tio_group[step] = 0.0;
    Tio_close[step] = 0.0;
    Tio_write[step] = MPI_Wtime();

    nbytes_written = write (fh, data, nbytes);
    if (nbytes_written == -1) {
        fprintf (stderr, "Error writing data to file '%s' in output_posix, errno=%d\n",
                filename, errno);
    } else if (nbytes_written != nbytes) {
        fprintf (stderr, "Error writing data to file '%s' in output_posix. "
                 "Expected to write %zu bytes but wrote only %zu bytes\n",
                filename, nbytes, nbytes_written);
    }

    Tio_write[step] = MPI_Wtime() - Tio_write[step];
    last_step = step;
    return 0;
}

int output_finalize (int rank)
{
    Tio_close[last_step] = MPI_Wtime();
    if (close(fh) == -1) {
        fprintf (stderr, "Error closing file '%s' in output_posix, errno=%d\n",
                filename, errno);
    }
    Tio_close[last_step] = MPI_Wtime() - Tio_close[last_step];
    return 0;
}
