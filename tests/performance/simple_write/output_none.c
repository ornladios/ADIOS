/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* NO output functions for write.c */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include "output.h"
#include "timing.h"


MPI_Comm iocomm;

int output_init(MPI_Comm comm, int bufsizeMB)
{
    return 0;
}

int output_define(int nx, int ny, int gnx, int gny, int offsx, int offsy)
{
    return 0;
}

int output_dump(char *filename, int step, void *data)
{
    double t1, t2;
    t1 = MPI_Wtime();
    t2 = MPI_Wtime();
    Tio_open[step] = t2-t1;
    t1 = MPI_Wtime();
    Tio_group[step] = t1-t2;
    t2 = MPI_Wtime();
    Tio_write[step] = t2-t1;
    t1 = MPI_Wtime();
    Tio_close[step] = t1-t2;
    return 0;
}

int output_finalize (int rank)
{
    return 0;
}
