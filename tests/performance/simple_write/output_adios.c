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
#include <fcntl.h>
#include <errno.h>

#include "output.h"
#include "timing.h"

#include "mpi.h"
#include "adios.h"

static int rank;
static int64_t gh;
static uint64_t groupsize;
static MPI_Comm iocomm;
static int file_per_process = 0;  // 0: method dependent, 1: N process writes to N files
static int streaming = 0;  // 0: separate file names, 1: append to data.bp 

char * DEFAULT_ADIOSMETHOD_NAME   = "MPI";
char * DEFAULT_ADIOSMETHOD_PARAMS = "";

int output_init(MPI_Comm comm, int bufsizeMB)
{
    char * wmethodname, *wmethodparams;
    MPI_Comm_rank (comm, &rank);
    adios_init_noxml(comm);
    adios_declare_group(&gh,"writer","",adios_flag_yes);
    adios_set_max_buffer_size (bufsizeMB);

    // Select output method
    wmethodname = getenv("ADIOSMETHOD");
    wmethodparams = getenv("ADIOSMETHOD_PARAMS");
    if (!wmethodname)
        wmethodname = DEFAULT_ADIOSMETHOD_NAME;
    if (!wmethodparams)
        wmethodparams = DEFAULT_ADIOSMETHOD_PARAMS;
    adios_select_method (gh, wmethodname, wmethodparams, "");

    groupsize = 0;

    if (!strcmp(wmethodname, "POSIX1")) 
    {
        file_per_process = 1;
        iocomm = MPI_COMM_SELF;
    } 
    else 
    {
        iocomm = comm;
    }

    if (!strcmp(wmethodname, "DATASPACES") ||  
        !strcmp(wmethodname, "DIMES")      || 
        !strcmp(wmethodname, "FLEXPATH")    ) 
    {
        streaming = 1;
    } else {
        streaming = 0;
    }
    return 0;
}

int output_define(int nx, int ny, int gnx, int gny, int offsx, int offsy)
{
    char ldims[256], gdims[256], offs[256];
    snprintf (ldims, sizeof(ldims), "%d,%d", nx, ny);
    snprintf (gdims, sizeof(gdims), "%d,%d", gnx, gny);
    snprintf (offs,  sizeof(offs),  "%d,%d", offsx, offsy);
    adios_define_var (gh, "xy", "", adios_double, ldims, gdims, offs);
    groupsize = nx*ny*sizeof(double);
    return 0;
}

int output_dump(char *filename, int step, void *data)
{
    int64_t fh;
    uint64_t tsize;
    double t1, t2;
    char fname[256], mode[2]="w";

    if (streaming) 
    {
        snprintf (fname, sizeof(fname), "data.bp",filename);
        mode[0]='a'; mode[1] = 0;
    } 
    else if (file_per_process) 
    {
        snprintf (fname, sizeof(fname), "%s_%d.bp",filename, rank);
        mode[0]='w'; mode[1] = 0;
    } 
    else 
    {
        snprintf (fname, sizeof(fname), "%s.bp",filename);
        mode[0]='w'; mode[1] = 0;
    }

    t1 = MPI_Wtime();
    adios_open (&fh, "writer", fname, mode, iocomm);
    t2 = MPI_Wtime();
    Tio_open[step] = t2-t1;
    adios_group_size (fh, groupsize, &tsize);
    t1 = MPI_Wtime();
    Tio_group[step] = t1-t2;
    adios_write (fh, "xy", data);
    t2 = MPI_Wtime();
    Tio_write[step] = t2-t1;
    adios_close (fh);
    t1 = MPI_Wtime();
    Tio_close[step] = t1-t2;
    return 0;
}

int output_finalize (int rank)
{
    adios_finalize (rank);
    return 0;
}
