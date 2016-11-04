/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"

const int  NX = 10;
const int  NY = 5;
const char diagfilename[] = "diag.bp";
const char ckptfilename[] = "ckpt.bp";

const MPI_Comm comm = MPI_COMM_WORLD;
int rank, size;

void write_diag (int step, double * p)
{
    int64_t     adios_handle;
    if (rank==0) printf("Timestep %d write diagnostics\n", step);

    if (step==1)
        adios_open (&adios_handle, "diagnostics", diagfilename, "w", comm);
    else
        adios_open (&adios_handle, "diagnostics", diagfilename, "a", comm);

    adios_write (adios_handle, "NY", &NY);
    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "pressure", p);

    adios_close (adios_handle);
}

void write_checkpoint (int step, double * p, double *t)
{
    int64_t     adios_handle;
    if (rank==0) printf("Checkpointing at step %d\n", step);
    adios_open (&adios_handle, "checkpoint", ckptfilename, "w", comm);
    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "NY", &NY);
    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "step", &step);
    adios_write (adios_handle, "temperature", t);
    adios_write (adios_handle, "pressure", p);
    adios_close (adios_handle);
}

int main (int argc, char ** argv) 
{

    int         i, it;
    double      t[NX];
    double      p[NY];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_init ("time_aggregation.xml", comm);


    for (it = 1; it <= 100; it++)
    {

        for (i = 0; i < NX; i++)
            t[i] = it*1000.0 + rank*NX + i;

        for (i = 0; i < NY; i++)
            p[i] = it*1000.0 + rank*NY + i;

        write_diag(it, p);

        if ( it%30 == 0) {
            write_checkpoint(it, p, t);
        }

        MPI_Barrier (comm);
        if (rank==0) printf("Timestep %d written\n", it);
    }

    MPI_Barrier (comm);
    adios_finalize (rank);
    MPI_Finalize ();
    return 0;
}
