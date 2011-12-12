/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write a global array from N processors using MPI_AMR.
 *
 * How to run: mpirun -np <N> adios_amr_write
 * Output: adios_amr_write.bp
 * ADIOS config file: adios_amr_write.xml
 *
*/
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i;
    int         NX = 16; 
    double      t[NX], t1, t2;
    MPI_Comm    comm = MPI_COMM_WORLD;

    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    for (i = 0; i < NX; i++)
    {
        t[i] = rank*NX + i;
    }

    strcpy (filename, "adios_amr_write.bp");

    adios_init ("adios_amr_write.xml");

    adios_open (&adios_handle, "temperature", filename, "w", &comm);

    adios_groupsize = 4 \
                    + 4 \
                    + 4 \
                    + 8 * (1) * (NX);
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "temperature", t);

    adios_close (adios_handle);

    adios_finalize (rank);

    MPI_Finalize ();

    return 0;
}
