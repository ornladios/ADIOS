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

struct double_complex
{
   double r;
   double i;
};

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, it;
    int         NX = 10; 
    double      t[NX];
    struct double_complex c[NX];

    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;
    MPI_Comm    comm=MPI_COMM_WORLD;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (comm, &size);

    adios_init ("stat.xml", comm);
    strcpy (filename, "adios_stat.bp");

    for (it =0; it < 13; it++) {
        for (i = 0; i < NX; i++) {
            t[i] = it*100.0 + NX * rank + i;
            c[i].r = it * 10 + i + 8;
            c[i].i = it * 10 + i - 5;
        }
        // Introduce Inf value here    
        t[0] = 1 / 0.0;    

        if (it==0)
            adios_open (&adios_handle, "temperature", filename, "w", comm);
        else
            adios_open (&adios_handle, "temperature", filename, "a", comm);

#include "gwrite_stat.ch"
        adios_close (adios_handle);
        MPI_Barrier (comm);
    }

    printf ("[%d]: adios_stat.bp written successfully\n", rank);
    MPI_Barrier (comm);
    //if (rank==0) printf("Finalize adios\n");
    adios_finalize (rank);

    //if (rank==0) printf("Finalize MPI\n");
    MPI_Finalize ();
    return 0;
}
