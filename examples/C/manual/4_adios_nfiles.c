/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 * ADIOS example from the User's manual
 *
 * Write N files from P processes by using the ADIOS MPI method and 
 * by grouping the processes into N groups.
 *
 */
#include <stdio.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size;
    int         NX = 10; 
    int         N = 3; /* number of files to write */
    double      t[NX];
    int         i;

    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;
    int         color, key;
    MPI_Comm    comm;
 
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    /* MPI_Comm_split partitions the world group into N disjoint subgroups, 
     * the processes are ranked in terms of the argument key. 
     * A new communicator comm is returned for this specific grid configuration
     */
    color = rank % N;
    key = rank / N;
    MPI_Comm_split (MPI_COMM_WORLD, color, key, &comm);

    for (i=0; i<NX; i++)
        t[i] = rank*NX + i;
            
    /* every P/N processes write into the same file 
     * there are N files generated. 
     */
    sprintf (filename, "restart_%5.5d.bp", color);
    adios_init ("config.xml", MPI_COMM_WORLD);
    adios_open (&adios_handle, "temperature", filename, "w", comm);
    #include "gwrite_temperature.ch"
    adios_close (adios_handle);
    adios_finalize (rank);
    MPI_Finalize ();
    return 0;
}

