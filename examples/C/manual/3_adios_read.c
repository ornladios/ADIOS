/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 * ADIOS example from the User's manual
 *
 * Read back data written by 2_adios_posix.
 * Use the same method and the same number of processes that were used
 * for writing the data with 2_adios_posix.
 *
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank;
    int         NX = 10;
    double      t[NX];
    char        result[1024], s[32];
    int         i;
    
    /* ADIOS variables declarations for matching gread_temperature.ch */
    uint64_t    adios_groupsize, adios_totalsize, adios_buf_size;
    int64_t     adios_handle;
    MPI_Comm    comm =  MPI_COMM_WORLD;
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    sprintf (filename, "restart.bp");
    adios_init ("config.xml", comm);
    adios_open (&adios_handle, "temperature", filename, "r", comm);
    #include "gread_temperature.ch"
    adios_close (adios_handle);
    adios_finalize (rank);
    MPI_Finalize ();

    sprintf(result, "rank=%d t=[%g", rank, t[0]);
    for (i=1; i<NX; i++) {
        sprintf (s, ",%g", t[i]);
        strcat (result, s);
    }
    printf("%s]\n", result);

    return 0;
}    

