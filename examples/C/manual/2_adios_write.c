/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 * ADIOS example from the User's manual
 *
 * Write a separate file from each process by using the POSIX method or
 * write into a large single file from all processes using the MPI method.
 * You need to change only the method in the config.xml and rerun the 
 * program (no recompile is needed)
 *
 * In case of POSIX method, the output files will have the process rank 
 * appended to the file name (e.g. restart.bp.15).
 *
 * 4_posix_read.c example can read the output data from the same number 
 * of processors and using the same method. 
 * 
 * Application of the example: 
 *    Checkpoint/restart files 
 *
 * Note: bpls utility and the generic read API can see only the array
 * written from one of the processors. You need to use global arrays to
 * make the data available for the utilities or for reading data from 
 * arbitrary number of processors. 
 *
 */
#include <stdio.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank;
    int         NX = 10;
    double      t[NX];
    int         i;
    
    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;
    MPI_Comm    comm =  MPI_COMM_WORLD;
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    for (i=0; i<NX; i++)
        t[i] = rank*NX + i;

    sprintf (filename, "restart.bp");
    adios_init ("config.xml", comm);
    adios_open (&adios_handle, "temperature", filename, "w", comm);
    #include "gwrite_temperature.ch"
    adios_close (adios_handle);
    adios_finalize (rank);
    MPI_Finalize ();
    return 0;
}

