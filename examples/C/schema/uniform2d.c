/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write variables along with an unstructured mesh. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include "mpi.h"
#include "adios.h"

//will work with 12 cores, which are arranged by npx=4, npy=3 (4x3)
char   npx_str[256];       // # of procs in x dim (string value)
char   npy_str[256];       // # of procs in y dim (string value)
int    npx;                // # of procs in x direction
int    npy;                // # of procs in y direction
int    nproc;               // # of total procs

void printUsage(char *prgname)
{
    printf("Usage: mpirun -np <N> %s <nx> <ny>\n"
           "    <nx> <ny>  2D decomposition values in each dimension of an 2D array\n"
           "         The product of these number must be equal the number of processes\n"
           "         e.g. for N=12 you may use  4 3\n"
        ,prgname);
}


int processArgs(int argc, char ** argv)
{
    if (argc < 3) {
        printUsage (argv[0]);
        return 1;
    }

    strncpy(npx_str, argv[1], sizeof(npx_str));
    strncpy(npy_str, argv[2], sizeof(npy_str));  

    npx = atoi(npx_str);
    npy = atoi(npy_str);

    if (npx*npy != nproc) {
        printf ("ERROR: Product of decomposition numbers in X and Y dimension %d != number of processes %d\n", npx*npy, nproc);
        printUsage(argv[0]);
        return 1;
    }

    return 0;
}



int main (int argc, char ** argv ) 
{
    MPI_Comm    comm = MPI_COMM_WORLD;
    int         rank;
    int         ndx, ndy;             // size of array per processor
    double * data;

    int         O1 = 0;               //origin in x direction
    int         O2 = 0;               //origin in y direction
    int         S1 = 1;               //spacing in x direction
    int         S2 = 2;               //spacing in y direction

    // Offsets and sizes
    int         offs_x, offs_y;       //offset in x and y direction
    int         nx_local, ny_local;   //local address
    int         nx_global, ny_global; //global address
    int         posx, posy;           // position index in the array
    int         i,j;
  
    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &nproc);

    if (processArgs(argc, argv)) {
        return 1;
    }

    //will work with each core writing ndx = 65, ndy = 129, (65*4,129*3) global
    ndx = 65;
    ndy = 129;

    //2D array with block,block decomposition
    posx = rank%npx;           // 1st dim
    posy = rank/npx;           // 2nd dim
    offs_x = posx * ndx;
    offs_y = posy * ndy;
    nx_local = ndx;
    ny_local = ndy;
    nx_global = npx * ndx;
    ny_global = npy * ndy;

    data = malloc (ndx * ndy * sizeof(double));
    for( i = 0; i < ndx; i++ )
        for( j = 0; j < ndy; j++)
            data[i*ndy + j] = 1.0*rank;
 
    adios_init ("uniform2d.xml", comm);
    adios_open (&adios_handle, "uniform2d", "uniform2d.bp", "w", comm);

    adios_groupsize = 7*sizeof(int) + 4*sizeof(double)\
	+ sizeof(double) * (nx_local*ny_local) ;

    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "nproc", &nproc);
    adios_write (adios_handle, "nx_global", &nx_global);
    adios_write (adios_handle, "ny_global", &ny_global);
    adios_write (adios_handle, "offs_x", &offs_x);
    adios_write (adios_handle, "offs_y", &offs_y);
    adios_write (adios_handle, "nx_local", &nx_local);
    adios_write (adios_handle, "ny_local", &ny_local);
    adios_write (adios_handle, "O1", &O1);
    adios_write (adios_handle, "O2", &O2);
    adios_write (adios_handle, "S1", &S1);
    adios_write (adios_handle, "S2", &S2);
    adios_write (adios_handle, "data", data);

    adios_close (adios_handle);

    MPI_Barrier (comm);

    free (data);
    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}
