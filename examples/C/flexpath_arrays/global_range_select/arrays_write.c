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

/*************************************************************/
/*          Example of writing arrays in ADIOS               */
/*                                                           */
/*        Similar example is manual/2_adios_write.c          */
/*************************************************************/
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, j, offset, size_y;
    int         NX = 5; 
    int         NY = 2;
    int         NZ = 2;
    double      t[NX*NY*NZ];
    MPI_Comm    comm = MPI_COMM_WORLD;

    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
    
    strcpy (filename, "arrays");
    adios_init ("arrays.xml", comm);
    
    int total = NX * NY * NZ * size;
    int test_scalar = rank * 1000;

    offset = rank*NY;
    size_y = size*NY;

    int start = 0;
    int myslice = NX * NY * NZ;

    for (i = 0; i<5; i++) {       
	for (j=0; j<NY*NX*NZ; j++) {       
	    t[j] = rank*myslice + (start + j);
	}

	int s;
	for (s = 0; s<size; s++) {
	    if (s == rank) {
		fprintf(stderr, "rank %d: step: %d [", rank, i);
		int z;
		for(z=0; z<NX*NY*NZ;z++){
		    fprintf(stderr, "%lf, ", t[z]);
		}
		fprintf(stderr, "]\n");
	    }
	    fprintf(stderr, "\n");
	    fflush(stderr);
	    MPI_Barrier(MPI_COMM_WORLD);
	}

	start += total;
        //prints the array.
	adios_open (&adios_handle, "temperature", filename, "w", comm);
	
	adios_write (adios_handle, "/scalar/dim/NX", &NX);
	adios_write (adios_handle, "/scalar/dim/NY", &NY);
	adios_write (adios_handle, "/scalar/dim/NZ", &NZ);
	adios_write (adios_handle, "test_scalar", &test_scalar);
	adios_write (adios_handle, "size", &size);
	adios_write (adios_handle, "rank", &rank);
	adios_write (adios_handle, "offset", &offset);
	adios_write (adios_handle, "size_y", &size_y);
	adios_write (adios_handle, "var_2d_array", t);
	
	adios_close (adios_handle);
    }

    adios_finalize (rank);
    MPI_Finalize ();
    return 0;
}
