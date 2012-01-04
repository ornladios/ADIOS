/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write a global array from N processors with gwrite
 *
 * How to run: mpirun -np <N> adios_global
 * Output: adios_global.bp
 * ADIOS config file: adios_global.xml
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
		t[i] = rank*NX + i;

	strcpy (filename, "adios_global.bp");


	adios_init ("adios_global.xml");
//    MPI_Barrier(MPI_COMM_WORLD);
//    t1=MPI_Wtime();

	adios_open (&adios_handle, "temperature", filename, "w", &comm);
	#include "gwrite_temperature.ch"
	adios_close (adios_handle);
    MPI_Barrier(MPI_COMM_WORLD);
    t2=MPI_Wtime();
/*
    if (rank == 0)
        printf ("%f\n", t2 - t1);
*/

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}
