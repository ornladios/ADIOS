/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include "mpi.h"
#include "public/adios.h"
int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank;
	int         NX = 10; 
	double      t[NX];

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	uint64_t    adios_groupsize, adios_totalsize;
	int64_t     adios_handle;
	MPI_Comm comm;

	int         color, key;
	int         size;
	int         i;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	/* MPI_Comm_split paritions the world group into disjointed 2 subgroups, 
	 * the processes are ranked in terms of the argument key  
	 *  a new communicator comm is returned for this specific grid configuration
	 */
	color = rank % 2;
	key = rank / 2;
	MPI_Comm_split (MPI_COMM_WORLD, color, key, &comm);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);
	for (i = 0; i < NX; i++)
		t [i] = 10*(key*2+color)+i;
	/* every P/2 processes write into the same file 
	 * there are 2 files generated. 
	 */
	sprintf (filename, "adios_global_%5.5d.bp", color);
	adios_init ("adios_global.xml", MPI_COMM_WORLD);
	adios_open (&adios_handle, "temperature", filename, "w", comm);
	#include "gwrite_temperature.ch"
	adios_close (adios_handle);
	adios_finalize (rank);
	MPI_Finalize ();
	return 0;
}
