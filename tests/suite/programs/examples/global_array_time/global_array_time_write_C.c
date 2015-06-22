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
int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank, size, i, it;
	int         NX = 10;
        // NY = 1 for testing purpose
	int         NY = 1; 
	double      t[NX];
	double      p[NY];

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	uint64_t    adios_groupsize, adios_totalsize;
	int64_t     adios_handle;
	MPI_Comm    comm=MPI_COMM_WORLD;
 
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (comm, &size);

	adios_init ("global_array_time_C.xml", comm);
    	strcpy (filename, "global_array_time_C.bp");
    	for (it =1; it <= 13; it++) {

        	for (i = 0; i < NX; i++)
            		t[i] = it*100.0 + rank*NX + i;

        	for (i = 0; i < NY; i++)
            		p[i] = it*1000.0 + rank*NY + i;
		
                if (it==1)
		    adios_open (&adios_handle, "restart", filename, "w", comm);
                else
		    adios_open (&adios_handle, "restart", filename, "a", comm);

        	#include "gwrite_restart.ch"
        	adios_close (adios_handle);
		MPI_Barrier (comm);
                //if (rank==0) printf("Timestep %d written\n", it+1);
 	}
	MPI_Barrier (comm);
        //if (rank==0) printf("Finalize adios\n");
    	adios_finalize (rank);

        //if (rank==0) printf("Finalize MPI\n");
    	MPI_Finalize ();
	return 0;
}
