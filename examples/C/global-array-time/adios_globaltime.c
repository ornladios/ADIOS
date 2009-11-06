#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank, size, i, it;
	int         NX = 10; 
	double      t[NX];

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	int         adios_err;
	uint64_t    adios_groupsize, adios_totalsize;
	int64_t     adios_handle;
	MPI_Comm    comm=MPI_COMM_WORLD;
 
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (comm, &size);

	adios_init ("adios_globaltime.xml");
    	strcpy (filename, "adios_globaltime.bp");
    	for (it =0; it < 13; it++) {

        	for (i = 0; i < NX; i++)
            		t[i] = it*100.0 + rank*NX + i;
		
                if (it==0)
		    adios_open (&adios_handle, "temperature", filename, "w");
                else
		    adios_open (&adios_handle, "temperature", filename, "a");

        	#include "gwrite_temperature.ch"
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
