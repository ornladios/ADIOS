/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mpi.h"
#include "adios.h"

// t1 is split time, t2 is write time, t3 is close time
double t1=0, t2=0, t3=0;
#define START(t) t -= MPI_Wtime()
#define STOP(t) t += MPI_Wtime()

#define print(...) fprintf (stderr, __VA_ARGS__); 
#define print0(...) if (!rank) fprintf (stderr, __VA_ARGS__); 

int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank, size, i, it;
	int         NX = 10;
    // NY = 1 for testing purpose
	int         NY = 1; 
	double      *t;
	double      *p;

    
	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	uint64_t    adios_groupsize, adios_totalsize;
	int64_t     adios_handle;
	MPI_Comm    comm=MPI_COMM_WORLD;

    if(argc > 1)
    {
        NX = atoi(argv[1]);
        NY = NX;
    }
    else
    {
        NX=10;
        NY=10;
    }

    t = (double*)malloc(sizeof(double)* NX);
    p = (double*)malloc(sizeof(double)* NY);

    
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (comm, &size);

	adios_init ("adios_globaltime.xml", comm);
    strcpy (filename, "adios_globaltime.bp");
    for (it =1; it <= 1; it++) {

        for (i = 0; i < NX; i++)
            t[i] = it*100.0 + rank*NX + i;

        for (i = 0; i < NY; i++)
            p[i] = it*1000.0 + rank*NY + i;
		
        if (it==1)
		    adios_open (&adios_handle, "restart", filename, "w", MPI_COMM_SELF);
        else
		    adios_open (&adios_handle, "restart", filename, "a", MPI_COMM_SELF);

        START(t1);
#include "gwrite_restart.ch"
        STOP(t1);
        START(t2);
        adios_close (adios_handle);
        STOP(t2);
		MPI_Barrier (comm);
        //if (rank==0) printf("Timestep %d written\n", it+1);
 	}
	MPI_Barrier (comm);
    //if (rank==0) printf("Finalize adios\n");
    print ("RESULTS: %i: split time: %f close time: %f\n", rank, t1, t2);

    adios_finalize (rank);

    //if (rank==0) printf("Finalize MPI\n");
    MPI_Finalize ();
	return 0;
}
