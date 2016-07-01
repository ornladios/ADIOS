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
	int         rank, size, i, j;
	int         NX = 8, NY = 16, GX, GY, OX, OY;
	double      t[NX*NY];
	MPI_Comm    comm = MPI_COMM_WORLD;

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	uint64_t    adios_groupsize, adios_totalsize;
	int64_t     adios_handle;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

    GX = NX;
    GY = size * NY;
    OX = 0;
    OY = rank * NY;

	for (i = 0; i < NX; i++)
        for (j = 0; j < NY; j++)
	        t[i * NY + j] = 1.0;
	        //t[i * NY + j] = 10 * rank + i * NY + j;
            //
    t[18] = 10.0;
    t[19] = 10.0;
 //   t[20] = 10.0;
    t[26] = 10.0;
    t[27] = 10.0;
//    t[28] = 10.0;
//    t[34] = 10.0;
//    t[35] = 10.0;
//    t[36] = 10.0;


//    t[18 + 64] = 10.0;
//    t[19 + 64] = 10.0;
//    t[20 + 64] = 10.0;
//    t[26 + 64] = 10.0;
//    t[27 + 64] = 10.0;
//    t[28 + 64] = 10.0;
//    t[34 + 64] = 10.0;
//    t[35 + 64] = 10.0;
    t[36 + 64] = 10.0;

	strcpy (filename, "adios_global.bp");

	adios_init ("adios_global.xml", comm);

	adios_open (&adios_handle, "temperature", filename, "w", comm);
    adios_groupsize = 6 * 4 + 8 * NX * NY; 
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "NY", &NY);
    adios_write (adios_handle, "GX", &GX);
    adios_write (adios_handle, "GY", &GY);
    adios_write (adios_handle, "OX", &OX);
    adios_write (adios_handle, "OY", &OY);
    adios_write (adios_handle, "temperature", t);

	adios_close (adios_handle);

    MPI_Barrier (comm);

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}
