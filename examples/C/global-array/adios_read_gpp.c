/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read a global array from N processors with gread
 *
 * How to run: mpirun -np <N> adios_global
 * Reads output from adios_global.c
 * ADIOS config file: adios_global.xml
 *
*/
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "public/adios_read.h"
int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank, size, i;
	int         NX = 4096; 
	double      t[NX];
	MPI_Comm    comm = MPI_COMM_WORLD;

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	ADIOS_FILE *fp;
	ADIOS_SELECTION *s;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

	strcpy (filename, "adios_global.bp");

    //adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "");

    fp = adios_read_open_file ("adios_global.bp", ADIOS_READ_METHOD_BP, comm);

    #include "gread_temperature.ch"

    adios_read_close (fp);

    // Verify data
	for (i = 0; i < NX; i++)
		if (t[i] != rank*NX + i)
        {
            fprintf (stderr, "Error detected\n");
        }

    
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);

	MPI_Finalize ();
	return 0;
}
