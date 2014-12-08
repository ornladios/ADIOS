/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>


#define _NOMPI

#include "adios.h"
#include "core/adios_logger.h"
#include "adios_read.h"

/*************************************************************/
/*          Example of writing arrays in ADIOS               */
/*                                                           */
/*        Similar example is manual/2_adios_write.c          */
/*************************************************************/
int main (int argc, char ** argv) 
{
	int         rank, i, j;
	MPI_Comm    comm = MPI_COMM_WORLD;
	int x= 0;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);

	// x = open("lock", O_CREAT|O_EXCL, 0666);
	// fprintf(stderr, "rank = %d\n", x);
	
		char        filename [256];
		int         NX = 10, NY = 100; 
		double      t[NX][NY];
		int         p[NX];

		uint64_t    adios_groupsize, adios_totalsize;
		int64_t     adios_handle;

		for (i = 0; i < NX; i++)
			for (j = 0; j< NY; j++)
				t[i][j] = rank * NX + i + j*(1.0/NY);

		for (i = 0; i < NX; i++)
			p[i] = rank * NX + i;

		strcpy (filename, "arrays.bp");
		adios_init ("arrays.xml", comm);
		adios_verbose_level = 4;
		adios_logger_open("stderr", rank);
		adios_open (&adios_handle, "arrays", filename, "w", comm);
#include "gwrite_arrays.ch"
		adios_close (adios_handle);

		MPI_Barrier (comm);

		sleep(1000);
		adios_finalize (rank);

		MPI_Finalize ();
	return 0;
}
