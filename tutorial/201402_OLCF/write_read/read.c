/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*************************************************************/
/*          Example of reading arrays in ADIOS               */
/*    which were written from the same number of processors  */
/*                                                           */
/*        Similar example is manual/2_adios_read.c           */
/*************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
//#define _NOMPI

#include "adios.h"
#include "core/adios_logger.h"

#include "adios_read.h"

int main (int argc, char ** argv) 
{
	char        filename [256];
	double      *t;
	int         *p;
	enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_XPMEM;
	ADIOS_SELECTION * sel;
	int         rank, i, j;
	MPI_Comm    comm = MPI_COMM_WORLD;
	int x= 0;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);

	fprintf(stderr,      "in the read method\n");

	adios_verbose_level = 4;
	adios_logger_open("stderr", -1);

	adios_read_init_method (method, comm, "verbose=4");


	strcpy (filename, "arrays.bp");
	ADIOS_FILE * f = adios_read_open_file (filename, method, comm);

	int         NX = 10, NY = 100; 

	fprintf(stderr, "NX = %d\tNY=%d\n", NX, NY);

	if(0)
	{
		/* Specify a selection that points to a specific writer's block */
		sel = adios_selection_writeblock (rank);

		/* First get the scalars to calculate the size of the arrays */
		// adios_schedule_read (f, sel, "NX", 0, 1, &NX);
		// adios_schedule_read (f, sel, "NY", 0, 1, &NY);
		// adios_perform_reads (f, 1);

		printf("rank=%d: NX=%d NY=%d\n", rank, NX, NY);

		/* Allocate space for the arrays */
		t = (double *) malloc (NX*NY*sizeof(double));
		p = (int *) malloc (NX*sizeof(int));

		/* Read the arrays */
		adios_schedule_read (f, sel, "var_double_2Darray", 0, 1, t);
		adios_schedule_read (f, sel, "var_int_1Darray", 0, 1, p);
		adios_perform_reads (f, 1);
	}
	else
	{
		uint64_t    start[] = {0,0};
		uint64_t count[] = {2,2};

		/* Allocate space for the arrays */
		t = (double *) malloc (NX*NY*sizeof(double));
		p = (int *) malloc (NX*sizeof(int));

		sel = adios_selection_boundingbox(2, start, count);
		adios_schedule_read(f, sel, "var_double_2Darray", 0, 1, t);
		adios_perform_reads(f, 1);
	}
	
	/* At this point, we have the data in memory */
	// printf("rank=%d: p = [%d", rank, p[0]);
	// for (i = 1; i < NX; i++)
	// 	printf(", %d", p[i]);
	// printf("]\n");
    
	// printf("rank=%d: t[5,*] = [%6.2f", rank, t[5*NY]);
	// for (j = 1; j < NY; j++)
	// 	printf(", %6.2f", t[5*NY+j]);
	// printf("]\n");

	free (t);
	free (p);

	adios_read_close (f);
	MPI_Barrier (comm);
	adios_read_finalize_method (method);
	MPI_Finalize ();


	return 0;
}
