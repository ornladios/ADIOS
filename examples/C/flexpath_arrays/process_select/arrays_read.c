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
#include "mpi.h"
#include "adios.h"
#include "public/adios_read.h"

int main (int argc, char ** argv) 
{
    /* application data structures */
    int         rank;
    int         NX, NY; 
    double      *t;
    int         *p;

    /* MPI and ADIOS data structures */
    MPI_Comm    comm = MPI_COMM_WORLD;

    /* MPI and ADIOS setup */
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    adios_read_init_method(ADIOS_READ_METHOD_FLEXPATH, comm, "");

    /* First read in the scalars to calculate the size of the arrays */

    /* get everything from single process - rank 0 for now*/
    ADIOS_SELECTION process_select;
    process_select.type=ADIOS_SELECTION_WRITEBLOCK;
    process_select.u.block.index = rank;

    /* read the size of arrays using local inq_var */

    /* Note: at this moment, timeout is not handled. It blocks until writer appears */
    ADIOS_FILE* afile = adios_read_open("arrays", ADIOS_READ_METHOD_FLEXPATH, comm, 
                                        ADIOS_LOCKMODE_NONE, 30.0);
    /* Read arrays for each time step */
    while(adios_errno != err_end_of_stream){

	ADIOS_VARINFO* nx_info = adios_inq_var( afile, "NX");
	if(nx_info->value) {
	    NX = *((int *)nx_info->value);
	}

	ADIOS_VARINFO* ny_info = adios_inq_var( afile, "NY");
	if(ny_info->value) {
	    NY = *((int *)ny_info->value);
	}
    
	/* Allocate space for the arrays */
	t = (double *) malloc (NX*NY*sizeof(double));    
	p = (int *) malloc (NX*sizeof(int));
	memset(t, 0, NX*NY*sizeof(double));
	memset(p, 0, NX*sizeof(int));
        /* schedule a read of the arrays */
        adios_schedule_read (afile, &process_select, "var_double_2Darray", 0, 1, t);
        adios_schedule_read (afile, &process_select, "var_int_1Darray", 0, 1, p);	
        /* commit request and retrieve data */
        adios_perform_reads (afile, 1);

        /* print result */
        printf("Results Rank=%d Step=%d p[] = [%d, %d,...] t[][] = [%.2f, %.2f]\n", 
                rank, afile->current_step, p[0], p[1], t[0], t[1]);
    
	/* block until next step is available (30 sec timeout unsupported) */
	adios_release_step(afile);
        adios_advance_step(afile, 0, 30);
	MPI_Barrier (comm);
	
	/* shutdown ADIOS and MPI */
    }
    adios_read_close(afile);	
    /* wait until all readers finish */
    adios_read_finalize_method(ADIOS_READ_METHOD_FLEXPATH);
    MPI_Finalize ();

    return 0;
}
