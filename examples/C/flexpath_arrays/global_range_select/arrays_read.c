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
    char        filename [256];
    int         rank, size, i, j;
    int         NX, NY; 
    double      *t;
    MPI_Comm    comm = MPI_COMM_WORLD;

    int         adios_err;
    int64_t     adios_handle, adios_buf_size;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);

    adios_read_init_method(ADIOS_READ_METHOD_FLEXPATH, comm, "");

    ADIOS_SELECTION global_range_select;
    //if(rank == 0){
    global_range_select.type=ADIOS_SELECTION_BOUNDINGBOX;
    global_range_select.u.bb.start = malloc(sizeof(int));
    global_range_select.u.bb.count = malloc(sizeof(int));
    (global_range_select.u.bb.start)[0] = 1;
    (global_range_select.u.bb.count)[0] = 1;
    (global_range_select.u.bb.start)[1] = 5;
    (global_range_select.u.bb.count)[1] = 5;
    global_range_select.u.bb.ndim = 2;
    //fprintf(stderr, "app got here\n");
    /* read the size of arrays using local inq_var */
    ADIOS_FILE* afile = adios_read_open_file("arrays", 
					     ADIOS_READ_METHOD_FLEXPATH, 
					     comm);
    
    int ii;
    while(adios_errno != err_end_of_stream){
      //for(ii = 0; ii<1000; ii++){   
	/* get a bounding box - rank 0 for now*/
	ADIOS_VARINFO* nx_info = adios_inq_var( afile, "NX");
	ADIOS_VARINFO* ny_info = adios_inq_var( afile, "NY");
	fprintf(stderr, "after inq var\n");
	if(nx_info->value) {
	    NX = *((int *)nx_info->value);
	}
	if(ny_info->value){
	    NY= *((int*)ny_info->value);
	}
    
	printf("\trank=%d: NX=%d\n", rank, NX);
	printf("\trank=%d: NY=%d\n", rank, NY);
    
	/* Allocate space for the arrays */
	t = (double *) malloc (5*1*sizeof(double));
	fprintf(stderr, "t %p\n", t);
      
	/* Read the arrays */
	fprintf(stderr, "example is calling schedule read\n");
	
	adios_schedule_read (afile, 
			     &global_range_select, 
			     "var_2d_array", 
			     0, 1, t);
        fprintf(stderr, "example is calling perform reads\n");
	adios_perform_reads (afile, 1);		
    
	//sleep(20);
    
	printf("rank=%d: t[0,5+x] = [%6.2f", rank, t[0]);
	for(j=1; j<5; j++) {
	    printf(", %6.2f", t[j]);
	}
	printf("]\n");
        adios_advance_step(afile, 0, 30);
	//MPI_Barrier (comm);
	//sleep(1);
    }
    //
    adios_read_close(afile);

    adios_read_finalize_method(ADIOS_READ_METHOD_FLEXPATH);

    MPI_Finalize ();

    return 0;
}
