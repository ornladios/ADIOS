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
    int         rank, j;
    int         NX, NY; 
    double      *t;
    MPI_Comm    comm = MPI_COMM_WORLD;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);

    adios_read_init_method(ADIOS_READ_METHOD_FLEXPATH, comm, "");

    ADIOS_SELECTION *global_range_select;

    ADIOS_SELECTION scalar_block_select;
    scalar_block_select.type = ADIOS_SELECTION_WRITEBLOCK;
    scalar_block_select.u.block.index = rank;

    /* schedule_read of a scalar. */    
    int test_scalar = -1;
    ADIOS_FILE* afile = adios_read_open("arrays", 
					ADIOS_READ_METHOD_FLEXPATH, 
					comm,
					ADIOS_LOCKMODE_NONE, 0.0);

    /* for(i=0; i<afile->nvars; i++){ */
    /* 	printf("var: %s\n", afile->var_namelist[i]); */
    /* } */
    
    int ii = 0;
    while(adios_errno != err_end_of_stream){       
        /* get a bounding box - rank 0 for now*/
        ADIOS_VARINFO *nx_info = adios_inq_var( afile, "scalar/dim/NX");
        ADIOS_VARINFO *ny_info = adios_inq_var( afile, "scalar/dim/NY");
	ADIOS_VARINFO *size_info = adios_inq_var( afile, "size");
	ADIOS_VARINFO *arry = adios_inq_var( afile, "var_2d_array");

	int nx_val = *((int*)nx_info->value);
	int ny_val = *((int*)ny_info->value);
	int size_val = *((int*)size_info->value);

	printf("nx: %d, ny: %d, size: %d\n", nx_val, ny_val, size_val);
	
	uint64_t xcount = arry->dims[0];
	uint64_t ycount = arry->dims[1];

	uint64_t starts[] = {0,0};
	uint64_t counts[] = {xcount, ycount};

	global_range_select = adios_selection_boundingbox(2, starts, counts);

	int nelem = xcount*ycount;

        if(nx_info->value) {
            NX = *((int *)nx_info->value);
        }
        if(ny_info->value){
            NY= *((int*)ny_info->value);
        }
    
	if(rank == 0){
	    int n;
	    printf("dims: [ ");
	    for(n=0; n<arry->ndim; n++){
		printf("%d ", (int)arry->dims[n]);
	    }
	    printf("]\n");
	}
    
        /* Allocate space for the arrays */
        int arr_size = sizeof(double) * nelem;
        t = (double *) malloc (arr_size);
        memset(t, 0, arr_size);
        //fprintf(stderr, "t %p\n", t);
      
        /* Read the arrays */        
        adios_schedule_read (afile, 
                             global_range_select, 
                             "var_2d_array", 
                             0, 1, t);
	adios_schedule_read (afile,
			     &scalar_block_select,
			     "test_scalar",
			     0, 1, &test_scalar);

        adios_perform_reads (afile, 1);                
    
        //sleep(20);
    
        printf("Rank=%d: test_scalar: %d step: %d, t[0,5+x] = [", rank, test_scalar, ii);
        for(j=0; j<nelem; j++) {
            printf(", %6.2f", t[j]);
        }
        printf("]\n\n");
        adios_release_step(afile);
        adios_advance_step(afile, 0, 30);
        ii++;
        //MPI_Barrier (comm);
        //sleep(1);
    }
    //
    adios_read_close(afile);

    adios_read_finalize_method(ADIOS_READ_METHOD_FLEXPATH);

    MPI_Finalize ();

    return 0;
}
