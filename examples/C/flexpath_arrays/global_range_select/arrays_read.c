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


void
slice(uint64_t length, uint64_t *s, uint64_t *e, int rank, int mpisize)
{
    uint64_t start = 0;
    uint64_t end = 0;
    uint64_t rem = length % mpisize;

    start = length/mpisize * rank;
    end = length/mpisize * (rank+1);
    *s = start;
    *e = end;
    
    /* If our MPI size is greater
       than the number of y dimensions,
       then read the whole thing. */
    if (mpisize > length) {
	*e = length;
	*s = 0;
	return;\
    }
    if (end > length) {
        end = length;
        *e = end;
        return;
    }
    if (rank == mpisize-1) {
        end += rem;
        *e = end;
    }
}


int main (int argc, char ** argv) 
{
    int         rank, size, j;
    int         NX, NY, NZ; 
    double      *t;
    MPI_Comm    comm = MPI_COMM_WORLD;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method(ADIOS_READ_METHOD_FLEXPATH, comm, "");

    ADIOS_SELECTION *global_range_select;

    ADIOS_SELECTION scalar_block_select;
    scalar_block_select.type = ADIOS_SELECTION_WRITEBLOCK;
    scalar_block_select.u.block.index = 0;

    /* schedule_read of a scalar. */    
    int test_scalar = -1;
    ADIOS_FILE* afile = adios_read_open("arrays", 
					ADIOS_READ_METHOD_FLEXPATH, 
					comm,
					ADIOS_LOCKMODE_NONE, 0.0);

    /*int i;*/
    /* for(i=0; i<afile->nvars; i++){ */
    /* 	printf("var: %s\n", afile->var_namelist[i]); */
    /* } */
    
    int ii = 0;
    while(adios_errno != err_end_of_stream){       
        /* get a bounding box - rank 0 for now*/
        ADIOS_VARINFO *nx_info = adios_inq_var(afile, "/scalar/dim/NX");
        ADIOS_VARINFO *ny_info = adios_inq_var(afile, "/scalar/dim/NY");
	ADIOS_VARINFO *nz_info = adios_inq_var(afile, "/scalar/dim/NZ");

	ADIOS_VARINFO *arry = adios_inq_var( afile, "var_2d_array");

	/*
	ADIOS_VARINFO *size_info = adios_inq_var( afile, "size");
	int nx_val = *((int*)nx_info->value);
	int ny_val = *((int*)ny_info->value);
	int size_val = *((int*)size_info->value);
	printf("nx: %d, ny: %d, size: %d\n", nx_val, ny_val, size_val);
	*/
	
	// slice array along y dimension
	uint64_t my_ystart, my_yend;
	slice(arry->dims[1], &my_ystart, &my_yend, rank, size);

	/* printf("rank: %d my_ystart: %d, my_yend: %d\n", */
	/*        rank, (int)my_ystart, (int)my_yend); */

	uint64_t xcount = arry->dims[0];
	uint64_t ycount = my_yend - my_ystart;
	uint64_t zcount = arry->dims[2];

	uint64_t starts[] = {0, my_ystart, 0};
	uint64_t counts[] = {xcount, ycount, zcount};

	/* printf("rank: %d starts: %d %d %d. counts: %d %d %d\n", */
	/*        rank, */
	/*        (int)starts[0], (int)starts[1], (int)starts[2], */
	/*        (int)counts[0], (int)counts[1], (int)counts[2]); */

	global_range_select = adios_selection_boundingbox(arry->ndim, starts, counts);

	int nelem = xcount*ycount*zcount;

        if(nx_info->value) {
            NX = *((int *)nx_info->value);
        }
        if(ny_info->value){
            NY= *((int*)ny_info->value);
        }
        if(nz_info->value){
            NZ= *((int*)nz_info->value);
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
