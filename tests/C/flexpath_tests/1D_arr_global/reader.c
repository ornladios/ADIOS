/**
 * reader.c
 *
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 *
 * Created on: Jul 1, 2013
 * Author: Magda Slawinska aka Magic Magg magg dot gatech at gmail.com
 *
 * This is a test for the FlexPath method based on examples/C/flexpath_arrays
 */

#include "mpi.h"
#include "adios.h"
#include "adios_read.h"

#include "misc.h"
#include "utils.h"
#include "test_common.h"
#include "cfg.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>


int main (int argc, char **argv){
	int rank =0, writer_size =0;
	int NX = 0;
	double *t = NULL;
	// this is an array we expect as a reference array
	double *t_ref = NULL;
	MPI_Comm comm = MPI_COMM_WORLD;
	diag_t diag = DIAG_OK;  // to store the diagnostic information
	struct test_info test_result = { TEST_PASSED, "1D_arr_global" };
	struct err_counts err = { 0, 0};
	struct adios_tsprt_opts adios_opts;

	GET_ENTRY_OPTIONS(adios_opts, "Runs readers. It is recommended to run as many readers as writers.");

	// adios read initialization
	MPI_Init( &argc, &argv);
	MPI_Comm_rank (comm, &rank);

        /*volatile int qur = 0;
        while(qur == 0) {  Wait for gdb  }
        MPI_Barrier(comm);*/

	// depending on the method
	SET_ERROR_IF_NOT_ZERO(adios_read_init_method(adios_opts.method, comm, adios_opts.adios_options), err.adios);
	RET_IF_ERROR(err.adios, rank);

	// I will be working with streams so the lock mode is necessary,
	// return immediately if the stream unavailable
	ADIOS_FILE *adios_handle = adios_read_open(FILE_NAME,adios_opts.method, comm, ADIOS_LOCKMODE_NONE, 0.0);
	if ( !adios_handle){
		p_error("Quitting ... (%d) %s\n", adios_errno, adios_errmsg());
		return DIAG_ERR;
	}

	// define portions of data how they will be read
	ADIOS_SELECTION *sel = NULL;
	ADIOS_VARINFO *avi = NULL;


	// read how many processors wrote that array
	avi = adios_inq_var (adios_handle, "size");
	if (!avi){
		p_error("rank %d: Quitting no inq var for writer size... (%d) %s\n", rank, adios_errno, adios_errmsg());
		CLOSE_ADIOS_READER(adios_handle, adios_opts.method);
		return DIAG_ERR;
	}
	writer_size = *((int*)avi->value);
	adios_free_varinfo(avi);
	avi = NULL;


	// read the size of the array per process
	avi = adios_inq_var (adios_handle, "NX");
	if (!avi){
		p_error("rank %d: Quitting ... (%d) %s\n", rank, adios_errno, adios_errmsg());
		CLOSE_ADIOS_READER(adios_handle, adios_opts.method)
		return DIAG_ERR;
	}

	// I expect a scalar that will tell me the size of an array
	assert(0 == avi->ndim);
	assert(adios_integer == avi->type);
	NX = *((int*)avi->value);
	// I don't need that variable any more
	adios_free_varinfo(avi);
	assert(NX_DIM == NX);
	avi = NULL;


	// this will define the slice that we want to read; each rank should
	// read its own slice written by a corresponding writer rank
	// the var_1d is a 2-dim ADIOS variable (because we decided to have
	// it as a global array); that's why we need a 2-dim ADIOS variable
	uint64_t count[2] = {1,0};
	count[1] = NX;
	uint64_t start[2] = {0,0};
	start[0] = rank % writer_size;

	sel = adios_selection_boundingbox(2,start, count);
	if( !sel ){
		p_error("rank %d: Quitting because no sel... (%d) %s\n", rank, adios_errno, adios_errmsg());
		diag = DIAG_ERR;
		goto close_adios;
	}

	// allocate the memory for the actual array to be read
	t = calloc(NX, sizeof(double));

	if (adios_schedule_read(adios_handle, sel, "var_1d_array",0,1,t) != 0){
		p_error("rank %d: Quitting schedule read for var array failed...(%d) %s\n", 
                                    rank, adios_errno, adios_errmsg());
		diag = DIAG_ERR;
		goto just_clean;
	}

	// not sure if this assumption is correct; difficult to find in the ADIOS sources
	if (adios_perform_reads(adios_handle, 1) != 0){
		p_error("rank %d: Couldn't perform the read...(%d) %s\n", rank, adios_errno, adios_errmsg());
		diag = DIAG_ERR;
		goto just_clean;
	}

	// make the reference array with reference values I expect to get
	t_ref = calloc(NX, sizeof(double));
	gen_1D_array(t_ref, NX, (rank % writer_size));

	// compare the values what you get with what you expect
	int i = 0;
        if(rank == 0)
        {
            printf("Reader received array: ");
        }
	for (i = 0; i < NX; ++i) {
                if(rank == 0)
                    printf("%f ", t[i]);
		if (t[i] != t_ref[i]) {
		    	test_failed(test_result.name, rank);
			test_result.result = TEST_FAILED;
			break;
		}
	}
        if(rank == 0)
            printf("\n");

	if (TEST_PASSED == test_result.result)
		test_passed(test_result.name, rank);
	else 
	    test_failed(test_result.name, rank);


	adios_release_step(adios_handle);
	// 0 - next available step, block for max 30 seconds until the next step
	// is available
	adios_advance_step(adios_handle, 0, 30);
	if (err_end_of_stream == adios_errno){
		printf("Rank %d: finished the step and moved on appropriately\n", rank);
	} else {
		printf("ERROR: adios_advance_step(); anyway Quitting ... Rank %d: (%d) %s\n", rank, adios_errno, adios_errmsg());
	}


just_clean:
	// clean everything
	adios_selection_delete(sel);
	sel = NULL;
	free(t);
	t = NULL;
	free(t_ref);
	t_ref = NULL;

close_adios:
	CLOSE_ADIOS_READER(adios_handle, adios_opts.method);

	if ((DIAG_OK == diag) && (TEST_PASSED == test_result.result)) {
		return 0;
	} else {
		return 1;
	}
}
