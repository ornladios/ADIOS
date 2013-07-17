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

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>


// for printing the values of the variable
#define STR_BUFFER_SIZE 3000

#define CLOSE_ADIOS \
	do { 														\
		adios_read_close(adios_handler);  						\
		adios_read_finalize_method(method);						\
		MPI_Finalize(); \
	} while (0)

#define JUST_CLEAN \
	do {								\
		adios_selection_delete(sel);	\
		sel = NULL;						\
		free(t);						\
		t = NULL;						\
	} while (0)

#define CLEAN_ON_ERROR_AND_CLOSE_ADIOS 	\
	do {								\
		JUST_CLEAN;						\
		CLOSE_ADIOS;					\
	} while (0)


int main (int argc, char **argv){
	char filename[256];
	int rank =0, size =0;
	int NX = 0;
	double *t = NULL;
	// this is an array we expect as a reference array
	double *t_ref = NULL;
	MPI_Comm comm = MPI_COMM_WORLD;
	int diag = 0;  // to store the diagnostic information; 0 is OK, -1 error
	int test_passed = TEST_PASSED;   // TEST_PASSED if test passed; TEST_FAILED if test failed
	enum ADIOS_READ_METHOD method = METHOD;


	if (1 < argc){
		usage(argv[0], "Runs readers. It is recommended to run as many readers as writers.");
		return 0;
	}
	// adios read initialization
	MPI_Init( &argc, &argv);
	MPI_Comm_rank (comm, &rank);

	// depending on the method
	if( (diag = adios_read_init_method( method, comm, ADIOS_OPTIONS)) != 0){
		printf("ERROR: Quitting ... Issues with initialization adios_read: (%d) %s\n", adios_errno, adios_errmsg());
		return -1;
	}

	// I will be working with streams so the lock mode is necessary,
	// return immediately if the stream unavailable
	ADIOS_FILE *adios_handler = adios_read_open(FILE_NAME,method, comm, ADIOS_LOCKMODE_NONE, 0.0);
	if ( !adios_handler){
		printf("ERROR: Quitting ... (%d) %s\n", adios_errno, adios_errmsg());
		return -1;
	}

	// define portions of data how they will be read
	ADIOS_SELECTION *sel = NULL;
	ADIOS_VARINFO *avi = NULL;


	int64_t adios_handle, adios_buf_size;
	// for storing the variables
	char buf[STR_BUFFER_SIZE];

	int step = 0;
	int n = 0;

	// read how many processors wrote that array
	avi = adios_inq_var (adios_handler, "size");
	if (!avi){
		printf("ERROR: rank %d: Quitting ... (%d) %s\n", rank, adios_errno, adios_errmsg());
		CLOSE_ADIOS;
		return -1;
	}
	size = *((int*)avi->value);
	adios_free_varinfo(avi);
	avi = NULL;

	// if I run the more readers than writers; just release
	// the excessive readers
	if (rank >= size){
		printf("INFO: rank %d: I am an excessive rank. Nothing to read ...\n", rank);
		CLOSE_ADIOS;
		return 0;
	}

	// read the size of the array
	avi = adios_inq_var (adios_handler, "NX");
	if (!avi){
		printf("ERROR: rank %d: Quitting ... (%d) %s\n", rank, adios_errno, adios_errmsg());
		CLOSE_ADIOS;
		return -1;
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
	start[0] = rank;

	sel = adios_selection_boundingbox(2,start, count);
	if( !sel ){
		printf("ERROR: rank %d: Quitting ... (%d) %s\n", rank, adios_errno, adios_errmsg());
		CLOSE_ADIOS;
		return -1;
	}

	// allocate the memory for the actual array to be read
	t = calloc(NX, sizeof(double));

	if (adios_schedule_read(adios_handler, sel, "var_1d_array",0,1,t) != 0){
		printf("ERROR: rank %d: Quitting ...(%d) %s\n", rank, adios_errno, adios_errmsg());
		CLEAN_ON_ERROR_AND_CLOSE_ADIOS;
		return -1;
	}

	// not sure if this assumption is correct; difficult to find in the ADIOS sources
	if (adios_perform_reads(adios_handler, 1) != 0){
		printf("ERROR: rank %d: Quitting ...(%d) %s\n", rank, adios_errno, adios_errmsg());
		CLEAN_ON_ERROR_AND_CLOSE_ADIOS;
		return -1;
	}

	// make the reference array with reference values I expect to get
	t_ref = calloc(NX, sizeof(double));
	gen_1D_array(t_ref, NX, rank);

	// compare the values what you get with what you expect
	int i = 0;
	for (i = 0; i < NX; ++i) {
		if (t[i] != t_ref[i]) {
			printf(
					"\nRank %d: for t[%d] TEST FAILED (expected %.1f, got %.1f)\n",
					rank, i, t_ref[i], t[i]);
			test_passed = TEST_FAILED;
			break;
		}
	}

	if (TEST_PASSED == test_passed)
		printf("\nRank %d: TEST PASSED\n", rank);

	// clean everything
	JUST_CLEAN;

/*#ifdef FLEXPATH_METHOD
	adios_release_step(adios_handler);
	// 0 - next available step, block for max 30 seconds until the next step
	// is available
	adios_advance_step(adios_handler, 0, 30);
	if (0 == adios_errno){
		printf("Rank %d: proceeding to the next step ...\n", rank);
	} else {
		printf("ERROR: adios_advance_step(); anyway Quitting ... Rank %d: (%d) %s\n", rank, adios_errno, adios_errmsg());
	}
#endif
*/
	CLOSE_ADIOS;

	return test_passed;
}
