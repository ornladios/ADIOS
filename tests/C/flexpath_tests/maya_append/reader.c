/**
 * reader.c
 *
 *  Created on: Aug 21, 2013
 *  Author: Magda Slawinska aka Magic Magg magg dot gatech at gmail.com
 */

#include <string.h>
#include "mpi.h"
#include "adios.h"
#include "adios_read.h"

#include "misc.h"
#include "test_common.h"
#include "cfg.h"

/**
 * For storing errors
 */
struct err_counts {
	int adios;        // counter for adios calls errors
	int test;		  // counter for comparisons errors
};


/**
 * breaks the loop if error count is positive
 *
 * @param err_count The variable that positive value indicates that ther is an error
 */
#define BREAK_IF_ERROR(err_count) \
	if ( err_count > 0) { \
		break; \
	}


/**
 * wrapper for scheduling adios reads; this macro assumes existence of
 * quite a few important variables; please take a look and be careful
 * how to use it
 *
 * @param path_str The path to the variable
 * @param out_buf  The output buffer
 */
#define READ_FULLPATH(attribute, grid_func_name, out_buf) \
	sprintf(fullpath, "%s%s", attribute, grid_func_name);  \
	TEST_ADIOS_ERROR__IF_NOT_ZERO(adios_schedule_read(adios_handle, sel, fullpath,0,10, out_buf), error_counts.adios);


/**
 * checks if the adios call returned not zero; requires declaration of error_count
 * int adios_err_count++
 * @param fn_call adios function call
 * @param (in/out) err_count incremented if the error observed
 */
#define TEST_ADIOS_ERROR__IF_NOT_ZERO(fn_call, err_count)             \
  do {                                                                  \
	 int _error_code = fn_call;                                         \
     if (_error_code != 0){                                             \
       p_error("rank %d: %s: (%d) %s\n", rank, #fn_call, adios_errno, adios_errmsg()) ;\
       err_count++;                                                   \
     }                                                                  \
  } while (0)

/**
 * assumes that err_count is defined
 * @param value_ref The reference value
 * @param value     The actual value
 * @param err_count The value of the error counter; if the error increases
 *                  the value will be increased
 */
#define TEST_INT_EQUAL(value_ref, value, err_count) \
	if (value != value_ref ){ \
		test_passed = TEST_FAILED; \
		p_test_failed("(expected=%d, got=%d)\n", value_ref, value); \
		err_count++; \
	}

/**
 * assumes that err_count is defined
 * @param value_ref The reference value
 * @param value     The actual value
 * @param err_count The value of the error counter; if the error increases
 *                  the value will be increased
 */
#define TEST_DOUBLE_EQUAL(value_ref, value, err_count) \
	if (value != value_ref ){ \
		test_passed = TEST_FAILED; \
		p_test_failed("(expected=%f, got=%f)\n", value_ref, value); \
		err_count++; \
	}


// for printing the values of the variable
#define STR_BUFFER_SIZE 100

int main (int argc, char **argv){
	char filename[256];
	int rank =0, size =0;
	MPI_Comm comm = MPI_COMM_WORLD;
	int test_passed = TEST_PASSED;   // TEST_PASSED if test passed; TEST_FAILED if test failed
	enum ADIOS_READ_METHOD method = METHOD;
	const char * TEST_NAME="maya_append";
	struct err_counts error_counts = {0, 0};

	if (1 < argc){
		usage(argv[0], "Runs readers as many as you want.");
		return 0;
	}

	// adios read initialization
	MPI_Init( &argc, &argv);
	MPI_Comm_rank (comm, &rank);

	// get the name of the file
	strcpy(filename, FILE_NAME);

	// depending on the method
	if( adios_read_init_method( method, comm, ADIOS_OPTIONS) != 0){
		p_error("Quitting ... Issues with initialization adios_read: (%d) %s\n", adios_errno, adios_errmsg());
		return PROGRAM_ERROR;
	}

	// I will be working with streams so the lock mode is necessary,
	// return immediately if the stream unavailable
	ADIOS_FILE *adios_handle = adios_read_open_file(filename, method, comm);
	if ( !adios_handle){
		p_error("Quitting ... (%d) %s\n", adios_errno, adios_errmsg());
		return PROGRAM_ERROR;
	}


	int i = 0;
	// I will only support reading TIMESTEP_COUNT integers for the level value
	int level[TIMESTEP_COUNT];
	int level_scalar[TIMESTEP_COUNT];
	int cctk_bbox[TIMESTEP_COUNT * 6];
	double *data = malloc(8 * 11* 12*13 * TIMESTEP_COUNT);

	memset(level, 0, sizeof(int) * TIMESTEP_COUNT);
	memset(cctk_bbox, 0, sizeof(int) * TIMESTEP_COUNT*6);
	memset(data, 0, sizeof(double) * 11* 12*13 * TIMESTEP_COUNT);



	char fullpath[STR_BUFFER_SIZE];

	// selection should be NULL or as a single variable
	// just say that you want to take different steps
	adios_schedule_read(adios_handle, NULL, "/level/maya_gf_var",0,10, level);
	adios_schedule_read(adios_handle, NULL, "/scalar/maya_gf_var", 0, 10, level_scalar);
	adios_schedule_read(adios_handle, NULL, "/cctk_bbox/maya_gf_var", 0, 10, cctk_bbox);

	// now I will try to read with a single variable but with multiple steps
	ADIOS_SELECTION *sel = NULL;
	uint64_t start_3D[] = {0, 0, 0};
	uint64_t count_3D[] = {11, 12, 13};
	sel = adios_selection_boundingbox(3, start_3D, count_3D );
	adios_schedule_read(adios_handle, sel, "/data/maya_gf_var", 0, 10, data);
	adios_selection_delete(sel);
	sel = NULL;


	TEST_ADIOS_ERROR__IF_NOT_ZERO(adios_perform_reads(adios_handle, 1), error_counts.adios);

	int j = 0;
	if (error_counts.adios){
		test_passed = TEST_FAILED;
	} else {
		// reference data
		int level_ref[TIMESTEP_COUNT];
		gen_1D_array(level_ref, TIMESTEP_COUNT, rank);

		int cctk_bbox_ref[6];
		gen_1D_array(cctk_bbox_ref, 6, rank);

		double * data_ref =  (double *)malloc(11*12*13*sizeof(double));
		gen_1D_array_double(data_ref, 11*12*13, rank);

		// compare with reference values
		for( i = 0; i < TIMESTEP_COUNT; ++i){
			TEST_INT_EQUAL(level_ref[i], level[i], error_counts.test);
			BREAK_IF_ERROR(error_counts.test);
			TEST_INT_EQUAL(level_ref[i], level_scalar[i], error_counts.test);
			BREAK_IF_ERROR(error_counts.test);

			for( j = 0; j < 6; ++j){
				TEST_INT_EQUAL(cctk_bbox_ref[j], cctk_bbox[i*6 + j], error_counts.test);
				BREAK_IF_ERROR(error_counts.test);
			}
			BREAK_IF_ERROR(error_counts.test);

			for( j = 0; j < 11 * 12 * 13; ++j){
				TEST_DOUBLE_EQUAL(data_ref[j], data[i * 11 * 12 * 13 +j], error_counts.test);
				BREAK_IF_ERROR(error_counts.test);
			}
			BREAK_IF_ERROR(error_counts.test);
		}

		free(data_ref);
		data_ref = NULL;
	}

	if (TEST_PASSED == test_passed)
		p_test_passed("%s: rank %d\n", TEST_NAME, rank);

	free(data);
	data = NULL;
	CLOSE_ADIOS(adios_handle, method);

	return test_passed;
}
