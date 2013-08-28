/**
 * utils.c
 *
 *  Created on: Jul 5, 2013
 *  Author: Magda Slawinska aka Magic Magg magg dot gatech at gmail.com
 *
 *  Utility functions.
 */
#include <stdio.h>
#include "misc.h"
#include "test_common.h"

/**
 * @param program_name The name of the program
 * @param program_desc The description of the program
 * @return 0 Returns always 0
 */
int usage(char *program_name, char *program_desc){
	printf("Usage: mpirun -np nprocs %s\n", program_name);
	printf("  %s\n", program_desc);
	printf("  METHOD: ");
#ifdef FLEXPATH_METHOD
	printf("FLEXPATH/ADIOS_READ_METHOD_FLEXPATH" );
#else`
	printf("MPI/ADIOS_READ_METHOD_BP");
#endif
	printf(" will be used.\n");
	printf("  nprocs Number of processors you want to use\n");
	printf("E.g. mpirun -np 5 %s\n", program_name);

	return 0;
}

/**
 * TODO this is not true; each rank generates the same sequence of numbers
 * Generates the 1D array of length arr_len, based on the provided rank
 * It assumes that the memory is allocated for p_arr
 *
 * 	 p_arr : rank 0: 1, 2, 3, 4, 5, ....
 *           rank 1: 10, 11, 12, 13, ...
 *           rank 2: 20, 21, 22, 23, ...
 *
 * The function does not check if the memory is overwritten
 *
 * @param p_arr The pointer to the array that will hold the values
 * @param arr_len The number of elements in the array
 * @param rank The rank for which I want to have the number generated
 *
 * @return -1 if the p_arr is NULL
 *          0 otherwise
 *          p_arr (with generated numbers)
 */
int gen_1D_array(int *p_arr, int arr_len, int rank){

	if (!p_arr){
		p_error("p_arr is NULL\n.");
		return -1;
	}
	int i = 0;
	for( i = 0; i < arr_len; i++){
		// the diversification with rank will not work for as the last one will be used
		// if more than 1 ranks will be used and the reader will signal the error
		// that's why I decided to go with the same value independent on
		// the rank
		// TODO but the above and below is something to explore
		//p_arr[i] = rank * 10 + i;
		p_arr[i] = i;
	}

	return 0;
}

/**
 * TODO each rank generates the same sequence of numbers
 * Generates the 1D array of length arr_len, based on the provided rank
 * It assumes that the memory is allocated for p_arr
 *
 * 	 p_arr : rank 0: 1, 2, 3, 4, 5, ....
 *           rank 1: 10, 11, 12, 13, ...
 *           rank 2: 20, 21, 22, 23, ...
 *
 * The function does not check if the memory is overwritten
 *
 * @param p_arr The pointer to the array that will hold the values
 * @param arr_len The number of elements in the array
 * @param rank The rank for which I want to have the number generated
 *
 * @return -1 if the p_arr is NULL
 *          0 otherwise
 *          p_arr (with generated numbers)
 */
int gen_1D_array_double(double *p_arr, int arr_len, int rank){

	if (!p_arr){
		p_error("p_arr is NULL\n.");
		return -1;
	}
	int i = 0;
	for( i = 0; i < arr_len; i++){
		// see comment in @gen_1D_array
		//p_arr[i] = rank * 10.0 + i;
		p_arr[i] = i * 1.0;
	}

	return 0;
}
