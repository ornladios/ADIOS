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
	printf("ADIOS_READ_METHOD_FLEXPATH" );
#else
	printf("ADIOS_READ_METHOD_BP+MPI");
#endif
	printf(" will be used.\n");
	printf("Please check " XML_ADIOS_INIT_FILENAME "\n");
	printf("  nprocs Number of processors you want to use\n");
	printf("E.g. mpirun -np 5 %s\n", program_name);

	return 0;
}

/**
 * Generates the 1D array of length arr_len, based on the provided rank
 *
 * 	 p_arr : rank 0: 1, 2, 3, 4, 5, ....
 *           rank 1: 10, 11, 12, 13, ...
 *           rank 2: 20, 21, 22, 23, ...
 *
 * The function does not check if arr_len does check if the memory is overwritten
 *
 * @param p_arr The pointer to the array that will hold the values
 * @param arr_len The number of elements in the array
 * @param rank The rank for which I want to have the number generated
 *
 * @return -1 if the p_arr is NULL
 *          0 otherwise
 *          p_arr (with generated numbers)
 */
int gen_1D_array(double *p_arr, int arr_len, int rank){

	if (!p_arr){
		printf("ERROR: p_arr is NULL\n.");
		return -1;
	}
	int i = 0;
	for( i = 0; i < arr_len; i++){
		p_arr[i] = rank * 10 + i;
	}

	return 0;
}
