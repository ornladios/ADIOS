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
#include "config.h"

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
#else
	printf("MPI/ADIOS_READ_METHOD_BP");
#endif
	printf(" will be used.\n");
	printf("  nprocs Number of processors you want to use\n");
	printf("E.g. mpirun -np 5 %s\n", program_name);

	return 0;
}

/**
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
int gen_1D_array(double *p_arr, int arr_len, int rank){

	if (!p_arr){
		p_error("p_arr is NULL\n.");
		return -1;
	}
	int i = 0;
	for( i = 0; i < arr_len; i++){
		p_arr[i] = rank * 10 + i;
	}

	return 0;
}

/**
 * Sets values the 1D array of length arr_len with a specified value
 * It assumes that the memory is allocated for p_arr
 *
 * The function does not check if the memory is overwritten
 *
 * @param p_arr The pointer to the array that will hold the values
 * @param arr_len The number of elements in the array
 * @param value The value to be set
 *
 * @return -1 if the p_arr is NULL
 *          0 otherwise
 *          p_arr (with generated numbers)
 */
int set_value(double *p_arr, int arr_len, double value){

	if (!p_arr){
		p_error("p_arr is NULL\n.");
		return -1;
	}
	int i = 0;
	for( i = 0; i < arr_len; i++){
		p_arr[i] = value;
	}

	return 0;
}

/**
 * Get the size of the data. The function assumes that the data will be
 * of a "double" type
 *
 * @param shape The array containing the shape values
 * @param shape_elem_count How many elements is in th the shape array
 * @param data_size (OUT) This value is computed based on the shape array
 *
 * @return -1 if the shape is null
 *          0 if the data_size contains the valid value
 */
int get_data_size(int *shape, int shape_elem_count, int* data_size){

	if (!shape){
		p_error("shape is NULL\n");
		return -1;
	}

	int i = 0;
	int my_size = 8; //sizeof(double);

	for(i = 0; i < shape_elem_count; ++i){
		my_size *= shape[i];
	}

	*data_size = my_size;

	return 0;
}

/**
 * Generates the maya var name based on the value of MAYA_GF_VAR_PFX macro and the number.
 * Now it generates: MAYA_GF_VAR_PFX0, MAYA_GF_VAR_PFX1, ...
 *
 * The function doesn't protect against too small arra
 * @param (in/out) buf  the buffer for holding the maya variable; it should be big enough
 *                to hold additional characters that will create the final
 *                name of the variable; it is cleaned
 * @param (in) buf_size The size of the buffer in characters
 * @param (in) number  The number that will be concatenated with the prefix;
 *                     the number should be positive
 *
 * @return 0 everything is ok
 *         -1 the number is negative
 */
int gen_maya_var_name(char *buf, int buf_size, int number){
	if (number < 0){
		return -1;
	}

	memset(buf, 0, buf_size);
	sprintf(buf, MAYA_GF_VAR_PFX "%d", number);

	return 0;
}

