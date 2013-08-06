/**
 * utils.h
 *
 *  Created on: Jul 5, 2013
 *  Author: Magda Slawinska aka Magic Magg magg dot gatech at gmail.com
 */

#ifndef UTILS_H_
#define UTILS_H_

extern int usage(char *program_name, char *program_desc);
extern int gen_1D_array(double *p_arr, int arr_len, int rank);
extern int set_value(double *p_arr, int arr_len, double value);

extern int get_data_size(int *shape, int shape_elem_count, int* data_size);
extern int get_maya_var_name(char *prefix, int number);

#endif /* UTILS_H_ */
