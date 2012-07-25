/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef HW_UTILS_H
#define HW_UTILS_H 1

#include "stdint.h"
#include "hdf5.h"
#include "adios_types.h"

enum lang_convention {
    USE_FORTRAN = 1,
    USE_C = 2
};

enum scalar_convention {
    USE_SCALAR = 1,
    USE_SINGLE_ELE_ARRAY = 2
};

enum verbose_level {
    NO_INFO = 0,
    LIST_INFO = 1,    
    DEBUG_INFO = 2   
};

/*
 * Initialization
 * initialize_bp2h5() sets config variables and allocate
 * internal read buffer.
 * It returns 0 if no error and -1 otherwise
 */
int initialize_bp2h5(enum lang_convention array_dim_order,
                     enum lang_convention var_str,
                     enum lang_convention attr_str_ds,
                     enum lang_convention attr_str_gp,
                     enum scalar_convention scalar_ds,
                     enum verbose_level verb
                     );

/*
 * Convert bp file fnamein to a h5 file fnameout
 */
int hw_makeh5(char *filename, char *filename1);

/*
 * String-typed dataset attribute is converted to hdf5 dataset attribute 
 * in C convention or Fortran convention. 
 */
void hw_string_attr_ds_internal ( hid_t parent_id, const char *name,const char *value);

/*
 * String-typed group attribute is converted to hdf5 dataset attribute 
 * in C convention or Fortran convention. 
 */
void hw_string_attr_gp_internal (hid_t parent_id, const char *name,const char *value);

/*
 * Write string-typed attribute in Fortran convention
 */
void hw_string_attr_f (hid_t parent_id, const char *name,const char *value);

/*
 * Write string-typed attribute in C convention
 */
void hw_string_attr_c (hid_t parent_id, const char *name,const char *value);

/*
 * Write an integer attribute as a scalar or a single-element-array
 */
void hw_scalar_attr(hid_t parent_id, const char *name, void *value,enum ADIOS_DATATYPES type);

/*
 * Write an integer attribute to an h5 file as a scalar
 *
 * Edited from the AVS example file src/hdf5/examp/write_struct.c
 */
void hw_scalar_attr_scalar ( hid_t parent_id, const char *name, void *value, enum ADIOS_DATATYPES type);

/*
 * Write an integer attribute to an h5 file as a single-element-array
 */
void hw_scalar_attr_array ( hid_t parent_id, const char *name, void *value, enum ADIOS_DATATYPES type);

/*
 * Write dataset to h5 file
 */
/*
 * Write dataset to h5 file
 */
void hw_dset (hid_t root_id, char * dirstr, char * name, void * data
             ,enum ADIOS_DATATYPES type, int rank
             ,uint64_t *dims
             ,uint64_t *global_dims
             ,uint64_t *offsets
             ,uint32_t time_index 
             );

/*
 * Write a scalar var to h5 file as a scalar
 */
void hw_scalar_as_scalar (hid_t parent_id, char * name, void * data
                 ,enum ADIOS_DATATYPES type, int append
                 );

/*
 * Write a scalar var to h5 file as a single-element array
 */
void hw_scalar_as_array (hid_t parent_id, char * name, void * data
                 ,enum ADIOS_DATATYPES type, int ndims, hsize_t * dims
                 ,int append
                 );

/*
 * write a scalar var to h5 file
 */
void hw_scalar (hid_t root_id, char * dirstr, char * name, void * data
               ,enum ADIOS_DATATYPES type, int append
               );

/*
 * hw_attr_str_ds() writes string-typed dataset attribute 
 * hw_attr_str_ds() finds the dataset to which the attribute is attached and
 * writes the attribute. 
 */
void hw_attr_str_ds (hid_t root_id, char * dirstr, char * aname, char * aval);

/*
 * hw_attr_num_ds() writes numeric-typed dataset attribute 
 * hw_attr_num_ds() finds the dataset to which the attribute is attached and
 * writes the attribute. 
 */
void hw_attr_num_ds(hid_t root_id, char *dirstr, char *aname, void *avalue, enum ADIOS_DATATYPES type);

/*
 * hw_attr_str_gp() writes string-typed group attribute 
 * hw_attr_str_gp() finds the group to which the attribute is attached and
 * writes the attribute. 
 */
void hw_attr_str_gp (hid_t root_id, char * dirstr, char * aname, char * aval);

/*
 * hw_attr_num_gp() writes numeric-typed group attribute 
 * hw_attr_num_gp() finds the group to which the attribute is attached and
 * writes the attribute. 
 */
void hw_attr_num_gp(hid_t root_id, char *dirstr, char *aname, void *avalue, enum ADIOS_DATATYPES type);

/*
 * Write an array as a dataset to an h5 file 
 */
void hw_dataset(hid_t parent_id, char* name, void* data,enum ADIOS_DATATYPES type, int ndims, hsize_t * dims);

/*
 * Maps bp datatypes to h5 datatypes 
 */
 int bp_getH5TypeId(enum ADIOS_DATATYPES type, hid_t* h5_type_id, void * val);

#endif
