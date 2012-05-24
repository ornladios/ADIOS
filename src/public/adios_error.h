/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef __ADIOS_ERROR_H_
#define __ADIOS_ERROR_H_

enum ADIOS_ERRCODES {
     err_no_error = 0
    ,err_no_memory
    ,err_MPI_open_error
    ,err_file_not_found_error
    ,err_invalid_file_pointer
    ,err_invalid_group
    ,err_invalid_group_struct
    ,err_invalid_varid
    ,err_invalid_varname
    ,err_corrupted_variable
    ,err_invalid_attrid
    ,err_invalid_attrname
    ,err_corrupted_attribute
    ,err_invalid_attribute_reference
    ,err_invalid_timestep
    ,err_no_data_at_timestep
    ,err_time_at_wrong_dimension
    ,err_invalid_read_method
    ,err_connection_failed
    ,err_out_of_bound
    ,err_operation_not_supported
    ,err_end_of_file    // stream: fopen() returns if reached end of stream
    ,err_too_many_files  // DART allows for using only a fixed number of different filenames
    ,err_unspecified
};

void adios_error (enum ADIOS_ERRCODES errno, char *fmt, ...);
void adios_error_at_line (enum ADIOS_ERRCODES errno, const char* filename, unsigned int linenum, char *fmt, ...);

const char* adios_get_last_errmsg (void);

#endif
