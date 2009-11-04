#ifndef __ADIOS_ERROR_H_
#define __ADIOS_ERROR_H_

enum ADIOS_ERRCODES {
     err_no_error = 0
    ,err_no_memory
    ,err_MPI_open_error
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
    ,err_out_of_bound
};

void error (enum ADIOS_ERRCODES errno, char *fmt, ...);

const char* adios_get_last_errmsg (void);

#endif
