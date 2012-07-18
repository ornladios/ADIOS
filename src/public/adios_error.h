/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef __ADIOS_ERROR_H_
#define __ADIOS_ERROR_H_

enum ADIOS_ERRCODES {
    err_no_error                        = 0,
    err_no_memory                       = -1,
    err_file_open_error                 = -2,
    err_file_not_found                  = -3,
    err_invalid_file_pointer            = -4,
    err_invalid_group                   = -5,
    err_invalid_group_struct            = -6,
    err_invalid_varid                   = -7,
    err_invalid_varname                 = -8,
    err_corrupted_variable              = -9,

    err_invalid_attrid                  = -10,
    err_invalid_attrname                = -11,
    err_corrupted_attribute             = -12,
    err_invalid_attribute_reference     = -13,
    err_invalid_timestep                = -14,
    err_no_data_at_timestep             = -15,
    err_time_at_wrong_dimension         = -16,
    err_invalid_read_method             = -17,
    err_connection_failed               = -18,
    err_out_of_bound                    = -19,

    err_operation_not_supported         = -20,
    err_end_of_stream                   = -21, /* stream: reached end of stream, 
                                                  returned by adios_read_open() or
                                                           by adios_advance_step()         */
    err_step_notready                   = -22, /* stream: tried to advance the step, 
                                                  but no new step has arrived yet          */
    err_step_disappeared                = -23, /* stream: tried to advance the step, 
                                                  but next step has already disappeared    */
    err_too_many_files                  = -24, /* some staging methods allow for using only
                                                  a fixed number of different filenames    */

    err_unknown_char                    = -30,

    // XML parsing errors
    err_invalid_host_language           = -50,
    err_global_dim_required             = -51,
    err_global_offset_required          = -52,
    err_invalid_method                  = -53,
    err_invalid_buffer_size             = -54,
    err_missing_config_file             = -55,
    err_expected_read_size_mismatch     = -56,
    err_allocating_buffer_size          = -57,
    err_invalid_xml_doc                 = -58,
    err_no_group_defined                = -59,
    err_no_method_defined               = -60,
    err_no_buffer_defined               = -61,
    err_missing_invalid_group           = -62,
    err_group_method_mismatch           = -63,

    // Write method errors
    err_invalid_file_mode               = -100,
    err_invalid_file_version            = -101,
    err_invalid_data                    = -102,
    err_buffer_overflow                 = -103,
    err_too_many_variables              = -104,
    err_invalid_write_method            = -105,

    //buffering errors
    err_invalid_buffer                  = -130,
    err_invalid_buffer_version          = -131,
    err_invalid_buffer_index            = -132,
    err_invalid_buffer_group            = -133,
    err_invalid_buffer_vars             = -134,
    err_invalid_buffer_attrs            = -135,
    
    

    //invalid argument to function
    err_invalid_argument                = -140,


    err_unspecified                     = -200
};

void adios_error (enum ADIOS_ERRCODES errcode, char *fmt, ...);
void adios_error_at_line (enum ADIOS_ERRCODES errcode, const char* filename, unsigned int linenum, char *fmt, ...);

const char* adios_get_last_errmsg (void);

#endif
