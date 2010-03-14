/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_INTERNALS_V2_H
#define ADIOS_INTERNALS_V2_H

#include <stdint.h>
#include <stdlib.h>

// need the enum for the transports
#include "adios_transport_hooks.h"
#include "adios_bp_v2.h"

// ADIOS file format functions
uint16_t adios_calc_var_overhead_v2 (struct adios_var_struct * v);
uint32_t adios_calc_attribute_overhead_v2 (struct adios_attribute_struct * a);
uint64_t adios_calc_overhead_v2 (struct adios_file_struct * fd);

int adios_write_version_v2 (char ** buffer
                           ,uint64_t * buffer_size
                           ,uint64_t * buffer_offset
                           );
int adios_write_process_group_header_v2 (struct adios_file_struct * fd
                                        ,uint64_t total_size
                                        );

// data is only there for sizing
uint64_t adios_write_var_header_v2 (struct adios_file_struct * fd
                                   ,struct adios_var_struct * v
                                   );
int adios_generate_var_characteristics_v2 (struct adios_file_struct * fd
                                          ,struct adios_var_struct * var
                                          );
uint16_t adios_write_var_characteristics_v2 (struct adios_file_struct * fd
                                            ,struct adios_var_struct * var
                                            );
int adios_write_var_payload_v2 (struct adios_file_struct * fd
                               ,struct adios_var_struct * var
                               );
int adios_write_attribute_v2 (struct adios_file_struct * fd
                             ,struct adios_attribute_struct * a
                             );
int adios_write_open_vars_v2 (struct adios_file_struct * fd);
int adios_write_close_vars_v2 (struct adios_file_struct * fd);
int adios_write_open_attributes_v2 (struct adios_file_struct * fd);
int adios_write_close_attributes_v2 (struct adios_file_struct * fd);
int adios_write_index_v2 (char ** buffer
                         ,uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,uint64_t index_start
                         ,struct adios_index_process_group_struct_v2 * pg_root
                         ,struct adios_index_var_struct_v2 * vars_root
                         ,struct adios_index_attribute_struct_v2 * attrs_root
                         );

void adios_build_index_v2 (struct adios_file_struct * fd
                       ,struct adios_index_process_group_struct_v2 ** pg_root
                       ,struct adios_index_var_struct_v2 ** vars_root
                       ,struct adios_index_attribute_struct_v2 ** attrs_root
                       );
void adios_merge_index_v2 (
                   struct adios_index_process_group_struct_v2 ** main_pg_root
                  ,struct adios_index_var_struct_v2 ** main_vars_root
                  ,struct adios_index_attribute_struct_v2 ** main_attrs_root
                  ,struct adios_index_process_group_struct_v2 * new_pg_root
                  ,struct adios_index_var_struct_v2 * new_vars_root
                  ,struct adios_index_attribute_struct_v2 * new_attrs_root
                  );
void adios_clear_index_v2 (struct adios_index_process_group_struct_v2 * pg_root
                          ,struct adios_index_var_struct_v2 * vars_root
                          ,struct adios_index_attribute_struct_v2 * attrs_root
                          );

uint64_t adios_get_type_size (enum ADIOS_DATATYPES type, void * var);
uint64_t adios_get_var_size (struct adios_var_struct * var
                            ,struct adios_group_struct * group, void * data
                            );

const char * adios_type_to_string_int (int type);
const char * adios_file_mode_to_string (int mode);

// the following are defined in adios_transport_hooks.c
void adios_init_transports (struct adios_transport_struct ** transports);
int adios_parse_method (const char * buf, enum ADIOS_IO_METHOD * method
                       ,int * requires_group_comm
                       );

/* some internal functions that adios_internals.c and adios_internals_mxml.c share */
int adios_int_is_var (const char * temp); // 1 == yes, 0 == no
int adios_int_is_num (char * temp); // 1 == yes, 0 == no

#endif
