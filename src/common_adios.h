/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef COMMON_ADIOS_H
#define COMMON_ADIOS_H

#include <stdint.h>
#include "adios_types.h"
#include "adios_internals.h"

/* Write functions for ADIOS
 *
 * Used by the write C api (adios.c) and Fortran api (adiosf.c)
 */

// Global setup using the XML file
int common_adios_init (const char * config);

// setup, but all XML file pieces will be provided by another series of calls
// yet to be worked out
// TODO
int common_adios_init_noxml (void);

int common_adios_finalize (int mype);

int common_adios_allocate_buffer (enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when
                                 ,uint64_t buffer_size);

// end user calls for each I/O operation
// modes = "r" = "read", "w" = "write", "a" = "append", "u" = "update"
int common_adios_open (int64_t * fd, const char * group_name, const char * name
               ,const char * mode, void * comm
               );

int common_adios_group_size (int64_t fd_p, uint64_t data_size
                     ,uint64_t * total_size
                     );

//int common_adios_write (int64_t fd_p, const char * name, void * var);
int common_adios_write (struct adios_file_struct * fd, struct adios_var_struct * v, void * var);

int common_adios_get_write_buffer (int64_t fd_p, const char * name
                           ,uint64_t * size
                           ,void ** buffer
                           );

int common_adios_read (int64_t fd_p, const char * name, void * buffer
               ,uint64_t buffer_size
               );

int common_adios_set_path (int64_t fd_p, const char * path);

int common_adios_set_path_var (int64_t fd_p, const char * path, const char * name);

int common_adios_end_iteration (void);

int common_adios_start_calculation (void);

int common_adios_stop_calculation (void);

int common_adios_close (int64_t fd_p);

// Generally internal use called when parsing the XML file
int common_adios_declare_group (int64_t * id, const char * name
                        ,const char * coordination_comm
                        ,const char * coordination_var
                        ,const char * time_index
                        );

int common_adios_define_var (int64_t group_id, const char * name
                     ,const char * path, int type
                     ,const char * dimensions
                     ,const char * global_dimensions
                     ,const char * local_offsets
                     );

int common_adios_define_attribute (int64_t group, const char * name
                           ,const char * path, enum ADIOS_DATATYPES type
                           ,const char * value, const char * var
                           );

int common_adios_select_method (int priority, const char * method
                        ,const char * parameters, const char * type
                        ,const char * base_path, int iters
                        );

int common_adios_select_method (int priority, const char * method
                        ,const char * parameters, const char * type
                        ,const char * base_path, int iters
                        );

#endif
