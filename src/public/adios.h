/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_H
#define ADIOS_H

#include "adios_mpi.h"
#include "adios_types.h"
#include "adios_error.h"
#include <stdint.h>

// ADIOS - Adaptable IO System

#ifdef __cplusplus
extern "C" {
#endif

/* Most functions return 0 if OK, and !=0 on error 
   which is the value of the variable 'adios_errno'.
   On error, one can use char * adios_get_last_errmsg() from adios_error.h
   to retrieve the error string of the last error. 

   exceptions: int64_t adios_define_var() returns a variable ID, 0 indicates an error
*/

// Global setup using the XML file
// Only processes of the provided communicator can later participate
// in any adios activity
int adios_init (const char * config, MPI_Comm comm);

int adios_finalize (int mype);

// end user calls for each I/O operation
// modes = "r" = "read", "w" = "write", "a" = "append", "u" = "update"
int adios_open (int64_t * fd, 
                const char * group_name, 
                const char * name,
                const char * mode, 
                MPI_Comm comm
               );

int adios_group_size (int64_t fd_p, 
                      uint64_t data_size,
                      uint64_t * total_size
                     ); 

int adios_write (int64_t fd_p, const char * name, void * var);

int adios_get_write_buffer (int64_t fd_p, 
                            const char * name,
                            uint64_t * size,
                            void ** buffer
                           );

int adios_read (int64_t fd_p, 
                const char * name, 
                void * buffer,
                uint64_t buffer_size
               );

int adios_set_path (int64_t fd_p, const char * path);

int adios_set_path_var (int64_t fd_p, const char * path, const char * name);

int adios_end_iteration (void);

int adios_start_calculation (void);

int adios_stop_calculation (void);

int adios_close (int64_t fd_p);

// ADIOS No-XML API's
int adios_init_noxml (MPI_Comm comm);

// To allocate ADIOS buffer
int adios_allocate_buffer (
        enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when,
        uint64_t buffer_size
        );

// To declare a ADIOS group
int adios_declare_group (int64_t * id, 
                         const char * name,
                         const char * time_index, 
                         enum ADIOS_FLAG stats
                        );

// To free a ADIOS group
int adios_free_group (int64_t id);

// To select a I/O method for a ADIOS group
int adios_select_method (int64_t group, 
                         const char * method,
                         const char * parameters,
                         const char * base_path
                        );

// To define a ADIOS variable
// Returns a variable ID, which can be used in adios_write_byid()
// 0 return value indicates an error
int64_t adios_define_var (int64_t group_id, 
                          const char * name,
                          const char * path,
                          enum ADIOS_DATATYPES type,
                          const char * dimensions,
                          const char * global_dimensions,
                          const char * local_offsets
                         );

int adios_define_attribute (int64_t group, 
                            const char * name,
                            const char * path, 
                            enum ADIOS_DATATYPES type,
                            const char * value, 
                            const char * var
                           );
/** This function does similar function as adios_write. It is, however, used
 * in the following scenario that
 * 1. numbers, instead of a variable, are used to annotate array dimensions, and
 * 2. a variable is written mutiple times on a processor (e.g., AMR codes)
 */
int adios_write_byid (int64_t fd_p, int64_t id, void * var);

/** Set the application's ID for adios_read_init()
 *  when using a staging method (DATASPACES, DIMES, NSSI or DATATAP).
 *  The ID should be unique for each application accessing the staging area
 *  IN:  id   a number unique for this application
 *  RETURN:       0 if accepted, <0 on error
 *  It is optional to use it before calling adios_init. Default is 1. 
 *  It has no effect for file based methods.
 *  Note: this function is defined both in adios.h and adios_read.h so that
 *  writing-only and reading-only applications can both use it.
 */ 
/*int adios_set_application_id (int id);*/

void adios_timing_write_xml (int64_t fd_p, const char* filename);


#ifdef __cplusplus
}
#endif


#endif
