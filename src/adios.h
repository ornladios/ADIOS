/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_H
#define ADIOS_H

#include "adios_types.h"
#include <stdint.h>

#ifdef _NOMPI
/* Sequential processes can use the library compiled with -D_NOMPI */
#   include "mpidummy.h"
#else
/* Parallel applications should use MPI to communicate file info and slices of data */
#   include "mpi.h"
#endif

// ADIOS - Adaptable IO System

#ifdef __cplusplus
extern "C" {
#endif

    // Global setup using the XML file
    int adios_init (const char * config);

    int adios_finalize (int mype);

    // end user calls for each I/O operation
    // modes = "r" = "read", "w" = "write", "a" = "append", "u" = "update"
    int adios_open (int64_t * fd, const char * group_name, const char * name
            ,const char * mode, void * comm
            );

    int adios_group_size (int64_t fd_p, uint64_t data_size
            ,uint64_t * total_size
            );

    int adios_write (int64_t fd_p, const char * name, void * var);

    int adios_get_write_buffer (int64_t fd_p, const char * name
            ,uint64_t * size
            ,void ** buffer
            );

    int adios_read (int64_t fd_p, const char * name, void * buffer
            ,uint64_t buffer_size
            );

    int adios_set_path (int64_t fd_p, const char * path);

    int adios_set_path_var (int64_t fd_p, const char * path, const char * name);

    int adios_end_iteration (void);

    int adios_start_calculation (void);

    int adios_stop_calculation (void);

    int adios_close (int64_t fd_p);

    // ADIOS No-XML API's
    int adios_init_noxml (void);

    // To allocate ADIOS buffer
    int adios_allocate_buffer (enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when
            ,uint64_t buffer_size);

    // To declare a ADIOS group
    int adios_declare_group (int64_t * id, const char * name
            ,const char * time_index
            ,enum ADIOS_FLAG stats
            );
    // To free a ADIOS group
    int adios_free_group (int64_t id);

    // To select a I/O method for a ADIOS group
    int adios_select_method (int64_t group, const char * method
            ,const char * parameters
            ,const char * base_path
            );

    // To define a ADIOS variable
    int adios_define_var (int64_t group_id, const char * name
            ,const char * path, int type
            ,const char * dimensions
            ,const char * global_dimensions
                     ,const char * local_offsets
                     );

int adios_define_attribute (int64_t group, const char * name
                           ,const char * path, enum ADIOS_DATATYPES type
                           ,const char * value, const char * var
                           );

/** Set the application's ID for adios_read_init()
 *  when using a staging method (DART, DIMES, NSSI or DATATAP).
 *  The ID should be unique for each application accessing the staging area
 *  IN:  id   a number unique for this application
 *  RETURN:       0 if accepted, <0 on error
 *  It is optional to use it before calling adios_init. Default is 1. 
 *  It has no effect for file based methods.
 *  Note: this function is defined both in adios.h and adios_read.h so that
 *  writing-only and reading-only applications can both use it.
 */ 
int adios_set_application_id (int id);


/********************************************/
/*           F O R T R A N  A P I           */
/********************************************/
/* In fortran, you do not need to include this header file.
   Just link the code with the adios library.
*/
/*
void FC_FUNC_(adiosf_init, ADIOSF_INIT) (const char * config, int* err, int config_size);
void FC_FUNC_(adiosf_init_local, ADIOSF_INIT_LOCAL) (int * err);
void FC_FUNC_(adiosf_finalize, ADIOSF_FINALIZE) (int * mype, int * err);
void FC_FUNC_(adiosf_allocate_buffer, ADIOSF_ALLOCATE_BUFFER) (int * err);
void FC_FUNC_(adiosf_open, ADIOSF_OPEN) (int64_t * fd, const char * group_name, const char * name
                 ,const char * mode, int * err
                 ,int group_name_size, int name_size, int mode_size
                 );
void FC_FUNC_(adiosf_group_size, ADIOSF_GROUP_SIZE) (int64_t * fd_p, int64_t * data_size
                       ,int64_t * total_size, void * comm, int * err
                       );
void FC_FUNC_(adiosf_write, ADIOSF_WRITE) (int64_t * fd_p, const char * name, void * var, int * err
                  ,int name_size
                  );
void FC_FUNC_(adiosf_get_write_buffer, ADIOSF_GET_WRITE_BUFFER) (int64_t * fd_p, const char * name
                             ,int64_t * size
                             ,void ** buffer, int * err, int name_size
                             );
void FC_FUNC_(adiosf_read, ADIOSF_READ) (int64_t * fd_p, const char * name, void * buffer
                 ,int64_t * buffer_size, int * err, int name_size
                 );
void FC_FUNC_(adiosf_set_path, ADIOSF_SET_PATH) (int64_t * fd_p, const char * path, int * err
                     ,int path_size
                     );
void FC_FUNC_(adiosf_set_path_var, ADIOSF_SET_PATH_VAR) (int64_t * fd_p, const char * path
                         ,const char * name, int * err
                         ,int path_size, int name_size
                         );
void FC_FUNC_(adiosf_end_iteration, ADIOSF_END_ITERATION) (int * err);
void FC_FUNC_(adiosf_start_calculation, ADIOSF_START_CALCULATION) (int * err);
void FC_FUNC_(adiosf_stop_calculation, ADIOSF_STOP_CALCULATION) (int * err);
void FC_FUNC_(adiosf_close, ADIOSF_CLOSE) (int64_t * fd_p, int * err);
void FC_FUNC_(adiosf_declare_group, ADIOSF_DECLARE_GROUP) (int64_t * id, const char * name
                          ,const char * coordination_comm
                          ,const char * coordination_var
                          ,const char * time_index, int * err
                          ,int name_size, int coordination_comm_size
                          ,int coordination_var_size, int time_index_size
                          );
void FC_FUNC_(adiosf_define_var, ADIOSF_DEFINE_VAR) (int64_t * group_id, const char * name
                       ,const char * path, int * type
                       ,const char * dimensions
                       ,const char * global_dimensions
                       ,const char * local_offsets, int * err
                       ,int name_size, int path_size, int dimensions_size
                       ,int global_dimensions_size, int local_offsets_size
                       );
void FC_FUNC_(adiosf_define_attribute, ADIOSF_DEFINE_ATTRIBUTE) (int64_t * group, const char * name
                             ,const char * path, int type, const char * value
                             ,const char * var, int * err
                             ,int name_size, int path_size, int value_size
                             ,int var_size
                             );
void FC_FUNC_(adiosf_select_method, ADIOSF_SELECT_METHOD) (int * priority, const char * method
                          ,const char * parameters, const char * type
                          ,const char * base_path, int * iters, int * err
                          ,int method_size
                          ,int parameters_size, int type_size
                          ,int base_path_size
                          );
void FC_FUNC_(adiosf_select_method, ADIOSF_SELECT_METHOD) (int * priority, const char * method
                          ,const char * parameters, const char * type
                          ,const char * base_path, int * iters, int * err
                          ,int method_size
                          ,int parameters_size, int type_size
                          ,int base_path_size
                          );
*/

#ifdef __cplusplus
}
#endif


#endif
