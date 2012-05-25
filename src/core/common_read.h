/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 *   Internal read API for C and Fortran read APIs
 */
#ifndef __COMMON_READ_H__
#define __COMMON_READ_H__

#include "public/adios_types.h"
#include "public/adios_read_v2.h"  /* C API's struct's are used here */

#ifdef _NOMPI
    /* Sequential processes can use the library compiled with -D_NOMPI */
#   include "mpidummy.h"
#else
    /* Parallel applications should use MPI to communicate file info and slices of data */
#   include "mpi.h"
#endif

#include <stdint.h>

int common_read_init_method (enum ADIOS_READ_METHOD method, 
                             MPI_Comm comm, 
                             const char * parameters);

int common_read_finalize_method(enum ADIOS_READ_METHOD method);

ADIOS_FILE * common_read_open_stream (const char * fname,
                                     enum ADIOS_READ_METHOD method,
                                     MPI_Comm comm,
                                     enum ADIOS_LOCKMODE lock_mode,
                                     int timeout_msec);

ADIOS_FILE * common_read_open_file   (const char * fname,
                                     enum ADIOS_READ_METHOD method,
                                     MPI_Comm comm);

int common_read_close (const ADIOS_FILE *fp);

int common_read_advance_step (const ADIOS_FILE *fp, int last, int wait_for_step);
void common_read_release_step (const ADIOS_FILE *fp);

ADIOS_VARINFO * common_read_inq_var (const ADIOS_FILE  *fp, const char * varname);
ADIOS_VARINFO * common_read_inq_var_byid (const ADIOS_FILE  *fp, int varid);
int common_read_inq_var_stat (const ADIOS_FILE *fp, const ADIOS_VARINFO * varinfo,
                             int per_step_stat, int per_block_stat);
int common_read_inq_var_blockinfo (const ADIOS_FILE *fp, const ADIOS_VARINFO * varinfo);
void common_read_free_varinfo (const ADIOS_VARINFO *vp);

int common_read_schedule_read (const ADIOS_FILE      * fp,
                               const ADIOS_SELECTION * sel,
                               const char            * varname,
                               int                     from_steps,
                               int                     nsteps,
                               void                  * data);

int common_read_schedule_read_byid (const ADIOS_FILE      * fp,
                                    const ADIOS_SELECTION * sel,
                                    int                     varid,
                                    int                     from_steps,
                                    int                     nsteps,
                                    void                  * data);

int common_read_perform_reads (const ADIOS_FILE *fp, int blocking);
int common_read_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk);
void common_read_free_chunk (const ADIOS_VARCHUNK *chunk);


int common_read_get_attr (const ADIOS_FILE            * fp,
                    const char            * attrname,
                    enum ADIOS_DATATYPES  * type,
                    int                   * size,
                    void                 ** data);

int common_read_get_attr_byid (const ADIOS_FILE  * fp, int attrid, enum ADIOS_DATATYPES * type, 
                         int * size, void ** data); 

const char * common_read_type_to_string (enum ADIOS_DATATYPES type);
int common_read_type_size(enum ADIOS_DATATYPES type, void *data);

int common_read_get_grouplist (const ADIOS_FILE  *fp, char **group_namelist);
int common_read_group_view (const ADIOS_FILE  *fp, int groupid);

// internal function to support version 1 time-dimension reads
int common_read_is_var_timed (const ADIOS_FILE *fp, int varid);

void common_read_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran);
void common_read_print_fileinfo (const ADIOS_FILE *fp);
#endif
