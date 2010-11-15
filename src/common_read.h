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

#include "adios_types.h"
#include "adios_read.h"  /* C API's struct's are used here */

#ifdef _NOMPI
    /* Sequential processes can use the library compiled with -D_NOMPI */
#   include "mpidummy.h"
#else
    /* Parallel applications should use MPI to communicate file info and slices of data */
#   include "mpi.h"
#endif

#include <stdint.h>

int common_read_set_read_method (enum ADIOS_READ_METHOD method);
ADIOS_FILE * common_read_fopen (const char * fname, MPI_Comm comm);
int common_read_fclose (ADIOS_FILE *fp);
void common_read_reset_dimension_order (ADIOS_FILE *fp, int is_fortran);
ADIOS_GROUP * common_read_gopen (ADIOS_FILE *fp, const char * grpname);
ADIOS_GROUP * common_read_gopen_byid (ADIOS_FILE *fp, int grpid);
int common_read_gclose (ADIOS_GROUP *gp);
ADIOS_VARINFO * common_read_inq_var (ADIOS_GROUP *gp, const char * varname);
ADIOS_VARINFO * common_read_inq_var_byid (ADIOS_GROUP *gp, int varid);
void common_read_free_varinfo (ADIOS_VARINFO *cp);
int64_t common_read_read_var (ADIOS_GROUP    * gp, 
                        const char     * varname, 
                        const uint64_t * start,
                        const uint64_t * count, 
                        void           * data);
int64_t common_read_read_local_var (ADIOS_GROUP    * gp,
                                    const char     * varname,
                                    int            idx,
                                    const uint64_t * start,
                                    const uint64_t * count,
                                    void           * data);
int64_t common_read_read_var_byid (ADIOS_GROUP * gp, int varid,
                             const uint64_t * start, const uint64_t * count, 
                             void * data);
int common_read_get_attr (ADIOS_GROUP           * gp,
                    const char            * attrname,
                    enum ADIOS_DATATYPES  * type,
                    int                   * size,
                    void                 ** data);

int common_read_get_attr_byid (ADIOS_GROUP * gp, int attrid, enum ADIOS_DATATYPES * type, 
                         int * size, void ** data); 

const char * common_read_type_to_string (enum ADIOS_DATATYPES type);
int common_read_type_size(enum ADIOS_DATATYPES type, void *data);

void common_read_print_groupinfo (ADIOS_GROUP *gp);
void common_read_print_fileinfo (ADIOS_FILE *fp);
#endif
