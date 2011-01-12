/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_READ_HOOKS_H
#define ADIOS_READ_HOOKS_H

#include "config.h"
#include <stdint.h>
#include <string.h>
#include "adios_read.h"

#define FORWARD_DECLARE(a) \
int adios_read_##a##_init (MPI_Comm comm); \
int adios_read_##a##_finalize (); \
int adios_read_##a##_fclose (ADIOS_FILE *fp); \
ADIOS_FILE * adios_read_##a##_fopen (const char * fname, MPI_Comm comm); \
int adios_read_##a##_fclose (ADIOS_FILE *fp); \
ADIOS_GROUP * adios_read_##a##_gopen (ADIOS_FILE *fp, const char * grpname); \
ADIOS_GROUP * adios_read_##a##_gopen_byid (ADIOS_FILE *fp, int grpid); \
int adios_read_##a##_gclose (ADIOS_GROUP *gp); \
ADIOS_VARINFO * adios_read_##a##_inq_var (ADIOS_GROUP *gp, const char * varname); \
ADIOS_VARINFO * adios_read_##a##_inq_var_byid (ADIOS_GROUP *gp, int varid); \
int64_t adios_read_##a##_read_var (ADIOS_GROUP * gp, const char * varname, const uint64_t * start, const uint64_t * count, void * data); \
int64_t adios_read_##a##_read_local_var (ADIOS_GROUP * gp, const char * varname, int idx, const uint64_t * start, const uint64_t * count, void * data); \
int64_t adios_read_##a##_read_var_byid (ADIOS_GROUP * gp, int varid, const uint64_t * start, const uint64_t * count, void * data); \
int adios_read_##a##_get_attr (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type, int * size, void ** data); \
int adios_read_##a##_get_attr_byid (ADIOS_GROUP * gp, int attrid, enum ADIOS_DATATYPES * type, int * size, void ** data); \
void adios_read_##a##_reset_dimension_order (ADIOS_FILE *fp, int is_fortran); \

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//// SETUP YOUR NEW READ METHODS BELOW (FOLLOW THE PATTERN):                  ////
//// 1. Add an entry to the adios_read.h/ADIOS_READ_METHOD                    ////
//// 2. Update the ADIOS_METHOD_COUNT                                         ////
//// 2. Add a FOWARD_DECLARE line (assuming standard naming)                  ////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

#define ADIOS_READ_METHOD_COUNT 7

// forward declare the functions (or dummies for internals use)
FORWARD_DECLARE(bp)
FORWARD_DECLARE(bp_subfile)
FORWARD_DECLARE(hdf5)
#if HAVE_DART
FORWARD_DECLARE(dart)
#endif
#if HAVE_DIMES
FORWARD_DECLARE(dimes)
#endif
#if HAVE_NSSI
FORWARD_DECLARE(nssi)
#endif
#if HAVE_NSSI
FORWARD_DECLARE(datatap)
#endif


typedef int (* ADIOS_INIT_FN) (MPI_Comm comm);
typedef int (* ADIOS_FINALIZE_FN) ();
typedef ADIOS_FILE * (* ADIOS_FOPEN_FN) (const char * fname, MPI_Comm comm);
typedef int (* ADIOS_FCLOSE_FN) (ADIOS_FILE *fp);
typedef ADIOS_GROUP * (* ADIOS_GOPEN_FN) (ADIOS_FILE *fp, const char * grpname);
typedef ADIOS_GROUP * (* ADIOS_GOPEN_BYID_FN) (ADIOS_FILE *fp, int grpid);
typedef int (* ADIOS_GCLOSE_FN) (ADIOS_GROUP *gp);
typedef ADIOS_VARINFO * (* ADIOS_INQ_VAR_FN) (ADIOS_GROUP *gp, const char * varname);
typedef ADIOS_VARINFO * (* ADIOS_INQ_VAR_BYID_FN) (ADIOS_GROUP *gp, int varid);
typedef int64_t (* ADIOS_READ_VAR_FN) (ADIOS_GROUP * gp, const char * varname, const uint64_t * start, const uint64_t * count, void * data);
typedef int64_t (* ADIOS_READ_LOCAL_VAR_FN) (ADIOS_GROUP * gp, const char * varname, int idx, const uint64_t * start, const uint64_t * count, void * data);
typedef int64_t (* ADIOS_READ_VAR_BYID_FN) (ADIOS_GROUP * gp, int varid, const uint64_t * start, const uint64_t * count, void * data);
typedef int (* ADIOS_GET_ATTR_FN) (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type, int * size, void ** data);
typedef int (* ADIOS_GET_ATTR_BYID_FN) (ADIOS_GROUP * gp, int attrid, enum ADIOS_DATATYPES * type, int * size, void ** data);
typedef void (* ADIOS_RESET_DIMENSION_ORDER_FN) ();

struct adios_read_hooks_struct
{
    ADIOS_INIT_FN           adios_init_fn;
    ADIOS_FINALIZE_FN       adios_finalize_fn;
    ADIOS_FOPEN_FN          adios_fopen_fn;
    ADIOS_FCLOSE_FN         adios_fclose_fn;
    ADIOS_GOPEN_FN          adios_gopen_fn;
    ADIOS_GOPEN_BYID_FN     adios_gopen_byid_fn;
    ADIOS_GCLOSE_FN         adios_gclose_fn;
    ADIOS_INQ_VAR_FN        adios_inq_var_fn;
    ADIOS_INQ_VAR_BYID_FN   adios_inq_var_byid_fn;
    ADIOS_READ_VAR_FN       adios_read_var_fn;
    ADIOS_READ_LOCAL_VAR_FN adios_read_local_var_fn;
    ADIOS_READ_VAR_BYID_FN  adios_read_var_byid_fn;
    ADIOS_GET_ATTR_FN       adios_get_attr_fn;
    ADIOS_GET_ATTR_BYID_FN  adios_get_attr_byid_fn;
    ADIOS_RESET_DIMENSION_ORDER_FN adios_reset_dimension_order_fn;
};

#undef FORWARD_DECLARE
#endif
