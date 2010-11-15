/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdlib.h>
#include <string.h>
#include "adios.h"
#include "bp_utils.h"
#include "bp_types.h"
#include "adios_read.h"
#include "common_read.h"
#include "adios_error.h"
#include "futils.h"
#define BYTE_ALIGN 8

#ifdef DMALLOC
#include "dmalloc.h"
#endif

const char *adios_errmsg ()
{
    return adios_get_last_errmsg();
}

int adios_set_read_method (enum ADIOS_READ_METHOD method)
{
    return common_read_set_read_method (method);
} 

int adios_read_init(MPI_Comm comm)
{
    return common_read_init (comm);
}

int adios_read_finalize()
{
    return common_read_finalize();
}

ADIOS_FILE * adios_fopen (const char * fname, MPI_Comm comm)
{
    return common_read_fopen (fname, comm);
}

int adios_fclose (ADIOS_FILE *fp) 
{
    return common_read_fclose (fp);
}


ADIOS_GROUP * adios_gopen (ADIOS_FILE *fp, const char * grpname)
{
    return common_read_gopen (fp, grpname);
}

ADIOS_GROUP * adios_gopen_byid (ADIOS_FILE *fp, int grpid)
{
    return common_read_gopen_byid (fp, grpid);
}
                   
int adios_gclose (ADIOS_GROUP *gp)
{
    return common_read_gclose (gp);
}


int adios_get_attr (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type,
                    int * size, void ** data)
{
    return common_read_get_attr (gp, attrname, type, size, data);
}

int adios_get_attr_byid (ADIOS_GROUP * gp, int attrid, 
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    return common_read_get_attr_byid (gp, attrid, type, size, data);
}


ADIOS_VARINFO * adios_inq_var (ADIOS_GROUP *gp, const char * varname) 
{
    return common_read_inq_var (gp, varname);
}

ADIOS_VARINFO * adios_inq_var_byid (ADIOS_GROUP *gp, int varid)
{
    return common_read_inq_var_byid (gp, varid);
}

void adios_free_varinfo (ADIOS_VARINFO *vp)
{
    common_read_free_varinfo (vp);
}

int64_t adios_read_var (ADIOS_GROUP * gp, const char * varname,
                        const uint64_t * start, const uint64_t * count,
                        void * data)
{
    return common_read_read_var (gp, varname, start, count, data);
}

int64_t adios_read_local_var (ADIOS_GROUP    * gp,
                              const char     * varname,
                              int            idx,
                              const uint64_t * start,
                              const uint64_t * count,
                              void           * data)
{
    return common_read_read_local_var (gp, varname, idx, start, count, data);
}

int64_t adios_read_var_byid (ADIOS_GROUP    * gp,
                             int              varid,
                             const uint64_t  * start,
                             const uint64_t  * count,
                             void           * data)
{
    return common_read_read_var_byid (gp, varid, start, count, data);
}

const char * adios_type_to_string (enum ADIOS_DATATYPES type)
{
    return common_read_type_to_string (type);
}

int adios_type_size(enum ADIOS_DATATYPES type, void *data)
{
    return common_read_type_size(type, data);
}


void adios_print_groupinfo (ADIOS_GROUP *gp) 
{
    common_read_print_groupinfo(gp);
}


void adios_print_fileinfo (ADIOS_FILE *fp) 
{
    common_read_print_fileinfo(fp);
}

