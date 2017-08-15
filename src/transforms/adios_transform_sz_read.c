#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "adios_logger.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef HAVE_SZ

#include "sz.h"

int adios_transform_sz_is_implemented (void) {return 1;}

int adios_transform_sz_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                 adios_transform_pg_read_request *pg_reqgroup)
{
    //log_debug("function: %s\n", __FUNCTION__);
    void *buf = malloc(pg_reqgroup->raw_var_length);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_sz_subrequest_completed(adios_transform_read_request *reqgroup,
                                                          adios_transform_pg_read_request *pg_reqgroup,
                                                          adios_transform_raw_read_request *completed_subreq)
{
    //log_debug("function: %s\n", __FUNCTION__);
    return NULL;
}



adios_datablock * adios_transform_sz_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                           adios_transform_pg_read_request *completed_pg_reqgroup)
{
    //log_debug("function: %s\n", __FUNCTION__);
    size_t raw_size = (size_t) completed_pg_reqgroup->raw_var_length;
    unsigned char *raw_buff = completed_pg_reqgroup->subreqs->data;

    // Get type info
    int dtype;
    switch (reqgroup->transinfo->orig_type)
    {
        case adios_double:
            dtype = SZ_DOUBLE;
            break;
        case adios_real:
            dtype = SZ_FLOAT;
            break;
        default:
            adios_error(err_transform_failure, "No supported data type\n");
            return NULL;
            break;
    }

    // Get dimension info
    int ndims = reqgroup->transinfo->orig_ndim;
    if (ndims > 5)
    {
        adios_error(err_transform_failure, "No more than 5 dimension is supported.\n");
        return NULL;
    }

    size_t r[5] = {0,0,0,0,0};
    int i = 0;
    for(i = 0; i < ndims; i++)
    {
        uint64_t dsize = (uint64_t)(completed_pg_reqgroup->orig_varblock->count[i]);
        r[ndims-i-1] = dsize;
    }

    void* orig_buff;
    orig_buff = SZ_decompress(dtype, raw_buff, raw_size, r[4], r[3], r[2], r[1], r[0]);

    if (dtype == SZ_FLOAT)
    {
        log_debug("%10s: %g %g %g %g %g ... \n", "out",
                  ((float*)orig_buff)[0], ((float*)orig_buff)[1], ((float*)orig_buff)[2], ((float*)orig_buff)[3], ((float*)orig_buff)[4]);
    }
    else if (dtype == SZ_DOUBLE)
    {
        log_debug("%10s: %g %g %g %g %g ... \n", "out",
                  ((double*)orig_buff)[0], ((double*)orig_buff)[1], ((double*)orig_buff)[2], ((double*)orig_buff)[3], ((double*)orig_buff)[4]);
    }
    //log_debug("=== SZ decompress ===\n");
    log_debug("%s: %d\n", "SZ dtype", dtype);
    log_debug("%s: %lu\n", "SZ raw_size", raw_size);
    /*
    log_debug("%s: %d %d %d %d %d ... %d %d %d %d %d\n", "SZ out_buff",
              raw_buff[0], raw_buff[1], raw_buff[2], raw_buff[3], raw_buff[4],
              raw_buff[raw_size-5], raw_buff[raw_size-4], raw_buff[raw_size-3], raw_buff[raw_size-2], raw_buff[raw_size-1]);
    int sum = 0;
    for (i=0; i<raw_size; i++)
    {
        sum += raw_buff[i];
    }
    log_debug("%s: %d\n", "SZ sum", sum);
     */
    log_debug("%s: %lu %lu %lu %lu %lu\n", "SZ dim", r[0], r[1], r[2], r[3], r[4]);
    //log_debug("=====================\n");

    return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, orig_buff);
}

// Do nothing for the full read request complete (typical)
adios_datablock * adios_transform_sz_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    //log_debug("function: %s\n", __FUNCTION__);
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(sz);

#endif
