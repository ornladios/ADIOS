#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "adios_logger.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef SZ

int adios_transform_sz_is_implemented (void) {return 1;}

int adios_transform_sz_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                       adios_transform_pg_read_request *pg_reqgroup)
{
    log_info("function: %s\n", __FUNCTION__);
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
    log_info("function: %s\n", __FUNCTION__);
    return NULL;
}



adios_datablock * adios_transform_sz_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                                 adios_transform_pg_read_request *completed_pg_reqgroup)
{
    log_info("function: %s\n", __FUNCTION__);
    uint64_t raw_size = (uint64_t)completed_pg_reqgroup->raw_var_length;
    void* raw_buff = completed_pg_reqgroup->subreqs->data;

    uint64_t orig_size = adios_get_type_size(reqgroup->transinfo->orig_type, "");
    int d = 0;
    for(d = 0; d < reqgroup->transinfo->orig_ndim; d++)
        orig_size *= (uint64_t)(completed_pg_reqgroup->orig_varblock->count[d]);

    void* orig_buff = malloc(orig_size);

    // Decompress into orig_buff

    return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, orig_buff);
}

// Do nothing for the full read request complete (typical)
adios_datablock * adios_transform_sz_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    log_info("function: %s\n", __FUNCTION__);
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(zlib);

#endif

