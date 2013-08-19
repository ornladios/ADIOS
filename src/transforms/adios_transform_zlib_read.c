#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "adios_logger.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef ZLIB

#include "zlib.h"

int adios_transform_zlib_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                    adios_transform_pg_read_request *pg_reqgroup)
{
    void *buf = malloc(pg_reqgroup->raw_var_length);
    if (!buf) {
        log_error("Out of memory during read allocating %llu bytes for variable ID %d under zlib transform\n", pg_reqgroup->raw_var_length, reqgroup->raw_varinfo->varid);
        return 0;
    }

    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 1;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_zlib_subrequest_completed(adios_transform_read_request *reqgroup,
                                                            adios_transform_pg_read_request *pg_reqgroup,
                                                            adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}



adios_datablock * adios_transform_zlib_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                             adios_transform_pg_read_request *completed_pg_reqgroup)
{
    uint64_t compressed_size = (uint64_t)completed_pg_reqgroup->raw_var_length;
    void *compressed_data = completed_pg_reqgroup->subreqs->data;

    // Read whether it was compressed
    char is_compressed = *(char*)reqgroup->transinfo->transform_metadata;

    // Compute original variable size
    uint64_t uncompressed_size = 1;
    int i;
    for(i = 0; i < reqgroup->transinfo->orig_ndim; i++)
        uncompressed_size *= completed_pg_reqgroup->orig_varblock->count[i];

    void *uncompressed_data = NULL;
    if (is_compressed) { // Data is compressed
        // Allocate the output buffer
        uncompressed_data = malloc(uncompressed_size);
        if (!uncompressed_data) {
            log_warn("Out of memory allocating read buffer of %llu bytes for zlib transform\n", uncompressed_size);
            return NULL;
        }

        int rtn = uncompress(uncompressed_data, &uncompressed_size, compressed_data, compressed_size);
        if (rtn != Z_OK)
            return NULL;
    } else { // Data is not compressed, so just copy it
        uncompressed_data = compressed_data;
        completed_pg_reqgroup->subreqs->data = NULL; // Give up ownership of the buffer
    }

    // Return the buffer as a datablock
    return adios_datablock_new(reqgroup->transinfo->orig_type,
                               completed_pg_reqgroup->timestep,
                               completed_pg_reqgroup->pg_bounds_sel,
                               uncompressed_data);
}

adios_datablock * adios_transform_zlib_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(zlib);

#endif

