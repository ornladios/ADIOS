#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "core/transforms/adios_transforms_common.h"
#include "core/transforms/adios_transforms_read.h"
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/adios_logger.h"
#include "core/adios_subvolume.h"
#include "public/adios_error.h"

#ifdef ALACRITY
#include "alacrity.h"
#include <stdbool.h>

int adios_transform_alacrity_generate_read_subrequests(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *pg_reqgroup) {

    assert(reqgroup);
    assert(pg_reqgroup);

    //int read_layer_buf_allowed = (reqgroup->orig_data == NULL);
    void *buf = malloc(pg_reqgroup->raw_var_length);
    adios_transform_read_subrequest *subreq = adios_transform_new_subreq_whole_pg(pg_reqgroup->raw_varblock, buf);
    adios_transform_pg_reqgroup_append_subreq(pg_reqgroup, subreq);

    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_alacrity_subrequest_completed(
                    adios_transform_read_reqgroup *reqgroup,
                    adios_transform_pg_reqgroup *pg_reqgroup,
                    adios_transform_read_subrequest *completed_subreq) {
    return NULL;
}

static void * decodeBuffer(void **buf_ptr, _Bool free_buf) {
    void * buf = *buf_ptr;

    ALPartitionData part;
    uint64_t out_size;
    memstream_t ms = memstreamInitReturn(buf);

    ALDeserializePartitionData(&part, &ms);
    memstreamDestroy(&ms, false);
    if (free_buf) {
        free(buf);
        *buf_ptr = NULL;
    }

    // Decode the ALACRITY data
    void *decode_buf = malloc(ALGetDecodeLength(&part));
    ALDecode(&part, decode_buf, &out_size);
    ALPartitionDataDestroy(&part);

    return decode_buf;
}

adios_datablock * adios_transform_alacrity_pg_reqgroup_completed(
        adios_transform_read_reqgroup *reqgroup,
        adios_transform_pg_reqgroup *completed_pg_reqgroup) {

    // Decode the data (NOTE: this will free encoded the given data buffer)
    void *decoded_data = decodeBuffer(&completed_pg_reqgroup->subreqs->data, true);

    return adios_datablock_new(reqgroup->transinfo->orig_type,
                               completed_pg_reqgroup->timestep,
                               completed_pg_reqgroup->pg_bounds_global,
                               decoded_data);
}

adios_datablock * adios_transform_alacrity_reqgroup_completed(
        adios_transform_read_reqgroup *completed_reqgroup) {
    return NULL;
}

#else
DECLARE_TRANSFORM_READ_METHOD_UNIMPL(alacrity)
#endif
