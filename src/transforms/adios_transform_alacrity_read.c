#include <stdint.h>
#include <assert.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_hooks_read.h"
#include "adios_subvolume.h"
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
ADIOS_VARCHUNK * adios_transform_alacrity_subrequest_completed(
                    adios_transform_read_reqgroup *reqgroup,
                    adios_transform_pg_reqgroup *pg_reqgroup,
                    adios_transform_read_subrequest *completed_subreq,
                    enum ADIOS_READ_RESULT_MODE mode) {
    return NULL;
}

static void * decodeBuffer(void *buf, _Bool free_buf) {
    ALPartitionData part;
    uint64_t out_size;
    memstream_t ms = memstreamInitReturn(buf);

    ALDeserializePartitionData(&part, &ms);
    memstreamDestroy(&ms, free_buf);

    // Decode the ALACRITY data
    void *decode_buf = malloc(ALGetDecodeLength(&part));
    ALDecode(&part, decode_buf, &out_size);
    ALPartitionDataDestroy(&part);

    return decode_buf;
}

ADIOS_VARCHUNK * adios_transform_alacrity_pg_reqgroup_completed(
        adios_transform_read_reqgroup *reqgroup,
        adios_transform_pg_reqgroup *completed_pg_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) {

    void *encoded_data, *decoded_data;
    encoded_data = completed_pg_reqgroup->subreqs->data; // The first (and only) subrequest data buffer
    decoded_data = decodeBuffer(encoded_data, true);
    completed_pg_reqgroup->subreqs->data = NULL; // This was free'd in decodeBuffer

    // If we are allowed to return a partial result, return this PG's data
    // Else, copy the PG's data into the global result buffer, and do nothing
    //   all data is present
    if (mode == adios_read_return_partial) {
        ADIOS_VARCHUNK *retchunk = (ADIOS_VARCHUNK *)malloc(sizeof(ADIOS_VARCHUNK));
        retchunk->varid = reqgroup->raw_varinfo->varid;
        retchunk->type = reqgroup->transinfo->orig_type;
        retchunk->sel = adios_copyspec_to_dst_selection(completed_pg_reqgroup->pg_intersection_to_global_copyspec);
        retchunk->data = decoded_data;
        return retchunk;
    } else {
        copy_subvolume_with_spec(reqgroup->orig_data,					// Copy TO original buffer
                                 decoded_data,							// Copy FROM decoded buffer
                                 completed_pg_reqgroup->pg_intersection_to_global_copyspec,	// Copy USING the PG-to-global copy spec
                                 reqgroup->transinfo->orig_type,		// Copy elements of the original type
                                 reqgroup->swap_endianness);			// Swap endianness if needed
        return NULL;
    }
}

ADIOS_VARCHUNK * adios_transform_alacrity_reqgroup_completed(
        adios_transform_read_reqgroup *completed_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) {

    if (mode == adios_read_return_complete) {
        ADIOS_VARCHUNK *retchunk = (ADIOS_VARCHUNK *)malloc(sizeof(ADIOS_VARCHUNK));
        retchunk->varid = completed_reqgroup->raw_varinfo->varid;
        retchunk->type = completed_reqgroup->transinfo->orig_type;
        retchunk->sel = copy_selection(completed_reqgroup->orig_sel);
        retchunk->data = completed_reqgroup->orig_data;
        return retchunk;
    } else {
        // Mode is either noreturn, in which case there should be no return,
        // or return_partial, in which case all data was already returned
        return NULL;
    }
}

#else
DECLARE_TRANSFORM_READ_METHOD_UNIMPL(alacrity)
#endif
