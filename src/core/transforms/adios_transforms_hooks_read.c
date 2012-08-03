/*
 * adios_transforms_hooks_read.c
 *
 *  Created on: Jul 24, 2012
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "util.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(none);
DECLARE_TRANSFORM_READ_METHOD(identity);
DECLARE_TRANSFORM_READ_METHOD(alacrity);

// Transform read method registry
adios_transform_read_method TRANSFORM_READ_METHODS[num_adios_transform_types];

void adios_transform_read_init() {
    static int adios_transforms_initialized = 0;
    if (adios_transforms_initialized)
        return;

    REGISTER_TRANSFORM_READ_METHOD(none, adios_transform_none);
    REGISTER_TRANSFORM_READ_METHOD(identity, adios_transform_identity);
    REGISTER_TRANSFORM_READ_METHOD(alacrity, adios_transform_alacrity);

    adios_transforms_initialized = 1;
}




// Delegate functions

static int get_system_endianness() {
    uint16_t word = 0x1234;
    return *(uint8_t*)(&word) == 0x12; // Returns 1 (big endian) iff the high byte comes first
}

// TODO: Implement
adios_transform_read_reqgroup * adios_transform_generate_read_reqgroup(const ADIOS_VARINFO *raw_varinfo, const ADIOS_TRANSINFO* transinfo, const ADIOS_FILE *fp,
                                                                       const ADIOS_SELECTION *sel, int from_steps, int nsteps, void *data) {
    // Declares
    adios_transform_read_reqgroup *new_reqgroup;
    int idx;
    int curblocks, startblock_idx, endblock_idx;
    int intersects;
    ADIOS_VARBLOCK *raw_vb, *orig_vb;
    adios_subvolume_copy_spec *pg_intersection_to_global_copyspec;

    enum ADIOS_FLAG swap_endianness = (fp->endianness == get_system_endianness()) ? adios_flag_no : adios_flag_yes;
    int to_steps = from_steps + nsteps;

    assert(is_transform_type_valid(transinfo->transform_type));
    assert(from_steps >= 0 && to_steps <= raw_varinfo->nsteps);

    if (sel->type != ADIOS_SELECTION_BOUNDINGBOX) {
        adios_error(err_operation_not_supported, "Reads of transformed variables using selection types other than bounding box are not supported.");
        assert(0);
    }

    // Allocate a new, empty request group
    new_reqgroup = adios_transform_new_read_reqgroup(fp, raw_varinfo, transinfo, sel, from_steps, nsteps, data, swap_endianness);

    // Find the block index for the start and end timestep
    curblocks = 0;
    for (idx = 0; idx < raw_varinfo->nsteps; idx++) {
        if (idx == from_steps) { startblock_idx = curblocks;		}	// Find the start block
        curblocks += raw_varinfo->nblocks[idx];
        if (idx == to_steps - 1) { endblock_idx = curblocks; break;	}	// Find the end block and stop
    }

    // Retrieve blockinfos, if they haven't been done retrieved
    if (!raw_varinfo->blockinfo)
        common_read_inq_var_blockinfo_raw(fp, raw_varinfo);
    if (!transinfo->orig_blockinfo)
        common_read_inq_trans_blockinfo(fp, raw_varinfo, transinfo);

    // Assemble read requests for each varblock
    pg_intersection_to_global_copyspec = NULL;
    for (idx = startblock_idx; idx != endblock_idx; idx++) {
        raw_vb = &raw_varinfo->blockinfo[idx];
        orig_vb = &transinfo->orig_blockinfo[idx];

        // Get a new copyspec if we need one
        if (!pg_intersection_to_global_copyspec)
            pg_intersection_to_global_copyspec = malloc(sizeof(adios_subvolume_copy_spec));

        // Find the intersection, if any
        intersects = adios_copyspec_init_from_selection_intersect(pg_intersection_to_global_copyspec, &sel->u.bb, orig_vb->count, orig_vb->start);
        if (intersects) {
            // Make a PG read request group, and fill it with some subrequests, and link it into the read reqgroup
            adios_transform_pg_reqgroup *new_pg_reqgroup;

            new_pg_reqgroup = adios_transform_new_pg_reqgroup(idx, orig_vb, raw_vb,
                                                              adios_copyspec_to_src_selection(pg_intersection_to_global_copyspec),
                                                              pg_intersection_to_global_copyspec);
            pg_intersection_to_global_copyspec = NULL;

            TRANSFORM_READ_METHODS[transinfo->transform_type].transform_generate_read_subrequests(new_reqgroup, new_pg_reqgroup);

            adios_transform_read_reqgroup_append_pg_reqgroup(new_reqgroup, new_pg_reqgroup);
        }
    }
    // If there is a leftover (uninitialized) copy spec, free it
    if (pg_intersection_to_global_copyspec)
        free(pg_intersection_to_global_copyspec);

    return new_reqgroup;
}

ADIOS_VARCHUNK * adios_transform_subrequest_completed(adios_transform_read_reqgroup *reqgroup,
                                                      adios_transform_pg_reqgroup *pg_reqgroup,
                                                      adios_transform_read_subrequest *completed_subreq,
                                                      enum ADIOS_READ_RESULT_MODE mode) {
    enum ADIOS_TRANSFORM_TYPE transform_type = reqgroup->transinfo->transform_type;
    assert(is_transform_type_valid(transform_type));
    return TRANSFORM_READ_METHODS[transform_type].transform_subrequest_completed(reqgroup, pg_reqgroup, completed_subreq, mode);
}

ADIOS_VARCHUNK * adios_transform_pg_reqgroup_completed(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *completed_pg_reqgroup,
                                                       enum ADIOS_READ_RESULT_MODE mode) {

    enum ADIOS_TRANSFORM_TYPE transform_type = reqgroup->transinfo->transform_type;
    assert(is_transform_type_valid(transform_type));
    return TRANSFORM_READ_METHODS[transform_type].transform_pg_reqgroup_completed(reqgroup, completed_pg_reqgroup, mode);
}

ADIOS_VARCHUNK * adios_transform_read_reqgroup_completed(adios_transform_read_reqgroup *completed_reqgroup,
                                                         enum ADIOS_READ_RESULT_MODE mode) {
    enum ADIOS_TRANSFORM_TYPE transform_type = completed_reqgroup->transinfo->transform_type;
    assert(is_transform_type_valid(transform_type));
    return TRANSFORM_READ_METHODS[transform_type].transform_reqgroup_completed(completed_reqgroup, mode);
}
