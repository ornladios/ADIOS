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
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_reqgroup.h"
#include "core/adios_subvolume.h"

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(none);
DECLARE_TRANSFORM_READ_METHOD(identity);
DECLARE_TRANSFORM_READ_METHOD(alacrity);
DECLARE_TRANSFORM_READ_METHOD(compress);
DECLARE_TRANSFORM_READ_METHOD(mloc);

// Transform read method registry
adios_transform_read_method TRANSFORM_READ_METHODS[num_adios_transform_types];

void adios_transform_read_init() {
    static int adios_transforms_initialized = 0;
    if (adios_transforms_initialized)
        return;

    REGISTER_TRANSFORM_READ_METHOD(none, adios_transform_none);
    REGISTER_TRANSFORM_READ_METHOD(identity, adios_transform_identity);
    REGISTER_TRANSFORM_READ_METHOD(alacrity, adios_transform_alacrity);
    REGISTER_TRANSFORM_READ_METHOD(compress, adios_transform_compress);
    REGISTER_TRANSFORM_READ_METHOD(mloc, adios_transform_mloc);

    adios_transforms_initialized = 1;
}


// Datablock management

adios_datablock * adios_datablock_new(
        enum ADIOS_DATATYPES elem_type,
        int timestep,
        const ADIOS_SELECTION *bounds,
        void *data) {

    assert(bounds);
    assert(data);
    return adios_datablock_new_ragged_offset(elem_type, timestep, bounds, 0, data);
}

/*
 * Note: only valid for bounding box selections (since there are no ragged
 * arrays for point selections, and no other selection types are supported
 * right now).
 */
adios_datablock * adios_datablock_new_ragged(
        enum ADIOS_DATATYPES elem_type,
        int timestep,
        const ADIOS_SELECTION *bounds,
        const uint64_t *ragged_offsets, void *data) {

    assert(bounds);
    assert(data);
    assert(bounds->type == ADIOS_SELECTION_BOUNDINGBOX);

    const uint64_t ragged_offset = ragged_offsets ?
            compute_ragged_array_offset(bounds->u.bb.ndim, ragged_offsets, bounds->u.bb.count) :
            0;

    return adios_datablock_new_ragged_offset(elem_type, timestep, bounds, ragged_offset, data);
}

adios_datablock * adios_datablock_new_ragged_offset(
        enum ADIOS_DATATYPES elem_type,
        int timestep,
        const ADIOS_SELECTION *bounds,
        uint64_t ragged_offset, void *data) {

    assert(bounds);
    assert(data);

    adios_datablock *datablock = malloc(sizeof(adios_datablock));

    datablock->elem_type = elem_type;
    datablock->bounds = copy_selection(bounds);
    datablock->timestep = timestep;
    datablock->ragged_offset = ragged_offset;
    datablock->data = data;

    return datablock;
}

#define MYFREE(p) {if (p) free(p); (p)=NULL;}
void adios_datablock_free(adios_datablock **datablock_ptr, int free_data) {
    adios_datablock *datablock = *datablock_ptr;
    if (datablock) {
        if (datablock->bounds)
            common_read_selection_delete(datablock->bounds);
        if (free_data)
            MYFREE(datablock->data);
    }
    MYFREE(*datablock_ptr);
}
#undef MYFREE

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
    int blockidx, timestep, timestep_blockidx;
    int curblocks, start_blockidx, end_blockidx;
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
    for (blockidx = 0; blockidx < raw_varinfo->nsteps; blockidx++) {
        // Find the start block
        if (blockidx == from_steps) {
            start_blockidx = curblocks;
        }
        curblocks += raw_varinfo->nblocks[blockidx];
        // Find the end block, then stop
        if (blockidx == to_steps - 1) {
            end_blockidx = curblocks;
            break;
        }
    }

    // Retrieve blockinfos, if they haven't been done retrieved
    if (!raw_varinfo->blockinfo)
        common_read_inq_var_blockinfo_raw(fp, raw_varinfo);
    if (!transinfo->orig_blockinfo)
        common_read_inq_trans_blockinfo(fp, raw_varinfo, transinfo);

    // Assemble read requests for each varblock
    pg_intersection_to_global_copyspec = NULL;
    blockidx = start_blockidx;
    timestep = from_steps;
    timestep_blockidx = 0;
    while (blockidx != end_blockidx) { //for (blockidx = startblock_idx; blockidx != endblock_idx; blockidx++) {
        raw_vb = &raw_varinfo->blockinfo[blockidx];
        orig_vb = &transinfo->orig_blockinfo[blockidx];

        // Find the intersection, if any
        pg_intersection_to_global_copyspec = malloc(sizeof(adios_subvolume_copy_spec));
        intersects = adios_copyspec_init_from_bb_intersection(pg_intersection_to_global_copyspec, &sel->u.bb, orig_vb->count, orig_vb->start);

        if (intersects) {
            // Make a PG read request group, and fill it with some subrequests, and link it into the read reqgroup
            adios_transform_pg_reqgroup *new_pg_reqgroup;
            ADIOS_SELECTION *intersection_pg_rel;
            ADIOS_SELECTION *intersection_orig_sel_rel;
            ADIOS_SELECTION *intersection_global;
            ADIOS_SELECTION *pg_bounds_global;

            intersection_pg_rel = adios_copyspec_to_src_selection(pg_intersection_to_global_copyspec);
            intersection_orig_sel_rel = adios_copyspec_to_dst_selection(pg_intersection_to_global_copyspec);
            // Derelativize from PG space to global space
            intersection_global = new_derelativized_selection(intersection_pg_rel, orig_vb->start);
            pg_bounds_global = varblock_to_bb(transinfo->orig_ndim, orig_vb);

            // Transfer ownership of pg_intersection_to_global_copyspec
            new_pg_reqgroup = adios_transform_new_pg_reqgroup(timestep, timestep_blockidx,
                                                              blockidx,
                                                              orig_vb, raw_vb,
                                                              intersection_pg_rel,
                                                              intersection_orig_sel_rel,
                                                              intersection_global,
                                                              pg_bounds_global,
                                                              pg_intersection_to_global_copyspec);
            pg_intersection_to_global_copyspec = NULL;

            TRANSFORM_READ_METHODS[transinfo->transform_type].transform_generate_read_subrequests(new_reqgroup, new_pg_reqgroup);

            adios_transform_read_reqgroup_append_pg_reqgroup(new_reqgroup, new_pg_reqgroup);
        } else {
            adios_copyspec_free(&pg_intersection_to_global_copyspec, 1);
        }

        // Increment block indexes
        blockidx++;
        timestep_blockidx++;
        if (timestep_blockidx == raw_varinfo->nblocks[timestep]) {
            timestep_blockidx = 0;
            timestep++;
        }
    }
    assert(!pg_intersection_to_global_copyspec);

    return new_reqgroup;
}

adios_datablock * adios_transform_subrequest_completed(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *pg_reqgroup,
                                                       adios_transform_read_subrequest *completed_subreq) {
    enum ADIOS_TRANSFORM_TYPE transform_type = reqgroup->transinfo->transform_type;
    assert(is_transform_type_valid(transform_type));
    return TRANSFORM_READ_METHODS[transform_type].transform_subrequest_completed(reqgroup, pg_reqgroup, completed_subreq);
}

adios_datablock * adios_transform_pg_reqgroup_completed(adios_transform_read_reqgroup *reqgroup,
                                                        adios_transform_pg_reqgroup *completed_pg_reqgroup) {

    enum ADIOS_TRANSFORM_TYPE transform_type = reqgroup->transinfo->transform_type;
    assert(is_transform_type_valid(transform_type));
    return TRANSFORM_READ_METHODS[transform_type].transform_pg_reqgroup_completed(reqgroup, completed_pg_reqgroup);
}

adios_datablock * adios_transform_read_reqgroup_completed(adios_transform_read_reqgroup *completed_reqgroup) {
    enum ADIOS_TRANSFORM_TYPE transform_type = completed_reqgroup->transinfo->transform_type;
    assert(is_transform_type_valid(transform_type));
    return TRANSFORM_READ_METHODS[transform_type].transform_reqgroup_completed(completed_reqgroup);
}
