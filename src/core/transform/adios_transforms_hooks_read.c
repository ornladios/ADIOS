/*
 * adios_transforms_hooks_read.c
 *
 *  Created on: Jul 24, 2012
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
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
    adios_subvolume_copy_spec *pg_to_global_copyspec;

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
    pg_to_global_copyspec = NULL;
    for (idx = startblock_idx; idx != endblock_idx; idx++) {
        raw_vb = &raw_varinfo->blockinfo[idx];
        orig_vb = &transinfo->orig_blockinfo[idx];

        // Get a new copyspec if we need one
        if (!pg_to_global_copyspec)
            pg_to_global_copyspec = malloc(sizeof(adios_subvolume_copy_spec));

        // Find the intersection, if any
        intersects = adios_selection_to_copy_spec(pg_to_global_copyspec, &sel->u.bb, orig_vb->count, orig_vb->start);
        if (intersects) {
            // Make a PG read request group, and fill it with some subrequests, and link it into the read reqgroup
            adios_transform_pg_reqgroup *new_pg_reqgroup;
            new_pg_reqgroup = adios_transform_new_pg_reqgroup(idx, orig_vb, raw_vb,
                                                              adios_copy_spec_to_src_selection(pg_to_global_copyspec),
                                                              pg_to_global_copyspec);
            pg_to_global_copyspec = NULL;

            TRANSFORM_READ_METHODS[transinfo->transform_type].transform_generate_read_subrequests(new_reqgroup, new_pg_reqgroup);

            adios_transform_read_reqgroup_append_pg_reqgroup(new_reqgroup, new_pg_reqgroup);
        }
    }
    // If there is a leftover (uninitialized) copy spec, free it
    if (pg_to_global_copyspec)
        free(pg_to_global_copyspec);

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

enum ADIOS_ERRCODES adios_transform_retrieve_subvolume(
        enum ADIOS_TRANSFORM_TYPE transform_type,
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *subv_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness) {

    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_READ_METHODS[transform_type].transform_retrieve_subvolume(var, global_out, time_index, subv_spec, read_state, read_delegate, swap_endianness);
}


// Implementation of the "identity" transform, which does nothing, but
// exercises the transform framework for testing.
// (also in adios_transforms_hooks_write.c)

enum ADIOS_ERRCODES adios_transform_identity_retrieve_subvolume(
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *copy_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness) {

    assert(var->characteristics[time_index].transform.transform_type == adios_transform_identity);

    // Read the variable in its entirety
    void *buf = read_delegate(var, time_index, 0, adios_transform_var_get_transformed_size(var, time_index), read_state);

    copy_subvolume_with_spec(global_out, buf, copy_spec, adios_transform_get_var_original_type(var), swap_endianness);

    return err_no_error;
}

int adios_transform_identity_generate_read_subrequests(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *pg_reqgroup) {

    assert(reqgroup);
    assert(pg_reqgroup);
    // TODO: Optimize by only reading some of the PG where possible
    void *buf = malloc(pg_reqgroup->raw_var_length);
    adios_transform_read_subrequest *subreq = adios_transform_new_subreq_whole_pg(pg_reqgroup->raw_varblock, buf);
    adios_transform_pg_reqgroup_append_subreq(pg_reqgroup, subreq);

    return 0;
}

// Do nothing for individual subrequest
ADIOS_VARCHUNK * adios_transform_identity_subrequest_completed(
                    adios_transform_read_reqgroup *reqgroup,
                    adios_transform_pg_reqgroup *pg_reqgroup,
                    adios_transform_read_subrequest *completed_subreq,
                    enum ADIOS_READ_RESULT_MODE mode) {
    return NULL;
}

ADIOS_VARCHUNK * adios_transform_identity_pg_reqgroup_completed(
        adios_transform_read_reqgroup *reqgroup,
        adios_transform_pg_reqgroup *completed_pg_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) {

    // If we are allowed to return a partial result, return this PG's data
    // Else, copy the PG's data into the global result buffer, and do nothing
    //   all data is present
    if (mode == adios_read_return_partial) {
        ADIOS_VARCHUNK *retchunk = (ADIOS_VARCHUNK *)malloc(sizeof(ADIOS_VARCHUNK));
        retchunk->varid = reqgroup->raw_varinfo->varid;
        retchunk->type = reqgroup->transinfo->orig_type;
        retchunk->sel = adios_copy_spec_to_dst_selection(completed_pg_reqgroup->pg_to_global_copyspec);
        retchunk->data = completed_pg_reqgroup->subreqs->data; // The first (and only) subrequest data buffer
        return retchunk;
    } else {
        copy_subvolume_with_spec(reqgroup->orig_data,					// Copy TO original buffer
                                 completed_pg_reqgroup->subreqs->data,	// Copy FROM buffer of first (and only) subrequest
                                 completed_pg_reqgroup->pg_to_global_copyspec,	// Copy USING the PG-to-global copy spec
                                 reqgroup->transinfo->orig_type,		// Copy elements of the original type
                                 reqgroup->swap_endianness);			// Swap endianness if needed
        return NULL;
    }
}

ADIOS_VARCHUNK * adios_transform_identity_reqgroup_completed(
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
