#include <assert.h>

#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "public/adios_error.h"
#include "public/adios_types.h"

#include "transforms/adios_transforms_common.h"
#include "transforms/adios_transforms_datablock.h"
#include "transforms/adios_transforms_hooks_read.h"
#include "transforms/adios_transforms_read.h"
#include "transforms/adios_transforms_util.h"

#define MYFREE(p) {free(p); (p)=NULL;}

enum ADIOS_TRANSFORM_REQGROUP_RESULT_MODE adios_transform_reqgroup_get_result_mode(adios_transform_read_reqgroup *reqgroup) {
    return reqgroup->orig_data != NULL ? FULL_RESULT_MODE : PARTIAL_RESULT_MODE;
}

// Delegate functions

static int get_system_endianness() {
    uint16_t word = 0x1234;
    return *(uint8_t*)(&word) == 0x12; // Returns 1 (big endian) iff the high byte comes first
}

/*
 * Determines the block indices corresponding to a start and end timestep.
 * Both the input start/end timesteps and the output start/end blockidx are lower bound inclusive, upper bound exclusive: [start, end)
 */
static void compute_blockidx_range(const ADIOS_VARINFO *raw_varinfo, int from_steps, int to_steps, int *start_blockidx, int *end_blockidx) {
    int blockidx;

    // Find the block index for the start and end timestep
    int curblocks = 0;
    for (blockidx = 0; blockidx < raw_varinfo->nsteps; blockidx++) {
        // Find the start block
        if (blockidx == from_steps) {
            *start_blockidx = curblocks;
        }
        curblocks += raw_varinfo->nblocks[blockidx];
        // Find the end block, then stop
        if (blockidx == to_steps - 1) {
            *end_blockidx = curblocks;
            break;
        }
    }
}

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

    // Precondition checking
    assert(is_transform_type_valid(transinfo->transform_type));
    assert(from_steps >= 0 && to_steps <= raw_varinfo->nsteps);

    if (sel->type != ADIOS_SELECTION_BOUNDINGBOX) {
        adios_error(err_operation_not_supported, "Reads of transformed variables using selection types other than bounding box are not supported.");
        assert(0);
    }

    // Compute the blockidx range, given the timesteps
    compute_blockidx_range(raw_varinfo, from_steps, to_steps, &start_blockidx, &end_blockidx);

    // Retrieve blockinfos, if they haven't been done retrieved
    if (!raw_varinfo->blockinfo)
        common_read_inq_var_blockinfo_raw(fp, raw_varinfo);
    if (!transinfo->orig_blockinfo)
        common_read_inq_trans_blockinfo(fp, raw_varinfo, transinfo);

    // Allocate a new, empty request group
    new_reqgroup = adios_transform_new_read_reqgroup(fp, raw_varinfo, transinfo, sel, from_steps, nsteps, data, swap_endianness);

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

            adios_transform_generate_read_subrequests(new_reqgroup, new_pg_reqgroup);

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

/*
 * Called whenever a subreq has been served by the read layer. Marks
 * all subreqs, pg_reqgroups and read_reqgroups as completed as necessary,
 * calls the appropriate hooks in the transform method, and returns an
 * adios_datablock if the transform method produces one.
 */
static adios_datablock * finish_subreq(
        adios_transform_read_reqgroup *reqgroup,
        adios_transform_pg_reqgroup *pg_reqgroup,
        adios_transform_read_subrequest *subreq) {

    adios_datablock *result, *tmp_result;

    // Mark the subrequest as complete
    assert(!subreq->completed && !pg_reqgroup->completed && !reqgroup->completed);
    adios_transform_subreq_mark_complete(reqgroup, pg_reqgroup, subreq);

    // Invoke all callbacks, depending on what completed, and
    // get at most one ADIOS_VARCHUNK to return
    result = adios_transform_subrequest_completed(reqgroup, pg_reqgroup, subreq);

    if (pg_reqgroup->completed) {
        tmp_result = adios_transform_pg_reqgroup_completed(reqgroup, pg_reqgroup);
        if (tmp_result) {
            assert(!result); // pg_reqgroup_completed returned a result, but subrequest_completed did as well
            result = tmp_result;
        }
    }

    if (reqgroup->completed) {
        tmp_result = adios_transform_read_reqgroup_completed(reqgroup);
        if (tmp_result) {
            assert(!result); // read_reqgroup_completed returned a result, but subrequest_completed or pg_reqgroup_completed did as well
            result = tmp_result;
        }
    }

    return result;
}

/*
 * Takes a datablock and applies its data to the user buffer for the given
 * read request group, then frees the given datablock. Assumes there is, in
 * fact, a user buffer (i.e., it is not NULL).
 *
 * Assumes that the datablock selection is of type bounding box.
 *
 * NOTE: also frees the data buffer within the datablock
 *
 * @return non-zero if some data in the datablock intersected the read
 *         request's selection, and was applied; returns 0 otherwise.
 */
static int apply_datablock_to_result_and_free(adios_datablock *datablock,
                                              adios_transform_read_reqgroup *reqgroup) {
    adios_subvolume_copy_spec *copyspec = malloc(sizeof(adios_subvolume_copy_spec));

    assert(datablock); assert(reqgroup);
    assert(reqgroup->orig_sel);
    assert(reqgroup->orig_data);
    assert(reqgroup->orig_sel->type == ADIOS_SELECTION_BOUNDINGBOX);
    assert(datablock->bounds->type == ADIOS_SELECTION_BOUNDINGBOX);
    assert(reqgroup->orig_sel->u.bb.ndim == datablock->bounds->u.bb.ndim);

    const int intersects =
        adios_copyspec_init_from_2bb_intersection(copyspec,
                                                  &reqgroup->orig_sel->u.bb,
                                                  &datablock->bounds->u.bb);
    if (intersects) {
        int rel_timestep = datablock->timestep - reqgroup->from_steps;
        void *timestep_data_slice =
                (char*)reqgroup->orig_data +
                rel_timestep * reqgroup->orig_sel_timestep_size;

        copy_subvolume_ragged_offset_with_spec(
                timestep_data_slice, datablock->data, copyspec,
                0, datablock->ragged_offset,
                datablock->elem_type, reqgroup->swap_endianness);
    }

    adios_copyspec_free(&copyspec, 1);
    adios_datablock_free(&datablock, 1);
    return intersects;
}

/*
 * Takes a datablock containing data potentially applicable to the given read
 * request group, identifies that data (if any), and returns it as an
 * ADIOS_VARCHUNK. Additionally, free the datablock.
 *
 * NOTE: This function transfers ownership of the ->data field of the datablock
 * and places it in the returned ADIOS_VARCHUNK (if an ADIOS_VARCHUNK is
 * returned).
 */
static ADIOS_VARCHUNK * apply_datablock_to_chunk_and_free(adios_datablock *result, adios_transform_read_reqgroup *reqgroup) {
    ADIOS_VARCHUNK *chunk;
    uint64_t *inter_goffset;
    uint64_t *inter_offset_within_result;
    uint64_t *inter_dims;

    assert(result); assert(reqgroup);
    assert(reqgroup->orig_sel);
    assert(reqgroup->orig_sel->type == ADIOS_SELECTION_BOUNDINGBOX);
    assert(reqgroup->orig_sel->u.bb.ndim == result->bounds->u.bb.ndim);
    assert(result->bounds->type == ADIOS_SELECTION_BOUNDINGBOX);

    const int ndim = result->bounds->u.bb.ndim;
    const int dimsize = ndim * sizeof(uint64_t);

    inter_goffset = malloc(dimsize);
    inter_offset_within_result = malloc(dimsize);
    inter_dims = malloc(dimsize);

    const int intersects =
            intersect_bb(&result->bounds->u.bb, &reqgroup->orig_sel->u.bb,
                         inter_goffset, inter_offset_within_result, NULL, inter_dims);

    if (intersects) {
        chunk = malloc(sizeof(ADIOS_VARCHUNK));

        // Compact the data within the datablock buffer so it is fully
        // contiguous according to the dimensions we will return to the user
        compact_subvolume_ragged_offset(result->data, ndim, inter_dims,
                                        result->bounds->u.bb.count, result->ragged_offset,
                                        inter_offset_within_result, result->elem_type);

        // Populate the chunk struct
        chunk->varid = reqgroup->raw_varinfo->varid;
        chunk->type = result->elem_type;

        // Transfer ownership of the data buffer
        chunk->data = result->data;
        result->data = NULL;

        // Transfer ownership of our global offset and dimension vectors
        chunk->sel = common_read_selection_boundingbox(ndim, inter_goffset, inter_dims);
        inter_goffset = NULL;
        inter_dims = NULL;
    } else {
        chunk = NULL;
    }

    MYFREE(inter_goffset);
    MYFREE(inter_offset_within_result);
    MYFREE(inter_dims);

    // Do free the data buffer; we removed that buffer from the datablock above if it was needed
    adios_datablock_free(&result, 1);
    return chunk;
}

static ADIOS_VARCHUNK * extract_chunk_from_finished_read_reqgroup(adios_transform_read_reqgroup *reqgroup) {
    assert(reqgroup);
    assert(reqgroup->completed);

    ADIOS_VARCHUNK *chunk = malloc(sizeof(ADIOS_VARCHUNK));
    chunk->varid = reqgroup->raw_varinfo->varid;
    chunk->type = reqgroup->transinfo->orig_type;

    // Transfer ownership of orig_data
    chunk->data = reqgroup->orig_data;
    reqgroup->orig_data = NULL;

    // Transfer ownership of orig_sel
    chunk->sel = (ADIOS_SELECTION*)reqgroup->orig_sel; // Remove const
    reqgroup->orig_sel = NULL;

    return chunk;
}

// Take an ADIOS_VARCHUNK that was just read and process it with the transform
// system. If it was part of a read request corresponding to a transformed
// variable, consume it, and possibly replacing it with a detransformed chunk.
// Otherwise, do nothing.
void adios_transform_process_read_chunk(adios_transform_read_reqgroup **reqgroups_head, ADIOS_VARCHUNK ** chunk) {
    adios_transform_read_reqgroup *reqgroup;
    adios_transform_pg_reqgroup *pg_reqgroup;
    adios_transform_read_subrequest *subreq;
    adios_datablock *result, *tmp_result;

    // Find the subrequest that matches the VARCHUNK that was just read (if any)
    int found = adios_transform_read_reqgroups_find_subreq(*reqgroups_head, *chunk, 1, &reqgroup, &pg_reqgroup, &subreq);

    // If no subrequest matches the VARCHUNK, it must correspond to a non-transformed variable.
    // In this case, return immediately and let it be processed as-is.
    if (!found)
        return;

    // Otherwise, this VARCHUNK corresponds to a subrequest.
    // Therefore, consume it, and perhaps replace it with a detransformed chunk.

    // Consume the chunk, as it will be passed to a transform method and should
    // not be processed by the caller.
    // (NOTE: Freeing this does not free the memory it points to)
    common_read_free_chunk(*chunk);
    *chunk = NULL;

    // Next, free any buffers held by the last-returned VARCHUNK, as they are now invalidated
    // by the user's call to check_reads (which in turn is the caller of this function)
    if (reqgroup->lent_varchunk && reqgroup->lent_varchunk->data)
        free(reqgroup->lent_varchunk->data);

    // Next, update the subreq that corresponds to this VARCHUNK as completed, retrieving any
    // produced result
    result = finish_subreq(reqgroup, pg_reqgroup, subreq);

    // Now, if a new adios_datablock is now available as a result of the above completed subreq,
    // apply it as a result for the user in a way appropriate to the current result mode
    if (result) {
        // Then, return data as appropriate depending on the return mode of this read operation
        //   PARTIAL: no user-allocated buffer is given for the full result, so results must be
        //            returned one VARCHUNK at a time.
        //   FULL: the user has supplied a buffer for full results, so patch relevant data from
        //         the returned VARCHUNK into this buffer.
        enum ADIOS_TRANSFORM_REQGROUP_RESULT_MODE result_mode = adios_transform_reqgroup_get_result_mode(reqgroup);
        switch (result_mode) {
        case PARTIAL_RESULT_MODE:
            // Apply this VARCHUNK
            *chunk = apply_datablock_to_chunk_and_free(result, reqgroup);

            reqgroup->lent_varchunk = *chunk;
            break;
        case FULL_RESULT_MODE:
            apply_datablock_to_result_and_free(result, reqgroup);

            // If the whole variable is now ready, return it as a VARCHUNK
            // Otherwise, return no chunk (NULL)
            if (reqgroup->completed) {
                *chunk = extract_chunk_from_finished_read_reqgroup(reqgroup);
            } else {
                assert(!*chunk); // No chunk to return, and *chunk is already NULL
            }
            break;
        }
    } else {
        assert(!*chunk); // No chunk to return, and *chunk is already NULL
    }

    // Free the read request group if it was completed
    if (reqgroup->completed) {
        adios_transform_read_reqgroups_remove(reqgroups_head, reqgroup);
        adios_transform_free_read_reqgroup(&reqgroup);
    }
}

/*
 * Process all read reqgroups, assuming they have been fully completed,
 * producing all required results based on the raw data read.
 * (This function is called after a blocking perform_reads completes)
 */
void adios_transform_process_all_reads(adios_transform_read_reqgroup **reqgroups_head) {
    // Mark all subrequests, PG request groups and read request groups
    // as completed, calling callbacks as needed
    adios_transform_read_reqgroup *reqgroup;
    adios_transform_pg_reqgroup *pg_reqgroup;
    adios_transform_read_subrequest *subreq;
    adios_datablock *result;

    // Complete each read reqgroup in turn
    while ((reqgroup = adios_transform_read_reqgroups_pop(reqgroups_head)) != NULL) {
        // Free leftover read request groups immediately, with no further processing
        if (reqgroup->completed) {
            adios_transform_free_read_reqgroup(&reqgroup);
            continue;
        }

        // Complete every child PG reqgroup
        for (pg_reqgroup = reqgroup->pg_reqgroups; pg_reqgroup; pg_reqgroup = pg_reqgroup->next) {
            // Skip completed PG reqgroups
            if (pg_reqgroup->completed) continue;

            // Complete every child subreq
            for (subreq = pg_reqgroup->subreqs; subreq; subreq = subreq->next) {
                // Skip completed subreqs
                if (subreq->completed) continue;

                // Mark the subreq as completed
                adios_transform_subreq_mark_complete(reqgroup, pg_reqgroup, subreq);
                assert(subreq->completed);

                // Make the required call to the transform method to apply the results
                result = adios_transform_subrequest_completed(reqgroup, pg_reqgroup, subreq);
                if (result) apply_datablock_to_result_and_free(result, reqgroup);
            }
            assert(pg_reqgroup->completed);

            // Make the required call to the transform method to apply the results
            result = adios_transform_pg_reqgroup_completed(reqgroup, pg_reqgroup);
            if (result) apply_datablock_to_result_and_free(result, reqgroup);
        }
        assert(reqgroup->completed);

        // Make the required call to the transform method to apply the results
        result = adios_transform_read_reqgroup_completed(reqgroup);
        if (result) apply_datablock_to_result_and_free(result, reqgroup);

        // Now that the read reqgroup has been processed, free it (which also frees all children)
        adios_transform_free_read_reqgroup(&reqgroup);
    }
}
