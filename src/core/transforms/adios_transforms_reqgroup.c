/*
 * adios_transforms_reqgroup.c
 *
 *  Created on: Jul 30, 2012
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_reqgroup.h"
#include "core/common_read.h"
#include "core/adios_subvolume.h"
#include "public/adios_selection.h"

// An adios_transform_read_reqgroup corresponds to a variable read request
// An adios_transform_pg_reqgroup corresponds to the portion of a read request intersecting a PG
// An adios_transform_read_subrequest corresponds to a single byte range read within a PG

// adios_transform_read_subrequest owns ->sel and ->data (both will be free'd)
// adios_transform_pg_reqgroup owns ->pg_selection and ->pg_intersection_to_global_copyspec (both will be free'd)
// adios_transform_read_reqgroup owns ->orig_sel, ->varinfo and ->transinfo (all will be free'd)

// Also, in all cases, ->transform_internal will be free'd if it is non-NULL
// Thus, it is best to keep only a single level malloc in transform_internal (no pointers to other mallocs)
// If this is not possible, you should free() malloced memory yourself as soon as it is no longer needed


#define MYFREE(p) {if (p) free(p); (p)=NULL;}

// Generic list manipulation
// Assumes the list node struct has a ->next field

#define LIST_APPEND(head, elem, elem_type) 		\
    if (!(head)) {								\
        (head) = (elem);						\
    } else {									\
        elem_type *_cur = (head);    	    	\
        while (_cur->next)						\
            _cur = _cur->next;					\
        _cur->next = (elem);					\
    }											\


#define LIST_REMOVE(head, elem, elem_type, removed) \
    LIST_REMOVE_PRED(head, (_cur == (elem)), elem_type, removed)

#define LIST_REMOVE_KEY(head, key, keyfield, elem_type, removed) \
    LIST_REMOVE_PRED(head, (_cur->keyfield == (key)), elem_type, removed)

#define LIST_REMOVE_PRED(head, pred, elem_type, removed)	\
    if (!(head)) {									\
        removed = 0;								\
    } else {										\
        elem_type *_prev = NULL;					\
        elem_type *_cur = (head);					\
        while (_cur) {								\
            if (pred)								\
                break;								\
            _prev = _cur;							\
            _cur = _cur->next;						\
        }											\
        if (!_cur) {								\
            removed = NULL;							\
        } else {									\
            if (!_prev) {							\
                (head) = (head)->next;				\
            } else {								\
                _prev->next = _cur->next;			\
            }										\
            _cur->next = NULL;						\
            removed = _cur;							\
        }											\
    }

static int common_adios_selection_equal(ADIOS_SELECTION *sel1, ADIOS_SELECTION *sel2) {
    if (sel1->type != sel2->type)
        return 0;

    switch (sel1->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb1 = &sel1->u.bb;
        ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb2 = &sel2->u.bb;
        return (bb1->ndim == bb2->ndim) &&
               memcmp(bb1->start, bb2->start, bb1->ndim * sizeof(uint64_t)) == 0 &&
               memcmp(bb1->count, bb2->count, bb1->ndim * sizeof(uint64_t)) == 0;
    }
    default:
        adios_error(err_operation_not_supported, "Selection types other than bounding box not supported in %s\n", __FUNCTION__);
        return 0;
    }
}

//
// adios_transform_read_subreq struct manipulation
//

adios_transform_read_subrequest * adios_transform_new_subreq(ADIOS_SELECTION *sel, void *data) {
    adios_transform_read_subrequest *new_subreq = malloc(sizeof(adios_transform_read_subrequest));
    new_subreq->sel = sel;
    new_subreq->data = data;
    new_subreq->completed = 0;
    new_subreq->transform_internal = 0;
    new_subreq->next = 0;
    return new_subreq;
}

adios_transform_read_subrequest * adios_transform_new_subreq_byte_segment(const ADIOS_VARBLOCK *raw_varblock, uint64_t start, uint64_t count, void *data) {
    ADIOS_SELECTION *sel;
    uint64_t *start_sel, *count_sel;

    // TODO: Move this bounding box construction to a separate function?
    start_sel = (uint64_t*)malloc(2 * sizeof(uint64_t));
    count_sel = (uint64_t*)malloc(2 * sizeof(uint64_t));
    // PG ID dim
    start_sel[0] = raw_varblock->start[0]; // PG ID
    count_sel[0] = 1;
    // Buffer dim
    start_sel[1] = start;
    count_sel[1] = count;

    // Transfer ownership of the our start/count vectors
    sel = common_read_selection_boundingbox(2, start_sel, count_sel);
    start_sel = count_sel = NULL;

    return adios_transform_new_subreq(sel, data);
}

adios_transform_read_subrequest * adios_transform_new_subreq_whole_pg(const ADIOS_VARBLOCK *raw_varblock, void *data) {
    return adios_transform_new_subreq_byte_segment(raw_varblock, 0, raw_varblock->count[1], data);
}

void adios_transform_subreq_mark_complete(adios_transform_read_reqgroup *parent_reqgroup, adios_transform_pg_reqgroup *parent_pg_reqgroup,
                                          adios_transform_read_subrequest *subreq) {
    if (subreq->completed)
        return;

    subreq->completed = 1;
    parent_pg_reqgroup->num_completed_subreqs++;

    if (parent_pg_reqgroup->num_completed_subreqs == parent_pg_reqgroup->num_subreqs) {
        parent_pg_reqgroup->completed = 1;
        parent_reqgroup->num_completed_pg_reqgroups++;

        if (parent_reqgroup->num_completed_pg_reqgroups == parent_reqgroup->num_pg_reqgroups) {
            parent_reqgroup->completed = 1;
        }
    }
}

// NOTE: MUST have removed the subrequest from the PG request group BEFORE calling this
void adios_transform_free_subreq(adios_transform_read_subrequest **subreq_ptr) {
    adios_transform_read_subrequest *subreq = *subreq_ptr;
    assert(!subreq->next); // Not a perfect check, but will catch many requests that are still linked

    // Free malloc'd resources
    common_read_selection_delete(subreq->sel);
    MYFREE(subreq->data);
    MYFREE(subreq->transform_internal);

    // Clear all data to 0's for safety
    memset(subreq, 0, sizeof(adios_transform_read_subrequest));

    // Free the entire struct and the user's pointer to it
    MYFREE(*subreq_ptr);
}

//
// adios_transform_pg_reqgroup struct manipulation
//

adios_transform_pg_reqgroup * adios_transform_new_pg_reqgroup(
        int timestep, int timestep_blockidx, int blockidx,
        const ADIOS_VARBLOCK *orig_varblock,
        const ADIOS_VARBLOCK *raw_varblock,
        const ADIOS_SELECTION *intersection_pg_rel,
        const ADIOS_SELECTION *intersection_orig_sel_rel,
        const ADIOS_SELECTION *intersection_global,
        const ADIOS_SELECTION *pg_bounds_global,
        adios_subvolume_copy_spec *pg_intersection_to_global_cs) {

    adios_transform_pg_reqgroup *new_pg_reqgroup;

    assert(orig_varblock); assert(pg_intersection_to_global_cs);
    assert(blockidx >= 0);

    new_pg_reqgroup = calloc(sizeof(adios_transform_pg_reqgroup), 1);
    new_pg_reqgroup->timestep = timestep;
    new_pg_reqgroup->timestep_blockidx = timestep_blockidx;
    new_pg_reqgroup->blockidx = blockidx;
    new_pg_reqgroup->raw_var_length = raw_varblock->count[1]; // TODO: Break out into helper function in transforms_common
    new_pg_reqgroup->raw_varblock = raw_varblock;
    new_pg_reqgroup->orig_varblock = orig_varblock;
    new_pg_reqgroup->intersection_pg_rel = intersection_pg_rel;
    new_pg_reqgroup->intersection_orig_sel_rel = intersection_orig_sel_rel;
    new_pg_reqgroup->intersection_global = intersection_global;
    new_pg_reqgroup->pg_bounds_global = pg_bounds_global;
    new_pg_reqgroup->pg_intersection_to_global_copyspec = pg_intersection_to_global_cs;
    // Other fields are 0'd

    return new_pg_reqgroup;
}

void adios_transform_pg_reqgroup_append_subreq(adios_transform_pg_reqgroup *pg_reqgroup, adios_transform_read_subrequest *subreq) {
    LIST_APPEND(pg_reqgroup->subreqs, subreq, adios_transform_read_subrequest);
    pg_reqgroup->num_subreqs++;
}

int adios_transform_pg_reqgroup_find_subreq(const adios_transform_pg_reqgroup *pg_reqgroup, const ADIOS_VARCHUNK *chunk,
                                            int skip_completed, adios_transform_read_subrequest **matching_subreq) {
    adios_transform_read_subrequest *cur;

    for (cur = pg_reqgroup->subreqs; cur; cur = cur->next) {
        if (skip_completed && cur->completed)
            continue;

        if (common_adios_selection_equal(cur->sel, chunk->sel))
            break;
    }

    // cur will be NULL if none matched
    *matching_subreq = cur;
    return cur != NULL;
}

int adios_transform_pg_reqgroup_remove_subreq(adios_transform_pg_reqgroup *pg_reqgroup, adios_transform_read_subrequest *subreq) {
    adios_transform_read_subrequest *removed;
    LIST_REMOVE(pg_reqgroup->subreqs, subreq, adios_transform_read_subrequest, removed);

    if (removed) pg_reqgroup->num_subreqs--;
    return removed != NULL;
}

adios_transform_read_subrequest * adios_transform_pg_reqgroup_pop_subreq(adios_transform_pg_reqgroup *pg_reqgroup) {
    adios_transform_read_subrequest *to_remove = pg_reqgroup->subreqs;
    if (adios_transform_pg_reqgroup_remove_subreq(pg_reqgroup, to_remove))
        return to_remove;
    else
        return NULL;
}

void adios_transform_free_pg_reqgroup(adios_transform_pg_reqgroup **pg_reqgroup_ptr) {
    adios_transform_pg_reqgroup *pg_reqgroup = *pg_reqgroup_ptr;
    adios_transform_read_subrequest *removed_subreq;

    assert(!pg_reqgroup->next);

    // Free any remaining subrequests
    while ((removed_subreq = adios_transform_pg_reqgroup_pop_subreq(pg_reqgroup)) != NULL) {
        adios_transform_free_subreq(&removed_subreq);
    }

    // Free malloc'd resources
    if (pg_reqgroup->intersection_pg_rel)
        common_read_selection_delete(pg_reqgroup->intersection_pg_rel);
    if (pg_reqgroup->intersection_orig_sel_rel)
        common_read_selection_delete(pg_reqgroup->intersection_orig_sel_rel);
    if (pg_reqgroup->intersection_global)
        common_read_selection_delete(pg_reqgroup->intersection_global);
    if (pg_reqgroup->pg_bounds_global)
        common_read_selection_delete(pg_reqgroup->pg_bounds_global);
    adios_copyspec_free(&pg_reqgroup->pg_intersection_to_global_copyspec, 1);
    MYFREE(pg_reqgroup->transform_internal);

    // Clear all data to 0's for safety
    memset(pg_reqgroup, 0, sizeof(adios_transform_pg_reqgroup));
    // Free the entire struct and the user's pointer to it
    MYFREE(*pg_reqgroup_ptr);
}

//
// adios_transform_read_reqgroup struct manipulation
//

adios_transform_read_reqgroup * adios_transform_new_read_reqgroup(
        const ADIOS_FILE *fp, const ADIOS_VARINFO *varinfo, const ADIOS_TRANSINFO *transinfo,
        const ADIOS_SELECTION *sel, int from_steps, int nsteps,
        void *data, enum ADIOS_FLAG swap_endianness) {

    adios_transform_read_reqgroup *new_reqgroup;
    assert(fp); assert(varinfo); assert(transinfo);
    assert(nsteps > 0);

    new_reqgroup = calloc(sizeof(adios_transform_read_reqgroup), 1);
    new_reqgroup->fp = fp;
    new_reqgroup->raw_varinfo = varinfo;
    new_reqgroup->transinfo = transinfo;

    new_reqgroup->from_steps = from_steps;
    new_reqgroup->nsteps = nsteps;
    new_reqgroup->orig_sel = copy_selection(sel);
    new_reqgroup->orig_data = data;
    new_reqgroup->swap_endianness = swap_endianness;

    new_reqgroup->orig_sel_timestep_size = compute_selection_size(sel) *
                                           common_read_type_size(transinfo->orig_type, NULL);

    // Other fields are 0'd

    return new_reqgroup;
}

void adios_transform_read_reqgroup_append_pg_reqgroup(adios_transform_read_reqgroup *reqgroup, adios_transform_pg_reqgroup *pg_reqgroup) {
    LIST_APPEND(reqgroup->pg_reqgroups, pg_reqgroup, adios_transform_pg_reqgroup);
    reqgroup->num_pg_reqgroups++;
}

int adios_transform_read_reqgroup_remove_pg_reqgroup(adios_transform_read_reqgroup *reqgroup, adios_transform_pg_reqgroup *pg_reqgroup) {
    adios_transform_pg_reqgroup *removed;
    LIST_REMOVE(reqgroup->pg_reqgroups, pg_reqgroup, adios_transform_pg_reqgroup, removed);

    if (removed) reqgroup->num_pg_reqgroups--;
    return removed != NULL;
}

adios_transform_pg_reqgroup * adios_transform_read_reqgroup_pop_pg_reqgroup(adios_transform_read_reqgroup *reqgroup) {
    adios_transform_pg_reqgroup *to_remove = reqgroup->pg_reqgroups;
    if (adios_transform_read_reqgroup_remove_pg_reqgroup(reqgroup, to_remove))
        return to_remove;
    else
        return NULL;
}

int adios_transform_read_reqgroup_find_subreq(const adios_transform_read_reqgroup *reqgroup, const ADIOS_VARCHUNK *chunk, int skip_completed,
                                              adios_transform_pg_reqgroup **matching_pg_reqgroup, adios_transform_read_subrequest **matching_subreq) {
    adios_transform_pg_reqgroup *cur;

    if (reqgroup->raw_varinfo->varid != chunk->varid)
        return 0;

    int found = 0;
    // Search all PG request groups
    for (cur = reqgroup->pg_reqgroups; cur; cur = cur->next) {
        // Skip completed PG reqgroups if required
        if (skip_completed && cur->completed)
            continue;

        // Skip PG reqgroups that are for other timesteps
        if (cur->timestep != chunk->from_steps)
            continue;

        // Delegate remaining search to the PG regroup
        found = adios_transform_pg_reqgroup_find_subreq(cur, chunk, skip_completed, matching_subreq);
        if (found)
            break;
    }

    // cur will be NULL if nothing matched
    *matching_pg_reqgroup = cur;
    return found;
}

int adios_transform_read_reqgroups_find_subreq(const adios_transform_read_reqgroup *reqgroup_head,
                                               const ADIOS_VARCHUNK *chunk, int skip_completed,
                                               adios_transform_read_reqgroup **matching_reqgroup,
                                               adios_transform_pg_reqgroup **matching_pg_reqgroup,
                                               adios_transform_read_subrequest **matching_subreq) {
    int found;
    adios_transform_read_reqgroup *cur;
    for (cur = reqgroup_head; cur; cur = cur->next) {
        found = adios_transform_read_reqgroup_find_subreq(cur, chunk, skip_completed, matching_pg_reqgroup, matching_subreq);
        if (found)
            break;
    }

    *matching_reqgroup = cur;
    return found;
}

void adios_transform_read_reqgroups_append(adios_transform_read_reqgroup **head, adios_transform_read_reqgroup *new_reqgroup) {
    LIST_APPEND(*head, new_reqgroup, adios_transform_read_reqgroup);
}

adios_transform_read_reqgroup * adios_transform_read_reqgroups_remove(adios_transform_read_reqgroup **head, adios_transform_read_reqgroup *reqgroup) {
    adios_transform_read_reqgroup *removed;
    LIST_REMOVE(*head, reqgroup, adios_transform_read_reqgroup, removed);
    return removed;
}

adios_transform_read_reqgroup * adios_transform_read_reqgroups_pop(adios_transform_read_reqgroup **head) {
    adios_transform_read_reqgroup *to_remove = *head;
    if (adios_transform_read_reqgroups_remove(head, to_remove))
        return to_remove;
    else
        return NULL;
}

void adios_transform_free_read_reqgroup(adios_transform_read_reqgroup **reqgroup_ptr) {
    adios_transform_read_reqgroup *reqgroup = *reqgroup_ptr;
    adios_transform_pg_reqgroup *removed_pg_reqgroup;

    assert(!reqgroup->next);

    // Free any remaining subrequests
    while ((removed_pg_reqgroup = adios_transform_read_reqgroup_pop_pg_reqgroup(reqgroup)) != NULL) {
        adios_transform_free_pg_reqgroup(&removed_pg_reqgroup);
    }

    // Free malloc'd resources

    // Free any data buffer lent to the user, but don't free the VARCHUNK; that
    // should have been done already by the user
    if (reqgroup->lent_varchunk)
        MYFREE(reqgroup->lent_varchunk->data);

    common_read_selection_delete((ADIOS_SELECTION*)reqgroup->orig_sel); // Remove const
    common_read_free_transinfo(reqgroup->raw_varinfo,
                               (ADIOS_TRANSINFO*)reqgroup->transinfo); // Remove const
    common_read_free_varinfo((ADIOS_VARINFO*)reqgroup->raw_varinfo); // Remove const
    MYFREE(reqgroup->transform_internal);

    // Clear all data to 0's for safety
    memset(reqgroup, 0, sizeof(adios_transform_read_reqgroup));
    // Free the entire struct and the user's pointer to it
    MYFREE(*reqgroup_ptr);
}


