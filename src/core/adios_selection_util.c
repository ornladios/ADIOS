/*
 * adios_selection_util.c
 *
 *  Created on: Jan 5, 2013
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "public/adios_error.h"
#include "public/adios_selection.h"
#include "adios_subvolume.h"
#include "adios_selection_util.h"

//
// NOTE: Intersection type guarantees:
// * The intersection of any two selections of the same type returns a third selection of that
//   same type
// * BB  + PTS  -> PTS
//

//
// One-on-one intersection functions
//
ADIOS_SELECTION * adios_selection_intersect_bb_bb(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb1,
                                                  const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb2) {
    const int ndim = bb1->ndim;
    uint64_t *new_start = malloc(ndim * sizeof(uint64_t));
    uint64_t *new_count = malloc(ndim * sizeof(uint64_t));

    assert(bb1->ndim == bb2->ndim);
    if (!new_start || !new_count) {
        adios_error(err_no_memory, "Cannot allocate memory for BOUNDINGBOX-BOUNDINGBOX selection intersection");
        return NULL;
    }

    if (intersect_bb(bb1, bb2, new_start, NULL, NULL, new_count)) {
        return common_read_selection_boundingbox(ndim, new_start, new_count);
    } else {
        free(new_start);
        free(new_count);
        return NULL;
    }
}

ADIOS_SELECTION * adios_selection_intersect_bb_pts(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb1,
                                                   const ADIOS_SELECTION_POINTS_STRUCT *pts2) {
    const int ndim = bb1->ndim;
    const uint64_t max_new_npts = pts2->npoints;

    uint64_t *new_pts = malloc(max_new_npts * ndim * sizeof(uint64_t));
    uint64_t i;
    int j;
    uint64_t *new_pts_ptr = new_pts;
    uint64_t *pts2_ptr;
    const uint64_t * const pts2_end_ptr = pts2->points + pts2->npoints * ndim;
    uint64_t new_npts = 0;

    assert(bb1->ndim == pts2->ndim);
    if (!new_pts) {
        adios_error(err_no_memory, "Cannot allocate memory for BOUNDINGBOX-POINTS selection intersection");
        return NULL;
    }

    // Check every pair of points for equality; whenever a shared point is found, output
    // it into the new point list
    for (pts2_ptr = pts2->points; pts2_ptr < pts2_end_ptr; pts2_ptr += ndim) {
        // Check each dimension component of the point for containment in the bounding box
        for (j = 0; j < ndim; j++)
            if (pts2_ptr[j] < bb1->start[j] ||
                pts2_ptr[j] >= bb1->start[j] + bb1->count[j])
                break;

        // Check whether any component was out of bounds; if so, skip this point; otherwise,
        // output the point
        if (j != ndim) {
            continue;
        } else {
            memcpy(new_pts_ptr, pts2_ptr, ndim * sizeof(uint64_t));
            new_pts_ptr += ndim;
            new_npts++;
        }
    }

    if (new_npts == 0) {
        free(new_pts);
        return NULL;
    } else {
        new_pts = (uint64_t*)realloc(new_pts, new_npts * ndim * sizeof(uint64_t));
        return common_read_selection_points(ndim, new_npts, new_pts);
    }
}

ADIOS_SELECTION * adios_selection_intersect_pts_pts(const ADIOS_SELECTION_POINTS_STRUCT *pts1,
                                                    const ADIOS_SELECTION_POINTS_STRUCT *pts2) {
    const int ndim = pts1->ndim;
    const uint64_t max_new_npts = pts1->npoints > pts2->npoints ? pts1->npoints : pts2->npoints;

    uint64_t *new_pts = malloc(max_new_npts * ndim * sizeof(uint64_t));
    uint64_t i, j;
    int k;
    uint64_t *new_pts_ptr = new_pts;
    uint64_t *pts1_ptr, *pts2_ptr;
    const uint64_t * const pts1_end_ptr = pts1->points + pts1->npoints * ndim;
    const uint64_t * const pts2_end_ptr = pts2->points + pts2->npoints * ndim;
    uint64_t new_npts = 0;

    assert(pts1->ndim == pts2->ndim);
    if (!new_pts) {
        adios_error(err_no_memory, "Cannot allocate memory for POINTS-POINTS selection intersection");
        return NULL;
    }

    // Check every pair of points for equality; whenever a shared point is found, output
    // it into the new point list
    for (pts1_ptr = pts1->points; pts1_ptr < pts1_end_ptr; pts1_ptr += ndim) {
        for (pts2_ptr = pts2->points; pts2_ptr < pts2_end_ptr; pts2_ptr += ndim) {
            // Check each dimension component of the pair of points for equality
            for (k = 0; k < ndim; k++)
                if (pts1_ptr[k] != pts2_ptr[k])
                    break;

            // Check whether any component was unequal; if so, skip this pair; otherwise,
            // output the shared point
            if (k != ndim) {
                continue;
            } else {
                memcpy(new_pts_ptr, pts1_ptr, ndim * sizeof(uint64_t));
                new_pts_ptr += ndim;
                new_npts++;
            }
        }
    }

    if (new_npts == 0) {
        free(new_pts);
        return NULL;
    } else {
        new_pts = (uint64_t*)realloc(new_pts, new_npts * sizeof(uint64_t));
        return common_read_selection_points(ndim, new_npts, new_pts);
    }
}

ADIOS_SELECTION * adios_selection_intersect_wb_wb(const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb1,
                                                  const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb2) {
    return wb1->index == wb2->index ?
           common_read_selection_writeblock(wb1->index):
           NULL;
}

//
// One-on-any intersection functions
//

// s2 can be selection type
inline static ADIOS_SELECTION * adios_selection_intersect_bb(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb1,
                                                             const ADIOS_SELECTION *s2) {
    switch (s2->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb2 = &s2->u.bb;
        return adios_selection_intersect_bb_bb(bb1, bb2);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *pts2 = &s2->u.points;
        return adios_selection_intersect_bb_pts(bb1, pts2);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Intersection of selection type BOUNDINGBOX and WRITEBLOCK not currently supported");
        return NULL;
    case ADIOS_SELECTION_AUTO:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Intersection of selection type BOUNDINGBOX and AUTO not currently supported");
        return NULL;
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", s2->type);
        return NULL;
    }
}

// s2 can be any selection type except boundingbox
inline static ADIOS_SELECTION * adios_selection_intersect_pts(const ADIOS_SELECTION_POINTS_STRUCT *pts1,
                                                              const ADIOS_SELECTION *s2) {
    switch (s2->type) {
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *pts2 = &s2->u.points;
        return adios_selection_intersect_pts_pts(pts1, pts2);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Intersection of selection type POINTS and WRITEBLOCK not currently supported");
        return NULL;
    case ADIOS_SELECTION_AUTO:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Intersection of selection type POINTS and AUTO not currently supported");
        return NULL;
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", s2->type);
        return NULL;
    }
}

// s2 can be any selection except boundingbox and points
inline static ADIOS_SELECTION * adios_selection_intersect_wb(const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb1,
                                                             const ADIOS_SELECTION *s2) {
    switch (s2->type) {
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb2 = &s2->u.block;
        return adios_selection_intersect_wb_wb(wb1, wb2);
    }
    case ADIOS_SELECTION_AUTO:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Intersection of selection type WRITEBLOCK and AUTO not currently supported");
        return NULL;
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", s2->type);
        return NULL;
    }
}

// s2 can only be auto
inline static ADIOS_SELECTION * adios_selection_intersect_auto(const ADIOS_SELECTION_AUTO_STRUCT *as1,
                                                               const ADIOS_SELECTION *s2) {
    return copy_selection(s2);
}

//
// Any-on-any intersection function
//

// The if statements impose a total order on the selection types, and call this function
// with arguments swapped if they are out of this order.
ADIOS_SELECTION * adios_selection_intersect(const ADIOS_SELECTION *s1, const ADIOS_SELECTION *s2) {
    switch (s1->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *bb1 = &s1->u.bb;
        return adios_selection_intersect_bb(bb1, s2);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *pts1 = &s1->u.points;
        if (s1->type == ADIOS_SELECTION_BOUNDINGBOX) {
            return adios_selection_intersect(s2, s1);
        } else {
            return adios_selection_intersect_pts(pts1, s2);
        }
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb1 = &s1->u.block;
        if (s1->type == ADIOS_SELECTION_BOUNDINGBOX ||
            s1->type == ADIOS_SELECTION_POINTS) {
            return adios_selection_intersect(s2, s1);
        } else {
            return adios_selection_intersect_wb(wb1, s2);
        }
    }
    case ADIOS_SELECTION_AUTO: {
        const ADIOS_SELECTION_AUTO_STRUCT *as1 = &s1->u.autosel;
        if (s1->type == ADIOS_SELECTION_BOUNDINGBOX ||
            s1->type == ADIOS_SELECTION_POINTS ||
            s1->type == ADIOS_SELECTION_WRITEBLOCK) {
            return adios_selection_intersect(s2, s1);
        } else {
            return adios_selection_intersect_auto(as1, s2);
        }
    }
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", s1->type);
        return NULL;
    }
}

