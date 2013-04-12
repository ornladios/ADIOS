/*
 * adios_patchdata.c
 *
 * Provides a main function, adios_patch_data, which copies all relevant data
 * from one buffer/selection to another buffer/selection. It supports various
 * combinations of source and destination selection types. In order to
 * minimize implementation, though, it classifies the current selection types
 * as follows:
 *
 * > Global geometric:
 *   > Bounding box
 *   > Points
 * > Local PG
 *   > Writeblock
 * > Other
 *   > Auto
 *
 * Patching is only supported between two selections within the same class.
 *
 *  Created on: Jan 15, 2013
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include "public/adios_error.h"
#include "public/adios_selection.h"
#include "adios_subvolume.h"
#include "adios_selection_util.h"
#include "adios_patchdata.h"

#define PATCH_UNIMPL(dsttype,srctype) \
    adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, \
                        "Patching of data from '%s' selection to '%s' selection not currently supported", \
                        srctype, dsttype);

// One-to-one patch functions
static uint64_t adios_patch_data_bb_to_bb(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb,
                                          void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb,
                                          enum ADIOS_DATATYPES datum_type,
                                          enum ADIOS_FLAG swap_endianness) {

    const int ndim = dst_bb->ndim;
    const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *inter_bb;
    uint64_t *inter_off_relative_to_dst;
    uint64_t *inter_off_relative_to_src;
    uint64_t volume;

    // Intersect the two bounding boxes
    const ADIOS_SELECTION *inter_sel = adios_selection_intersect_bb_bb(dst_bb, src_bb);

    // If there is no intersection, stop now, nothing to do
    if (!inter_sel)
        return 0;

    // (this is guaranteed by the selection intersection code; this is just to check for bugs)
    assert(inter_sel->type == ADIOS_SELECTION_BOUNDINGBOX);
    inter_bb = &inter_sel->u.bb;

    // Compute the offset of the intersection bounding box within each of
    // the source and destination bounding boxes
    assert(dst_bb->ndim == src_bb->ndim);
    inter_off_relative_to_dst = malloc(ndim * sizeof(uint64_t));
    inter_off_relative_to_src = malloc(ndim * sizeof(uint64_t));
    vector_sub(ndim, inter_off_relative_to_dst, inter_bb->start, dst_bb->start);
    vector_sub(ndim, inter_off_relative_to_src, inter_bb->start, src_bb->start);

    // Perform a subvolume memcpy
    copy_subvolume_ragged_offset(
        dst, src, dst_bb->ndim, inter_bb->count,
        dst_bb->count, inter_off_relative_to_dst, dst_ragged_offset,
        src_bb->count, inter_off_relative_to_src, src_ragged_offset,
        datum_type, swap_endianness);

    // Compute the number of elements copied
    volume = compute_volume(ndim, inter_bb->count);

    // Cleanup
    free(inter_off_relative_to_dst);
    free(inter_off_relative_to_src);
    common_read_selection_delete(inter_sel);

    return volume;
}

// Whenever we are patching between a point selection and bounding box, we
// will always iterate over the point selection, check each points for containment
// within the bounding box, and compute byte offsets for the point within the point list
// and bounding box, regardless of whether the point selection is on the source or
// destination buffer. Therefore, we include a helper function with a boolean flag
// to switch the copy direction, and just branch for a few lines during the copy
// operation. This simplifies the code a lot, and reduces the LoC in this file (although
// this comment makes up for a good bit of that savings).
static uint64_t adios_patch_data_bb_pts_helper(void *dst, uint64_t dst_ragged_offset, void *src, uint64_t src_ragged_offset,
                                               const ADIOS_SELECTION_POINTS_STRUCT *pts,
                                               const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb,
                                               _Bool isDestPoints, enum ADIOS_DATATYPES datum_type,
                                               enum ADIOS_FLAG swap_endianness) {
    const int ndim = pts->ndim;
    uint64_t i;
    int j;
    uint64_t pts_copied = 0;
    uint64_t byte_offset_in_bb_buffer, byte_offset_in_pt_buffer;
    const uint64_t *cur_pt;
    uint64_t *bb_byte_strides = malloc(sizeof(uint64_t) * ndim);
    uint64_t *pt_relative_to_bb = malloc(sizeof(uint64_t) * ndim);

    // Compute the strides into the source bounding box array
    int typelen = adios_get_type_size(datum_type, NULL);
    uint64_t bb_volume = typelen;
    for (j = ndim - 1; j >= 0; j--) {
        bb_byte_strides[j] = bb_volume;
        bb_volume *= bb->count[j];
    }

    // Check that the selection dimensions are compatible
    assert(pts->ndim == bb->ndim);

    // Check each point; if it's in the bounding box, perform a copy
    for (i = 0; i < pts->npoints; i++) {
        cur_pt = &pts->points[i * ndim];

        for (j = 0; j < ndim; j++) {
            // If the point's coordinate in some dimension is outside the bounding box
            if (cur_pt[j] < bb->start[j] ||
                cur_pt[j] >= bb->start[j] + bb->count[j]) {
                break;
            }
        }

        // If the point is within the bounding box
        if (j == ndim) {
            vector_sub(ndim, pt_relative_to_bb, cur_pt, bb->start);

            byte_offset_in_bb_buffer = 0;
            for (j = 0; j < ndim; j++)
                byte_offset_in_bb_buffer += pt_relative_to_bb[j] * bb_byte_strides[j];

            byte_offset_in_pt_buffer = i * typelen;

            if (isDestPoints) {
                assert(byte_offset_in_pt_buffer >= dst_ragged_offset);
                assert(byte_offset_in_bb_buffer >= src_ragged_offset);
                memcpy((char*)dst + byte_offset_in_pt_buffer - dst_ragged_offset, (char*)src + byte_offset_in_bb_buffer - src_ragged_offset, typelen);
            } else {
                assert(byte_offset_in_bb_buffer >= dst_ragged_offset);
                assert(byte_offset_in_pt_buffer >= src_ragged_offset);
                memcpy((char*)dst + byte_offset_in_bb_buffer - dst_ragged_offset, (char*)src + byte_offset_in_pt_buffer - src_ragged_offset, typelen);
            }
            pts_copied++;
        }
    }

    free(bb_byte_strides);
    free(pt_relative_to_bb);

    return pts_copied;
}

static uint64_t adios_patch_data_pts_to_bb(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb,
                                           void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION_POINTS_STRUCT *src_pts,
                                           enum ADIOS_DATATYPES datum_type,
                                           enum ADIOS_FLAG swap_endianness) {
    return adios_patch_data_bb_pts_helper(dst, dst_ragged_offset, src, src_ragged_offset, src_pts, dst_bb, false, datum_type, swap_endianness);
}

static uint64_t adios_patch_data_bb_to_pts(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_POINTS_STRUCT *dst_pts,
                                           void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb,
                                           enum ADIOS_DATATYPES datum_type,
                                           enum ADIOS_FLAG swap_endianness) {
    return adios_patch_data_bb_pts_helper(dst, dst_ragged_offset, src, src_ragged_offset, dst_pts, src_bb, true, datum_type, swap_endianness);
}

static uint64_t adios_patch_data_pts_to_pts(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_POINTS_STRUCT *dst_pts,
                                            void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION_POINTS_STRUCT *src_pts,
                                            enum ADIOS_DATATYPES datum_type,
                                            enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("points","points");
    return 0;
}

static uint64_t adios_patch_data_wb_to_wb(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb,
                                          void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb,
                                          enum ADIOS_DATATYPES datum_type,
                                          enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("writeblock","writeblock");
    return 0;
}

static uint64_t adios_patch_data_auto_to_auto(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_AUTO_STRUCT *dst_auto,
                                              void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION_AUTO_STRUCT *src_auto,
                                              enum ADIOS_DATATYPES datum_type,
                                              enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("auto","auto");
    return 0;
}


// One-to-any patch functions

static uint64_t adios_patch_data_to_bb(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb,
                                       void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION *src_sel,
                                       enum ADIOS_DATATYPES datum_type,
                                       enum ADIOS_FLAG swap_endianness) {
    switch (src_sel->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb = &src_sel->u.bb;
        return adios_patch_data_bb_to_bb(dst, dst_ragged_offset, dst_bb, src, src_ragged_offset, src_bb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *src_pts = &src_sel->u.points;
        return adios_patch_data_pts_to_bb(dst, dst_ragged_offset, dst_bb, src, src_ragged_offset, src_pts, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    case ADIOS_SELECTION_AUTO:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Incompatible selection types %d, %d were used while patching decoded "
                                                                      "transformed data into the user buffer (this is an error in the current "
                                                                      "transform plugin)", src_sel->type, ADIOS_SELECTION_BOUNDINGBOX);
        return 0;
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", src_sel->type);
        return 0;
    }

}

static uint64_t adios_patch_data_to_pts(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_POINTS_STRUCT *dst_pts,
                                        void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION *src_sel,
                                        enum ADIOS_DATATYPES datum_type,
                                        enum ADIOS_FLAG swap_endianness) {
    switch (src_sel->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb = &src_sel->u.bb;
        return adios_patch_data_bb_to_pts(dst, dst_ragged_offset, dst_pts, src, src_ragged_offset, src_bb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *src_pts = &src_sel->u.points;
        return adios_patch_data_pts_to_pts(dst, dst_ragged_offset, dst_pts, src, src_ragged_offset, src_pts, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    case ADIOS_SELECTION_AUTO:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Incompatible selection types %d, %d were used while patching decoded "
                                                                      "transformed data into the user buffer (this is an error in the current "
                                                                      "transform plugin)", src_sel->type, ADIOS_SELECTION_POINTS);
        return 0;
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", src_sel->type);
        return 0;
    }

}

static uint64_t adios_patch_data_to_wb(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb,
                                       void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION *src_sel,
                                       enum ADIOS_DATATYPES datum_type,
                                       enum ADIOS_FLAG swap_endianness) {
    switch (src_sel->type) {
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb = &src_sel->u.block;
        return adios_patch_data_wb_to_wb(dst, dst_ragged_offset, dst_wb, src, src_ragged_offset, src_wb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_BOUNDINGBOX:
    case ADIOS_SELECTION_POINTS:
    case ADIOS_SELECTION_AUTO:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Incompatible selection types %d, %d were used while patching decoded "
                                                                      "transformed data into the user buffer (this is an error in the current "
                                                                      "transform plugin)", src_sel->type, ADIOS_SELECTION_WRITEBLOCK);
        return 0;
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", src_sel->type);
        return 0;
    }
}

static uint64_t adios_patch_data_to_auto(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION_AUTO_STRUCT *dst_auto,
                                         void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION *src_sel,
                                         enum ADIOS_DATATYPES datum_type,
                                         enum ADIOS_FLAG swap_endianness) {
    switch (src_sel->type) {
    case ADIOS_SELECTION_AUTO:
    {
        const ADIOS_SELECTION_AUTO_STRUCT *src_auto = &src_sel->u.autosel;
        return adios_patch_data_auto_to_auto(dst, dst_ragged_offset, dst_auto, src, src_ragged_offset, src_auto, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_BOUNDINGBOX:
    case ADIOS_SELECTION_POINTS:
    case ADIOS_SELECTION_WRITEBLOCK:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Incompatible selection types %d, %d were used while patching decoded "
                                                                      "transformed data into the user buffer (this is an error in the current "
                                                                      "transform plugin)", src_sel->type, ADIOS_SELECTION_AUTO);
        return 0;
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", src_sel->type);
        return 0;
    }
}

//
// Any-on-any patch function
//

uint64_t adios_patch_data(void *dst, uint64_t dst_ragged_offset, const ADIOS_SELECTION *dst_sel,
                          void *src, uint64_t src_ragged_offset, const ADIOS_SELECTION *src_sel,
                          enum ADIOS_DATATYPES datum_type,
                          enum ADIOS_FLAG swap_endianness) {
    switch (dst_sel->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb = &dst_sel->u.bb;
        return adios_patch_data_to_bb(dst, dst_ragged_offset, dst_bb, src, src_ragged_offset, src_sel, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *dst_pts = &dst_sel->u.points;
        return adios_patch_data_to_pts(dst, dst_ragged_offset, dst_pts, src, src_ragged_offset, src_sel, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb = &dst_sel->u.block;
        return adios_patch_data_to_wb(dst, dst_ragged_offset, dst_wb, src, src_ragged_offset, src_sel, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_AUTO:
    {
        const ADIOS_SELECTION_AUTO_STRUCT *dst_auto = &dst_sel->u.autosel;
        return adios_patch_data_to_auto(dst, dst_ragged_offset, dst_auto, src, src_ragged_offset, src_sel, datum_type, swap_endianness);
    }
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", dst_sel->type);
        return 0;
    }
}
