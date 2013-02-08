/*
 * adios_patchdata.c
 *
 *  Created on: Jan 15, 2013
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
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

static inline void vec_sub(int ndim, uint64_t *dst, const uint64_t *a, const uint64_t *b) {
    while (ndim--)
        *dst++ = *a++ - *b++;
}

// One-to-one patch functions
static uint64_t adios_patch_data_bb_to_bb(void *dst, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb,
                                          void *src, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb,
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

    // Ensure the intersection is actually a bounding box
    // (this is guaranteed by the selection intersection code; this is just to check for bugs)
    assert(inter_sel->type == ADIOS_SELECTION_BOUNDINGBOX);
    inter_bb = &inter_sel->u.bb;

    // Compute the offset of the intersection bounding box within each of
    // the source and destination bounding boxes
    assert(dst_bb->ndim == src_bb->ndim);
    inter_off_relative_to_dst = malloc(ndim * sizeof(uint64_t));
    inter_off_relative_to_src = malloc(ndim * sizeof(uint64_t));
    vec_sub(ndim, inter_off_relative_to_dst, inter_bb->start, dst_bb->start);
    vec_sub(ndim, inter_off_relative_to_src, inter_bb->start, src_bb->start);

    // Perform a subvolume memcpy
    copy_subvolume(dst, src, dst_bb->ndim, inter_bb->count,
                   dst_bb->count, inter_off_relative_to_dst,
                   src_bb->count, inter_off_relative_to_src,
                   datum_type, swap_endianness);

    // Compute the number of elements copied
    volume = compute_volume(ndim, inter_bb->count);

    // Cleanup
    free(inter_off_relative_to_dst);
    free(inter_off_relative_to_src);
    common_read_selection_delete(inter_sel);

    return volume;
}

static uint64_t adios_patch_data_pts_to_bb(void *dst, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb,
                                           void *src, const ADIOS_SELECTION_POINTS_STRUCT *src_pts,
                                           enum ADIOS_DATATYPES datum_type,
                                           enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("bounding box","points");
    return 0;
}

static uint64_t adios_patch_data_wb_to_bb(void *dst, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb,
                                          void *src, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb,
                                          enum ADIOS_DATATYPES datum_type,
                                          enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("bounding box","writeblock");
    return 0;
}

static uint64_t adios_patch_data_auto_to_bb(void *dst, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb,
                                            void *src, const ADIOS_SELECTION_AUTO_STRUCT *src_auto,
                                            enum ADIOS_DATATYPES datum_type,
                                            enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("bounding box","auto");
    return 0;
}

static uint64_t adios_patch_data_bb_to_pts(void *dst, const ADIOS_SELECTION_POINTS_STRUCT *dst_pts,
                                           void *src, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb,
                                           enum ADIOS_DATATYPES datum_type,
                                           enum ADIOS_FLAG swap_endianness) {
    const int ndim = dst_pts->ndim;
    uint64_t i;
    int j;
    uint64_t pts_copied = 0;
    uint64_t byte_offset_in_src;
    const uint64_t *cur_pt;
    uint64_t *src_strides = malloc(sizeof(uint64_t) * ndim);
    uint64_t *pt_relative_to_src = malloc(sizeof(uint64_t) * ndim);

    // Compute the strides into the source bounding box array
    int typelen = adios_get_type_size(datum_type, NULL);
    uint64_t src_volume = typelen;
    for (j = ndim - 1; j >= 0; j--) {
        src_strides[j] = src_volume;
        src_volume *= src_bb->count[j];
    }

    // Check that the selection dimensions are compatible
    assert(dst_pts->ndim == src_bb->ndim);

    // Check each point in the destination; if it's in the source bounding box, copy it over
    for (i = 0; i < dst_pts->npoints; i++) {
        cur_pt = &dst_pts->points[i * ndim];

        for (j = 0; j < ndim; j++) {
            // If the point's coordinate in some dimension is outside the bounding box
            if (cur_pt[j] < src_bb->start[j] ||
                cur_pt[j] >= src_bb->start[j] + src_bb->count[j]) {
                break;
            }
        }

        // If the point is within the bounding box
        if (j == ndim) {
            vec_sub(ndim, pt_relative_to_src, cur_pt, src_bb->start);

            byte_offset_in_src = 0;
            for (j = 0; j < ndim; j++)
                byte_offset_in_src += pt_relative_to_src[j] * src_strides[j];

            memcpy((char*)dst + i * typelen, (char*)src + byte_offset_in_src, typelen);
            pts_copied++;

            printf("Copied into point at index %llu!\n", i);
        }
    }

    free(src_strides);
    free(pt_relative_to_src);

    printf("Copied %llu points!\n", pts_copied);

    return pts_copied;
}

static uint64_t adios_patch_data_pts_to_pts(void *dst, const ADIOS_SELECTION_POINTS_STRUCT *dst_pts,
                                            void *src, const ADIOS_SELECTION_POINTS_STRUCT *src_pts,
                                            enum ADIOS_DATATYPES datum_type,
                                            enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("points","points");
    return 0;
}

static uint64_t adios_patch_data_wb_to_pts(void *dst, const ADIOS_SELECTION_POINTS_STRUCT *dst_pts,
                                           void *src, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb,
                                           enum ADIOS_DATATYPES datum_type,
                                           enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("points","writeblock");
    return 0;
}

static uint64_t adios_patch_data_auto_to_pts(void *dst, const ADIOS_SELECTION_POINTS_STRUCT *dst_pts,
                                             void *src, const ADIOS_SELECTION_AUTO_STRUCT *src_auto,
                                             enum ADIOS_DATATYPES datum_type,
                                             enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("points","auto");
    return 0;
}

static uint64_t adios_patch_data_bb_to_wb(void *dst, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb,
                                          void *src, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb,
                                          enum ADIOS_DATATYPES datum_type,
                                          enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("writeblock","bounding box");
    return 0;
}

static uint64_t adios_patch_data_pts_to_wb(void *dst, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb,
                                           void *src, const ADIOS_SELECTION_POINTS_STRUCT *src_pts,
                                           enum ADIOS_DATATYPES datum_type,
                                           enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("writeblock","points");
    return 0;
}

static uint64_t adios_patch_data_wb_to_wb(void *dst, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb,
                                          void *src, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb,
                                          enum ADIOS_DATATYPES datum_type,
                                          enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("writeblock","writeblock");
    return 0;
}

static uint64_t adios_patch_data_auto_to_wb(void *dst, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb,
                                            void *src, const ADIOS_SELECTION_AUTO_STRUCT *src_auto,
                                            enum ADIOS_DATATYPES datum_type,
                                            enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("writeblock","auto");
    return 0;
}

static uint64_t adios_patch_data_bb_to_auto(void *dst, const ADIOS_SELECTION_AUTO_STRUCT *dst_auto,
                                            void *src, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb,
                                            enum ADIOS_DATATYPES datum_type,
                                            enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("auto","bounding box");
    return 0;
}

static uint64_t adios_patch_data_pts_to_auto(void *dst, const ADIOS_SELECTION_AUTO_STRUCT *dst_auto,
                                             void *src, const ADIOS_SELECTION_POINTS_STRUCT *src_pts,
                                             enum ADIOS_DATATYPES datum_type,
                                             enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("auto","points");
    return 0;
}

static uint64_t adios_patch_data_wb_to_auto(void *dst, const ADIOS_SELECTION_AUTO_STRUCT *dst_auto,
                                            void *src, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb,
                                            enum ADIOS_DATATYPES datum_type,
                                            enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("auto","writeblock");
    return 0;
}

static uint64_t adios_patch_data_auto_to_auto(void *dst, const ADIOS_SELECTION_AUTO_STRUCT *dst_auto,
                                              void *src, const ADIOS_SELECTION_AUTO_STRUCT *src_auto,
                                              enum ADIOS_DATATYPES datum_type,
                                              enum ADIOS_FLAG swap_endianness) {
    PATCH_UNIMPL("auto","auto");
    return 0;
}


// One-to-any patch functions

static uint64_t adios_patch_data_to_bb(void *dst, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb,
                                       void *src, const ADIOS_SELECTION *src_sel,
                                       enum ADIOS_DATATYPES datum_type,
                                       enum ADIOS_FLAG swap_endianness) {
    switch (src_sel->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb = &src_sel->u.bb;
        return adios_patch_data_bb_to_bb(dst, dst_bb, src, src_bb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *src_pts = &src_sel->u.points;
        return adios_patch_data_pts_to_bb(dst, dst_bb, src, src_pts, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb = &src_sel->u.block;
        return adios_patch_data_wb_to_bb(dst, dst_bb, src, src_wb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_AUTO:
    {
        const ADIOS_SELECTION_AUTO_STRUCT *src_auto = &src_sel->u.autosel;
        return adios_patch_data_auto_to_bb(dst, dst_bb, src, src_auto, datum_type, swap_endianness);
    }
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", src_sel->type);
        return 0;
    }

}

static uint64_t adios_patch_data_to_pts(void *dst, const ADIOS_SELECTION_POINTS_STRUCT *dst_pts,
                                        void *src, const ADIOS_SELECTION *src_sel,
                                        enum ADIOS_DATATYPES datum_type,
                                        enum ADIOS_FLAG swap_endianness) {
    switch (src_sel->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb = &src_sel->u.bb;
        return adios_patch_data_bb_to_pts(dst, dst_pts, src, src_bb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *src_pts = &src_sel->u.points;
        return adios_patch_data_pts_to_pts(dst, dst_pts, src, src_pts, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb = &src_sel->u.block;
        return adios_patch_data_wb_to_pts(dst, dst_pts, src, src_wb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_AUTO:
    {
        const ADIOS_SELECTION_AUTO_STRUCT *src_auto = &src_sel->u.autosel;
        return adios_patch_data_auto_to_pts(dst, dst_pts, src, src_auto, datum_type, swap_endianness);
    }
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", src_sel->type);
        return 0;
    }

}

static uint64_t adios_patch_data_to_wb(void *dst, const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb,
                                       void *src, const ADIOS_SELECTION *src_sel,
                                       enum ADIOS_DATATYPES datum_type,
                                       enum ADIOS_FLAG swap_endianness) {
    switch (src_sel->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb = &src_sel->u.bb;
        return adios_patch_data_bb_to_wb(dst, dst_wb, src, src_bb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *src_pts = &src_sel->u.points;
        return adios_patch_data_pts_to_wb(dst, dst_wb, src, src_pts, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb = &src_sel->u.block;
        return adios_patch_data_wb_to_wb(dst, dst_wb, src, src_wb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_AUTO:
    {
        const ADIOS_SELECTION_AUTO_STRUCT *src_auto = &src_sel->u.autosel;
        return adios_patch_data_auto_to_wb(dst, dst_wb, src, src_auto, datum_type, swap_endianness);
    }
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", src_sel->type);
        return 0;
    }
}

static uint64_t adios_patch_data_to_auto(void *dst, const ADIOS_SELECTION_AUTO_STRUCT *dst_auto,
                                         void *src, const ADIOS_SELECTION *src_sel,
                                         enum ADIOS_DATATYPES datum_type,
                                         enum ADIOS_FLAG swap_endianness) {
    switch (src_sel->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *src_bb = &src_sel->u.bb;
        return adios_patch_data_bb_to_auto(dst, dst_auto, src, src_bb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *src_pts = &src_sel->u.points;
        return adios_patch_data_pts_to_auto(dst, dst_auto, src, src_pts, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *src_wb = &src_sel->u.block;
        return adios_patch_data_wb_to_auto(dst, dst_auto, src, src_wb, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_AUTO:
    {
        const ADIOS_SELECTION_AUTO_STRUCT *src_auto = &src_sel->u.autosel;
        return adios_patch_data_auto_to_auto(dst, dst_auto, src, src_auto, datum_type, swap_endianness);
    }
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", src_sel->type);
        return 0;
    }
}

//
// Any-on-any patch function
//

uint64_t adios_patch_data(void *dst, const ADIOS_SELECTION *dst_sel,
                          void *src, const ADIOS_SELECTION *src_sel,
                          enum ADIOS_DATATYPES datum_type,
                          enum ADIOS_FLAG swap_endianness) {
    switch (dst_sel->type) {
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *dst_bb = &dst_sel->u.bb;
        return adios_patch_data_to_bb(dst, dst_bb, src, src_sel, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_POINTS:
    {
        const ADIOS_SELECTION_POINTS_STRUCT *dst_pts = &dst_sel->u.points;
        return adios_patch_data_to_pts(dst, dst_pts, src, src_sel, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        const ADIOS_SELECTION_WRITEBLOCK_STRUCT *dst_wb = &dst_sel->u.block;
        return adios_patch_data_to_wb(dst, dst_wb, src, src_sel, datum_type, swap_endianness);
    }
    case ADIOS_SELECTION_AUTO:
    {
        const ADIOS_SELECTION_AUTO_STRUCT *dst_auto = &dst_sel->u.autosel;
        return adios_patch_data_to_auto(dst, dst_auto, src, src_sel, datum_type, swap_endianness);
    }
    default:
        adios_error_at_line(err_invalid_argument, __FILE__, __LINE__, "Unknown selection type %d", dst_sel->type);
        return 0;
    }
}
