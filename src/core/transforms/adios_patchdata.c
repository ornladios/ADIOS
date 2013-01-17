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
    const ADIOS_SELECTION *inter_sel = adios_selection_intersect_bb_bb(dst_bb, src_bb);
    const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *inter_bb;
    uint64_t *inter_off_relative_to_dst;
    uint64_t *inter_off_relative_to_src;
    uint64_t volume;

    // If there is no intersection, stop now, nothing to do
    if (!inter_sel)
        return 0;

    assert(inter_sel->type == ADIOS_SELECTION_BOUNDINGBOX);
    inter_bb = &inter_sel->u.bb;

    assert(dst_bb->ndim == src_bb->ndim);
    inter_off_relative_to_dst = malloc(ndim * sizeof(uint64_t));
    inter_off_relative_to_src = malloc(ndim * sizeof(uint64_t));

    vec_sub(ndim, inter_off_relative_to_dst, inter_bb->start, dst_bb->start);
    vec_sub(ndim, inter_off_relative_to_src, inter_bb->start, src_bb->start);

    copy_subvolume(dst, src, dst_bb->ndim, inter_bb->count,
                   dst_bb->count, inter_off_relative_to_dst,
                   src_bb->count, inter_off_relative_to_src,
                   datum_type, swap_endianness);

    volume = compute_volume(ndim, inter_bb->count);

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
    PATCH_UNIMPL("points","bounding box");
    return 0;
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
