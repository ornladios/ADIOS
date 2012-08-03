/*
 * adios_subvolume.c
 *
 *  Created on: Jul 25, 2012
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "util.h"
#include "adios_subvolume.h"

/*
 * Copies a given subvolume from 'src' to 'dst'. It recursively copies a
 * hyperplane of progressively lower dimensions, until it reaches the lowest
 * dimension, in which case it simply performs a memcpy.
 *
 * At any given call level to this function, 'src' and 'dst' point to the first
 * element of a hyperplane of the subvolume with dimension 'ndim'.
 *
 * @param dst the destination buffer
 * @param src the source buffer
 * @param ndim the number of dimensions
 * @param next_subv_dim the list of dimensions for the subvolume, starting with
 *        the slowest-changing
 * @param next_dst_stride the list of destination buffer strides for each
 *        dimension
 * @param next_src_stride the list of source buffer strides for each dimension
 * @param type the datatype of the elements in dst/src (NOTE: this is ONLY for
 *        the purpose of endianness swapping; strides/dimensions are in bytes)
 * @param swap_endianness if true, swap the endianness of each element
 */
static void copy_subvolume_helper(char *dst, const char *src,
                                  int ndim, const uint64_t *next_subv_dim,
                                  const uint64_t *next_dst_stride, const uint64_t *next_src_stride,
                                  enum ADIOS_DATATYPES buftype, int swap_endianness) {
    if (ndim == 1) {
        memcpy(dst, src, *next_subv_dim);
        if (swap_endianness) {
            change_endianness(dst, *next_subv_dim, buftype);
        }
    } else {
        int i;
        for (i = 0; i < *next_subv_dim; i++) {
            copy_subvolume_helper(dst, src, ndim - 1,
                                  next_subv_dim + 1, next_dst_stride + 1, next_src_stride + 1,
                                  buftype, swap_endianness);

            src += *next_src_stride;
            dst += *next_dst_stride;
        }
    }
}

/*
 * copy_subvolume delegates to the copy_subvolume_helper function, with the
 * following parameter translations:
 * 1) Element size (i.e. 8 bytes for double, etc.) is considered an additional,
 *    fastest-varying dimension
 * 2) All "contiguous" fastest-varying dimensions are collapsed into a single
 *    fastest-varying dimension
 */
void copy_subvolume(void *dst, const void *src, int ndim, const uint64_t *subv_dims,
                    const uint64_t *dst_dims, const uint64_t *dst_subv_offsets,
                    const uint64_t *src_dims, const uint64_t *src_subv_offsets,
                    enum ADIOS_DATATYPES datum_type, enum ADIOS_FLAG swap_endianness) {
    copy_subvolume_ragged(dst, src, ndim, subv_dims,
                          dst_dims, dst_subv_offsets, NULL,
                          src_dims, src_subv_offsets, NULL,
                          datum_type, swap_endianness);
}

/*
 *
 * {src,dst}_ragged_offsets are offset vectors that represent the position of
 * the first element in ragged arrays {src,dst}. In order words, the argument
 * {src,dst} logically points into the complete array to the element at
 * {src,dst}_ragged_offsets, and the data continues from there in the normal
 * order.
 *
 * {src,dst}_ragged_offsets may be NULL, in which case it is considered to be a
 * zero vector.
 */
void copy_subvolume_ragged(void *dst, const void *src, int ndim, const uint64_t *subv_dims,
                           const uint64_t *dst_dims, const uint64_t *dst_subv_offsets,
                           const uint64_t *dst_ragged_offsets,
                           const uint64_t *src_dims, const uint64_t *src_subv_offsets,
                           const uint64_t *src_ragged_offsets,
                           enum ADIOS_DATATYPES datum_type, enum ADIOS_FLAG swap_endianness) {

    const uint64_t src_ragged_offset =
            src_ragged_offsets ?
            compute_ragged_array_offset(ndim, src_ragged_offsets,
                                        src_dims, datum_type) : 0;
    const uint64_t dst_ragged_offset =
            dst_ragged_offsets ?
            compute_ragged_array_offset(ndim, dst_ragged_offsets,
                                        dst_dims, datum_type) : 0;

    copy_subvolume_ragged_offset(dst, src, ndim, subv_dims,
                                 dst_dims, dst_subv_offsets,
                                 dst_ragged_offsets,
                                 src_dims, src_subv_offsets,
                                 src_ragged_offsets,
                                 datum_type, swap_endianness);
}

/*
 * {src,dst}_ragged_offset are the byte offsets of the start of the
 * corresponding arrays relative to the logical start of the complete arrays.
 */
void copy_subvolume_ragged_offset(void *dst, const void *src, int ndim, const uint64_t *subv_dims,
                                  const uint64_t *dst_dims, const uint64_t *dst_subv_offsets,
                                  uint64_t dst_ragged_offset,
                                  const uint64_t *src_dims, const uint64_t *src_subv_offsets,
                                  uint64_t src_ragged_offset,
                                  enum ADIOS_DATATYPES datum_type, enum ADIOS_FLAG swap_endianness) {

    int i;
    int last_noncovering_dim = 0; // Index of dimension that starts a contiguous block
    uint64_t src_strides[32];
    uint64_t dst_strides[32];
    const int type_size = adios_get_type_size(datum_type, NULL);

    // Find the last dimension for which the subvolume, source and destination
    // spaces do not exactly match (i.e. non-zero offset or unequal lengths).
    for (i = 0; i < ndim; i++) {
        // If the subvolume, source and destination do not exactly match along
        // this dimension, it is marked as incomplete
        if (src_subv_offsets[i] != 0 ||
            dst_subv_offsets[i] != 0 ||
            subv_dims[i] != src_dims[i] ||
            subv_dims[i] != dst_dims[i]) {

            last_noncovering_dim = i;
        }
    }

    // Calculate the volume (number of bytes) of the region subtended by the
    // contiguous dimensions (which start with the last non-covering dimension)
    uint64_t contig_dims_volume = 1;
    for (i = last_noncovering_dim; i < ndim; i++) {
        contig_dims_volume *= subv_dims[i];
    }
    // Add the element size as a new "dimension", to convert from element-space
    // to byte-space
    contig_dims_volume *= type_size; // Assumes non-string type

    // Compute strides for the dimensions
    uint64_t src_volume = type_size;
    uint64_t dst_volume = type_size;
    for (i = ndim - 1; i >= 0; i--) {
        src_strides[i] = src_volume;
        dst_strides[i] = dst_volume;

        src_volume *= src_dims[i];
        dst_volume *= dst_dims[i];
    }

    // Compute the starting offsets for src and dst. In theory we could skip
    // the contiguous dimensions other than the first one, but just to be
    // safe/simple we calculate them all.
    uint64_t src_offset = 0, dst_offset = 0;
    for (i = 0; i < ndim; i++) {
        src_offset += src_subv_offsets[i] * src_strides[i];
        dst_offset += dst_subv_offsets[i] * dst_strides[i];
    }
    // Incorporate ragged offsets
    src_offset -= src_ragged_offset;
    dst_offset -= dst_ragged_offset;

    // Save the old value of the first contiguous dimension, then replace it
    // with a "collapsed" dimension consolidating all contiguous dimensions
    // into one
    // We "cheat" a bit by removing the const modifier. This is OK because we
    // carefully put back the original value before returning. This is the only
    // argument modification we do in this function besides filling 'dst'.
    uint64_t first_contig_dim_value_old = subv_dims[last_noncovering_dim];
    ((uint64_t*)subv_dims)[last_noncovering_dim] = contig_dims_volume;

    //printf(">>> copy_subvolume is using %d contiguous dimensions...\n", ndim - first_contig_dim);

    // Finally, delegate to the recursive worker function
    copy_subvolume_helper((char*)dst + dst_offset,			/* Offset dst buffer to the first element */
                          (char*)src + src_offset,			/* Offset src buffer to the first element */
                          last_noncovering_dim + 1,			/* Reduce dimensions to the non-contiguous dimensions, plus 1 for the collapsed contiguous dimensions */
                          subv_dims,						/* Subvolume dimensions (modified to collapse all contiguous dimensions) */
                          dst_strides,						/* dst buffer dimension strides */
                          src_strides,						/* src buffer dimension strides */
                          datum_type,						/* The datatype of the buffer elements */
                          swap_endianness == adios_flag_yes	/* Whether to swap endianness */
                          );

    // Restore the old first contiguous dimension value
    ((uint64_t*)subv_dims)[last_noncovering_dim] = first_contig_dim_value_old;
}

void copy_subvolume_with_spec(void *dst, const void *src,
                              const adios_subvolume_copy_spec *copy_spec,
                              enum ADIOS_DATATYPES datum_type,
                              enum ADIOS_FLAG swap_endianness) {
    copy_subvolume(dst, src, copy_spec->ndim, copy_spec->subv_dims,
                   copy_spec->dst_dims, copy_spec->dst_subv_offsets,
                   copy_spec->src_dims, copy_spec->src_subv_offsets,
                   datum_type, swap_endianness);
}

static int intersect_segments(uint64_t start1, uint64_t len1, uint64_t start2, uint64_t len2,
                              uint64_t *inter_start, uint64_t *inter_len) {
    int end1, end2;
    int inter_end;

    // Swap the segments if segment 1 starts after segment 2, to simplify calculations
    if (start1 > start2) {
        int tmps = start1;
        int tmpl = len1;
        start1 = start2;
        len1 = len2;
        start2 = tmps;
        len2 = tmpl;
    }

    end1 = start1 + len1;
    end2 = start2 + len2;

    // If segment 1 ends before segment 2 starts, no intersection
    if (end1 <= start2)
        return 0;

    *inter_start = start2;						// Intersection starts at the beginning of the later segment
    inter_end = end1 <= end2 ? end1 : end2;	// Intersection ends at the earlier of the two segment ends
    *inter_len = inter_end - *inter_start;

    return 1;
}

int intersect_volumes(int ndim,
                      const uint64_t *offset1, const uint64_t *dims1,
                      const uint64_t *offset2, const uint64_t *dims2,
                      uint64_t *inter_offset,
                      uint64_t *inter_offset_rel1, uint64_t *inter_offset_rel2,
                      uint64_t *inter_dims) {
    // For every dimension, find the intersection in that dimension
    // If ever the volumes are disjoint in some dimension, stop immediately
    // and report no intersection
    int dim;
    uint64_t tmp_inter_offset;
    for (dim = 0; dim < ndim; dim++) {
        const int intersects = intersect_segments(*offset1, *dims1,
                                                  *offset2, *dims2,
                                                  &tmp_inter_offset, inter_dims);
        if (!intersects)
            return 0;

        // Calculate/store offsets as the user requests, and then advance to
        // the next dimension on these as well
        if (inter_offset)
            *inter_offset++ = tmp_inter_offset;
        if (inter_offset_rel1)
            *inter_offset_rel1++ = tmp_inter_offset - *offset1;
        if (inter_offset_rel2)
            *inter_offset_rel2++ = tmp_inter_offset - *offset2;

        // Advance the other arrays to the next dimension, as well
        // NOTE: this must be done after offset calculations, as offset[12] are
        // accessed
        offset1++; dims1++;
        offset2++; dims2++;
        inter_dims++;
    }

    // If we have a non-null intersection in every dimension, the entire
    // volumes intersect, so return true to indicate this
    return 1;
}

uint64_t compute_ragged_array_offset(int ndim, const uint64_t *start_offset, const uint64_t *overall_dims, enum ADIOS_DATATYPES elem_type) {
    int dim;
    uint64_t ragged_off = 0;
    uint64_t volume_so_far = adios_get_type_size(elem_type, NULL);
    for (dim = ndim - 1; dim >= 0; dim--) {
        ragged_off += start_offset[dim] * volume_so_far;
        volume_so_far *= overall_dims[dim];
    }
    return ragged_off;
}
