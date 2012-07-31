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
#include "public/adios_selection.h"

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
 *
 * copy_data_space calls copy_data_space helper with the following paramter
 * translations:
 * 1) Element size (i.e. 8 bytes for double, etc.) is considered an additional,
 *    fastest-varying dimension
 * 2) All "contiguous" fastest-varying dimensions are collapsed into a single
 *    fastest-varying dimension
 */
void copy_subvolume(void *dst, const void *src, int ndim, const uint64_t *subv_dims,
                    const uint64_t *dst_dims, const uint64_t *dst_subv_offsets,
                    const uint64_t *src_dims, const uint64_t *src_subv_offsets,
                    enum ADIOS_DATATYPES datum_type, enum ADIOS_FLAG swap_endianness) {

    int i;
    int first_contig_dim = 0; // Index of dimension that starts a contiguous block
    uint64_t src_strides[32];
    uint64_t dst_strides[32];
    int type_size = adios_get_type_size(datum_type, NULL);

    // Find the last dimension for which the subvolume, source and destination
    // spaces do not exactly match (i.e. non-zero offset or unequal lengths).
    for (i = 0; i < ndim; i++) {
        // If the subvolume, source and destination do not exactly match along
        // this dimension, it is marked as incomplete
        if (src_subv_offsets[i] != 0 ||
            dst_subv_offsets[i] != 0 ||
            subv_dims[i] != src_dims[i] ||
            subv_dims[i] != dst_dims[i]) {

            first_contig_dim = i;
        }
    }

    // Calculate the volume (number of bytes) of the region subtended by the
    // contiguous dimensions.
    uint64_t contig_dims_volume = 1;
    for (i = first_contig_dim; i < ndim; i++) {
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

    // Save the old value of the first contiguous dimension, then replace it
    // with a "collapsed" dimension consolidating all contiguous dimensions
    // into one
    // We "cheat" a bit by removing the const modifier. This is OK because we
    // carefully put back the original value before returning. This is the only
    // modification we do in this function besides filling 'dst'.
    uint64_t first_contig_dim_value_old = subv_dims[first_contig_dim];
    ((uint64_t*)subv_dims)[first_contig_dim] = contig_dims_volume;

    //printf(">>> copy_subvolume is using %d contiguous dimensions...\n", ndim - first_contig_dim);

    // Finally, delegate to the recursive worker function
    copy_subvolume_helper((char*)dst + dst_offset,			/* Offset dst buffer to the first element */
                          (char*)src + src_offset,			/* Offset src buffer to the first element */
                          first_contig_dim + 1,				/* Reduce dimensions to the non-contiguous dimensions, plus 1 for the collapsed contiguous dimensions */
                          subv_dims,						/* Subvolume dimensions (modified to collapse all contiguous dimensions) */
                          dst_strides,						/* dst buffer dimension strides */
                          src_strides,						/* src buffer dimension strides */
                          datum_type,						/* The datatype of the buffer elements */
                          swap_endianness == adios_flag_yes	/* Whether to swap endianness */
                          );

    // Restore the old first contiguous dimension value
    ((uint64_t*)subv_dims)[first_contig_dim] = first_contig_dim_value_old;
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

int is_subvolume_src_covered(const adios_subvolume_copy_spec *subv_spec) {
    int dim;
    for (dim = 0; dim < subv_spec->ndim; dim++) {
        if (subv_spec->src_subv_offsets[dim] != 0 ||
            subv_spec->src_dims[dim] != subv_spec->subv_dims[dim])
            return 0;
    }
    return 1;
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


int adios_selection_to_copy_spec(adios_subvolume_copy_spec *copy_spec,
                                 const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *sel,
                                 const uint64_t *ldims, const uint64_t *offsets) {
    int dim;
    const int ndim = sel->ndim;
    const int dimsize = sizeof(uint64_t) * ndim;

    copy_spec->ndim = ndim;
    copy_spec->subv_dims = malloc(dimsize);
    copy_spec->dst_dims = malloc(dimsize);
    copy_spec->dst_subv_offsets = malloc(dimsize);
    copy_spec->src_dims = malloc(dimsize);
    copy_spec->src_subv_offsets = malloc(dimsize);

    //printf(">>> Intersecting local bb dim1 start:count %llu:%llu with sel bb start:count %llu:%llu\n", offsets[0], ldims[0], sel->start[0], sel->count[0]);

    // Find the intersection between the selection box and the variable bounds
    for (dim = 0; dim < ndim; dim++) {
        const uint64_t sel_dim_goffset = sel->start[dim];
        const uint64_t sel_dim_len = sel->count[dim];
        const uint64_t varbox_dim_goffset = offsets[dim];
        const uint64_t varbox_dim_len = ldims[dim];
        uint64_t intersection_goffset, intersection_len;

        int intersects = intersect_segments(
                sel_dim_goffset, sel_dim_len,
                varbox_dim_goffset, varbox_dim_len,
                &intersection_goffset, &intersection_len);

        if (!intersects)
            return 0;

        // Initialize parameters for this dimension based on the intersection
        copy_spec->subv_dims[dim] = intersection_len;
        copy_spec->src_dims[dim] = varbox_dim_len;
        copy_spec->src_subv_offsets[dim] = intersection_goffset - varbox_dim_goffset;
        copy_spec->dst_dims[dim] = sel_dim_len;
        copy_spec->dst_subv_offsets[dim] = intersection_goffset - sel_dim_goffset;
    }

    return 1;
}

static void * bufdup(const void *buf, uint64_t elem_size, uint64_t count) {
    const uint64_t len = elem_size * count;
    void *newbuf = malloc(len);
    memcpy(newbuf, buf, len);
    return newbuf;
}

// Extracts a selection corresponding to the subvolume within the source buffer
ADIOS_SELECTION * adios_copy_spec_to_src_selection(adios_subvolume_copy_spec *copy_spec) {
    return common_read_selection_boundingbox(copy_spec->ndim,
                                             bufdup(copy_spec->src_subv_offsets, sizeof(uint64_t), copy_spec->ndim),
                                             bufdup(copy_spec->subv_dims, sizeof(uint64_t), copy_spec->ndim));
}

// Extracts a selection corresponding to the subvolume within the destination buffer
ADIOS_SELECTION * adios_copy_spec_to_dst_selection(adios_subvolume_copy_spec *copy_spec) {
    return common_read_selection_boundingbox(copy_spec->ndim,
                                             bufdup(copy_spec->dst_subv_offsets, sizeof(uint64_t), copy_spec->ndim),
                                             bufdup(copy_spec->subv_dims, sizeof(uint64_t), copy_spec->ndim));
}


void adios_subvolume_copy_spec_init(adios_subvolume_copy_spec *copy_spec,
                                    int ndim, uint64_t *subv_dims,
                                    uint64_t *dst_dims, uint64_t *dst_subv_offsets,
                                    uint64_t *src_dims, uint64_t *src_subv_offsets) {
    copy_spec->ndim = ndim;
    copy_spec->subv_dims = subv_dims;
    copy_spec->dst_dims = dst_dims;
    copy_spec->dst_subv_offsets = dst_subv_offsets;
    copy_spec->src_dims = src_dims;
    copy_spec->src_subv_offsets = src_subv_offsets;
}

#define MY_FREE(x) if (x) free(x); x = 0;
void adios_subvolume_copy_spec_destroy(adios_subvolume_copy_spec *copy_spec, int free_buffers) {
    if (free_buffers) {
        MY_FREE(copy_spec->subv_dims);
        MY_FREE(copy_spec->dst_dims);
        MY_FREE(copy_spec->dst_subv_offsets);
        MY_FREE(copy_spec->src_dims);
        MY_FREE(copy_spec->src_subv_offsets);
    }
    memset(copy_spec, 0, sizeof(adios_subvolume_copy_spec));
}
#undef MY_FREE
