/*
 * adios_subvolume.h
 *
 *  Created on: Jul 25, 2012
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_SUBVOLUME_H_
#define ADIOS_SUBVOLUME_H_

#include <stdint.h>
#include "public/adios_types.h"
#include "public/adios_selection.h"

typedef struct {
    int ndim;
    uint64_t *subv_dims;
    uint64_t *dst_dims, *dst_subv_offsets;
    uint64_t *src_dims, *src_subv_offsets;
} adios_subvolume_copy_spec;

/*
 * Copies a multi-dimensional subvolume from one buffer to another.
 *
 * 'src' and 'dst' are assumed to be buffers laid out in C ordering (the
 * first dimension is the slowest-changing). 'dst' is assumed to have
 * sufficient size and dimensions to accommodate the incoming subvolume.
 *
 * @param dst the destination buffer
 * @param src the source buffer
 * @param ndim the number of dimensions in the source and destination space
 * @param subv_dims the dimensions of the subvolume to copy
 * @param dst_dims the dimensions of the entire destination buffer
 * @param dst_subv_offsets the offsets at which to write the subvolume
 * @param src_dims the dimensions of the entire source buffer
 * @param src_subv_offsets the offsets from which to read the subvolume
 * @param datum_type the ADIOS datatype of the elements in the source and
 *        destination buffers
 */
void copy_subvolume(void *dst, const void *src, int ndim, const uint64_t *subv_dims,
                    const uint64_t *dst_dims, const uint64_t *dst_subv_offsets,
                    const uint64_t *src_dims, const uint64_t *src_subv_offsets,
                    enum ADIOS_DATATYPES datum_type,
                    enum ADIOS_FLAG swap_endianness);

void copy_subvolume_with_spec(void *dst, const void *src,
                              const const adios_subvolume_copy_spec *copy_spec,
                              enum ADIOS_DATATYPES datum_type,
                              enum ADIOS_FLAG swap_endianness);

/*
 * Converts a variable bounding box and a selection box into a subvolume
 * copy spec. The provided subvolume copy spec must be uninitialized, and must
 * be properly destroyed after use (including freeing its buffers, as they are
 * malloc()ed by this function.
 *
 * @param copy_spec the subvolume copy spec to initialize (it must not have
 *        been previously initialized, or must have been subsequently
 *        destroyed)
 * @param sel the bounding box selection
 * @param ldims the dimensions of the variable local bounding box
 * @param offsets the global offsets of the variable's local bounding box
 * @return non-zero if the selection box intersects the variable's bounding
 *         box, and zero if the selection box does not inersect the variable's
 *         bounding box. If zero is returned, the subvolume copy spec's
 *         contents will be undefined.
 */
int adios_selection_to_copy_spec(adios_subvolume_copy_spec *copy_spec,
                                 const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *sel,
                                 const uint64_t *ldims, const uint64_t *offsets);

// Extracts a selection corresponding to the subvolume within the source buffer
ADIOS_SELECTION * adios_copy_spec_to_src_selection(adios_subvolume_copy_spec *copy_spec);

// Extracts a selection corresponding to the subvolume within the destination buffer
ADIOS_SELECTION * adios_copy_spec_to_dst_selection(adios_subvolume_copy_spec *copy_spec);

void adios_subvolume_copy_spec_init(adios_subvolume_copy_spec *copy_spec,
                                    int ndim, uint64_t *subv_dims,
                                    uint64_t *dst_dims, uint64_t *dst_subv_offsets,
                                    uint64_t *src_dims, uint64_t *src_subv_offsets);

void adios_subvolume_copy_spec_destroy(adios_subvolume_copy_spec *copy_spec, int free_buffers);

int is_subvolume_src_covering(adios_subvolume_copy_spec subv_spec);

#endif /* ADIOS_SUBVOLUME_H_ */
