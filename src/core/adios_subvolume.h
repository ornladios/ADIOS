/*
 * adios_subvolume.h
 *
 * Utility functions for manipulating subvolumes of multi-dimensional arrays.
 *
 *  Created on: Jul 25, 2012
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_SUBVOLUME_H_
#define ADIOS_SUBVOLUME_H_

#include <stdint.h>
#include "public/adios_types.h"
#include "public/adios_selection.h"
#include "core/adios_copyspec.h"

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

/*
 * The same as copy_subvolume, with the addition of optional ragged src/dst
 * arrays. These arrays are ragged iff the pointer supplied does not point to
 * the logical (0,0,...,0) element of the corresponding array, but instead
 * points to some element (r1,r2,...,rn) with some ri != 0. In this case,
 * the corresponding {src,dst}_ragged_offsets designates the element pointed to
 * by the corresponding pointer.
 */
void copy_subvolume_ragged(void *dst, const void *src, int ndim, const uint64_t *subv_dims,
                           const uint64_t *dst_dims, const uint64_t *dst_subv_offsets,
                           const uint64_t *dst_ragged_offsets,
                           const uint64_t *src_dims, const uint64_t *src_subv_offsets,
                           const uint64_t *src_ragged_offsets,
                           enum ADIOS_DATATYPES datum_type, enum ADIOS_FLAG swap_endianness);

/*
 * Same as copy_subvolume_ragged, but takes a scalar byte offset for ragged
 * arrays instead of an array of element offsets.
 */
void copy_subvolume_ragged_offset(void *dst, const void *src, int ndim, const uint64_t *subv_dims,
                                  const uint64_t *dst_dims, const uint64_t *dst_subv_offsets,
                                  uint64_t dst_ragged_offset,
                                  const uint64_t *src_dims, const uint64_t *src_subv_offsets,
                                  uint64_t src_ragged_offset,
                                  enum ADIOS_DATATYPES datum_type, enum ADIOS_FLAG swap_endianness);

void copy_subvolume_with_spec(void *dst, const void *src,
                              const const adios_subvolume_copy_spec *copy_spec,
                              enum ADIOS_DATATYPES datum_type,
                              enum ADIOS_FLAG swap_endianness);

/*
 * Calculates the intersection, if any, between two volumes. For each volume,
 * dimensions and global offsets must be specified. If the volumes do
 * intersect, the size dimensions of the intersection are returned, as well as
 * the offset of the intersection in three forms: global, and relative to each
 * of the two volumes.
 *
 * All global offsets (offset1, offset2, inter_offset) are defined relative to
 * the same global coordinate space (typically a global array).
 *
 * All buffer arguments (everything but ndim) must be arrays of uint64_t of
 * length at least ndim (or NULL, in the case of the optional output arguments;
 * see below).
 *
 * If the volumes intersect, this function will return a non-zero value, and
 * the output arguments (inter_offset, inter_offset_rel1, inter_offset_rel2,
 * and inter_dims) will be populated with the global offset, offset relative
 * to volume 1's offset, offset relative to volume 2's offset, and dimensions
 * of the intersection volume, respectively. An appropriate buffer must be
 * supplied for inter_dims, but NULL may be supplied for any of the other
 * three output parameters if that information is not desired.
 *
 * If the volumes are disjoint, this function will return 0, and the content of
 * the output arguments is undefined.
 *
 * @param offset1 the global offset of volume 1
 * @param dims1 the dimensions of volume 1
 * @param offset2 the global offset of volume 2
 * @param dims2 the dimensions of volume 2
 * @param inter_offset a buffer to hold the offset of the intersection
 *        volume, or NULL if this information isn't required.
 * @param inter_offset_rel1 a buffer to hold the offset of the intersection
 *        volume relative to offset1, or NULL if this information isn't
 *        required.
 * @param inter_offset_rel2 a buffer to hold the offset of the intersection
 *        volume relative to offset2, or NULL if this information isn't
 *        required.
 * @param inter_dims a buffer to hold the dimensions of the intersection volume
 */
int intersect_volumes(int ndim,
                      const uint64_t *offset1, const uint64_t *dims1,
                      const uint64_t *offset2, const uint64_t *dims2,
                      uint64_t *inter_offset,
                      uint64_t *inter_offset_rel1, uint64_t *inter_offset_rel2,
                      uint64_t *inter_dims);

/*
 * Computes the byte offset of the beginning of a ragged multidimensional
 * volume array relative to the beginning of the corresponding complete volume
 * array.
 *
 * @param ndim number of dimensions of the volume
 * @param start_offset the offsets of the start of the ragged array
 * @param overall_dims the dimensions of the complete array
 * @param elem_type the datatype of each element in the array
 */
uint64_t compute_ragged_array_offset(int ndim, const uint64_t *start_offset, const uint64_t *overall_dims, enum ADIOS_DATATYPES elem_type);

#endif /* ADIOS_SUBVOLUME_H_ */
