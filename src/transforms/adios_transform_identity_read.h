/*
 * adios_transform_identity_read.h
 *
 *  Created on: Apr 12, 2013
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_TRANSFORM_IDENTITY_READ_H_
#define ADIOS_TRANSFORM_IDENTITY_READ_H_

// This function is shared with APLOD
/*
 * Computes the linearized element start and end offsets of a selection (intersect_sel) within
 * a PG (bounded by pgbb). If all elements between these two offsets is read, and then
 * returned as a datablock with ragged offset equal to the start offset returned here,
 * all data within intersect_sel will be present and results patching will work.
 * @param intersect_sel the intersecting portion of some global read selection with a PG. It must
 *        be fully contained in the PG bounds, as specified by pgbb.
 * @param pgbb the bounds of the PG
 * @param start_off_ptr a pointer to where the start offset will be stored
 * @param end_off_ptr a pointer to where the end offset will be stored
 */
void compute_sieving_offsets_for_pg_selection(const ADIOS_SELECTION *intersect_sel,
                                              const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *pgbb,
                                              uint64_t *start_off_ptr, uint64_t *end_off_ptr);

#endif /* ADIOS_TRANSFORM_IDENTITY_READ_H_ */
