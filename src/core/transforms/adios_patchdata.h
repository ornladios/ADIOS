/*
 * adios_patchdata.h
 *
 *  Created on: Jan 15, 2013
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_PATCHDATA_H_
#define ADIOS_PATCHDATA_H_

#include <public/adios_types.h>
#include <public/adios_selection.h>

/*
 * Copies data from a source buffer to a destination buffer, where the buffers each contain
 * data to fill given ADIOS_SELECTIONs. The data to be copied is the intersection of the source
 * and destination buffers, as defined by their selections.
 */
uint64_t adios_patch_data(void *dst, const ADIOS_SELECTION *dst_sel,
                          void *src, const ADIOS_SELECTION *src_sel,
                          enum ADIOS_DATATYPES datum_type,
                          enum ADIOS_FLAG swap_endianness);

#endif /* ADIOS_PATCHDATA_H_ */
