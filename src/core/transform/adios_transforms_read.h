/*
 * Contains read-specific code for handling variable transforms in ADIOS
 *
 *  Created on: Jun 27, 2012
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_TRANSFORMS_READ_H_
#define ADIOS_TRANSFORMS_READ_H_

#include "public/adios_error.h"
#include "public/adios_types.h"
#include "adios_transforms_common.h"
#include "adios_subvolume.h"

/*
 * A read delegate function that can be called by
 * adios_transform_retrieve_transformed_data_subvolume to read a specified
 * chunk of a variable's data. The data is contained in the returned buffer.
 * Do not free() this buffer; it is managed by the read method.
 *
 * slice_offset is relative to the variable's payload.
 */
typedef void * (*adios_transform_var_read_delegate)(struct adios_index_var_struct_v1 *v,
                                                    int time_index,
                                                    uint64_t slice_offset,
                                                    uint64_t slice_size,
                                                    void *state);

/*
 * Reads and detransforms the given subvolume of the given variable.
 * @param var the variable to read from
 * @param out a data buffer into which the detransformed data should be stored
 * @param out_size an output variable for the the number of bytes written to "out"
 * @param read_state opaque state data to be passed to the read delegate
 * @param read_delegate a read delegate function that can be called to read chunks
 *        of data
 * @param sel the subvolume selected to be retrieved
 */
enum ADIOS_ERRCODES adios_transform_retrieve_transformed_data_subvolume(
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *copy_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness);

#endif /* ADIOS_TRANSFORMS_WRITE_H_ */
