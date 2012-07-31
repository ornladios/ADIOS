
#ifdef COMPRESS

#include <stdint.h>
#include <assert.h>

#include "adios_logger.h"
#include "adios_transforms_hooks.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_write.h"
#include "public/adios_error.h"

#include "compress.h"

uint64_t adios_transform_compress_calc_vars_transformed_size(uint64_t orig_size, int num_vars) {
    return EXPAND_SIZE(orig_size);
}

// *apply function is in adios_transform_compress_write.c, since it requires
// access to functions from adios_internals.c, and therefore cannot link in
// read-only libraries.

typedef enum {
    DISJOINT,			// Bounding box and variable PG are completely disjoint
    PARTIALLY_COVERED,	// Bounding box only partially covers variable PG
    COVERED,			// Bounding box completely covers variable PG
    UNKNOWN				// Some unknown selection method
} bb_intersetion_type;

enum ADIOS_ERRCODES adios_transform_compress_retrieve_subvolume(
                    struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
                    const adios_subvolume_copy_spec *copy_spec,
                    void *read_state, adios_transform_var_read_delegate read_delegate,
                    enum ADIOS_FLAG swap_endianness) {

    assert(var->characteristics[time_index].transform.transform_type == adios_transform_compress);
    // struct adios_index_characteristic_dims_struct_v1 *orig_dims =
            // &var->characteristics[time_index].pre_transform_dimensions;

 
    return 0;
}






#else

DEFINE_FNS_UNIMPL(compress);

#endif
