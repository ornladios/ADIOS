
#include <assert.h>

#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "public/adios_error.h"
#include "public/adios_types.h"

#include "adios_transforms_hooks.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_util.h"

enum ADIOS_ERRCODES adios_transform_retrieve_transformed_data_subvolume(
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *copy_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness) {
    assert(var->characteristics_count > 0);

    enum ADIOS_TRANSFORM_TYPE transform_type = var->characteristics[0].transform.transform_type;
    assert(transform_type != adios_transform_none);

    return adios_transform_retrieve_subvolume(transform_type, var, global_out, time_index, copy_spec, read_state, read_delegate, swap_endianness);
}
