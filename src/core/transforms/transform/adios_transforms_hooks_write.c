
#include <stdint.h>
#include <assert.h>

#include "adios_logger.h"
#include "adios_internals.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"
#include "public/adios_selection.h"

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(none);
DECLARE_TRANSFORM_WRITE_METHOD(identity);
DECLARE_TRANSFORM_WRITE_METHOD(alacrity);

// Transform write method registry
adios_transform_write_method TRANSFORM_WRITE_METHODS[num_adios_transform_types];

void adios_transform_init() {
    static int adios_transforms_initialized = 0;
    if (adios_transforms_initialized)
        return;

    REGISTER_TRANSFORM_WRITE_METHOD(none, adios_transform_none);
    REGISTER_TRANSFORM_WRITE_METHOD(identity, adios_transform_identity);
    REGISTER_TRANSFORM_WRITE_METHOD(alacrity, adios_transform_alacrity);

    adios_transforms_initialized = 1;
}

// Delegate functions

uint16_t adios_transform_get_metadata_size(enum ADIOS_TRANSFORM_TYPE transform_type) {
    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_WRITE_METHODS[transform_type].transform_get_metadata_size();
}

uint64_t adios_transform_calc_vars_transformed_size(enum ADIOS_TRANSFORM_TYPE transform_type, uint64_t orig_size, int num_vars) {
    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_WRITE_METHODS[transform_type].transform_calc_vars_transformed_size(orig_size, num_vars);
}

int adios_transform_apply(
        enum ADIOS_TRANSFORM_TYPE transform_type,
        struct adios_file_struct *fd, struct adios_var_struct *var,
        uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer) {

    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_WRITE_METHODS[transform_type].transform_apply(fd, var, transformed_len, use_shared_buffer, wrote_to_shared_buffer);
}

// Implementation of the "identity" transform, which does nothing, but
// exercises the transform framework for testing.
// (the other identity method hooks are in adios_transforms_hooks.c)

uint64_t adios_transform_identity_calc_vars_transformed_size(uint64_t orig_size, int num_vars) {
    return orig_size;
}

uint16_t adios_transform_identity_get_metadata_size() {
    return 0;
}

int adios_transform_identity_apply(struct adios_file_struct *fd, struct adios_var_struct *var, uint64_t *transformed_len,
                                   int use_shared_buffer, int *wrote_to_shared_buffer) {
    // Just use what is already in var->data; size remains the same, and no
    // shared buffer is used
    *transformed_len = adios_transform_get_pre_transform_var_size(fd->group, var);
    *wrote_to_shared_buffer = 0;
    return 1;
}
