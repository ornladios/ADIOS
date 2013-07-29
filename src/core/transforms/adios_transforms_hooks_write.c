
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
DECLARE_TRANSFORM_WRITE_METHOD(zlib);
DECLARE_TRANSFORM_WRITE_METHOD(bzip2);
DECLARE_TRANSFORM_WRITE_METHOD(szip);
DECLARE_TRANSFORM_WRITE_METHOD(isobar);
DECLARE_TRANSFORM_WRITE_METHOD(aplod);
DECLARE_TRANSFORM_WRITE_METHOD(alacrity);

// Transform write method registry
adios_transform_write_method TRANSFORM_WRITE_METHODS[num_adios_transform_types];

void adios_transform_init() {
    static int adios_transforms_initialized = 0;
    if (adios_transforms_initialized)
        return;

    REGISTER_TRANSFORM_WRITE_METHOD(none, adios_transform_none);
    REGISTER_TRANSFORM_WRITE_METHOD(identity, adios_transform_identity);
    REGISTER_TRANSFORM_WRITE_METHOD(zlib, adios_transform_zlib);
    REGISTER_TRANSFORM_WRITE_METHOD(bzip2, adios_transform_bzip2);
    REGISTER_TRANSFORM_WRITE_METHOD(szip, adios_transform_szip);
    REGISTER_TRANSFORM_WRITE_METHOD(isobar, adios_transform_isobar);
    REGISTER_TRANSFORM_WRITE_METHOD(aplod, adios_transform_aplod);
    REGISTER_TRANSFORM_WRITE_METHOD(alacrity, adios_transform_alacrity);

    adios_transforms_initialized = 1;
}

// Delegate functions

uint16_t adios_transform_get_metadata_size(struct adios_transform_spec *transform_spec) {
    assert(transform_spec->transform_type >= adios_transform_none && transform_spec->transform_type < num_adios_transform_types);
    return TRANSFORM_WRITE_METHODS[transform_spec->transform_type].transform_get_metadata_size(transform_spec);
}

uint64_t adios_transform_calc_vars_transformed_size(enum ADIOS_TRANSFORM_TYPE transform_type, uint64_t orig_size, int num_vars) {
    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_WRITE_METHODS[transform_type].transform_calc_vars_transformed_size(transform_type, orig_size, num_vars);
}

int adios_transform_apply(
        struct adios_file_struct *fd, struct adios_var_struct *var,
        uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer) {

    assert(var->transform_type >= adios_transform_none && var->transform_type < num_adios_transform_types);
    return TRANSFORM_WRITE_METHODS[var->transform_type].transform_apply(fd, var, transformed_len, use_shared_buffer, wrote_to_shared_buffer);
}
