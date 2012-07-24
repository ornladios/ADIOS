
#include <stdint.h>
#include <assert.h>

#include "adios_logger.h"
#include "adios_internals.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks.h"
#include "public/adios_selection.h"

/*
 * The transform method registry, containing all necessary function pointers
 */
extern adios_transform_method TRANSFORM_METHODS[num_adios_transform_types];
extern int adios_transforms_initialized;

// We just define the write version of the init function here, since it
// references the *apply functions, which are not linked with read
// libraries
void adios_transform_init() {
    if (adios_transforms_initialized)
        return;

    ASSIGN_FNS(none, adios_transform_none);
    ASSIGN_FNS(identity, adios_transform_identity);
    ASSIGN_FNS(alacrity, adios_transform_alacrity);

    adios_transforms_initialized = 1;
}

uint16_t adios_transform_get_metadata_size(enum ADIOS_TRANSFORM_TYPE transform_type) {
    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_METHODS[transform_type].transform_get_metadata_size();
}

int adios_transform_apply(
        enum ADIOS_TRANSFORM_TYPE transform_type,
        struct adios_file_struct *fd, struct adios_var_struct *var,
        uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer) {

    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_METHODS[transform_type].transform_apply(fd, var, transformed_len, use_shared_buffer, wrote_to_shared_buffer);
}

// Implementation of the "identity" transform, which does nothing, but
// exercises the transform framework for testing.
// (the other identity method hooks are in adios_transforms_hooks.c)

uint16_t adios_transform_identity_get_metadata_size() {
    return 24;
}

int adios_transform_identity_apply(struct adios_file_struct *fd, struct adios_var_struct *var, uint64_t *transformed_len,
                                   int use_shared_buffer, int *wrote_to_shared_buffer) {
    // Just use what is already in var->data; size remains the same, and no
    // shared buffer is used
    *transformed_len = adios_transform_get_pre_transform_var_size(fd->group, var);
    *wrote_to_shared_buffer = 0;
    return 1;
}
