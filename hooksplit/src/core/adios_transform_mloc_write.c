#ifdef MLOC

#include <stdint.h>
#include <assert.h>
#include "adios_logger.h"
#include "adios_transforms_hooks.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_write.h"
#include "mloc.h"

int adios_transform_mloc_apply(struct adios_file_struct *fd, struct adios_var_struct *var, uint64_t *transformed_len, int *use_shared_buffer) {
    // Assume this function is only called for MLOC transform type
    assert(var->transform_type == adios_transform_mloc);

    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input = var->data;

    
    return 1;
}

#endif
