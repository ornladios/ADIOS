#include <stdint.h>
#include <assert.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"

#ifdef MLOC
#include "mloc.h"

uint16_t adios_transform_mloc_get_metadata_size() { return 0; }

#define MAX_POSSIBLE_BINS 65536
static const int MAX_PART_METADATA_SIZE =
    sizeof(uint64_t) +
    sizeof(unsigned short int) +
    MAX_POSSIBLE_BINS * ( sizeof(unsigned short int) +  2 * sizeof(uint64_t) + sizeof(unsigned char) );


uint64_t adios_transform_mloc_calc_vars_transformed_size(uint64_t orig_size, int num_vars) {
    return num_vars * MAX_PART_METADATA_SIZE +	// For the metadata
           orig_size * 5/4;						// For the index + data
}

int adios_transform_mloc_apply(struct adios_file_struct *fd, struct adios_var_struct *var, uint64_t *transformed_len, int *use_shared_buffer) {
    // Assume this function is only called for MLOC transform type
    assert(var->transform_type == adios_transform_mloc);

    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input = var->data;


    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(mloc)

#endif
