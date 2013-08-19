#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"

#ifdef ZLIB

#include "zlib.h"

uint16_t adios_transform_zlib_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    return sizeof(char);    // metadata: compression success flag (char)
}

uint64_t adios_transform_zlib_calc_vars_transformed_size(enum ADIOS_TRANSFORM_TYPE type, uint64_t orig_size, int num_vars)
{
    return orig_size;
}

int adios_transform_zlib_apply(struct adios_file_struct *fd,
                               struct adios_var_struct *var,
                               uint64_t *transformed_len,
                               int use_shared_buffer,
                               int *wrote_to_shared_buffer)
{
    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input_buff= var->data;

    // Parse the parameter
    int compress_level = Z_DEFAULT_COMPRESSION;
    if (var->transform_spec->param_count > 0) {
        compress_level = atoi(var->transform_spec->params[0].key);
        if (compress_level < 1 || compress_level > 9)
            compress_level = Z_DEFAULT_COMPRESSION;
    }

    // Determine whether to use the ADIOS shared buffer
    uint64_t output_size = input_size; // for compression, at most the original data size
    void *output_buff = NULL;

    // If shared buffer is permitted, use it
    // Else, fall back to allocating a new buffer
    if (use_shared_buffer) {
        // Reserve space
        if (shared_buffer_reserve(fd, output_size)) {
            // Set to write directly to the shared buffer
            *wrote_to_shared_buffer = 1;
            output_buff = fd->buffer + fd->offset;
        }
    } else {
        *wrote_to_shared_buffer = 0;
        output_buff = malloc(output_size);
        // If out of memory, output_buff is NULL
    }

    // Ensure we actually got a buffer
    if (!output_buff) {
        log_error("Out of memory allocating %llu bytes for %s for zlib transform\n", output_size, var->name);
        return 0;
    }

    // Apply compression
    int rtn = compress2(output_buff, &output_size, input_buff, input_size, compress_level);
    char compression_successful = (rtn == Z_OK);

    // If compression failed, just copy the original data
    if (!compression_successful) {
        memcpy(output_buff, input_buff, input_size);
        output_size = input_size;
    }

    // Wrap up, depending on buffer mode
    if (use_shared_buffer) {
        shared_buffer_mark_written(fd, output_size);
    } else {
        var->data = output_buff;
        var->data_size = output_size;
        var->free_data = adios_flag_yes;
    }

    // Store whether the compression was actually applied
    *(char*)var->transform_metadata = compression_successful;

    // Return the size of the data buffer
    *transformed_len = output_size;
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(zlib)

#endif

