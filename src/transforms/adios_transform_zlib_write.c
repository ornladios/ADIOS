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

int compress_zlib_pre_allocated(const void* input_data, const uint64_t input_len,
                                void* output_data, uint64_t* output_len, int compress_level)
{
    assert(input_data != NULL && input_len > 0 && output_data != NULL && output_len != NULL && *output_len > 0);

    uLongf temp = *output_len;
    int zerr = compress2((Bytef*)(output_data), &temp,
                            (Bytef*)input_data, (uLongf)input_len,
                            compress_level);

    *output_len = (uint64_t)temp;

    if(zerr != Z_OK)
    {
        // printf("zlib compress2 error %d\n", zerr);
        return -1;
    }
    return 0;
}

uint16_t adios_transform_zlib_get_metadata_size()
{
    return 0;
}

uint64_t adios_transform_zlib_calc_vars_transformed_size(uint64_t orig_size, int num_vars)
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
    const void *input_buff = var->data;

    // parse the compressiong parameter
    int compress_level = Z_DEFAULT_COMPRESSION;
    if(var->transform_type_param
        && strlen(var->transform_type_param) > 0)
    {
        compress_level = atoi(var->transform_type_param);
        if(compress_level > 9 || compress_level < 1)
        {
            compress_level = Z_DEFAULT_COMPRESSION;
        }
    }

    // decide the output buffer
    uint64_t output_size = adios_transform_zlib_calc_vars_transformed_size(input_size, 1);
    void* output_buff = NULL;

    if (use_shared_buffer) {
        // If shared buffer is permitted, serialize to there
        assert(shared_buffer_reserve(fd, output_size));

        // Write directly to the shared buffer
        output_buff = fd->buffer + fd->offset;
    } else { // Else, fall back to var->data memory allocation
        output_buff = malloc(output_size);
        assert(output_buff);
    }
    *wrote_to_shared_buffer = use_shared_buffer;

    // compress it
    int rtn = compress_zlib_pre_allocated(input_buff, input_size, output_buff, &output_size, compress_level);

    if(0 != rtn 					// compression failed for some reason, then just copy the buffer
        || output_size > input_size)  // or size after compression is even larger (not likely to happen since compression lib will return non-zero in this case)
    {
        return 0;
    }

    // Wrap up, depending on buffer mode
    if (*wrote_to_shared_buffer) {
        shared_buffer_mark_written(fd, output_size);
    } else {
        var->data = output_buff;
        var->data_size = output_size;
        var->free_data = adios_flag_yes;
    }

    *transformed_len = output_size; // Return the size of the data buffer
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(zlib)

#endif

