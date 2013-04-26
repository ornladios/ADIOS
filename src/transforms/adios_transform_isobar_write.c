#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"

#ifdef ISOBAR

#include "isobar.h"

#define ELEMENT_BYTES	8

int compress_isobar_pre_allocated(const void* input_data, const uint64_t input_len,
                                  void* output_data, uint64_t* output_len, int compress_level)
{
    assert(input_data != NULL && input_len > 0 && output_data != NULL && output_len != NULL && *output_len > 0);

    enum ISOBAR_status status;
    struct isobar_stream i_strm;

    status = isobarDeflateInit(&i_strm, ELEMENT_BYTES, compress_level);
    if(status != ISOBAR_SUCCESS)
    {
        return -1;
    }

    i_strm.next_in = (void*) input_data;
    i_strm.avail_in = input_len;

    status = isobarDeflateAnalysis(&i_strm);
    if(status != ISOBAR_SUCCESS)
    {
        return -1;
    }

    i_strm.next_out = (void*) output_data;
    i_strm.avail_out = input_len;

    status = isobarDeflate(&i_strm, ISOBAR_FINISH);
    if(status != ISOBAR_SUCCESS)
    {
        return -1;
    }

    status = isobarDeflateEnd(&i_strm);
    if(status != ISOBAR_SUCCESS)
    {
        return -1;
    }

    *output_len = input_len - i_strm.avail_out;
    return 0;
}

uint16_t adios_transform_isobar_get_metadata_size()
{
    return 0;
}

uint64_t adios_transform_isobar_calc_vars_transformed_size(uint64_t orig_size, int num_vars)
{
    return orig_size;
}

int adios_transform_isobar_apply(struct adios_file_struct *fd,
                                 struct adios_var_struct *var,
                                 uint64_t *transformed_len,
                                 int use_shared_buffer,
                                 int *wrote_to_shared_buffer)
{
    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input_buff = var->data;

    // parse the compressiong parameter
    int compress_level = ISOBAR_SPEED;
    if(var->transform_type_param
        && strlen(var->transform_type_param) > 0)
    {
        compress_level = atoi(var->transform_type_param);
        if(compress_level > 9 || compress_level < 1)
        {
            compress_level = ISOBAR_SPEED;
        }
    }

    // decide the output buffer
    uint64_t output_size = adios_transform_isobar_calc_vars_transformed_size(input_size, 1);
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
    int rtn = compress_isobar_pre_allocated(input_buff, input_size, output_buff, &output_size, compress_level);

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

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(isobar)

#endif

