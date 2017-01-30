/*
 * adios_transform_lz4_write.c
 *
 * 	Author: Rene Widera
 * 	Contact: r.widera@hzdr.de
 *
 *  This code based on `adios_transform_zlib_write.c`.
 */

#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>

#include "core/adios_logger.h"
#include "core/transforms/adios_transforms_common.h"
#include "core/transforms/adios_transforms_write.h"
#include "core/transforms/adios_transforms_hooks_write.h"
#include "core/transforms/adios_transforms_util.h"

#ifdef LZ4

#include <lz4.h>
#include "adios_transform_lz4_common.h"


/** compress data with LZ4
 * 
 * @param stream [in|out] LZ4 stream object
 * @param input_data [in] data to compress
 * @param input_len size of the input data (in byte)
 * @param output_data [out] buffer for the compressed output
 * @param max_output_len size of the output buffer (in byte)
 * @param compressed_size [out] number of bytes after compression
 * @param compress_level LZ4 compression level
 * @return zero (0) if something goes wrong else non zero
 */
int adios_transform_lz4_compress_lz4_compress(LZ4_stream_t* stream,
                                              const char* input_data,
                                              const adiosLz4Size_t input_len,
                                              char* output_data,
                                              const adiosLz4Size_t max_output_len,
                                              adiosLz4Size_t* compressed_size,
                                              int compress_level)
{
    assert(stream != NULL &&
           input_data != NULL &&
           input_len > 0 &&
           output_data != NULL &&
           max_output_len > 0 &&
           compressed_size != NULL);


    adiosLz4Size_t result = LZ4_compress_fast_continue(stream,
                                                       input_data,
                                                       output_data,
                                                       input_len,
                                                       max_output_len,
                                                       compress_level);
    *compressed_size = 0;
    if (result > 0)
    {
        *compressed_size = result;
    }

    return result <= 0;
}

uint16_t adios_transform_lz4_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    /* number of chunks with max size and compressed size of the last chunk 
     * if both are zero the data are uncompressed
     */
    return (sizeof (adiosLz4Size_t) + sizeof (adiosLz4Size_t));
}

/** calculate the maximum data overhead for non compressible data
 * 
 * @param input_size data size in byte of the input data
 * @param n_full_chunks [out] number of full chunks needed to compress the data (NULL is allowed)
 * @param max_size_last_chunks [out] maximum size in bytes needed to compress the last not full chunk (NULL is allowed)
 * @return overhead of bytes for non compressible data
 */
uint64_t adios_transform_lz4_max_overhead(const uint64_t input_size, uint64_t* n_full_chunks, uint64_t* max_size_last_chunks)
{
    const uint64_t max_input_per_chunk = ADIOS_LZ4_MAX_INPUT_SIZE;
    const uint64_t num_full_chunks = input_size / max_input_per_chunk;
    const uint64_t size_last_chunk = input_size % max_input_per_chunk;

    const uint64_t max_output_per_chunk = LZ4_compressBound(ADIOS_LZ4_MAX_INPUT_SIZE);
    const uint64_t max_output_last_chunk = LZ4_compressBound(size_last_chunk);

    if (n_full_chunks)
        *n_full_chunks = num_full_chunks;
    if (max_size_last_chunks)
        *max_size_last_chunks = max_output_last_chunk;

    const uint64_t full_size = max_output_per_chunk * num_full_chunks + max_output_last_chunk;
    return full_size - input_size;
}

void adios_transform_lz4_transformed_size_growth(
                                                 const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
                                                 uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
    struct adios_var_struct *var_no_const = (struct adios_var_struct *) var;
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var_no_const);

    const uint64_t max_overhead = adios_transform_lz4_max_overhead(input_size, NULL, NULL);

    *constant_factor = max_overhead;
}

int adios_transform_lz4_apply(struct adios_file_struct *fd,
                              struct adios_var_struct *var,
                              uint64_t *transformed_len,
                              int use_shared_buffer,
                              int *wrote_to_shared_buffer)
{
    // Assume this function is only called for LZ4 transform type
    assert(var->transform_type == adios_transform_lz4);

    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var);
    const char *input_buff = (char*) var->data;

    int compress_level = 1;
    /* input size under this bound (in byte) would not compressed */
    uint64_t threshold_size = 128;

    int num_param = var->transform_spec->param_count;
    int p;
    for (p = 0; p < num_param; ++p)
    {
        /* Which zfp mode to use */
        const struct adios_transform_spec_kv_pair * const param = &var->transform_spec->params[p];
        if (strcmp(param->key, "lvl") == 0)
        {
            compress_level = atoi(param->value);
            if (compress_level < 1)
                compress_level = 1;
        }
        else if (strcmp(param->key, "threshold") == 0)
        {
            threshold_size = atoi(param->value);
            if (threshold_size < 128)
                threshold_size = 128;
        }
        else
        {
            adios_error(err_invalid_argument, "An unknown LZ4 compression mode '%s' was specified for variable %s. "
                        "Available choices are: lvl, threshold.\n",
                        param->key, var->name);
            return 0;
        }
    }
    // number of full chunks
    uint64_t num_chunks = 0;
    // maximum size for the last not full chunk
    uint64_t max_size_last_chunk = 0;
    const uint64_t max_overhead = adios_transform_lz4_max_overhead(input_size, &num_chunks, &max_size_last_chunk);

    // maximum size (in byte) for the full output buffer
    uint64_t output_size = input_size + max_overhead;
    char* output_buff = NULL;

    if (use_shared_buffer) // If shared buffer is permitted, serialize to there
    {
        *wrote_to_shared_buffer = 1;
        if (!shared_buffer_reserve(fd, output_size))
        {
            log_error("Out of memory allocating %" PRIu64 " bytes for %s for lz4 transform\n", output_size, var->name);
            return 0;
        }

        // Write directly to the shared buffer
        output_buff = (char*) (fd->buffer + fd->offset);
    }
    else // Else, fall back to var->adata memory allocation
    {
        *wrote_to_shared_buffer = 0;
        output_buff = (char*) malloc(output_size);
        if (!output_buff)
        {
            log_error("Out of memory allocating %" PRIu64 " bytes for %s for lz4 transform\n", output_size, var->name);
            return 0;
        }
    }

    // temporary data for lz4 stream compression
    LZ4_stream_t lz4Stream_body = {0};
    LZ4_stream_t* lz4Stream = &lz4Stream_body;

    uint64_t actual_output_size = 0;
    uint64_t input_offset = 0;

    int disable_compression = 0;
    if (input_size < threshold_size)
    {
        /* disable compression */
        disable_compression = 1;
    }

    adiosLz4Size_t compressed_size_last_chunk = 0;
    uint64_t chunk = 0;
    for (; (chunk < num_chunks || input_offset < input_size) && !disable_compression; ++chunk)
    {
        adiosLz4Size_t max_intput_size = ADIOS_LZ4_MAX_INPUT_SIZE;
        // handle the last not full chunk
        if (chunk >= num_chunks)
            max_intput_size = input_size - input_offset;
        adiosLz4Size_t max_chunk_size = LZ4_compressBound(max_intput_size);

        const char* in_ptr = input_buff + input_offset;
        char* out_ptr = output_buff + actual_output_size;

        adiosLz4Size_t compressed_size = 0;
        disable_compression = adios_transform_lz4_compress_lz4_compress(lz4Stream, in_ptr,
                                                                        max_intput_size,
                                                                        out_ptr, max_chunk_size,
                                                                        &compressed_size,
                                                                        compress_level);

        if (chunk >= num_chunks)
            compressed_size_last_chunk = compressed_size;

        /* add size to written output data */
        input_offset += (uint64_t) max_intput_size;
        actual_output_size += (int64_t) compressed_size;
    }

    if (!disable_compression)
        assert(input_offset == input_size);

    if (disable_compression)
    {
        memcpy(output_buff, input_buff, input_size);
        actual_output_size = input_size;
    }

    // Wrap up, depending on buffer mode
    if (use_shared_buffer)
    {
        shared_buffer_mark_written(fd, actual_output_size);
    }
    else
    {
        var->adata = output_buff;
        var->data_size = actual_output_size;
        var->free_data = adios_flag_yes;
    }

    if (var->transform_metadata && var->transform_metadata_len > 0)
    {
        adiosLz4Size_t n_chunks = num_chunks;

        if (disable_compression)
        {
            n_chunks = 0;
            compressed_size_last_chunk = 0;
        }
        /* store number of full chunks */
        memcpy((char*) var->transform_metadata, &n_chunks, sizeof (adiosLz4Size_t));
        /* size (in byte) of the last not full chunk */
        memcpy((char*) var->transform_metadata + sizeof (adiosLz4Size_t), &compressed_size_last_chunk, sizeof (adiosLz4Size_t));
    }

    // return the size of the data buffer
    *transformed_len = actual_output_size; 
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(lz4)

#endif
