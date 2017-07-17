/*
 * adios_transform_blosc_write.c
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

#ifdef BLOSC

#include <blosc.h>
#include "adios_transform_blosc_common.h"


int adios_transform_blosc_compress( 
                                const void* input_data,
                                const adiosBloscSize_t input_len,
                                void* output_data,
                                const adiosBloscSize_t max_output_len,
                                adiosBloscSize_t* compressed_size,
                                int compress_level,
                                int shuffle, int typesize
                                )
{
    assert( input_data != NULL && input_len > 0 && output_data != NULL && max_output_len > 0 && compressed_size != NULL);

    
    adiosBloscSize_t result = blosc_compress( 
        compress_level, shuffle, typesize, 
        input_len, input_data,
        output_data, 
        max_output_len
    );
    *compressed_size = 0;
    if( result > 0)
    {
        *compressed_size = result;
    }
    
    return result <= 0;
}

uint16_t adios_transform_blosc_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    // number of chunks with max size and compressed size of the last chunk
    return (sizeof(adiosBloscSize_t) + sizeof(adiosBloscSize_t));
}

uint64_t calculate_max_overhead(const uint64_t input_size, uint64_t* n_full_chunks, uint64_t* max_size_last_chunks)
{   
    const uint64_t max_input_per_chunk = ADIOS_BLOSC_MAX_INPUT_SIZE;
    const uint64_t num_full_chunks = input_size / max_input_per_chunk;
    const uint64_t size_last_chunk = input_size % max_input_per_chunk;
    
    const uint64_t max_output_per_chunk = ADIOS_BLOSC_MAX_INPUT_SIZE + ADIOS_BLOSC_MAX_OVERHEAD;
    const uint64_t max_output_last_chunk = size_last_chunk + ADIOS_BLOSC_MAX_OVERHEAD;
    
    if(n_full_chunks)
        *n_full_chunks = num_full_chunks;
    if(max_size_last_chunks)
        *max_size_last_chunks = max_output_last_chunk;
    
    const uint64_t full_size = max_output_per_chunk * num_full_chunks + max_output_last_chunk;
    return full_size - input_size;
}

void adios_transform_blosc_transformed_size_growth(
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
    struct adios_var_struct *var_no_const = (struct adios_var_struct *)var;
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var_no_const);

    const uint64_t max_overhead = calculate_max_overhead(input_size, NULL, NULL);
    
    *constant_factor = max_overhead;
}

int adios_transform_blosc_apply(struct adios_file_struct *fd,
                                struct adios_var_struct *var,
                                uint64_t *transformed_len,
                                int use_shared_buffer,
                                int *wrote_to_shared_buffer)
{
    // Assume this function is only called for Blosc transform type
    assert(var->transform_type == adios_transform_blosc);

    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var);
    const char *input_buff= (char*)var->data;
    
    uint64_t num_chunks = 0;
    uint64_t max_size_last_chunk = 0;
    const uint64_t max_overhead = calculate_max_overhead(input_size, &num_chunks, &max_size_last_chunk);
   
    int compress_level = 1;
    int shuffle = BLOSC_NOSHUFFLE;
    int num_threads = 1;
    char compressor[32];
    strcpy( compressor, "memcpy" );
    /* input size under this bound would not compressed */
    uint64_t threshold_size = 128;
    
    int typesize = adios_get_type_size( var->pre_transform_type, NULL );
    if( typesize > ADIOS_BLOSC_MAX_TYPESIZE )
        typesize = 1;

    int num_param = var->transform_spec->param_count;
    int p;
    for( p = 0; p < num_param; ++p )
    {
        const struct adios_transform_spec_kv_pair* const param = &var->transform_spec->params[p];
        if( strcmp( param->key, "lvl") == 0 )
        {
            compress_level = atoi( param->value );
            if( compress_level < 1 || compress_level > 9 )
            {
                if( compress_level < 1 )
                    compress_level = 1;
                if( compress_level > 9 )
                    compress_level = 9;
                log_warn("Blosc: invalid compression level %s, switch to lvl %i\n",param->value, compress_level);
            }
        }
        else if( strcmp( param->key, "threshold") == 0 )
        {
            threshold_size = atoi( param->value );
            if (threshold_size < 128)
                threshold_size = 128;
        }
        else if( strcmp( param->key, "shuffle") == 0 )
        {
            if (strcmp(param->value, "byte") == 0) 
               shuffle = BLOSC_SHUFFLE;
            else if (strcmp(param->value, "bit") == 0)
               shuffle = BLOSC_BITSHUFFLE;
            else if (strcmp(param->value, "no") == 0) 
               shuffle = BLOSC_NOSHUFFLE;
            else
            {
                log_warn("Blosc: unknown shuffle %s, disable shuffle\n",param->value);
            }
        }
        else if( strcmp(param->key, "threads") == 0 )
        {
            num_threads = atoi( param->value );
            if( num_threads < 1 || num_threads > 8)
                num_threads = 1;
        }
        else if( strcmp( param->key, "compressor") == 0 )
        {
            strcpy( compressor, param->value );
            if( strcmp( compressor, "memcpy") != 0 && blosc_set_compressor( compressor ) < 0 )
            {
                log_warn("Blosc: unknown compressor %s, switch to memcpy\n",param->value);
                strcpy( compressor, "memcpy" );
            }
        }
        else
        {
            adios_error(err_invalid_argument, "An unknown Blosc compression parameter '%s' was specified for variable %s.\n",
                        param->key, var->name);
            return 0;
        }
    }

    
    // decide the output buffer
    uint64_t output_size = input_size + max_overhead; // for compression, at most the original data size
    char* output_buff = NULL;

    if (use_shared_buffer)    // If shared buffer is permitted, serialize to there
    {
        *wrote_to_shared_buffer = 1;
        if (!shared_buffer_reserve(fd, output_size))
        {
            log_error("Out of memory allocating %" PRIu64 " bytes for %s for blosc transform\n", output_size, var->name);
            return 0;
        }

        // Write directly to the shared buffer
        output_buff = (char*)( fd->buffer + fd->offset );
    }
    else    // Else, fall back to var->adata memory allocation
    {
        *wrote_to_shared_buffer = 0;
        output_buff = (char*)malloc(output_size);
        if (!output_buff)
        {
            log_error("Out of memory allocating %" PRIu64 " bytes for %s for blosc transform\n", output_size, var->name);
            return 0;
        }
    }

    uint64_t actual_output_size = 0;
    uint64_t input_offset = 0;
    int compress_failed = 0;
    
    if( input_size < threshold_size )
    {
        /* disable compression */
        compress_failed = 1;
    }
    
    if( strcmp( compressor, "memcpy") == 0 )
    {
        char* envvar = getenv("BLOSC_COMPRESSOR");
        if (envvar == NULL)
            compress_failed = 1;
    }
        
    blosc_set_nthreads(num_threads);
    adiosBloscSize_t compressed_size_last_chunk = 0;
    uint64_t chunk = 0;
    for(; (chunk < num_chunks || input_offset < input_size) && !compress_failed; ++chunk)
    {
        adiosBloscSize_t max_intput_size = ADIOS_BLOSC_MAX_INPUT_SIZE;
        // handle the last not full chunk
        if (chunk >= num_chunks)
            max_intput_size = input_size - input_offset;
        adiosBloscSize_t max_chunk_size = max_intput_size + ADIOS_BLOSC_MAX_OVERHEAD;
       
        const char* in_ptr = input_buff + input_offset;
        char* out_ptr = output_buff + actual_output_size;

        adiosBloscSize_t compressed_size = 0;
        compress_failed = adios_transform_blosc_compress(in_ptr, max_intput_size,
                                          out_ptr, max_chunk_size,&compressed_size,
                                          compress_level, shuffle, typesize );
        
        if (chunk >= num_chunks)
            compressed_size_last_chunk = compressed_size;

        /* add size to written output data */
        input_offset += (uint64_t) max_intput_size;
        actual_output_size += (int64_t) compressed_size;
    }
    
    if (!compress_failed)
        assert(input_offset == input_size);
    
    if(compress_failed)
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

    // copy the metadata, simply the original size before compression
    if(var->transform_metadata && var->transform_metadata_len > 0)
    {
        adiosBloscSize_t n_chunks = num_chunks;
        
        if(compress_failed == 1)
        {
            n_chunks = 0;
            compressed_size_last_chunk = 0;
        }
        /* store size of the chunk behind the bit mask */
        memcpy((char*)var->transform_metadata, &n_chunks, sizeof(adiosBloscSize_t));
        /* store the bit mask */
        memcpy((char*)var->transform_metadata + sizeof(adiosBloscSize_t), &compressed_size_last_chunk, sizeof(adiosBloscSize_t));
    }

    *transformed_len = actual_output_size; // Return the size of the data buffer
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(blosc)

#endif
