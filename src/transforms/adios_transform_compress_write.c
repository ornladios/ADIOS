#include <stdint.h>
#include <assert.h>
#include <limits.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"

#ifdef COMPRESS
#include "compress.h"

const struct
{
    const char *alias;				// A possible name for a compres method
    enum COMPRESS_TYPE type;	// The corresponding COMPRESS_TYPE
} COMPRESS_TYPE_NAMES[] = {
     { "unknown"	, compress_type_unknown }

    ,{ "none"		, compress_type_none }
    ,{ "no"			, compress_type_none }
    ,{ "raw"		, compress_type_none }
    ,{ ""			, compress_type_none }

    ,{ "zlib"	, compress_type_zlib }
    ,{ "zip"	, compress_type_zlib }

    ,{ "bzlib2"	, compress_type_bzlib2 }
    ,{ "bzip"	, compress_type_bzlib2 }

    ,{ "szlib"	, compress_type_szip }
    ,{ "szip"	, compress_type_szip }
};

const int NUM_COMPRESS_TYPE_NAMES = sizeof(COMPRESS_TYPE_NAMES)/sizeof(COMPRESS_TYPE_NAMES[0]);


uint16_t adios_transform_compress_get_metadata_size() { return sizeof(enum COMPRESS_TYPE); }

uint64_t adios_transform_compress_calc_vars_transformed_size(uint64_t orig_size, int num_vars) {
    return EXPAND_SIZE(orig_size);
}

int adios_transform_compress_apply(struct adios_file_struct *fd, struct adios_var_struct *var,
                                uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer) {
    // Assume this function is only called for COMPRESS transform type
    assert(var->transform_type == adios_transform_compress);

    // Get the input data and data length
    const uint64_t input_size_64 = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input_buff= var->data;

    if(input_size_64 > UINT_MAX) // too large, do not support 64 bit now
    {
        log_error("data size too long %llu bytes for %s for compress transform\n", input_size_64, var->name);
        return 0;
    }

    const uint32_t input_size = (uint32_t)input_size_64;

    // parse the parameter to know the compressiong type
    enum COMPRESS_TYPE compress_type_flag = compress_type_zlib; // using zlib by default

    if(var->transform_type_param && var->transform_type_param[0] != '\0')
    {
        int i = 0;
        for (i = 0; i < NUM_COMPRESS_TYPE_NAMES; i++)
        {
            if (strcasecmp(var->transform_type_param, COMPRESS_TYPE_NAMES[i].alias) == 0)
            {
                compress_type_flag = COMPRESS_TYPE_NAMES[i].type;
                break;
            }
        }
    }

    // decide the output buffer
    uint32_t output_size = EXPAND_SIZE(input_size);
    void* output_buff = NULL;

    uint64_t mem_allowed = 0;
    if (use_shared_buffer)
    {	// If shared buffer is permitted, serialize to there
        if (!shared_buffer_reserve(fd, output_size))
        {
            log_error("Out of memory allocating %llu bytes for %s for compress transform\n", output_size, var->name);
            return 0;
        }

        // Write directly to the shared buffer
        output_buff = fd->buffer + fd->offset;
    }
    else
    {					// Else, fall back to var->data memory allocation

        mem_allowed = adios_method_buffer_alloc(output_size);
        // If we aren't allowed enough memory, or the malloc fails...
        if (mem_allowed != output_size ||
            !(output_buff = malloc(output_size)) )
        {
            adios_method_buffer_free(mem_allowed);
            log_error("Out of memory allocating %llu bytes for %s for ALACRITY transform\n", output_size, var->name);
            return 0;
        }
    }

    // compress it
    int rtn = 0;

    switch(compress_type_flag)
    {
        case compress_type_zlib:
            rtn = compress_zlib_pre_allocated(input_buff, input_size, output_buff, &output_size);
        break;

        case compress_type_bzlib2:
            rtn = compress_bzlib2_pre_allocated(input_buff, input_size, output_buff, &output_size);
        break;

        case compress_type_szip:
            // rtn = compress_szip_pre_allocated(input_buff, input_size, output_buff, &output_size);
            rtn = compress_bzlib2_pre_allocated(input_buff, input_size, output_buff, &output_size);
        break;

        default: // default: do not do compression, just copy the buffer
		{
            memcpy(output_buff, input_buff, input_size);
            output_size = input_size;
            compress_type_flag = compress_type_none;
            rtn = 0;
		}
        break;

    }

    if(0 != rtn 					// compression failed for some reason, then just copy the buffer
        || output_size > input_size)  // or size after compression is even larger
    {
        memcpy(output_buff, input_buff, input_size);
        output_size = input_size;
        compress_type_flag = compress_type_none;
    }

    // Wrap up, depending on buffer mode
    if (use_shared_buffer)
    {
        shared_buffer_mark_written(fd, output_size);
    }
    else
    {
        var->data = output_buff;
        var->data_size = output_size;
        var->free_data = adios_flag_yes;
    }

    // copy the metadata, simply the compress type
    if(var->transform_metadata_len)
    {
        memcpy(var->transform_metadata, &compress_type_flag, var->transform_metadata_len);
    }
	
	printf("adios_transform_compress_apply %d %d\n", input_size, output_size);


    *transformed_len = output_size; // Return the size of the data buffer
    return 1;
}

// uint64_t adios_transform_compress_vars_size(uint64_t orig_size, int num_vars) {
    // return num_vars * sizeof(double) +	// For the metadata
           // orig_size * 5/4;						// For the index + data
// }

#endif
