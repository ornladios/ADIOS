#include <stdint.h>
#include <assert.h>
#include <limits.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"

#ifdef ALACRITY

#include "alacrity.h"
#include "alacrity-serialization-debug.h"

#define TEST_SIZE(s) (1.75 * s)

static int is_digit_str(char* input_str)
{
	return 1;
}

uint16_t adios_transform_alacrity_get_metadata_size() 
{ 
	return (3 * sizeof(uint64_t));
}

uint64_t adios_transform_alacrity_calc_vars_transformed_size(uint64_t orig_size, int num_vars) 
{
    return TEST_SIZE(orig_size);
}

int adios_transform_alacrity_apply(struct adios_file_struct *fd, 
								struct adios_var_struct *var,
								uint64_t *transformed_len, 
								int use_shared_buffer, 
								int *wrote_to_shared_buffer)
{
    // Assume this function is only called for ALACRITY transform type
    assert(var->transform_type == adios_transform_alacrity);

    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input_buff= var->data;

    ALEncoderConfig config;
    uint32_t numElements = 0;

    if (var->pre_transform_type == adios_real) {
        assert (sizeof (DATATYPE_FLOAT32) == sizeof (adios_real));
        ALEncoderConfigure(&config, 16, DATATYPE_FLOAT32, ALInvertedIndex);
        numElements = input_size / sizeof (float);
    } else if (var->pre_transform_type == adios_double) {
        assert (sizeof (DATATYPE_FLOAT64) == sizeof (adios_double));
        ALEncoderConfigure(&config, 16, DATATYPE_FLOAT64, ALInvertedIndex);
        numElements = input_size / sizeof (double);
    } else {
        // What should I do to ensure that the transform actually does nothing
        // Check Error
        log_error("Can index only real datatypes. \n");

        return 0;
    }

	// parse the parameter relating to sigbits, with index compression
    // Read all ALACRITY parameters here
	if(var->transform_type_param)
	{
        char transform_param [1024];
        char *transform_param_ptr       = 0;
        uint16_t transform_param_length = 0;

        char transform_param_option [256];

        uint16_t idx = 0;

        strcpy (transform_param, var->transform_type_param);
        transform_param_ptr     = transform_param;
        transform_param_length  = strlen (transform_param);

        // Change all the delimiters to a space
        while (idx < transform_param_length) { 
            if (transform_param [idx] == ':') { 
                transform_param [idx] = ' ';       
            }
            idx ++;
        } 

        // For each option in the transform parameter,
        // get its key and possibly its value.
        idx = 0;
        while (idx < transform_param_length) {
            // Get the first key, value pair
            sscanf (transform_param_ptr, "%s", transform_param_option);

            // Advance the pointer
            idx += strlen (transform_param_option) + 1;  
            transform_param_ptr = transform_param + idx;

            // Get the key
            char *key = strtok (transform_param_option, "="); 

            if (strcmp (key, "indexForm") == 0) {
                char *value = strtok (NULL, "=");
                if (strcmp (value, "ALCompressedInvertedIndex") == 0) { 
                    config.indexForm = ALCompressedInvertedIndex;
                } else if (strcmp (value, "ALInvertedIndex") == 0) {
                    config.indexForm = ALInvertedIndex;    
                }

            } else if (strcmp (key, "sigBits") == 0) {
                char *value = strtok (NULL, "=");
                int significantBits = 0;
                sscanf (value, "%d", &(significantBits));
                config.significantBits = significantBits;
            } else {
                printf ("Option %s not found. \n", key);
            }
        }        
	}

    // decide the output buffer
    uint64_t output_size = 0;
    void* output_buff = NULL;
    ALPartitionData output_partition;

    uint64_t mem_allowed = 0;

    uint8_t indexed = false;

    if (config.indexForm == ALInvertedIndex) {
        ALEncode (&config, input_buff, numElements, &output_partition);
        output_size = ALGetPartitionDataSize (&output_partition);
        indexed = true;
    } else if (config.indexForm == ALCompressedInvertedIndex) {
        config.indexForm = ALInvertedIndex;

        ALBinLookupTable binLookupTable; 
        ALBuildBinLayout (&config, input_buff, numElements, &binLookupTable, &output_partition);

        ALBuildInvertedIndexFromLayout (&config, input_buff, numElements, &binLookupTable, &output_partition);
        ALConvertIndexForm (&output_partition.metadata, &output_partition.index, ALCompressedInvertedIndex);

        output_size = ALGetPartitionDataSize (&output_partition);
        indexed = true;
    } 
    
    if (!indexed) {
		*wrote_to_shared_buffer = 0;
		output_buff = 0;
        log_error("Error on indexing %s. Reverting to original data\n", var->name);

        return 0;
    }

    if (use_shared_buffer)
    {	
		*wrote_to_shared_buffer = 1;
		// If shared buffer is permitted, serialize to there
        if (!shared_buffer_reserve(fd, output_size))
        {
            log_error("Out of memory allocating %llu bytes for %s for ALACRITY transform\n", output_size, var->name);
            return 0;
        }

        // Write directly to the shared buffer
        output_buff = fd->buffer + fd->offset;
    } else {
		*wrote_to_shared_buffer = 0;
        output_buff = malloc (output_size);
    }

    // Else, fall back to var->data memory allocation
    memstream_t ms;
    memstreamInit (&ms, output_buff);

    if (ms.buf == 0) {
        log_error("Out of memory allocating %llu bytes for %s for ALACRITY transform\n", output_size, var->name);
        ALPartitionDataDestroy (&output_partition);
        return 0;
    }

    ALSerializePartitionData (&output_partition, &ms);

    // Check this for later. What do you intend to add in the metadata
    if(var->transform_metadata && var->transform_metadata_len > 0) {
        ((uint64_t * ) (var->transform_metadata)) [0] = ALGetMetadataSize (& (output_partition.metadata));
        ((uint64_t * ) (var->transform_metadata)) [1] = ALGetIndexSize (& (output_partition.index), & (output_partition.metadata));
        ((uint64_t * ) (var->transform_metadata)) [2] = ALGetDataSize  (& (output_partition.data), & (output_partition.metadata));
    }

    ALPartitionDataDestroy (&output_partition);
    memstreamDestroy (&ms, false);

    // if ( actual_output_size > input_size)  // or size after compression is even larger (not likely to happen since compression lib will return non-zero in this case)
    // {
    //     memcpy(output_buff, input_buff, input_size);
    //     actual_output_size = input_size;
    // }

    // Wrap up, depending on buffer mode
    if (*wrote_to_shared_buffer)
    {
        shared_buffer_mark_written(fd, output_size);
    }
    else
    {
        var->data = output_buff;
        var->data_size = output_size;
        var->free_data = adios_flag_yes;
    }

	
	// printf("adios_transform_compress_apply compress %d input_size %d actual_output_size %d\n", 
			// rtn, input_size, actual_output_size);
	
    *transformed_len = output_size; // Return the size of the data buffer
	
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(alacrity)

#endif

