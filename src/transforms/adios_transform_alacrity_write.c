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

uint16_t adios_transform_alacrity_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    return (3 * sizeof(uint64_t));
}

uint64_t adios_transform_alacrity_calc_vars_transformed_size(struct adios_transform_spec *transform_spec, uint64_t orig_size, int num_vars)
{
    return (uint64_t)(1.75 * orig_size);
}

int adios_transform_alacrity_apply(struct adios_file_struct *fd,
                                   struct adios_var_struct *var,
                                   uint64_t *transformed_len,
                                   int use_shared_buffer,
                                   int *wrote_to_shared_buffer)
{
    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input_buff = var->data;

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
        log_error("Can index only real datatypes. \n");
        return 0;
    }

    // Longest parameter-parsing code ever.
    // parse the parameter relating to sigbits, with index compression
    // Read all ALACRITY parameters here
    if(var->transform_type_param) {
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

    if (config.indexForm == ALInvertedIndex) {
        ALEncode (&config, input_buff, numElements, &output_partition);
    } else if (config.indexForm == ALCompressedInvertedIndex) {
        config.indexForm = ALInvertedIndex;
        ALEncode (&config, input_buff, numElements, &output_partition);
        ALConvertIndexForm (&output_partition.metadata, &output_partition.index, ALCompressedInvertedIndex);
    } else {
        log_error("Error on indexing %s\n", var->name);
        return 0;
    }

    output_size = ALGetPartitionDataSize (&output_partition);

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

    memstream_t ms = memstreamInitReturn (output_buff);
    ALSerializePartitionData (&output_partition, &ms);

    // Check this for later. What do you intend to add in the metadata
    if(var->transform_metadata && var->transform_metadata_len > 0) {
        ((uint64_t * ) (var->transform_metadata)) [0] = ALGetMetadataSize (& (output_partition.metadata));
        ((uint64_t * ) (var->transform_metadata)) [1] = ALGetIndexSize (& (output_partition.index), & (output_partition.metadata));
        ((uint64_t * ) (var->transform_metadata)) [2] = ALGetDataSize  (& (output_partition.data), & (output_partition.metadata));
    }

    assert(output_size == ((char*)ms.ptr - (char*)ms.buf)); // Make sure we computed the output size right

    ALPartitionDataDestroy (&output_partition);
    memstreamDestroy(&ms, false);

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

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(alacrity)

#endif

