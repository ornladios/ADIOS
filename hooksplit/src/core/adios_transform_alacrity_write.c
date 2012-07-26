#ifdef ALACRITY

#include <stdint.h>
#include <assert.h>
#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_write.h"
#include "adios_transforms_util.h"
#include "alacrity.h"

uint16_t adios_transform_alacrity_get_metadata_size() { return 2 * sizeof(uint64_t); }

#define MAX_POSSIBLE_BINS 65536
static const int MAX_PART_METADATA_SIZE =
    sizeof(uint64_t) +
    sizeof(unsigned short int) +
    MAX_POSSIBLE_BINS * ( sizeof(unsigned short int) +  2 * sizeof(uint64_t) + sizeof(unsigned char) );


uint64_t adios_transform_alacrity_calc_vars_transformed_size(uint64_t orig_size, int num_vars) {
    return num_vars * MAX_PART_METADATA_SIZE +	// For the metadata
           orig_size * 5/4;						// For the index + data
}

int adios_transform_alacrity_apply(struct adios_file_struct *fd, struct adios_var_struct *var, uint64_t *transformed_len,
                                   int use_shared_buffer, int *wrote_to_shared_buffer) {

    memstream_t payload_ms, transform_meta_ms;
    uint64_t data_offset, index_offset;
    ALEncoderConfig econfig;
    ALPartitionData part;

    // Assume this function is only called for ALACRITY transform type
    assert(var->transform_type == adios_transform_alacrity);

    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input = var->data;

    // Transform the data using the ALACRITY encoder
    econfig.elementSize = adios_get_type_size(var->pre_transform_type, NULL);
    econfig.significantBytes = 2;
    econfig.indexForm = ALInvertedIndex;

    if (econfig.elementSize != 2 && econfig.elementSize != 4 && econfig.elementSize != 8) {
        log_error("Error: Cannot apply indexing method \"alacrity\" to variable %s: datatype %d is unsupported\n",
                  var->name, var->pre_transform_type);
        return 0;
    }

    ALEncode(&econfig, input, input_size / adios_get_type_size(var->pre_transform_type, NULL), &part);

    // Allocate an output buffer, as a chunk of the shared buffer if possible
    const uint64_t output_size = ALGetPartitionDataSize(&part);;
    void *output;
    if (use_shared_buffer) {	// If shared buffer is permitted, serialize to there
        *wrote_to_shared_buffer = 1;

        if (!shared_buffer_reserve(fd, output_size)) {
            log_error("Out of memory allocating %llu bytes for %s for ALACRITY transform\n", output_size, var->name);
            return 0;
        }

        // Write directly to the shared buffer
        output = fd->buffer + fd->offset;
    } else {					// Else, fall back to var->data memory allocation
        *wrote_to_shared_buffer = 0;

        /*
        uint64_t mem_allowed = adios_method_buffer_alloc(output_size);
        // If we aren't allowed enough memory, or the malloc fails...
        if (mem_allowed != output_size ||
            !(output = malloc(output_size)) )
        {
            adios_method_buffer_free(mem_allowed);
            log_error("Out of memory allocating %llu bytes for %s for ALACRITY transform\n", output_size, var->name);
            return 0;
        }*/
        output = malloc(output_size);
        if (!output)
        {
            log_error("Out of memory allocating %llu bytes for %s for ALACRITY transform\n", output_size, var->name);
            return 0;
        }
    }

    // Serialize ALACRITY data
    payload_ms = memstreamInitReturn(output);
    transform_meta_ms = memstreamInitReturn(var->transform_metadata);

    ALSerializeMetadata(&part.metadata, &payload_ms);
    data_offset = memstreamGetPosition(&payload_ms);
    ALSerializeData(&part.data, &part.metadata, &payload_ms);
    index_offset = memstreamGetPosition(&payload_ms);
    ALSerializeIndex(&part.index, &part.metadata, &payload_ms);
    ALPartitionDataDestroy(&part);

    // Write out the transform metadata
    memstreamAppendUint64(&transform_meta_ms, data_offset);
    memstreamAppendUint64(&transform_meta_ms, index_offset);

    // Make sure we wrote at most what we said we would
    const uint64_t actual_output_size = memstreamGetPosition(&payload_ms);
    assert(actual_output_size <= output_size);

    // Wrap up, depending on buffer mode
    if (use_shared_buffer) {
        shared_buffer_mark_written(fd, actual_output_size);
    } else {
        var->data = payload_ms.buf;
        var->data_size = actual_output_size;
        var->free_data = adios_flag_yes;
    }

    // Clean up memstream cursors
    memstreamDestroy(&transform_meta_ms, false);
    memstreamDestroy(&payload_ms, false);

    // Return
    *transformed_len = output_size; // Return the size of the data buffer
    return 1;
}

#endif
