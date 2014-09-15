#include <stdint.h>
#include <assert.h>
#include <limits.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"
#include "adios_transforms_util.h"

#ifdef ALACRITY

#include "alacrity.h"
#include "alacrity-serialization-debug.h"

uint64_t adios_get_type_size(enum ADIOS_DATATYPES type, void *var);

uint16_t adios_transform_alacrity_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    return (3 * sizeof(uint64_t));
}

static ALEncoderConfig parse_configuration(const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec) {
    ALEncoderConfig config;

    if (var->pre_transform_type == adios_real) {
        ALEncoderConfigure(&config, 16, DATATYPE_FLOAT32, ALInvertedIndex);
    } else if (var->pre_transform_type == adios_double) {
        ALEncoderConfigure(&config, 16, DATATYPE_FLOAT64, ALInvertedIndex);
    } else {
        log_error("Can index only real datatypes. \n");
    }

    int i;
    for (i = 0; i < transform_spec->param_count; i++) {
        const struct adios_transform_spec_kv_pair * const param = &transform_spec->params[i];
        if (strcmp(param->key, "indexForm") == 0) {
            if (strcmp(param->value, "ALCompressedInvertedIndex") == 0)
                config.indexForm = ALCompressedInvertedIndex;
            else if (strcmp (param->value, "ALInvertedIndex") == 0)
                config.indexForm = ALInvertedIndex;
        } else if (strcmp(param->key, "sigBits") == 0) {
            config.significantBits = atoi(param->value);
        }
    }

    return config;
}

/*
uint64_t ALGetBinLayoutSize(const ALBinLayout *binLayout, uint8_t significantBits) {
    const uint8_t sigbytes = (significantBits + 0x07) >> 3;
    return    sizeof(binLayout->numBins) +                      // Number of bins
            (binLayout->numBins) * sizeof(bin_offset_t) +                   // Bin Values
            (binLayout->numBins + 1) * sizeof(bin_offset_t);    // The Bin Offsets
}
 */
static void add_bin_layout_size_growth(
		ALEncoderConfig *config,
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
	const int sigbytes = (config->significantBits + 0x07) >> 3;

	*constant_factor += sizeof(bin_id_t) + // Number of bins
						sizeof(bin_offset_t); // Ending bin offset

	// Computing an upper bound on the part of the metadata that depends on number of bins is hard.
	// Since this is data dependent, we must assume worst-case data, which results in one bin per value
	// However, there is always a maximum number of bins possible, based on the sigbits used, so after
	// we reach that many bins, we don't need to keep adding this factor. Thus, we use a "capped linear"
	// factor that scales linearly (one bin per value) until we reach the "cap" (the max number of bins,
	// already set in adios_transform_alacrity_transformed_size_growth).
	const int metadata_bytes_per_bin = 2 * sizeof(bin_offset_t);
	const int datatype_size = adios_get_type_size(var->pre_transform_type, NULL);
	const double metadata_bytes_per_data_bytes = metadata_bytes_per_bin / datatype_size;

	*capped_linear_factor += metadata_bytes_per_data_bytes; // For every value under the cap, add the corresponding number of metadata bytes
}

/*
uint64_t ALGetIndexMetadataSize(const ALIndexMetadata *indexMeta, const ALMetadata *metadata) {
    uint64_t size = sizeof(indexMeta->indexForm);
    switch (indexMeta->indexForm) {
    case ALCompressionIndex:
    case ALInvertedIndex:
        break;
    case ALCompressedInvertedIndex:
        size += (metadata->binLayout.numBins + 1) * sizeof(indexMeta->u.ciim.indexBinStartOffsets[0]);
        break;
    case ALCompressedHybridInvertedIndex:
    case ALCompressedSkipInvertedIndex:
    case ALCompressedMixInvertedIndex:
    case ALCompressedExpansionII:
    	size += (metadata->binLayout.numBins + 1) * sizeof(indexMeta->u.ciim.indexBinStartOffsets[0]);
    	break;
    default:
        abort();
    }

    return size;
}
*/
static void add_index_size_growth(
		ALEncoderConfig *config,
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
	const int datatype_size = adios_get_type_size(var->pre_transform_type, NULL);
	int metadata_bytes_per_bin = 0;

	*constant_factor += sizeof(ALIndexForm); // Index form

	switch (config->indexForm) {
    case ALCompressionIndex:
    	fprintf(stderr, "ALCompressionIndex is UNSUPPORTED in the ALACRITY transform plugin (it should not be possible to reach this error; something has gone very, very wrong)\n");
    	abort();
    	break;
    case ALInvertedIndex:
    	// No extra metadata for this index forms
    	*linear_factor += (double)sizeof(rid_t) / datatype_size; // RIDs
    	break;
    case ALCompressedInvertedIndex:
    	break;
    	*constant_factor += sizeof(uint64_t) + // Ending index bin offset
    	                    sizeof(uint64_t); // PFOR-Delta adds this much overhead to the first chunk of compressed RIDs. If there
    	                                      // are >1 chunks in a bin, assume the compression makes up for the additional overhead
    	metadata_bytes_per_bin += sizeof(uint64_t); // Index bin offsets
    	*linear_factor += (double)sizeof(rid_t) / datatype_size; // RIDs (assume 1x compression ratio as worst case)
    	break;
    case ALCompressedHybridInvertedIndex:
    case ALCompressedSkipInvertedIndex:
    case ALCompressedMixInvertedIndex:
    case ALCompressedExpansionII:
    	// All same as above, but there is a higher upper bound linear factor
    	*constant_factor += sizeof(uint64_t) +
    	                    sizeof(uint64_t);
    	metadata_bytes_per_bin += sizeof(uint64_t);
    	*linear_factor += (double)(2 * sizeof(rid_t)) / datatype_size; // RIDs (2x due to the hybrid method)
    	break;
    default:
        break;
	}

	if (metadata_bytes_per_bin > 0) {
		const double metadata_bytes_per_data_bytes = metadata_bytes_per_bin / datatype_size;
		*capped_linear_factor += metadata_bytes_per_data_bytes; // For every value under the cap, add the corresponding number of metadata bytes
	}
}

static void add_data_size_growth(
		ALEncoderConfig *config,
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
	const int datatype_size = adios_get_type_size(var->pre_transform_type, NULL);
	const int insigbytes = ((datatype_size << 3) - config->significantBits + 0x07) >> 3;

	*linear_factor += (double)insigbytes / datatype_size; // insigbytes per data value
}

/*
uint64_t ALGetMetadataSize(const ALMetadata *metadata) {
    return  sizeof(global_rid_t) +                                                  // Global RID offset
            sizeof(partition_length_t) +                                            // Partition length
            ALGetBinLayoutSize(&metadata->binLayout, metadata->significantBits) +   // Bin layout metadata
            ALGetIndexMetadataSize(&metadata->indexMeta, metadata) +                // Index metadata
            sizeof(char) +                                                          // Significant bytes
            sizeof(char) +                                                          // Element size
            sizeof(ALDatatype) +                                                    // Datatype
            sizeof(char);                                                           // Endianness
}
*/
void adios_transform_alacrity_transformed_size_growth(
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
	ALEncoderConfig config = parse_configuration(var, transform_spec);

	// ALACRITY metadata (not counting bin layout metadata and index-specific metadata, added in below)
	*constant_factor =
			sizeof(global_rid_t) +													// Global RID offset
			sizeof(partition_length_t) +                                            // Partition length
			sizeof(char) +                                                          // Significant bytes
			sizeof(char) +                                                          // Element size
			sizeof(ALDatatype) +                                                    // Datatype
			sizeof(char);															// Endianness

	*linear_factor = 0; // We totally re-encode the data, so start from scratch (i.e., we aren't adding to existing data, we're replacing it)

	// Some of the following metadata upper bounds are proportional to the number of bins,
	// which is worst-case linear with number of values, but only up to the max possible bins.
	// These bounds will all use the same capped linear cap (bytes corresponding to number of
	// values equal to max possible bins), so compute it here, and let each function add to
	// the capped linear factor assuming this cap is in place.
	const uint64_t max_possible_bins = (1ULL << config.significantBits);
	const int datatype_size = adios_get_type_size(var->pre_transform_type, NULL);
	*capped_linear_cap = max_possible_bins * datatype_size; // There can be at most one bin per value, so stop adding this factor after (max bins) * (bytes per value) bytes
	*capped_linear_factor = 0;

#ifdef ALACRITY_DEBUG
#define PRINT_FACTORS(msg) fprintf(stderr, "%s: const=%llu, lin=%lf, lincap=%lf->%llu\n", (msg), *constant_factor, *linear_factor, *capped_linear_factor, *capped_linear_cap);
	fprintf(stderr, "ALACRITY growth compute info: datatype size: %d, maxbins: %llu\n", datatype_size, max_possible_bins);
	PRINT_FACTORS("after metadata");
#else
#define PRINT_FACTORS(msg) ((void)0)
#endif

	add_bin_layout_size_growth(&config, var, transform_spec, constant_factor, linear_factor, capped_linear_factor, capped_linear_cap);
	PRINT_FACTORS("after bin layout");
	add_index_size_growth(&config, var, transform_spec, constant_factor, linear_factor, capped_linear_factor, capped_linear_cap);
	PRINT_FACTORS("after index");
	add_data_size_growth(&config, var, transform_spec, constant_factor, linear_factor, capped_linear_factor, capped_linear_cap);
	PRINT_FACTORS("after data");

	// Phew!
}

int adios_transform_alacrity_apply(struct adios_file_struct *fd,
                                   struct adios_var_struct *var,
                                   uint64_t *transformed_len,
                                   int use_shared_buffer,
                                   int *wrote_to_shared_buffer)
{
    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var);
    const void *input_buff = var->data;

    // Determine the ALACRITY encoder configuration to use
    ALEncoderConfig config = parse_configuration(var, var->transform_spec);

    uint32_t numElements = 0;
    numElements = input_size / adios_get_type_size(var->pre_transform_type, NULL);

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

    output_size = ALGetPartitionDataSize(&output_partition);

#ifdef ALACRITY_DEBUG
    fprintf(stderr, "ALGetPartitionDataSize: %llu\n", output_size);
#endif

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

