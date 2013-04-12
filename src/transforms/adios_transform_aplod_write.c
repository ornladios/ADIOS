#include <stdint.h>
#include <assert.h>
#include <limits.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"

#ifdef APLOD

#include "aplod.h"

#define TEST_SIZE(s,numComponents) (s + sizeof (numComponents) + numComponents * sizeof (int32_t))

static int is_digit_str(char* input_str)
{
	return 1;
}

int compress_aplod_pre_allocated(const void* input_data, const uint64_t input_len, 
								void* output_data, uint64_t* output_len, int compress_level)
{
	assert(input_data != NULL && input_len > 0 && output_data != NULL && output_len != NULL && *output_len > 0);
	
	// uLongf temp = *output_len;	
	// int zerr = compress2((Bytef*)(output_data), &temp, 
							// (Bytef*)input_data, (uLongf)input_len, 
							// compress_level);
	
	// *output_len = (uint64_t)temp;
	
	// if(zerr != Z_OK)
	// {
		// printf("aplod compress2 error %d\n", zerr);
		// return -1;
	// }
	return 0;
}

uint16_t adios_transform_aplod_get_metadata_size() 
{ 
    // Write the component vector here
    // No more than 8 components per variable
	return (sizeof (uint64_t) + sizeof (int8_t) + 8 * sizeof(int32_t));
}

uint64_t adios_transform_aplod_calc_vars_transformed_size(uint64_t orig_size, int num_vars) 
{
    return TEST_SIZE(orig_size, sizeof (long double));
}

int adios_transform_aplod_apply(struct adios_file_struct *fd, 
								struct adios_var_struct *var,
								uint64_t *transformed_len, 
								int use_shared_buffer, 
								int *wrote_to_shared_buffer)
{
    // Assume this function is only called for COMPRESS transform type
    assert(var->transform_type == adios_transform_aplod);

    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input_buff= var->data;
 
    // max size supported is long double
    int32_t componentVector [4];
    int8_t numComponents = 0;
    int32_t componentTotals = 0; 

    // printf ("[%s] byte = %d, real = %d, double = %d, this = %d\n", var->name, adios_byte, adios_real, adios_double, var->pre_transform_type);
	// parse the aplod component vector parameters
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

        // Get the key
        char *key = strtok (transform_param, "=");

        if (strcmp (key, "CV") == 0) {

            char *value = strtok (key, ",");
            int32_t componentID = 0;

            while (value) {
                int32_t component = atoi (value);
                if (component <= 0) {
                    numComponents = 0;
                    break ;
                }

                componentVector [componentID] = component;
                componentTotals += component;
                componentID ++;
            }
        }
	}

    if ((numComponents == 0) || (componentTotals != bp_get_type_size (var->pre_transform_type, ""))) {
        if (var->pre_transform_type == adios_double) {
            componentVector [0] = 2;
            componentVector [1] = 2;
            componentVector [2] = 2;
            componentVector [3] = 2;
            numComponents = 4;
        } else if (var->pre_transform_type == adios_real) {
            componentVector [0] = 2;
            componentVector [1] = 2;
            numComponents = 2;
        }
    }

    // decide the output buffer
    uint64_t output_size = TEST_SIZE(input_size, numComponents);
    void* output_buff = NULL;

    uint64_t mem_allowed = 0;
    if (use_shared_buffer)
    {	
		*wrote_to_shared_buffer = 1;
		// If shared buffer is permitted, serialize to there
        if (!shared_buffer_reserve(fd, output_size))
        {
            log_error("Out of memory allocating %llu bytes for %s for aplod transform\n", output_size, var->name);
            return 0;
        }

        // Write directly to the shared buffer
        output_buff = fd->buffer + fd->offset;
    }
    else	// Else, fall back to var->data memory allocation
    {
		*wrote_to_shared_buffer = 0;
		output_buff = malloc(output_size);
        if (!output_buff)
        {
            log_error("Out of memory allocating %llu bytes for %s for aplod transform\n", output_size, var->name);
            return 0;
        }
    }

    // printf ("Numcomponents: %hu\n", numComponents);

    // APLOD specific code - Start
    uint32_t numElements = input_size / bp_get_type_size (var->pre_transform_type, "");

    APLODConfig_t *config = APLODConfigure (componentVector, numComponents);

    config->blockLengthElts = numElements;
    APLODShuffleComponents (config, numElements, 0, numComponents, input_buff, output_buff);

    // void *decompressed_buff = (void *) calloc (input_size, sizeof (char));
    // APLODReconstructComponents  (config, numElements, 0, numComponents, 0, 0, decompressed_buff, output_buff);
    // assert (memcmp (input_buff, decompressed_buff, numElements * sizeof (double)) == 0);
    // free (decompressed_buff);
 
    // APLOD specific code - End 

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

    // Do I copy the PLODHandle_t object as the metadata or do I serialize it into the buffer as well
    if(var->transform_metadata && var->transform_metadata_len > 0)
    {
        memcpy (var->transform_metadata, &input_size, sizeof(uint64_t));
        memcpy (var->transform_metadata + sizeof (uint64_t), &numComponents, sizeof (numComponents));
        memcpy (var->transform_metadata + sizeof (uint64_t) + sizeof (numComponents), componentVector, numComponents * sizeof (int32_t));

    }

    free (config->byteVector);
    free (config->byteVectorPS);
    free (config);

	// printf("adios_transform_compress_apply compress %d input_size %d actual_output_size %d\n", 
			// rtn, input_size, actual_output_size);
	
    *transformed_len = output_size; // Return the size of the data buffer
	
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(aplod)

#endif

