/*
 * adios_transform_zfp_write.c
 *
 *  Created on: June 30, 2016
 *      Author: Eric Suchyta
 */

#ifdef BZIP2

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>

#include "core/adios_logger.h"
#include "core/transforms/adios_transforms_common.h"
#include "core/transforms/adios_transforms_write.h"
#include "core/transforms/adios_transforms_hooks_write.h"
#include "core/transforms/adios_transforms_util.h"
#include "core/adios_internals.h"

#include "zfp.h"



uint16_t adios_transform_zfp_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    return (sizeof(uint64_t) + sizeof(char));    // metadata: original data size (uint64_t) + compression succ flag (char)
}

void adios_transform_zfp_transformed_size_growth(
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
	// Do nothing (defaults to "no transform effect on data size")
}

int adios_transform_zfp_apply(struct adios_file_struct *fd,
                                struct adios_var_struct *var,
                                uint64_t *transformed_len,
                                int use_shared_buffer,
                                int *wrote_to_shared_buffer)
{

	void* outbuffer;	// What to send to ADIOS
	zfp_type type;		// Map adios_type into zfp_type
	uint choice; 


	/* adios to zfp datatype */
	if (var->pre_transform_type == adios_double) type = zfp_type_double;
	else if (var->pre_transform_type == adios_float) type = zfp_type_float;
	else zpf_error("A datatype ZFP does not understand was given for compression. Understood types are adios_double, adios_float.");

	/* dimensionality */
	uint ndims = (uint) count_dimensions(var->pre_transform_dimensions);
	int dims = get_dimension(var->pre_transform_dimensions, ndims);
	for (i=0; i<ndims; i++) dims[i] = var->pre_transform_dimensions[i]->rank;

	/* Which zfp mode to use */
	for (i=0; i<var->transform_spec->param_count; i++) 
	{
		const struct adios_transform_spec_kv_pair* const param = &var->transform_spec->params[i];
		if (strcmp(param->key, "accuracy") == 0) choice = 0;
		else if (strcmp(param->key, "precision") == 0) choice = 1;
		else if (strcmp(param->key, "rate") == 0) choice = 2;
		else zfp_error("An unknown compression mode was specified for zfp: %s. Availble choices are: accuracy, precision, rate.", param->key)
	}

	/* do compression */
	int status = _zfp_comp(param->value, var->data, ndims, dims, type, choice, 0, field, zfp, stream, outbuffer);
        var->adata = outbuffer;



    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var);
    const void *input_buff = var->data;

    // parse the compressiong parameter
    /* Pre-specparse code
    if(var->transform_type_param
        && strlen(var->transform_type_param) > 0
        && is_digit_str(var->transform_type_param))
    {
        compress_level = atoi(var->transform_type_param);
        if(compress_level > 9 || compress_level < 1)
        {
            compress_level = 9;
        }
    }
    */
    int compress_level = 9;
    if (var->transform_spec->param_count > 0) {
        compress_level = atoi(var->transform_spec->params[0].key);
        if (compress_level < 1 || compress_level > 9)
            compress_level = 9;
    }


    // decide the output buffer
    uint64_t output_size = input_size; //adios_transform_bzip2_calc_vars_transformed_size(adios_transform_bzip2, input_size, 1);
    void* output_buff = NULL;

    if (use_shared_buffer)    // If shared buffer is permitted, serialize to there
    {
        *wrote_to_shared_buffer = 1;
        if (!shared_buffer_reserve(fd, output_size))
        {
            log_error("Out of memory allocating %" PRIu64 " bytes for %s for bzip2 transform\n", output_size, var->name);
            return 0;
        }

        // Write directly to the shared buffer
        output_buff = fd->buffer + fd->offset;
    }
    else    // Else, fall back to var->adata memory allocation
    {
        *wrote_to_shared_buffer = 0;
        output_buff = malloc(output_size);
        if (!output_buff)
        {
            log_error("Out of memory allocating %" PRIu64 " bytes for %s for bzip2 transform\n", output_size, var->name);
            return 0;
        }
    }

    uint64_t actual_output_size = output_size;
    char compress_ok = 1;

    int rtn = compress_bzip2_pre_allocated(input_buff, input_size, output_buff, &actual_output_size, compress_level);

    if(0 != rtn                     // compression failed for some reason, then just copy the buffer
        || actual_output_size > input_size)  // or size after compression is even larger (not likely to happen since compression lib will return non-zero in this case)
    {
        // printf("compression failed, fall back to memory copy\n");
        memcpy(output_buff, input_buff, input_size);
        actual_output_size = input_size;
        compress_ok = 0;    // succ sign set to 0
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
        memcpy((char*)var->transform_metadata, &input_size, sizeof(uint64_t));
        memcpy((char*)var->transform_metadata + sizeof(uint64_t), &compress_ok, sizeof(char));
    }

    *transformed_len = actual_output_size; // Return the size of the data buffer

    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(zfp)

#endif

