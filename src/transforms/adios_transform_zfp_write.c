/*
 * adios_transform_zfp_write.c
 *
 *  Created on: June 30, 2016
 *      Author: Eric Suchyta
 */

#ifdef ZFP

#include <stdint.h>	// uint64_t
#include <stdio.h> 	// NULL, sprintf
#include <stdlib.h>	// NULL, malloc, free
#include <string.h>	// memcpy, strcmp, strlen

#include "core/transforms/adios_transforms_common.h"
#include "core/transforms/adios_transforms_write.h"
#include "core/transforms/adios_transforms_hooks_write.h"
#include "core/transforms/adios_transforms_util.h"
#include "core/adios_internals.h" 	// count_dimensions

#include "zfp.h"
#include "adios_tranform_zfp_common.h"



uint16_t adios_transform_zfp_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    return (2*sizeof(uint64_t));    // metadata: original data size (uint64_t) + compressed size (uint64_t)
}

void adios_transform_zfp_transformed_size_growth(
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
	// Do nothing (defaults to "no transform effect on data size")
}


int adios_transform_zfp_apply(struct adios_file_struct *fd, struct adios_var_struct *var, uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer)
{

	void* outbuffer;							// What to send to ADIOS
	uint64_t outsize;							// size of output buffer
	uint64_t insize = adios_transform_get_pre_transform_var_size(var); 	// size of input buffer


	int success; 			// Did compression succeed?
	struct zfp_buffer* zbuff;	// Handle zfp streaming
	zbuff->name = var->name


	/* adios to zfp datatype */
	if (var->pre_transform_type == adios_double) 
	{
		zbuff->type = zfp_type_double;
	}
	else if (var->pre_transform_type == adios_float) 
	{
		zbuff->type = zfp_type_float;
	}
	else 
	{
		sprintf(zbuff->msg, "A datatype ZFP does not understand was given for compression. Understood types are adios_double, adios_float.")
		zfp_error(zbuff);
		return 0;
	}

	/* dimensionality */
	zbuff->ndims = (uint) count_dimensions(var->pre_transform_dimensions);
	get_dimension(var->pre_transform_dimensions, zbuff);


	/* make sure the user only gives the sensible number of key:values -- 1. */
	if (var->transform_spec->param_count == 0)
	{
		sprintf(zbuff->errmsg, "No compression mode specified. Choose from: accuracy, precision, rate");
		zfp_error(zbuff);
		return 0;
	}
	else if (var->transform_spec->param_count > 1)
	{
		sprintf(zbuff->errmsg, "Too many parameters specified. You can only give one key:value, the compression mode and it's tolerance.");
		zfp_error(zbuff);
		return 0;
	}
	else if (var->transform_spec->param_count < 0)
	{
		sprintf(zbuff->errmsg, "Negative number of parameters interpretted. This shouldn't happen.");
		zfp_error(zbuff);
		return 0;
	}


	/* Which zfp mode to use */
	const struct adios_transform_spec_kv_pair* const param = &var->transform_spec->params[0];
	if (strcmp(param->key, "accuracy") == 0) 
	{
		zbuff->mode = 0;
	}
	else if (strcmp(param->key, "precision") == 0)
	{
		zbuff->mode = 1;
	}
	else if (strcmp(param->key, "rate") == 0)
	{
		zbuff->mode = 2;
	}
	else 
	{
		sprintf(zbuff->errmsg, "An unknown compression mode was specified for zfp: %s. Availble choices are: accuracy, precision, rate.", param->key)
		zfp_error(zbuff);
		return 0;
	}
	zbuff->ctol = param->value;


	/* do compression */
	success = zfp_compression(zbuff, var->array, outbuffer, outsize, use_shared_buffer, fd);


	/* What do do if compresssion fails. For now, just give up */
	if(!success)
	{
		return 0;
		// memcpy(output_buff, input_buff, input_size);
		// actual_output_size = input_size;
		// compress_ok = 0; 
	}


	*wrote_to_shared_buffer = use_shared_buffer;
	if (*wrote_to_shared_buffer) 
	{
		shared_buffer_mark_written(fd, outsize);
	} 
	else 
	{
		var->adata = outbuffer;
		var->data_size = outsize;
		var->free_data = adios_flag_yes;
	}


	/* Copy the transform metadata: 
	 * original size before compression
	 * compressed size 
	 */
	if(var->transform_metadata && var->transform_metadata_len > 0)
	{
		memcpy((char*)var->transform_metadata, &insize, sizeof(uint64_t));
		memcpy((char*)var->transform_metadata + sizeof(uint64_t), &outsize, sizeof(uint64_t));
	}
	*transformed_len = outsize; // Return the size of the data buffer

    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(zfp)

#endif

