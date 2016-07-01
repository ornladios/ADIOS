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


int adios_transform_zfp_apply(struct adios_file_struct *fd, struct adios_var_struct *var, uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer)
{

	void* outbuffer;	// What to send to ADIOS
	int outsize;		// size of output buffer
	zfp_type type;		// Map adios_type into zfp_type
	uint choice; 		// zfp mode
	int success; 		// Did compression succeed?


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
	int status = _zfp_compress(param->value, var->data, ndims, dims, type, choice, 0, field, zfp, stream, outbuffer);
	success = zfp_compression(param->value, var->array, ndims, dims, type, choice, var->name, outbuffer, outsize, fd, use_shared_buffer, fd);


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

