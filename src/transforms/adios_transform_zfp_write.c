/*
 * adios_transform_zfp_write.c
 *
 * 	Author: Eric Suchyta
 * 	Contact: eric.d.suchyta@gmail.com
 */

#ifdef ZFP


/* general C stuff */
#include <stdint.h>	// uint64_t
#include <stdio.h> 	// NULL, sprintf
#include <stdlib.h>	// NULL, malloc, free
#include <string.h>	// memcpy, strcmp, strlen


/* Were in the template included from ADIOS. Not necessarily sure if they're all strictly needed. */
#include "core/transforms/adios_transforms_common.h"
#include "core/transforms/adios_transforms_write.h"
#include "core/transforms/adios_transforms_hooks_write.h"
#include "core/transforms/adios_transforms_util.h"


/* Extra ADIOS headers that weren't added in the template */
//#include "core/adios_internals.h" 	// count_dimensions


/* ZFP specific */
#include "adios_transform_zfp_common.h"
#include "zfp.h"


/* see zfp_metadata in adios_transform_zfp_common.h */
uint16_t adios_transform_zfp_get_metadata_size(struct adios_transform_spec *transform_spec)
{
	return (2*sizeof(uint64_t) + sizeof(uint) + 2*ZFP_STRSIZE);
}


/* Template says: 'Doing nothing defaults to "no transform effect on data size"'. 
 * I think this means a blank function equates to I don't need to grow the array to do transform */
void adios_transform_zfp_transformed_size_growth(const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
	return;
}


static void get_dims(const struct adios_dimension_struct* d, struct zfp_buffer* zbuff, struct adios_var_struct* var)
{
	int i;	
	zbuff->dims = (uint*) malloc(zbuff->ndims*sizeof(uint));
	while (d)
	{
		zbuff->dims[i] = (uint) adios_get_dimension_space_size(var, (struct adios_dimension_struct*) d);
		d = d->next;
	}
	return;
}


/* Does the main compression work */
int adios_transform_zfp_apply(struct adios_file_struct *fd, struct adios_var_struct *var, 
		uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer)
{

	int success; 			// Did (some part of) compression succeed?
	void* outbuffer = NULL;		// What to send to ADIOS
	uint64_t outsize;		// size of output buffer

	struct zfp_buffer* zbuff = (struct zfp_buffer*) malloc(sizeof(struct zfp_buffer));	// Handle zfp streaming
	uint64_t insize = adios_transform_get_pre_transform_var_size(var); 			// size of input buffer

	/* adios to zfp datatype */
	strcpy(zbuff->name, var->name);
	memset(zbuff->msg, '\0', ZFP_STRSIZE);
	success = zfp_get_datatype(zbuff, var->pre_transform_type);
	if (!success)
	{
		return 0;
	}


	/* dimensionality */
	zbuff->ndims = (uint) count_dimensions(var->pre_transform_dimensions);
	//get_dimensions(var->pre_transform_dimensions, zbuff, var);

	//struct adios_dimension_struct* d = adios_get_dimension_space_size(var, var->pre_transform_dimensions);
	struct adios_dimension_struct* d = var->pre_transform_dimensions;
	get_dims(d, zbuff, var);


	/* make sure the user only gives the sensible number of key:values -- 1. */
	if (var->transform_spec->param_count == 0)
	{
		sprintf(zbuff->msg, "No compression mode specified. Choose from: accuracy, precision, rate");
		zfp_error(zbuff);
		return 0;
	}
	else if (var->transform_spec->param_count > 1)
	{
		sprintf(zbuff->msg, "Too many parameters specified. You can only give one key:value, the compression mode and it's tolerance.");
		zfp_error(zbuff);
		return 0;
	}
	else if (var->transform_spec->param_count < 0)
	{
		sprintf(zbuff->msg, "Negative number of parameters interpretted. This shouldn't happen.");
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
		sprintf(zbuff->msg, "An unknown compression mode was specified for zfp: %s. Availble choices are: accuracy, precision, rate.", param->key);
		zfp_error(zbuff);
		return 0;
	}
	strcpy(zbuff->ctol, param->value);


	/* do compression */
	if (use_shared_buffer)
	{
		*wrote_to_shared_buffer = 1;
	}
	else
	{
		*wrote_to_shared_buffer = 0;
	}
	//*wrote_to_shared_buffer = use_shared_buffer;
	success = zfp_compression(zbuff, var->data, &outbuffer, &outsize, use_shared_buffer, fd);


	/* What do do if compresssion fails. For now, just give up. Maybe eventually use raw data. */
	if(!success)
	{
		return 0;
	}

	
	/* Write the data */
	if (*wrote_to_shared_buffer) 
	{
		shared_buffer_mark_written(fd, outsize);
		//shared_buffer_mark_written(fd, zbuff->buffsize);
	} 
	else 
	{
		var->adata = outbuffer;
		var->data_size = outsize;
		var->free_data = adios_flag_yes;
	}


	/* Write the transform metadata */
	char* pos = (char*)var->transform_metadata;
	size_t offset = 0;
	
	if(var->transform_metadata && var->transform_metadata_len > 0)
	{
		zfp_write_metadata_var(pos, &insize, sizeof(uint64_t), &offset);
		zfp_write_metadata_var(pos, &outsize, sizeof(uint64_t), &offset);
		zfp_write_metadata_var(pos, &zbuff->mode, sizeof(uint), &offset);
		zfp_write_metadata_var(pos, zbuff->ctol, ZFP_STRSIZE, &offset);
		zfp_write_metadata_var(pos, zbuff->name, ZFP_STRSIZE, &offset);
	}

	printf("buffsize: %u\n", zbuff->buffsize);
	printf("outsize: %u\n", outsize);

	*transformed_len = outsize; // Return the size of the data buffer
	//*transformed_len = zbuff->buffsize; // Return the size of the data buffer

	/* clean up */
	free(zbuff);

	return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(zfp)

#endif

