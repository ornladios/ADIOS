#ifdef COMPRESS

#include <stdint.h>
#include <assert.h>
#include <limits.h>

#include "adios_logger.h"
#include "adios_transforms_hooks.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_write.h"
#include "compress.h"

int adios_transform_compress_apply(struct adios_file_struct *fd, struct adios_var_struct *var, uint64_t *transformed_len, int *use_shared_buffer) {
    // Assume this function is only called for COMPRESS transform type
    assert(var->transform_type == adios_transform_compress 
			|| var->transform_type == adios_transform_compress_zlib 
			|| var->transform_type == adios_transform_compress_bzlib2 
			|| var->transform_type == adios_transform_compress_szip);

    // Get the input data and data length
    const uint64_t input_size_64 = adios_transform_get_pre_transform_var_size(fd->group, var);
    const void *input_buff= var->data;
	
	if(input_size_64 > UINT_MAX) // too large, do not support 64 bit now
	{
		log_error("data size too long %llu bytes for %s for compress transform\n", input_size_64, var->name);       
		return 0;
	}
	
	const uint32_t input_size = (uint32_t)input_size_64;
	
	// compress it
		
	uint32_t output_size = EXPAND_SIZE(input_size);
	void* output_buff = NULL;
	
	uint64_t mem_allowed = 0;
	 
	if (*use_shared_buffer) 
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
	
	int rtn = 0;
	// uint32_t output_size_save = output_size;
	
	switch(var->transform_type)
	{
		case adios_transform_compress: 
		case adios_transform_compress_zlib: 
			rtn = compress_zlib_pre_allocated(input_buff, input_size, output_buff, &output_size);
		break;
		
		case adios_transform_compress_bzlib2:
			rtn = compress_bzlib2_pre_allocated(input_buff, input_size, output_buff, &output_size);
		break;
		
		case adios_transform_compress_szip:
			rtn = compress_zlib_pre_allocated(input_buff, input_size, output_buff, &output_size);
		break;
		
		default:
			rtn = compress_zlib_pre_allocated(input_buff, input_size, output_buff, &output_size);
		break;
	
	}
	
	if(0 != rtn)
	{
		if (!(*use_shared_buffer))
		{			
			adios_method_buffer_free(mem_allowed);
			if(output_buff)
			{
				free(output_buff);
				output_buff = NULL;
			}
		}
        log_error("Compress failed for %s for compress transform\n", var->name);
        return 0;
	}
    
	// Wrap up, depending on buffer mode
    if (*use_shared_buffer) 
	{
        shared_buffer_mark_written(fd, output_size);
    } 
	else 
	{
        var->data = output_buff;
        var->data_size = output_size;
        var->free_data = adios_flag_yes;
    }


    *transformed_len = output_size; // Return the size of the data buffer
    return 1;
}

#endif
