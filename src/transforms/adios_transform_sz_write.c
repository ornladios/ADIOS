/*
 * adios_transform_sz_write.c
 *
 * 	Author: Jong Choi
 * 	Contact: choij@ornl.gov
 */

#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <sys/time.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_hooks_write.h"
#include "adios_transforms_util.h"

#include "sz.h"

#ifdef HAVE_SZ

typedef struct
{
    int r[5];
} sz_info_t;

uint16_t adios_transform_sz_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    log_info("function: %s\n", __FUNCTION__);
    return 0; // Set amount of transform-internal metadata space to allocate
}

void adios_transform_sz_transformed_size_growth(
		const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
		uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
    log_info("function: %s\n", __FUNCTION__);
}

int adios_transform_sz_apply(struct adios_file_struct *fd,
                                   struct adios_var_struct *var,
                                   uint64_t *transformed_len,
                                   int use_shared_buffer,
                                   int *wrote_to_shared_buffer)
{
    log_info("function: %s\n", __FUNCTION__);
    log_info("use_shared_buffer: %d\n", use_shared_buffer);
    
    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var);
    const void *input_buff = var->data;

    sz_params sz;
    sz.dataEndianType = LITTLE_ENDIAN_DATA;
    //sz.sysEndianType = BIG_ENDIAN_DATA;
    sz.sol_ID = SZ;
    sz.layers = 1;
    sz.sampleDistance = 50;
    sz.quantization_intervals = 0;
    sz.predThreshold = 0.98;
    sz.offset = 0;
    sz.szMode = SZ_DEFAULT_COMPRESSION;
    sz.gzipMode = 1;
    sz.errorBoundMode = REL;
    sz.absErrBound = 1E-6;
    sz.relBoundRatio = 1E-5;

    SZ_Init_Params(&sz);

    unsigned char *bytes;
    int outsize;
    int r[5] = {0,0,0,0,0};

    // Get type info
    int dtype;
    switch (var->pre_transform_type)
    {
        case adios_double:
            dtype = SZ_DOUBLE;
            break;
        case adios_real:
            dtype = SZ_FLOAT;
            break;
        default:
            adios_error(err_transform_failure, "No supported data type\n");
            return NULL;
            break;
    }
    
    // Get dimension info
    struct adios_dimension_struct* d = var->pre_transform_dimensions;
    int ndims = (uint) count_dimensions(d);
    printf("ndims: %d\n", ndims);
    if (ndims > 5)
    {
        adios_error(err_transform_failure, "No more than 5 dimension is supported.\n");
        return NULL;
    }

    int i, ii;
    for (i=0; i<ndims; i++)
    {
        uint dsize = (uint) adios_get_dim_value(&d->dimension);
        if (fd->group->adios_host_language_fortran == adios_flag_yes)
            ii = ndims - 1 - i;
        else
            ii = i;
        r[5-ndims+ii] = dsize;
        d = d->next;
    }
    printf("r: %d %d %d %d %d\n", r[0], r[1], r[2], r[3], r[4]);
    
    bytes = SZ_compress(dtype, (void *) input_buff, &outsize, r[0], r[1], r[2], r[3], r[4]);
    printf("SZ outsize: %d\n", outsize);
        
    // Output
    uint64_t output_size = outsize/* Compute how much output size we need */;
    void* output_buff;

    if (use_shared_buffer) {
        // If shared buffer is permitted, serialize to there
        assert(shared_buffer_reserve(fd, output_size));

        // Write directly to the shared buffer
        output_buff = fd->buffer + fd->offset;
        memcpy(output_buff, bytes, outsize);
    } else { // Else, fall back to var->adata memory allocation
        output_buff = bytes;
        //assert(output_buff);
    }
    *wrote_to_shared_buffer = use_shared_buffer;

    // Do transform from input_buff into output_buff, and update output_size to the true output size

    // Wrap up, depending on buffer mode
    if (*wrote_to_shared_buffer) {
        shared_buffer_mark_written(fd, output_size);
    } else {
        var->adata = output_buff;
        var->data_size = output_size;
        var->free_data = adios_flag_yes;
    }

    *transformed_len = output_size; // Return the size of the data buffer
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(sz)

#endif

