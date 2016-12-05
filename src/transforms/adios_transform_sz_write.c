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
    //log_debug("function: %s\n", __FUNCTION__);
    return 0; // Set amount of transform-internal metadata space to allocate
}

void adios_transform_sz_transformed_size_growth(const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
                                                uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
    //log_debug("function: %s\n", __FUNCTION__);
}

int adios_transform_sz_apply(struct adios_file_struct *fd,
                             struct adios_var_struct *var,
                             uint64_t *transformed_len,
                             int use_shared_buffer,
                             int *wrote_to_shared_buffer)
{
    //log_debug("function: %s\n", __FUNCTION__);
    //log_debug("use_shared_buffer: %d\n", use_shared_buffer);
    
    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var);
    const void *input_buff = var->data;
    
    sz_params sz = {};
    sz.dataEndianType = LITTLE_ENDIAN_DATA;
    sz.sysEndianType = LITTLE_ENDIAN_DATA;
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
    
    int  errorBoundMode = REL;
    double absErrBound = 1E-6;
    double relBoundRatio = 1E-5;
    
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
    //log_debug("ndims: %d\n", ndims);
    if (ndims > 5)
    {
        adios_error(err_transform_failure, "No more than 5 dimension is supported.\n");
        return NULL;
    }
    
    int i = 0, ii = 0;
    for (i=0; i<ndims; i++)
    {
        uint dsize = (uint) adios_get_dim_value(&d->dimension);
        if (fd->group->adios_host_language_fortran == adios_flag_yes)
            ii = ndims - 1 - i;
        else
            ii = i;
        r[ii] = dsize;
        d = d->next;
    }
    
    /* SZ parameters */
    struct adios_transform_spec_kv_pair* param;
    for (i=0; i<var->transform_spec->param_count; i++)
    {
        param = &(var->transform_spec->params[i]);
        if (strcmp(param->key, "errorboundmode") == 0)
        {
            errorBoundMode = atoi(param->value);
        }
        else if (strcmp(param->key, "abserrbound") == 0)
        {
            absErrBound = atof(param->value);
        }
        else if (strcmp(param->key, "relboundratio") == 0)
        {
            relBoundRatio = atof(param->value);
        }
        else
        {
            log_warn("An unknown SZ parameter: %s\n", param->key);
        }
    }
    
    
    bytes = SZ_compress_args (dtype, (void *) input_buff, &outsize,
                              errorBoundMode, absErrBound, relBoundRatio,
                              r[4], r[3], r[2], r[1], r[0]);
    
    unsigned char *raw_buff = (unsigned char*) bytes;
    int raw_size = outsize;
    //log_debug("=== SZ compress ===\n");
    log_debug("%s: %d\n", "SZ dtype", dtype);
    log_debug("%s: %d\n", "SZ out_size", raw_size);
    /*
    log_debug("%10s: %d %d %d %d %d ... %d %d %d %d %d\n", "out_buff",
              raw_buff[0], raw_buff[1], raw_buff[2], raw_buff[3], raw_buff[4],
              raw_buff[raw_size-5], raw_buff[raw_size-4], raw_buff[raw_size-3], raw_buff[raw_size-2], raw_buff[raw_size-1]);
    int sum = 0;
    for (i=0; i<raw_size; i++)
    {
        sum += raw_buff[i];
    }
    log_debug("%10s: %d\n", "sum", sum);
     */
    log_debug("%s: %d %d %d %d %d\n", "SZ dim", r[0], r[1], r[2], r[3], r[4]);
    //log_debug("===================\n");
    
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
