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

#ifdef HAVE_SZ

#include "sz.h"

typedef unsigned int uint;

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

    sz_params sz;
    memset(&sz, 0, sizeof(sz_params));
    sz.max_quant_intervals = 65536;
    sz.quantization_intervals = 0;
    sz.dataEndianType = LITTLE_ENDIAN_DATA;
    sz.sysEndianType = LITTLE_ENDIAN_DATA;
    sz.sol_ID = SZ;
    sz.layers = 1;
    sz.sampleDistance = 100;
    sz.predThreshold = 0.99;
    sz.offset = 0;
    sz.szMode = SZ_BEST_COMPRESSION; //SZ_BEST_SPEED; //SZ_BEST_COMPRESSION;
    sz.gzipMode = 1;
    sz.errorBoundMode = ABS;
    sz.absErrBound = 1E-4;
    sz.relBoundRatio = 1E-3;
    sz.pw_relBoundRatio = 1E-5;
    sz.segment_size = 32;

    unsigned char *bytes;
    int outsize;
    int r[5] = {0,0,0,0,0};

    /* SZ parameters */
    int use_configfile = 0;
    char *sz_configfile = NULL;
    struct adios_transform_spec_kv_pair* param;
    int i = 0;
    if (adios_verbose_level>7) log_debug("param_count: %d\n", var->transform_spec->param_count);
    for (i=0; i<var->transform_spec->param_count; i++)
    {
        param = &(var->transform_spec->params[i]);
        if (adios_verbose_level>7) log_debug("param: %s %s\n", param->key, param->value);
        if (strcmp(param->key, "init") == 0)
        {
            use_configfile = 1;
            sz_configfile = strdup(param->value);
        }
        else if (strcmp(param->key, "max_quant_intervals") == 0)
        {
            sz.max_quant_intervals = atoi(param->value);
        }
        else if (strcmp(param->key, "quantization_intervals") == 0)
        {
            sz.quantization_intervals = atoi(param->value);
        }
        else if (strcmp(param->key, "dataEndianType") == 0)
        {
            sz.dataEndianType = atoi(param->value);
        }
        else if (strcmp(param->key, "sysEndianType") == 0)
        {
            sz.sysEndianType = atoi(param->value);
        }
        else if (strcmp(param->key, "sol_ID") == 0)
        {
            sz.sol_ID = atoi(param->value);
        }
        else if (strcmp(param->key, "layers") == 0)
        {
            sz.layers = atoi(param->value);
        }
        else if (strcmp(param->key, "sampleDistance") == 0)
        {
            sz.sampleDistance = atoi(param->value);
        }
        else if (strcmp(param->key, "predThreshold") == 0)
        {
            sz.predThreshold = atof(param->value);
        }
        else if (strcmp(param->key, "offset") == 0)
        {
            sz.offset = atoi(param->value);
        }
        else if (strcmp(param->key, "szMode") == 0)
        {
            int szMode = SZ_BEST_SPEED;
            if (strcmp(param->value, "SZ_BEST_SPEED") == 0)
            {
              szMode = SZ_BEST_SPEED;
            }
            else if (strcmp(param->value, "SZ_BEST_COMPRESSION") == 0)
            {
              szMode = SZ_BEST_COMPRESSION;
            }
            else if (strcmp(param->value, "SZ_DEFAULT_COMPRESSION") == 0)
            {
              szMode = SZ_DEFAULT_COMPRESSION;
            }
            else
            {
              log_warn("An unknown szMode: %s\n", param->value);
            }
            sz.szMode = szMode;
        }
        else if (strcmp(param->key, "gzipMode") == 0)
        {
            sz.gzipMode = atoi(param->value);
        }
        else if (strcmp(param->key, "errorBoundMode") == 0)
        {
            int errorBoundMode = ABS;
            if (strcmp(param->value, "ABS") == 0)
            {
              errorBoundMode = ABS;
            }
            else if (strcmp(param->value, "REL") == 0)
            {
              errorBoundMode = REL;
            }
            else if (strcmp(param->value, "ABS_AND_REL") == 0)
            {
              errorBoundMode = ABS_AND_REL;
            }
            else if (strcmp(param->value, "ABS_OR_REL") == 0)
            {
              errorBoundMode = ABS_OR_REL;
            }
            else if (strcmp(param->value, "PW_REL") == 0)
            {
              errorBoundMode = PW_REL;
            }
            else
            {
              log_warn("An unknown errorBoundMode: %s\n", param->value);
            }
            sz.errorBoundMode = errorBoundMode;
        }
        else if (strcmp(param->key, "absErrBound") == 0)
        {
            sz.absErrBound = atof(param->value);
        }
        else if (strcmp(param->key, "relBoundRatio") == 0)
        {
            sz.relBoundRatio = atof(param->value);
        }
        else if (strcmp(param->key, "pw_relBoundRatio") == 0)
        {
            sz.pw_relBoundRatio = atof(param->value);
        }
        else if (strcmp(param->key, "segment_size") == 0)
        {
            sz.segment_size = atoi(param->value);
        }
        else if (!strcmp(param->key, "abs") || !strcmp(param->key, "absolute") || !strcmp(param->key, "accuracy"))
        {
            sz.errorBoundMode = ABS;
            sz.absErrBound = atof(param->value);
        }
        else if (!strcmp(param->key, "rel") || !strcmp(param->key, "relative"))
        {
            sz.errorBoundMode = REL;
            sz.relBoundRatio = atof(param->value);
        }
        else
        {
            log_warn("An unknown SZ parameter: %s\n", param->key);
        }
    }

    if (use_configfile)
    {
        log_debug("%s: %s\n", "SZ config", sz_configfile);
        SZ_Init(sz_configfile);
        free(sz_configfile);
    }
    else
    {
        if (adios_verbose_level>7)
        {
            log_debug("%s: %d\n", "sz.max_quant_intervals", sz.max_quant_intervals);
            log_debug("%s: %d\n", "sz.quantization_intervals", sz.quantization_intervals);
            log_debug("%s: %d\n", "sz.dataEndianType", sz.dataEndianType);
            log_debug("%s: %d\n", "sz.sysEndianType", sz.sysEndianType);
            log_debug("%s: %d\n", "sz.sol_ID", sz.sol_ID);
            log_debug("%s: %d\n", "sz.layers", sz.layers);
            log_debug("%s: %d\n", "sz.sampleDistance", sz.sampleDistance);
            log_debug("%s: %g\n", "sz.predThreshold", sz.predThreshold);
            log_debug("%s: %d\n", "sz.offset", sz.offset);
            log_debug("%s: %d\n", "sz.szMode", sz.szMode);
            log_debug("%s: %d\n", "sz.gzipMode", sz.gzipMode);
            log_debug("%s: %d\n", "sz.errorBoundMode", sz.errorBoundMode);
            log_debug("%s: %g\n", "sz.absErrBound", sz.absErrBound);
            log_debug("%s: %g\n", "sz.relBoundRatio", sz.relBoundRatio);
            log_debug("%s: %g\n", "sz.pw_relBoundRatio", sz.pw_relBoundRatio);
            log_debug("%s: %d\n", "sz.segment_size", sz.segment_size);
        }
        SZ_Init_Params(&sz);
    }

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
            return -1;
            break;
    }

    // Get dimension info
    struct adios_dimension_struct* d = var->pre_transform_dimensions;
    int ndims = (uint) count_dimensions(d);
    //log_debug("ndims: %d\n", ndims);
    if (ndims > 5)
    {
        adios_error(err_transform_failure, "No more than 5 dimension is supported.\n");
        return -1;
    }

    // r[0] is the fastest changing dimension and r[4] is the lowest changing dimension
    // In C, r[0] is the last dimension. In Fortran, r[0] is the first dimension
    for (i=0; i<ndims; i++)
    {
        uint dsize = (uint) adios_get_dim_value(&d->dimension);
        if (fd->group->adios_host_language_fortran == adios_flag_yes)
            r[i] = dsize;
        else
            r[ndims-i-1] = dsize;
        d = d->next;
    }

    bytes = SZ_compress (dtype, (void *) input_buff, &outsize,
                         r[4], r[3], r[2], r[1], r[0]);

    unsigned char *raw_buff = (unsigned char*) bytes;

    int raw_size = outsize;
    //log_debug("=== SZ compress ===\n");
    log_debug("%s: %d\n", "SZ dtype", dtype);
    log_debug("%s: %d\n", "SZ out_size", raw_size);
    /*
    log_debug("%s: %d %d %d %d %d ... %d %d %d %d %d\n", "SZ out_buff",
              raw_buff[0], raw_buff[1], raw_buff[2], raw_buff[3], raw_buff[4],
              raw_buff[raw_size-5], raw_buff[raw_size-4], raw_buff[raw_size-3], raw_buff[raw_size-2], raw_buff[raw_size-1]);
    int sum = 0;
    for (i=0; i<raw_size; i++)
    {
        sum += raw_buff[i];
    }
    log_debug("%s: %d\n", "SZ sum", sum);
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
    SZ_Finalize();
    return 1;
}

#else

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(sz)

#endif
