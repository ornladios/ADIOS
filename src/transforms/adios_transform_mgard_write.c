/*
 * adios_transform_mgard_write.c
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

#ifdef HAVE_MGARD

#include "mgard_capi.h"

typedef unsigned int uint;

uint16_t adios_transform_mgard_get_metadata_size(struct adios_transform_spec *transform_spec)
{
    //log_debug("function: %s\n", __FUNCTION__);
    return 0; // Set amount of transform-internal metadata space to allocate
}

void adios_transform_mgard_transformed_size_growth(const struct adios_var_struct *var, const struct adios_transform_spec *transform_spec,
                                                uint64_t *constant_factor, double *linear_factor, double *capped_linear_factor, uint64_t *capped_linear_cap)
{
    //log_debug("function: %s\n", __FUNCTION__);
}

int adios_transform_mgard_apply(struct adios_file_struct *fd,
                             struct adios_var_struct *var,
                             uint64_t *transformed_len,
                             int use_shared_buffer,
                             int *wrote_to_shared_buffer)
{
    //log_debug("function: %s\n", __FUNCTION__);
    //log_debug("use_shared_buffer: %d\n", use_shared_buffer);

    int iflag = 1; //0 -> float, 1 -> double
    int nrow, ncol;
    double tol;
    int out_size;

    // Get type info
    switch (var->pre_transform_type)
    {
        case adios_double:
            iflag = 1;
            break;
        case adios_real:
            iflag = 0;
            break;
        default:
            adios_error(err_transform_failure, "No supported data type\n");
            return -1;
            break;
    }

    // Get the input data and data length
    const uint64_t input_size = adios_transform_get_pre_transform_var_size(var);
    const void *input_buff = var->data;

    // Get dimension info
    struct adios_dimension_struct* d = var->pre_transform_dimensions;
    int ndims = (uint) count_dimensions(d);
    //log_debug("ndims: %d\n", ndims);
    if (ndims != 2)
    {
        adios_error(err_transform_failure, "Support only 2 dimension.\n");
        return -1;
    }

    if (fd->group->adios_host_language_fortran == adios_flag_yes)
    {
        nrow = (int) adios_get_dim_value(&d->dimension);
        ncol = (int) adios_get_dim_value(&d->next->dimension);
    }
    else
    {
        ncol = (int) adios_get_dim_value(&d->dimension);
        nrow = (int) adios_get_dim_value(&d->next->dimension);
    }
    log_debug("nrow,ncol: %d,%d\n", ncol, nrow);

    struct adios_transform_spec_kv_pair* param;
    int i = 0;
    if (adios_verbose_level>7) log_debug("param_count: %d\n", var->transform_spec->param_count);
    for (i=0; i<var->transform_spec->param_count; i++)
    {
        param = &(var->transform_spec->params[i]);
        if (adios_verbose_level>7) log_debug("param: %s %s\n", param->key, param->value);

        if (strcmp(param->key, "tol") == 0)
        {
            tol = atof(param->value);
        }
    }
    log_debug("tol: %g\n", tol);

    unsigned char* mgard_comp_buff;
    mgard_comp_buff = mgard_compress(iflag, input_buff, &out_size,  nrow,  ncol, &tol );
    log_debug("out_size: %d\n", out_size);

    // Output
    uint64_t output_size = (uint64_t) out_size/* Compute how much output size we need */;
    void* output_buff;

    if (use_shared_buffer) {
        // If shared buffer is permitted, serialize to there
        assert(shared_buffer_reserve(fd, output_size));

        // Write directly to the shared buffer
        output_buff = fd->buffer + fd->offset;
        memcpy(output_buff, mgard_comp_buff, (size_t)output_size);
    } else { // Else, fall back to var->adata memory allocation
        output_buff = mgard_comp_buff;
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

DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(mgard)

#endif
