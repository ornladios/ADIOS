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
        //log_debug("%s: Fortran dimension enabled\n", "MGARD");
        ncol = (int) adios_get_dim_value(&d->dimension);
        nrow = (int) adios_get_dim_value(&d->next->dimension);
    }
    else
    {
        //log_debug("%s: C dimension enabled\n", "MGARD");
        nrow = (int) adios_get_dim_value(&d->dimension);
        ncol = (int) adios_get_dim_value(&d->next->dimension);
    }

    struct adios_transform_spec_kv_pair* param;
    int i = 0;
    if (adios_verbose_level>7) log_debug("param_count: %d\n", var->transform_spec->param_count);
    for (i=0; i<var->transform_spec->param_count; i++)
    {
        param = &(var->transform_spec->params[i]);
        if (adios_verbose_level>7) log_debug("param: %s %s\n", param->key, param->value);

        if (!strcmp(param->key, "tol") || !strcmp(param->key, "accuracy"))
        {
            tol = atof(param->value);
        }
    }

    unsigned char* mgard_comp_buff;
    mgard_comp_buff = mgard_compress(iflag, input_buff, &out_size,  nrow,  ncol, &tol );
    //out_size = 15671;
    //mgard_comp_buff = malloc(out_size);
    //memset(mgard_comp_buff, 0xFF, out_size);

    log_debug("%s: %d,%d\n", "MGARD now,ncol", nrow, ncol);
    log_debug("%s: %g\n", "MGARD tol", tol);
    log_debug("%s: %d\n", "MGARD out_size", out_size);
    log_debug("%s: %p\n", "MGARD output buffer", mgard_comp_buff);

    FILE *qfile;
    char fname[80];
    sprintf(fname, "input_type%d_%dx%d.dat", iflag, nrow, ncol);
    log_debug("%s: %s\n", "MGARD save", fname);
    qfile = fopen (fname, "wb");
    fwrite (input_buff, 1, nrow*ncol*8, qfile);
    fclose(qfile);

    double *v = (double*) input_buff;
    log_debug("%s: %g %g %g %g %g ... %g %g %g %g %g\n", "MGARD input_buff",
              v[0], v[1], v[2], v[3], v[4],
              v[nrow*ncol-5], v[nrow*ncol-4], v[nrow*ncol-3], v[nrow*ncol-2], v[nrow*ncol-1]);

    // Output
    uint64_t output_size = (uint64_t) out_size/* Compute how much output size we need */;
    void* output_buff;

    log_debug("%s: %ld\n", "MGARD output_size", output_size);
    log_debug("%s: %d\n", "MGARD use_shared_buffer", use_shared_buffer);
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
