#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "adios_logger.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef HAVE_MGARD

#include "mgard_capi.h"

int adios_transform_mgard_is_implemented (void) {return 1;}

int adios_transform_mgard_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                 adios_transform_pg_read_request *pg_reqgroup)
{
    //log_debug("function: %s\n", __FUNCTION__);
    void *buf = malloc(pg_reqgroup->raw_var_length);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_mgard_subrequest_completed(adios_transform_read_request *reqgroup,
                                                          adios_transform_pg_read_request *pg_reqgroup,
                                                          adios_transform_raw_read_request *completed_subreq)
{
    //log_debug("function: %s\n", __FUNCTION__);
    return NULL;
}



adios_datablock * adios_transform_mgard_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                           adios_transform_pg_read_request *completed_pg_reqgroup)
{
    //log_debug("function: %s\n", __FUNCTION__);
    int iflag = 1; //0 -> float, 1 -> double
    int nrow, ncol;
    int raw_size = (int) completed_pg_reqgroup->raw_var_length;
    unsigned char *raw_buff = completed_pg_reqgroup->subreqs->data;

    // Get type info
    switch (reqgroup->transinfo->orig_type)
    {
        case adios_double:
            iflag = 1;
            break;
        case adios_real:
            iflag = 0;
            break;
        default:
            adios_error(err_transform_failure, "No supported data type\n");
            return NULL;
            break;
    }

    void* mgard_out_buff;
    // Get dimension info
    int ndims = reqgroup->transinfo->orig_ndim;
    if (ndims != 2)
    {
        mgard_out_buff = (unsigned char *) malloc (raw_size);
        memcpy(mgard_out_buff, raw_buff, raw_size);
    }
    else
    {
        nrow = (int) (completed_pg_reqgroup->orig_varblock->count[0]);
        ncol = (int) (completed_pg_reqgroup->orig_varblock->count[1]);

        log_debug("%s: %d,%d\n", "MGARD now,ncol", nrow, ncol);
        log_debug("%s: %d\n", "MGARD out_size", raw_size);

        mgard_out_buff = mgard_decompress(iflag, raw_buff, raw_size,  nrow,  ncol); 
    }

    return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, mgard_out_buff);
}

// Do nothing for the full read request complete (typical)
adios_datablock * adios_transform_mgard_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    //log_debug("function: %s\n", __FUNCTION__);
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(mgard);

#endif
