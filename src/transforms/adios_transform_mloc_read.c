#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "util.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef MLOC
#include "mloc.h"


int adios_transform_mloc_generate_read_subrequests(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *pg_reqgroup) {

    assert(reqgroup);
    assert(pg_reqgroup);

    return 0;
}

// Do nothing for individual subrequest
ADIOS_VARCHUNK * adios_transform_mloc_subrequest_completed(
                    adios_transform_read_reqgroup *reqgroup,
                    adios_transform_pg_reqgroup *pg_reqgroup,
                    adios_transform_read_subrequest *completed_subreq,
                    enum ADIOS_READ_RESULT_MODE mode) {
    return NULL;
}



ADIOS_VARCHUNK * adios_transform_mloc_pg_reqgroup_completed(
        adios_transform_read_reqgroup *reqgroup,
        adios_transform_pg_reqgroup *completed_pg_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) {

    return NULL;
}

ADIOS_VARCHUNK * adios_transform_mloc_reqgroup_completed(
        adios_transform_read_reqgroup *completed_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) {

    return NULL;
}


#else
DECLARE_TRANSFORM_READ_METHOD_UNIMPL(mloc);

#endif
