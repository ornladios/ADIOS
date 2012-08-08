/*
 * adios_transform_identity_read.c
 *
 * An implementation of the "identity" transform, which does nothing, but
 * exercises the transform framework for testing.
 *
 *  Created on: Jul 31, 2012
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "util.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

// Implementation of the "identity" transform, which does nothing, but
// exercises the transform framework for testing.

int adios_transform_identity_generate_read_subrequests(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *pg_reqgroup) {

    assert(reqgroup);
    assert(pg_reqgroup);
    // TODO: Optimize by only reading some of the PG where possible
    void *buf = malloc(pg_reqgroup->raw_var_length);
    adios_transform_read_subrequest *subreq = adios_transform_new_subreq_whole_pg(pg_reqgroup->raw_varblock, buf);
    adios_transform_pg_reqgroup_append_subreq(pg_reqgroup, subreq);

    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_identity_subrequest_completed(
                    adios_transform_read_reqgroup *reqgroup,
                    adios_transform_pg_reqgroup *pg_reqgroup,
                    adios_transform_read_subrequest *completed_subreq) {
    return NULL;
}

adios_datablock * adios_transform_identity_pg_reqgroup_completed(
        adios_transform_read_reqgroup *reqgroup,
        adios_transform_pg_reqgroup *completed_pg_reqgroup) {

    // Transfer ownership of the data buffer
    void *pg_data = completed_pg_reqgroup->subreqs->data;
    completed_pg_reqgroup->subreqs->data = NULL;

    return adios_datablock_new(reqgroup->transinfo->orig_type,
                               completed_pg_reqgroup->pg_bounds_global,
                               pg_data);
}

adios_datablock * adios_transform_identity_reqgroup_completed(
        adios_transform_read_reqgroup *completed_reqgroup) {
    return NULL;
}
