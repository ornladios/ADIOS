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

enum ADIOS_ERRCODES adios_transform_identity_retrieve_subvolume(
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *copy_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness) {

    assert(var->characteristics[time_index].transform.transform_type == adios_transform_identity);

    // Read the variable in its entirety
    void *buf = read_delegate(var, time_index, 0, adios_transform_var_get_transformed_size(var, time_index), read_state);

    copy_subvolume_with_spec(global_out, buf, copy_spec, adios_transform_get_var_original_type(var), swap_endianness);

    return err_no_error;
}

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
ADIOS_VARCHUNK * adios_transform_identity_subrequest_completed(
                    adios_transform_read_reqgroup *reqgroup,
                    adios_transform_pg_reqgroup *pg_reqgroup,
                    adios_transform_read_subrequest *completed_subreq,
                    enum ADIOS_READ_RESULT_MODE mode) {
    return NULL;
}

ADIOS_VARCHUNK * adios_transform_identity_pg_reqgroup_completed(
        adios_transform_read_reqgroup *reqgroup,
        adios_transform_pg_reqgroup *completed_pg_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) {

    // If we are allowed to return a partial result, return this PG's data
    // Else, copy the PG's data into the global result buffer, and do nothing
    //   all data is present
    if (mode == adios_read_return_partial) {
        ADIOS_VARCHUNK *retchunk = (ADIOS_VARCHUNK *)malloc(sizeof(ADIOS_VARCHUNK));
        retchunk->varid = reqgroup->raw_varinfo->varid;
        retchunk->type = reqgroup->transinfo->orig_type;
        retchunk->sel = adios_copy_spec_to_dst_selection(completed_pg_reqgroup->pg_intersection_to_global_copyspec);
        retchunk->data = completed_pg_reqgroup->subreqs->data; // The first (and only) subrequest data buffer
        return retchunk;
    } else {
        copy_subvolume_with_spec(reqgroup->orig_data,					// Copy TO original buffer
                                 completed_pg_reqgroup->subreqs->data,	// Copy FROM buffer of first (and only) subrequest
                                 completed_pg_reqgroup->pg_intersection_to_global_copyspec,	// Copy USING the PG-to-global copy spec
                                 reqgroup->transinfo->orig_type,		// Copy elements of the original type
                                 reqgroup->swap_endianness);			// Swap endianness if needed
        return NULL;
    }
}

ADIOS_VARCHUNK * adios_transform_identity_reqgroup_completed(
        adios_transform_read_reqgroup *completed_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) {

    if (mode == adios_read_return_complete) {
        ADIOS_VARCHUNK *retchunk = (ADIOS_VARCHUNK *)malloc(sizeof(ADIOS_VARCHUNK));
        retchunk->varid = completed_reqgroup->raw_varinfo->varid;
        retchunk->type = completed_reqgroup->transinfo->orig_type;
        retchunk->sel = copy_selection(completed_reqgroup->orig_sel);
        retchunk->data = completed_reqgroup->orig_data;
        return retchunk;
    } else {
        // Mode is either noreturn, in which case there should be no return,
        // or return_partial, in which case all data was already returned
        return NULL;
    }
}
