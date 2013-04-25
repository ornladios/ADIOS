/*
 * Contains read-specific code for handling variable transforms in ADIOS
 *
 *  Created on: Jun 27, 2012
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_TRANSFORMS_READ_H_
#define ADIOS_TRANSFORMS_READ_H_

#include "public/adios_error.h"
#include "public/adios_types.h"
#include "adios_subvolume.h"
#include "transforms/adios_transforms_common.h"
#include "transforms/adios_transforms_reqgroup.h"

enum ADIOS_TRANSFORM_REQGROUP_RESULT_MODE {
    FULL_RESULT_MODE,
    PARTIAL_RESULT_MODE
};

enum ADIOS_TRANSFORM_REQGROUP_RESULT_MODE adios_transform_read_request_get_mode(adios_transform_read_request *reqgroup);

// Delegation functions
adios_transform_read_request * adios_transform_generate_read_reqgroup(const ADIOS_VARINFO *vi, const ADIOS_TRANSINFO* ti, const ADIOS_FILE *fp,
                                                                       const ADIOS_SELECTION *sel, int from_steps, int nsteps, const char *param, void *data);

/*
 * Processes a VARCHUNK just returned by the read layer against the given list of outstanding transform
 * read requests.
 */
void adios_transform_process_read_chunk(adios_transform_read_request **reqgroups_head, ADIOS_VARCHUNK ** chunk);

/*
 * Processes all data after a blocking read that has serviced all read requests,
 * completing all transform read requests.
 */
void adios_transform_process_all_reads(adios_transform_read_request **reqgroups_head);

#endif /* ADIOS_TRANSFORMS_READ_H_ */
