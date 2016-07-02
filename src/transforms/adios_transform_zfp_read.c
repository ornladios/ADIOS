/*
 * adios_transform_zfp_read.c
 *
 *  Created on: June 30, 2016
 *      Author: Eric Suchyta
 */
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>

#include "util.h"
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_reqgroup.h"
#include "core/adios_internals.h" // adios_get_type_size()
#include "adios_logger.h"

#ifdef ZFP

#include "zfp.h"
#include "adios_tranform_zfp_common.h"


int adios_transform_zfp_is_implemented (void) {return 1;}


// Do the default
int adios_transform_zfp_generate_read_subrequests(adios_transform_read_request *reqgroup, adios_transform_pg_read_request *pg_reqgroup)
{
    void *buf = malloc(pg_reqgroup->raw_var_length);
    assert(buf);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}


// Do nothing for individual subrequest
adios_datablock * adios_transform_bzip2_subrequest_completed(adios_transform_read_request *reqgroup, adios_transform_pg_read_request *pg_reqgroup, adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}


adios_datablock * adios_transform_bzip2_pg_reqgroup_completed(adios_transform_read_request *reqgroup, adios_transform_pg_read_request *completed_pg_reqgroup)
{
	
	/* Get the transform metadata */
	struct zfp_metadata metadata metadata = read_metadata(completed_pg_reqgroup);


	/* Get the data native to ADIOS (as opposed to the metadata which only the tranform plugin knows about) */
	uint64_t csize = (uint64_t)completed_pg_reqgroup->raw_var_length;
	uint64_t usize = adios_get_type_size(reqgroup->transinfo->orig_type, "");
	int d = 0;
	for(d = 0; d < reqgroup->transinfo->orig_ndim; d++)
	{
		usize *= (uint64_t)(completed_pg_reqgroup->orig_varblock->count[d]);
	}

	
	/* Do the metadata and ADIOS agree? */
	if(metadata.csize != csize)
	{
		log_warn("variable %s: Metadata thinks compressed size is %" PRIu64 \
				"bytes. ADIOS thinks compressed size is %" PRIu64 \
				"bytes. Likely corruption.\n", metadata.name, metadata.csize, csize);
	}

	if(metadata.usize != usize)
	{
		log_warn("variable %s: Metadata thinks uncompressed size is %" PRIu64 \
				"bytes. ADIOS thinks uncompressed size is %" PRIu64 \
				"bytes. Likely corruption.\n", metadata.name, metadata.usize, usize);
	}


	void* cdata = completed_pg_reqgroup->subreqs->data;
	void* udata = malloc(usize);
	if(!udata)
	{
		return NULL;
	}

       	/* possibly add check for successful compression eventually */ 
	int success = zfp_decompression(zbuff, var->array, outbuffer, *outsize, use_shared_buffer, fd);
        if(rtn != 0)
        {
            return NULL;
        }
	
	return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, udata);
}

adios_datablock * adios_transform_bzip2_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(bzip2);

#endif

