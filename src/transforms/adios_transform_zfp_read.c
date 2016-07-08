/*
 * adios_transform_zfp_read.c
 *
 * 	Author: Eric Suchyta
 * 	Contact: eric.d.suchyta@gmail.com
 */

#ifdef ZFP


/* general C stuff */
#include <stdint.h>	// uint64_t
#include <stdio.h> 	// NULL, sprintf
#include <stdlib.h>	// NULL, malloc, free
#include <string.h>	// memcpy, strcmp, strlen  
#include <assert.h>


/* Were in the template included from ADIOS. Not necessarily sure if they're all strictly needed. */
#include "util.h"
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_reqgroup.h"


/* Extra ADIOS headers that weren't added in the template */
#include "core/adios_internals.h" 	// adios_get_type_size()


/* ZFP specific */
#include "adios_transform_zfp_common.h"
#include "zfp.h"


/* ZFP is installed */
int adios_transform_zfp_is_implemented (void) {return 1;}


/* Kept the default. I think this is piecing together how to read a "block" from smaller subrequests*/
int adios_transform_zfp_generate_read_subrequests(adios_transform_read_request *reqgroup, adios_transform_pg_read_request *pg_reqgroup)
{
    void *buf = malloc(pg_reqgroup->raw_var_length);
    assert(buf);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}


/* Kept the default. Template says "Do nothing for individual subrequest" */
adios_datablock * adios_transform_zfp_subrequest_completed(adios_transform_read_request *reqgroup, 
		adios_transform_pg_read_request *pg_reqgroup, adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}


/* The main work is being done here. I think this happens when a "block" of data gets returned. */
adios_datablock * adios_transform_zfp_pg_reqgroup_completed(adios_transform_read_request *reqgroup, 
		adios_transform_pg_read_request *completed_pg_reqgroup)
{

	/* Set up ZFP */
	int i;
	int success;		// was (a piece of) the decompression okay


	struct zfp_buffer* zbuff = (struct zfp_buffer*) malloc(sizeof(struct zfp_buffer));		// Handle zfp streaming
	struct zfp_metadata* metadata = (struct zfp_metadata*) malloc(sizeof(struct zfp_metadata));	// allocate metadata
	metadata = zfp_read_metadata(metadata, completed_pg_reqgroup);

	void* cdata = completed_pg_reqgroup->subreqs->data;			// get the compressed data
	void* udata;								// decompress into this

	
	/* Get the transform metadata */
	strcpy(zbuff->name, metadata->name);


	/* Get the data native to ADIOS (as opposed to the metadata which only the tranform plugin knows about) */
	uint64_t csize = (uint64_t)completed_pg_reqgroup->raw_var_length;
	uint64_t usize = adios_get_type_size(reqgroup->transinfo->orig_type, "");
	int d = 0;
	for(d = 0; d < reqgroup->transinfo->orig_ndim; d++)
	{
		usize *= (uint64_t)(completed_pg_reqgroup->orig_varblock->count[d]);
	}


	/* Allocate the array we'll store the uncompressed data in */
	udata = malloc(usize);			// allocate space for uncompressed data
	if(!udata)
	{
		sprintf(zbuff->msg, "Ran out of memory allocating uncompressed buffer\n");
		zfp_error(zbuff);
		return NULL;
	}

	
	/* Do the metadata and ADIOS agree? */
	if (metadata->csize != csize)
	{
		sprintf(zbuff->msg, "Metadata thinks compressed size is %" PRIu64 \
				"bytes. ADIOS thinks compressed size is %" PRIu64 \
				"bytes. Likely corruption.\n", metadata->csize, csize);
		zfp_warn(zbuff);
	}

	if (metadata->usize != usize)
	{
		sprintf(zbuff->msg, "Metadata thinks uncompressed size is %" PRIu64 \
				"bytes. ADIOS thinks uncompressed size is %" PRIu64 \
				"bytes. Likely corruption.\n", metadata->usize, usize);
		zfp_warn(zbuff);
	}


	/* zfp datatype */
	success = zfp_get_datatype(zbuff, reqgroup->transinfo->orig_type);
	if (!success)
	{
		return NULL;
	}


	/* dimensionality */
	zbuff->ndims = (uint) reqgroup->transinfo->orig_ndim;
	zbuff->dims = (uint*) malloc(zbuff->ndims*sizeof(uint));
	for (i=0; i<zbuff->ndims; i++)
	{
		zbuff->dims[i] = (uint) reqgroup->transinfo->orig_dims[i];
	}

	
	/* mode */
	zbuff->mode = metadata->cmode;
	strcpy(zbuff->ctol, metadata->ctol);


       	/* possibly add check for successful decompression eventually */ 
	success = zfp_decompression(zbuff, udata, cdata);
        if(!success)
        {
            return NULL;
        }

	free(zbuff);
	free(metadata);
	
	return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, udata);
}


/* Kept the default. Template says "Do nothing for the full read request complete (typical)" */
adios_datablock * adios_transform_zfp_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(zfp);

#endif

