/*
 * adios_transform_zfp_read.c
 *
 * 	Author: Eric Suchyta
 * 	Contact: eric.d.suchyta@gmail.com
 */

#include "core/transforms/adios_transforms_hooks_read.h"

#ifdef ZFP

/* general C stuff */
#include <stdint.h>	// uint64_t
#include <stdio.h> 	// NULL, sprintf
#include <stdlib.h>	// NULL, malloc, free
#include <string.h>	// memcpy, strcmp, strlen  
#include <assert.h>


/* Were in the template included from ADIOS. Not necessarily sure if they're all strictly needed. */
#include "util.h"
#include "core/transforms/adios_transforms_reqgroup.h"


/* Extra ADIOS headers that weren't added in the template */
#include "core/adios_internals.h" 	// adios_get_type_size()


/* ZFP specific */
#include "adios_transform_zfp_common.h"
#include "zfp.h"



/* Basically giving an address and a size */
static void* zfp_read_metadata_var(const void* pos, size_t size, size_t* offset)
{
    *offset += size;
    return ((void*) pos + *offset - size);
}


/* Read string by iterating over byes */
static void read_metastring(char s[ZFP_STRSIZE], const void* pos, size_t* offset)
{
    int i;
    for (i=0; i<ZFP_STRSIZE; i++)
    {
        s[i] = *((char*)zfp_read_metadata_var(pos, 1, offset));
    }
    return;
}


/* Read each memory location and cast to the correct type */
static struct zfp_metadata* zfp_read_metadata(struct zfp_metadata* metadata, adios_transform_pg_read_request *completed_pg_reqgroup)
{
    const void* pos = completed_pg_reqgroup->transform_metadata;
    size_t offset = 0;

    metadata->usize = *((uint64_t*)zfp_read_metadata_var(pos, sizeof(uint64_t), &offset));
    metadata->csize = *((uint64_t*)zfp_read_metadata_var(pos, sizeof(uint64_t), &offset));
    metadata->cmode = *((uint*)zfp_read_metadata_var(pos, sizeof(uint), &offset));
    read_metastring(metadata->ctol, pos, &offset);
    read_metastring(metadata->name, pos, &offset);

    return metadata;
}



/* This is called in the main transform-level function.
 * In a nutshell: decompress array, using (undoing) the settings in the other args. Connect to the ADIOS buffer.
 */
static int zfp_decompression(struct zfp_buffer* zbuff, void* uarray, void* carray)
{
    zfp_initialize(uarray, zbuff);
    if (zbuff->error)
    {
        return 0;
    }

    /* The buffersize and compressed size aren't the same, this check doesn't make a lot of sense.
    if (zbuff->buffsize != csize)
    {
        log_warn("ZFP thinks compressed size should be %u" \
                "bytes. ADIOS thinks compressed size should be %" PRIu64 \
                "bytes. Likely corruption.\n", zbuff->buffsize, csize);
    }
    */


    zfp_streaming(zbuff, carray, 1, NULL);
    if (zbuff->error)
    {
        return 0;
    }
    /* Possibly add a check for output size later, if it matches what ADIOS though it should be.
    if (zusize != usize)
    {
        log_warn("ZFP thinks uncompressed size is %u" \
                "bytes. ADIOS thinks uncompressed size is %" PRIu64 \
                "bytes. Likely corruption.\n", zusize, usize);
    }
    */

    return 1;
}

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

	int i;
	int success;	// was (a piece of) the decompression okay


	/* Get the transform metadata */
	struct zfp_metadata* metadata = (struct zfp_metadata*) malloc(sizeof(struct zfp_metadata));	// allocate metadata
	metadata = zfp_read_metadata(metadata, completed_pg_reqgroup);


	/* Set up ZFP */
	void* cdata = completed_pg_reqgroup->subreqs->data;						// get the compressed data
	struct zfp_buffer* zbuff = (struct zfp_buffer*) malloc(sizeof(struct zfp_buffer));		// Handle zfp streaming
	init_zfp_buffer(zbuff, metadata->name);

	
	/* Get the data native to ADIOS (as opposed to the metadata which only the tranform plugin knows about) */
	uint64_t csize = (uint64_t)completed_pg_reqgroup->raw_var_length;
	uint64_t usize = adios_get_type_size(reqgroup->transinfo->orig_type, "");
	int d = 0;
	for(d = 0; d < reqgroup->transinfo->orig_ndim; d++)
	{
		usize *= (uint64_t)(completed_pg_reqgroup->orig_varblock->count[d]);
	}


	/* Do the metadata and ADIOS agree? */
	if (metadata->csize != csize)
	{
		log_warn("zfp processing variable %s: Metadata thinks compressed size is %" PRIu64 \
				"bytes. ADIOS thinks compressed size is %" PRIu64 \
				"bytes. Likely corruption.\n", zbuff->name, metadata->csize, csize);
	}

	if (metadata->usize != usize)
	{
		log_warn("zfp processing variable %s: Metadata thinks uncompressed size is %" PRIu64 \
				"bytes. ADIOS thinks uncompressed size is %" PRIu64 \
				"bytes. Likely corruption.\n", zbuff->name, metadata->usize, usize);
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


	/* Allocate the array we'll store the uncompressed data in */
	void* udata;
	udata = malloc(usize);
	if(!udata)
	{
		adios_error(err_no_memory, "Ran out of memory allocating uncompressed "
		        "buffer for ZFP transformation.\n");
		return NULL;
	}


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

