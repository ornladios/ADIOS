#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "util.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef COMPRESS
#include "compress.h"


int adios_transform_compress_generate_read_subrequests(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *pg_reqgroup) 
{

    assert(reqgroup);
    assert(pg_reqgroup);
	
	void *buf = malloc(pg_reqgroup->raw_var_length);
    adios_transform_read_subrequest *subreq = adios_transform_new_subreq_whole_pg(pg_reqgroup->raw_varblock, buf);
    adios_transform_pg_reqgroup_append_subreq(pg_reqgroup, subreq);

    return 0;
}

// Do nothing for individual subrequest
ADIOS_VARCHUNK * adios_transform_compress_subrequest_completed(
                    adios_transform_read_reqgroup *reqgroup,
                    adios_transform_pg_reqgroup *pg_reqgroup,
                    adios_transform_read_subrequest *completed_subreq,
                    enum ADIOS_READ_RESULT_MODE mode) 
{

    return NULL;
}



ADIOS_VARCHUNK * adios_transform_compress_pg_reqgroup_completed(
        adios_transform_read_reqgroup *reqgroup,
        adios_transform_pg_reqgroup *completed_pg_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) 
{
	uint32_t compressed_len = (uint32_t)completed_pg_reqgroup->raw_var_length;
	void* compressed_buff = completed_pg_reqgroup->subreqs->data;
	
	uint32_t decompressed_len = adios_get_type_size(reqgroup->raw_varinfo->type, "");
	int d = 0;
	for(d = 0; d < reqgroup->raw_varinfo->ndim; d++)
	{
		decompressed_len *= (uint32_t)(completed_pg_reqgroup->orig_varblock->count[d]);
	}
	
	void* decompressed_buff = malloc(decompressed_len);	
	
	enum COMPRESS_TYPE compress_type_flag = 0;
	int rtn = 0;
	
	switch(compress_type_flag)
    {
        case compress_type_zlib:
            rtn = decompress_zlib_pre_allocated(compressed_buff, compressed_len, decompressed_buff, &decompressed_len);
        break;

        case compress_type_bzlib2:
            rtn = decompress_bzlib2_pre_allocated(compressed_buff, compressed_len, decompressed_buff, &decompressed_len);
        break;

        case compress_type_szip:
            // rtn = compress_szip_pre_allocated(input_buff, input_size, output_buff, &output_size);
            rtn = decompress_bzlib2_pre_allocated(compressed_buff, compressed_len, decompressed_buff, &decompressed_len);
        break;

        default: // default: do not do compression, just copy the buffer
		{
            // memcpy(output_buff, input_buff, input_size);
            // output_size = input_size;
            // compress_type_flag = compress_type_none;
            rtn = -1;
		}
        break;

    }
	
	if(rtn != 0)
	{
		if(decompressed_buff)
		{
			free(decompressed_buff);
			decompressed_buff = NULL;
		}
		return NULL;
	}
	
	if (mode == adios_read_return_partial) {
        ADIOS_VARCHUNK *retchunk = (ADIOS_VARCHUNK *)malloc(sizeof(ADIOS_VARCHUNK));
        retchunk->varid = reqgroup->raw_varinfo->varid;
        retchunk->type = reqgroup->transinfo->orig_type;
        retchunk->sel = adios_copy_spec_to_dst_selection(completed_pg_reqgroup->pg_intersection_to_global_copyspec);
        retchunk->data = decompressed_buff;
        return retchunk;
    } else {
        copy_subvolume_with_spec(reqgroup->orig_data,					// Copy TO original buffer
                                 decompressed_buff,							// Copy FROM decoded buffer
                                 completed_pg_reqgroup->pg_intersection_to_global_copyspec,	// Copy USING the PG-to-global copy spec
                                 reqgroup->transinfo->orig_type,		// Copy elements of the original type
                                 reqgroup->swap_endianness);			// Swap endianness if needed
        return NULL;
    }
	
    return NULL;
}

ADIOS_VARCHUNK * adios_transform_compress_reqgroup_completed(
        adios_transform_read_reqgroup *completed_reqgroup,
        enum ADIOS_READ_RESULT_MODE mode) 
{
	
	if (mode == adios_read_return_complete) 
	{
        ADIOS_VARCHUNK *retchunk = (ADIOS_VARCHUNK *)malloc(sizeof(ADIOS_VARCHUNK));
        retchunk->varid = completed_reqgroup->raw_varinfo->varid;
        retchunk->type = completed_reqgroup->transinfo->orig_type;
        retchunk->sel = copy_selection(completed_reqgroup->orig_sel);
        retchunk->data = completed_reqgroup->orig_data;
        return retchunk;
    } 
	else 
	{
        // Mode is either noreturn, in which case there should be no return,
        // or return_partial, in which case all data was already returned
        return NULL;
    }
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(compress);

#endif
