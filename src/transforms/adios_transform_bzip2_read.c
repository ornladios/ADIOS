#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include "util.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef BZIP2

#include "bzlib.h"

int decompress_bzip2_pre_allocated(const void* input_data, const uint64_t input_len, 
									void* output_data, uint64_t* output_len)
{
	// bzip2 only support input size of 32 bit integer

	assert(input_data != NULL 
			&& input_len > 0 && input_len <= UINT_MAX 
			&& output_data != NULL 
			&& output_len != NULL && *output_len > 0 && *output_len < UINT_MAX);

	unsigned int input_len_32 = (unsigned int)input_len;
	unsigned int output_len_32 = (unsigned int)(*output_len);
	
	int bz_rtn = BZ2_bzBuffToBuffDecompress((char*)output_data, &output_len_32,
												(char*)input_data, input_len_32,
												0, 0 );
	
	if(bz_rtn != BZ_OK)
	{
		printf("BZ2_bzBuffToBuffDecompress error %d\n", bz_rtn);
		return -1;
	}
	
	*output_len = output_len_32;
	return 0;
}

int adios_transform_bzip2_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                       adios_transform_pg_read_request *pg_reqgroup) 
{
	assert(reqgroup && pg_reqgroup);

	void *buf = malloc(pg_reqgroup->raw_var_length);

	// printf("[adios_transform_compress_generate_read_subrequests] raw_var_length %d %d %d %d %d\n", 
			// pg_reqgroup->raw_var_length, pg_reqgroup->raw_varblock->start[0], pg_reqgroup->raw_varblock->count[0], 
			// pg_reqgroup->raw_varblock->start[1], pg_reqgroup->raw_varblock->count[1]);

	// adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new(pg_reqgroup->raw_varblock, buf);
    // adios_transform_raw_read_request_append(pg_reqgroup, subreq);
	
	adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup->raw_varblock, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);

	return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_bzip2_subrequest_completed(adios_transform_read_request *reqgroup,
																adios_transform_pg_read_request *pg_reqgroup,
																adios_transform_raw_read_request *completed_subreq) 
{
    return NULL;
}



adios_datablock * adios_transform_bzip2_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
																adios_transform_pg_read_request *completed_pg_reqgroup) 
{
	uint64_t compressed_len = (uint64_t)completed_pg_reqgroup->raw_var_length;
	void* compressed_buff = completed_pg_reqgroup->subreqs->data;
	
	uint64_t decompressed_len_test = adios_get_type_size(reqgroup->raw_varinfo->type, "");
	int d = 0;
	for(d = 0; d < reqgroup->raw_varinfo->ndim; d++)
	{
		decompressed_len_test *= (uint64_t)(completed_pg_reqgroup->orig_varblock->count[d]);
	}
	
	// retrieve the original buffer length from metadata
	uint64_t decompressed_len = *((uint64_t*)(reqgroup->transinfo->transform_metadata));
	void* decompressed_buff = malloc(decompressed_len);	
	
	int rtn = decompress_bzip2_pre_allocated(compressed_buff, compressed_len, decompressed_buff, &decompressed_len);
	if(rtn != 0)
	{
		if(decompressed_buff)
		{
			free(decompressed_buff);
			decompressed_buff = NULL;
		}
		return NULL;
	}
	
	return adios_datablock_new(reqgroup->transinfo->orig_type,								
                               completed_pg_reqgroup->timestep,
                               completed_pg_reqgroup->pg_bounds_sel,
                               decompressed_buff);
}

adios_datablock * adios_transform_bzip2_reqgroup_completed(adios_transform_read_request *completed_reqgroup) 
{
	return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(bzip2);

#endif

