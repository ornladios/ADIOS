#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef ZLIB

#include "zlib.h"

int decompress_zlib_pre_allocated(const void* input_data, 
									const uint64_t input_len,
                                    void* output_data, 
									uint64_t* output_len)
{
    assert(input_data != NULL && input_len > 0 && output_data != NULL && output_len != NULL && *output_len > 0);

    uLongf dest_temp = *output_len;

    int z_rtn = uncompress((Bytef*)output_data, &dest_temp, (Bytef*)input_data, input_len);
    if(z_rtn != Z_OK)
    {
        // printf("zlib uncompress error %d\n", z_rtn);
        return -1;
    }

    *output_len = (uint64_t)dest_temp;

    return 0;
}

int adios_transform_zlib_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                    adios_transform_pg_read_request *pg_reqgroup)
{
    void *buf = malloc(pg_reqgroup->raw_var_length);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_zlib_subrequest_completed(adios_transform_read_request *reqgroup,
                                                            adios_transform_pg_read_request *pg_reqgroup,
                                                            adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}



adios_datablock * adios_transform_zlib_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                             adios_transform_pg_read_request *completed_pg_reqgroup)
{
    uint64_t raw_size = (uint64_t)completed_pg_reqgroup->raw_var_length;
    void* raw_buff = completed_pg_reqgroup->subreqs->data;
	
	uint64_t orig_size_meta = *((uint64_t*)reqgroup->transinfo->transform_metadata);
	char compress_succ = *((char*)(reqgroup->transinfo->transform_metadata + sizeof(uint64_t)));

    uint64_t orig_size = adios_get_type_size(reqgroup->transinfo->orig_type, "");
    int d = 0;
    for(d = 0; d < reqgroup->transinfo->orig_ndim; d++)
	{
        orig_size *= (uint64_t)(completed_pg_reqgroup->orig_varblock->count[d]);
	}
	
	if(orig_size_meta != orig_size)
	{
		printf("possible wrong data size or corrupted metadata\n");
	}

    void* orig_buff = malloc(orig_size);
	
	if(compress_succ == 1)	// compression is successful
	{
		int rtn = decompress_zlib_pre_allocated(raw_buff, raw_size, orig_buff, &orig_size);
		if(0 != rtn)
		{
			return NULL;
		}
	}
	else	// just copy the buffer since data is not compressed
	{
		memcpy(orig_buff, raw_buff, raw_size);
	}
	
    return adios_datablock_new(reqgroup->transinfo->orig_type,
                               completed_pg_reqgroup->timestep,
                               completed_pg_reqgroup->pg_bounds_sel,
                               orig_buff);
}

adios_datablock * adios_transform_zlib_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(zlib);

#endif

