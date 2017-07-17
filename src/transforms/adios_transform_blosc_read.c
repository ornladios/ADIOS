/*
 * adios_transform_blosc_read.c
 *
 * 	Author: Rene Widera
 * 	Contact: r.widera@hzdr.de
 *
 *  This code based on `adios_transform_zlib_read.c`.
 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "core/util.h"
#include "core/adios_logger.h"
#include "core/transforms/adios_transforms_hooks_read.h"
#include "core/transforms/adios_transforms_reqgroup.h"
#include "core/adios_internals.h" // adios_get_type_size()
#include "core/adios_endianness.h"

#ifdef BLOSC

#include "blosc.h"
#include "adios_transform_blosc_common.h"

int adios_transform_blosc_is_implemented (void) {return 1;}

int adios_transform_blosc_decompress( 
                                  const void* input_data,
                                  void* output_data,
                                  const adiosBloscSize_t max_output_len,
                                  adiosBloscSize_t* decoded_bytes)
{
    assert(input_data != NULL &&  output_data != NULL && max_output_len > 0 && decoded_bytes != 0);

    
    adiosBloscSize_t result = blosc_decompress(input_data, output_data, max_output_len);

    *decoded_bytes = 0;
    if( result > 0)
    {
        *decoded_bytes = result;
    }

    return result <= 0;
}

int adios_transform_blosc_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                    adios_transform_pg_read_request *pg_reqgroup)
{
    void *buf = malloc(pg_reqgroup->raw_var_length);
    assert(buf);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_blosc_subrequest_completed(adios_transform_read_request *reqgroup,
                                                            adios_transform_pg_read_request *pg_reqgroup,
                                                            adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}

/* convert a 4 byte size from the blosc meta data to adiosBloscSize_t
 * 
 * big and little endian is observed
 */
static adiosBloscSize_t adios_transform_blosc_endian_convert(const char* const input)
{
    int32_t result = *((int32_t*)input);
    /* check for big and little endian */
    int check = 1;  
    /* first byte of big enian is zero */
    if ( *((char *)&check) != 1) {
        /* big endian need swap */
        swap_32_ptr(&result);
    }
    return result;
}

adios_datablock * adios_transform_blosc_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                             adios_transform_pg_read_request *completed_pg_reqgroup)
{
    uint64_t input_size = (uint64_t)completed_pg_reqgroup->raw_var_length;
    char* input_buff = (char*)(completed_pg_reqgroup->subreqs->data);

    adiosBloscSize_t num_chunks = *((adiosBloscSize_t*)completed_pg_reqgroup->transform_metadata);
    adiosBloscSize_t compressed_size_last_chunk = *((adiosBloscSize_t*)(completed_pg_reqgroup->transform_metadata + sizeof(adiosBloscSize_t)));

    int is_compressed = 1;
    if( num_chunks == 0 && compressed_size_last_chunk == 0 )
        is_compressed = 0;
    
    uint64_t uncompressed_size = adios_get_type_size(reqgroup->transinfo->orig_type, "");
    int d = 0;
    for(d = 0; d < reqgroup->transinfo->orig_ndim; d++)
    {
        uncompressed_size *= (uint64_t)(completed_pg_reqgroup->orig_varblock->count[d]);
    }

    char* output_buff = (char*)malloc(uncompressed_size);
    if(!output_buff)
    {
        return NULL;
    }
    
   
    uint64_t input_offset = 0;
    uint64_t actual_output_size = 0;

    adiosBloscSize_t chunk = 0;
    for(; (chunk < num_chunks || input_offset < input_size) && is_compressed; ++chunk)
    {
        /* move over the size of the compressed data */
        const char* in_ptr = input_buff + input_offset;
        
        /** read the size of the compress block from the blosc meta data
         * 
         * blosc meta data format (all little endian):
         *   - 1 byte blosc format version
         *   - 1 byte blosclz format version
         *   - 1 byte flags
         *   - 1 byte typesize
         *   - 4 byte uncompressed data size
         *   - 4 byte block size
         *   - 4 byte compressed data size
         * 
         * we need only the compressed size ( source address + 12 byte)
         */
        adiosBloscSize_t max_input_size = adios_transform_blosc_endian_convert(in_ptr+12);

        char* out_ptr = output_buff + actual_output_size;
        adiosBloscSize_t max_output_size = ADIOS_BLOSC_MAX_INPUT_SIZE;
        // handle the last not full chunk
        if (chunk >= num_chunks)
        {
            max_output_size = uncompressed_size - actual_output_size;
        }
        
        adiosBloscSize_t max_chunk_size = max_output_size + BLOSC_MAX_OVERHEAD;
        adiosBloscSize_t decoded_bytes = 0;
        int rtn = adios_transform_blosc_decompress( in_ptr, 
                                                   out_ptr, max_output_size, &decoded_bytes);
        /* abort if decompression failed */
        if( 0 != rtn )
        {
            return NULL;
        }

        actual_output_size += (uint64_t) decoded_bytes;
        input_offset += (uint64_t) max_input_size;
    }

    if(!is_compressed)
    {
        memcpy(output_buff, input_buff, input_size);
        actual_output_size = input_size;
        input_offset += input_size;
    }

    assert(actual_output_size == uncompressed_size);
    assert(input_offset == input_size);

    return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, output_buff);
}

adios_datablock * adios_transform_blosc_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(blosc);

#endif
