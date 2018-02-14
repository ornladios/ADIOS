/*
 * adios_transform_lz4_read.c
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

#ifdef LZ4

#include "lz4.h"
#include "adios_transform_lz4_common.h"

int adios_transform_lz4_is_implemented(void)
{
    return 1;
}

/** decompress lz4 compressed data
 * 
 * @param stream [in|out] LZ4 decode stream object
 * @param input_data [in] data to decompress
 * @param input_len size of the input data (in byte)
 * @param output_data [out] buffer for the decompressed output
 * @param max_output_len size of the output buffer (in byte)
 * @param decoded_byte [out] number of decompressed bytes from the input buffer (can be smaller than input_len)
 * @return zero (0) if something goes wrong else non zero
 */
int adios_transform_lz4_decompress(LZ4_streamDecode_t* stream,
                                   const char* input_data,
                                   const adiosLz4Size_t input_len,
                                   char* output_data,
                                   const adiosLz4Size_t max_output_len,
                                   adiosLz4Size_t* decoded_byte)
{
    assert(stream != NULL &&
           input_data != NULL &&
           input_len > 0 &&
           output_data != NULL &&
           max_output_len > 0);


    adiosLz4Size_t result = LZ4_decompress_fast_continue(stream, input_data, output_data, max_output_len);

    *decoded_byte = 0;
    if (result > 0)
    {
        *decoded_byte = result;
    }

    return result <= 0;
}

int adios_transform_lz4_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                  adios_transform_pg_read_request *pg_reqgroup)
{
    void *buf = malloc(pg_reqgroup->raw_var_length);
    assert(buf);
    adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_whole_pg(pg_reqgroup, buf);
    adios_transform_raw_read_request_append(pg_reqgroup, subreq);
    return 0;
}

/* Do nothing for individual subrequest */
adios_datablock * adios_transform_lz4_subrequest_completed(adios_transform_read_request *reqgroup,
                                                           adios_transform_pg_read_request *pg_reqgroup,
                                                           adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}

adios_datablock * adios_transform_lz4_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                            adios_transform_pg_read_request *completed_pg_reqgroup)
{
    uint64_t input_size = (uint64_t) completed_pg_reqgroup->raw_var_length;
    char* input_buff = (char*) (completed_pg_reqgroup->subreqs->data);

    // empty chunk in process group
    if(completed_pg_reqgroup->transform_metadata == NULL)
        return NULL;

    // number of full chunks
    adiosLz4Size_t num_chunks = *((adiosLz4Size_t*) completed_pg_reqgroup->transform_metadata);
    // compressed size of the last chunk
    adiosLz4Size_t compressed_size_last_chunk = *((adiosLz4Size_t*) (completed_pg_reqgroup->transform_metadata + sizeof (adiosLz4Size_t)));

    // flag to mark compressed data (0 means not compressed)
    int is_compressed = 1;
    if (num_chunks == 0 && compressed_size_last_chunk == 0)
        is_compressed = 0;

    /* size of data after compression */
    uint64_t uncompressed_size = adios_get_type_size(reqgroup->transinfo->orig_type, "");
    int d = 0;
    for (d = 0; d < reqgroup->transinfo->orig_ndim; d++)
    {
        uncompressed_size *= (uint64_t) (completed_pg_reqgroup->orig_varblock->count[d]);
    }

    char* output_buff = (char*) malloc(uncompressed_size);
    if (!output_buff)
    {
        log_error("Out of memory allocating %" PRIu64 " bytes in lz4 transform (read)\n", uncompressed_size);
        return NULL;
    }

    LZ4_streamDecode_t lz4StreamDecode_body = {0};
    LZ4_streamDecode_t* lz4StreamDecode = &lz4StreamDecode_body;

    /* offset inside the input buffer */
    uint64_t input_offset = 0;
    /* count decompressed data size */
    uint64_t actual_output_size = 0;

    adiosLz4Size_t chunk = 0;
    for (; (chunk < num_chunks || input_offset < input_size) && is_compressed; ++chunk)
    { 
        const char* in_ptr = input_buff + input_offset;
        char* out_ptr = output_buff + actual_output_size;

        adiosLz4Size_t max_output_size = ADIOS_LZ4_MAX_INPUT_SIZE;
        // handle the last not full chunk
        if (chunk >= num_chunks)
            max_output_size = uncompressed_size - actual_output_size;

        adiosLz4Size_t max_chunk_size = LZ4_compressBound(max_output_size);
        adiosLz4Size_t decoded_bytes = 0;
        int rtn = adios_transform_lz4_decompress(lz4StreamDecode, in_ptr, max_chunk_size,
                                                 out_ptr, max_output_size, &decoded_bytes);
        /* abort if decompression failed */
        if (0 != rtn)
        {
            return NULL;
        }
        actual_output_size += (uint64_t) max_output_size;
        input_offset += (uint64_t) decoded_bytes;
    }

    if (is_compressed == 0)
    {
        memcpy(output_buff, input_buff, input_size);
        actual_output_size = input_size;
        input_offset += input_size;
    }
    
    assert(input_offset == input_size);
    assert(actual_output_size == uncompressed_size);

    return adios_datablock_new_whole_pg(reqgroup, completed_pg_reqgroup, output_buff);
}

adios_datablock * adios_transform_lz4_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(lz4);

#endif
