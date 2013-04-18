#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "util.h"
#include "adios_transforms_hooks_read.h"
#include "adios_transforms_reqgroup.h"

#ifdef APLOD

#include <stdint.h>

#include "aplod.h"
#include "adios_transform_identity_read.h"

#define MAX_COMPONENTS 16

typedef struct {
    int numComponents;
    int32_t components[MAX_COMPONENTS];
} aplod_meta_t;

void parse_aplod_meta(char *transform_metadata, aplod_meta_t *metaout) {
    transform_metadata += sizeof (uint64_t);

    metaout->numComponents = *(int8_t*)transform_metadata;
    transform_metadata += sizeof(int8_t);

    memcpy(metaout->components, transform_metadata, metaout->numComponents * sizeof(int32_t));
    transform_metadata += metaout->numComponents * sizeof(int32_t);
}

typedef struct {
    uint64_t numElements;
    uint64_t startOff;
    char *outputBuf;
} aplod_read_meta_t;

int adios_transform_aplod_generate_read_subrequests(adios_transform_read_request *reqgroup,
                                                    adios_transform_pg_read_request *pg_reqgroup)
{
    uint64_t start_off, end_off; // Start/end byte offsets to read between
    compute_sieving_offsets_for_pg_selection(pg_reqgroup->pg_intersection_sel, &pg_reqgroup->pg_bounds_sel->u.bb, &start_off, &end_off);

    aplod_meta_t aplodmeta;
    parse_aplod_meta(reqgroup->transinfo->transform_metadata, &aplodmeta);

    int typelen = adios_get_type_size(reqgroup->transinfo->orig_type, NULL);
    uint64_t elemcount = end_off - start_off;
    char *buf = malloc(elemcount * typelen);

    int numByteColsDone = 0;
    int i;

    // Read the element start_pos..end_pos sieving segment from each column
    for (i = 0; i < aplodmeta.numComponents; i++) {
        adios_transform_raw_read_request *subreq = adios_transform_raw_read_request_new_byte_segment(
                pg_reqgroup,
                numByteColsDone * elemcount + start_off * aplodmeta.components[i],
                numByteColsDone * elemcount + end_off * aplodmeta.components[i],
                buf + elemcount * numByteColsDone);
        adios_transform_raw_read_request_append(pg_reqgroup, subreq);
        numByteColsDone += aplodmeta.components[i];
    }

    aplod_read_meta_t *arm = (aplod_read_meta_t*)malloc(sizeof(aplod_read_meta_t));
    *arm = (aplod_read_meta_t){ .numElements = elemcount, .outputBuf = buf, .startOff = start_off };
    pg_reqgroup->transform_internal = arm; // Store it here to be safe, since each of the subreqs has a different piece of it

    return 0;
}

// Do nothing for individual subrequest
adios_datablock * adios_transform_aplod_subrequest_completed(adios_transform_read_request *reqgroup,
                                                            adios_transform_pg_read_request *pg_reqgroup,
                                                            adios_transform_raw_read_request *completed_subreq)
{
    return NULL;
}



adios_datablock * adios_transform_aplod_pg_reqgroup_completed(adios_transform_read_request *reqgroup,
                                                             adios_transform_pg_read_request *completed_pg_reqgroup)
{
    uint32_t elementSize = adios_get_type_size(reqgroup->transinfo->orig_type, "");

    aplod_read_meta_t *arm = (aplod_read_meta_t *)completed_pg_reqgroup->transform_internal;
    void *compressed_buff = arm->outputBuf;

    uint32_t numElements = arm->numElements;
    uint64_t decompressed_len = numElements * elementSize;
    void* decompressed_buff = malloc (decompressed_len);

    int8_t numComponents = 0;
    int32_t *componentVector = 0;

    aplod_meta_t aplodmeta;
    parse_aplod_meta(reqgroup->transinfo->transform_metadata, &aplodmeta);

    APLODConfig_t *config = APLODConfigure (aplodmeta.components, aplodmeta.numComponents);
    config->blockLengthElts = numElements;

    APLODReconstructComponents  (config,
                                    numElements,
                                    0,
                                    aplodmeta.numComponents,
                                    0,
                                    0,
                                    decompressed_buff,
                                    compressed_buff
                                );

    free (config->byteVector);
    free (config->byteVectorPS);
    free (config);

    // Clear the buffer pointers for all raw read requests, because they all point
    // to the same buffer, and would be free'd by the framework if we didn't clear here
    adios_transform_raw_read_request *rrr = completed_pg_reqgroup->subreqs;
    for (; rrr; rrr = rrr->next)
        rrr->data = NULL;

    // Free the actual buffer, and the special metadata container we added
    free (arm->outputBuf);
    free (arm);

    return adios_datablock_new_ragged_offset(reqgroup->transinfo->orig_type,
                                             completed_pg_reqgroup->timestep,
                                             completed_pg_reqgroup->pg_bounds_sel,
                                             arm->startOff,
                                             decompressed_buff);
}

adios_datablock * adios_transform_aplod_reqgroup_completed(adios_transform_read_request *completed_reqgroup)
{
    return NULL;
}


#else

DECLARE_TRANSFORM_READ_METHOD_UNIMPL(aplod);

#endif

