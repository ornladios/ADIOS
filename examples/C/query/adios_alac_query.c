/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios_read_ext.h"
#include <alacrity.h>


typedef struct {
    // Shared state
    uint64_t partition_capacity;
    uint64_t *meta_offsets;
    uint64_t *data_offsets;
    uint64_t *index_offsets;

    // Read state

    // Write state
    uint64_t cur_meta_offset;
    uint64_t cur_data_offset;
    uint64_t cur_index_offset;
//    _Bool has_last_partition_been_written;
} adios_alac_store_state;


typedef struct {
  _Bool global_meta_loaded;
  ALGlobalMetadata global_meta;
  uint64_t cur_partition;
  adios_alac_store_state *impl_state;
} adios_alac_store;

#define PARTITION_CAPACITY 100;
ALError alacInitStore(adios_alac_store *store, const char *basename, const char *access_mode, _Bool legacyFormat) {

	store->global_meta_loaded = false;

	store->global_meta.num_partitions = 0;
	store->global_meta.total_elements = 0;
	store->global_meta.partition_size = 0;

	store->cur_partition = 0;
    // Set up our implementation-specific parameters
    adios_alac_store_state *params = (adios_alac_store_state*)malloc(sizeof(adios_alac_store_state));
    store->impl_state = params;

    //*******This could be changed **************//
    params->partition_capacity = PARTITION_CAPACITY ;
    params->meta_offsets = malloc(params->partition_capacity * sizeof(uint64_t));
    params->data_offsets = malloc(params->partition_capacity * sizeof(uint64_t));
    params->index_offsets = malloc(params->partition_capacity * sizeof(uint64_t));

    params->cur_meta_offset = 0;
    params->cur_data_offset = 0;
    params->cur_index_offset = 0;

    return ALErrorNone;
}


void alacQueryEngineInit(ALQueryEngine *qe, ALStore *store) {
    qe->store = store;
    qe->gmeta = store->global_meta;
    qe->metadatas = NULL;
}


void alacLoadGlobalMeta(ALQueryEngine *qe, ALStore *store){
//TODO
}

int main (int argc, char ** argv) 
{
    char        filename [256];
    char        varName[256];
    int         rank, size, i, j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;
    uint64_t start[2], count[2], bytes_read = 0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=3");

    strcpy(filename, "adios_alac.bp");
    double lb = 0, hb  =0;
    strcpy(varName, "temp");

    ADIOS_FILE * f = adios_read_open (filename, method, comm, ADIOS_LOCKMODE_NONE, 0);
    if (f == NULL) {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    adios_alac_store alStore;     alacInitStore(&alStore);
    ALGlobalMetadata *gmeta = &alStore.global_meta;
    ALQueryEngine qe;             alacQueryEngineInit(&qe, &alStore);
    adios_alac_store_state *params = (adios_alac_store_state*)alStore.impl_state;
    ALUnivariateQuery query;


    /***** Begin to Load Global Meta DO WE NEED THIS??*********/
    /**** Load global meta from adios file metadata header ***/
    /* questions:
     * 1) how to load the adios metadata?
     * 2) how to write alacrity global metadata as alacrity plugin only writes its partition metadata, index size, and low order bytes into PG
     */
/*
    fread(&gmeta->total_elements, sizeof(uint64_t), 1, params->metadatafp);
    fread(&gmeta->partition_size, sizeof(uint64_t), 1, params->metadatafp);
    fread(&gmeta->num_partitions, sizeof(uint64_t), 1, params->metadatafp);

    if (gmeta->num_partitions + 1 > params->partition_capacity)
        extendOffsetBuffers(params, gmeta->num_partitions + 1, false);

    fread(params->meta_offsets, sizeof(uint64_t), gmeta->num_partitions + 1, params->metadatafp);
    fread(params->data_offsets, sizeof(uint64_t), gmeta->num_partitions + 1, params->metadatafp);
    fread(params->index_offsets, sizeof(uint64_t), gmeta->num_partitions + 1, params->metadatafp);

    if (!ferror(params->metadatafp)) {
        store->global_meta_loaded = true;
        return ALErrorNone;
    } else {
        return ALErrorSomething; // TODO: Error code
    }*/

    /***** End of Loading Global Meta *********/

    ALQueryEngineStartUnivariateDoubleQuery(&qe, lb, hb, VALUE_RETRIEVAL_QUERY_TYPE, &query);

    /*********** doQuery ***************
     *
     * 1. Open partition  [locate offsets of meta, data, and index for the partition]
     * 2. Read Partition Meta from file => meta
       3. find touched bins:  are_bins_touched = findBinRange1C(meta, uniquery, &start_bin, &end_bin);
       4. read index of touched bins: ALPartitionStoreReadIndexBins
       5. read dataBin
     */

    adios_read_set_transforms_enabled(f, 0); // load transformed 1-D byte array
    ADIOS_VARINFO * v = adios_inq_var (f, varName);
    ADIOS_TRANSINFO *ti = adios_inq_transinfo(f, v);
    adios_inq_trans_metadata(f, v, ti);
    // ensure its alacrity plugin & get metaData, index and data Size
    uint64_t metaSize = 0, indexSize = 0,  dataSize =0;
    if (ti->transform_type == adios_transform_alacrity && ti->transform_metadata_len == 3 * sizeof(uint64_t) ) {
    	metaSize = ((uint64_t * ) (ti->transform_metadata)) [0];
    	indexSize = ((uint64_t * ) (ti->transform_metadata)) [1];
    	dataSize = ((uint64_t * ) (ti->transform_metadata)) [2];
    }else {
    	// some thing is wrong
    	printf(".....\n");
    	adios_read_close (f);
		MPI_Barrier (comm);
		adios_read_finalize_method (method); //TODO: verify that do we need this?
		MPI_Finalize ();
    	return ;
    }

    // transformed data has 1 dimension,
    // 1. load partition Metadata
    // TODO: how to ensure the element type is byte?
    int ndim = 1;
    int startStep = 0, numStep = 1;
    sel = adios_selection_boundingbox (ndim,  0, &metaSize);
    ALMetadata partitionMeta ;
    memstream_t ms = memstreamInitReturn(malloc(metaSize));
    adios_schedule_read (f, sel, varName, startStep, numStep, ms.buf);
    adios_perform_reads (f, 1);
    ALDeserializeMetadata(&partitionMeta, &ms);
    memstreamDestroy(&ms, true);

    //2. find touched bin
    bin_id_t low_bin , hi_bin;
    _Bool are_bins_touched = findBinRange1C(partitionMeta, query, &low_bin, &hi_bin);

    //3. load index size
    uint64_t indexOffset = metaSize ; // metadata size is index start offset

    const char insigbytes = insigBytesCeil(partitionMeta);
    const uint64_t first_bin_off = indexOffset +  ALGetIndexBinOffset(partitionMeta, low_bin);
    const uint64_t last_bin_off = indexOffset + ALGetIndexBinOffset(partitionMeta, hi_bin);
    const uint64_t bin_read_len = last_bin_off - first_bin_off;
    sel = adios_selection_boundingbox (ndim,  &first_bin_off, &bin_read_len);
    void *indexData = NULL;
    adios_schedule_read (f, sel, varName, startStep, numStep, indexData);
    adios_perform_reads (f, 1);
    ALIndex index = (ALIndex) indexData;

    ALUnivariateQueryResult result;
    const ALBinLayout * const bl = &(partitionMeta.binLayout);
    uint64_t resultCount;
    ALData resultData;

    resultCount = bl->binStartOffsets[hi_bin] - bl->binStartOffsets[low_bin];

    if (query->queryType == REGION_RETRIEVAL_INDEX_ONLY_QUERY_TYPE) {
        		// No data, no candidate checks, just return the results as-is
        		data = NULL;
            	populateQueryResult(*result, resultData, index, resultCount);
	} else {
		// TODO
	}

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);
    MPI_Finalize ();
    return 0;
}
