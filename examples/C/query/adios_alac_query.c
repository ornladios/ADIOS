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
    ALQueryEngine qe;
    ALUnivariateQuery query;
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

    int totalPG =  v->sum_nblocks;

    uint64_t * metaSizes = (uint64_t *) malloc(sizeof(uint64_t) * totalPG);
    uint64_t * indexSizes = (uint64_t *) malloc(sizeof(uint64_t) * totalPG);
    uint64_t * dataSizes = (uint64_t *) malloc(sizeof(uint64_t) * totalPG);
    // Will assume user providing PG range for us
    // Right now, use all PGs
    for (int i = 0; i < totalPG ; i ++){
    	ADIOS_TRANSINFO_TRANSMETA * tmeta = ti->transform_metadatas[i];
    	assert(tmeta->length == 24);
    	uint64_t * threeData = (uint64_t *) tmeta->content;
    	metaSizes[i] = threeData[0];
    	indexSizes[i] = threeData[1];
    	dataSizes[i] = threeData[2];

    	uint64_t metaSize = metaSizes[i];

        // transformed data has 1 dimension,
        // 1. load partition Metadata
        int ndim = 1;
        int startStep = 0, numStep = 1; // TODO: what time step we will have?
        sel = adios_selection_writeblock_bounded(i, 0,metaSize);
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
        sel = adios_selection_writeblock_bounded(i, metaSize,bin_read_len);
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

    }

    free(metaSizes);
    free(indexSizes);
    free(dataSizes);


    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);
    MPI_Finalize ();
    return 0;
}
