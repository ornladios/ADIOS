/*
 * test_query_alac_internal.c
 *
 *  Created on: Jun 22, 2014
 *      Author: xczou
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
//#include "mpi.h"
#include "adios_read.h"
#include "adios_read_ext.h"
#include "./query/query_alac.c"

void printALMetadata(ALMetadata *pm){
	printf("global offset id %" PRIu64" \n", pm->globalOffset);
	printf("partition length %" PRIu32" \n" , pm->partitionLength);
	printf("significant bits %d \n" , pm->significantBits);
	printf("num bins %" PRIu32" \n", pm->binLayout.numBins);
	printf("index form %d \n", pm->indexMeta.indexForm);

}

void test_adios_alac_meta_read(ADIOS_FILE *fp, ADIOS_VARINFO *v){
	int i = adios_inq_var_blockinfo(fp, v);
	if (i != 0){
		printf("error from adios_inq_var_blockinfo \n");
		return ;
	}
	ADIOS_VARTRANSFORM *tfv = adios_inq_var_transform(fp, v);
	ADIOS_TRANSFORM_METADATA * tmetas = tfv->transform_metadatas;
	ADIOS_VARBLOCK * bi = v->blockinfo;
	ADIOS_SELECTION * sel;
	// TODO: bi includes all block ids spanning mutliple timesteps?
	uint64_t metaSize;
	int startStep =0 /*TODO*/, numStep =1;
	for(i=0; i< v->sum_nblocks; i++){
		ADIOS_TRANSFORM_METADATA tmeta = tmetas[i]; // i = blockId
		uint64_t * threeData = (uint64_t *) tmeta.content;
		assert(tmeta.length == 24);
		metaSize = threeData[0];
		ALMetadata metadata;
		readPartitionMeta(i, metaSize,fp, v
							,startStep,numStep,&metadata);

		printf("partitionMeta in block [%d]\n", i);
		printf("meta size: %" PRIu64 ", index size: %" PRIu64
							", data size: %" PRIu64 "\n"
							, threeData[0], threeData[1], threeData[2]);



		printALMetadata(&metadata);
//		ALPrintMetaData(metadata)
		free(metadata.binLayout.binStartOffsets);
		free(metadata.binLayout.binValues);
		if (metadata.indexMeta.indexForm == ALCompressedInvertedIndex)
			free(metadata.indexMeta.u.ciim.indexBinStartOffsets);
//		printPartitionMeta(&partitionMeta);
	}
}


/*
 * RUN with one processor is enough
 * sample usage:
 * ./test_query_alac_internal ../*.bp rdm
 */
int main (int argc, char ** argv)
{

    char        filename [256];
    int         i, j, datasize, if_any;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel1, * sel2;
    ADIOS_VARCHUNK * chunk = 0;
    double * data = NULL;
    uint64_t start[2], count[2], npoints, * points;

    MPI_Init (&argc, &argv);

    if (argc < 3 ){
    	printf(" usage: %s {input bp file} {variable name} \n", argv[0]);
    	return 1;
    }

    adios_read_init_method (method, comm, NULL);

    ADIOS_FILE * fp = adios_read_open_file (argv[1], method, comm);

    if ( fp == NULL){
    	printf(" can not open file %s \n", argv[1]);
    	return 1;
    }

    char varName[256];
    strcpy(varName, argv[2]);

    ADIOS_VARINFO* v = adios_inq_var(fp, varName);

    data_view_t  dv = PHYSICAL_DATA_VIEW;
    adios_read_set_data_view(fp, dv);

    //====================start to test ==================//
    test_adios_alac_meta_read(fp, v);

    adios_read_close (fp);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);


    MPI_Finalize ();
    return 0;
}

