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

void printALMetadata(ALMetadata *pm);
void printBB(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *olBB);


void printALMetadata(ALMetadata *pm){
	printf("global offset id %" PRIu64" \n", pm->globalOffset);
	printf("partition length %" PRIu32" \n" , pm->partitionLength);
	printf("significant bits %d \n" , pm->significantBits);
	printf("num bins %" PRIu32" \n", pm->binLayout.numBins);
        int i;
        printf("bin values:\n");
        for (i = 0; i < pm->binLayout.numBins; i++) {
	    printf("%x \t", pm->binLayout.binValues[i]);
        }
        printf("\n");
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
		free(metadata.binLayout.binStartOffsets);
		free(metadata.binLayout.binValues);
		if (metadata.indexMeta.indexForm == ALCompressedInvertedIndex)
			free(metadata.indexMeta.u.ciim.indexBinStartOffsets);
	}
}


void printBB(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *olBB)
{
	int i ;
    for(i=0; i < olBB->ndim; i++){
			printf("[%" PRIu64 ": %"PRIu64"], " , olBB->start[i], olBB->start[i] + olBB->count[i] -1 );
	}
    printf("\n");
}

void checkRidConversion(rid_t rid, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *pgBB, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *userBB, int Corder)
{
    rid_t newRid = 0;
    if (ridConversionWithCheck(rid, pgBB->start, pgBB->count
				, userBB->start, userBB->count, userBB->ndim, &newRid, Corder)) {
			printf("rid %"PRIu32" is converted to %" PRIu32 " \n"
					,rid, newRid);
		}else{
			printf("rid is not in user selection \n");
		}
}

void testContainingAndRidConversion(ADIOS_FILE *fp, ADIOS_VARINFO *v){

	//******************* Get Selection BoundingBox first ***************//
	int from_step =0;
	int nsteps =1, i = 0,j =0;
	/*uint64_t start[3] = {45, 10,0};
	uint64_t count[3] = {65, 40, 60};*/
	//********this is boundary case, in which user selection is the entire domain *******//
	uint64_t start[3] = {0, 0,0};
	uint64_t count[3] = {128, 64, 64};
	ADIOS_SELECTION *userSelection =  adios_selection_boundingbox(3, start, count);
	ADIOS_PG_INTERSECTIONS* intersectedPGs2 = adios_find_intersecting_pgs(
				fp, v->varid, userSelection, from_step, nsteps);
	int totalnpg = intersectedPGs2->npg;
	ADIOS_PG_INTERSECTION *  PGs = intersectedPGs2->intersections;

	printf("user's input selection box: ");
	assert(userSelection->type == ADIOS_SELECTION_BOUNDINGBOX );
	const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *userBB = &(userSelection->u.bb);
    printBB(userBB);

    ADIOS_PG_INTERSECTION pg;
	for(j= 0; j < totalnpg; j ++){
		pg = PGs[j];
		printf("PG[%d], timestep[%d], PG id in TS[%d]\n", pg.blockidx, pg.timestep, pg.blockidx_in_timestep);
		const ADIOS_SELECTION * pgSelBox = pg.pg_bounds_sel;
		assert (pgSelBox->type == ADIOS_SELECTION_BOUNDINGBOX );
		printf("touched PG bounding box : ");
		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *pgBB = &(pgSelBox->u.bb);
		printBB(pgBB);


		ADIOS_SELECTION * intersectedBox = pg.intersection_sel;
		assert (intersectedBox->type == ADIOS_SELECTION_BOUNDINGBOX );
		printf("overlapped/intersected bounding box: ");
		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *interBB = &(intersectedBox->u.bb);
		printBB(interBB);

		printf("does user input bounding box contain the PG box: %d \n", boxEqual(pgBB, interBB));

		rid_t rid = 0;
		int Corder =1;
		checkRidConversion(rid, pgBB, userBB, Corder);
		rid = 13;
		checkRidConversion(rid, pgBB, userBB, Corder);
	}

	adios_selection_delete(userSelection);
	adios_free_pg_intersections(&intersectedPGs2);

}

ADIOS_ALAC_BITMAP * mockBitmap(uint64_t length, uint64_t numSetBits, uint64_t realElmSize){
	ADIOS_ALAC_BITMAP * b = (ADIOS_ALAC_BITMAP *) malloc(sizeof(ADIOS_ALAC_BITMAP));
	b->length = length;
	b->bits = (uint64_t *) calloc(b->length, sizeof(uint64_t));
	b->numSetBits = numSetBits;
	b->realElmSize = realElmSize;
	initLastConvRid(b);
	return b;
}
void  testBitsCal(){
	ADIOS_ALAC_BITMAP *b = mockBitmap(10, 0, 630) ;
	uint64_t tbits = calSetBitsNum(b);
	assert(tbits == 0);
	b->bits[0] = 3;
	tbits = calSetBitsNum(b);
	assert(tbits == 2);
	FreeALACBITMAP(b);
}

void testResolveBoundary(){
	ADIOS_QUERY q;
	q.predicateValue = "0.651";
	q.predicateOp = ADIOS_LT;
	double lb , hb;
	resolveQueryBoundary(&q, &hb, &lb);
	assert(hb == 0.651);
	assert(lb == DBL_MIN);
}

void testAlacBitmapConversion(){
	ADIOS_ALAC_BITMAP *b = mockBitmap(1024, 0, 65536) ;
	void * mem = NULL;
	mem = convertALACBitmapTomemstream(b);
	ADIOS_ALAC_BITMAP rb ;
	convertMemstreamToALACBitmap(mem, &rb);
	assert(rb.length == b->length);
	assert(rb.numSetBits == b->numSetBits);
	assert(rb.realElmSize == b->realElmSize );
	assert(rb.lastConvRid == b->lastConvRid);
	assert(memcmp(rb.bits, b->bits, sizeof(uint64_t)* (b->length)) == 0);
	FREE(mem);
	FreeALACBITMAP(b);
	FREE(rb.bits);
}


void test_adios_read_alac_index(ADIOS_FILE *fp, ADIOS_VARINFO *v){
    int blockId = 0;
    ADIOS_VARTRANSFORM *tfv = adios_inq_var_transform(fp, v);
    ADIOS_TRANSFORM_METADATA *tmetas = tfv->transform_metadatas;
    ADIOS_TRANSFORM_METADATA tmeta = tmetas[blockId];
    assert(tmeta.length == 24);
    uint64_t *threeData = (uint64_t*)tmeta.content;
    printf("meta size: %" PRIu64 ", index size: %" PRIu64
			", data size: %" PRIu64 "\n"
			, threeData[0], threeData[1], threeData[2]);

    uint64_t metaSize = threeData[0];
    uint64_t indexSize = threeData[1];
    uint64_t dataSize  = threeData[2];
    int startStep = 0, numStep = 1;

    ALMetadata partitionMeta;
	const ALBinLayout * const bl = &(partitionMeta.binLayout);
    readPartitionMeta(blockId, metaSize,fp, v ,startStep,numStep,&partitionMeta);

    bin_id_t low_bin =0 , hi_bin = 1;
    uint64_t indexStartPos = metaSize + dataSize; // meta + data + index

    const uint64_t first_bin_off = ALGetIndexBinOffset( &partitionMeta, low_bin);
	const uint64_t last_bin_off = ALGetIndexBinOffset( &partitionMeta, hi_bin);
	const uint64_t bin_read_len = last_bin_off - first_bin_off;

    void *indexData =  (void *) malloc(sizeof(char)*bin_read_len);

    ADIOS_SELECTION *sel = adios_selection_writeblock_bounded(blockId, indexStartPos, bin_read_len, 0);
	adios_schedule_read_byid(fp, sel, v->varid, startStep, numStep, indexData);
	adios_perform_reads(fp, 1);
	adios_selection_delete(sel);

	ALIndex index = (ALIndex) indexData;
	rid_t * idx = (rid_t *) index;
	uint64_t resultCount = bl->binStartOffsets[hi_bin] - bl->binStartOffsets[low_bin];
	uint64_t ni;
	printf("index from bin %" PRIu32 " to bin %"PRIu32 " [%" PRIu64 "]\n "
			, low_bin, hi_bin, resultCount);
	for ( ni = 0; ni < resultCount; ni++) {
		rid_t rid_val = idx[ni];
		assert(rid_val <=65535);
		assert(rid_val >= 0 );

//		printf("%"PRIu32 ",", rid_val);
	}
//	printf("\n");

	FREE(index);
	adios_free_var_transform(tfv);
}


void test_adios_query_alac_retrieval_points3d(){
	create_lookup(set_bit_count, set_bit_position);
	ADIOS_ALAC_BITMAP *b = mockBitmap(10, 0, 630) ;
	b->bits[0] = 6;     // 0110 .... | 0 ->        rid: 1, 2
	b->bits[1] = 7;     // 111 ....  |  64 ->      rid: 64 65 66
	b->bits[2] = 8;     // 0001 .... |  128 ->     rid: 131
	b->bits[3] = 255;   // 1111 1111 ...  | 192 -> rid: 192 193 194 195 196 197 198 199
	b->bits[4] = 100;   // 0010 0110 ...  | 256 -> rid: 258 256 257
	b->bits[5] = 130;   // 0100 0001 ...  | 320 -> rid: 321 328

	uint64_t retrieval_size =9, i = 0;
	uint64_t start1[] = {0, 0, 0};
	uint64_t count1[] = {64, 32,32};
	int ndim = 3;
	int Corder =1;
	ADIOS_SELECTION* box1 = adios_selection_boundingbox(ndim, start1, count1);
	int retrieveTimes =2, p = 0;
	switch (box1->type) {
		case ADIOS_SELECTION_BOUNDINGBOX: {
			ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(box1->u.bb);

			for (p = 0; p < retrieveTimes ; p ++ ){
				uint64_t dataSize = retrieval_size * (bb->ndim);
				uint64_t* points = (uint64_t*) (malloc(dataSize * sizeof(uint64_t)));
				adios_query_alac_retrieval_pointsNd( b,  retrieval_size, bb , points /*OUT*/ , Corder);
				printf("[%d] retrieved points: ", p);
				for(i = 0; i < retrieval_size ; i ++){
					printf("[%"PRIu64", %"PRIu64",%"PRIu64"] , ", points[i*3], points[i*3+1], points[i*3+2]);
				}
				printf("\n");
				FREE(points);
			}
		}
	}

	adios_selection_delete(box1);
	FreeALACBITMAP(b);
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
    printf("//====================test_adios_alac_meta_read ==================//\n");
    test_adios_alac_meta_read(fp, v);
    printf("///////////////////////////////////////\n");

    printf("//====================testContainingAndRidConversion==================//\n");
//    testContainingAndRidConversion(fp, v);
    printf("///////////////////////////////////////\n");

    printf("//====================testAlacBitmapConversion==================//\n");
    testAlacBitmapConversion();
    printf("// test is passed \n");
    printf("///////////////////////////////////////\n");


    printf("//====================testBitsCal==================//\n");
    testBitsCal();
    printf("// test is passed \n");
    printf("///////////////////////////////////////\n");

    printf("//====================testResolveBoundary==================//\n");
    testResolveBoundary();
    printf("// test is passed \n");
    printf("///////////////////////////////////////\n");


    printf("//====================test_adios_read_alac_index==================//\n");
    test_adios_read_alac_index(fp,v);
    printf("// test is passed \n");
    printf("///////////////////////////////////////\n");


    printf("//====================test_adios_query_alac_retrieval_points3d==================//\n");
    test_adios_query_alac_retrieval_points3d();
    printf("// test is passed \n");
    printf("///////////////////////////////////////\n");



    adios_read_close (fp);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);


    MPI_Finalize ();
    return 0;
}

