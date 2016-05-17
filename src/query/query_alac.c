/*
 * query_alacrity.c
 *
 *  Created on: Jun 1, 2014
 *      Author: xczou
 */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "public/adios_error.h"
#include "public/adios_read_ext.h"
#include "public/adios_query.h"
#include "public/adios_selection.h"
#include "transforms/adios_transform_alacrity_common.h"
#include "core/a2sel.h"
#include "core/common_read.h"
#include "common_query.h"
#include "query_utils.h"
#include <alacrity.h>


#ifdef ALACRITY
#include "alacrity.h"

/************Uncompressed bitmap*********/
/***********This is an internal data structure for multi-variate constraints query processing ******/
typedef struct{
	uint64_t *bits;       // uint64_t array holds the bits
	uint64_t length;      // the bits array size
	uint64_t numSetBits;  // actual element size of this bits array represents ( the number of set bits(1s) in the bits array)
	uint64_t realElmSize; // the number of elements the bits are presenting
	                      // if number of elements is 100, then, the length = 2 ( 100 / 64 + 100%64 ), and realElmSize = 100
	uint64_t lastConvRid; // it is only used when the bit map is converted to RIDs, it indicates the bit position (local RID) of last RID conversion
	                      // if the last converting bit array is ....(234) 0001 0001 0000 0000 , then lastConvRid = 238 indicates the left RID 242(0001) has
	                      // not been converted
} ADIOS_ALAC_BITMAP;

#define FreeALACBITMAP(b) { FREE(b->bits); FREE(b) }


#define MAX_DIMS 32

#define PFORDELTA_CHUNK_SIZE 128 // the chunk size of PFORDELTA 

// This three variables are used for timing the transformer layer,
// the reason of being global variables is because calling transformer layer APIs scatter
// I have to sum every place that calls transformer layer
double ckSpatialCheck = 0.0;
double ckSetBitMap = 0.0; // timing is in the macro in which the pair of "#ifdef" and "#endif" can not be used
#ifdef BREAKDOWN
	double preparationTime = 0.0, preparationStart = 0;
	double metaTotal = 0.0, metaStart=0.0;// timing metadata read
	double idxTotal = 0.0, idxStart =0.0;  // timing index read
	uint64_t idxReadSize = 0, dataReadSize =0; // in KB
	double dataTotal = 0.0, dataStart =0.0; // timing low-order byte read
	double procTotal = 0.0, procStart = 0.0; // timing for porc_write_block
	double findPGTotal = 0.0, findPGStart = 0.0; // timing for porc_write_block
	double candidateCheckTotal = 0.0, candidateCheckStart = 0.0;
	double ckRids = 0.0, ckReconstituteStart= 0.0, ckReconstitute=0.0;
	uint64_t ckTotalElm = 0;
	double decodeTotal = 0.0, wholeDecodeTotal=0.0, decodeStart = 0.0, pgCTotal = 0.0, pgPartialTotal = 0.0; // timing for porc_write_block
    double setRidTotal = 0.0, setRidStart = 0.0;
    double alacPartitionMetaTotal=0.0, alacPartitionMetaStart = 0.0;
    int numTouchedPGs = 0;
    bin_id_t numTouchedBins =0;

#endif

/**** Funcs. that are internal funcs. ********/

uint64_t * convertALACBitmapTomemstream( ADIOS_ALAC_BITMAP * b);

void convertMemstreamToALACBitmap( void *mem , ADIOS_ALAC_BITMAP * bout /*OUT*/);

ADIOS_ALAC_BITMAP * adios_alac_process(ADIOS_QUERY* q, int timeStep,
		bool estimate);

ADIOS_ALAC_BITMAP * adios_alac_bitsOp(ADIOS_ALAC_BITMAP * op1,
		ADIOS_ALAC_BITMAP * op2, enum ADIOS_CLAUSE_OP_MODE operator);

uint64_t calSetBitsNum(ADIOS_ALAC_BITMAP *b);

// since we could not set initial lastConvRid to 0 ( 0 represents the first RID in the bitmap )
// we need extra effort to set lastConvRid to initial state
void initLastConvRid (ADIOS_ALAC_BITMAP *b);

bool isLastConvRidInit(const ADIOS_ALAC_BITMAP *b);

static inline int coordinateConversionWithCheck(uint64_t * coordinates, const  int dim
		, const  uint64_t *srcstart, const  uint64_t *deststart, const  uint64_t *destend);

static inline bool ridConversionWithCheck(const rid_t rid/*relative to local src selectoin*/
		, const uint64_t *srcstart, const uint64_t *srccount, const uint64_t *deststart, const uint64_t *destcount,
		const int dim, rid_t *relativeRid , const int Corder);

static rid_t ridConversionWithoutCheck(const rid_t rid/*relative to local src selectoin*/,
		const uint64_t *srcstart, const uint64_t *srccount, const uint64_t *deststart, const uint64_t *destcount,
		const int dim, const int Corder);

static inline void ridToCoordinates(const int dim, const int Corder, const rid_t rid, const uint64_t * const dimSize , uint64_t * const coordinates /*OUT*/);

static void init_lookup();

bool boxEqual(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *pgBB, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *interBB);

void readPartitionMeta(int blockId, uint64_t metaSize, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, ALMetadata *pm /*OUT*/);

void readTransformedElms(ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, int blockId, uint64_t start_elem, uint64_t num_elems, int is_timestep_relative, void * outputData/*out*/);

void readIndexData(int blockId, uint64_t offsetSize /*in bytes*/
		,uint64_t length /*in bytes*/, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, void *idxBytes /*OUT*/);

void readLowOrderBytes(int blockId, uint64_t offsetSize /*in bytes*/
		,uint64_t length /*in bytes*/, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, void *idxBytes /*OUT*/);

void resolveQueryBoundary(ADIOS_QUERY *adiosQuery, double *hb, double *lb);

static inline void setRidToBits(bool isPGCovered ,uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int ndim, rid_t * idx, uint64_t totalRids,
		ADIOS_ALAC_BITMAP *alacResultBitmap /*IN&OUT*/, int Corder);


void proc_write_block(int blockId, bool isPGCovered, ADIOS_VARTRANSFORM *ti, ADIOS_QUERY * adiosQuery, int startStep, bool estimate
		, ALUnivariateQuery * alacQuery , double lb , double hb
		,uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount
		, ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/ , int Corder);

void setBitsinBitMap(rid_t rid, ADIOS_ALAC_BITMAP * alacResultBitmap);

void setBitsinBitMap(rid_t rid, ADIOS_ALAC_BITMAP * alacResultBitmap){
	uint32_t word = (uint32_t) (rid >> 6);
	if (word > alacResultBitmap->length){
		printf("what a hell\n");
	}
	assert(word <= alacResultBitmap->length);
	alacResultBitmap->bits[word]
			|= (1LL << (rid & 0x3F));
}

#define BITNSLOTS64(nb) ((nb + 64 - 1) / 64)

static uint8_t bits_in_char[256] = {
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
		B6(0), B6(1), B6(1), B6(2)};

static int isInitialized = 0; // 0: the lookup table is not initialized ; 1 : the lookup table is initialized
static unsigned char set_bit_count[65536];
static unsigned char set_bit_position[65536][16];


/**** END -- Funcs. that are internal funcs. ********/

void initLastConvRid (ADIOS_ALAC_BITMAP *b){
	b->lastConvRid = b->realElmSize + 1;
}

bool isLastConvRidInit(const ADIOS_ALAC_BITMAP *b){
	return (b->lastConvRid == (b->realElmSize +1));
}




/*since we are set data view is physical view, the element size is 1 byte */
//this is a basic read function for other read functions , such as readIndex, readLowOrderBytes, readMeta.
void readTransformedElms(ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, int blockId, uint64_t start_elem, uint64_t num_elems, int is_timestep_relative, void * outputData/*out*/){
	ADIOS_SELECTION *sel = adios_selection_writeblock_bounded(blockId, start_elem, num_elems, is_timestep_relative);
	common_read_schedule_read_byid(fp, sel, vi->varid, startStep, numStep, NULL, outputData);
	common_read_perform_reads(fp, 1);
	// adios_selection_writeblock_bounded internally malloc data for adios_selection
	// so I need to free it before the next usage
	a2sel_free(sel);
}

void readBlockData(int gBlockId /*global block id */, ADIOS_QUERY * adiosQuery, int startStep,
		ADIOS_VARINFO * varInfo, uint64_t dataElmNum, void ** data){
	adios_read_set_data_view(adiosQuery->file, LOGICAL_DATA_VIEW); // switch to the transform view,
	int dataElmSize = common_read_type_size(varInfo->type, NULL); // data element size, in bytes
	char * blockData = (char*) (*data);
	blockData = (char *) malloc(sizeof(char) * dataElmSize * dataElmNum);
	ADIOS_SELECTION *sel = adios_selection_writeblock_bounded(gBlockId, 0, dataElmNum, 0); // entire PG selection
	common_read_schedule_read_byid(adiosQuery->file, sel, varInfo->varid, startStep, 1, NULL, blockData);
}
void readIndexData(int blockId, uint64_t offsetSize /*in bytes*/
		,uint64_t length /*in bytes*/, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, void *idxBytes /*OUT*/){
	const int is_streaming = fp->is_streaming;
	readTransformedElms(fp, vi, startStep, numStep, blockId, offsetSize, length, is_streaming, idxBytes);
}

void readLowOrderBytes(int blockId/*if file, global block id; if streaming, block id relative to current step*/, uint64_t offsetSize /*in bytes*/
		,uint64_t length /*in bytes*/, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, void *idxBytes /*OUT*/){
	const int is_streaming = fp->is_streaming;
	readTransformedElms(fp, vi, startStep, numStep, blockId, offsetSize,length, is_streaming, idxBytes);
}

void readPartitionMeta( int blockId /*if file, global block id; if streaming, block id relative to current step*/, uint64_t metaSize, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, ALMetadata *pm /*OUT*/){
	const int is_streaming = fp->is_streaming;
	memstream_t ms = memstreamInitReturn(malloc(metaSize));
	uint64_t metaStartPos= 0;
	readTransformedElms(fp, vi, startStep,numStep,blockId, metaStartPos, metaSize, is_streaming, ms.buf);
	ALDeserializeMetadata(pm, &ms);
	memstreamDestroy(&ms, true);

}

bool boxEqual(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *pgBB, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *interBB){
	if (pgBB->ndim != interBB->ndim)
		return false;
	 // if these two boxes are not equal, then,
	// it means two boxes are intersecting / partial overlapping
	int k = 0;
	for(k = 0 ; k < pgBB->ndim ; k ++){
		if (pgBB->count[k]!=interBB->count[k] ||
				pgBB->start[k] != interBB->start[k]){
			return false;
		}
	}
	return true;
}


static void init_lookup() {
	memset(set_bit_count, 0, 256);
	int i = 0, j;
	for (i = 0; i < 65536; i++) {
//		set_bit_count[i] = __builtin_popcount(i); // total bit 1 for value
		set_bit_count[i] =  bits_in_char [i & 0xff]
						   +  bits_in_char [(i >>  8) & 0xff]
						   +  bits_in_char [(i >> 16) & 0xff]
						   +  bits_in_char [(i >> 24) & 0xff]
						   ;
		unsigned short int temp = i;
		int counter = set_bit_count[i] - 1;
		for (j = 15; j >= 0; j--) {
			unsigned int temp1 = temp >> j & 0x0001;
			if (temp1 == 1) {
				set_bit_position[i][counter--] = j;
			}
		}

	}
}




uint64_t calSetBitsNum(ADIOS_ALAC_BITMAP *b) {

	uint64_t total = 0, count = 0, i = 0;
	for (; i < b->length; i++) {
		  count = bits_in_char[b->bits[i] & 0xff]
				+ bits_in_char[(b->bits[i] >> 8) & 0xff]
				+ bits_in_char[(b->bits[i] >> 16) & 0xff]
				+ bits_in_char[(b->bits[i] >> 24) & 0xff]
				+ bits_in_char[(b->bits[i] >> 32) & 0xff]
				+ bits_in_char[(b->bits[i] >> 40) & 0xff]
				+ bits_in_char[(b->bits[i] >> 48) & 0xff]
				+ bits_in_char[(b->bits[i] >> 56) & 0xff];
		total += count;
	}
	return total;
}

// Supports bitmap AND or OR operation.
// The space of second operand is freed, and the first operand serves as the results of the operation
ADIOS_ALAC_BITMAP * adios_alac_bitsOp(ADIOS_ALAC_BITMAP * op1,
		ADIOS_ALAC_BITMAP * op2, enum ADIOS_CLAUSE_OP_MODE operator) {
	int64_t i = 0;
	if (operator == ADIOS_QUERY_OP_AND) {
		for (i = 0; i < op1->length; i++) {
			op1->bits[i] &= op2->bits[i];
		}

		FreeALACBITMAP(op2);
//		free(op2->bits);
//		op1->elmSize = calSetBitsNum(op1);

	} else if (operator == ADIOS_QUERY_OP_OR) {

		for (i = 0; i < op1->length; i++) {
			op1->bits[i] |= op2->bits[i];
		}

		FreeALACBITMAP(op2);
	//		free(op2->bits);
//		op1->elmSize = calSetBitsNum(op1);

	} else {
		printf("Operator[%d] is not surpported now \n ", operator);
	}
	return op1;
}


static inline int coordinateConversionWithCheck(uint64_t * coordinates, const  int dim, const  uint64_t *srcstart, const  uint64_t *deststart, const  uint64_t *destend){

	int i = 0;
	for (i = 0; i < dim; i++) {
		coordinates[i] += (srcstart[i] /*global coordinate*/ );
		if ( coordinates[i] < deststart[i] || coordinates[i] > destend[i]){
			return (i+1) * -1;
		}
	}

	/*change coordinate to the destination box*/
	for (i = 0; i < dim; i++) {
		coordinates[i] -=  deststart[i];
	}
	return 1;
}


/* Give a rid that is relative to a src region
 * return a rid that is relative to dest selection box
 * Assume all the start & count array has slowest dimension at first position
 * NOTE: some redundant codes are for performance sake
 */
static inline bool ridConversionWithCheck(const rid_t rid/*relative to local src selectoin*/, const uint64_t *srcstart, const uint64_t *srccount, const uint64_t *deststart, const uint64_t *destcount,
		const int dim, rid_t *relativeRid , const int Corder ){

	int i = 0;
	uint64_t coordinates[MAX_DIMS];
	uint64_t destend[MAX_DIMS];

	for(i=0; i < dim ; i ++){
		destend[i] = deststart[i] + destcount[i] -1;
	}

	ridToCoordinates(dim, Corder, rid, srccount, coordinates);

	if (coordinateConversionWithCheck(coordinates, dim, srcstart, deststart, destend) < 0) {
			return false;
	}

	rid_t tmpRelativeRid = 0;
	if (Corder) {
		for (i = 0; i < dim; i++){
			tmpRelativeRid *= destcount[i];
			tmpRelativeRid += coordinates[i];
		}
	} else {
		for (i = dim-1; i >= 0; i--){
			tmpRelativeRid *= destcount[i];
			tmpRelativeRid += coordinates[i];
		}
	}

	*relativeRid = tmpRelativeRid;
	return true;
}


/*
 * TODO:
 * handle with c-order and fortran order
 */
static inline void ridToCoordinates(const int dim, const int Corder, const rid_t rid, const uint64_t * const dimSize , uint64_t * const coordinates /*OUT*/) {
	int i;
	rid_t remain = rid;
	if (Corder) {	// c-order
		for (i = dim - 1; i >= 0; i--) {
			coordinates[i] = remain % dimSize[i];
			remain /= dimSize[i];
		}
	} else { //f-order
		for (i = 0; i < dim; i++) {
			coordinates[i] = remain % dimSize[i];
			remain /= dimSize[i];
		}
	}
}

/*
 * Corder: 1 => c-order , 0=> fortran-order
		count={2,5, 7}
		coord= {1,3,4}
		c-order :    1 * 5 * 7 + 3 * 7 + 4 = 35 + 21 + 4 =
		expected RID: 4 + 3*7 + 1*35 = 4 + 21 + 35 = 60
		relativeRID = 0	(0)
		relativeRID *= 2	(0)
		relativeRID += 1	(1)
		relativeRID *= 5	(5)
		relativeRID += 3	(8)
		relativeRID *= 7	(56)
		relativeRID += 4	(60)

	( (0*2 + 1 ) *5 + 3 )* 7 + 4 = 1*5*7 + 3*7 + 4
 */

static inline rid_t ridConversionWithoutCheck(const rid_t rid/*relative to local src selectoin*/,
		const uint64_t *srcstart, const uint64_t *srccount, const uint64_t *deststart, const uint64_t *destcount,
		const int dim, const int Corder){

	int i;
	uint64_t coordinates[MAX_DIMS];

	// First convert the RID in the src box to coordinates in the dest box
	ridToCoordinates(dim, Corder, rid, srccount, coordinates);
	for (i = 0; i < dim; i++)
		coordinates[i] = coordinates[i] + srcstart[i] - deststart[i];

	// Then, convert the coordinates in the dest box back to an RID
	uint64_t relativeRid = 0;
	if (Corder) {
		for (i = 0; i < dim; i++){
			relativeRid *= destcount[i];
			relativeRid += coordinates[i];
		}
	} else {
		for (i = dim-1; i >= 0; i--){
			relativeRid *= destcount[i];
			relativeRid += coordinates[i];
		}
	}

	return relativeRid;
}

/*
 * usage: this MACRO tightly depends on the usage context,
 * which requires the context to declare the exact variables
  */
#define CHECK_ELEMENT(code) { \
		for(el= 0; el < totalElm; el ++){              \
			if ((code)){                               \
				rid = decodeRids[el];                  \
				double tmps = dclock();                 \
				bool flag = ridConversionWithCheck(rid, srcstart,srccount, deststart,destcount, dim, &newRid, Corder); \
				ckSpatialCheck += (dclock()-tmps);     \
				if (flag){  \
					double tmps2 = dclock();           \
					setBitsinBitMap(newRid, alacResultBitmap);       \
					alacResultBitmap->numSetBits ++;   \
					ckSetBitMap += (dclock() - tmps2);   \
				}                                      \
			}                                          \
		}                                              \
}

// printf("%"PRIu32"->%"PRIu32",", rid, newRid);

/*
 * NOTE : relies on  a lots of variables:
 * isPGCovered, el, srcstart, srcount, deststart, destcount, ndim, newRid, totalElm, op, Corder
 */
#define CHECK_NODECODE_ELEMENT(code) {                         \
		if (isPGCovered){                             \
			for(el= 0; el < totalElm; el ++){         \
				if (code) {                           \
					newRid = ridConversionWithoutCheck(el, srcstart,srccount, deststart,destcount, ndim,Corder);   \
					setBitsinBitMap(newRid, alacResultBitmap);    \
					alacResultBitmap->numSetBits ++;              \
				}                                                 \
			}                                                     \
		}else {                                                   \
			for(el= 0; el < totalElm; el ++){                     \
				if (code) {                                       \
					if (ridConversionWithCheck(el,	srcstart,srccount, deststart,destcount, ndim, &newRid,Corder)){  \
						setBitsinBitMap(newRid, alacResultBitmap); \
						alacResultBitmap->numSetBits ++;           \
					}                                              \
				}                                                  \
			}                                                      \
		}                                                          \
}

/*
 * NOTE: for ADIOS_NE (!=), we still use == to do the candidate check,
 * since we first treat != to be =, and at end, we will flip the bits
 */
#define CHECK_GENERIC_DATA(data, FUNC) {                       \
	switch(op){                                          \
		case(ADIOS_LT):                                  \
				FUNC((data)[el] < hb);            \
				break;                                   \
		case (ADIOS_LTEQ):                               \
				FUNC((data)[el] <= hb);           \
				break;                                   \
		case (ADIOS_GT):                                 \
				FUNC((data)[el] > lb);            \
				break;                                   \
		case (ADIOS_GTEQ):                               \
				FUNC((data)[el] >= lb);           \
				break;                                   \
		case (ADIOS_EQ):                                 \
				FUNC((data)[el] ==  lb);          \
				break;                                   \
		case (ADIOS_NE):                                 \
				FUNC((data)[el] == lb);           \
				break;                                   \
		default :                                        \
			printf("unknown predicate mode %d \n", op);  \
	}                                                    \
}



/*
 * return : 1 => yes, there is unsupported data type
 *          0 => no
 */
int checkUnsupportedDataType(enum ADIOS_DATATYPES dataType) {

	return !( (dataType == adios_real ) ||
			 (dataType == adios_double) ) ;
}


void reconstituteElementsUint64T(const ALMetadata *meta, bin_id_t start_bin, bin_id_t end_bin,
		                         char *start_bin_input, char *start_bin_output) {

	const ALBinLayout * const bl = &meta->binLayout;
	int insigbytes = alacrity_util_insigBytesCeil(meta);
	int insigbits = (meta->elementSize << 3) - meta->significantBits;

	uint64_t reconst_elem;
	uint64_t high_mask;
	uint64_t *next_output = (uint64_t*)start_bin_output;
	char *next_input = (char*)start_bin_input;

	// First reconstitute the values // Iterate through the bins backwards
	bin_id_t bin = start_bin;
	for (; bin < end_bin; bin++) {
		const bin_offset_t off_start = bl->binStartOffsets[bin];
		const bin_offset_t off_end = bl->binStartOffsets[bin + 1];

		high_mask = ((uint64_t)bl->binValues[bin]) << insigbits;
		// Iterate through the elements backwards
		bin_offset_t off = off_start;
		for (; off < off_end; off++) {
			GET_BUFFER_ELEMENT(reconst_elem, next_input, insigbytes);
			reconst_elem |= high_mask;
			*next_output++ = reconst_elem;
			next_input += insigbytes;
		}
	}
}

/*
 * reassemble the low-order byte & bin header value to the original value
 * for lower-order bins ( [start_bin , end_bin) )
 * I copied the implementation from the alacrity in order to profile the performance
 */
void reconstituteData2(const ALMetadata *meta, bin_id_t start_bin, bin_id_t end_bin,
		                         char *start_bin_input, char *start_bin_output){
	switch (meta->elementSize) {
	case sizeof(uint64_t):
		reconstituteElementsUint64T( meta, start_bin, end_bin, start_bin_input, start_bin_output);
		break;
	case sizeof(uint32_t):
		 //TODO
		break;
	case sizeof(uint16_t):
		//TODO
		break;
	case sizeof(uint8_t):
		//TODO:
		break;
	default:
		eprintf("Unsupported element size %d in %s\n", (int)meta->elementSize, __FUNCTION__);
		assert(false);
	}
	return ;
}


int adios_alac_check_candidate(ALMetadata *partitionMeta, bin_id_t startBin, bin_id_t endBin , double hb, double lb
		, uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount, int dim
		, ADIOS_QUERY * adiosQuery , const char *inputCurPtr /*index bytes of entire PG*/
		, bool decoded  /*true: need decoding */ , char * lowOrderBytes /*low order bytes of from startBin to endBin*/
		, enum ADIOS_DATATYPES dataType
		,ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/, int Corder){

#ifdef BREAKDOWN
	candidateCheckStart = dclock();
#endif

	const ALBinLayout * bl = &(partitionMeta->binLayout);
	/*assert(adiosQuery->_sel->type == ADIOS_SELECTION_BOUNDINGBOX);
	const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(adiosQuery->_sel->u.bb);
	uint64_t * destcount = bb->count;  uint64_t * deststart = bb->start; //region dimension of Selection box
	uint64_t *srcstart = pgBB->start; uint64_t *srccount = pgBB->count; //PG region dimension
	int dim = pgBB->ndim;*/
	uint32_t * decodeRids;
	bin_offset_t totalElm = bl->binStartOffsets[endBin] - bl->binStartOffsets[startBin];
	if (decoded){
		const uint64_t *compBinStartOffs = partitionMeta->indexMeta.u.ciim.indexBinStartOffsets;
		uint64_t  binCompressedLen = compBinStartOffs[endBin] - compBinStartOffs[startBin];
		decodeRids= (uint32_t*) malloc(sizeof(uint32_t)*totalElm);
		ALDecompressRIDs( inputCurPtr, binCompressedLen, decodeRids, &totalElm);
		assert(totalElm== bl->binStartOffsets[endBin] - bl->binStartOffsets[startBin]);
	}else{
		decodeRids = (uint32_t *) inputCurPtr;
	}

#ifdef BREAKDOWN
	ckTotalElm += totalElm;
	ckRids += (dclock()- candidateCheckStart);
	ckReconstituteStart = dclock();
#endif
	char * data = (char *) malloc(partitionMeta->elementSize*totalElm); // recovered data in bytes
	// reconstituteData(partitionMeta, startBin, endBin, lowOrderBytes, data); //this takes quite amount of time when the data size is large

	reconstituteData2(partitionMeta, startBin, endBin, lowOrderBytes, data);

#ifdef BREAKDOWN
	ckReconstitute += (dclock() - ckReconstituteStart);
#endif
	//Following variables are needed for the micros definition
	rid_t newRid;
	enum ADIOS_PREDICATE_MODE op = adiosQuery->predicateOp;
	bin_offset_t el = 0;   rid_t rid ;

	switch(dataType){
		case adios_unsigned_byte:
			CHECK_GENERIC_DATA((unsigned char  *)data, CHECK_ELEMENT );
			break;
		case adios_byte:
			CHECK_GENERIC_DATA((signed char  *)data, CHECK_ELEMENT );
			break;
		case adios_string:
			CHECK_GENERIC_DATA((char *)data, CHECK_ELEMENT );
			break;
		case adios_unsigned_short:
			CHECK_GENERIC_DATA((unsigned short  *)data, CHECK_ELEMENT );
			break;
		case adios_short:
			CHECK_GENERIC_DATA((signed short  *)data, CHECK_ELEMENT );
			break;
		case adios_unsigned_integer:
			CHECK_GENERIC_DATA((unsigned int *)data, CHECK_ELEMENT );
			break;
		case adios_integer:
			CHECK_GENERIC_DATA((signed int *)data, CHECK_ELEMENT );
			break;
		case adios_unsigned_long:
			CHECK_GENERIC_DATA((unsigned long long  *)data, CHECK_ELEMENT );
			break;
		case adios_long:
			CHECK_GENERIC_DATA((signed long long  *)data, CHECK_ELEMENT );
			break;
		case adios_real:
			CHECK_GENERIC_DATA((float *)data, CHECK_ELEMENT );
			break;
		case adios_double:
			CHECK_GENERIC_DATA((double *)data, CHECK_ELEMENT );
			break;
	}


	free(data);
	if(decoded)
		free(decodeRids);

#ifdef BREAKDOWN
	candidateCheckTotal += (dclock() - candidateCheckStart);
#endif

	return 0;

}

/*
 * go through every data value in the array, check whether the value satisfies the query condition
 * return: 1 ==> error, 0 ==> no error
 */
int literallyCheckData(void *data, uint64_t totalElm , enum ADIOS_DATATYPES dataType, ADIOS_QUERY * adiosQuery, double hb, double lb
		, uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount, int ndim, bool isPGCovered
		,ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/,int Corder){
	uint32_t newRid ;
	uint64_t el = 0;
	enum ADIOS_PREDICATE_MODE op = adiosQuery->predicateOp;
//	GENERIC_DATA_CONV(data, CHECK_NODECODE_ELEMENT );
	switch(dataType){
		case adios_unsigned_byte:
			CHECK_GENERIC_DATA((unsigned char  *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_byte:
			CHECK_GENERIC_DATA((signed char  *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_string:
			CHECK_GENERIC_DATA((char *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_unsigned_short:
			CHECK_GENERIC_DATA((unsigned short  *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_short:
			CHECK_GENERIC_DATA((signed short  *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_unsigned_integer:
			CHECK_GENERIC_DATA((unsigned int *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_integer:
			CHECK_GENERIC_DATA((signed int *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_unsigned_long:
			CHECK_GENERIC_DATA((unsigned long long  *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_long:
			CHECK_GENERIC_DATA((signed long long  *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_real:
			CHECK_GENERIC_DATA((float *)data, CHECK_NODECODE_ELEMENT );
			break;
		case adios_double:
			CHECK_GENERIC_DATA((double *)data, CHECK_NODECODE_ELEMENT );
			break;
	}

	free(data);
	return 0;

}
/*for each RID, set bit as 1 at corresponding position */
/*check the rid whether it belongs to user selection box
 if the current PG is not fully covered by user selection box*/
void setRidToBits(bool isPGCovered
		,uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int ndim, rid_t * idx, uint64_t totalRids,
		ADIOS_ALAC_BITMAP *alacResultBitmap /*IN&OUT*/, int Corder){
	uint32_t newRid ;
	if (isPGCovered){ // if PG is fully covered by the user bounding box, we dont need to check the new location(position) of converted point
		uint64_t ni;
		for ( ni = 0; ni < totalRids; ni++) {
			rid_t rid_val = idx[ni];
			newRid = ridConversionWithoutCheck(rid_val, srcstart,srccount, deststart,destcount, ndim, Corder);
#ifdef RIDBUG
			printf("%"PRIu32"->%"PRIu32",", rid_val, newRid);
#endif
			setBitsinBitMap(newRid, alacResultBitmap);
		}
		alacResultBitmap->numSetBits += totalRids;
	}else {
		uint64_t ni;
		for ( ni = 0; ni < totalRids; ni++) {
			rid_t rid_val = idx[ni];
			if (ridConversionWithCheck(rid_val,	srcstart,srccount, deststart,destcount, ndim, &newRid, Corder)){
#ifdef RIDBUG
			printf("%"PRIu32"->%"PRIu32",", rid_val, newRid);
#endif
				setBitsinBitMap(newRid, alacResultBitmap);
				alacResultBitmap->numSetBits ++;
			}
		}
	}
}



void resolveQueryBoundary(ADIOS_QUERY *adiosQuery, double *hb, double *lb)
{
	(*hb)= DBL_MAX;
	(*lb)= -DBL_MAX;
    if(adiosQuery->predicateOp == ADIOS_LT || adiosQuery->predicateOp == ADIOS_LTEQ){
    	*hb = atof(adiosQuery->predicateValue);
    }else if(adiosQuery->predicateOp == ADIOS_GT || adiosQuery->predicateOp == ADIOS_GTEQ){
    	*lb = atof(adiosQuery->predicateValue);
	}else if(adiosQuery->predicateOp == ADIOS_EQ){
		//following two cases are tricky to ALACRITY
		*hb = atof(adiosQuery->predicateValue);
		*lb = *hb;
	}else if(adiosQuery->predicateOp == ADIOS_NE){
		//TODO: flip the bits once the evaluation is done
		*hb = atof(adiosQuery->predicateValue);
		*lb = *hb;
	}else{
			printf("Unsupported predicate type[%d] \n", adiosQuery->predicateOp);
	}

}

/*
 * read low order bytes from 'low_bin'  bin to 'hi_bin'
 * allocated the data buffer in this function, requiring caller to free the readData
 */
char * readLowDataAmongBins(ALMetadata *partitionMeta, bin_id_t low_bin , bin_id_t hi_bin
		, uint64_t lowDataByteStartPos,  ADIOS_QUERY * adiosQuery
		, int blockId, int startStep, int numStep){
#ifdef BREAKDOWN
	dataStart = dclock();
#endif

	const uint64_t bin_read_len = ALGetDataBinOffset( partitionMeta, hi_bin) - ALGetDataBinOffset( partitionMeta, low_bin); // in bytes
	uint64_t lowDataBinOffset = lowDataByteStartPos+ ALGetDataBinOffset( partitionMeta, low_bin); /*element offset*/;
	// low order bytes from low_bin to hi_bin, ITS NOT entire low order byte
	char * readData= (char *) calloc(bin_read_len, sizeof(char));
	readLowOrderBytes(blockId,lowDataBinOffset, bin_read_len
	         , adiosQuery->file, adiosQuery->varinfo, startStep, numStep, (void *)readData);

#ifdef BREAKDOWN
	dataTotal = dataTotal + (dclock()- dataStart);
	dataReadSize += (bin_read_len/1024) ; // in KB
#endif
	return readData;
}


char * readIndexAmongBins(ALMetadata *partitionMeta, bin_id_t low_bin , bin_id_t hi_bin
		, uint64_t indexStartPos,  ADIOS_QUERY * adiosQuery
		, int blockId, int startStep, int numStep){

#ifdef BREAKDOWN
	idxStart = dclock();
#endif

	const uint64_t bin_read_len = ALGetIndexBinOffset( partitionMeta, hi_bin) - ALGetIndexBinOffset( partitionMeta, low_bin); // in bytes
	uint64_t lowDataBinOffset = indexStartPos+ ALGetIndexBinOffset( partitionMeta, low_bin); /*element offset*/;
	// low order bytes from low_bin to hi_bin, ITS NOT entire low order byte
	char * readData= (char *) calloc(bin_read_len, sizeof(char));
	readLowOrderBytes(blockId,lowDataBinOffset, bin_read_len
	         , adiosQuery->file, adiosQuery->varinfo, startStep, numStep, (void *)readData);
#ifdef BREAKDOWN
	idxTotal = idxTotal + (dclock()- idxStart);
	idxReadSize += (bin_read_len/1024); // in KB
#endif

	return readData;
}


void set_consecutive_bits(uint64_t relativeRid, uint64_t ones , ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/) {
	uint64_t * bitmap = alacResultBitmap->bits;
	uint64_t slot,  pass_bits, leftMask , pass_mask;
	uint64_t allOnesMask = (uint64_t) (~0); // must have casting, otherwise, will only generate 32 1s
	// starting point: relativeRid,  number of 1 bits: srccount[ndim-1]
	slot =  (relativeRid >> 6);
	pass_bits = relativeRid & 0x3F;

	/*************LITTLE ENDIAN : Significant bit goes to high address (right side)*******/
	/*       0000  | 111111111              */
	/*   pass_bits | */
	pass_mask = allOnesMask << pass_bits;
	if ( ones + pass_bits <= 64 ) { // the number of bits to be set in the same slot
		/*     111111111111111 |  0000
		 *                     | 64-ones-pass_bits
		 */
		leftMask = allOnesMask >> (64 - pass_bits - ones);
		bitmap[slot] |= (leftMask & pass_mask);
	}else { // if it has more than one slot bits to be set
		bitmap[slot] |= pass_mask ; // the first slot
		uint64_t leftOnes = ones - (64 - pass_bits);
		uint64_t spanSlot = leftOnes >> 6 ; // left ones span the number of slots
		uint64_t  remainOnes = leftOnes & 0x3F; // remaining ones
		uint64_t s = 1;
		for (s = 1; s <= spanSlot; ++ s) {
			bitmap[slot + s] |= allOnesMask ; // slot+s in generally doesn't exceed the length of bitmap
		}
		if ( remainOnes > 0) {
			 /*
			  *       1111111 | 0000
			  *    remain_ones|
			  */
			uint64_t remainMask = allOnesMask >> (64- remainOnes) ;
			bitmap[slot + spanSlot +1 ] |= remainMask;
		}

	}

	alacResultBitmap->numSetBits += ones;
}


/*  desstart     ,      destcount        ,,      srcstart (x)        , srccount
 *   (2, 3)      ,       (5, 6)           ,      (4, 5)          ,   (2, 3 )       => (4 + 2 -1, 5+3-1) = (5, 7)
 *    (2,3)___________________
 *         |__|__|__|__|__|__|           (2,3) 0 | 1 | 2 | 3 | 4 | 5 |
 *         |__|__|__|__|__|__|           (3,3) 6 | 7 | 8 | 9 | 10 | 11 |
 *         |__|_x|__|__|__|__|           (4,3) 12| 13| 14
 *         |__|__|__|__|__|__|
 *         |__|__|__|_y|__|__|
 *         |__|__|__|__|__|__|
 *
 *  // 1. calculate the first relative rids (4,5)
 *         (4-2) * 6 + (5-3)* 1 -1 = 2* 6 + 2  = 14
 *  // 2. calculate the next consecutive rids    ( 5,  5)  (5-2)*6 + (5-3)*1 = 18 + 2 = 20, so the delta = 6 (destcount[ndim-1])
 *  // 3. locate the bit start position using relativeRid in the bitmap, and set consecutive 3 consecutive 1s ( 3= srccount[ndim-1])
 *  // NOTE that length argument is not used at this point
 */
void setBits_for_wholePG_1_2_dim(uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int ndim, ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/) {

	uint64_t relativeRid = 0;
	uint64_t count = 1;
	int dim = 0;
	//1.calculate the first relativeRids (4,5)
	while ( dim < ndim ) {
		relativeRid += (srcstart[ndim - dim -1 ] - deststart[ndim - dim -1 ]) * count;
		count *= destcount[ndim - dim -1 ];
		dim ++;
	}
	//2. delta between any two consecutive rids
	uint64_t delta = destcount[ndim -1 ] ; //count / destcount[0];
	//2.2 times of how many relativeRids
	uint64_t times = 1;
	dim = 0;
	while (dim < ndim - 1) {
		times *= srccount[dim];
		dim ++;
	}
	int j ;
	uint64_t  ones = srccount[ndim -1 ]; // the number of 1s needs to be reflected in the bitmap

	for ( j = 1; j <= times; ++ j) {
		set_consecutive_bits(relativeRid, ones, alacResultBitmap);
		relativeRid += delta;
	}

}


/*
 *
0                                  1407
0,0,0                          0,0,1407
0,1,0                          0,1,1407
...
0,1079                         0,1079,1407

              relativeRid, delta, times, ones
uint64_t zero [4] = {0, 1408, 50000, 400};
 *
0,0,0     0,0,399
0,1,0     0,1,399      --- 400 -> delta
...
0,249,0   0,249,399

1,0,0     1,0,399   --- 1080*1408 = 1520640 => gap
1,1,0     1,1,399
..
1,249,0   1,249,399

...

199,0,0   199,0,399
199,1,0	  199,1,399
..
199,249,0  199,249,399
 *
 */
void setBits_for_wholePG_3_dim(uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int ndim, ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/) {

	uint64_t relativeRid = 0, tmpRid = 0;
	uint64_t count = 1;
	int dim = 0;
	//1.calculate the first relativeRids
	while ( dim < ndim ) {
		relativeRid += (srcstart[ndim - dim -1 ] - deststart[ndim - dim -1 ]) * count;
		count *= destcount[ndim - dim -1 ];
		dim ++;
	}
	//2. delta between any two consecutive rids
	uint64_t delta = destcount[ndim -1 ] , gap = destcount[1]*destcount[2];

	int i = 0, j = 0 ;
	uint64_t  ones = srccount[ndim -1 ]; // the number of 1s needs to be reflected in the bitmap

	for (i=0; i < srccount[0]; i++) {
	   tmpRid = relativeRid;
	   for (j=0; j < srccount[1]; j ++ ) {
		   set_consecutive_bits(tmpRid,ones, alacResultBitmap);
		   tmpRid += destcount[2];
	   }
	   relativeRid += gap; // move to another dimension, need to calculate the jump again
	}
}


uint64_t getRealBitNum(ADIOS_ALAC_BITMAP * alacResultBitmap) {

	uint64_t kk = 0, count = 0, pcount = 0;
	uint64_t * p_bitmap = alacResultBitmap->bits;
	for (; kk < alacResultBitmap->length; kk++) {

			count =  bits_in_char [p_bitmap[kk] & 0xff]
								   +  bits_in_char [(p_bitmap[kk] >>  8) & 0xff]
								   +  bits_in_char [(p_bitmap[kk] >> 16) & 0xff]
								   +  bits_in_char [(p_bitmap[kk] >> 24) & 0xff]
								   +  bits_in_char [(p_bitmap[kk] >> 32) & 0xff]
								   +  bits_in_char [(p_bitmap[kk] >> 40) & 0xff]
								   +  bits_in_char [(p_bitmap[kk] >> 48) & 0xff]
								   +  bits_in_char [(p_bitmap[kk] >> 56) & 0xff]
								   ;
		 pcount += count;
	}
	return pcount;
}



void deltaBitSetsWithPGCoveredCorderIn3D(
		 uint32_t *input /*decoded detla list, with first element as frame of reference and the rest values are detlas*/
		, uint32_t inputLength
		, uint64_t *srcstart //PG region dimension
		, uint64_t *srccount
		, uint64_t *srcboundary //[ srccount[2], srccount[2]*srccount[1], srccount[2]*srccount[1]*srccount[0]
		, uint64_t *deststart //region dimension of Selection box
		, uint64_t *destcount
		, uint64_t *destboundary // [destcount[2], destcount[2]*destcount[1], destcount[2]*destcount[1]*destcount[0]
		, int ndim
		,ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/
		) {

	uint32_t x = 0;

	for (x = 0; x < inputLength ; x += PFORDELTA_CHUNK_SIZE) { // perform for every pfordetal chunk (128 elements)
		uint32_t realSize = x + PFORDELTA_CHUNK_SIZE >= inputLength ? inputLength - x : PFORDELTA_CHUNK_SIZE ;
		uint32_t first_rid_val = input[x];
		int i;
		uint64_t deltaCoordPG[ndim], coordRegion[ndim]; // coordinate in pg context

		// First convert the RID in the src box to first_coordinates in the dest box
		ridToCoordinates(ndim, 1, first_rid_val, srccount, deltaCoordPG);
		for (i = 0; i < ndim; i++){
			coordRegion[i] = deltaCoordPG[i] + srcstart[i] - deststart[i];
		}

		// Then, convert the coordinates in the dest box back to an RID
		uint32_t relativeRid = 0;
		for (i = 0; i < ndim; i++){
			relativeRid *= destcount[i];
			relativeRid  += coordRegion[i];
		}

		uint32_t word = (uint32_t) (relativeRid >> 6);
		alacResultBitmap->bits[word]
					|= (1LL << (relativeRid & 0x3F));

		uint32_t j = x+1, r = 0; // skip the first frame of the reference NOTE: j = x + 1
		//results[r++] = relativeRid;
		uint32_t nextRelativeRid = 0;
		for ( j = x + 1; j < x + realSize ; ++ j ) {
			uint32_t delta = input[j];
			uint64_t nextDim = delta + deltaCoordPG[ 2 ];
			//update first_coordinate with detal
			if (  nextDim < srcboundary[ 2]  ) { // have to be ndim (=3) specific , srcboundary[2]= srccount[2]
				nextRelativeRid = relativeRid + delta;
				deltaCoordPG[2 ] = nextDim; // along with the coordinates
			}else  if (deltaCoordPG[1] * srccount[2] + nextDim < srcboundary[1]) { // srcboundary[1] = srccount[2]*srccount[1]
				uint64_t inc = nextDim / srccount[2];  // increment dimension
				uint64_t co2 = nextDim % srccount[2]; // it is the new dim number
				if ( co2 > deltaCoordPG[2]) {
					nextRelativeRid = relativeRid +  inc * destboundary[2] + (co2 - deltaCoordPG[2]) ;
				}else {
					nextRelativeRid = relativeRid +  inc * destboundary[2] - (deltaCoordPG[2] - co2) ;
				}
				deltaCoordPG[2] = co2;  // along with the coordinates
			    deltaCoordPG[1] += inc;
			}else {
				uint64_t nextDim2 = (deltaCoordPG[1] * srcboundary[2] + nextDim);
				uint64_t inc0 =  nextDim2 / srcboundary[1];  // this is the delta dimsion on the highest dimension
				uint64_t remain0 = nextDim2 % srcboundary[1];
				uint64_t coord1 = remain0 / srcboundary[2];  //they are the new dim number
				uint64_t coord2 = remain0 % srcboundary[2];
				if (coord1 > deltaCoordPG[1]) { // we deal with unsign data type, it doesn't have sign
					if ( coord2 > deltaCoordPG[2]) {
						nextRelativeRid = relativeRid + inc0 * destboundary[1] +  (coord1 - deltaCoordPG[1]) * destboundary[2] + (coord2 - deltaCoordPG[2]) ;
					}else {
						nextRelativeRid = relativeRid + inc0 * destboundary[1] +  (coord1 - deltaCoordPG[1]) * destboundary[2] - (deltaCoordPG[2] - coord2) ;
					}

				}else {
					if ( coord2 > deltaCoordPG[2]) {
						nextRelativeRid = relativeRid + inc0 * destboundary[1] -  (deltaCoordPG[1] - coord1 ) * destboundary[2] + (coord2 - deltaCoordPG[2]) ;
					}else {
						nextRelativeRid = relativeRid + inc0 * destboundary[1] -  (deltaCoordPG[1] - coord1 ) * destboundary[2] - (deltaCoordPG[2] - coord2) ;
					}
				}
				deltaCoordPG[0] +=inc0;
				deltaCoordPG[1] = coord1;
				deltaCoordPG[2] = coord2;
			}

			relativeRid = nextRelativeRid;  // update the relativeRid

			word = (uint32_t) (relativeRid >> 6);
			alacResultBitmap->bits[word] |= (1LL << (relativeRid & 0x3F));
			//results[r++] = relativeRid;

		}
	}
}





void deltaBitSetsWithPGPartialCoveredCorderIn3D(
		 uint32_t *input /*decoded detla list, with first element as frame of reference and the rest values are detlas*/
		, uint32_t inputLength
		, uint64_t *srcstart //PG region dimension
		, uint64_t *srccount
		, uint64_t *srcboundary //[ srccount[2], srccount[2]*srccount[1], srccount[2]*srccount[1]*srccount[0]
		, uint64_t *deststart //region dimension of Selection box
		, uint64_t *destcount
		, uint64_t *destboundary // [destcount[2], destcount[2]*destcount[1], destcount[2]*destcount[1]*destcount[0]
		, int ndim
		,ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/
		) {
	uint32_t x = 0;
		uint64_t destend [ndim] ;
		destend[0] = deststart[0] + destcount[0] -1 ; destend[1] = deststart[1] + destcount[1] -1 ; destend[2] = deststart[2] + destcount[2] -1 ;
	for (x = 0; x < inputLength ; x += PFORDELTA_CHUNK_SIZE) {
		uint32_t realSize = x + PFORDELTA_CHUNK_SIZE >= inputLength ? inputLength - x : PFORDELTA_CHUNK_SIZE ;

		uint32_t first_rid_val = input[x];
		int i;
		uint64_t deltaCoordPG[ndim], coordRegion[ndim], tmpCoord[ndim]; // coordinate in pg context

		// First convert the RID in the src box to first_coordinates in the dest box
		ridToCoordinates(ndim, 1, first_rid_val, srccount, deltaCoordPG);
		for (i = 0; i < ndim; i++){
			coordRegion[i] = deltaCoordPG[i] + srcstart[i] - deststart[i];
		}
		// the first rid is not in the bounding box, then, the entire chunk is not as well
		if ( coordRegion[0] < deststart[0] || coordRegion[0] > destend[0] ) continue  ;
		if ( coordRegion[1] < deststart[1] || coordRegion[1] > destend[1] ) continue  ;
		if ( coordRegion[2] < deststart[2] || coordRegion[2] > destend[2] ) continue  ;

		// Then, convert the coordinates in the dest box back to an RID
		uint32_t relativeRid = 0;
		for (i = 0; i < ndim; i++){
			relativeRid *= destcount[i];
			relativeRid  += coordRegion[i];
		}

		uint32_t word = (uint32_t) (relativeRid >> 6);
		alacResultBitmap->bits[word]
					|= (1LL << (relativeRid & 0x3F));

		uint32_t j = x + 1, r = 0; // skip the first frame of the reference
		//results[r++] = relativeRid;
		uint32_t nextRelativeRid = 0;
		for ( j = x + 1; j < x + realSize ; ++ j ) {
			uint32_t delta = input[j];
			uint64_t nextDim = delta + deltaCoordPG[ 2 ];
			//update first_coordinate with detal
			if (  nextDim < srcboundary[ 2]  ) { // have to be ndim (=3) specific , srcboundary[2]= srccount[2]
				nextRelativeRid = relativeRid + delta;
				deltaCoordPG[2 ] = nextDim; // along with the coordinates
			}else  if (deltaCoordPG[1] * srccount[2] + nextDim < srcboundary[1]) { // srcboundary[1] = srccount[2]*srccount[1]
				uint64_t inc = nextDim / srccount[2];  // increment dimension
				uint64_t co2 = nextDim % srccount[2]; // it is the new dim number
				if ( co2 > deltaCoordPG[2]) {
					nextRelativeRid = relativeRid +  inc * destboundary[2] + (co2 - deltaCoordPG[2]) ;
				}else {
					nextRelativeRid = relativeRid +  inc * destboundary[2] - (deltaCoordPG[2] - co2) ;
				}
				deltaCoordPG[2] = co2;  // along with the coordinates
			    deltaCoordPG[1] += inc;
			}else {
				uint64_t nextDim2 = (deltaCoordPG[1] * srcboundary[2] + nextDim);
				uint64_t inc0 =  nextDim2 / srcboundary[1];  // this is the delta dimsion on the highest dimension
				uint64_t remain0 = nextDim2 % srcboundary[1];
				uint64_t coord1 = remain0 / srcboundary[2];  //they are the new dim number
				uint64_t coord2 = remain0 % srcboundary[2];
				if (coord1 > deltaCoordPG[1]) { // we deal with unsign data type, it doesn't have sign
					if ( coord2 > deltaCoordPG[2]) {
						nextRelativeRid = relativeRid + inc0 * destboundary[1] +  (coord1 - deltaCoordPG[1]) * destboundary[2] + (coord2 - deltaCoordPG[2]) ;
					}else {
						nextRelativeRid = relativeRid + inc0 * destboundary[1] +  (coord1 - deltaCoordPG[1]) * destboundary[2] - (deltaCoordPG[2] - coord2) ;
					}

				}else {
					if ( coord2 > deltaCoordPG[2]) {
						nextRelativeRid = relativeRid + inc0 * destboundary[1] -  (deltaCoordPG[1] - coord1 ) * destboundary[2] + (coord2 - deltaCoordPG[2]) ;
					}else {
						nextRelativeRid = relativeRid + inc0 * destboundary[1] -  (deltaCoordPG[1] - coord1 ) * destboundary[2] - (deltaCoordPG[2] - coord2) ;
					}
				}
				deltaCoordPG[0] +=inc0;
				deltaCoordPG[1] = coord1;
				deltaCoordPG[2] = coord2;
			}

			relativeRid = nextRelativeRid;  // update the relativeRid

			tmpCoord[0] = deltaCoordPG[0] + srcstart[0] - deststart[0];
//			tmpCoord[0] = deltaCoordPG[0] + srcstart[0] ;
			if ( tmpCoord[0] <  deststart[0] || tmpCoord[0] > destend[0]) continue ;
			tmpCoord[1] = deltaCoordPG[1] + srcstart[1] - deststart[1];
//			tmpCoord[1] = deltaCoordPG[1] + srcstart[1] ;
			if ( tmpCoord[1] <  deststart[1] || tmpCoord[1] > destend[1]) continue ;
			tmpCoord[2] = deltaCoordPG[2] + srcstart[2] - deststart[2];
//			tmpCoord[2] = deltaCoordPG[2] + srcstart[2] ;
			if ( tmpCoord[2] <  deststart[2] || tmpCoord[2] > destend[2]) continue ;

			word = (uint32_t) (relativeRid >> 6);
			alacResultBitmap->bits[word] |= (1LL << (relativeRid & 0x3F));
//			results[r++] = relativeRid;

		}
	}

}



void proc_write_block(int gBlockId /*its a global block id*/, bool isPGCovered, ADIOS_VARTRANSFORM *ti, ADIOS_QUERY * adiosQuery, int startStep, bool estimate
		, ALUnivariateQuery * alacQuery , double lb , double hb
		,uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount
		, ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/ , int Corder){

#ifdef RIDBUG
	printf("PG [%d], RIDs relative to PG converted to RIDs relative to output BoundingBox: ", gBlockId);
#endif

	int numStep = 1; // only deal with one timestep
//	uint64_t metaSize, indexSize, dataSize;
	ADIOS_VARINFO * varInfo = adiosQuery->varinfo;
	int ndim = varInfo->ndim;
	ADIOS_TRANSFORM_METADATA * tmetas = ti->transform_metadatas;
	ADIOS_TRANSFORM_METADATA tmeta = tmetas[gBlockId];
	//	assert(tmeta->length == 24);

	adios_transform_alacrity_metadata *alac_metadata = (adios_transform_alacrity_metadata *) malloc(sizeof(adios_transform_alacrity_metadata));

#ifdef BREAKDOWN
	metaStart = dclock();
#endif

	//the metadata order is alacrity meta, alacrity index, alacrity LoB, raw data (original data)
	read_alacrity_transform_metadata(tmeta.length, tmeta.content, alac_metadata);
	//It now assumes LoB along with OD (original data)
#ifdef BREAKDOWN
	metaTotal = metaTotal + (dclock()- metaStart);
	alacPartitionMetaStart = dclock();
#endif

	//TODO: offset of each PG should be included
//	printf("PG[%d] has meta size[ %" PRIu64 "], index size[ %" PRIu64 "], and data size[ %" PRIu64 "] \n", blockId, metaSize,  indexSize, dataSize);

	// transformed data has 1 dimension,
	// 1. load partition Metadata
	// NOTE: One ALACRITY PG Data is written in the below format:  [meta data] | [low order bytes data] | [ index data]
	ALMetadata partitionMeta;
	readPartitionMeta(gBlockId, alac_metadata->meta_size,adiosQuery->file, varInfo
					,startStep,numStep,&partitionMeta);

#ifdef BREAKDOWN
	//printf("read partition meta %d with size %"PRIu64"\n", gBlockId, alac_metadata->meta_size);
	alacPartitionMetaTotal += dclock() - alacPartitionMetaStart;
#endif

	const uint8_t insigbytes = alacrity_util_insigBytesCeil(&partitionMeta);

	//2. find touched bin
	bin_id_t low_bin, hi_bin;
	_Bool are_bins_touched = findBinRange1C(&partitionMeta, alacQuery, &low_bin, &hi_bin);

	if (!are_bins_touched) {
		FREE(alac_metadata);
		return ;
	}

	const ALBinLayout * bl = &(partitionMeta.binLayout);


#ifdef BREAKDOWN
	numTouchedBins  += (hi_bin - low_bin);

	printf("****Touched PG[%d]: the # of total bins[%"PRIu32"] & touched bins [%"PRIu32"], each touched bin has element #: [", gBlockId, bl->numBins , (hi_bin- low_bin));

	bin_id_t bin ;
	for ( bin = low_bin; bin < hi_bin; bin++) {
		printf("%"PRIu32",", bl->binStartOffsets[bin+1] - bl->binStartOffsets[bin]);
	}
	printf("]\n");
#endif

	bool fullyContained = false;
	if (hi_bin - low_bin == bl->numBins) {

		// with regards to the 3th & 4th arguments,  they are 0 or 1, false or true. You want per-block statistics, so call with (fp, vi, 0, 1) arguments
		common_read_inq_var_stat(adiosQuery->file, varInfo, 0, 1); // retrieve the statistics structure per-block

#ifdef BREAKDOWN
		printf("****Touched PG[%d] has all bins retrieved \n", gBlockId);
		printf("****Touched PG[%d] is fully covered by the selection region [%d] \n", gBlockId, isPGCovered);
		printf("****Var info has statistics: [%s]  \n", varInfo->statistics != NULL ? "yes" : "no");
#endif

		if (isPGCovered && varInfo->statistics && varInfo->statistics->blocks) {
				ADIOS_VARSTAT *stat = varInfo->statistics;
				double min =  *(double*)stat->blocks->mins[gBlockId];
				double max =  *(double*)stat->blocks->maxs[gBlockId];
				fullyContained =  (lb <= min && hb >= max) ;
		#ifdef BREAKDOWN
			    printf("*****PG [%d] is fully contained [%s], min [%f], max [%f], query low bound [%f], and query high bound [%f] \n ", gBlockId, fullyContained ? "true": "false", min, max, lb, hb);
		#endif
		}
	}

	//3. load index size
	uint64_t indexStartPos =  alac_metadata->index_offset ;
	char * index = readIndexAmongBins(&partitionMeta
								, low_bin,  hi_bin, indexStartPos, adiosQuery, gBlockId, startStep, numStep);
	char * input_index = index;
	ALIndex* indexPtr = &index;


	//TODO: distinguish the offset btw two bins for compressed and uncompressed index
	// is the offset byte-level or element-level?
	if (partitionMeta.indexMeta.indexForm == ALInvertedIndex) {
		// indexes are inverted indexes that are not compressed,  we build bitmaps for each rid;
		//element offset, instead of byte element
		uint64_t resultCount = bl->binStartOffsets[hi_bin] - bl->binStartOffsets[low_bin];
		if (estimate || fullyContained) {

if ( fullyContained ) {
#ifdef BREAKDOWN
		    printf("*****PG [%d] is fully contained in both spatial and value selections, we set all bits at once without doing ALACRITY process\n ", gBlockId);
#endif
}

			 // fullyContained: since this PG is fully contained in the spatial selection as well as in the value constraint, we now just set the entire bitmap, without doing the ALACRITY process
			setRidToBits(isPGCovered, srcstart, srccount, deststart, destcount, ndim
					, (rid_t *)index, resultCount, alacResultBitmap,Corder);
		} else {
			uint64_t lowByteStartPos2 = alac_metadata->lob_offset;
			rid_t * decodedRid = (rid_t *) index;

			// It touches at least 3 bins, so, we need to check RIDs that are in first and last bins
			if (hi_bin - low_bin > 2) {
				bin_offset_t lowBinElm = bl->binStartOffsets[low_bin + 1] - bl->binStartOffsets[low_bin];
				uint64_t hiBinElm = bl->binStartOffsets[hi_bin] - bl->binStartOffsets[hi_bin-1];

				// low boundary bin
				char * lobLowBin= readLowDataAmongBins(&partitionMeta
													,low_bin, low_bin+1,  lowByteStartPos2, adiosQuery, gBlockId, startStep, numStep);
				adios_alac_check_candidate(&partitionMeta, low_bin, low_bin+1 , hb, lb
						, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , (char*) decodedRid /*index bytes of entire PG*/
						, false  /*don't need decoding*/ , lobLowBin, varInfo->type
						,alacResultBitmap /*OUT*/,Corder);
				decodedRid += lowBinElm;
				FREE(lobLowBin);

				uint64_t innerElm = resultCount- lowBinElm - hiBinElm;

#ifdef BREAKDOWN
setRidStart = dclock();
#endif
				setRidToBits(isPGCovered, srcstart, srccount, deststart, destcount, ndim
											, decodedRid, innerElm, alacResultBitmap,Corder);

#ifdef BREAKDOWN
setRidTotal += (dclock() - setRidStart);
#endif
				decodedRid  +=  innerElm;

				// high boundary bin
				char * lobHighBin= readLowDataAmongBins(&partitionMeta
								  					,hi_bin-1, hi_bin,  lowByteStartPos2, adiosQuery, gBlockId, startStep, numStep);

				adios_alac_check_candidate(&partitionMeta, hi_bin-1, hi_bin , hb, lb
						, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , (char *)decodedRid
						, false , lobHighBin, varInfo->type
						,alacResultBitmap /*OUT*/,Corder);
				FREE(lobHighBin);

			} else { // for 1 or 2 bins touched, we need to check all RIDs

				char * lowOrderBytes2  = readLowDataAmongBins(&partitionMeta
													,low_bin, hi_bin,  lowByteStartPos2, adiosQuery, gBlockId, startStep, numStep);
				char *lowOrderPtr2 = lowOrderBytes2; // temporary pointer
				adios_alac_check_candidate(&partitionMeta, low_bin, hi_bin  , hb, lb
						, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , (char*)decodedRid
						, false , lowOrderPtr2, varInfo->type
						,alacResultBitmap /*OUT*/,Corder);
				FREE(lowOrderBytes2);
			}



		}

	}else if (partitionMeta.indexMeta.indexForm == ALCompressedInvertedIndex) {
		const uint64_t *compBinStartOffs = partitionMeta.indexMeta.u.ciim.indexBinStartOffsets;
		uint64_t binCompressedLen;
		const char *inputCurPtr = input_index;

		if (estimate || fullyContained) {
			 // fullyContained: since this PG is fully contained in the spatial selection as well as in the value constraint, we now just set the entire bitmap, without doing the ALACRITY process
			// Now compress each bin in turn
if ( fullyContained ) {
#ifdef BREAKDOWN
			printf("*****PG [%d] is fully contained in both spatial and value selections, we set all bits at once without doing ALACRITY process\n ", gBlockId);
#endif
}

#ifdef BREAKDOWN
double deS = dclock();
uint64_t beforeBits = alacResultBitmap->numSetBits;
#endif

		if ( ndim <= 2 ) {
				setBits_for_wholePG_1_2_dim(srcstart, srccount /*PG region dimension*/ ,
						deststart, destcount /*region dimension of Selection box*/
						, ndim, alacResultBitmap);
		}else if ( ndim == 3) {
				setBits_for_wholePG_3_dim(srcstart, srccount /*PG region dimension*/ ,
						deststart, destcount /*region dimension of Selection box*/
						, ndim, alacResultBitmap);
		}else {  // the setBits_for_wholePG is hard to be generalized to any dimension.
				 // basically, if there is a ndim, we have to have ndim-1 for-loops inside of the function
				bin_id_t bin ;
				for ( bin = low_bin; bin < hi_bin; bin++) {
					binCompressedLen = compBinStartOffs[bin + 1] - compBinStartOffs[bin];
					uint32_t decodedElm = ALDecompressRIDtoSelBox(isPGCovered , inputCurPtr, binCompressedLen
							, srcstart, srccount /*PG region dimension*/ , deststart, destcount /*region dimension of Selection box*/
							, ndim , (void **)&(alacResultBitmap->bits));
					inputCurPtr += binCompressedLen;
					alacResultBitmap->numSetBits += decodedElm;
				}
		}
#ifdef BREAKDOWN
//uint64_t realBit = getRealBitNum(alacResultBitmap);
//if(ndim == 3)
//printf("BitSet block_id: %d  srcstart: %"PRIu64" %"PRIu64" %"PRIu64" srccount: %"PRIu64" %"PRIu64" %"PRIu64" deststart: %"PRIu64" %"PRIu64" %"PRIu64" destcount: [%"PRIu64" %"PRIu64" %"PRIu64"] bits_before: %"PRIu64" bits_after: %"PRIu64" real_bit: %"PRIu64" exist: %s bits_add: %"PRIu64"\n"
//			, gBlockId, srcstart[0],  srcstart[1], srcstart[2], srccount[0], srccount[1], srccount[2], deststart[0], deststart[1], deststart[2], destcount[0], destcount[1], destcount[2], beforeBits, alacResultBitmap->numSetBits, realBit, realBit == alacResultBitmap->numSetBits ? "true" : "false", alacResultBitmap->numSetBits - beforeBits );

double dt= dclock() - deS;
wholeDecodeTotal += dt;
decodeTotal += dt;
#endif

		}else{
			// element count of touched bins, which is also the count of low order data
			uint64_t lowByteStartPos = alac_metadata->lob_offset;


			// It touches at least 3 bins, so, we need to check RIDs that are in first and last bins
			if (hi_bin - low_bin > 2) {
				// low boundary bin, compressed byte offset
				binCompressedLen = compBinStartOffs[low_bin + 1] - compBinStartOffs[low_bin];

				// low boundary bin
				char * lobLowBin= readLowDataAmongBins(&partitionMeta
													,low_bin, low_bin+1,  lowByteStartPos, adiosQuery, gBlockId, startStep, numStep);
				adios_alac_check_candidate(&partitionMeta, low_bin, low_bin+1 , hb, lb
						, srcstart, srccount, deststart, destcount, ndim,  adiosQuery  , inputCurPtr /*index bytes of entire PG*/
						, true  /*need decoding*/ , lobLowBin /*it points to the start of `low_bin` */, varInfo->type
						,alacResultBitmap /*OUT*/,Corder);
				inputCurPtr += binCompressedLen;
				FREE(lobLowBin);


				bin_id_t innerlowBin = low_bin + 1;
				bin_id_t innerHiBin = hi_bin -1;
				bin_id_t bin;


				if ( ndim == 3 && Corder == 1) {

						uint64_t * destboundary = (uint64_t *) malloc(sizeof(uint64_t)*ndim);
						destboundary[0]  = destcount[0] * destcount[1]* destcount[2]; destboundary[1] = destcount [1]* destcount[2]; destboundary[2] = destcount[2];
						uint64_t * srcboundary = (uint64_t *) malloc(sizeof(uint64_t)*ndim);
						srcboundary[0] = srccount[0] * srccount[1]*srccount[2]; srcboundary [1] = srccount[1]*srccount[2]; srcboundary[2] = srccount[2];
						bin_offset_t binLen, totalBinLen = 0;
						if (isPGCovered) {

#ifdef BREAKDOWN
double pgC= dclock();
#endif
							for ( bin = innerlowBin; bin < innerHiBin; bin++) {
								binCompressedLen = compBinStartOffs[bin + 1] - compBinStartOffs[bin];
								binLen = bl->binStartOffsets[bin + 1] - bl->binStartOffsets[bin];
								rid_t *output_index = malloc(binLen * sizeof(rid_t));
								rid_t *outputCurPtr = output_index;
								uint64_t beforeBits = alacResultBitmap->numSetBits;
								ALDecompressDeltas(inputCurPtr, binCompressedLen, outputCurPtr, &binLen);
								deltaBitSetsWithPGCoveredCorderIn3D(outputCurPtr, binLen, srcstart, srccount
								,srcboundary, deststart, destcount, destboundary, ndim, alacResultBitmap);
								alacResultBitmap->numSetBits += binLen;
								totalBinLen += binLen;
								inputCurPtr += binCompressedLen;
								free(output_index);

							}

#ifdef BREAKDOWN
				double pgcTmpTotal = (dclock() - pgC);
				pgCTotal += pgcTmpTotal;
				decodeTotal += pgcTmpTotal;
#endif
						}else {

#ifdef BREAKDOWN
				double pgPartialC= dclock();
#endif
								for ( bin = innerlowBin; bin < innerHiBin; bin++) {
									binCompressedLen = compBinStartOffs[bin + 1] - compBinStartOffs[bin];
									binLen = bl->binStartOffsets[bin + 1] - bl->binStartOffsets[bin];
									rid_t *output_index = malloc(binLen * sizeof(rid_t));
									rid_t *outputCurPtr = output_index;

								ALDecompressDeltas(inputCurPtr, binCompressedLen, outputCurPtr, &binLen);
								deltaBitSetsWithPGPartialCoveredCorderIn3D(outputCurPtr, binLen, srcstart, srccount
								,srcboundary, deststart, destcount, destboundary, ndim, alacResultBitmap);
									alacResultBitmap->numSetBits += binLen;
									inputCurPtr += binCompressedLen;
									free(output_index);
								}
#ifdef BREAKDOWN
				double pgTmpPartialTotal = (dclock() - pgPartialC);
				pgPartialTotal += pgTmpPartialTotal ;
				decodeTotal += pgTmpPartialTotal;
#endif
					}
					free(destboundary);
					free(srcboundary);
				} else {
								// Now compress each bin in turn
								for ( bin = innerlowBin; bin < innerHiBin; bin++) {
									binCompressedLen = compBinStartOffs[bin + 1] - compBinStartOffs[bin];

#ifdef BREAKDOWN
decodeStart= dclock();
#endif
									uint32_t decodedElm = ALDecompressRIDtoSelBox(isPGCovered
											, inputCurPtr, binCompressedLen
											, srcstart, srccount //PG region dimension
											, deststart, destcount //region dimension of Selection box
											, ndim , (void**)&(alacResultBitmap->bits));
#ifdef BREAKDOWN
decodeTotal += dclock() - decodeStart;
#endif
									alacResultBitmap->numSetBits += decodedElm;
									inputCurPtr += binCompressedLen;
								}

				}



				// high boundary bin
				binCompressedLen = compBinStartOffs[hi_bin]- compBinStartOffs[hi_bin-1];
				char * lobHighBin= readLowDataAmongBins(&partitionMeta
												  					,hi_bin-1, hi_bin,  lowByteStartPos, adiosQuery, gBlockId, startStep, numStep);

				adios_alac_check_candidate(&partitionMeta, hi_bin-1, hi_bin , hb, lb
						, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , inputCurPtr /*index bytes of entire PG*/
						, true  /*need decoding*/ , lobHighBin /*low order bytes of entire PG*/, varInfo->type
						,alacResultBitmap /*OUT*/,Corder);
				inputCurPtr += binCompressedLen;
				FREE(lobHighBin);

			} else { // for 1 or 2 bins touched, we need to check all RIDs
				char * lowOrderBytes  = readLowDataAmongBins(&partitionMeta
								, low_bin,  hi_bin, lowByteStartPos , adiosQuery, gBlockId, startStep, numStep);
				char *lowOrderPtr = lowOrderBytes; // temporary pointer
				adios_alac_check_candidate(&partitionMeta, low_bin, hi_bin , hb, lb
						, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , inputCurPtr /*index bytes of entire PG*/
						, true , lowOrderPtr, varInfo->type
						,alacResultBitmap /*OUT*/,Corder);
				FREE(lowOrderBytes);
			}


		}
	} else {
		printf("index form %d in alacrity is not supported", partitionMeta.indexMeta.indexForm);
		exit(EXIT_FAILURE);
	}
	FREE(input_index);

	FREE(alac_metadata);

#ifdef RIDBUG
	printf("\n");
#endif
}





/*
 * PG selection, [startPG, endPG)
 */
ADIOS_ALAC_BITMAP* adios_alac_uniengine(ADIOS_QUERY * adiosQuery, int timeStep, bool estimate) {


	if (checkUnsupportedDataType(adiosQuery->varinfo->type)){
		printf("unsupported data type [%d] at this point \n", adiosQuery->varinfo->type);
		exit(EXIT_FAILURE);
	}

	double lb , hb ;
    resolveQueryBoundary(adiosQuery, &hb, &lb); // query constraints
	//printf("constraint: %s\n", adiosQuery->condition);

	ADIOS_VARINFO * varInfo = adiosQuery->varinfo;

#ifdef BREAKDOWN
	printf("**************************\n");
	printf("process %s constraint \n", adiosQuery->condition);
	// initialization
	preparationTime = 0;
	numTouchedBins =0;
	numTouchedPGs = 0;
	findPGTotal =0;
	procTotal = 0;
	metaTotal =0;
	alacPartitionMetaTotal =0;
	idxTotal =0;
	idxReadSize = 0;
	dataTotal =0;
	dataReadSize = 0;
	candidateCheckTotal = 0;
	ckRids = 0; ckReconstitute = 0; ckSpatialCheck = 0; ckSetBitMap = 0;

	decodeTotal =0;
	wholeDecodeTotal=0;
	pgCTotal = 0;
	pgPartialTotal = 0;
	setRidTotal =0;
	preparationStart = dclock();
#endif
	adios_read_set_data_view(adiosQuery->file, LOGICAL_DATA_VIEW); // switch to the transform view,
 	ADIOS_VARTRANSFORM *ti = adios_inq_var_transform(adiosQuery->file, varInfo); // this func. will fill the blockinfo field


	int startStep = timeStep, numStep = 1;
	uint64_t totalElm = adiosQuery->rawDataSize; // no matter bounding box or writeblock selection, the rawDataSize has been calculated in the common query layer
	ADIOS_ALAC_BITMAP *alacResultBitmap =  (ADIOS_ALAC_BITMAP *) malloc(sizeof(ADIOS_ALAC_BITMAP ));
	alacResultBitmap->length = BITNSLOTS64(totalElm);
	//initially, no 1s at all
	alacResultBitmap->bits = (uint64_t *) calloc( alacResultBitmap->length , sizeof(uint64_t));
	alacResultBitmap->numSetBits = 0;
	alacResultBitmap->realElmSize = totalElm;

	// check the transform type, if it is alacrity,
	ALQueryEngine qe;
	ALUnivariateQuery alacQuery;
	if (ti->transform_type == adios_get_transform_type_by_uid("ncsu-alacrity")) { // if it is not alacrity type, do not initialize this query engine
		if ( adiosQuery->varinfo->type == adios_double){
			ALQueryEngineStartUnivariateDoubleQuery(&qe, lb, hb, REGION_RETRIEVAL_INDEX_ONLY_QUERY_TYPE, &alacQuery);
		}else if (adiosQuery->varinfo->type == adios_real){
			ALQueryEngineStartUnivariateFloatQuery(&qe, lb, hb, REGION_RETRIEVAL_INDEX_ONLY_QUERY_TYPE, &alacQuery);
		}
	}

	const int Corder = !futils_is_called_from_fortran(); // Use the dimension order of the caller; the common read layer will also mimic this order

#ifdef BREAKDOWN
	preparationTime += (dclock() - preparationStart) ;
#endif
	/*********** doQuery ***************
	 *
	 * 1. Open partition  [locate offsets of meta, data, and index for the partition]
	 * 2. Read Partition Meta from file => meta
	 3. find touched bins:  are_bins_touched = findBinRange1C(meta, uniquery, &start_bin, &end_bin);
	 4. read index of touched bins: ALPartitionStoreReadIndexBins
	 5. read dataBin
	 */
	uint64_t* deststart ;  uint64_t* destcount ;// current variables selection box
	uint64_t * srcstart;  	uint64_t * srccount; // PG's bounding box is the global bounding box
	int ndim = varInfo->ndim;

    if (adiosQuery->sel == NULL) {
         static uint64_t ZERO[32] = { 0 };
		 destcount = varInfo->dims;   deststart = ZERO;

         if (varInfo->blockinfo == NULL) {
        	 adios_read_set_data_view(adiosQuery->file, LOGICAL_DATA_VIEW);
             common_read_inq_var_blockinfo(adiosQuery->file, varInfo);
         }

		adios_read_set_data_view(adiosQuery->file, PHYSICAL_DATA_VIEW);
		ADIOS_VARTRANSFORM *ti = adios_inq_var_transform(adiosQuery->file, varInfo);

		int totalPG = varInfo->nblocks[timeStep];
		int blockId, j;

		// global block id
		int startBlkId = 0;
		int blkIter;
		for (blkIter = 0; blkIter < timeStep; blkIter++) {
			startBlkId += varInfo->nblocks[blkIter];
		}

		for (j = startBlkId; j < startBlkId + totalPG; j++) {
			// false: PG selection box is intersecting with variable's selection box
			// true : PG selection box is fully contained within variable's selection box
                        // isPGCovered is alwasy true for selection is NULL
			bool isPGCovered = true;

			srcstart = varInfo->blockinfo[j].start;
			srccount = varInfo->blockinfo[j].count;

			if (ti->transform_type == adios_get_transform_type_by_uid("ncsu-alacrity")) {
				proc_write_block(j,isPGCovered,ti, adiosQuery,startStep,estimate,&alacQuery,lb,hb
						,srcstart, srccount, deststart, destcount,alacResultBitmap	,Corder);
			}else {
				char * blockData  = NULL;
				uint64_t totalElm = 1;
				int t = 0;
				for(t=0; t < ndim; t++){
					totalElm *= srccount[t];
				}
				readBlockData(j, adiosQuery, startStep, varInfo,totalElm, (void**)&blockData );
				literallyCheckData(blockData, totalElm, varInfo->type,adiosQuery,hb, lb
						, srcstart, srccount, deststart, destcount, ndim, isPGCovered
						, alacResultBitmap,Corder);
			}
		}

	}else if (adiosQuery->sel->type == ADIOS_SELECTION_BOUNDINGBOX) {

		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(adiosQuery->sel->u.bb);
		destcount = bb->count;   deststart = bb->start;
#ifdef RIDBUG
		printf("Output Bounding box, start: [");
		int t = 0;
		for(t  = 0 ; t < ndim ; t++){
			printf("%"PRIu64", ", deststart[t]);
		}
		printf("], count: [");
		for(t  = 0 ; t < ndim ; t++){
			printf("%"PRIu64", ", destcount[t]);
		}
		printf("]\n");

#endif


		adios_read_set_data_view(adiosQuery->file, PHYSICAL_DATA_VIEW);
		ADIOS_VARTRANSFORM *ti = adios_inq_var_transform(adiosQuery->file, varInfo);

#ifdef BREAKDOWN
	findPGStart = dclock();
#endif
		ADIOS_PG_INTERSECTIONS* intersectedPGs = adios_find_intersecting_pgs( adiosQuery->file, varInfo->varid, adiosQuery->sel, timeStep, numStep);

#ifdef BREAKDOWN
	findPGTotal += dclock() - findPGStart;
	numTouchedPGs = intersectedPGs->npg;
#endif
		int totalPG = intersectedPGs->npg;
		int blockId, j;

		ADIOS_PG_INTERSECTION *  PGs = intersectedPGs->intersections;
		for (j = 0; j < totalPG; j++) {
			ADIOS_PG_INTERSECTION pg = PGs[j];
			ADIOS_SELECTION * interSelBox = pg.intersection_sel;
			// false: PG selection box is intersecting with variable's selection box
			// true : PG selection box is fully contained within variable's selection box
			bool isPGCovered = false;
			ADIOS_SELECTION * pgSelBox = pg.pg_bounds_sel;
			assert(pgSelBox->type == ADIOS_SELECTION_BOUNDINGBOX );
			assert(pgSelBox->type == interSelBox->type );
			const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *pgBB = &(pgSelBox->u.bb);
			const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *interBB = &(interSelBox->u.bb);
			isPGCovered = boxEqual(pgBB, interBB);

			srcstart = pgBB->start;
			srccount = pgBB->count;
			// blockId = pg.blockidx; // pg.blockidx is relative block id           //	blockId = pg.blockidx_in_timestep ;
			if (ti->transform_type == adios_get_transform_type_by_uid("ncsu-alacrity")) {
#ifdef BREAKDOWN
	  procStart = dclock();
#endif
				proc_write_block(pg.blockidx,isPGCovered,ti, adiosQuery,startStep,estimate,&alacQuery,lb,hb
						,srcstart, srccount, deststart, destcount,alacResultBitmap,Corder	);
#ifdef BREAKDOWN
	  procTotal = procTotal + (dclock() - procStart);
#endif

			}else {
				char * blockData  = NULL;
				uint64_t totalElm = 1;
				int t = 0;
				for(t=0; t < ndim; t++){
					totalElm *= pgBB->count[t];
				}
				readBlockData(pg.blockidx, adiosQuery, startStep, varInfo,totalElm, (void**)&blockData );
				literallyCheckData(blockData, totalElm, varInfo->type,adiosQuery,hb, lb
						, srcstart, srccount, deststart, destcount, ndim, isPGCovered
						, alacResultBitmap,Corder);
			}
		}
		adios_free_pg_intersections(&intersectedPGs);

	} else if (adiosQuery->sel->type == ADIOS_SELECTION_WRITEBLOCK){

		const ADIOS_SELECTION_WRITEBLOCK_STRUCT *writeBlock = &(adiosQuery->sel->u.block);
		int blockId= writeBlock->index; //relative block id
		int globalBlockId = query_utils_getGlobalWriteBlockId(blockId, startStep, varInfo);
		common_read_inq_var_blockinfo(adiosQuery->file, varInfo);
		ADIOS_VARBLOCK block = varInfo->blockinfo[globalBlockId];
		// since user supplies the query with block id, in this case, the start and destination(querying) bounding box are the block itself
		srcstart = block.start; srccount = block.count;
		deststart = block.start;  destcount = block.count;
		bool isPGCovered = true; // because they are same bounding boxes, they are fully contained

		if (ti->transform_type == adios_get_transform_type_by_uid("ncsu-alacrity")){
			adios_read_set_data_view(adiosQuery->file, PHYSICAL_DATA_VIEW); // switch to the transform view,
			ADIOS_VARTRANSFORM *ti = adios_inq_var_transform(adiosQuery->file, varInfo); // this func. will fill the blockinfo field

			proc_write_block(globalBlockId,isPGCovered,ti, adiosQuery,startStep,estimate,&alacQuery,lb,hb
					,srcstart, srccount, deststart, destcount,alacResultBitmap,Corder	);
		}else {
			char * blockData  = NULL;
			readBlockData(globalBlockId, adiosQuery, startStep, varInfo,adiosQuery->rawDataSize, (void**) &blockData );
			literallyCheckData(blockData, adiosQuery->rawDataSize,  varInfo->type,adiosQuery,hb, lb
					, srcstart, srccount, deststart, destcount, ndim, isPGCovered
					, alacResultBitmap,Corder);
		}

	} else if (adiosQuery->sel->type == ADIOS_SELECTION_POINTS){
		// TODO: at this point, this type of querying is took careful by the common query layer
	} else {
		printf("not supported selection typed in alacrity \n");
		exit(EXIT_FAILURE);
	}


	if (adiosQuery->predicateOp == ADIOS_NE){ // we flip the bits, since we treat it "=" before
		uint64_t i  = 0;
		for (; i < alacResultBitmap->length; i++) {
			alacResultBitmap->bits[i]  = ~(alacResultBitmap->bits[i]);
		}
	}
	// NOTE: this is for correctness of reading info later. We switch back to ensure the end-use are not affected
	adios_read_set_data_view(adiosQuery->file, LOGICAL_DATA_VIEW);

#ifdef BREAKDOWN
	printf("Preparation time: %f \n", preparationTime);
	printf("Total # of touched PGs: %d \n", numTouchedPGs);
	printf("Sum touched bin # in all touched PGs: %"PRIu32"\n", numTouchedBins);
	printf("Find total # of PGs touched by the bounding box in the query: %f \n", findPGTotal);
	printf("One PG processing (Proc write block) time : %f \n", procTotal);
	printf("The following time should approximately sum up to one PG processing time \n");
	printf("===>Total time of reading alacrity metadata: %f \n", metaTotal);
	printf("===>Total time of reading partition meta in Alacrity metadata: %f \n", alacPartitionMetaTotal);
	printf("===>Total index read time : %f \n", idxTotal);
	printf("===>Total index read size (in KB): %"PRIu64"\n", idxReadSize);
	printf("===>Total low-order bytes read time : %f \n", dataTotal);
	printf("===>Total low-order bytes read size (in KB): %"PRIu64"\n", dataReadSize);
	printf("===>Candidate check total time: %f (below times should sum up to this)\n", candidateCheckTotal);
	printf("======>The number of total elements/Rids examined in Candidate check: %"PRIu64"\n", ckTotalElm);
	printf("======>Rids Decoded (if compressed) in Candidate check %f \n", ckRids);
	printf("======>Reconstitute data in Candidate check %f \n", ckReconstitute);
	printf("======>Spatial check (whether in the BB) in Candidate check %f \n", ckSpatialCheck);
	printf("======>Set bitmaps after pass the spatial check in Candidate check %f \n", ckSetBitMap);
	printf("===>Decode compressed RIDs to bitmap (optional, only applied to compressed indexes): %f %f %f %f %f \n", decodeTotal, wholeDecodeTotal,  pgCTotal , pgPartialTotal, decodeTotal - wholeDecodeTotal - pgCTotal - pgPartialTotal);
	printf("===>Set every RID in the bitmap (optional, only applied to uncompressed indexes): %f \n", setRidTotal);

#endif
	return alacResultBitmap;
}

static void dump_bitmap_word(FILE *stream, uint64_t word, int nbits) {
	int i;
	char bitstr[65];
	for (i = 0; i < nbits; ++i) {
		const int bit = (word >> i) & 1;
		bitstr[i] = bit ? '1' : '0';
	}
	bitstr[nbits] = 0;
	fprintf(stream, "%s", bitstr);

}

static void dump_bitmap(FILE *stream, ADIOS_ALAC_BITMAP *bitmap) {
	uint64_t i;
	for (i = 0; i < bitmap->length - 1; ++i)
		dump_bitmap_word(stream, bitmap->bits[i], 64);

	if (bitmap->realElmSize > 0) {
		// Dump only the valid bits of the last word
		dump_bitmap_word(
			stream,								// output stream
			bitmap->bits[bitmap->length - 1],	// the last word
			(bitmap->realElmSize - 1) % 64 + 1	// the number of valid bits in the last word
		);
	}
}

/*
 * This is an internal function processing the expression tree
 */
ADIOS_ALAC_BITMAP * adios_alac_process(ADIOS_QUERY* q, int timestep,
		bool estimate) {

	ADIOS_ALAC_BITMAP *result;
	//LEAF NODE
	ADIOS_ALAC_BITMAP * rbitmap, *lbitmap;
	if (q ->left == NULL && q->right == NULL) {
		result = adios_alac_uniengine(q, timestep, estimate);

#ifdef ADIOS_ALAC_QUERY_DEBUG
		fprintf(stderr, "constraint %s = ", q->condition);
		dump_bitmap(stderr, result);
		fprintf(stderr, "\n");
#endif
		return result;
	}

	if (q->left) {
		lbitmap = adios_alac_process((ADIOS_QUERY*) q->left, timestep, estimate);
		result = q->left;
	}

	if (q->right) {
		rbitmap = adios_alac_process((ADIOS_QUERY*) q->right, timestep, estimate);
		result = q->right;
	}

	if (q->left && q->right) {
	    result = adios_alac_bitsOp(lbitmap, rbitmap, q->combineOp);
	}

#ifdef ADIOS_ALAC_QUERY_DEBUG
	fprintf(stderr, "op %s = ", (q->combineOp == ADIOS_QUERY_OP_AND ? "AND" : "OR"));
	dump_bitmap(stderr, result);
	fprintf(stderr, "\n");
#endif
	return result;
}



/*
 * it was previously adios_query_alac_retrieval_pointsNd which returns coordinates within query bounding box,
 * now, for the minmax branch, it returns offsets within the query bounding box
 */
void adios_query_alac_retrieval_offset( ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size
		, ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb , uint64_t *points /*OUT*/, int Corder ){
	int d = 0, td;
	uint64_t start_pos = 0;
	if (!isLastConvRidInit(b))
		start_pos =	b->lastConvRid  / 64;

	uint64_t * p_bitmap = b->bits;
	uint64_t pidx = 0, off = start_pos, retrieveCount = 0;
	uint64_t reconstct_rid;
	uint64_t coordinates[MAX_DIMS];
#ifdef RIDBUG
	printf("reconstrcuted RID relative to timestep : ");
#endif
	while (off < b->length ){
		uint16_t * temp = (uint16_t *) &(p_bitmap[off]); // 2 bytes (unsigned short int)  = 16 bits
		uint64_t offset_long_int = off * 64; // original index offset ; // 4 bytes (unsigned long int )= 64 bit
		uint64_t offset;
		int j , m;
		for (j = 0; j < 4; j++) {
			offset = offset_long_int + j * 16; // here, 16 is used because temp is 16bits (unsigned short int) pointer
			// set_bit_count for each 2 bytes, the number of 1
			/*
			 * *******|               64 bits                 | => final_result_bitmap []
			 * *******| 16 bits | 16 bits | 16 bits | 16 bits | => temp[]
			 */

			for (m = 0; m < set_bit_count[temp[j]]  ; m++) {
				reconstct_rid = offset+ set_bit_position[temp[j]][m];
#ifdef RIDBUG
				printf("%"PRIu64", ", reconstct_rid);
#endif
				if (!isLastConvRidInit(b)) {
					if (reconstct_rid > b->lastConvRid) { // skip the RIDs in the 16-bits part
						points[pidx++] = reconstct_rid; // this is the change, reconstct_rid is offset
						/*ridToCoordinates(bb->ndim, Corder, reconstct_rid, bb->count, coordinates );
						for (d = 0; d < bb->ndim; d ++){
							points[pidx++] = bb->start[d] + coordinates[d];
						}*/
						retrieveCount++;
					}
				}else { // lastConvRid == realElmSize+1 represents the initial state
					// in which the RID recovering has not started yet
					points[pidx++] = reconstct_rid; // this is the change, reconstct_rid is offset

					/*ridToCoordinates(bb->ndim, Corder, reconstct_rid, bb->count, coordinates);
					for (d = 0; d < bb->ndim; d ++){
						points[pidx++] = bb->start[d] + coordinates[d];
					}*/
					retrieveCount++;
				}

				if (retrieveCount == retrieval_size){
					b->lastConvRid = reconstct_rid; 					 // updated the status

#ifdef RIDBUG
		printf("\n");
#endif
					return ;
				}
			}
		}
		off ++;
	}

#ifdef RIDBUG
		printf("\n");
#endif
}


//void adios_query_alac_init() {} // not used

void adios_query_alac_retrieval_pointsNd( ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size
		, ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb , uint64_t *points /*OUT*/, int Corder ){
	int d = 0, td;
	uint64_t start_pos = 0;
	if (!isLastConvRidInit(b))
		start_pos =	b->lastConvRid  / 64;

	uint64_t * p_bitmap = b->bits;
	uint64_t pidx = 0, off = start_pos, retrieveCount = 0;
	uint64_t reconstct_rid;
	uint64_t coordinates[MAX_DIMS];
#ifdef RIDBUG
	printf("reconstrcuted RID relative to timestep : ");
#endif
	while (off < b->length ){
		uint16_t * temp = (uint16_t *) &(p_bitmap[off]); // 2 bytes (unsigned short int)  = 16 bits
		uint64_t offset_long_int = off * 64; // original index offset ; // 4 bytes (unsigned long int )= 64 bit
		uint64_t offset;
		int j , m;
		for (j = 0; j < 4; j++) {
			offset = offset_long_int + j * 16; // here, 16 is used because temp is 16bits (unsigned short int) pointer
			// set_bit_count for each 2 bytes, the number of 1
			/*
			 * *******|               64 bits                 | => final_result_bitmap []
			 * *******| 16 bits | 16 bits | 16 bits | 16 bits | => temp[]
			 */

			for (m = 0; m < set_bit_count[temp[j]]  ; m++) {
				reconstct_rid = offset+ set_bit_position[temp[j]][m];
#ifdef RIDBUG
				printf("%"PRIu64", ", reconstct_rid);
#endif
				if (!isLastConvRidInit(b)) {
					if (reconstct_rid > b->lastConvRid) { // skip the RIDs in the 16-bits part
						ridToCoordinates(bb->ndim, Corder, reconstct_rid, bb->count, coordinates );
						for (d = 0; d < bb->ndim; d ++){
							points[pidx++] = bb->start[d] + coordinates[d];
						}
						retrieveCount++;
					}
				}else { // lastConvRid == realElmSize+1 represents the initial state
					// in which the RID recovering has not started yet
					ridToCoordinates(bb->ndim, Corder, reconstct_rid, bb->count, coordinates);
					for (d = 0; d < bb->ndim; d ++){
						points[pidx++] = bb->start[d] + coordinates[d];
					}
					retrieveCount++;
				}

				if (retrieveCount == retrieval_size){
					b->lastConvRid = reconstct_rid; 					 // updated the status

#ifdef RIDBUG
		printf("\n");
#endif
					return ;
				}
			}
		}
		off ++;
	}

#ifdef RIDBUG
		printf("\n");
#endif
}

/*
 * it was adios_query_build_results_boundingbox, i copy it for the minmax branch change, which returns the offets within bounding box, rather than the coordinates
 */
static ADIOS_SELECTION * adios_query_build_offsets_boundingbox(ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size,
                                                               ADIOS_SELECTION *box, int Corder)
{

    ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(box->u.bb);
    int ndim = 1 ; // offset is 1 dimensional data
    const uint64_t dataSize = retrieval_size * (ndim);
    uint64_t *points = (uint64_t *)(malloc(dataSize * sizeof(uint64_t)));
    adios_query_alac_retrieval_offset(b,retrieval_size, bb, points,Corder);

    return a2sel_points(ndim, retrieval_size, points, box, 1);
}


static ADIOS_SELECTION * adios_query_build_results_boundingbox(ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size,
                                                               ADIOS_SELECTION *box, int Corder)
{
    ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(box->u.bb);
	const uint64_t dataSize = retrieval_size * (bb->ndim);
	uint64_t *points = (uint64_t *)(malloc(dataSize * sizeof(uint64_t)));
	adios_query_alac_retrieval_pointsNd(b,retrieval_size, bb, points,Corder);

	return a2sel_points(bb->ndim, retrieval_size, points, box, 1);
}


void adios_query_alac_build_results(
		uint64_t retrieval_size, ADIOS_SELECTION* outputBoundry, ADIOS_ALAC_BITMAP *b,
		ADIOS_VARINFO *varinfo,  ADIOS_QUERY_RESULT* adiosQueryResult /*ADIOS_SELECTION ** queryResult*/ , int Corder)
{

	//last bounding box / points supplied by user
	switch (outputBoundry->type) {
	case ADIOS_SELECTION_BOUNDINGBOX: {
		ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(outputBoundry->u.bb);

#ifdef ALACRITY_EVAL_1D_POINTS
		adiosQueryResult->selections = adios_query_build_offsets_boundingbox(b, retrieval_size, outputBoundry, Corder);
#else
		adiosQueryResult->selections = adios_query_build_results_boundingbox(b, retrieval_size, outputBoundry, Corder);
#endif
        adiosQueryResult ->nselections = 1 ;
        adiosQueryResult -> npoints = retrieval_size ;
        adiosQueryResult ->method_used = ADIOS_QUERY_METHOD_ALACRITY;
		break;
	}
	case ADIOS_SELECTION_WRITEBLOCK : {
		adios_error(err_unspecified,
					"Detected writeblock selection as query output selection internally, "
					"but the writeblock selection should have already been mapped to a bounding box! (at %s:%s)",
					__FILE__, __LINE__);
		break;
	}
	case ADIOS_SELECTION_POINTS: {
		const ADIOS_SELECTION_POINTS_STRUCT *points =
				&(outputBoundry->u.points);
		uint64_t arraySize = retrieval_size * (points->ndim);
		uint64_t* pointArray =
				(uint64_t*) (malloc(arraySize * sizeof(uint64_t)));
		//TODO:
		break;
	}
	default:
		printf("Error: Type of selection is not supported!");
	}
}

int adios_query_alac_can_evaluate(ADIOS_QUERY* q)
{
	int is_alac = 0;
	if (!q->left && !q->right) {
		// If this is a query leaf node, we support ALACRITY query iff
		// that variable is transformed using the ALACRITY indexing method
		ADIOS_VARTRANSFORM *vartrans = adios_inq_var_transform(q->file, q->varinfo);
		is_alac = (vartrans->transform_type == adios_get_transform_type_by_uid("ncsu-alacrity"));
		adios_free_var_transform(vartrans);
	} else {
		// Else, this is an internal node, and we support ALACRITY query if
		// any descendent node supports ALACRITY (since ALACRITY query is
		// capable of processing indexed and non-indexed variables in the
		// same query)
		if (q->left)
			is_alac |= adios_query_alac_can_evaluate((ADIOS_QUERY *)q->left);
		if (q->right)
			is_alac |= adios_query_alac_can_evaluate((ADIOS_QUERY *)q->right);
	}
    return is_alac;
}



int64_t adios_query_alac_estimate(ADIOS_QUERY* q, int timestep) {
	ADIOS_ALAC_BITMAP* b = adios_alac_process(q, timestep, true);
	return calSetBitsNum(b);
}

/* memory format:
 * | length | numSetBits | realNumSize | lastConvRid  |bits ....
 */
uint64_t * convertALACBitmapTomemstream( ADIOS_ALAC_BITMAP * b ){
	int metaLen = 4; //----- 1 /*b->length*/ + 1 /*b->numSetBits*/ + 1 /*b->realNumSize*/ + 1 /*b->lastConvRid*/----//
	uint64_t totalSize = b->length + metaLen;
	uint64_t * headPtr  = (uint64_t *) calloc (totalSize, sizeof(uint64_t));
	if (headPtr == NULL){
		printf("%s failed to allocat %"PRIu64 " bytes memory \n", __FUNCTION__, sizeof(uint64_t)*totalSize);
		return NULL;
	}
	uint64_t * movePtr = headPtr;
	movePtr[0] = b->length;   movePtr[1] = b->numSetBits; movePtr[2] = b->realElmSize; movePtr[3] = b->lastConvRid;
	memcpy(movePtr+ metaLen, b->bits, sizeof(uint64_t) * b->length);
	return headPtr;
}

/*
 * do memory copy for the bits, so that the `mem` could be free later
 */
void convertMemstreamToALACBitmap( void *mem , ADIOS_ALAC_BITMAP * bout /*OUT*/){

	uint64_t * ptr  = (uint64_t *) mem;
	bout->length = ptr[0];
	bout->numSetBits = ptr[1];
	bout->realElmSize = ptr[2];
	bout->lastConvRid = ptr[3];
	bout->bits = (uint64_t *) malloc(sizeof(uint64_t)*(bout->length));
	memcpy(bout->bits, ptr+4, sizeof(uint64_t)*(bout->length));
}

int adios_query_alac_evaluate(ADIOS_QUERY* q,
                   int timestep,
			       uint64_t batchSize, // limited by maxResult
			       ADIOS_SELECTION* outputBoundry,
				   //ADIOS_SELECTION** result
				   	ADIOS_QUERY_RESULT* adiosQueryResult)
{
	double alacStart = 0, bimapToCoordinates= 0;

#ifdef BREAKDOWN
	alacStart = dclock();
#endif
	if (!isInitialized){ // if this is the very first time of calling the queries, we initialize the lookup tables
		init_lookup();
		isInitialized= 1;
	}
	const int absoluteTimestep = adios_get_actual_timestep(q, timestep);

	ADIOS_ALAC_BITMAP* b;
	if (q->onTimeStep != absoluteTimestep) { // if this is the first call to evaluate the query for a new timestep
		b = adios_alac_process(q, timestep, q->estimate);
		initLastConvRid(b);
		q->maxResultsDesired =  calSetBitsNum(b);
		q->resultsReadSoFar = 0;
		q->onTimeStep = absoluteTimestep;
	} else { //convert void* _internal to ADIOS_ALAC_BITMAP
		b = (ADIOS_ALAC_BITMAP*)malloc(sizeof(ADIOS_ALAC_BITMAP));
		convertMemstreamToALACBitmap(q->queryInternal, b);
	}
	uint64_t retrievalSize = q->maxResultsDesired - q->resultsReadSoFar;
	if (retrievalSize <= 0) {
		//(*queryResult) = NULL;
		adiosQueryResult ->status = ADIOS_QUERY_NO_MORE_RESULTS;
		FreeALACBITMAP(b);
		//q->onTimeStep = NO_EVAL_BEFORE;
		//printf(":: ==> no more results to fetch\n");
		return 0;
	}
	if (retrievalSize > batchSize) {
		retrievalSize = batchSize;
	}

	const int Corder = !futils_is_called_from_fortran(); // Use the dimension order of the caller; the common read layer will also mimic this order
#ifdef BREAKDOWN
	bimapToCoordinates = dclock();
#endif

	adios_query_alac_build_results(retrievalSize, outputBoundry, b, q->varinfo, adiosQueryResult, Corder);

#ifdef BREAKDOWN
	printf("Time of converting the final resultant bitmap to coordinates: %f \n", dclock() - bimapToCoordinates );
#endif


	// b->lastConvRid is updated in the above func., so the bitmap serializing function has to wait until the above function is finished
	q->queryInternal = convertALACBitmapTomemstream(b);

	/*if (q->_maxResultDesired >= 0) {
		FREE(b->bits); // these data is copied to q->_queryInternal
	}
	FREE(b); // NOTE: only free the structure*/
	FreeALACBITMAP(b);
	q->resultsReadSoFar += retrievalSize;

#ifdef BREAKDOWN
	printf("Total Time [alac plugin + adios] : %f \n", dclock() - alacStart);
#endif

	int moreResults =  q->resultsReadSoFar < q->maxResultsDesired;
	adiosQueryResult ->status = moreResults == 1 ? ADIOS_QUERY_HAS_MORE_RESULTS : ADIOS_QUERY_NO_MORE_RESULTS;
	return moreResults ;
}

int adios_query_alac_free_one_node(ADIOS_QUERY* query){
	if (query == NULL) {
		return 0;
	}

	// Currently, no ALAC-specific information to free...
	// Note: the common layer frees non-plugin-specific stuff
	// like query->predicateValue and query->dataSlice

	return 1;
}

int adios_query_alac_free(ADIOS_QUERY* query) {

	// free the tree in a bottom-to-up manner
	if (query->left == NULL && query->right == NULL) {
		return adios_query_alac_free_one_node(query);
	}else if  (query->right){
		return adios_query_alac_free(query->right);
	}else if (query->left) {
		return adios_query_alac_free(query->left);
	}

	return 1;
}

void adios_query_alac_finalize() { /* there is nothing to finalize */ }


#endif
