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
#include "public/adios_read_ext.h"
#include "public/adios_query.h"
#include "common_query.h"
//#include <alacrity.h>

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

int coordinateConversionWithCheck(uint64_t * coordinates, const  int dim
		, const  uint64_t *srcstart, const  uint64_t *deststart, const  uint64_t *destend);

bool ridConversionWithCheck(rid_t rid/*relative to local src selectoin*/
		, uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int dim, rid_t *relativeRid  );

uint64_t ridConversionWithoutCheck(uint64_t rid/*relative to local src selectoin*/,
		uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int dim);

void create_lookup(unsigned char set_bit_count[],
		unsigned char set_bit_position[][16]);


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

void setRidToBits(bool isPGCovered ,uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int ndim, rid_t * idx, uint64_t totalRids,
		ADIOS_ALAC_BITMAP *alacResultBitmap /*IN&OUT*/);


void adios_query_alac_retrieval_points2d( ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size
		, ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb , uint64_t *points /*OUT*/ );

void adios_query_alac_retrieval_points3d( ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size
		, ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb , uint64_t *points /*OUT*/ );


void proc_write_block(int blockId, bool isPGCovered, ADIOS_VARTRANSFORM *ti, ADIOS_QUERY * adiosQuery, int startStep, bool estimate
		, ALUnivariateQuery * alacQuery , double lb , double hb
		,uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount
		, ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/ );


inline void setBitsinBitMap(rid_t rid, ADIOS_ALAC_BITMAP * alacResultBitmap){

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

unsigned char set_bit_count[65536];
unsigned char set_bit_position[65536][16];

/**** END -- Funcs. that are internal funcs. ********/

void initLastConvRid (ADIOS_ALAC_BITMAP *b){
	b->lastConvRid = b->realElmSize + 1;
}

bool isLastConvRidInit(const ADIOS_ALAC_BITMAP *b){
	return (b->lastConvRid == (b->realElmSize +1));
}




/*since we are set data view is physical view, the element size is 1 byte */

void readTransformedElms(ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, int blockId, uint64_t start_elem, uint64_t num_elems, int is_timestep_relative, void * outputData/*out*/){
	ADIOS_SELECTION *sel = adios_selection_writeblock_bounded(blockId, start_elem, num_elems, is_timestep_relative);
	adios_schedule_read_byid(fp, sel, vi->varid, startStep, numStep, outputData);
	adios_perform_reads(fp, 1);
	// adios_selection_writeblock_bounded internally malloc data for adios_selection
	// so I need to free it before the next usage
	adios_selection_delete(sel);

}

void readBlockData(int blockId , ADIOS_QUERY * adiosQuery, int startStep,
		ADIOS_VARINFO * varInfo, uint64_t dataElmNum, void ** data){
	adios_read_set_data_view(adiosQuery->_f, LOGICAL_DATA_VIEW); // switch to the transform view,
	int dataElmSize = adios_type_size(varInfo->type, NULL); // data element size, in bytes
	char * blockData = (char*) (*data);
	blockData = (char *) malloc(sizeof(char) * dataElmSize * dataElmNum);
	ADIOS_SELECTION *sel = adios_selection_writeblock_bounded(blockId, 0, dataElmNum, 1); // entire PG selection
	adios_schedule_read_byid(adiosQuery->_f, sel, varInfo->varid, startStep, 1, blockData);
}
void readIndexData(int blockId, uint64_t offsetSize /*in bytes*/
		,uint64_t length /*in bytes*/, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, void *idxBytes /*OUT*/){

	readTransformedElms(fp, vi, startStep, numStep, blockId, offsetSize, length, 0, idxBytes);
}

void readLowOrderBytes(int blockId, uint64_t offsetSize /*in bytes*/
		,uint64_t length /*in bytes*/, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, void *idxBytes /*OUT*/){
	readTransformedElms(fp, vi, startStep, numStep, blockId, offsetSize,length, 0, idxBytes);
}

void readPartitionMeta( int blockId, uint64_t metaSize, ADIOS_FILE* fp,ADIOS_VARINFO* vi
		, int startStep, int numStep
		, ALMetadata *pm /*OUT*/){
	memstream_t ms = memstreamInitReturn(malloc(metaSize));
	uint64_t metaStartPos= 0;
	readTransformedElms(fp, vi, startStep,numStep,blockId, metaStartPos, metaSize, 0, ms.buf);
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


void create_lookup(unsigned char set_bit_count[],
		unsigned char set_bit_position[][16]) {
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
			op1->bits[i] ^= op2->bits[i];
		}

		FreeALACBITMAP(op2);
	//		free(op2->bits);
//		op1->elmSize = calSetBitsNum(op1);

	} else {
		printf("Operator[%d] is not surpported now \n ", operator);
	}
	return op1;
}


int coordinateConversionWithCheck(uint64_t * coordinates, const  int dim, const  uint64_t *srcstart, const  uint64_t *deststart, const  uint64_t *destend){

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
 */
bool ridConversionWithCheck(rid_t rid/*relative to local src selectoin*/, uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int dim, rid_t *relativeRid  ){

	int i = 0;
	*relativeRid = 0;
	if (dim == 3) {
		uint64_t coordinates[3]= {0}, destend[3]={0};
		for(i = 0; i < dim; i ++){
				destend[i] = deststart[i] + destcount[i] -1;
		}
		coordinates[0] = rid / (srccount[1] * srccount[2]) ;
		coordinates[1] = (rid % (srccount[1] * srccount[2])) / srccount[2];
		coordinates[2] = (rid % (srccount[1] * srccount[2])) % srccount[2] ;

		if (coordinateConversionWithCheck(coordinates, dim, srcstart, deststart, destend) < 0){
			return false;
		}

		*relativeRid = coordinates[2] + coordinates[1] * destcount[2] + coordinates[0]* destcount[1] * destcount[2];

	}

	if (dim == 2){
		uint64_t coordinates[2]= {0}, destend[2]={0};
		for(i = 0; i < dim; i ++){
			destend[i] = deststart[i] + destcount[i] -1;
		}

		coordinates[0] = rid / (srccount[1]);
		coordinates[1] = rid % (srccount[1] );
		if (coordinateConversionWithCheck(coordinates, dim, srcstart, deststart, destend) < 0){
			return false;
		}

		*relativeRid = coordinates[1] + coordinates[0] * destcount[1] ;

	}

	return true;
}


uint64_t ridConversionWithoutCheck(uint64_t rid/*relative to local src selectoin*/,
		uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount,
		int dim){

	uint64_t relativeRid = 0;
	int * coordinates = (int *) malloc(sizeof(int) * dim); // coordinate of current PG
	if (dim == 3) {
		coordinates[0] = rid / (srccount[1] * srccount[2]);
		coordinates[1] = (rid % (srccount[1] * srccount[2])) / srccount[2];
		coordinates[2] = (rid % (srccount[1] * srccount[2])) % srccount[2] ;
		relativeRid = coordinates[2] + coordinates[1] * destcount[2] + coordinates[0]* destcount[1] * destcount[2];
	}

	if (dim == 2){
		coordinates[0] = rid / (srccount[1]);
		coordinates[1] = rid % (srccount[1] );
		relativeRid = coordinates[1] + coordinates[0] * destcount[1] ;
	}

	free(coordinates);
	return relativeRid;
}

/*
 * usage: this MACRO tightly depends on the usage context,
 * which requires the context declare the exact variables
 */
#define CHECK_ELEMENT(code) { \
		for(el= 0; el < totalElm; el ++){              \
			if ((code)){                               \
				rid = decodeRids[el];                  \
				if (ridConversionWithCheck(rid, srcstart,srccount, deststart,destcount, dim, &newRid)){  \
					setBitsinBitMap(newRid, alacResultBitmap);       \
					alacResultBitmap->numSetBits ++;   \
				}                                      \
			}                                          \
		}                                              \
}

/*
 * NOTE : relies on  a lots of variables:
 * isPGCovered, el, srcstart, srcount, deststart, destcount, ndim, newRid, totalElm, op
 */
#define CHECK_NODECODE_ELEMENT(code) {                         \
		if (isPGCovered){                             \
			for(el= 0; el < totalElm; el ++){         \
				if (code) {                           \
					newRid = ridConversionWithoutCheck(el, srcstart,srccount, deststart,destcount, ndim);   \
					setBitsinBitMap(newRid, alacResultBitmap);    \
					alacResultBitmap->numSetBits ++;              \
				}                                                 \
			}                                                     \
		}else {                                                   \
			for(el= 0; el < totalElm; el ++){                     \
				if (code) {                                       \
					if (ridConversionWithCheck(el,	srcstart,srccount, deststart,destcount, ndim, &newRid)){  \
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
 * TODO : GOT compile errors using following codes
 */
/*#define GENERIC_DATA_CONV(data, FUNC)  {                                             \
	switch(dataType){                                                                \
		case adios_unsigned_byte:                                                    \
			CHECK_GENERIC_DATA((unsigned char  *)data, (FUNC) );           \
			break;                                                                   \
		case adios_byte:                                                             \
			CHECK_GENERIC_DATA((signed char  *)data, (FUNC) );             \
			break;                                                                   \
		case adios_string:                                                           \
			CHECK_GENERIC_DATA((char *)data, (FUNC) );                   \
			break;                                                                   \
		case adios_unsigned_short:                                                   \
			CHECK_GENERIC_DATA((unsigned short  *)data, (FUNC) );          \
			break;                                                                   \
		case adios_short:                                                            \
			CHECK_GENERIC_DATA((signed short  *)data, (FUNC) );            \
			break;                                                                   \
		case adios_unsigned_integer:                                                 \
			CHECK_GENERIC_DATA((unsigned int *)data, (FUNC) );            \
			break;                                                                   \
		case adios_integer:                                                          \
			CHECK_GENERIC_DATA((signed int *)data, (FUNC) );              \
			break;                                                                   \
		case adios_unsigned_long:                                                    \
			CHECK_GENERIC_DATA((unsigned long long  *)data, (FUNC) );       \
			break;                                                                   \
		case adios_long:                                                             \
			CHECK_GENERIC_DATA((signed long long  *)data, (FUNC) );         \
			break;                                                                   \
		case adios_real:                                                             \
			CHECK_GENERIC_DATA((float *)data, (FUNC) );                  \
			break;                                                                   \
		case adios_double:                                                           \
			CHECK_GENERIC_DATA((double *)data, (FUNC) );                 \
			break;                                                                   \
	}                                                                                \
}*/


/*
 * return : 1 => yes, there is unsupported data type
 *          0 => no
 */
int checkUnsupportedDataType(enum ADIOS_DATATYPES dataType) {

	return (dataType == adios_long_double) ||
			(dataType == adios_complex)  ||
			(dataType ==  adios_double_complex) ;
}
int adios_alac_check_candidate(ALMetadata *partitionMeta, bin_id_t startBin, bin_id_t endBin , double hb, double lb
		, uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount, int dim
		, ADIOS_QUERY * adiosQuery , const char *inputCurPtr /*index bytes of entire PG*/
		, bool decoded  /*true: need decoding */ , char * lowOrderBytes /*low order bytes of from startBin to endBin*/
		, enum ADIOS_DATATYPES dataType
		,ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/){
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

	char * data = (char *) malloc(partitionMeta->elementSize*totalElm); // recovered data in bytes
	reconstituteData(partitionMeta, startBin, endBin,
										 lowOrderBytes, data);
	//Following variables are needed for the micros definition
	rid_t newRid;
	enum ADIOS_PREDICATE_MODE op = adiosQuery->_op;
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
	return 0;

}

/*
 * go through every data value in the array, check whether the value satisfies the query condition
 * return: 1 ==> error, 0 ==> no error
 */
int literallyCheckData(void *data, uint64_t totalElm , enum ADIOS_DATATYPES dataType, ADIOS_QUERY * adiosQuery, double hb, double lb
		, uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount, int ndim, bool isPGCovered
		,ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/){
	uint32_t newRid ;
	uint64_t el = 0;
	enum ADIOS_PREDICATE_MODE op = adiosQuery->_op;
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
		ADIOS_ALAC_BITMAP *alacResultBitmap /*IN&OUT*/){
	uint32_t newRid ;
	if (isPGCovered){ // if PG is fully covered by the user bounding box, we dont need to check the new location(position) of converted point
		uint64_t ni;
		for ( ni = 0; ni < totalRids; ni++) {
			rid_t rid_val = idx[ni];
			newRid = ridConversionWithoutCheck(rid_val, srcstart,srccount, deststart,destcount, ndim);
			setBitsinBitMap(newRid, alacResultBitmap);
		}
		alacResultBitmap->numSetBits += totalRids;
	}else {
		uint64_t ni;
		for ( ni = 0; ni < totalRids; ni++) {
			rid_t rid_val = idx[ni];
			if (ridConversionWithCheck(rid_val,	srcstart,srccount, deststart,destcount, ndim, &newRid)){
				setBitsinBitMap(newRid, alacResultBitmap);
				alacResultBitmap->numSetBits ++;
			}
		}
	}
}



void resolveQueryBoundary(ADIOS_QUERY *adiosQuery, double *hb, double *lb)
{
	(*hb)= DBL_MAX;
	(*lb)= DBL_MIN;
    if(adiosQuery->_op == ADIOS_LT || adiosQuery->_op == ADIOS_LTEQ){
    	*hb = atof(adiosQuery->_value);
    }else if(adiosQuery->_op == ADIOS_GT || adiosQuery->_op == ADIOS_GTEQ){
    	*lb = atof(adiosQuery->_value);
	}else if(adiosQuery->_op == ADIOS_EQ){
		//following two cases are tricky to ALACRITY
		*hb = atof(adiosQuery->_value);
		*lb = *hb;
	}else if(adiosQuery->_op == ADIOS_NE){
		//TODO: flip the bits once the evaluation is done
		*hb = atof(adiosQuery->_value);
		*lb = *hb;
	}else{
			printf("Unsupported predicate type[%d] \n", adiosQuery->_op);
	}

}

/*
 * read low order bytes from 'low_bin'  bin to 'hi_bin'
 * allocated the data buffer in this function, requiring caller to free the readData
 */
char * readLowDataAmongBins(ALMetadata *partitionMeta, bin_id_t low_bin , bin_id_t hi_bin
		, uint64_t lowDataByteStartPos,  ADIOS_QUERY * adiosQuery
		, int blockId, int startStep, int numStep){
	const uint64_t bin_read_len = ALGetDataBinOffset( partitionMeta, hi_bin) - ALGetDataBinOffset( partitionMeta, low_bin); // in bytes
	uint64_t lowDataBinOffset = lowDataByteStartPos+ ALGetDataBinOffset( partitionMeta, low_bin); /*element offset*/;
	// low order bytes from low_bin to hi_bin, ITS NOT entire low order byte
	char * readData= (char *) calloc(bin_read_len, sizeof(char));
	readLowOrderBytes(blockId,lowDataBinOffset, bin_read_len
	         , adiosQuery->_f, adiosQuery->_var, startStep, numStep, (void *)readData);
	return readData;
}


char * readIndexAmongBins(ALMetadata *partitionMeta, bin_id_t low_bin , bin_id_t hi_bin
		, uint64_t lowDataByteStartPos,  ADIOS_QUERY * adiosQuery
		, int blockId, int startStep, int numStep){
	const uint64_t bin_read_len = ALGetIndexBinOffset( partitionMeta, hi_bin) - ALGetIndexBinOffset( partitionMeta, low_bin); // in bytes
	uint64_t lowDataBinOffset = lowDataByteStartPos+ ALGetIndexBinOffset( partitionMeta, low_bin); /*element offset*/;
	// low order bytes from low_bin to hi_bin, ITS NOT entire low order byte
	char * readData= (char *) calloc(bin_read_len, sizeof(char));
	readLowOrderBytes(blockId,lowDataBinOffset, bin_read_len
	         , adiosQuery->_f, adiosQuery->_var, startStep, numStep, (void *)readData);
	return readData;
}


void proc_write_block(int blockId, bool isPGCovered, ADIOS_VARTRANSFORM *ti, ADIOS_QUERY * adiosQuery, int startStep, bool estimate
		, ALUnivariateQuery * alacQuery , double lb , double hb
		,uint64_t *srcstart, uint64_t *srccount, uint64_t *deststart, uint64_t *destcount
		, ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/ ){
	int numStep = 1; // only deal with one timestep
	uint64_t metaSize, indexSize, dataSize;
	ADIOS_VARINFO * varInfo = adiosQuery->_var;
	int ndim = varInfo->ndim;
	ADIOS_TRANSFORM_METADATA * tmetas = ti->transform_metadatas;
	ADIOS_TRANSFORM_METADATA tmeta = tmetas[blockId];
	//	assert(tmeta->length == 24);
	uint64_t * threeData = (uint64_t *) tmeta.content;
	metaSize   = threeData[0];
	indexSize  = threeData[1];
	dataSize   = threeData[2];

	//TODO: offset of each PG should be included
	printf("PG[%d] has meta size[ %" PRIu64 "], index size[ %" PRIu64 "], and data size[ %" PRIu64 "] \n",
			blockId, metaSize,  indexSize, dataSize);

	// transformed data has 1 dimension,
	// 1. load partition Metadata
	// NOTE: One ALACRITY PG Data is written in the below format:  [meta data] | [low order bytes data] | [ index data]
	ALMetadata partitionMeta;
	readPartitionMeta(blockId, metaSize,adiosQuery->_f, varInfo
					,startStep,numStep,&partitionMeta);
	const uint8_t insigbytes = insigBytesCeil(&partitionMeta);

	//2. find touched bin
	bin_id_t low_bin, hi_bin;
	_Bool are_bins_touched = findBinRange1C(&partitionMeta, alacQuery, &low_bin,
			&hi_bin);

	if (are_bins_touched) {

		//3. load index size
		uint64_t indexStartPos = metaSize + dataSize;
		char * index = readIndexAmongBins(&partitionMeta
									, low_bin,  hi_bin, indexStartPos, adiosQuery, blockId, startStep, numStep);
		char * input_index = index;
		const ALBinLayout * bl = &(partitionMeta.binLayout);
		ALIndex* indexPtr = &index;
		//TODO: distinguish the offset btw two bins for compressed and uncompressed index
		// is the offset byte-level or element-level?
		if (partitionMeta.indexMeta.indexForm == ALInvertedIndex) {
			// indexes are inverted indexes that are not compressed,  we build bitmaps for each rid;
			//element offset, instead of byte element
			uint64_t resultCount = bl->binStartOffsets[hi_bin] - bl->binStartOffsets[low_bin];

			if (estimate) {
				setRidToBits(isPGCovered, srcstart, srccount, deststart, destcount, ndim
						, (rid_t *)index, resultCount, alacResultBitmap);

			} else {
				uint64_t lowByteStartPos2 = metaSize;
				char * lowOrderBytes2  = readLowDataAmongBins(&partitionMeta
										,low_bin, hi_bin,  lowByteStartPos2, adiosQuery, blockId, startStep, numStep);
				char *lowOrderPtr2 = lowOrderBytes2; // temporary pointer
				rid_t * decodedRid = (rid_t *) index;

				// It touches at least 3 bins, so, we need to check RIDs that are in first and last bins
				if (hi_bin - low_bin > 2) {
					bin_offset_t lowBinElm = bl->binStartOffsets[low_bin + 1] - bl->binStartOffsets[low_bin];
					uint64_t hiBinElm = bl->binStartOffsets[hi_bin] - bl->binStartOffsets[hi_bin-1];
					// low boundary bin
					adios_alac_check_candidate(&partitionMeta, low_bin, low_bin+1 , hb, lb
							, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , (char*) decodedRid /*index bytes of entire PG*/
							, false  /*don't need decoding*/ , lowOrderPtr2, varInfo->type
							,alacResultBitmap /*OUT*/);
					decodedRid += lowBinElm;

					uint64_t innerElm = resultCount- lowBinElm - hiBinElm;
					setRidToBits(isPGCovered, srcstart, srccount, deststart, destcount, ndim
												, decodedRid, innerElm, alacResultBitmap);
					decodedRid  +=  innerElm;

					lowOrderPtr2 += (( bl->binStartOffsets[hi_bin-1] - bl->binStartOffsets[low_bin]) * insigbytes);
					// high boundary bin
					adios_alac_check_candidate(&partitionMeta, hi_bin-1, hi_bin , hb, lb
							, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , (char *)decodedRid
							, false , lowOrderPtr2, varInfo->type
							,alacResultBitmap /*OUT*/);

				} else { // for 1 or 2 bins touched, we need to check all RIDs
					adios_alac_check_candidate(&partitionMeta, low_bin, hi_bin  , hb, lb
							, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , (char*)decodedRid
							, false , lowOrderPtr2, varInfo->type
							,alacResultBitmap /*OUT*/);
				}

				FREE(lowOrderBytes2);

			}

		}else if (partitionMeta.indexMeta.indexForm == ALCompressedInvertedIndex) {
			const uint64_t *compBinStartOffs = partitionMeta.indexMeta.u.ciim.indexBinStartOffsets;
			uint64_t binCompressedLen;
			const char *inputCurPtr = input_index;

			if (estimate) {
				// Now compress each bin in turn
				bin_id_t bin ;
				for ( bin = low_bin; bin < hi_bin; bin++) {
					binCompressedLen = compBinStartOffs[bin + 1] - compBinStartOffs[bin];
					uint32_t decodedElm = ALDecompressRIDtoSelBox(isPGCovered , inputCurPtr, binCompressedLen
							, srcstart, srccount /*PG region dimension*/ , deststart, destcount /*region dimension of Selection box*/
							, ndim , &(alacResultBitmap->bits));
					inputCurPtr += binCompressedLen;
					alacResultBitmap->numSetBits += decodedElm;
				}

			}else{
				// element count of touched bins, which is also the count of low order data
				uint64_t lowByteStartPos = metaSize;
				char * lowOrderBytes  = readLowDataAmongBins(&partitionMeta
						, low_bin,  hi_bin, lowByteStartPos , adiosQuery, blockId, startStep, numStep);
				char *lowOrderPtr = lowOrderBytes; // temporary pointer

				// It touches at least 3 bins, so, we need to check RIDs that are in first and last bins
				if (hi_bin - low_bin > 2) {
					// low boundary bin, compressed byte offset
					binCompressedLen = compBinStartOffs[low_bin + 1] - compBinStartOffs[low_bin];

					adios_alac_check_candidate(&partitionMeta, low_bin, low_bin+1 , hb, lb
							, srcstart, srccount, deststart, destcount, ndim,  adiosQuery  , inputCurPtr /*index bytes of entire PG*/
							, true  /*need decoding*/ , lowOrderPtr /*it points to the start of `low_bin` */, varInfo->type
							,alacResultBitmap /*OUT*/);
					inputCurPtr += binCompressedLen;

					bin_id_t innerlowBin = low_bin + 1;
					bin_id_t innerHiBin = hi_bin -1;
					bin_id_t bin;
					// Now compress each bin in turn
					for ( bin = innerlowBin; bin < innerHiBin; bin++) {
						binCompressedLen = compBinStartOffs[bin + 1] - compBinStartOffs[bin];

						uint32_t decodedElm = ALDecompressRIDtoSelBox(isPGCovered
								, inputCurPtr, binCompressedLen
								, srcstart, srccount //PG region dimension
								, deststart, destcount //region dimension of Selection box
								, ndim , &(alacResultBitmap->bits));
						alacResultBitmap->numSetBits += decodedElm;
						inputCurPtr += binCompressedLen;
					}

					// high boundary bin
					binCompressedLen = compBinStartOffs[hi_bin]- compBinStartOffs[hi_bin-1];
					// point to the low order byte of one bin before hi_bin
					lowOrderPtr += (( bl->binStartOffsets[hi_bin-1] - bl->binStartOffsets[low_bin]) * insigbytes);
					adios_alac_check_candidate(&partitionMeta, hi_bin-1, hi_bin , hb, lb
							, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , inputCurPtr /*index bytes of entire PG*/
							, true  /*need decoding*/ , lowOrderPtr /*low order bytes of entire PG*/, varInfo->type
							,alacResultBitmap /*OUT*/);
					inputCurPtr += binCompressedLen;

				} else { // for 1 or 2 bins touched, we need to check all RIDs

					adios_alac_check_candidate(&partitionMeta, low_bin, hi_bin , hb, lb
							, srcstart, srccount, deststart, destcount, ndim,  adiosQuery , inputCurPtr /*index bytes of entire PG*/
							, true , lowOrderPtr, varInfo->type
							,alacResultBitmap /*OUT*/);
				}

				FREE(lowOrderBytes);
			}
		} else {
			printf("index form %d in alacrity is not supported", partitionMeta.indexMeta.indexForm);
			exit(EXIT_FAILURE);
		}
		FREE(input_index);
	}else {
		printf("there is no touched bin for constraint \n");
	}
}



/*
 * PG selection, [startPG, endPG)
 */
ADIOS_ALAC_BITMAP* adios_alac_uniengine(ADIOS_QUERY * adiosQuery, int timeStep, bool estimate) {

	if (checkUnsupportedDataType(adiosQuery->_var->type)){
		printf("unsupported data type [%d] at this point \n", adiosQuery->_var->type);
		exit(EXIT_FAILURE);
	}

	double lb , hb ;
    resolveQueryBoundary(adiosQuery, &hb, &lb); // query constraints
    printf("%s\n", adiosQuery->_condition);

	ADIOS_VARINFO * varInfo = adiosQuery->_var;

	adios_read_set_data_view(adiosQuery->_f, LOGICAL_DATA_VIEW); // switch to the transform view,
 	ADIOS_VARTRANSFORM *ti = adios_inq_var_transform(adiosQuery->_f, varInfo); // this func. will fill the blockinfo field
	int startStep = timeStep, numStep = 1;
	uint64_t totalElm = adiosQuery->_rawDataSize; // no matter bounding box or writeblock selection, the rawDataSize has been calculated in the common query layer
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
		ALQueryEngineStartUnivariateDoubleQuery(&qe, lb, hb, REGION_RETRIEVAL_INDEX_ONLY_QUERY_TYPE, &alacQuery);
	}

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
	if (adiosQuery->_sel->type == ADIOS_SELECTION_BOUNDINGBOX) {
		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(adiosQuery->_sel->u.bb);
		destcount = bb->count;   deststart = bb->start;
		adios_read_set_data_view(adiosQuery->_f, PHYSICAL_DATA_VIEW);
		ADIOS_VARTRANSFORM *ti = adios_inq_var_transform(adiosQuery->_f, varInfo);
		ADIOS_PG_INTERSECTIONS* intersectedPGs = adios_find_intersecting_pgs( adiosQuery->_f, varInfo->varid, adiosQuery->_sel, timeStep, numStep);
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
			blockId = pg.blockidx_in_timestep ;

			if (ti->transform_type == adios_get_transform_type_by_uid("ncsu-alacrity")) {
				proc_write_block(blockId,isPGCovered,ti, adiosQuery,startStep,estimate,&alacQuery,lb,hb
						,srcstart, srccount, deststart, destcount,alacResultBitmap	);
			}else {
				char * blockData  = NULL;
				uint64_t totalElm = 1;
				int t = 0;
				for(t=0; t < ndim; t++){
					totalElm *= pgBB->count[t];
				}
				readBlockData(blockId, adiosQuery, startStep, varInfo,totalElm, (void**)&blockData );
				literallyCheckData(blockData, totalElm, varInfo->type,adiosQuery,hb, lb
						, srcstart, srccount, deststart, destcount, ndim, isPGCovered
						, alacResultBitmap);
			}
		}
		adios_free_pg_intersections(&intersectedPGs);
	}
	else if (adiosQuery->_sel->type == ADIOS_SELECTION_WRITEBLOCK){
		const ADIOS_SELECTION_WRITEBLOCK_STRUCT *writeBlock = &(adiosQuery->_sel->u.block);
		int blockId= writeBlock->index;
		int globalBlockId = getGlobalWriteBlockId(blockId, startStep, varInfo);
		adios_inq_var_blockinfo(adiosQuery->_f, varInfo);
		ADIOS_VARBLOCK block = varInfo->blockinfo[globalBlockId];
		// since user supplies the query with block id, in this case, the start and destination(querying) bounding box are the block itself
		srcstart = block.start; srccount = block.count;
		deststart = block.start;  destcount = block.count;
		bool isPGCovered = true; // because they are same bounding boxes, they are fully contained

		if (ti->transform_type == adios_get_transform_type_by_uid("ncsu-alacrity")){
			adios_read_set_data_view(adiosQuery->_f, PHYSICAL_DATA_VIEW); // switch to the transform view,
			ADIOS_VARTRANSFORM *ti = adios_inq_var_transform(adiosQuery->_f, varInfo); // this func. will fill the blockinfo field
			proc_write_block(blockId,isPGCovered,ti, adiosQuery,startStep,estimate,&alacQuery,lb,hb
					,srcstart, srccount, deststart, destcount,alacResultBitmap	);
		}else {
			char * blockData  = NULL;
			readBlockData(blockId, adiosQuery, startStep, varInfo,adiosQuery->_rawDataSize, (void**) &blockData );
			literallyCheckData(blockData, adiosQuery->_rawDataSize,  varInfo->type,adiosQuery,hb, lb
					, srcstart, srccount, deststart, destcount, ndim, isPGCovered
					, alacResultBitmap);
		}

	} else if (adiosQuery->_sel->type == ADIOS_SELECTION_POINTS){
		// TODO: at this point, this type of querying is took careful by the common query layer
	} else {
		printf("not supported selection typed in alacrity \n");
		exit(EXIT_FAILURE);
	}


	if (adiosQuery->_op == ADIOS_NE){ // we flip the bits, since we treat it "=" before
		uint64_t i  = 0;
		for (; i < alacResultBitmap->length; i++) {
			alacResultBitmap->bits[i]  = ~(alacResultBitmap->bits[i]);
		}
	}
	// NOTE: this is for correctness of reading info later. We switch back to ensure the end-use are not affected
	adios_read_set_data_view(adiosQuery->_f, LOGICAL_DATA_VIEW);
	return alacResultBitmap;
}

/*
 * This is an internal function processing the expression tree
 */
ADIOS_ALAC_BITMAP * adios_alac_process(ADIOS_QUERY* q, int timestep,
		bool estimate) {

	//LEAF NODE
	ADIOS_ALAC_BITMAP * rbitmap, *lbitmap;
	if (q ->_left == NULL && q->_right == NULL) {
		return adios_alac_uniengine(q, timestep, estimate);
	}

	if (q->_left)
		lbitmap = adios_alac_process((ADIOS_QUERY*) q->_left, timestep, estimate);

	if (q->_right)
		rbitmap = adios_alac_process((ADIOS_QUERY*) q->_right, timestep, estimate);


	return adios_alac_bitsOp(lbitmap, rbitmap, q->_leftToRightOp );
}


void adios_query_alac_init_method() {}


void adios_query_alac_retrieval_points2d(
		ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size
		, ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb
		, uint64_t *points /*OUT*/ ){

	uint64_t start_pos = 0;
	if (!isLastConvRidInit(b))
		start_pos =	b->lastConvRid  / 64;

	uint64_t * p_bitmap = b->bits;
	uint64_t pidx = 0, off = start_pos, retrieveCount = 0;
	uint64_t reconstct_rid;
	while (off <= b->length ){
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
				if (!isLastConvRidInit(b)) {
					if (reconstct_rid > b->lastConvRid) { // skip the RIDs in the 16-bits part
//						printf("recovered RID %"PRIu64"\n", reconstct_rid);
						points[pidx++] = bb->start[0] + reconstct_rid / bb->count[1];
						points[pidx++] = bb->start[1] + reconstct_rid % bb->count[1];
						retrieveCount++;
					}
				}else { // lastConvRid == realElmSize+1 represents the initial state
					// in which the RID recovering has not started yet
//					printf("recovered RID %"PRIu64"\n", reconstct_rid);
					points[pidx++] = bb->start[0] + reconstct_rid / bb->count[1];
					points[pidx++] = bb->start[1] + reconstct_rid % bb->count[1];
					retrieveCount++;
				}

				if (retrieveCount == retrieval_size){
					b->lastConvRid = reconstct_rid; 					 // updated the status
					return ;
				}
			}
		}
		off ++;
	}
}

void adios_query_alac_retrieval_points3d( ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size
		, ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb , uint64_t *points /*OUT*/ ){

	uint64_t start_pos = 0;
	if (!isLastConvRidInit(b))
		start_pos =	b->lastConvRid  / 64;

	uint64_t * p_bitmap = b->bits;
	uint64_t pidx = 0, off = start_pos, retrieveCount = 0;
	uint64_t reconstct_rid;
	while (off <= b->length ){
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
				if (!isLastConvRidInit(b)) {
					if (reconstct_rid > b->lastConvRid) { // skip the RIDs in the 16-bits part
//						printf("recovered RID %"PRIu64"\n", reconstct_rid);
						points[pidx++] = bb->start[0] + reconstct_rid / (bb->count[1] * bb->count[2]) ;
						points[pidx++] = bb->start[1] + (reconstct_rid % (bb->count[1] * bb->count[2])) / bb->count[2];
						points[pidx++] = bb->start[2] + (reconstct_rid % (bb->count[1] * bb->count[2])) % bb->count[2] ;
						retrieveCount++;
					}
				}else { // lastConvRid == realElmSize+1 represents the initial state
					// in which the RID recovering has not started yet
//					printf("recovered RID %"PRIu64"\n", reconstct_rid);
					points[pidx++] = bb->start[0] + reconstct_rid / (bb->count[1] * bb->count[2]) ;
					points[pidx++] = bb->start[1] + (reconstct_rid % (bb->count[1] * bb->count[2])) / bb->count[2];
					points[pidx++] = bb->start[2] + (reconstct_rid % (bb->count[1] * bb->count[2])) % bb->count[2] ;
					retrieveCount++;
				}

				if (retrieveCount == retrieval_size){
					b->lastConvRid = reconstct_rid; 					 // updated the status
					return ;
				}
			}
		}
		off ++;
	}

}

void adios_query_alac_build_results(
		uint64_t retrieval_size, ADIOS_SELECTION* outputBoundry, ADIOS_ALAC_BITMAP *b
		, ADIOS_SELECTION ** queryResult){

	//last bounding box / points supplied by user
	switch (outputBoundry->type) {
	case ADIOS_SELECTION_BOUNDINGBOX: {
		ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(outputBoundry->u.bb);

		uint64_t dataSize = retrieval_size * (bb->ndim);
		uint64_t* points = (uint64_t*) (malloc(dataSize * sizeof(uint64_t)));
		if ( bb->ndim == 3){
			adios_query_alac_retrieval_points3d(b,retrieval_size, bb, points);
		}else if (bb -> ndim == 2){
			adios_query_alac_retrieval_points2d(b,retrieval_size, bb, points);
		}
		*queryResult = adios_selection_points(bb->ndim, retrieval_size, points);
	}
		break;
	case ADIOS_SELECTION_POINTS: {
		const ADIOS_SELECTION_POINTS_STRUCT *points =
				&(outputBoundry->u.points);
		uint64_t arraySize = retrieval_size * (points->ndim);
		uint64_t* pointArray =
				(uint64_t*) (malloc(arraySize * sizeof(uint64_t)));
		//TODO:
	}
		break;
	default:
		printf("Error: Type of selection is not supported!");
	}
}


int64_t adios_query_alac_estimate_method(ADIOS_QUERY* q) {
	ADIOS_ALAC_BITMAP* b = adios_alac_process(q, gCurrentTimeStep, true);
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

int  adios_query_alac_get_selection_method(ADIOS_QUERY* q,
			       uint64_t batchSize, // limited by maxResult
			       ADIOS_SELECTION* outputBoundry,
			       ADIOS_SELECTION** queryResult) {
	// first time, we have to evaluate it
	ADIOS_ALAC_BITMAP* b ;
	if (q->_onTimeStep == NO_EVAL_BEFORE ) { // negative number is not evaluated
		create_lookup(set_bit_count, set_bit_position);
		b = adios_alac_process(q, gCurrentTimeStep, false);
		initLastConvRid(b);
		q->_maxResultDesired =  calSetBitsNum(b);
		q->_lastRead = 0;
		q->_queryInternal = convertALACBitmapTomemstream(b);
		q->_onTimeStep = gCurrentTimeStep;
	}else { //convert void* _internal to ADIOS_ALAC_BITMAP
		b = (ADIOS_ALAC_BITMAP*) malloc(sizeof(ADIOS_ALAC_BITMAP));
		convertMemstreamToALACBitmap(q->_queryInternal, b);
	}
	uint64_t retrievalSize = q->_maxResultDesired - q->_lastRead;
	if (retrievalSize <= 0) {
		(*queryResult) = NULL;
		FreeALACBITMAP(b);
		q->_onTimeStep = NO_EVAL_BEFORE;
		printf(":: ==> no more results to fetch\n");
		return 0;
	}
	if (retrievalSize > batchSize) {
			retrievalSize = batchSize;
	}

	adios_query_alac_build_results(retrievalSize, outputBoundry,b,queryResult);

	/*if (q->_maxResultDesired >= 0) {
		FREE(b->bits); // these data is copied to q->_queryInternal
	}
	FREE(b); // NOTE: only free the structure*/
	FreeALACBITMAP(b);
	q->_lastRead += retrievalSize;
	if (q->_lastRead == q->_maxResultDesired) {
		q->_onTimeStep = NO_EVAL_BEFORE;
		return 0;
	}
	return 1;
}

int adios_query_alac_free_one_node(ADIOS_QUERY* query){
	if (query == NULL) {
		return 0;
	}

	//TODO: confirm, SHOULD WE DO free here?
	//ADIOS_VARINFO* v = adios_inq_var(f, varName);
	adios_free_varinfo(query->_var);

	FREE(query->_condition);
	FREE(query->_dataSlice);
	FREE(query->_value);

	//TODO: confirm: adios_selection_delete(query->_sel);
	//RIGHT NOW, user will free this box


	//fastbit_selection_free(query->_queryInternal);
	FREE(query);
	return 1;

}

int adios_query_alac_free_method(ADIOS_QUERY* query) {

	// free the tree in a bottom-to-up manner
	if (query->_left == NULL && query->_right == NULL) {
		return adios_query_alac_free_one_node(query);
	}else if  (query->_right){
		return adios_query_alac_free_method(query->_right);
	}else if (query->_left) {
		return adios_query_alac_free_method(query->_left);
	}

	return 1;
}

void adios_query_alac_clean_method() { }


#endif
