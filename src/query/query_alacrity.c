/*
 * query_alacrity.c
 *
 *  Created on: Jun 1, 2014
 *      Author: xczou
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "mpi.h"
#include "adios_read_ext.h"
#include "adios_query.h"
//#include <alacrity.h>
#include "alacrity.h"

/**** Funcs. that are internal funcs. ********/

ADIOS_ALAC_BITMAP * adios_alac_process(ADIOS_QUERY* q, int timeStep,
		bool estimate);

ADIOS_ALAC_BITMAP * adios_alac_bitsOp(ADIOS_ALAC_BITMAP * op1,
		ADIOS_ALAC_BITMAP * op2, enum ADIOS_CLAUSE_OP_MODE operator);

uint64_t calSetBitsNum(ADIOS_ALAC_BITMAP *b1);

int coordinateConversionWithCheck(int * coordinates, const  int dim
		, const  int *srcstart, const  int *deststart, const  int *destend);

bool ridConversionWithCheck(rid_t rid/*relative to local src selectoin*/
		, int *srcstart, int *srccount, int *deststart, int *destcount,
		int dim, rid_t *relativeRid  );

void create_lookup(unsigned char set_bit_count[],
		unsigned char set_bit_position[][16]);

#define BITNSLOTS64(nb) ((nb + 64 - 1) / 64)
/**** END -- Funcs. that are internal funcs. ********/

void create_lookup(unsigned char set_bit_count[],
		unsigned char set_bit_position[][16]) {
	memset(set_bit_count, 0, 256);
	for (int i = 0; i < 65536; i++) {
//		set_bit_count[i] = __builtin_popcount(i); // total bit 1 for value
		set_bit_count[i] =  bits_in_char [i & 0xff]
						   +  bits_in_char [(i >>  8) & 0xff]
						   +  bits_in_char [(i >> 16) & 0xff]
						   +  bits_in_char [(i >> 24) & 0xff]
						   ;
		unsigned short int temp = i;
		int counter = set_bit_count[i] - 1;
		for (int j = 15; j >= 0; j--) {
			unsigned int temp1 = temp >> j & 0x0001;
			if (temp1 == 1) {
				set_bit_position[i][counter--] = j;
			}
		}

	}
}


static uint8_t bits_in_char[256] = {
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
		B6(0), B6(1), B6(1), B6(2)};


unsigned char set_bit_count[65536];
	unsigned char set_bit_position[65536][16];


uint64_t calSetBitsNum(ADIOS_ALAC_BITMAP *b1) {

	uint64_t total = 0;
	int kk = 0;
	for (; kk < b1->length; kk++) {
		uint64_t count = bits_in_char[b1->bits[kk] & 0xff]
				+ bits_in_char[(b1->bits[kk] >> 8) & 0xff]
				+ bits_in_char[(b1->bits[kk] >> 16) & 0xff]
				+ bits_in_char[(b1->bits[kk] >> 24) & 0xff]
				+ bits_in_char[(b1->bits[kk] >> 32) & 0xff]
				+ bits_in_char[(b1->bits[kk] >> 40) & 0xff]
				+ bits_in_char[(b1->bits[kk] >> 48) & 0xff]
				+ bits_in_char[(b1->bits[kk] >> 56) & 0xff];
		total += count;
	}
	return total;
}

// Supports bitmap AND or OR operation.
// The space of second operand is freed, and the first operand serves as the results of the operation
ADIOS_ALAC_BITMAP * adios_alac_bitsOp(ADIOS_ALAC_BITMAP * op1,
		ADIOS_ALAC_BITMAP * op2, enum ADIOS_CLAUSE_OP_MODE operator) {
	if (operator == ADIOS_QUERY_OP_AND) {
		for (uint64_t i = 0; i < op1->length; i++) {
			op1[i] &= op2[i];
		}

		free(op2->bits);
//		op1->elmSize = calSetBitsNum(op1);

	} else if (operator == ADIOS_QUERY_OP_OR) {

		for (uint64_t i = 0; i < op1->length; i++) {
			op1[i] ^= op2[i];
		}

		free(op2->bits);
//		op1->elmSize = calSetBitsNum(op1);

	} else {
		printf("Operator[%d] is not surpported now \n ", operator);
	}
	return op1;
}


int coordinateConversionWithCheck(int * coordinates, const  int dim, const  int *srcstart, const  int *deststart, const  int *destend){

	for (int i = 0; i < dim; i++) {
		coordinates[i] += (srcstart[i] /*global coordinate*/ );
		if ( coordinates[i] < deststart[i] || coordinates[i] > destend[i]){
			return (i+1) * -1;
		}
	}

	/*change coordinate to the destination box*/
	for (int i = 0; i < dim; i++) {
		coordinates[i] -=  deststart[i];
	}
	return 1;
}


/* Give a rid that is relative to a src region
 * return a rid that is relative to dest selection box
 * Assume all the start & count array has slowest dimension at first position
 */
bool ridConversionWithCheck(rid_t rid/*relative to local src selectoin*/, int *srcstart, int *srccount, int *deststart, int *destcount,
		int dim, rid_t *relativeRid  ){

	*relativeRid = 0;
	int * coordinates = (int *) malloc(sizeof(int) * dim); // coordinate of current PG
	int * destend = (int*) malloc(sizeof(int)*dim); // coordinate of ending points on the destination box
	for(int i = 0; i < dim; i ++){
		destend[i] = deststart[i] + destcount[i] -1;
	}
		if (dim == 3) {
			coordinates[0] = rid / (srccount[1] * srccount[2]) ;
			coordinates[1] = (rid % (srccount[1] * srccount[2])) / srccount[2];
			coordinates[2] = (rid % (srccount[1] * srccount[2])) % srccount[2] ;

			if (coordinateConversionWithCheck(coordinates, dim, srcstart, deststart, destend) < 0){
				free(coordinates);
				free(destend);
				return false;
			}

			*relativeRid = coordinates[2] + coordinates[1] * destcount[2] + coordinates[0]* destcount[1] * destcount[2];

		}

		if (dim == 2){

			coordinates[0] = rid / (srccount[1]);
			coordinates[1] = rid % (srccount[1] );

			if (coordinateConversionWitCheck(coordinates, dim, srcstart, deststart, destend) < 0){
				free(coordinates);
				free(destend);
				return false;
			}

			*relativeRid = coordinates[1] + coordinates[0] * destcount[1] ;

		}

	free(destend);
	free(coordinates);
	return true;
}


void adios_alac_check_candidate(ALMetadata *meta, bin_id_t startBin, bin_id_t endBin
		, double hb, double lb
		, int *srcstart, int *srccount //PG region dimension
		, int *deststart, int *destcount //region dimension of Selection box
		, int dim
		, const char *inputCurPtr /*index bytes of entire PG*/
		, bool decoded  /*true: need decoding */
		, char * lowOrderBytes /*low order bytes of entire PG*/
		,ADIOS_ALAC_BITMAP * alacResultBitmap /*OUT*/){
	ALMetadata partitionMeta = *meta;
	const ALBinLayout * const bl = &(partitionMeta.binLayout);
	uint64_t  binCompressedLen = bl->binStartOffsets[endBin]
										- bl->binStartOffsets[startBin];
	uint32_t * decodeRids;
	if (decoded){
		decodeRids= (uint32_t*) malloc(sizeof(uint32_t)*binCompressedLen);
		ALDecompressRIDs( inputCurPtr, binCompressedLen,
			decodeRids, &binCompressedLen);
	}else{
		decodeRids = (uint32_t *) inputCurPtr;
	}
	bin_offset_t loworderElm = bl->binStartOffsets[endBin] - bl->binStartOffsets[startBin];
	char * loadLowBytes = (char *) malloc(sizeof(partitionMeta.elementSize)*loworderElm);

		reconstituteData(partitionMeta, startBin, endBin,
										 lowOrderBytes, loadLowBytes);
		rid_t rid, newRid;
		  switch (partitionMeta.elementSize) {
			case sizeof(uint64_t):
			   uint64_t * lowData64 = (uint64_t *)loadLowBytes;
				for(bin_offset_t el= 0; el < loworderElm; el ++){
					if (lowData64[el] <= hb && lowData64[el] >= lb){
						rid = decodeRids[el];
						if (ridConversionWithCheck(rid,
								srcstart,srccount, deststart,destcount, dim, &newRid)){
							uint32_t word = (uint32_t) (rid >> 6);
							alacResultBitmap.bits[word]
									|= (1LL << (rid & 0x3F));
							alacResultBitmap->elmSize ++;
						}

					}
				}
				break;
			case sizeof(uint32_t):

				uint32_t * lowData32 = (uint32_t *)loadLowBytes;
						for(bin_offset_t el= 0; el < loworderElm; el ++){
							if (lowData32[el] <= hb && lowData32[el] >= lb){
								rid = decodeRids[el];
								if (ridConversionWithCheck(rid,
										srcstart,srccount, deststart,destcount, dim, &newRid)){
									uint32_t word = (uint32_t) (rid >> 6);
									alacResultBitmap.bits[word]
											|= (1LL << (rid & 0x3F));
									alacResultBitmap->elmSize ++;
								}
							}
						}
				break;
			case sizeof(uint16_t):
					uint16_t * lowData16 = (uint16_t *)loadLowBytes;
							for(bin_offset_t el= 0; el < loworderElm; el ++){
								if (lowData16[el] <= hb && lowData16[el] >= lb){
									rid = decodeRids[el];
									if (ridConversionWithCheck(rid,
											srcstart,srccount, deststart,destcount, dim, &newRid)){
										uint32_t word = (uint32_t) (rid >> 6);
										alacResultBitmap.bits[word]
												|= (1LL << (rid & 0x3F));
										alacResultBitmap->elmSize ++;
									}
								}
							}
				break;
			case sizeof(uint8_t):
				uint8_t * lowData8 = (uint8_t *)loadLowBytes;
						for(bin_offset_t el= 0; el < loworderElm; el ++){
							if (lowData8[el] <= hb && lowData8[el] >= lb){
								rid = decodeRids[el];
								if (ridConversionWithCheck(rid,
										srcstart,srccount, deststart,destcount, dim, &newRid)){
									uint32_t word = (uint32_t) (rid >> 6);
									alacResultBitmap.bits[word]
											|= (1LL << (rid & 0x3F));
									alacResultBitmap->elmSize ++;
								}
							}
						}
				break;
		   default:
				eprintf("Unsupported element size %d in %s\n", partitionMeta->elementSize, __FUNCTION__);
				assert(false);
				return 0;
			}

		free(loadLowBytes);
		if(decoded)
			free(decodeRids);

}

/*
 * PG selection, [startPG, endPG)
 */
ADIOS_ALAC_BITMAP* adios_alac_uniengine(ADIOS_QUERY * adiosQuery, int timeStep, bool estimate) {

	double lb = DBL_MIN;
	double hb = DBL_MAX;

	if (adiosQuery->_op == ADIOS_LT || adiosQuery->_op == ADIOS_LTEQ) {
		hb = memcpy(&hb, adiosQuery->_value, sizeof(double));
	} else if (adiosQuery->_op == ADIOS_GT || adiosQuery->_op == ADIOS_GTEQ) {
		lb = memcpy(&lb, adiosQuery->_value, sizeof(double));

	} else if (adiosQuery->_op == ADIOS_EQ) { //following two cases are tricky to ALACRITY
		//TODO: Verify
		hb = memcpy(&hb, adiosQuery->_value, sizeof(double));
		lb = hb;
	} else if (adiosQuery->_op == ADIOS_NE) {
		//TODO
	} else {
		printf("Unsupported predicate type[%d] \n", adiosQuery->_op);
	}


//	adios_find_intersecting_pgs();
	uint64_t ridOffset = 0;
	uint64_t* deststart ; // current variables selection box
	uint64_t* destcount ;
	uint64_t * srcstart;  // PG's bounding box is the global bounding box
	uint64_t * srccount;
	if (adiosQuery->_sel.type == ADIOS_SELECTION_BOUNDINGBOX) {
		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(adiosQuery->_sel.u.bb);
		destcount = bb->count;
		deststart = bb->start;
		ridOffset = deststart[0];
		// if for estimate, it will stop here
		// otherwise, continue searching for the precise results
		// TODO: verify this
		for (int j = 1; adiosQuery->_var.dims; j++) {
			ridOffset += (deststart[j] * deststart[j - 1]);
		}
	} else {
		printf("not supported selection typed in alacrity \n");
		exit(EXIT_FAILURE);
	}

	ALQueryEngine qe;
	ALUnivariateQuery query;
	ALQueryEngineStartUnivariateDoubleQuery(&qe, lb, hb,
			REGION_RETRIEVAL_INDEX_ONLY_QUERY_TYPE, &query);

	/*********** doQuery ***************
	 *
	 * 1. Open partition  [locate offsets of meta, data, and index for the partition]
	 * 2. Read Partition Meta from file => meta
	 3. find touched bins:  are_bins_touched = findBinRange1C(meta, uniquery, &start_bin, &end_bin);
	 4. read index of touched bins: ALPartitionStoreReadIndexBins
	 5. read dataBin
	 */

	//adios_read_set_transforms_enabled(f, 0); // load transformed 1-D byte array
	ADIOS_VARINFO * v = adiosQuery->_var;

	ADIOS_VARTRANSFORM *ti = adios_inq_var_transform(adiosQuery->_f, v);

	int startStep = timeStep, numStep = 1;

	ADIOS_PG_INTERSECTIONS* intersectedPGs = adios_find_intersecting_pgs(
			adiosQuery->_f, v->varid, adiosQuery->_sel, timeStep, numStep);
	int totalPG = intersectedPGs->npg;
	uint64_t totalElm = adiosQuery->_rawDataSize;

	ADIOS_ALAC_BITMAP alacResultBitmap;

	alacResultBitmap.length = BITNSLOTS64(totalElm);
	alacResultBitmap.bits = (uint64_t *) malloc(
			alacResultBitmap.length * sizeof(uint64_t));
	alacResultBitmap.elmSize = 0;


	uint64_t * metaSizes = (uint64_t *) malloc(sizeof(uint64_t) * totalPG);
	uint64_t * indexSizes = (uint64_t *) malloc(sizeof(uint64_t) * totalPG);
	uint64_t * dataSizes = (uint64_t *) malloc(sizeof(uint64_t) * totalPG);


	int j = 0, blockId;
	ADIOS_SELECTION * sel;
	ADIOS_PG_INTERSECTION *  PGs = *(intersectedPGs.intersections);
//	for (i = startPG; i < endPG; i++) {
	for (j = 0; j < totalPG; j++) {
		ADIOS_PG_INTERSECTION pg = PGs[j];
		//TODO: check we should use pg.pg_bounds_sel and pg.blockidx???
		ADIOS_SELECTION * interSelBox = pg.intersection_sel;
		// false: PG selection box is intersecting with variable's selection box
		// true : PG selection box is fully contained within variable's selection box
		bool isPGContained = false;
		ADIOS_SELECTION * pgSelBox = pg.pg_bounds_sel;



		if (pgSelBox.type == ADIOS_SELECTION_BOUNDINGBOX &&
				pgSelBox.type == interSelBox.type ){
			const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *pgBB = &(pgSelBox.u.bb);
			const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *interBB = &(interSelBox.u.bb);
			if (pgBB.ndim == interBB.ndim ){
				isPGContained = true;
				 // if these two boxes are not equal, then,
				// it means two boxes are intersecting / partial overlapping
				for(int k = 0 ; k < pgBB.ndim ; k ++){
					if (pgBB.count[k]!=interBB.count[k] ||
							pgBB.start[k] != interBB.start[k]){
						isPGContained =  false;
					}
				}
			}

		}


		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *pgBB = &(pgSelBox.u.bb);

		srcstart = pgBB.start;
		srccount = pgBB.count;
		blockId = pg.blockidx_in_timestep ;
		ADIOS_TRANSFORM_METADATA * tmetas = ti->transform_metadatas;
		ADIOS_TRANSFORM_METADATA tmeta = tmetas[blockId];
		//	assert(tmeta->length == 24);
		uint64_t * threeData = (uint64_t *) tmeta.content;
		metaSizes[j] = threeData[0];
		indexSizes[j] = threeData[1];
		dataSizes[j] = threeData[2];

		//TODO: offset of each PG should be included

		printf("PG[%d] has meta size[%ul], index size[%ul], and data size[%ul] \n",
				blockId, metaSizes[j], indexSizes[j], dataSizes[j]);

		uint64_t metaSize = metaSizes[j];

		// transformed data has 1 dimension,
		// 1. load partition Metadata
		int ndim = pgBB.ndim;

		sel = adios_selection_writeblock_bounded(blockId, 0, metaSize, 1);
		ALMetadata partitionMeta;
		memstream_t ms = memstreamInitReturn(malloc(metaSize));
		adios_schedule_read(adiosQuery->_f, sel, adiosQuery->_value, startStep,
				numStep, ms.buf);
		adios_perform_reads(adiosQuery->_f, 1);
		ALDeserializeMetadata(&partitionMeta, &ms);
		memstreamDestroy(&ms, true);

		//2. find touched bin
		bin_id_t low_bin, hi_bin;
		_Bool are_bins_touched = findBinRange1C(partitionMeta, query, &low_bin,
				&hi_bin);

		if (are_bins_touched) {

			//3. load index size
			uint64_t indexOffset = metaSize; // metadata size is index start offset

			const char insigbytes = insigBytesCeil(&partitionMeta);
			const uint64_t first_bin_off = indexOffset + ALGetIndexBinOffset(
					&partitionMeta, low_bin); // offset (in BYTES) of first touched bin
			const uint64_t last_bin_off = indexOffset + ALGetIndexBinOffset(
					&partitionMeta, hi_bin); // offset (in BYTES) of first touched bin
			const uint64_t bin_read_len = last_bin_off - first_bin_off; // in bytes
			sel = adios_selection_writeblock_bounded(blockId, metaSize, bin_read_len,
					1);
			void *indexData = NULL;
			adios_schedule_read_byId(adiosQuery->_f, sel,
					adiosQuery->_var.varid, startStep, numStep, indexData);
			adios_perform_reads(adiosQuery->_f, 1);
			ALIndex index = (ALIndex) indexData;
			const ALBinLayout * const bl = &(partitionMeta.binLayout);
			const ALIndexForm indexForm = partitionMeta.indexMeta.indexForm;
			ALIndex* indexPtr = &index;
			const char * input_index = (char*)*indexPtr;
			const char *inputCurPtr = input_index;

			//TODO: distinguish the offset btw two bins for compressed and uncompressed index
			// is the offset byte-level or element-level?
			if (indexForm == ALCompressedInvertedIndex) {
				const uint64_t *compBinStartOffs =
						partitionMeta.indexMeta.u.ciim.indexBinStartOffsets;
				uint64_t binCompressedLen;

				if (estimate) {

					// Now compress each bin in turn
					for (bin_id_t bin = low_bin; bin < hi_bin; bin++) {
						binCompressedLen = compBinStartOffs[bin + 1]
								- compBinStartOffs[bin];
						uint32_t decodedElm = ALDecompressRIDtoSelBox(isPGContained
								, inputCurPtr, binCompressedLen
								, srcstart, srccount //PG region dimension
								, deststart, destcount //region dimension of Selection box
								, ndim , &alacResultBitmap.bits);
						inputCurPtr += binCompressedLen;
						alacResultBitmap.elmSize += decodedElm;
					}

				}else{

					// Read low-order bytes
					ADIOS_SELECTION * sel1;
					char * lowOrderBytes = (char *) malloc(dataSizes[j]);
					sel1 = adios_selection_writeblock_bounded(blockId, metaSize+indexSizes[j], dataSizes[j], 1);
					adios_schedule_read(adiosQuery->_f, sel1, adiosQuery->_value, startStep,
							numStep, lowOrderBytes);
					adios_perform_reads(adiosQuery->_f, 1);

					// It touches at least 3 bins, so, we need to check RIDs that are in first and last bins
					if (hi_bin - low_bin > 2) {
						// low boundary bin, compressed byte offset
						binCompressedLen = compBinStartOffs[low_bin + 1]
															- compBinStartOffs[low_bin];

						adios_alac_check_candidate(&partitionMeta, low_bin, low_bin+1
								, hb, lb
								, srcstart, srccount //PG region dimension
								, deststart, destcount //region dimension of Selection box
								, ndim
								, inputCurPtr /*index bytes of entire PG*/
								, true  /*need decoding*/
								, lowOrderBytes /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);
						inputCurPtr += binCompressedLen;

						bin_id_t innerlowBin = low_bin + 1;
						bin_id_t innerHiBin = hi_bin -1;
						// Now compress each bin in turn
						for (bin_id_t bin = innerlowBin; bin < innerHiBin; bin++) {
							binCompressedLen = compBinStartOffs[bin + 1]
									- compBinStartOffs[bin];

							uint32_t decodedElm = ALDecompressRIDtoSelBox(isPGContained
									, inputCurPtr, binCompressedLen
									, srcstart, srccount //PG region dimension
									, deststart, destcount //region dimension of Selection box
									, ndim , &alacResultBitmap.bits);
							alacResultBitmap.elmSize += decodedElm;
							inputCurPtr += binCompressedLen;
						}


						// high boundary bin
						binCompressedLen = compBinStartOffs[hi_bin]
															- compBinStartOffs[hi_bin-1];

						adios_alac_check_candidate(&partitionMeta, hi_bin-1, hi_bin
								, hb, lb
								, srcstart, srccount //PG region dimension
								, deststart, destcount //region dimension of Selection box
								, ndim
								, inputCurPtr /*index bytes of entire PG*/
								, true  /*need decoding*/
								, lowOrderBytes /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);
						inputCurPtr += binCompressedLen;

					} else { // for 1 or 2 bins touched, we need to check all RIDs

						adios_alac_check_candidate(&partitionMeta, low_bin, hi_bin
								, hb, lb
								, srcstart, srccount //PG region dimension
								, deststart, destcount //region dimension of Selection box
								, ndim
								, inputCurPtr /*index bytes of entire PG*/
								, true
								, lowOrderBytes /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);

					}
				}
			}  else if (indexForm == ALInvertedIndex) {
				// indexes are inverted indexes that are not compressed
				// we build bitmaps for each rid;
				// element offset, instead of byte element
				uint64_t resultCount = bl->binStartOffsets[hi_bin] - bl->binStartOffsets[low_bin];

				if (estimate) {
					rid_t * idx = (rid_t *) index;
					uint32_t newRid ;
					/*for each index, set bit as 1 at corresponding position */
					for (uint64_t ni = 0; ni < resultCount; ni++) {
						rid_t rid_val = idx[ni];
						if (ridConversionWithCheck(rid_val,
								srcstart,srccount, deststart,destcount, pgBB.ndim, &newRid)){
							uint32_t word = (uint32_t) (rid_val >> 6);
							alacResultBitmap.bits[word]
									|= (1LL << (rid_val & 0x3F));
						}

					}
					alacResultBitmap.elmSize += resultCount;
				} else {
					rid_t newRid;
					ADIOS_SELECTION * sel2;
					char * lowOrderBytes2 = (char *) malloc(dataSizes[j]);
					sel2 = adios_selection_writeblock_bounded(blockId, metaSize+indexSizes[j], dataSizes[j], 1);
					adios_schedule_read(adiosQuery->_f, sel2, adiosQuery->_value, startStep,
							numStep, lowOrderBytes2);
					adios_perform_reads(adiosQuery->_f, 1);

					rid_t * decodedRid = (rid_t *) index;
					char * iptPtr = (char*)decodedRid;
					/*uint64_t firstBinElmBegin = 0, firstBinElmEnd = 0,
							lastBinElmBegin = 0, lastBinElmEnd = 0;*/
					// It touches at least 3 bins, so, we need to check RIDs that are in first and last bins
					if (hi_bin - low_bin > 2) {
						bin_offset_t lowBinElm = bl->binStartOffsets[low_bin + 1] - bl->binStartOffsets[low_bin];
						uint64_t hiBinElm = bl->binStartOffsets[hi_bin] - bl->binStartOffsets[hi_bin-1];
						// low boundary bin
						adios_alac_check_candidate(&partitionMeta, low_bin, low_bin+1
								, hb, lb
								, srcstart, srccount //PG region dimension
								, deststart, destcount //region dimension of Selection box
								, ndim
								, iptPtr /*index bytes of entire PG*/
								, false  /*don't need decoding*/
								, lowOrderBytes2 /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);

						iptPtr += (partitionMeta.elementSize * lowBinElm);

						uint64_t innerElm = resultCount- lowBinElm - hiBinElm;
						rid_t * innerIdx = (rid_t *)  iptPtr;
						/*for each index, set bit as 1 at corresponding position */
						for (uint64_t ni = 0; ni < innerElm; ni++) {
							rid_t rid_val = innerIdx[ni];
							if (ridConversionWithCheck(rid_val,
									srcstart,srccount, deststart,destcount, pgBB.ndim, &newRid)){
								uint32_t word = (uint32_t) (rid_val >> 6);
								alacResultBitmap.bits[word]
										|= (1LL << (rid_val & 0x3F));
							}

						}
						innerIdx  +=  innerElm;
						alacResultBitmap.elmSize += innerElm;
						// high boundary bin
						iptPtr = (int *) innerIdx;
						adios_alac_check_candidate(&partitionMeta, hi_bin-1, hi_bin
								, hb, lb
								, srcstart, srccount //PG region dimension
								, deststart, destcount //region dimension of Selection box
								, ndim
								, iptPtr /*index bytes of entire PG*/
								, false
								, lowOrderBytes2 /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);

					} else { // for 1 or 2 bins touched, we need to check all RIDs
						adios_alac_check_candidate(&partitionMeta, low_bin, hi_bin
								, hb, lb
								, srcstart, srccount //PG region dimension
								, deststart, destcount //region dimension of Selection box
								, ndim
								, iptPtr /*index bytes of entire PG*/
								, false
								, lowOrderBytes2 /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);
					}

				}

			} else {
				printf("index form %d in alacrity is not supported", indexForm);
				exit(EXIT_FAILURE);
			}
			FREE(input_index);
		}

	}

	//TODO: free intersectedPGs


	FREE(metaSizes);
	FREE(indexSizes);
	FREE(dataSizes);
	return &(alacResultBitmap);

}

/*
 * This is an internal function processing the expression tree
 */
ADIOS_ALAC_BITMAP * adios_alac_process(ADIOS_QUERY* q, int timestep,
		bool estimate) {

	//LEAF NODE
	ADIOS_ALAC_BITMAP * rbitmap, lbitmap;
	if (q ->_left == NULL && q->_right == NULL) {
		return adios_alac_uniengine(q, timestep, estimate);
	} else if (q->_right) {
		rbitmap = adios_alac_process((ADIOS_QUERY*) q->_right, timestep, estimate);
	} else if (q->_left) {
		lbitmap = adios_alac_process((ADIOS_QUERY*) q->_left, timestep, estimate);
	}

	return adios_alac_bitsOp(rbitmap, lbitmap, q->_leftToRightOp );


}


void adios_query_alac_init_method() {}


void adios_query_alac_retrieval_points2d(
		ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size
		, uint64_t lastRetrievalPos
		, ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb
		, uint64_t **points /*OUT*/
		){

	//TODO: indicates the left bits in the bit array
	uint64_t start_pos = lastRetrievalPos / 64; // each array element has 64 bits
	uint64_t end_pos = retrieval_size / 64 ;
	if ( retrieval_size % 64 > 0)
		end_pos ++;
	uint64_t remain = lastRetrievalPos % 64, count = 0;
	uint64_t * p_bitmap = b->bits;
	uint64_t pidx = 0;

	for (uint64_t off= start_pos; off <= end_pos ; off++) {
		uint64_t offset_long_int = off * 64; // original index offset
		// 2 bytes (unsigned short int)  = 16 bits
		// 4 bytes (unsigned long int )= 64 bit
		uint16_t * temp = (uint16_t *) &p_bitmap[off];
		uint64_t offset;
		for (int j = 0; j < 4; j++) {
			offset = offset_long_int + j * 16; // here, 16 is used because temp is 16bits (unsigned short int) pointer
			// set_bit_count for each 2 bytes, the number of 1
			/*
			 * *******|               64 bits                 | => final_result_bitmap []
			 * *******| 16 bits | 16 bits | 16 bits | 16 bits | => temp[]
			 */
			for (int m = 0; m < set_bit_count[temp[j]] && (count <= retrieval_size) ; m++) {
				uint64_t reconstct_rid = offset+ set_bit_position[temp[j]][m];
				if (count > remain ){
					(*points)[pidx++] = bb->start[0] + reconstct_rid / bb->count[1];
					(*points)[pidx++] = bb->start[1] + reconstct_rid % bb->count[1];
				}
				count ++;
			}
		}
	}

}

void adios_query_alac_retrieval_points3d(
		ADIOS_ALAC_BITMAP *b, uint64_t retrieval_size
		, uint64_t lastRetrievalPos
		, ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb
		, uint64_t **points /*OUT*/
		){

	//TODO: indicates the left bits in the bit array
	uint64_t start_pos = lastRetrievalPos / 64; // each array element has 64 bits
	uint64_t end_pos = retrieval_size / 64 ;
	if ( retrieval_size % 64 > 0)
		end_pos ++;
	uint64_t remain = lastRetrievalPos % 64, count = 0;
	uint64_t * p_bitmap = b->bits;
	uint64_t pidx = 0;

	for (uint64_t off= start_pos; off <= end_pos ; off++) {
		uint64_t offset_long_int = off * 64; // original index offset
		// 2 bytes (unsigned short int)  = 16 bits
		// 4 bytes (unsigned long int )= 64 bit
		uint16_t * temp = (uint16_t *) &p_bitmap[off];
		uint64_t offset;
		for (int j = 0; j < 4; j++) {
			offset = offset_long_int + j * 16; // here, 16 is used because temp is 16bits (unsigned short int) pointer
			// set_bit_count for each 2 bytes, the number of 1
			/*
			 * *******|               64 bits                 | => final_result_bitmap []
			 * *******| 16 bits | 16 bits | 16 bits | 16 bits | => temp[]
			 */
			for (int m = 0; m < set_bit_count[temp[j]] && (count <= retrieval_size) ; m++) {
				uint64_t reconstct_rid = offset+ set_bit_position[temp[j]][m];
				if (count > remain ){
					(*points)[pidx++] = bb->start[0] + reconstct_rid / (bb->count[1] * bb->count[2]) ;
					(*points)[pidx++] = bb->start[1] + (reconstct_rid % (bb->count[1] * bb->count[2])) / bb->count[2];
					(*points)[pidx++] = bb->start[2] + (reconstct_rid % (bb->count[1] * bb->count[2])) % bb->count[2] ;
				}
				count ++;
			}
		}
	}

}

void adios_query_alac_build_results(
		uint64_t retrieval_size, ADIOS_QUERY * q, ADIOS_ALAC_BITMAP *b
		, ADIOS_SELECTION ** queryResult){

	//last bounding box / points supplied by user
	ADIOS_SELECTION * outputBoundry = q->_sel;
	uint64_t offset = q->_lastRead;
	switch (outputBoundry->type) {
	case ADIOS_SELECTION_BOUNDINGBOX: {
		const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(outputBoundry->u.bb);

		uint64_t dataSize = retrieval_size * (bb->ndim);
		uint64_t* points = (uint64_t*) (malloc(dataSize * sizeof(uint64_t)));
		if ( bb->ndim == 3){
			adios_query_alac_retrieval_points3d(b,retrieval_size,offset, queryResult, bb);
		}else if (bb -> ndim == 2){
			adios_query_alac_retrieval_points2d(b,retrieval_size,offset, queryResult, bb);
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
	ADIOS_ALAC_BITMAP* b = adios_alac_process(q, false);
	return calSetBitsNum(b);
}

/* memory format:
 * | length | elmSize | bits ....
 */
void convertALACBitmapTomemstream( ADIOS_ALAC_BITMAP * b, void ** mem /*OUT*/){

	uint64_t totalSize = b->length + 1 /*b->length*/ + 1 /*b->elmSize*/;
	uint64_t * ptr  = (uint64_t *) malloc(sizeof(uint64_t)* totalSize);
	(*mem) = ptr;
	ptr[0] = b->length;   ptr[1] = b->elmSize;
	ptr += 2;
	memcpy(ptr, b->bits, sizeof(uint64_t)*(totalSize-2));

}

void convertMemstreamToALACBitmap( void *mem , ADIOS_ALAC_BITMAP ** bout /*OUT*/){

	ADIOS_ALAC_BITMAP * b = *bout;
	b = (ADIOS_ALAC_BITMAP *) malloc(sizeof(ADIOS_ALAC_BITMAP));
	uint64_t * ptr  = (uint64_t *) mem;
	b ->length = ptr[0];
	b -> elmSize = ptr[1];
	b-> bits = (ptr+2);
}



int  adios_query_get_selection(ADIOS_QUERY* q,
			       uint64_t batchSize, // limited by maxResult
			       ADIOS_SELECTION* outputBoundry,
			       ADIOS_SELECTION** queryResult) {
	// first time, we have to evaluate it
	ADIOS_ALAC_BITMAP* b ;
	if (q->_maxResultDesired < 0) { // negative number is not evaluated
		create_lookup(set_bit_count, set_bit_position);
		b = adios_alac_process(q, true);
		q->_maxResultDesired =  calSetBitsNum(b);
		q->_lastRead = 0;
		//convert ADIOS_ALAC_BITMAP to  " void* _internal"
		convertALACBitmapTomemstream(b, &(q->_queryInternal));
	}else { //convert void* _internal to ADIOS_ALAC_BITMAP
		convertMemstreamToALACBitmap(q->_queryInternal, &b);
	}
	uint64_t retrievalSize = q->_maxResultDesired - q->_lastRead;
	if (retrievalSize == 0) {
		(*queryResult) = NULL;
		printf(":: ==> no more results to fetch\n");
		return 0;
	}
	if (retrievalSize > batchSize) {
			retrievalSize = batchSize;
	}

	adios_query_alac_build_results(retrievalSize,q,b,queryResult);

	if (q->_maxResultDesired < 0) { // negative number is not evaluated
		free(b->bits); // these data is copied to q->_queryInternal
	}
	free(b); // NOTE: only free the structure
	q->_lastRead += retrievalSize;
	if (q->_lastRead == q->_maxResultDesired) {
		return 0;
	}
	return 1;
}


int adios_query_alac_free_method(ADIOS_QUERY* query) {

	//TODO:
	if (query == NULL) {
		return;
	}

	printf(":: free %s\n", query->_condition);
	free(query->_value);
	free(query->_dataSlice);
	free(query->_condition);

	//adios_selection_delete(query->_sel);
	adios_free_varinfo(query->_var);

	//fastbit_selection_free(query->_queryInternal);
	free(query);

}

void adios_query_alac_clean_method() {
	//  fastbit_iapi_free_all();
	//  fastbit_cleanup();
}
