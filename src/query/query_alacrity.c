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

#define BITNSLOTS64(nb) ((nb + 64 - 1) / 64)
/**** END -- Funcs. that are internal funcs. ********/

static uint8_t bits_in_char[256] = {
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
		B6(0), B6(1), B6(1), B6(2)};

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
		op1->elmSize = calSetBitsNum(op1);

	} else if (operator == ADIOS_QUERY_OP_OR) {

		for (uint64_t i = 0; i < op1->length; i++) {
			op1[i] ^= op2[i];
		}
		op1->elmSize = calSetBitsNum(op1);

	} else {
		printf("Operator[%d] is not surpported now \n ", operator);
	}
}




void adios_alac_check_candidate(ALMetadata *meta, bin_id_t startBin, bin_id_t endBin
		, double hb, double lb
		, uint64_t ridOffset , uint64_t pgRidOffset
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
	bin_offset_t loworderElm = bl->binValues[endBin] - bl->binValues[startBin];
	char * loadLowBytes = (char *) malloc(sizeof(partitionMeta.elementSize)*loworderElm);

		reconstituteData(partitionMeta, startBin, endBin,
										 lowOrderBytes, loadLowBytes);
		uint32_t rid;
		  switch (partitionMeta.elementSize) {
			case sizeof(uint64_t):
			   uint64_t * lowData64 = (uint64_t *)loadLowBytes;
				for(bin_offset_t el= 0; el < loworderElm; el ++){
					if (lowData64[el] <= hb && lowData64[el] >= lb){
						rid = decodeRids[el];
						rid = rid+ pgRidOffset -  ridOffset;
						uint32_t word = (uint32_t) (rid >> 6);
						alacResultBitmap->bits[word]
								|= (1LL << (rid & 0x3F));
						alacResultBitmap->elmSize ++;
					}
				}
				break;
			case sizeof(uint32_t):

				uint32_t * lowData32 = (uint32_t *)loadLowBytes;
						for(bin_offset_t el= 0; el < loworderElm; el ++){
							if (lowData32[el] <= hb && lowData32[el] >= lb){
								rid = decodeRids[el];
								rid = rid+ pgRidOffset -  ridOffset;
								uint32_t word = (uint32_t) (rid >> 6);
								alacResultBitmap->bits[word]
										|= (1LL << (rid & 0x3F));
								alacResultBitmap->elmSize ++;
							}
						}
				break;
			case sizeof(uint16_t):


					uint16_t * lowData16 = (uint16_t *)loadLowBytes;
							for(bin_offset_t el= 0; el < loworderElm; el ++){
								if (lowData16[el] <= hb && lowData16[el] >= lb){
									rid = decodeRids[el];
									rid = rid+ pgRidOffset -  ridOffset;
									uint32_t word = (uint32_t) (rid >> 6);
									alacResultBitmap->bits[word]
											|= (1LL << (rid & 0x3F));
									alacResultBitmap->elmSize ++;
								}
							}


				break;
			case sizeof(uint8_t):

				uint8_t * lowData8 = (uint8_t *)loadLowBytes;
						for(bin_offset_t el= 0; el < loworderElm; el ++){
							if (lowData8[el] <= hb && lowData8[el] >= lb){
								rid = decodeRids[el];
								rid = rid+ pgRidOffset -  ridOffset;
								uint32_t word = (uint32_t) (rid >> 6);
								alacResultBitmap->bits[word]
										|= (1LL << (rid & 0x3F));
								alacResultBitmap->elmSize ++;
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
	uint64_t resultCount;
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


		uint64_t pgRidOffset = srcstart[0];
		// if for estimate, it will stop here
		// otherwise, continue searching for the precise results
		// TODO: verify this
		for (int j = 1; j < pgBB->ndim; j++) {
			pgRidOffset += (srcstart[j] * srccount[j - 1]);
		}



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
		int ndim = 1;

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
			//			adiosQuery->_var.
			adios_schedule_read_byId(adiosQuery->_f, sel,
					adiosQuery->_var.varid, startStep, numStep, indexData);
			adios_perform_reads(adiosQuery->_f, 1);
			ALIndex index = (ALIndex) indexData;
			const ALBinLayout * const bl = &(partitionMeta.binLayout);
			resultCount = bl->binStartOffsets[hi_bin]
					- bl->binStartOffsets[low_bin];

			const ALIndexForm indexForm = partitionMeta.indexMeta.indexForm;
			ALIndex* indexPtr = &index;
			const bin_offset_t low_bin_elem_off = bl->binStartOffsets[low_bin]; // offset in elements
			const bin_offset_t hi_bin_elem_off = bl->binStartOffsets[hi_bin];
			const bin_offset_t outputCount = hi_bin_elem_off - low_bin_elem_off;

			const char * input_index = (char*) *indexPtr;
			const char *inputCurPtr = input_index;
			const uint64_t *compBinStartOffs =
					partitionMeta.indexMeta.u.ciim.indexBinStartOffsets;
			uint64_t binCompressedLen;

			if (indexForm == ALCompressedInvertedIndex) {

				if (estimate) {

					// Now compress each bin in turn
					for (bin_id_t bin = low_bin; bin < hi_bin; bin++) {
						binCompressedLen = compBinStartOffs[bin + 1]
								- compBinStartOffs[bin];

						ALDecompressRIDtoSelBox(isPGContained
								, inputCurPtr
								, binCompressedLen
								, srcstart //PG region dimension
								, srccount
								, deststart //region dimension of Selection box
								, destcount
								, pgBB.ndim
								, &alacResultBitmap.bits);

						/*ALRLEDecompressRIDs(inputCurPtr, binCompressedLen,
								&tmp_bitmap);*/
						inputCurPtr += binCompressedLen;
					}
				}else{

					// Read low-order bytes
					ADIOS_SELECTION * sel1;
					char * lowOrderBytes = (char *) malloc(dataSizes[j]);
					sel1 = adios_selection_writeblock_bounded(blockId, metaSize+indexSizes[j], dataSizes[j], 1);
					adios_schedule_read(adiosQuery->_f, sel1, adiosQuery->_value, startStep,
							numStep, lowOrderBytes);
					adios_perform_reads(adiosQuery->_f, 1);


					/*uint64_t firstBinElmBegin = 0, firstBinElmEnd = 0,
							lastBinElmBegin = 0, lastBinElmEnd = 0;*/
					// It touches at least 3 bins, so, we need to check RIDs that are in first and last bins
					if (hi_bin - low_bin > 2) {
						// low boundary bin
						binCompressedLen = compBinStartOffs[low_bin + 1]
															- compBinStartOffs[low_bin];

						adios_alac_check_candidate(&partitionMeta, low_bin, low_bin+1
								, hb, lb
								, ridOffset , pgRidOffset
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

							ALDecompressRIDtoSelBox(isPGContained
									, inputCurPtr
									, binCompressedLen
									, srcstart //PG region dimension
									, srccount
									, deststart //region dimension of Selection box
									, destcount
									, pgBB.ndim
									, &alacResultBitmap.bits);
							inputCurPtr += binCompressedLen;
						}


						// high boundary bin
						binCompressedLen = compBinStartOffs[hi_bin]
															- compBinStartOffs[hi_bin-1];

						adios_alac_check_candidate(&partitionMeta, hi_bin-1, hi_bin
								, hb, lb
								, ridOffset , pgRidOffset
								, inputCurPtr /*index bytes of entire PG*/
								, true  /*need decoding*/
								, lowOrderBytes /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);
						inputCurPtr += binCompressedLen;

					} else { // for 1 or 2 bins touched, we need to check all RIDs

						adios_alac_check_candidate(&partitionMeta, low_bin, hi_bin
								, hb, lb
								, ridOffset , pgRidOffset
								, inputCurPtr /*index bytes of entire PG*/
								, true
								, lowOrderBytes /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);

					}


				}
			}  else if (indexForm == ALInvertedIndex) {
				// indexes are inverted indexes that are not compressed
				// we build bitmaps for each rid;
				if (estimate) {
					rid_t * idx = (rid_t *) index;
					uint64_t ni = 0;
					/*for each index, set bit as 1 at corresponding position */
					for (; ni < resultCount; ni++) {
						rid_t rid_val = idx[ni];
						rid_val = rid_val + pgRidOffset -  ridOffset;
						uint32_t word = (uint32_t) (rid_val >> 6);
						alacResultBitmap.bits[word]
								|= (1LL << (rid_val & 0x3F));
					}
					alacResultBitmap.elmSize += resultCount;
				} else {

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

						binCompressedLen = bl->binStartOffsets[low_bin + 1]
								- bl->binStartOffsets[low_bin];

						uint64_t hiBinElm = bl->binStartOffsets[hi_bin]
								- bl->binStartOffsets[hi_bin-1];

						// low boundary bin
						adios_alac_check_candidate(&partitionMeta, low_bin, low_bin+1
								, hb, lb
								, ridOffset , pgRidOffset
								, iptPtr /*index bytes of entire PG*/
								, false  /*don't need decoding*/
								, lowOrderBytes2 /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);

						iptPtr += binCompressedLen;

						rid_t * innerIdx = (rid_t *)  iptPtr;
						uint64_t ni = 0;
						/*for each index, set bit as 1 at corresponding position */
						for (; ni < (resultCount- binCompressedLen - hiBinElm); ni++) {
							rid_t rid_val = innerIdx[ni];
							rid_val = rid_val + pgRidOffset -  ridOffset;
							uint32_t word = (uint32_t) (rid_val >> 6);
							alacResultBitmap.bits[word]
									|= (1LL << (rid_val & 0x3F));

						}
						innerIdx  +=   (resultCount- binCompressedLen - hiBinElm);
						alacResultBitmap.elmSize += (resultCount- binCompressedLen - hiBinElm);
						// high boundary bin
						binCompressedLen = compBinStartOffs[hi_bin]
															- compBinStartOffs[hi_bin-1];

						iptPtr = (int *) innerIdx;
						adios_alac_check_candidate(&partitionMeta, hi_bin-1, hi_bin
								, hb, lb
								, ridOffset , pgRidOffset
								, iptPtr /*index bytes of entire PG*/
								, false
								, lowOrderBytes2 /*low order bytes of entire PG*/
								,&alacResultBitmap /*OUT*/);

					} else { // for 1 or 2 bins touched, we need to check all RIDs

						adios_alac_check_candidate(&partitionMeta, low_bin, hi_bin
								, hb, lb
								, ridOffset , pgRidOffset
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
	if (q ->_left == NULL && q->_right == NULL) {
		// TODO: using existing query selection to select start and end PGs
		int startPG = 0;
		int endPG = q->_var->sum_nblocks;
		// TODO: check the (endPG - startPG)
		return adios_alac_uniengine(q, timestep, estimate);

	} else if (q->_right) {
		return adios_alac_process((ADIOS_QUERY*) q->_right, timestep, estimate);
	} else if (q->_left) {
		return adios_alac_process((ADIOS_QUERY*) q->_left, timestep, estimate);
	}
}

void adios_query_alac_init_method() {

}

int64_t adios_query_alac_estimate_method(ADIOS_QUERY* q) {
	ADIOS_ALAC_BITMAP* b = adios_alac_process(q, false);
	return calSetBitsNum(b);
}

int64_t adios_query_alac_evaluate_method(ADIOS_QUERY* q, int timeStep,
		uint64_t _maxResult) {

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
