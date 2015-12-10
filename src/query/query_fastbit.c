#include "public/adios_query.h"
   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/common_read.h"
#include "core/adios_logger.h"
#include "core/futils.h"
#include "fastbit_adios.h"
#include "common_query.h"
#include <iapi.h>
#include <math.h>

#define BITARRAY
#define INT_BIT 32

#define BITMASK(b) (1 << ((b) % INT_BIT))
#define BITSLOT(b) ((b) / INT_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + INT_BIT - 1) / INT_BIT)


int64_t getPosInBox(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* sel, int n, uint64_t* spatialCoordinates, int fortran_order);
int64_t getPosInVariable(const ADIOS_VARINFO* v, int n, uint64_t* spatialCoordinates, int fortran_order);

void getHandleFromBlockAtLeafQuery(int timeStep, int blockIdx, ADIOS_FILE* idxFile, ADIOS_QUERY* q, uint64_t blockSize);

FastBitSelectionHandle createHandle(ADIOS_QUERY* q, const char* registeredArrayName);
int mEvaluateTestFullRangeFancyQueryOnWhole(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep);

static inline void scan_r2(int dim, const uint64_t* const dimSize, uint64_t pos,  uint64_t * const result, uint64_t* const scan_start, int d, uint64_t curr0, uint64_t base)
{
  uint64_t i=0;

  if (d == dim-1) {
    result[d] = pos - curr0;
    return;
  }

  for (i=scan_start[d]; i<dimSize[d]; i++) {      
    uint64_t curr = i*base+curr0;
    uint64_t next = curr+base;

    if (pos >= next) {
      continue;
    }
      
    result[d] = i;
    if (i > scan_start[d]) {
      scan_start[d+1] = 0;      
    }
    
    if (d+1 == dim -2) {
      uint64_t k = 0;
      for (k=scan_start[dim-2]; k<dimSize[dim-2]; k++) {
	uint64_t curr1 = k*dimSize[dim-1]+curr;
	uint64_t next1 = curr1+dimSize[dim-1];
	if (pos >= next1) {
	  continue;
	}
	result[dim-2] = k;
	result[dim-1] = pos-curr1;
	return;
      }
    } else if (d+1 == dim-1)  { // save a recursive step
      result[dim-1] = pos-curr; 
      return;
    } else if (d+1 < dim-2) {
      uint64_t base2=base/dimSize[d+1];
      scan_r2(dim, dimSize, pos, result, scan_start, d+1, curr, base2);
    } 
    break;
  }  

}


static inline void scan3d(const uint64_t* const dimSize, uint64_t pos,  uint64_t * const result, uint64_t* const scan_start)
{

  int i,j,k;

  int base=dimSize[1]*dimSize[2];

  for (i=scan_start[0]; i<dimSize[0]; i++) {
    int curr = i*base;
    int next = curr+base;

    if (pos >= next) {
      continue;
    }
      
    result[0] = i;
    if (i > scan_start[0]) {
      scan_start[1] = 0;
    }

    for (j=scan_start[1]; j<dimSize[1]; j++) {
      int curr2 = curr+j*dimSize[2];
      int next2 = curr2+dimSize[2];
      if (pos >= next2) {
	continue;
      }
      result[1] = j;
      result[2] = pos-curr2;
      return;
    }
  }
}

static inline void scan2d(const uint64_t* const dimSize, uint64_t pos,  uint64_t * const result, uint64_t* const scan_start)
{

  int i,j;

  int base=dimSize[1];

  for (i=scan_start[0]; i<dimSize[0]; i++) {
    int curr = i*base;
    int next = curr+base;

    if (pos >= next) {
      continue;
    }
      
    result[0] = i;
    result[1] = pos-curr;
    return;
  }
}


static void simple_scan(int64_t* coordinateArray, uint64_t start, uint64_t count, 
			ADIOS_VARINFO* v, uint64_t* blockstart, uint64_t* blockcount,const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb)
{
  if (v->ndim == 1) { // to do: check
      uint64_t i = 0;
      uint64_t pos = 0; // pos in variable
      for (i=0; i<count; i++) {
	   pos = coordinateArray[i+start] + blockstart[0];
	   if (bb == NULL) {
	     coordinateArray[i+start] = pos;
	   } else {
	     coordinateArray[i+start] = pos - bb->start[0];
	   }
      }
      return;
  }

  if (v->ndim == 2) {
    uint64_t scan_start[2] = {0,0};
    uint64_t result[2] = {0,0};
    uint64_t pos = 0; // pos in variable
    uint64_t i=0;
    for (i=0; i<count; i++) {
         pos = coordinateArray[i+start];     
	 scan2d(blockcount, pos, result, scan_start) ;      
	 scan_start[0] = result[0];
	 
	 result[0] += blockstart[0];
	 result[1] += blockstart[1];
	 if (bb != NULL) {
	   coordinateArray[i+start] = getPosInBox(bb, bb->ndim, result, 0); 
	 } else {
	   coordinateArray[i+start] = getPosInVariable(v, v->ndim, result, 0); 
	 }
    }
    return;
  }

  if (v->ndim == 3) {
    uint64_t scan_start[3] = {0,0,0}; //{blockstart[0],blockstart[1], blockstart[2]};
    uint64_t result[3] = {0,0,0};
    uint64_t pos = 0; // pos in variable
    uint64_t i=0;
    for (i=0; i<count; i++) {
         pos = coordinateArray[i+start];     
	 scan3d(blockcount, pos, result, scan_start) ;      
	 scan_start[0] = result[0];
	 scan_start[1] = result[1];
	 
	 result[0] += blockstart[0];
	 result[1] += blockstart[1];
	 result[2] += blockstart[2];
	 
	 if (bb != NULL) { 
	   coordinateArray[i+start] = getPosInBox(bb, bb->ndim, result, 0); 
	 } else {
	   coordinateArray[i+start] = getPosInVariable(v, v->ndim, result, 0); 
	 }

	 //printf(" %llu, %llu pos=%llu, (%llu, %llu, %llu) => [%lld] = %lld  \n", i, start, pos, result[0], result[1], result[2], i+start, coordinateArray[i+start]);
    }
  }
}

static void arraymath_assign(uint64_t* left, uint64_t* right, int size)
{
  int i;
  for (i=0; i<size; i++) {
    left[i] = right[i];
  }
}

static void arraymath_pluseq(uint64_t* left, uint64_t* right, int size)
{
  int i;
  for (i=0; i<size; i++) {
    left[i] += right[i];
  }
}

static void arraymath_init(uint64_t* left,  int size)
{
  int i;
  for (i=0; i<size; i++) {
    left[i] = 0;
  }  
}


static void recursive_scan(int64_t* coordinateArray, uint64_t start, uint64_t count, ADIOS_VARINFO* v, 
			   uint64_t* blockstart, uint64_t* blockcount,const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb)
{

  uint64_t scan_start[v->ndim]; 
  uint64_t result[v->ndim];

  arraymath_init(scan_start, v->ndim);
  arraymath_init(result, v->ndim);

  uint64_t i=0; 
  uint64_t pos = 0;

  uint64_t startBase=1;
  for (i=1; i<v->ndim; i++) {
    startBase *= blockcount[i];
  }

  for (i=0; i<count; i++) {
      pos = coordinateArray[i+start];     
      scan_r2(v->ndim, blockcount,  pos, result, scan_start, 0, 0, startBase) ;      

      arraymath_assign(scan_start, result, v->ndim);	 
      arraymath_pluseq(result, blockstart, v->ndim);

      if (bb != NULL) { 
	coordinateArray[i+start] = getPosInBox(bb, bb->ndim, result, 0); 
      } else {
	coordinateArray[i+start] = getPosInVariable(v, v->ndim, result, 0); 
      }
  } 
}

static void fastscan(int64_t* coordinateArray, uint64_t start, uint64_t count, ADIOS_VARINFO* v, uint64_t* blockstart, uint64_t* blockcount,
		     const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb)
{
  if (v->ndim <= 3) {
    simple_scan(coordinateArray, start, count, v, blockstart, blockcount, bb);
  } else {
    recursive_scan(coordinateArray, start, count, v, blockstart, blockcount, bb);
  }
}

//
// map pos back to coordinates
//
static inline void posToSpace(uint64_t pos, const int isFortranClient,  
			     const uint64_t * const dimSize , uint64_t * const coordinates, const int dim,  const uint64_t* const start) {
  uint32_t i;
  uint64_t curr = pos;

  if (isFortranClient) {  
    for (i = 0; i < dim; i++) {
       uint64_t temp0 = curr/dimSize[i];
       uint64_t temp1 = curr-temp0*dimSize[i];

       coordinates[i] = temp1 + start[i];
       curr = temp0;
       //coordinates[i] = curr % dimSize[i]+start[i];
       //curr /= dimSize[i];
    }
  } else {    
    // note using i from dim to 1 b/c i is uint32, will not be -1 for 0--.
    // use unsigned to save compute time on div
    for (i = dim ; i >= 1; i--) {
       uint64_t temp0 = curr/dimSize[i-1];
       uint64_t temp1 = curr-temp0*dimSize[i-1];
       coordinates[i-1] = temp1 + start[i-1];
       curr = temp0;
    }
  }
}


uint32_t* bitarray_create(uint64_t max) 
{
  //uint32_t* result = malloc(BITNSLOTS(max) * sizeof(uint32_t));
  //memset(result, 0, BITNSLOTS(max));
  //uint32_t* result = (uint32_t *)calloc(0, BITNSLOTS(max) * sizeof(uint32_t));
  uint32_t* result = (uint32_t *)calloc(BITNSLOTS(max),  sizeof(uint32_t));
  return result;
}

void bitarray_unsetbit(uint32_t* array, uint64_t pos)
{
  BITCLEAR(array, pos);
}
void bitarray_setbit(uint32_t* array, uint64_t pos)
{
  BITSET(array, pos);
}

static const unsigned char BitsSetTable256[256] = 
  {
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
  };

//unsigned int v; // count the number of bits set in 32-bit value v
uint32_t getSetBits(uint32_t v) 
{
  uint32_t c = 0; // c is the total bits set in v

  // Option 1:
  c = BitsSetTable256[v & 0xff] + 
    BitsSetTable256[(v >> 8) & 0xff] + 
    BitsSetTable256[(v >> 16) & 0xff] + 
    BitsSetTable256[v >> 24]; 
  
  /*
  // Option 2:
  unsigned char * p = (unsigned char *) &v;
  c = BitsSetTable256[p[0]] + 
  BitsSetTable256[p[1]] + 
  BitsSetTable256[p[2]] +
  BitsSetTable256[p[3]];
  */
  return c;
}

uint64_t  bitarray_countHits(uint32_t* bitarray, uint64_t size) {
  struct timespec startT;
  casestudyLogger_getRealtime(&startT);

  uint64_t  i=0;
  uint64_t countHitser = 0;
  for (i=0; i<size; i++) {    
    countHitser += getSetBits(bitarray[i]);
  }

  //casestudyLogger_writeout(&startT, "count bit array hits ");
  casestudyLogger_frame_writeout(&startT, "count bit");
  return countHitser;
}


uint64_t  bitarray_getHits(uint32_t* bitarray, uint64_t* result, uint64_t size, uint64_t start, uint64_t count)
{
  //  printf("\n==> start=%ld count = %ld \n", start, count);

  if (bitarray == NULL) {
    return -1;
  }

  uint64_t i, j;
  uint64_t hitCounter = 0;    
  uint64_t collected = 0;

  for (i=0; i<size; i++) {
    int currCount = getSetBits(bitarray[i]);
    if (hitCounter + currCount > start) {
       for (j=0; j<INT_BIT; j++) {
	    uint64_t bitpos = i*INT_BIT+j;
	    if (BITTEST(bitarray, bitpos)) {
	      hitCounter += 1;
	      if (hitCounter > start) {
		  result[collected] = bitpos;
		  collected ++;
		  //printf("bit at %ld is  hit: %ld. ref: [ %ld th]\n", bitpos, hitCounter, collected);
	      }
	      if (collected  == count) {
		return collected;
	      }
	    }
       }
      break; 
    } 
    hitCounter += currCount;    
  }

  hitCounter += collected;
  //printf(" hit counter check: %ld\n", hitCounter);

  uint64_t rock = i+1;
  uint64_t lastrock = 0;

  for (i=rock; i<size; i++) {
    int currCount = getSetBits(bitarray[i]);
    if (collected + currCount > count) {
      lastrock = i;
      break;
    } 
    //collected += currCount;
    for (j=0; j<INT_BIT; j++) {
         uint64_t bitpos = i*INT_BIT + j;
	 if (BITTEST(bitarray, bitpos)) {
	   result[collected] = bitpos;
	   collected++;
	 }
    }
  }

  hitCounter += collected;
  //  printf(" hit counter check: %ld, collected: %ld\n", hitCounter, collected);

  if (lastrock == 0) {
    //printf("done. \n");
    return collected;
  }

  // check the last rock
  rock = i; 
  for (i=0; i<INT_BIT; i++) {
    uint64_t bitpos = rock*INT_BIT + i;
    if (BITTEST(bitarray, bitpos)) {
        result[collected] = bitpos;
        collected ++; 
	//printf("bit at %ld is  hit: %ld. ref: %ld th\n ", bitpos, hitCounter+collected, collected);
    }
    if (collected == count) {
      break;
    }
  }

  return collected;
}
 

uint64_t  bitarray_countHitsNative(uint64_t max, uint32_t* bitarray1) {
  uint64_t countHitser = 0;				
  uint64_t i;

  // printf("array size= %ld \n", BITNSLOTS(max));
  for (i=0; i<max; i++) {
    if (BITTEST(bitarray1, i)) {
      //printf(" ==> bit at: %ld is set. \n",i);
      countHitser ++;
    }
  }

  //printf("%ld\n", countHitser);
  return countHitser;
}

void bitarray_or(uint32_t* bitarray0, uint32_t* bitarray1, uint32_t* bitarray3, uint64_t arraysize)
{
  uint64_t i;
  for (i=0; i<arraysize; i++) {
    bitarray3[i] = bitarray0[i] | bitarray1[i];
  }  
}

void bitarray_and(uint32_t* bitarray0, uint32_t* bitarray1, uint32_t* bitarray3, uint64_t arraysize)
{
  uint64_t i;
  for (i=0; i<arraysize; i++) {
    bitarray3[i] = bitarray0[i] & bitarray1[i];
  }  
}











typedef struct {
  double* _keys; 
  int64_t* _offsets; 
  uint32_t* _bms;
  char* _arrayName;
  char* _bmsVarName;
  //void* _rawData;
  FastBitSelectionHandle _handle;
  
  ADIOS_FILE* _idxFile;
} FASTBIT_INTERNAL;
 

ADIOS_QUERY* getFirstLeaf(ADIOS_QUERY* q);
void getHandle(int timeStep, int blockIdx, ADIOS_FILE* idxFile, ADIOS_QUERY* q, uint64_t blockSize);

//#define __READ_BMS_AS_NEEDED__

/** A simple reader to be used by FastBit for index reconstruction.  In
    this simple case, the first argument is the whole array storing all the
    serialized bitmaps.  This first argument can be used to point to a data
    structure pointing to any complex object type necassary.
*/
//
// this static function is from fastbit example/tiapi.c
//
static int adios_bmreader(void *ctx, uint64_t start,uint64_t count, uint32_t *buf)
{
  struct timespec startT;
  casestudyLogger_getRealtime(&startT);

#ifdef __READ_BMS_AS_NEEDED__
  FASTBIT_INTERNAL* itn = (FASTBIT_INTERNAL*)ctx;
  // read bms from start to count:
  ADIOS_VARINFO * bmsV = common_read_inq_var (itn->_idxFile, itn->_bmsVarName);
  uint64_t bms_start[] = {start};
  uint64_t bms_count[] = {count};

  ADIOS_SELECTION* bmsSel = common_read_selection_boundingbox(bmsV->ndim, bms_start, bms_count);
  // idx file has one timestep
  common_read_schedule_read(itn->_idxFile, bmsSel, itn->_bmsVarName, 0, 1, NULL, buf);
  common_read_perform_reads(itn->_idxFile,1);
  common_read_free_varinfo(bmsV);
  common_read_selection_delete(bmsSel);

  //casestudyLogger_bms_writeout(&startT, 
  casestudyLogger_bms_writeout(&startT, "bmreader_adv visited ");
  return 0;
#else
  const uint32_t *bms = (uint32_t*)ctx + start;
  unsigned j;
  for (j = 0; j < count; ++ j) {
    buf[j] = bms[j];
  }

  //casestudyLogger_bms_writeout(&startT, "bmreader visited ");
  casestudyLogger_bms_writeout(&startT, "bmreader visited ");
  return 0;
#endif

}


void create_fastbit_internal_idxFile (ADIOS_QUERY* q, ADIOS_FILE* idxFile) 
{
  if (q->queryInternal == NULL) {
     FASTBIT_INTERNAL* internal = malloc(sizeof(FASTBIT_INTERNAL));
     internal->_keys = NULL;
     internal->_offsets = NULL;
     internal->_bms = NULL;
     internal->_arrayName = NULL;
     internal->_bmsVarName = NULL;
     internal->_handle = NULL;

     internal->_idxFile = idxFile;
     q->queryInternal = internal;
  }
  
  if (q->left != NULL) {
    create_fastbit_internal_idxFile(q->left, idxFile);
  } 
  if (q->right != NULL) {
    create_fastbit_internal_idxFile(q->right, idxFile);
  }
}

void create_fastbit_internal (ADIOS_QUERY* q) 
{
  if (q->queryInternal == NULL) {
      ADIOS_QUERY* leaf = getFirstLeaf(q);
      ADIOS_FILE* f = leaf->file;
      const char* basefileName = f->path;

      MPI_Comm comm_dummy = MPI_COMM_SELF;
      
      struct timespec idxStartT;
      casestudyLogger_getRealtime(&idxStartT);

      ADIOS_FILE* idxFile = fastbit_adios_util_getFastbitIndexFileToRead(basefileName, comm_dummy);
      casestudyLogger_idx_writeout(&idxStartT, "index file loaded ");

      create_fastbit_internal_idxFile (q, idxFile);
  }
}


void clear_fastbit_internal(ADIOS_QUERY* query) 
{
  FASTBIT_INTERNAL* s = (FASTBIT_INTERNAL*)(query->queryInternal);

  if (s == NULL) {
    return;
  }

  if (s->_keys != NULL) {
    free(s->_keys); 
  }
  if (s->_offsets != NULL) {
    free(s->_offsets);     
  }
  if (s->_bms != NULL) {
    free(s->_bms); 
  }

  if (s->_bmsVarName != NULL) {
    free(s->_bmsVarName);
    s->_bmsVarName = NULL;
  }

  s->_keys = NULL; s->_bms = NULL; s->_offsets = NULL;
  
  fastbit_iapi_free_array_by_addr(query->dataSlice);

  if (s->_arrayName != NULL) {
    fastbit_iapi_free_array(s->_arrayName);
    free(s->_arrayName);
    s->_arrayName = NULL;
  }
  
  if (query->hasParent == 0) {
    fastbit_selection_free(s->_handle); 
    s->_handle == NULL;
  }   

  //s = NULL;
  //free(s);
}


void clear_fastbit_internal_recursive(ADIOS_QUERY* query) 
{
  clear_fastbit_internal(query);

  /*
  if (query->hasParent != 0) {
    if (((FASTBIT_INTERNAL*)(query->queryInternal))->_handle != NULL) {
      fastbit_selection_free( ((FASTBIT_INTERNAL*)(query->queryInternal))->_handle); 
    }
  }
*/
 
  if (query->left != NULL) {
    clear_fastbit_internal_recursive(query->left);
  }
  if (query->right != NULL) {
    clear_fastbit_internal_recursive(query->right);
  }
}

void clean_fastbit_trace(ADIOS_QUERY* q) 
{

  if (q == NULL) {
    return;
  }

  clear_fastbit_internal_recursive(q);
  /*
  if (q->left == NULL) {
      free (((FASTBIT_INTERNAL*)(q->queryInternal))->_bms); ((FASTBIT_INTERNAL*)(q->queryInternal))->_bms = NULL; 
      free (((FASTBIT_INTERNAL*)(q->queryInternal))->_keys); ((FASTBIT_INTERNAL*)(q->queryInternal))->_keys = NULL; 
      free (((FASTBIT_INTERNAL*)(q->queryInternal))->_offsets); ((FASTBIT_INTERNAL*)(q->queryInternal))->_offsets = NULL; 
      free (((FASTBIT_INTERNAL*)(q->queryInternal))->_offsets); ((FASTBIT_INTERNAL*)(q->queryInternal))->_offsets = NULL; 
      fastbit_selection_free(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
      return;
  }

  ADIOS_QUERY* lq = (ADIOS_QUERY*)(q->left);
  ADIOS_QUERY* rq = (ADIOS_QUERY*)(q->right);

  clean_fastbit_trace(lq);
  clean_fastbit_trace(rq);

  */
}


void assertValue(char* input, char* endptr) {
  if (*endptr != '\0')  
    if ((errno == ERANGE) || (errno != 0) || (endptr == input) || (*endptr != '\0')) {
      //perror("strtol");
        printf("FastbitQueryEvaluation Exit due to :invalid type value: %s\n", input);
	exit(EXIT_FAILURE);
    }
}

void assert(void* ptr, const char* notes) 
{
  if (ptr == NULL) {
    printf("::Error allocating memory: %s\n", notes);
    exit(EXIT_FAILURE);
  }
}

void setQueryInternal(ADIOS_QUERY* q, FastBitCompareType compareOp, FastBitDataType dataType,   uint64_t dataSize, const char* arrayName) 
{
  /*
  char* endptr;
  
  if (dataType == FastBitDataTypeDouble) {
    double vv = strtod(q->_value, &endptr);
    assertValue(q->_value, endptr);
    q->_queryInternal = fastbit_selection_create(dataType, q->_dataSlice, dataSize, compareOp, &vv);
  } else if ((dataType == FastBitDataTypeInt) || (dataType == FastBitDataTypeLong) || (dataType == FastBitDataTypeUInt) || (dataType == FastBitDataTypeULong)) {    
    //long vv = strtol(q->_value, &endptr, 10);
    float vv = strtof(q->_value, &endptr);
    //so this is not a much to use a long value??
    assertValue(q->_value, endptr);
    q->_queryInternal = fastbit_selection_create(dataType, q->_dataSlice, dataSize, compareOp, &vv);
  } else if (dataType == FastBitDataTypeFloat) {
    float vv = strtof(q->_value, &endptr);
    assertValue(q->_value, endptr);
    q->_queryInternal = fastbit_selection_create(dataType, q->_dataSlice, dataSize, compareOp, &vv);
  } else {
    q->_queryInternal = NULL;
  }
  */

 
  fastbit_iapi_free_array_by_addr(q->dataSlice);
  fastbit_iapi_register_array(arrayName, dataType, q->dataSlice, dataSize);

  char* endptr;
  double vv = strtod(q->predicateValue, &endptr);
  FastBitSelectionHandle h = fastbit_selection_osr(arrayName, compareOp, vv);
  fastbit_adios_util_checkNotNull(h, arrayName);
  ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = h;
  //fastbit_adios_util_checkNotNull(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, arrayName);
}




void setCombinedQueryInternal(ADIOS_QUERY* q) {
    ADIOS_QUERY* left = (ADIOS_QUERY*)(q->left);
    ADIOS_QUERY* right = (ADIOS_QUERY*)(q->right);
    
    FASTBIT_INTERNAL* leftHandle = ((FASTBIT_INTERNAL*)left->queryInternal)->_handle;
    FASTBIT_INTERNAL* rightHandle = ((FASTBIT_INTERNAL*)right->queryInternal)->_handle;

    if ((leftHandle == NULL) && (rightHandle == NULL)) {
      log_debug("both fastbit handles are null, do nothing\n");
      // do nothing
    } else if (rightHandle == NULL) {
      log_debug("right fastbit handles are null\n");
      if (q->combineOp == ADIOS_QUERY_OP_OR) {            
	((FASTBIT_INTERNAL*)(q->queryInternal))->_handle =  leftHandle;
      }
    } else if (leftHandle == NULL) {
      log_debug("left fastbit handles are null\n");
      if (q->combineOp == ADIOS_QUERY_OP_OR) {            
	((FASTBIT_INTERNAL*)(q->queryInternal))->_handle =  rightHandle;
      }

    } else {
      if (q->combineOp == ADIOS_QUERY_OP_AND) {            
	((FASTBIT_INTERNAL*)(q->queryInternal))->_handle =  fastbit_selection_combine(leftHandle, FastBitCombineAnd, rightHandle);
      } else {
	((FASTBIT_INTERNAL*)(q->queryInternal))->_handle =  fastbit_selection_combine(leftHandle, FastBitCombineOr, rightHandle);
      }
    }
/*
    if (q->combineOp == ADIOS_QUERY_OP_AND) {      
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle =  fastbit_selection_combine(((FASTBIT_INTERNAL*)left->queryInternal)->_handle, FastBitCombineAnd, 
										    ((FASTBIT_INTERNAL*)right->queryInternal)->_handle);
    } else {
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle =  fastbit_selection_combine(((FASTBIT_INTERNAL*)left->queryInternal)->_handle, FastBitCombineOr, 
										    ((FASTBIT_INTERNAL*)right->queryInternal)->_handle);
    }    
*/
}

void getCoordinateFromPoints(uint64_t pos, const ADIOS_SELECTION_POINTS_STRUCT* sel, uint64_t* coordinates) 
{
  int i=0; 
  for(i=0; i<sel->ndim; i++) {
    (coordinates)[i] = sel->points[pos*(sel->ndim)+i];
  }
}

/*
void getCoordinateFromBlock(uint64_t pos, const ADIOS_VARBLOCK* sel, int n, uint64_t* coordinates, int blockDim) 
{  
  int fortran_order = futils_is_called_from_fortran(); 
  int dimToSlice = n-1;

  if (fortran_order == 1) {
    dimToSlice = blockDim  - n; // n is from 1-(bb->ndim)
  }

  //  log_debug("getCoordinateFromBlock: pos = %lld, n=%d, fortran_order? %d, dimToSlice %d\n", pos, n, fortran_order, dimToSlice);
  if (n == 1) {
      coordinates[n-1] = pos + sel->start[dimToSlice];
      //log_debug("     coordinate [%d]=%lld\n", n-1, coordinates[n-1]);
      return ;
  } 

  uint64_t lastDimSize= sel->count[dimToSlice];     
  uint64_t res  = pos % lastDimSize;  // time consuming 25%

  coordinates[n-1] = res + sel->start[dimToSlice];
  uint64_t stepUp = (pos - res)/lastDimSize;

  //log_debug("      coordinate [%d]=%lld\n", n-1, coordinates[n-1]);
  getCoordinateFromBlock(stepUp, sel, n-1, coordinates, blockDim);  
}
*/
//check point coordinates
//offset in the bounding box needs to be taken account
 /*
void getCoordinateFromBox(uint64_t pos, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* sel, int n, uint64_t* coordinates) 
{
  int fortran_order = futils_is_called_from_fortran();
  int dimToSlice = n-1;

  if (fortran_order == 1) {
    dimToSlice = sel->ndim  - n; // n is from 1-(bb->ndim)
  }

  //log_debug("pos = %lld, n=%d \n", pos, n);
  if (n == 1) {
      coordinates[n-1] = pos + sel->start[dimToSlice];
      log_debug("                                       coordinate[%d]=%lld\n", n-1, coordinates[n-1]);
      return ;
  } 
 
  uint64_t lastDimSize= sel->count[dimToSlice];     

  uint64_t res  = pos % lastDimSize;
  log_debug("       \t\t\t\t   pos = %lld lastDim = %lld, dimToSlice=%d n=%d res=%lld, start=%lld ", pos, lastDimSize, dimToSlice, n, res, sel->start[dimToSlice]);
  //log_debug("      lastDim = %lld, res=%d \n", lastDimSize, res);
  coordinates[n-1] = res + sel->start[dimToSlice];
  uint64_t stepUp = (pos - res)/lastDimSize;

  log_debug("      coordinate[%d]=%lld\n", n-1, coordinates[n-1]);

  getCoordinateFromBox(stepUp, sel, n-1, coordinates);  
}
 */
void getCoordinateFromVariable(uint64_t pos, const ADIOS_VARINFO* var, int n, uint64_t* coordinates) 
{
  if (n == 1) {
      coordinates[n-1] = pos;
      return ;
  } 
 
  uint64_t lastDimSize= var->dims[n-1];
  uint64_t res  = pos % lastDimSize;

  coordinates[n-1] = res;
  uint64_t stepUp = (pos - res)/lastDimSize;

  getCoordinateFromVariable(stepUp, var, n-1, coordinates);  
}


/* static ADIOS_VARINFO* getAdiosVariable(ADIOS_FILE* f, const char* varName)
{
  ADIOS_VARINFO * v = common_read_inq_var(f, varName);

  if (v != NULL) {
    log_debug("  found variable [%s] in file\n", varName);
    return v;
  } 

  return NULL;
}*/



//
// Initialization is private business, Finalize is called by ADIOS when read is finalized
//
static int is_method_initialized = 0;
static void adios_query_fastbit_init() 
{
  const char* conffile = 0;
#ifdef DEBUG
  int msglvl = 200;
#else
  int msglvl = 0;
#endif
  if (!is_method_initialized) {
      fastbit_init(conffile);
      fastbit_set_verbose_level(msglvl);
      log_debug("[fastbit has initialized with verbosity level = %d]\n", msglvl);
      is_method_initialized = 1;
  }
}

int adios_query_fastbit_finalize()
{
  if (is_method_initialized) {
      fastbit_iapi_free_all();
      fastbit_cleanup();
      is_method_initialized = 0;
  }
  return 0;
}

 

int64_t getPosInBox(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* sel, int n, uint64_t* spatialCoordinates, int fortran_order) 
{
  if (sel->ndim <= 0) {
    return -1;
  }

  //int fortran_order = futils_is_called_from_fortran();

  int i=0;
  if (n == sel->ndim) { // check validation once
    for (i=0; i<sel->ndim; i++) {
      int matchingBoxDim = i;
      if (fortran_order == 1) {
	matchingBoxDim = sel->ndim - 1 - i;
      }
      uint64_t min = sel->start[matchingBoxDim];
      uint64_t max = sel->count[matchingBoxDim] + min;
      if (spatialCoordinates[i] >= max) {
	return -1;
      }
      if (spatialCoordinates[i] < min) {
	return -1;
      } 
    }
  }

  // spatial Coordinate is valid

  if (fortran_order == 1) {
    int matchingBoxDim = sel->ndim  - n;
     if (n == 1) {
       //log_debug("n=1, c[0]=%lld ,  start[%d]=%lld", spatialCoordinates[0], matchingBoxDim, sel->start[matchingBoxDim]);
       return spatialCoordinates[0] - sel->start[matchingBoxDim];
     }
     //log_debug("n=%d, c[n-1]=%lld, start/count[%d]= %lld/%lld \n", n, spatialCoordinates[n-1], matchingBoxDim, sel->start[matchingBoxDim], sel->count[matchingBoxDim]);
     return (spatialCoordinates[n-1]- sel->start[matchingBoxDim]) + sel->count[matchingBoxDim] * getPosInBox(sel, n-1, spatialCoordinates, fortran_order);
  } else {
     if (n == 1) {
        return spatialCoordinates[0]-sel->start[0];    
     }
     return (spatialCoordinates[n-1]-sel->start[n-1]) + sel->count[n-1]*getPosInBox(sel, n-1, spatialCoordinates, fortran_order);
  }
}

int64_t getPosInVariable(const ADIOS_VARINFO* v, int n, uint64_t* spatialCoordinates, int fortran_order) 
{
  if (v->ndim <= 0) {
    return -1;
  }

  //  log_debug("getPosInVariables() v->dim[0]=%d sp[%d]=%lld\n", v->dims[0], n-1, spatialCoordinates[n-1]);

  if (n == 1) {
      return spatialCoordinates[0];
  }

  if (fortran_order == 1) {
    int matchingBoxDim = v->ndim  - n;
    return (spatialCoordinates[n-1]) + v->dims[matchingBoxDim] * getPosInVariable(v, n-1, spatialCoordinates, fortran_order);
  } else {
    return  spatialCoordinates[n-1] + v->dims[n-1]*getPosInVariable(v, n-1, spatialCoordinates, fortran_order); 
  }
}


int isSameRegion(const ADIOS_VARINFO* v, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb, int timestep)
{
  if (v->sum_nblocks > 1) {
    return 0;
  }

  if (v->nsteps > 1) { // non streaming case
    if (v->nblocks[timestep] > 1)  // bb = block for this timestep
      return 0;
  }

  if (bb == NULL) {
    return 1;
  }

  int i=0;
  for (i=0; i<v->ndim; i++) {
    if (bb->start[i] != 0) {
      return 0;
    }      
    if (v->dims[i] != bb->count[i]) {
      return 0;
    }
  }
  return 1;
}

//
// absBlockIdx is needed as this is what blockinfo[] takes to get measure of a block (in order to coordinates from a given 1-d pos)
//
int64_t getRelativeIdx(uint64_t currPosInBlock,  const ADIOS_VARINFO* v, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb, int absBlockIdx, int timestep)
{
  ADIOS_VARBLOCK* blockSel = &(v->blockinfo[absBlockIdx]);

  // only calculates when more than one block presented
    int isFortranClient = futils_is_called_from_fortran();
    if (bb == NULL) {
        uint64_t spatialCoordinates[v->ndim];
	posToSpace(currPosInBlock, isFortranClient, blockSel->count, spatialCoordinates, v->ndim, blockSel->start);
	return getPosInVariable(v, v->ndim, spatialCoordinates, isFortranClient);      
    } else {      
       uint64_t spatialCoordinates[bb->ndim];        
       posToSpace(currPosInBlock, isFortranClient, blockSel->count, spatialCoordinates, bb->ndim, blockSel->start);
       return getPosInBox(bb, bb->ndim, spatialCoordinates, isFortranClient);
    }
}

/*
int64_t getRelativeIdxInBoundingBox(uint64_t currPosInBlock, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb, const ADIOS_VARBLOCK* blockSel)
{
    uint64_t spatialCoordinates[bb->ndim];    
    //getCoordinateFromBlock(currPosInBlock, blockSel, bb->ndim, spatialCoordinates, bb->ndim);
    
    int isFortranClient = futils_is_called_from_fortran();
    posToSpace(currPosInBlock, isFortranClient, blockSel->count, spatialCoordinates, bb->ndim, blockSel->start);
    return getPosInBox(bb, bb->ndim, spatialCoordinates, isFortranClient);
}

int64_t getRelativeIdxInVariable(uint64_t currPosInBlock, const ADIOS_VARINFO* v, const ADIOS_VARBLOCK* blockSel)
{
    uint64_t spatialCoordinates[v->ndim];
    //getCoordinateFromBlock(currPosInBlock, blockSel, v->ndim, spatialCoordinates, v->ndim);
    int isFortranClient = futils_is_called_from_fortran();
    posToSpace(currPosInBlock, isFortranClient, blockSel->count, spatialCoordinates, v->ndim, blockSel->start);

    return getPosInVariable(v, v->ndim, spatialCoordinates, isFortranClient);
}
*/
//
// the calculated block start/end  are relative to given timestep
//
int getBlockCoveredOnBB(int* blockStart, int* blockEnd, uint64_t* sumBlocksBeforeThisTimeStep, ADIOS_QUERY* q, int timeStep)
{
    ADIOS_SELECTION* sel = q->sel;
    ADIOS_VARINFO* v = q->varinfo;

    int i=0;
    const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = NULL;
    if (sel == NULL) {
        *blockStart=0;
	if (v->nsteps == 1) { // could be streaming open
	  *blockEnd = v->nblocks[0]-1;
	} else {
	  *blockEnd = v->nblocks[timeStep] -1;      
	}
	log_debug(" got blockStart = %d, blockEnd = %d\n", *blockStart, *blockEnd);
    } else {
        bb = &(sel->u.bb);      
	if (v->blockinfo == NULL) {
	  common_read_inq_var_blockinfo(q->file, v);
	}
	*blockStart = fastbit_adios_util_getRelativeBlockNumForPoint(v,bb->start,timeStep);
	uint64_t end[v->ndim];
	for (i=0; i<v->ndim; i++) {
	  end[i] = bb->start[i]+bb->count[i]-1;
	}
	*blockEnd = fastbit_adios_util_getRelativeBlockNumForPoint(v, end, timeStep);
	log_debug(" figured blockStart = %d, blockEnd = %d\n", *blockStart, *blockEnd);
    }

    if ((*blockStart < 0) || (*blockEnd < 0) || (*blockStart > *blockEnd)) {
        adios_error (err_invalid_query_value, "Query processing failed. Unable to continue using index. vid=%d, dim[0]=%lld\n", v->varid, v->dims[0]);
	return -1;
    }
    casestudyLogger_setPrefix(" computed block ids to scan");

    if (v->nsteps > 1) { // non streaming case
      for (i=0;i<timeStep; i++) {
	*sumBlocksBeforeThisTimeStep += v->nblocks[i];
      }
    }

    return 0;
} // getBlockCoveredOnBB


void getIncrement(uint64_t absBlockIdx, uint64_t* increment, int packSize, int timeStep, ADIOS_VARINFO* v)
{
  int i=0;
  for (i=0; i<packSize; i++) {
    uint64_t blockSize = fastbit_adios_util_getBlockSize(v, timeStep, absBlockIdx+i);
    if (i == 0) {
      increment[i] = blockSize;
    } else {
      increment[i] = blockSize + increment[i-1];
    }
  }
}

void locateBlockFromPack(uint64_t currPosInPack, uint64_t* absBlockIdx, uint64_t* posInBlock, uint64_t* packIncrementByBlock, int packSize)
{
  if (currPosInPack > packIncrementByBlock[packSize-1]) {
      *absBlockIdx=-1;
      *posInBlock = -1;
      printf(" ERROR: %llu th element will not be found  in pack starting at block:%llu", currPosInPack, *absBlockIdx);
  }

  uint64_t packStarts = *absBlockIdx;
  int i=0;
  for (i=0; i<packSize; i++) {    
    if (currPosInPack < packIncrementByBlock[i]) {
      *absBlockIdx = i + packStarts;
      *posInBlock = packIncrementByBlock[i] - currPosInPack;
      break;
    }
  }
}

int recursiveIdentify(uint64_t* packIncrementByBlock, uint64_t currPosInPack, int start, int end)
{
  //printf(".... start=%d end = %d, currpos=%llu \n", start, end, currPosInPack);
  if (end - start == 1) {
    if (currPosInPack > packIncrementByBlock[start]) {
      //printf("     ===> %llu\n", end);
      return end;
    } else {
      //printf("     ===> %llu\n", start);
      return start;
    }
  }

  int pos = (start+end)/2;

  if (currPosInPack > packIncrementByBlock[pos]) {
    return recursiveIdentify(packIncrementByBlock, currPosInPack, pos, end);
  } else {
    return recursiveIdentify(packIncrementByBlock, currPosInPack, start, pos);
  }

}
int identifyBlockInPack(uint64_t* packIncrementByBlock, uint64_t currPosInPack, int packSize)
{
#ifdef SIZEN
  int i=0;
  for (i=0; i<packSize; i++) {
    if (currPosInPack < packIncrementByBlock[i]) {
      return i;
    }
  }
#else
  return recursiveIdentify(packIncrementByBlock, currPosInPack, 0, packSize);
#endif
  return -1;
}
// coordinateArray is hits within the whole pack
// this function sorts the hits with actual blocks
// say a pack has N blocks, relative block ids are from 0 to (N-1), 
// since hits are sorted, so we figure out start-a count-a rests in block-a
// etc
void sort(uint64_t* coordinateArray, uint64_t arraySize,  uint64_t* packIncrementByBlock, uint64_t* starts, uint64_t* counts, int* relativeBlockIds, int packSize)
{
  int i=0;
  int whichBlock;
  uint64_t arrayCounter = 0;

  for (i=0; i<packSize; i++) {
    relativeBlockIds[i]  = -1;
  }
  
  if (arraySize == 0) {
    return;
  }

  relativeBlockIds[0] = identifyBlockInPack(packIncrementByBlock, coordinateArray[0], packSize);
  starts[0] = 0;

  int lastBlock = identifyBlockInPack(packIncrementByBlock, coordinateArray[arraySize -1], packSize);

  if (relativeBlockIds[0] == lastBlock) {
      counts[0] = arraySize;
      return;
  }
  
  int bcounter =0 ; 
  while (true) {
    for (arrayCounter=starts[bcounter]+1; arrayCounter<arraySize; arrayCounter++) {
      //int bid = recursiveIdentify(packIncrementByBlock, coordinateArray[arrayCounter], relativeBlockIds[0], packSize);
	int bid = identifyBlockInPack(packIncrementByBlock, coordinateArray[arrayCounter], packSize);
        //int bid = coordinateArray[arrayCounter]/packIncrementByBlock[0];

      if (bid > relativeBlockIds[bcounter]) {
	 counts[bcounter] = arrayCounter-starts[bcounter];
	 bcounter++;
	 starts[bcounter] = arrayCounter;
	 relativeBlockIds[bcounter] = bid;
	 break;
      }
    }

    if (relativeBlockIds[bcounter] == lastBlock) {
      counts[bcounter] = arraySize - starts[bcounter];
      break;
    }
  }
}


void loopThrough(uint64_t dimMultiplier, uint64_t accumulatedPos, int currDim, ADIOS_VARINFO* var, uint64_t* regionStart, uint64_t* regionCount, uint32_t* bbSlice)
{
  int i = 0;
  if (currDim == var->ndim-1) {
    for (i=0; i<regionCount[currDim]; i++) {
      uint64_t pos = accumulatedPos + (i+regionStart[currDim]);
      bitarray_setbit(bbSlice, pos);
    }
    return;
  } 

  for (i=0; i<regionCount[currDim]; i++) {
      uint64_t m = dimMultiplier/var->dims[currDim];
      uint64_t pos = accumulatedPos + (i+regionStart[currDim]) * m;
      loopThrough(m, pos, currDim+1, var, regionStart, regionCount, bbSlice);
  }
}

//
// create bit slice that mark 1 if within BB, 0 elsewise
//
uint32_t* bitarray_create_markBB(uint64_t knownSize, ADIOS_VARINFO* var, uint64_t* regionStart, uint64_t* regionCount)
{
  uint32_t* bbSlice = bitarray_create(knownSize);

  if (var->ndim == 1) {
    int i=0;
    for (i = 0; i<regionCount[0]; i++) {
      uint64_t pos = (i+regionStart[0]);
      bitarray_setbit(bbSlice, pos);
    }
    return bbSlice;
  }

  if (var->ndim == 2) {
    int i=0, j=0;
    for (i = 0; i<regionCount[0]; i++) {
      for (j=0; j<regionCount[1]; j++) {
	uint64_t pos = (i+regionStart[0]) * var->dims[1] + (j+regionStart[1]);
	bitarray_setbit(bbSlice, pos);
      }
    }
    return bbSlice;
  }

  /*
  if (bb->ndim == 3) {
    int x=0, y=0, z=0;
    for (x=0; x<bb->count[0]; x++) {
      for (y=0; y<bb->count[1]; y++) {
	for (z=0; z<bb->count[2]; z++) {
	  uint64_t pos = (x+bb->start[0]) * var->dims[1]*var->dims[2] + (y+bb->start[1]) * var->dims[2] + (z+bb->start[2]);
	  bitarray_setbit(bbSlice, pos);
	}
      }     
    }
    return bbSlice;
  }
  */
  int i=0;
  uint64_t dimOneToN=1;
  for (i=1; i<var->ndim; i++) {
    dimOneToN *= var->dims[i];
  }

  for (i = 0; i<regionCount[0]; i++) {
    uint64_t pos = (i+regionStart[0]) * dimOneToN;
    loopThrough(dimOneToN, pos, 1, var, regionStart, regionCount, bbSlice);
  }

  return bbSlice;
}


void loopThroughAndAssign(uint64_t varDimOneToN, uint64_t bbDimOneToN, uint64_t rpos_sofar, uint64_t pos_sofar, int currDim, ADIOS_VARINFO* var, 
			  uint64_t* regionStart, uint64_t* regionCount,  uint32_t* bitSlice, uint32_t* resultSlice, uint64_t knownSize)
{
  int i = 0;
  if (currDim == var->ndim-1) {
    for (i=0; i<regionCount[currDim]; i++) {
      uint64_t rpos = rpos_sofar + (i+regionStart[currDim]);
      uint64_t pos =  pos_sofar + i;

      if ((rpos < knownSize) && BITTEST(resultSlice, rpos)) {
	bitarray_setbit(bitSlice, pos);
      }
    }
    return;
  } 

  for (i=0; i<regionCount[currDim]; i++) {
      uint64_t mv = varDimOneToN/(var->dims[currDim]);
      uint64_t mb = bbDimOneToN/(regionCount[currDim]);
      uint64_t rpos = rpos_sofar + (i+regionStart[currDim]) * mv;
      uint64_t pos  = pos_sofar + (i) * mb;
      loopThroughAndAssign(mv, mb, rpos, pos, currDim+1, var, regionStart, regionCount, bitSlice, resultSlice, knownSize);
  }
}

void bitMapBack(uint32_t* bitSlice, uint32_t* resultSlice, ADIOS_VARINFO* var, uint64_t* regionStart, uint64_t* regionCount, uint64_t knownSize)
{
  int i=0;
  if (var->ndim == 1) {
    for (i=0; i<regionCount[0]; i++) {
      uint64_t rpos = i+regionStart[0];
      uint64_t pos = i;
      if ((rpos < knownSize) && BITTEST(resultSlice, rpos)) {
	bitarray_setbit(bitSlice, pos);
      }
    }
    return;
  } 

  uint64_t varDimOneToN=1;
  uint64_t bbDimOneToN = 1;
  for (i=1; i<var->ndim; i++) {
    varDimOneToN *= var->dims[i];
    bbDimOneToN *= regionCount[i];    
  }

  for (i = 0; i<regionCount[0]; i++) {
    uint64_t rpos = (i+regionStart[0]) * varDimOneToN;
    uint64_t pos = i * bbDimOneToN;
    loopThroughAndAssign(varDimOneToN, bbDimOneToN, rpos, pos, 1, var, regionStart, regionCount, bitSlice, resultSlice, knownSize);
  }
}

 //
 //
void setBitArray(ADIOS_VARINFO* var, uint32_t* bitSlice, int64_t* coordinateArray, uint64_t count, uint64_t adjustment,
		 uint64_t* regionStart, uint64_t* regionCount, uint64_t eleStarts, uint64_t eleEnds)
{
  uint64_t knownSize = eleEnds+1; // at most this many elements. + 1 is due to C arrays starts at 0. 

  uint32_t* hitsSlice = bitarray_create(knownSize);
  uint32_t* resultSlice = bitarray_create(knownSize);

  uint32_t* bbSlice = bitarray_create_markBB(knownSize, var, regionStart, regionCount);
 
  int k=0;
  for (k=0; k<count; k++) {
    int64_t currPosInVar = coordinateArray[k] + adjustment;
    if (currPosInVar < knownSize) {
      bitarray_setbit(hitsSlice, currPosInVar);
    }
  }

  bitarray_and(hitsSlice, bbSlice,  resultSlice, BITNSLOTS(knownSize));      
 
  // ok
  int n = bitarray_countHits(resultSlice, BITNSLOTS(knownSize));      

  
  free(bbSlice);
  free(hitsSlice);

  bitMapBack(bitSlice, resultSlice, var, regionStart, regionCount, knownSize);

  free(resultSlice);
}


void setBitArray0(ADIOS_VARINFO* var, uint32_t* bitSlice, int64_t* coordinateArray, uint64_t count, uint64_t adjustment,
		 const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb, uint64_t eleStarts, uint64_t eleEnds)
{
  int isFortranClient = 0;

  int k=0;
  for (k=0; k<count; k++) {
    int64_t currPosInBlock = coordinateArray[k] + adjustment;
      
    if ((currPosInBlock >= eleStarts) && (currPosInBlock <= eleEnds)) {	  
        uint64_t spatialCoordinates[bb->ndim];
	getCoordinateFromVariable(currPosInBlock, var, var->ndim, spatialCoordinates);
	//posToSpace(currPosInBlock, isFortranClient, bb->count, spatialCoordinates, bb->ndim, bb->start);
	int64_t p = getPosInBox(bb, bb->ndim, spatialCoordinates, isFortranClient);
	if (p >= 0) {
	  bitarray_setbit(bitSlice, p);
	} 
    }
  }
}

void checkHits(ADIOS_VARINFO* v, ADIOS_QUERY* q, uint64_t boxStart, uint64_t* regionStart, uint64_t* regionCount, uint64_t eleStarts, uint64_t eleEnds)
{
      uint64_t count = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 	
      //uint64_t  coordinateArray[count];				
      uint64_t* coordinateArray = malloc(count*sizeof(uint64_t));
      fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, count, 0);      
      
      // set bits
      uint32_t* bitSlice = bitarray_create(q->rawDataSize);
      setBitArray(v, bitSlice, coordinateArray, count, boxStart, regionStart, regionCount, eleStarts, eleEnds);

      free(q->dataSlice);
      q->dataSlice = bitSlice;
      
      fastbit_iapi_free_array_by_addr(q->dataSlice);
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = 0;

      casestudyLogger_setPrefix(" summarized evaluation for bb");  
      free(coordinateArray);
}

int mEvaluateBBRangeFancyQueryOnWhole(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep, uint64_t* regionStart, uint64_t* regionCount)
{
  char bitsArrayName[60];
  sprintf(bitsArrayName, "%ld_%d", fastbit_adios_getCurrentTimeMillis(), timeStep);

  uint64_t s = 0;
  uint64_t startRef = 0;

  ADIOS_VARINFO* v = getFirstLeaf(q)->varinfo;

  uint64_t eleStarts = 0;
  uint64_t eleEnds = q->rawDataSize;

  uint64_t totalEle = 1;

  if ((regionStart == NULL) || (regionCount == NULL)) {
    return mEvaluateTestFullRangeFancyQueryOnWhole(idxFile, q, timeStep);
  } 

  eleStarts = getPosInVariable(v, v->ndim, regionStart, 0);
  uint64_t end[v->ndim];
  int coversAll = 1;
  int i = 0;
  for (i=0; i<v->ndim; i++) {
    if (regionCount[i] > 0) {
      end[i] = regionStart[i]+regionCount[i]-1;
    } else {
      end[i] = regionStart[i];
    }

    totalEle *= v->dims[i];
    if (v->dims[i] > regionCount[i]) {
      coversAll = 0;
    }
  }
  
  if (coversAll == 1) {
    return mEvaluateTestFullRangeFancyQueryOnWhole(idxFile, q, timeStep);
  }  
  eleEnds = getPosInVariable(v, v->ndim, end, 0);      
  

  ADIOS_VARINFO* packVar = common_read_inq_var (idxFile, "elements");
  uint64_t recommended_index_ele = 0;
  common_read_schedule_read (idxFile, NULL, "elements",           0, 1, NULL, &recommended_index_ele);
  common_read_perform_reads(idxFile, 1);

  //printf("dataSize from: %llu to %llu, elements = %lu\n", eleStarts, eleEnds, recommended_index_ele);

  uint64_t start[v->ndim], count[v->ndim];
  uint64_t split = totalEle/recommended_index_ele;
  uint32_t* bitSlice = NULL;

  int64_t eleBoxStarts = -1; 
  if (split == 0) {
      // index is on the whole timestep
      getHandle(timeStep, 0, idxFile,  q,  totalEle);
      checkHits(v, q, 0, regionStart, regionCount, eleStarts, eleEnds); 
      return 0;
  } else {
      int boxCounter = 0;
      while (startRef < v->dims[0]) {
	uint64_t count[v->ndim];
	uint64_t boxSize = 1;
	
	start[0] = startRef;
	startRef += v->dims[0]/split;
	if (startRef >= v->dims[0]) {
	  startRef = v->dims[0];
	}
	count[0] = startRef - start[0];
	boxSize *=count[0];
	for (s=1; s < v->ndim; s++) {
	  start[s] = 0;
	  count[s] = v->dims[s];
	  boxSize *= count[s];
	}

	if ((regionCount[0]+regionStart[0] >= start[0]) && (start[0]+count[0] > regionStart[0])) 
	{
	  if (eleBoxStarts == -1) {
	    eleBoxStarts = getPosInVariable(v, v->ndim, start, 0);
	  }
	  getHandle(timeStep, start[0], idxFile, q, boxSize);	  
	
	  // has to call evaluate before calling either register_selection_as_bit_array() or extend_bit_array()
	  uint64_t countMe = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 
	  
	  if (boxCounter == 0) {
	    fastbit_iapi_register_selection_as_bit_array(bitsArrayName, ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
	  } else {
	    fastbit_iapi_extend_bit_array_with_selection(bitsArrayName, ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
	  }
	  clean_fastbit_trace(q);	  
	}
      
	boxCounter++;
      }
  }

  FastBitSelectionHandle h  =  fastbit_selection_osr(bitsArrayName, FastBitCompareGreater, 0.0);//createHandle(q, bitsArrayName);

  fastbit_adios_util_checkNotNull(h, bitsArrayName);    
  ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = h;

  checkHits(v, q, eleBoxStarts, regionStart, regionCount, eleStarts, eleEnds); 
}



void checkHitsDefault(ADIOS_QUERY* q)
{
  uint64_t resultCount = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 	

  //uint64_t  coordinateArray[resultCount];
  uint64_t* coordinateArray = malloc(resultCount*sizeof(uint64_t));
  fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, resultCount, 0);      
  casestudyLogger_setPrefix(" got coordinates bb");

  uint32_t* bitSlice = bitarray_create(q->rawDataSize);

  int k;
  for (k=0; k<resultCount; k++) {
    int64_t currPosInBlock = coordinateArray[k];
    if (currPosInBlock >= 0) {
      bitarray_setbit(bitSlice, currPosInBlock);
    }
  } 
  
  free(q->dataSlice);
  q->dataSlice = bitSlice;
  
  fastbit_iapi_free_array_by_addr(q->dataSlice);
  ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = 0;

  casestudyLogger_setPrefix(" summarized evaluation for bb");  
  free(coordinateArray);
}
 //
 // 
 //

int mEvaluateTestFullRangeFancyQueryOnWhole(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{
  char bitsArrayName[60];
  sprintf(bitsArrayName, "%ld_%d", fastbit_adios_getCurrentTimeMillis(), timeStep);

  uint64_t s = 0;
  uint64_t startRef = 0;

  ADIOS_VARINFO* packVar = common_read_inq_var (idxFile, "elements");
  uint64_t recommended_index_ele = 0;
  common_read_schedule_read (idxFile, NULL, "elements",           0, 1, NULL, &recommended_index_ele);
  common_read_perform_reads(idxFile, 1);

  uint64_t dataSize = q->rawDataSize;

  uint64_t start[getFirstLeaf(q)->varinfo->ndim], count[getFirstLeaf(q)->varinfo->ndim];
  uint64_t split = q->rawDataSize/recommended_index_ele;
  //uint32_t* bitSlice = bitarray_create(q->rawDataSize);
  uint32_t* bitSlice = NULL;
  ADIOS_VARINFO* v = getFirstLeaf(q)->varinfo;

  if (split == 0) {
      // index is on the whole timestep
      getHandle(timeStep, 0, idxFile,  q,  dataSize);

      checkHitsDefault(q);
      return 0;
  } else {
      int boxCounter = 0;
      while (startRef < v->dims[0]) {
	uint64_t count[v->ndim];
	uint64_t boxSize = 1;
	
	start[0] = startRef;
	startRef += v->dims[0]/split;
	if (startRef >= v->dims[0]) {
	  startRef = v->dims[0];
	}
	count[0] = startRef - start[0];
	boxSize *=count[0];
	for (s=1; s < v->ndim; s++) {
	  start[s] = 0;
	  count[s] = v->dims[s];
	  boxSize *= count[s];
	}
	
	//getHandleFromBlockAtLeafQuery(timeStep, start[0], idxFile,  q,  boxSize);
	//
	getHandle(timeStep, start[0], idxFile, q, boxSize);
	
	// has to call evaluate before calling either register_selection_as_bit_array() or extend_bit_array()
	uint64_t countMe = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 
	//printf("  boxCounter=%d, countme=%ld \n", boxCounter, countMe);
	if (boxCounter == 0) {
	  fastbit_iapi_register_selection_as_bit_array(bitsArrayName, ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
	} else {
	  fastbit_iapi_extend_bit_array_with_selection(bitsArrayName, ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
	}

	clean_fastbit_trace(q);
	boxCounter++;
      }
  }

  FastBitSelectionHandle h  =  fastbit_selection_osr(bitsArrayName, FastBitCompareGreater, 0.0);//createHandle(q, bitsArrayName);
  fastbit_adios_util_checkNotNull(h, bitsArrayName);    
  ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = h;

  checkHitsDefault(q);
  return 0;
}


 // with fancy query, on leaf node
int mEvaluateTestFullRange(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{

  char bitsArrayName[60];
  sprintf(bitsArrayName, "%ld_%d_%d", fastbit_adios_getCurrentTimeMillis(), q->varinfo->varid, timeStep);

  uint64_t s = 0;
  uint64_t dataSize = 1;
  uint64_t start[q->varinfo->ndim], count[q->varinfo->ndim];
  for (s=0; s<q->varinfo->ndim; s++) {
    dataSize *= q->varinfo->dims[s];
  }

  uint64_t startRef = 0;

  ADIOS_VARINFO* packVar = common_read_inq_var (idxFile, "elements");
  uint64_t recommended_index_ele = 0;
  common_read_schedule_read (idxFile, NULL, "elements",           0, 1, NULL, &recommended_index_ele);
  common_read_perform_reads(idxFile, 1);

  //printf("dataSize = %llu, elements = %lu\n", dataSize, recommended_index_ele);

  uint64_t split = dataSize/recommended_index_ele;
  //uint32_t* bitSlice = bitarray_create(q->rawDataSize);
  uint32_t* bitSlice = NULL;

  if (split == 0) {
      // index is on the whole timestep
      getHandleFromBlockAtLeafQuery(timeStep, 0, idxFile,  q,  q->rawDataSize);
      uint64_t count = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 	
      //uint64_t  coordinateArray[count];
      uint64_t* coordinateArray = malloc(count*sizeof(uint64_t));
      fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, count, 0);      
      
      int k=0;
      // set bits
      bitSlice = bitarray_create(q->rawDataSize);
      for (k=0; k<count; k++) {
	int64_t currPosInBlock = coordinateArray[k];
	if (currPosInBlock >= 0) {
	  bitarray_setbit(bitSlice, currPosInBlock);
	}
      }

      free(q->dataSlice);
      q->dataSlice = bitSlice;
      
      fastbit_iapi_free_array_by_addr(q->dataSlice);
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = 0;
      
      casestudyLogger_setPrefix(" summarized evaluation for bb");  
      free(coordinateArray);
      return 0;
  } else {
      int boxCounter = 0;
      while (startRef < q->varinfo->dims[0]) {
	uint64_t count[q->varinfo->ndim];
	uint64_t boxSize = 1;
	
	start[0] = startRef;
	startRef += q->varinfo->dims[0]/split;
	if (startRef >= q->varinfo->dims[0]) {
	  startRef = q->varinfo->dims[0];
	}
	count[0] = startRef - start[0];
	boxSize *=count[0];
	for (s=1; s < q->varinfo->ndim; s++) {
	  start[s] = 0;
	  count[s] = q->varinfo->dims[s];
	  boxSize *= count[s];
	}
	
	getHandleFromBlockAtLeafQuery(timeStep, start[0], idxFile,  q,  boxSize);
	/*
	// index is on box identified by start/count
	uint64_t resultCount = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 	
	int64_t  coordinateArray[resultCount];
	fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, resultCount, 0);      
	
	int64_t startPosInVar = getPosInVariable(q->varinfo, q->varinfo->ndim, start, 0);
	printf("result count=%lu, startPos=%lld, (%llu, %llu, %llu) \n", resultCount, startPosInVar, start[0], start[1], start[2]);
	
	int k=0;
	// set bits
	for (k=0; k<resultCount; k++) {
	  int64_t currPosInBlock = coordinateArray[k]+startPosInVar;
	  if (currPosInBlock >= 0) {
	    bitarray_setbit(bitSlice, currPosInBlock);
	  }
	} 
	*/
	
	// has to call evaluate before calling either register_selection_as_bit_array() or extend_bit_array()
	uint64_t countMe = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 

	if (boxCounter == 0) {
	  fastbit_iapi_register_selection_as_bit_array(bitsArrayName, ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
	} else {
	  fastbit_iapi_extend_bit_array_with_selection(bitsArrayName, ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
	}

	clean_fastbit_trace(q);
      
	boxCounter++;
      }
  }

  FastBitSelectionHandle h  =  fastbit_selection_osr(bitsArrayName, FastBitCompareGreater, 0.0);//createHandle(q, bitsArrayName);

  fastbit_adios_util_checkNotNull(h, bitsArrayName);    
  ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = h;

  uint64_t resultCount = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 	
  //printf("resultCount = %llu\n", resultCount);
  //  uint64_t  coordinateArray[resultCount];
  uint64_t* coordinateArray = malloc(resultCount*sizeof(uint64_t));
  fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, resultCount, 0);      
  casestudyLogger_setPrefix(" got coordinates bb");

  if (resultCount > q->rawDataSize/2) {
      int k;
      bitSlice = bitarray_create(q->rawDataSize);

      for (k=0; k<resultCount; k++) {
	int64_t currPosInBlock = coordinateArray[k];
	if (currPosInBlock >= 0) {
	  bitarray_setbit(bitSlice, currPosInBlock);
	}
      } 
      
  } else {
    bitSlice = bitarray_create(q->rawDataSize);
    int k;
    for (k=0; k<resultCount; k++) {
      int64_t currPosInBlock = coordinateArray[k];
      if (currPosInBlock >= 0) {
	bitarray_setbit(bitSlice, currPosInBlock);
      }
    } 
  }

  free(coordinateArray);
  free(q->dataSlice);
  q->dataSlice = bitSlice;
  
  fastbit_iapi_free_array_by_addr(q->dataSlice);
  
  ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = 0;
  
  casestudyLogger_setPrefix(" summarized evaluation for bb");
  
}


//
// idx can be based on (m)ultiple blocks
//
int mEvaluateWithIdxOnBBoxWithBitArrayOnVar(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{
    int blockStart=0; // relative
    int blockEnd = 0; // relative

    uint64_t sumBlocksBeforeThisTimeStep = 0;

    int i=0;

    if (getBlockCoveredOnBB(&blockStart, &blockEnd, &sumBlocksBeforeThisTimeStep, q, timeStep) != 0) {
      return -1;
    }

    uint32_t* bitSlice = bitarray_create(q->rawDataSize);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);      
    char bitsArrayName[60];
    sprintf(bitsArrayName, "%ld_%d_%d_%d", fastbit_adios_getCurrentTimeMillis(), q->varinfo->varid, timeStep, rank);

    uint64_t currBlockIdx = blockStart;

    casestudyLogger_setPrefix(" start scanning");
    char casestudyLoggerPrefix[30];
    
    ADIOS_VARINFO* packVar = common_read_inq_var (idxFile, "pack");
    int packSize = -1;
    common_read_schedule_read (idxFile, NULL, "pack",           0, 1, NULL, &packSize);
    common_read_perform_reads(idxFile, 1);
        
    uint64_t virtualBlockStart = blockStart/packSize;
    uint64_t virtualBlockEnds = blockEnd/packSize;
    
    if (q->file->is_streaming == 1) {
      if (packSize > q->varinfo->sum_nblocks) {
	packSize = q->varinfo->sum_nblocks;
      }
    } else {
      if (packSize > q->varinfo->nblocks[timeStep]) {
	packSize = q->varinfo->nblocks[timeStep];
      }
    }

      
    const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = NULL;
    if (q->sel != NULL) {
      bb = &(q->sel->u.bb);      
    }	

    int k=0;
    uint64_t absBlockIdx = 0;

    for (currBlockIdx= virtualBlockStart; currBlockIdx <= virtualBlockEnds; currBlockIdx++) {
      //getHandle(timeStep, currBlockIdx*packSize, idxFile, q);	      
	absBlockIdx = currBlockIdx*packSize+sumBlocksBeforeThisTimeStep;	

	uint64_t packIncrementByBlock[packSize];
	if (q->file->is_streaming == 1) {
	  getIncrement(absBlockIdx, packIncrementByBlock, packSize, -1, q->varinfo);
	} else {
	  getIncrement(absBlockIdx, packIncrementByBlock, packSize, timeStep, q->varinfo);
	}

        getHandleFromBlockAtLeafQuery(timeStep, currBlockIdx*packSize, idxFile,  q, packIncrementByBlock[packSize-1]);

	if (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle == 0) {
	    log_warn(" Unable to construct fastbit query with NULL. Use _no_o idx method \n");
	    return -1;
	}
	
	sprintf(casestudyLoggerPrefix, "block:%llu", currBlockIdx);
	casestudyLogger_setPrefix(casestudyLoggerPrefix);
	
	struct timespec evalStartT; casestudyLogger_getRealtime(&evalStartT);	
	uint64_t count = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 	

	//uint64_t  coordinateArray[count];
	uint64_t* coordinateArray = malloc(count*sizeof(uint64_t));
	fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, count, 0);      

	casestudyLogger_pro_writeout(&evalStartT, "fastbitevaluated");	
	struct timespec frameStartT; casestudyLogger_getRealtime(&frameStartT);

 	clear_fastbit_internal_recursive(q);

	if (isSameRegion(q->varinfo, bb, timeStep)) { 
	    // bb covers all 
	    // no need to do mapping , to do: check
	} else {
	    uint64_t posInBlock;
	    int isFortranClient = futils_is_called_from_fortran();
	    if (isFortranClient) { // to do: check
	      for (k=0; k<count; k++) {
		uint64_t currPosInPack = coordinateArray[k];		
	        locateBlockFromPack(currPosInPack, &absBlockIdx, &posInBlock, packIncrementByBlock, packSize);
		coordinateArray[k] =  getRelativeIdx(posInBlock, q->varinfo, bb, absBlockIdx, timeStep); 	
	      }
	    } else { // fastscan on c-clients
	      uint64_t starts[packSize]; uint64_t counts[packSize]; int relativeBlockIds[packSize];
	      sort(coordinateArray, count, packIncrementByBlock, starts, counts, relativeBlockIds, packSize);

	      for (k=0; k<packSize; k++) {		
		if (relativeBlockIds[k] >= 0) {
		    uint64_t currIdx = absBlockIdx + relativeBlockIds[k];		  
		    //printf("starts[%d] = %llu counts[%d] = %llu, relativeIdx=%lld, hits=%llu  \n", k, starts[k], k, counts[k], relativeBlockIds[k], count);
		    ADIOS_VARBLOCK* blockSel = &(q->varinfo->blockinfo[currIdx]);
		    if (relativeBlockIds[k] > 0) {
		      int m = 0; 
		      for (m =starts[k]; m<starts[k]+counts[k]; m++) {
			//printf(" adjust coordinateArray[%ld],  %lld - %lld \n", m, coordinateArray[m], packIncrementByBlock[relativeBlockIds[k]-1]);
			coordinateArray[m] = coordinateArray[m] - packIncrementByBlock[relativeBlockIds[k]-1];
			
		      }
		    }
		    fastscan(coordinateArray, starts[k], counts[k], q->varinfo, blockSel->start, blockSel->count, bb);	
		} else {
		    break;
		}
	      }
	    }
	}
	// set bits
	for (k=0; k<count; k++) {
	  int64_t currPosInSel = coordinateArray[k];
	  if (currPosInSel >= 0) {
	    bitarray_setbit(bitSlice, currPosInSel);
	  }
	}
	
	//casestudyLogger_frame_writeout(&frameStartT, "block processed");
	casestudyLogger_frame_writeout(&frameStartT, "block processed");
	log_debug("----\n");
	free(coordinateArray);
    }
    
    casestudyLogger_setPrefix(" blocksProcessedIndividually!");

    //return fastbit_selection_create(dataType, dataOfInterest, dataSize, compareOp, &vv);

    free(q->dataSlice);
    q->dataSlice = bitSlice;

    fastbit_iapi_free_array_by_addr(q->dataSlice);

    ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = 0;
 
    casestudyLogger_setPrefix(" summarized evaluation for bb");
    //fastbit_adios_util_checkNotNull(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, bitsArrayName);
}
//
// for index that based on one block 
//
int evaluateWithIdxOnBBoxWithBitArrayOnVar0(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{
    int blockStart=0; // relative
    int blockEnd = 0; // relative
    int i=0;

    ADIOS_SELECTION* sel = q->sel;
    ADIOS_VARINFO* v = q->varinfo;

    const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = NULL;
    if (sel == NULL) {
        blockStart=0;
	if (v->nsteps == 1) { // could be streaming open
	  blockEnd = v->nblocks[0]-1;
	} else {
	  blockEnd = v->nblocks[timeStep] -1;      
	}
	log_debug(" got blockStart = %d, blockEnd = %d\n", blockStart, blockEnd);
    } else {
        bb = &(sel->u.bb);      
	if (v->blockinfo == NULL) {
	  common_read_inq_var_blockinfo(q->file, v);
	}
	blockStart = fastbit_adios_util_getRelativeBlockNumForPoint(v,bb->start,timeStep);
	uint64_t end[v->ndim];
	for (i=0; i<v->ndim; i++) {
	  end[i] = bb->start[i]+bb->count[i]-1;
	}
	blockEnd = fastbit_adios_util_getRelativeBlockNumForPoint(v, end, timeStep);
	log_debug(" figured blockStart = %d, blockEnd = %d\n", blockStart, blockEnd);
    }

    if ((blockStart < 0) || (blockEnd < 0) || (blockStart > blockEnd)) {
        adios_error (err_invalid_query_value, "Query processing failed. Unable to continue using index. vid=%d, dim[0]=%lld\n", v->varid, v->dims[0]);
	return -1;
    }

    casestudyLogger_setPrefix(" computed block ids to scan");

    uint32_t* bitSlice = bitarray_create(q->rawDataSize);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);      
    char bitsArrayName[60];
    sprintf(bitsArrayName, "%ld_%d_%d_%d", fastbit_adios_getCurrentTimeMillis(), v->varid, timeStep, rank);

    uint64_t currBlockIdx = blockStart;

    uint64_t sumBlocksBeforeThisTimeStep = 0;
    if (v->nsteps > 1) { // non streaming case
      for (i=0;i<timeStep; i++) {
	sumBlocksBeforeThisTimeStep += v->nblocks[i];
      }
    }

    casestudyLogger_setPrefix(" start scanning");
    char casestudyLoggerPrefix[30];
    
    for (currBlockIdx=blockStart; currBlockIdx <= blockEnd; currBlockIdx++) {
      getHandle(timeStep, currBlockIdx, idxFile, q, 0);	      
      if (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle == 0) {
	log_warn(" Unable to construct fastbit query with NULL. Use _no_o idx method \n");
	return -1;
      }

      sprintf(casestudyLoggerPrefix, "block:%llu", currBlockIdx);
      casestudyLogger_setPrefix(casestudyLoggerPrefix);

      struct timespec evalStartT;
      casestudyLogger_getRealtime(&evalStartT);

      uint64_t count = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 

      //uint64_t  coordinateArray[count];
      uint64_t* coordinateArray = malloc(count*sizeof(uint64_t));
      fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, count, 0);      

      casestudyLogger_pro_writeout(&evalStartT, "fastbitevaluated");

      struct timespec frameStartT;
      casestudyLogger_getRealtime(&frameStartT);
      clear_fastbit_internal_recursive(q);

      int absBlockIdx = currBlockIdx+sumBlocksBeforeThisTimeStep;
      int k=0;

      int isDefaultSelection = isSameRegion(v, bb, timeStep);
      /*
      for (k=0; k<count; k++) {
	uint64_t currPosInBlock = coordinateArray[k];
	int64_t currPos = currPosInBlock;
	
	if (!isDefaultSelection) {
	    currPos = getRelativeIdx(currPosInBlock, v, bb, absBlockIdx, timeStep);
	}

	log_debug("%lld th in block[%d],   =>  in actual box %lld  \n", currPosInBlock, absBlockIdx, currPos);
	if (currPos >= 0) {
            #ifdef BITARRAY
	    bitarray_setbit(bitSlice, currPos);
            #else
	    bitSlice[currPos] = 1;
            #endif
	}
      }
      */
      if (isDefaultSelection) {
	// no need to do mapping , to do: check
      } else {
	int isFortranClient = futils_is_called_from_fortran();
	if (isFortranClient) { // to do: check
	    for (k=0; k<count; k++) {
	      uint64_t currPosInBlock = coordinateArray[k];
	      coordinateArray[k] =  getRelativeIdx(currPosInBlock, v, bb, absBlockIdx, timeStep); 	
	    }
	} else { // fastscan on c-clients
	    ADIOS_VARBLOCK* blockSel = &(v->blockinfo[absBlockIdx]);
	    fastscan(coordinateArray, 0, count, v, blockSel->start, blockSel->count, bb);	
	}
      }

      // set bits
      for (k=0; k<count; k++) {
	int64_t currPosInBlock = coordinateArray[k];
	if (currPosInBlock >= 0) {
	  bitarray_setbit(bitSlice, currPosInBlock);
	}
      }
      
      //casestudyLogger_frame_writeout(&frameStartT, "block processed");
      casestudyLogger_pro_writeout(&frameStartT, "block processed");
      log_debug("----\n");
      free(coordinateArray);
    }

    casestudyLogger_setPrefix(" blocksProcessedIndividually!");

    //return fastbit_selection_create(dataType, dataOfInterest, dataSize, compareOp, &vv);

    free(q->dataSlice);
    q->dataSlice = bitSlice;

    fastbit_iapi_free_array_by_addr(q->dataSlice);

    ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = 0;
 
    casestudyLogger_setPrefix(" summarized evaluation for bb");
    //fastbit_adios_util_checkNotNull(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, bitsArrayName);
}


int isSameBB(ADIOS_QUERY* q, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb)
{
  if (q->left == NULL) {
    if (q->sel != NULL) {
      const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* cbb = &(q->sel->u.bb);
      int i=0;
      for (i=0; i<bb->ndim; i++) {
	if (bb->start[i] != cbb->start[i]) {
	  return 0;
	}
	if (bb->count[i] != cbb->count[i]) {
	  return 0;
	}
      }						
      return 1;
    }
    return 0;  // skip furthur investigation
  }
  if (isSameBB(q->left, bb) && isSameBB(q->right, bb)) {
    return 1;
  }
  
  return 0;
}
// roughly
int  isBasedOnSameBB(ADIOS_QUERY* q) {
  if (q->left == NULL) { // leaf query, ok
    return 1;
  }

  if (getFirstLeaf(q)->sel == NULL) {
    return 0;
  }

  const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb = &(getFirstLeaf(q)->sel->u.bb);    

  if (isSameBB(q->left, bb) && isSameBB(q->right, bb)) {
    return 1;
  } 
  return 0;
}
  //
  // in a general case, left query can have a different bb than right query
  // thus the index alighment is different. e.g. left ([0,0]-[N,N]), right ([N, N] - [2N, 2N])
  // index starts at 0 for left but something different likely for right query.
  // therefore get handle on each leaf and combine is not right, as handles would be covering  more than the actual data size.
  // and each leaf will vary
  //
// we identify the case when bb is the same through all the leaves,
// then we do composite query in fastbit to get best efficiency
int evaluateWithIdxOnBoundingBoxWithBitArray(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{  
  
  if (isBasedOnSameBB(q) > 0) {
    if (getFirstLeaf(q)->sel != NULL) {
      const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb = &(getFirstLeaf(q)->sel->u.bb);    
      mEvaluateBBRangeFancyQueryOnWhole(idxFile, q, timeStep, bb->start, bb->count);
    } else {
      mEvaluateBBRangeFancyQueryOnWhole(idxFile, q, timeStep, NULL, NULL);
    }
    return 0;
  }


  //ADIOS_SELECTION* sel = q->sel;
  ADIOS_VARINFO* v = q->varinfo;

  create_fastbit_internal(q);
  casestudyLogger_setPrefix(" created fastbit internal ");

  if (v == NULL) {
    ADIOS_QUERY* left = (ADIOS_QUERY*)(q->left);
    ADIOS_QUERY* right = (ADIOS_QUERY*)(q->right);

    if (evaluateWithIdxOnBoundingBoxWithBitArray(idxFile, left, timeStep) < 0) {
      return -1;
    }
    if (evaluateWithIdxOnBoundingBoxWithBitArray(idxFile, right, timeStep) < 0) {
      return -1;
    }

    casestudyLogger_setPrefix(" merge bitarray ");
    free(q->dataSlice);
    q->dataSlice = bitarray_create(q->rawDataSize);
    if (q->combineOp == ADIOS_QUERY_OP_OR) {            
      bitarray_or((uint32_t*)(left->dataSlice), (uint32_t*)(right->dataSlice), (uint32_t*)(q->dataSlice), BITNSLOTS(q->rawDataSize));      
    } else {
      bitarray_and((uint32_t*)(left->dataSlice), (uint32_t*)(right->dataSlice), (uint32_t*)(q->dataSlice), BITNSLOTS(q->rawDataSize));      
    }
    free(left->dataSlice);
    free(right->dataSlice);
    left->dataSlice = 0;
    right->dataSlice = 0;
    return 0;
  } else {
    // is a leaf
    if (q->rawDataSize == 0) {
      return 0;
    }
    //return mEvaluateWithIdxOnBBoxWithBitArrayOnVar(idxFile, q, timeStep);
    //!!return mEvaluateTestFullRange(idxFile, q, timeStep);
    if (getFirstLeaf(q)->sel != NULL) {
      const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb = &(getFirstLeaf(q)->sel->u.bb);    
      mEvaluateBBRangeFancyQueryOnWhole(idxFile, q, timeStep, bb->start, bb->count);
    } else {
      mEvaluateBBRangeFancyQueryOnWhole(idxFile, q, timeStep, NULL, NULL);
    }
    return 0;

  }

}


int evaluateWithIdxOnBlockWithBitArray(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{  
  
  ADIOS_VARINFO* v = q->varinfo;

  create_fastbit_internal(q);
  casestudyLogger_setPrefix(" created fastbit internal ");

  if (v == NULL) {
    ADIOS_QUERY* left = (ADIOS_QUERY*)(q->left);
    ADIOS_QUERY* right = (ADIOS_QUERY*)(q->right);

    if (evaluateWithIdxOnBlockWithBitArray(idxFile, left, timeStep) < 0) {
      return -1;
    }
    if (evaluateWithIdxOnBlockWithBitArray(idxFile, right, timeStep) < 0) {
      return -1;
    }

    casestudyLogger_setPrefix(" merge bitarray ");
    free(q->dataSlice);
    q->dataSlice = bitarray_create(q->rawDataSize);
    if (q->combineOp == ADIOS_QUERY_OP_OR) {            
      bitarray_or((uint32_t*)(left->dataSlice), (uint32_t*)(right->dataSlice), (uint32_t*)(q->dataSlice), BITNSLOTS(q->rawDataSize));      
    } else {
      bitarray_and((uint32_t*)(left->dataSlice), (uint32_t*)(right->dataSlice), (uint32_t*)(q->dataSlice), BITNSLOTS(q->rawDataSize));      
    }
    free(left->dataSlice);
    free(right->dataSlice);
    left->dataSlice = 0;
    right->dataSlice = 0;
    return 0;
  } else {
    // is a leaf
    if (q->rawDataSize == 0) {
      return 0;
    }

    if (q->sel != NULL) {
       const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb = &(q->sel->u.block);
       
       int absBlockCounter = wb->index;
       if (!q->file->is_streaming) 
	 absBlockCounter = query_utils_getGlobalWriteBlockId(wb->index, timeStep, v);
       
       if (v->blockinfo == NULL) {
	 common_read_inq_var_blockinfo(q->file, v);
       }
       
       ADIOS_VARBLOCK* blockSel = &(v->blockinfo[absBlockCounter]);
       
       mEvaluateBBRangeFancyQueryOnWhole(idxFile, q, timeStep, blockSel->start, blockSel->count);           
       return 0;
    }
    return -1;
  }

}

int evaluateWithIdxOnBoundingBox(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{  
  
  ADIOS_SELECTION* sel = q->sel;
  ADIOS_VARINFO* v = q->varinfo;

  create_fastbit_internal(q);
  casestudyLogger_setPrefix(" created fastbit internal ");

  if (v == NULL) {
    ADIOS_QUERY* left = (ADIOS_QUERY*)(q->left);
    ADIOS_QUERY* right = (ADIOS_QUERY*)(q->right);

    if (evaluateWithIdxOnBoundingBox(idxFile, left, timeStep) < 0) {
      return -1;
    }
    if (evaluateWithIdxOnBoundingBox(idxFile, right, timeStep) < 0) {
      return -1;
    }
#ifdef BITARRAY     
    casestudyLogger_setPrefix(" merge bitarray ");
    free(q->dataSlice);
    q->dataSlice = bitarray_create(q->rawDataSize);
    if (q->combineOp == ADIOS_QUERY_OP_OR) {            
      bitarray_or((uint32_t*)(left->dataSlice), (uint32_t*)(right->dataSlice), (uint32_t*)(q->dataSlice), BITNSLOTS(q->rawDataSize));      
    } else {
      bitarray_and((uint32_t*)(left->dataSlice), (uint32_t*)(right->dataSlice), (uint32_t*)(q->dataSlice), BITNSLOTS(q->rawDataSize));      
    }
    free(left->dataSlice);
    free(right->dataSlice);
    left->dataSlice = 0;
    right->dataSlice = 0;
#else
    setCombinedQueryInternal(q);
#endif
    return 0;
  } else {
    // is a leaf
    if (q->rawDataSize == 0) {
      return 0;
    }

    int blockStart=0; // relative
    int blockEnd = 0; // relative
    int i=0;

    const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = NULL;
    if (sel == NULL) {
        blockStart=0;
	if (v->nsteps == 1) { // could be streaming open
	  blockEnd = v->nblocks[0]-1;
	} else {
	  blockEnd = v->nblocks[timeStep] -1;      
	}
	log_debug(" got blockStart = %d, blockEnd = %d\n", blockStart, blockEnd);
    } else {
        bb = &(sel->u.bb);      
	if (v->blockinfo == NULL) {
	  common_read_inq_var_blockinfo(q->file, v);
	}
	blockStart = fastbit_adios_util_getRelativeBlockNumForPoint(v,bb->start,timeStep);
	uint64_t end[v->ndim];
	for (i=0; i<v->ndim; i++) {
	  end[i] = bb->start[i]+bb->count[i]-1;
	}
	blockEnd = fastbit_adios_util_getRelativeBlockNumForPoint(v, end, timeStep);
	log_debug(" figured blockStart = %d, blockEnd = %d\n", blockStart, blockEnd);
    }

    if ((blockStart < 0) || (blockEnd < 0) || (blockStart > blockEnd)) {
        adios_error (err_invalid_query_value, "Query processing failed. Unable to continue using index. vid=%d, dim[0]=%lld\n", v->varid, v->dims[0]);
	return -1;
    }

  casestudyLogger_setPrefix(" computed block ids to scan");
#ifdef FANCY_QUERY
#else  
    #ifdef BITARRAY
    uint32_t* bitSlice = bitarray_create(q->rawDataSize);
    #else
    uint16_t* bitSlice = malloc((q->rawDataSize)* sizeof(uint16_t));

    for (i=0; i<q->rawDataSize; i++) {
      bitSlice[i] = 0;
    }
    #endif
#endif

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);      
    char bitsArrayName[60];
    sprintf(bitsArrayName, "%ld_%d_%d_%d", fastbit_adios_getCurrentTimeMillis(), v->varid, timeStep, rank);

    uint64_t currBlockIdx = blockStart;

    uint64_t sumBlocksBeforeThisTimeStep = 0;
    if (v->nsteps > 1) { // non streaming case
      for (i=0;i<timeStep; i++) {
	sumBlocksBeforeThisTimeStep += v->nblocks[i];
      }
    }

    casestudyLogger_setPrefix(" start scanning");
    char casestudyLoggerPrefix[30];

    for (currBlockIdx=blockStart; currBlockIdx <= blockEnd; currBlockIdx++) {
      getHandle(timeStep, currBlockIdx, idxFile, q, 0);	      
      if (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle == 0) {
	log_warn(" Unable to construct fastbit query with NULL. Use _no_o idx method \n");
	return -1;
      }

      sprintf(casestudyLoggerPrefix, "block:%llu", currBlockIdx);
      casestudyLogger_setPrefix(casestudyLoggerPrefix);

      struct timespec evalStartT;
      casestudyLogger_getRealtime(&evalStartT);

      uint64_t count = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 
      
#ifdef FANCY_QUERY
      //printf("\nNot that this only works if the region covered exactly by blocks.\n");
      if (currBlockIdx == blockStart) {
	fastbit_iapi_register_selection_as_bit_array(bitsArrayName, ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
      } else {
	fastbit_iapi_extend_bit_array_with_selection(bitsArrayName, ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
      }
#else
      //i = currBlockIdx-blockStart;
      //uint64_t  coordinateArray[count];
      uint64_t* coordinateArray = malloc(count*sizeof(uint64_t));
      fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, count, 0);      

      casestudyLogger_pro_writeout(&evalStartT, "fastbitevaluated");

      struct timespec frameStartT;
      casestudyLogger_getRealtime(&frameStartT);
      //fastbit_selection_free(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
      ////fastbit_iapi_free_array_by_addr(q->_dataSlice); // if attached index    
      //fastbit_iapi_free_array_by_addr(q->dataSlice); // if attached index       
      clear_fastbit_internal_recursive(q);

      int absBlockIdx = currBlockIdx+sumBlocksBeforeThisTimeStep;
      int k=0;

      int isDefaultSelection = isSameRegion(v, bb, timeStep);

      for (k=0; k<count; k++) {
	uint64_t currPosInBlock = coordinateArray[k];
	int64_t currPos = currPosInBlock;
	
	if (!isDefaultSelection) {
	    currPos = getRelativeIdx(currPosInBlock, v, bb, absBlockIdx, timeStep);
	}

	log_debug("%lld th in block[%d],   =>  in actual box %lld  \n", currPosInBlock, absBlockIdx, currPos);
	if (currPos >= 0) {
            #ifdef BITARRAY
	    bitarray_setbit(bitSlice, currPos);
            #else
	    bitSlice[currPos] = 1;
            #endif
	}
      }
      casestudyLogger_frame_writeout(&frameStartT, "block processed");
#endif
      //casestudyLogger_setPrefix(" block processed");
      log_debug("----\n");
      free(coordinateArray);
    }

    casestudyLogger_setPrefix(" blocksProcessedIndividually ");

    //return fastbit_selection_create(dataType, dataOfInterest, dataSize, compareOp, &vv);

#ifdef FANCY_QUERY
    //
    // wrong result
    //fastbit_iapi_register_array(bitsArrayName, FastBitDataTypeBitRaw, q->dataSlice, q->rawDataSize);
    //FastBitSelectionHandle h = fastbit_selection_osr(bitsArrayName, FastBitCompareGreater, 0.1);
    //
    FastBitSelectionHandle h  = fastbit_selection_osr(bitsArrayName, FastBitCompareGreater, 0.0);
    fastbit_adios_util_checkNotNull(h, bitsArrayName);    
    ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = h;

#else
    free(q->dataSlice);
    q->dataSlice = bitSlice;

    fastbit_iapi_free_array_by_addr(q->dataSlice);

  #ifdef BITARRAY
    ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = 0;
  #else
    fastbit_iapi_register_array(bitsArrayName, FastBitDataTypeUShort, q->dataSlice, q->rawDataSize);
    FastBitSelectionHandle h = fastbit_selection_osr(bitsArrayName, FastBitCompareGreater, 0);
    fastbit_adios_util_checkNotNull(h, bitsArrayName);    
    ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = h;
  #endif

#endif    
 
    casestudyLogger_setPrefix(" summarized evaluation for bb");
    //fastbit_adios_util_checkNotNull(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, bitsArrayName);
  }
}

FastBitSelectionHandle createHandle(ADIOS_QUERY* q, const char* registeredArrayName)
{
  FastBitDataType  dataType = fastbit_adios_util_getFastbitDataType(q->varinfo->type);
  FastBitCompareType compareOp = fastbit_adios_util_getFastbitCompareType(q->predicateOp);

  uint64_t dataSize = q->rawDataSize;

  char* endptr;
  double vv = strtod(q->predicateValue, &endptr);

  if (dataType == FastBitDataTypeDouble) {
    //double vv = strtod(q->predicateValue, &endptr);
    assertValue(q->predicateValue, endptr);
    return fastbit_selection_osr(registeredArrayName, compareOp, vv);
  } else if ((dataType == FastBitDataTypeInt) || (dataType == FastBitDataTypeLong) || (dataType == FastBitDataTypeUInt) || (dataType == FastBitDataTypeULong)) {
    //long vv = strtol(q->predicateValue, &endptr, 10);    
    assertValue(q->predicateValue, endptr);
    return fastbit_selection_osr(registeredArrayName, compareOp, (double)vv);
  } else if (dataType == FastBitDataTypeFloat) {
    //float vv = strtof(q->predicateValue, &endptr);
    assertValue(q->predicateValue, endptr);
    return fastbit_selection_osr(registeredArrayName, compareOp, (double)vv);
  } else {
    return  NULL;
  }
}


void getHandleFromBlockAtLeafQuery(int timeStep, int blockIdx, ADIOS_FILE* idxFile, ADIOS_QUERY* q, uint64_t blockSize) 
{
  //double *keys = NULL; int64_t*offsets= NULL; uint32_t *bms = NULL;  
    uint64_t nk, no, nb;
    
    ADIOS_VARINFO* v = q->varinfo;
    
    ADIOS_FILE* dataFile = q->file;

    /*
    // read data from dataFile
    ADIOS_SELECTION* box = common_read_selection_writeblock(blockIdx);
    common_read_inq_var_blockinfo(dataFile, v);
    */
    if (v->blockinfo == NULL) {
      common_read_inq_var_blockinfo(dataFile, v);
    }

    //uint64_t blockSize = 0;
    if (blockSize == 0) {
      if (q->file->is_streaming == 1) {
        blockSize = fastbit_adios_util_getBlockSize(v, -1, blockIdx);
      } else {
        blockSize = fastbit_adios_util_getBlockSize(v, timeStep, blockIdx);
      }
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);      

    char blockDataName[40+strlen(q->condition)];
    //sprintf(blockDataName, "%d-%s-%d-%d-%ld-%d", v->varid, q->condition, timeStep, blockIdx, fastbit_adios_getCurrentTimeMillis(), rank);
    sprintf(blockDataName, "_%d_%d_%d_%ld_%d", v->varid,timeStep, blockIdx, fastbit_adios_getCurrentTimeMillis(), rank);

    FASTBIT_INTERNAL* itn = (FASTBIT_INTERNAL*)(q->queryInternal);

    struct timespec idxStartT;
    casestudyLogger_getRealtime(&idxStartT);

#ifdef __READ_BMS_AS_NEEDED__
    log_debug("FastBit will get bms as needed, not at the same time as key/offsets\n");
    if (fastbit_adios_util_readNoBMSFromIndexFile(idxFile, v, timeStep, blockIdx, &(itn->_keys), &nk, &(itn->_offsets), &no, &(itn->_bmsVarName)) < 0)     
#else
    if (fastbit_adios_util_readFromIndexFile(idxFile, v, timeStep, blockIdx, &(itn->_keys), &nk, &(itn->_offsets), &no, &(itn->_bms), &nb) < 0) 
#endif
    {
      // no idx for this variable, read from file:
      free(q->dataSlice);
      q->dataSlice = malloc(common_read_type_size(v->type, v->value)*blockSize);
      
      ADIOS_SELECTION* box = common_read_selection_writeblock(blockIdx);   
      common_read_inq_var_blockinfo(dataFile, v);        

      int errorCode = 0;
      if (dataFile->is_streaming == 1) {
	errorCode = common_read_schedule_read_byid(dataFile, box, v->varid, 0, 1, NULL, q->dataSlice);
      } else {
	errorCode = common_read_schedule_read_byid(dataFile, box, v->varid, timeStep, 1, NULL, q->dataSlice);
      }

      if (errorCode != 0) {
          log_error("      %s:%d  schedule read error code = %d adios_error=%d \n", __func__, __LINE__, errorCode, adios_errno);
          //return errorCode;
	  return;
      }
      common_read_perform_reads(dataFile,1);

      FastBitDataType  dataType = fastbit_adios_util_getFastbitDataType(q->varinfo->type);
      FastBitCompareType compareOp = fastbit_adios_util_getFastbitCompareType(q->predicateOp);

      setQueryInternal(q, compareOp, dataType, blockSize, blockDataName);
      common_read_selection_delete(box);
      return;
    }

    //casestudyLogger_idx_writeout(&idxStartT, "index file visited ");
    casestudyLogger_idx_writeout(&idxStartT, "index file visited ");

    //int err = fastbit_iapi_register_array(blockDataName, fastbit_adios_util_getFastbitDataType(v->type), q->_dataSlice, blockSize);
    uint64_t nv = blockSize;
    itn->_arrayName = malloc(strlen(blockDataName)+2);
    sprintf(itn->_arrayName, "%s", blockDataName);
    
    int ierr = fastbit_iapi_register_array_index_only(itn->_arrayName, fastbit_adios_util_getFastbitDataType(v->type), &nv, 1 , 
#ifdef __READ_BMS_AS_NEEDED__
						      itn->_keys, nk, itn->_offsets, no, itn, adios_bmreader);
#else
						      itn->_keys, nk, itn->_offsets, no, itn->_bms, adios_bmreader);
#endif
      /*
    if (ierr != 0) {
      log_error(" registering array failed. fastbit err code = %lld\n", ierr);
    }
    //printData(bms, adios_unsigned_integer, nb);

    ierr = fastbit_iapi_attach_index (blockDataName, keys, nk, offsets, no, bms, mybmreader);
      */
    if (ierr < 0) {
      log_debug(" reattaching index failed. %s, fastbit err code = %d\n", blockDataName, ierr);
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = NULL;
      //result = ierr;
    } else {
      FastBitSelectionHandle h = createHandle(q, blockDataName); //fastbit_selection_osr(blockDataName, getFastbitCompareType(q->_op), q->_value);
      fastbit_adios_util_checkNotNull(h, blockDataName);
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = h;      
      //fastbit_adios_util_checkNotNull(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, blockDataName);
    }
}



void printQueryData(ADIOS_QUERY* q, FastBitDataType dataType, int timeStep) {
  uint64_t dataSize = q->rawDataSize;
  int j;
  int batchSize = 31;
  log_debug ("::\t %s At timestep: %d datasize=%llu \n\t\t   raw data:  [", q->condition, timeStep, dataSize);
  for (j = 0; j < dataSize; j++) {
    if ((j < batchSize) || ((dataSize -j) < batchSize)) {
      if ((j % 10) == 0) {
	log_debug(" \n\t\t           ");
      }
      if (dataType == FastBitDataTypeDouble) {
	log_debug(" %lg ", ((double *)(q->dataSlice))[j]);
      } else if (dataType == FastBitDataTypeFloat) {
	log_debug(" %g ", ((float *)(q->dataSlice))[j]);
      } else if (dataType == FastBitDataTypeUInt) {
	log_debug("%d ", ((uint32_t  *)(q->dataSlice))[j]);
      } else if (dataType == FastBitDataTypeULong) {
	log_debug("%lld ", ((uint64_t  *)(q->dataSlice))[j]);
      } else {
	//log_debug("\t%g ", ((uint32_t *)(q->_dataSlice))[j]);
	log_debug(" *  ");
      }
    } else {
      if (j == batchSize) {
	log_debug(" ... ");
      }
      continue;
      //break;
    }
  }
  log_debug ("]\n");
}


// q is a leaf
int readWithTimeStepNoIdx(ADIOS_QUERY* q, int timeStep) {
  if (q->rawDataSize == 0) {
    return 0;
  }
  FastBitDataType  dataType = fastbit_adios_util_getFastbitDataType(q->varinfo->type);
  FastBitCompareType compareOp = fastbit_adios_util_getFastbitCompareType(q->predicateOp);

  int ts = timeStep;
  if (q->file->is_streaming) 
      ts = 0; // if file is opened as stream, the actual 'actual timestep' is always 0

  if (q->dataSlice == NULL) {
    q->dataSlice = malloc(common_read_type_size(q->varinfo->type, q->varinfo->value)*q->rawDataSize);
  }

  int errorCode = common_read_schedule_read_byid (q->file, q->sel, q->varinfo->varid, 
                                                  ts, 1, NULL,  q->dataSlice);
  log_debug("      schedule read error code = %d adios_error=%d \n", errorCode, adios_errno);
  if (errorCode != 0) {
    return errorCode;
  }

  common_read_perform_reads (q->file, 1); // return 0 regardless whether data is valid, so donnt need to check return value     
  log_debug("      perform read got error code = %d adios_errno=%d\n", errorCode, adios_errno);
  if (adios_errno != 0) {
    return -1;
  }
  
  printQueryData(q, dataType, timeStep);
  uint64_t dataSize = q->rawDataSize;
  //common_free_varinfo(q->_var);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);      

  char datasetName[strlen(q->condition) + 40];
  //sprintf(datasetName, "noidx_%s_%d_%ld_%d", q->condition, timeStep, fastbit_adios_getCurrentTimeMillis(),rank);  
  sprintf(datasetName, "_noidx_%d_%ld_%d", timeStep, fastbit_adios_getCurrentTimeMillis(),rank);  
  setQueryInternal(q, compareOp, dataType, dataSize, datasetName);

  return 0;
}

/*
int readWithTimeStep(ADIOS_QUERY* q, int timeStep, ADIOS_FILE* idxFile) {
  if (idxFile == NULL) {
    return readWithTimeStepNoIdx(q, timeStep);
  } else {
    return readWithTimeStepWithIdx(q, timeStep, idxFile);
  }
}
*/

int prepareData(ADIOS_QUERY* q, int timeStep) 
{
  if (q->onTimeStep == timeStep) {
    log_debug("::\t query data has been read for timestep: %d\n", timeStep);
    return 0;
  } 

  if (q->varinfo != NULL) {
    if (q->queryInternal != NULL) {
      //fastbit_selection_free(q->queryInternal);
      clear_fastbit_internal(q);
    }
    
    if (q->queryInternal == NULL) {
      q->queryInternal = malloc(sizeof(FASTBIT_INTERNAL));
    }

    int errorCode = readWithTimeStepNoIdx(q, timeStep);
    if (errorCode != 0) {
      return errorCode;
    }
  } else {
    int errorCode1 = prepareData((ADIOS_QUERY*)(q->left), timeStep);
    if (errorCode1 != 0) {
      return errorCode1;
    }

    int errorCode2 = prepareData((ADIOS_QUERY*)(q->right), timeStep);
    if (errorCode2 != 0) {
      return errorCode2;
    }

    setCombinedQueryInternal(q);
  }
  q->onTimeStep = timeStep;
  q->maxResultsDesired = 0;
  q->resultsReadSoFar = 0;

  return 0;
}



//
// selections in q are all block(same or different)
//
/*
FastBitSelectionHandle blockSelectionFastbitHandle(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{
  if (q->_var != NULL) {
    const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb = &(q->_sel->u.block);
    //int absBlockCounter = query_utils_getGlobalWriteBlockId(wb->index, timeStep, q->_var);    
    getHandleFromBlockAtLeafQuery(timeStep, wb->index, idxFile, q);
  } else {
    FastBitSelectionHandle left = blockSelectionFastbitHandle(idxFile, q->_left, timeStep);
    FastBitSelectionHandle right = blockSelectionFastbitHandle(idxFile, q->_right, timeStep);

    if (q->_leftToRightOp == ADIOS_QUERY_OP_AND) {
      return fastbit_selection_combine(left, FastBitCombineAnd, right);
    } else {
      return fastbit_selection_combine(left, FastBitCombineOr, right);
    }        
  }
}
*/
void blockSelectionFastbitHandle(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{
  if (q->queryInternal == NULL) {
    q->queryInternal = malloc(sizeof(FASTBIT_INTERNAL));
  }
  
  if (q->varinfo != NULL) {
    const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb = &(q->sel->u.block);
    //int absBlockCounter = query_utils_getGlobalWriteBlockId(wb->index, timeStep, q->_var);    
    getHandleFromBlockAtLeafQuery(timeStep, wb->index, idxFile, q, 0);
  } else {
    blockSelectionFastbitHandle(idxFile, q->left, timeStep);
    blockSelectionFastbitHandle(idxFile, q->right, timeStep);

    setCombinedQueryInternal(q);
  }
}


void getHandle(int timeStep, int blockIdx, ADIOS_FILE* idxFile, ADIOS_QUERY* q, uint64_t blockSize) 
{
  ADIOS_FILE* dataFile = q->file;
  ADIOS_VARINFO* v = q->varinfo;

  FastBitSelectionHandle result = NULL;
  if (v == NULL) {
    getHandle(timeStep, blockIdx, idxFile, q->left, blockSize);
    getHandle(timeStep, blockIdx, idxFile, q->right, blockSize);

    setCombinedQueryInternal(q);
    /*
    ADIOS_QUERY* left  = (ADIOS_QUERY*)(q->_left);
    ADIOS_QUERY* right = (ADIOS_QUERY*)(q->_right);

    if (q->_leftToRightOp == ADIOS_QUERY_OP_AND) {
      q->_queryInternal = fastbit_selection_combine(left->_queryInternal, FastBitCombineAnd, right->_queryInternal);
    } else {
      q->_queryInternal = fastbit_selection_combine(left->_queryInternal, FastBitCombineOr, right->_queryInternal);
      } */   
  } else {
    getHandleFromBlockAtLeafQuery(timeStep, blockIdx, idxFile, q, blockSize);
  }     
}

int64_t  applyIndexIfExists (ADIOS_QUERY* q, int timeStep) 
{
  if (q->onTimeStep == timeStep) {
    //log_debug("::\t query index data has been read for timestep: %d\n", gCurrentTimeStep);
    return 0;
  } 

  // call fastbit_estimate_num_hits(selection)
  ADIOS_QUERY* leaf = getFirstLeaf(q);
  ADIOS_FILE* f = leaf->file;
  const char* basefileName = f->path;

  int64_t result = -1;

  ADIOS_FILE* idxFile = ((FASTBIT_INTERNAL*)(q->queryInternal))->_idxFile;

  //MPI_Comm comm_dummy = MPI_COMM_SELF;
  //ADIOS_FILE* idxFile = fastbit_adios_util_getFastbitIndexFileToRead(basefileName, comm_dummy);
    
  if (idxFile != NULL) {
    //casestudyLogger_starts("idxEval");
    //clear_fastbit_internal(q);
      if ((leaf->sel == NULL) || (leaf->sel->type == ADIOS_SELECTION_BOUNDINGBOX)) {
#ifdef BITARRAY     
	  if (evaluateWithIdxOnBoundingBoxWithBitArray(idxFile,  q, timeStep) >= 0) {
#else
	  if (evaluateWithIdxOnBoundingBox(idxFile,  q, timeStep) >= 0) {
#endif
	     casestudyLogger_bms_print();
	     casestudyLogger_idx_print();
	     casestudyLogger_pro_print();
	     casestudyLogger_frame_print();
	     casestudyLogger_setPrefix(" preparedBoundingBox ");
	     if (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle != 0) {
	        result = fastbit_selection_estimate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);	
		casestudyLogger_setPrefix(" estimateDone ");
	     } else {	        
	        result = bitarray_countHits(q->dataSlice, BITNSLOTS(q->rawDataSize));
		casestudyLogger_setPrefix(" estimateDoneBitArray ");
	     }
	  }
      } else if (leaf->sel->type == ADIOS_SELECTION_WRITEBLOCK) {
#ifdef BITARRAY
	    evaluateWithIdxOnBlockWithBitArray(idxFile, q, timeStep);
#else
	    blockSelectionFastbitHandle(idxFile, q, timeStep);
#endif
	    //result= fastbit_selection_estimate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);       	  
	    if (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle != 0) {
	      result = fastbit_selection_estimate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);	
	      casestudyLogger_setPrefix(" estimateDone ");
	    } else {	        
	      result = bitarray_countHits(q->dataSlice, BITNSLOTS(q->rawDataSize));
	      casestudyLogger_setPrefix(" estimateDoneBitArray ");
	    }
      } 

      log_debug("idx evaluated with result=%" PRId64 "\n", result);
      if (result > -1) {
	q->onTimeStep = timeStep;
	q->maxResultsDesired = 0;
	q->resultsReadSoFar = 0;

	//return result;
      } // otherwise, use no idx method
      //common_read_close(idxFile);      
      //casestudyLogger_ends("idxEval");
  }
  
  
  //free(idxFile); //causes crash
  return result;
}


int adios_query_fastbit_can_evaluate(ADIOS_QUERY* q) 
{
    /* Return 1 if fastbit index file exists.
       Even though Fastbit library works without index
       we return 0 here, so that other query methods may pick up
       this query.
       If no method is found, the common layer may still call 
       the fastbit method (as default) to evaluate the query.
     */
   ADIOS_QUERY* leaf = getFirstLeaf(q);
   return fastbit_adios_util_FastbitIndexFileExists (leaf->file->path);
 }

/*
int assertTimeStepValidWithQuery(ADIOS_QUERY* q)
{
  ADIOS_QUERY* leaf = getFirstLeaf(q);
  if (leaf->varinfo == NULL) {
    log_debug("No variable on query leaf. Can not continue. Exiting.\n");
    return -1;
  }

  int timestep = gCurrentTimeStep;

  if (leaf->file->is_streaming) {
    int currentFileStep = leaf->file->current_step;
    if (timestep != currentFileStep) {
      adios_query_set_timestep(currentFileStep);
    }
    return 0;
  }

  if (leaf->varinfo->nsteps <= timestep) {
    log_debug("timestep %d is more than variables limit: %d, can not evaluate.\n", timestep, leaf->varinfo->nsteps);
    return -1;
  }
  return 0;
}
*/
int64_t adios_query_fastbit_estimate(ADIOS_QUERY* q, int incomingTimestep) 
{
  if (q == NULL) {
    return -1;
  }

  //int timeStep = gCurrentTimeStep;

  adios_query_fastbit_init();
  /*  if (assertTimeStepValidWithQuery(q) != 0) {
    return -1;
    }*/

  create_fastbit_internal(q);

  int timeStep = adios_get_actual_timestep(q, incomingTimestep);
  int64_t estimate = applyIndexIfExists(q, timeStep);
  if (estimate > 0) {
    return estimate;
  } else if (estimate == 0) { // estimated was called before
    return fastbit_selection_estimate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
  }

  //
  // no idx estimation
  //
  int errorCode = prepareData(q, timeStep);
  if (errorCode != 0) {
    return -1;
  }
  return fastbit_selection_estimate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
}
 

int64_t call_fastbit_evaluate(ADIOS_QUERY* q, int timeStep, uint64_t _maxResult) 
{
  create_fastbit_internal(q);
  
  log_debug("raw data size = %" PRId64 ", %s\n", q->rawDataSize, q->condition);
    
  int64_t estimate = applyIndexIfExists(q, timeStep);
  if (estimate < 0) {   // use no idx
    log_debug(" use no idx on timestep %d\n", timeStep);
    int errorCode = prepareData(q, timeStep);
    if (errorCode != 0) {
      return -1;
    }
  }

  if (q->maxResultsDesired > 0) { // evaluated already
      if (_maxResult <= 0) { // stay put
	  return q->maxResultsDesired;
      } 

      if (q->maxResultsDesired > _maxResult) {
         q->maxResultsDesired = _maxResult;
	 return q->maxResultsDesired;
     }
     log_debug(":: user required more results. will evaluate again. \n");
  }

  //casestudyLogger_starts("evaluatingFastbit");

#ifdef BITARRAY
  if (q->queryInternal == 0) {
#else
  if ((q->queryInternal == 0) || (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle == 0)) {
#endif
    //log_error(">>  Unable to use fastbit to evaluate NULL query.\n"); 
    log_debug("query is NULL, result is NULL.");
    return 0;
  }

  int64_t numHits = -1;
  if (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle != 0) {
    numHits = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 
  } else {
    numHits = bitarray_countHits(q->dataSlice, BITNSLOTS(q->rawDataSize));
  }

  log_debug(":: ==> fastbit_evaluate() num of hits found for [%s] = %lld, at timestep %d \n", q->condition, numHits, timeStep);  

  //casestudyLogger_ends("evaluatingFastbit");

  if (numHits < 0) {
    return 0;
  }

  if (numHits <= _maxResult) {
    // take it as max
    q->maxResultsDesired = numHits; 
  } else if (_maxResult > 0)  {
    q->maxResultsDesired = _maxResult;
  } else {
    q->maxResultsDesired = numHits; 
  }

  return q->maxResultsDesired;
}


void printOneSpatialCoordinate(int dim, uint64_t* spatialCoordinates)
{
      int k;
      log_debug(" spatial = [");
      for (k=0; k<dim; k++) {
	log_debug("%lld ", spatialCoordinates[k]);
      }
      log_debug("]\n");

}

 void fillUp(int dimSize, uint64_t* spatialCoordinates, uint64_t i, uint64_t* pointArray, int fortran_order) 
{
  int k=0;

  //  int fortran_order = futils_is_called_from_fortran();

  for (k = 0; k < dimSize; k++) {	  
    uint64_t idx = i * dimSize + k;
    //log_debug(" points[%d] = %lld ", idx, spatialCoordinates[k]);

    if (fortran_order == 1) {
      pointArray[idx] = spatialCoordinates[dimSize-1-k];
      //printf(" points[%lld] = %lld dimSize=%d k=%d sp=%lld", idx, pointArray[idx], dimSize, k, spatialCoordinates[k]); 
    } else {
      pointArray[idx] = spatialCoordinates[k];
    }
    

    //printf("   points[%d] = %lld ", idx, spatialCoordinates[k]);
    //pointArray[idx] = spatialCoordinates[k];
  }	
  log_debug("\n");
}

ADIOS_SELECTION* getSpatialCoordinatesDefault(ADIOS_VARINFO* var, uint64_t* coordinates, uint64_t retrivalSize, int timeStep)
{
  uint64_t arraySize = retrivalSize * (var->ndim);
  uint64_t* pointArray = (uint64_t*) (malloc(arraySize  * sizeof(uint64_t)));

  int isFortranClient = futils_is_called_from_fortran();

  int i;
  for (i=0; i<retrivalSize; i++) {
    uint64_t spatialCoordinates[var->ndim];
    getCoordinateFromVariable(coordinates[i], var, var->ndim, spatialCoordinates);
    
    fillUp(var->ndim, spatialCoordinates, i, pointArray, isFortranClient);
  }
  ADIOS_SELECTION* result =  common_read_selection_points(var->ndim, retrivalSize, pointArray);
  //free(pointArray); // user has to free this
  return result;
}

ADIOS_SELECTION* getSpatialCoordinates(ADIOS_SELECTION* outputBoundary, uint64_t* coordinates, uint64_t retrivalSize, ADIOS_VARINFO* v, int timeStep)
{
  int k = 0;
  uint64_t i=0;
  int isFortranClient = futils_is_called_from_fortran();

  switch (outputBoundary->type) {
  case  ADIOS_SELECTION_BOUNDINGBOX:    
    {
      const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(outputBoundary->u.bb);

      uint64_t arraySize = retrivalSize * (bb->ndim);
      uint64_t* pointArray = (uint64_t*) (malloc(arraySize  * sizeof(uint64_t)));
      
      for (i=0; i<retrivalSize; i++) {
	   uint64_t spatialCoordinates[bb->ndim];
	   //getCoordinateFromBox(coordinates[i], bb, bb->ndim, spatialCoordinates);
	   posToSpace(coordinates[i], isFortranClient, bb->count, spatialCoordinates, bb->ndim, bb->start);
	   //fillUp(bb->ndim, spatialCoordinates, i, pointArray, isFortranClient);

	   fillUp(bb->ndim, spatialCoordinates, i, pointArray, 0); // already fortran coordinates from posToSpace(.. isFortranClient ..)
      }
      ADIOS_SELECTION* result =  common_read_selection_points(bb->ndim, retrivalSize, pointArray);    
      //free(pointArray); // user has to free this
      return result;
      break;
    }
  case ADIOS_SELECTION_POINTS:
    {
      const ADIOS_SELECTION_POINTS_STRUCT *points = &(outputBoundary->u.points);	      

      uint64_t arraySize = retrivalSize * (points->ndim);
      uint64_t* pointArray = (uint64_t*) (malloc(arraySize  * sizeof(uint64_t)));
      
      for (i=0; i<retrivalSize; i++) {	
	uint64_t spatialCoordinates[points->ndim];
	getCoordinateFromPoints(coordinates[i], points, spatialCoordinates);
	fillUp(points->ndim, spatialCoordinates, i, pointArray, isFortranClient);
      }
      ADIOS_SELECTION* result = common_read_selection_points(points->ndim, retrivalSize, pointArray);	      
      //free(pointArray); // user has to free this
      return result;
      //printOneSpatialCoordinate(points->ndim, spatialCoordinates);      
      
      break;
    }
    //  if it is blocks, should use bounding box to retrive the coordinates
  case ADIOS_SELECTION_WRITEBLOCK:
    {
      const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb = &(outputBoundary->u.block);
      
      uint64_t arraySize = retrivalSize * (v->ndim);
      uint64_t* pointArray = (uint64_t*) (malloc(arraySize  * sizeof(uint64_t)));
      
      for (i=0; i<retrivalSize; i++) {
	   uint64_t spatialCoordinates[v->ndim];
	   //create bb from block;
	   int absBlockCounter = query_utils_getGlobalWriteBlockId(wb->index, timeStep, v);
	   //getCoordinateFromBlock(coordinates[i], &(v->blockinfo[absBlockCounter]), v->ndim, spatialCoordinates, v->ndim);
	   ADIOS_VARBLOCK* blockSel = &(v->blockinfo[absBlockCounter]);
	   posToSpace(coordinates[i], isFortranClient, blockSel->count, spatialCoordinates, v->ndim, blockSel->start);

	   //fillUp(v->ndim, spatialCoordinates, i, pointArray, isFortranClient); 
	   fillUp(v->ndim, spatialCoordinates, i, pointArray, 0); // already fortran coordinates from posToSpace(.. isFortranClient ..)
      }
      ADIOS_SELECTION* result = common_read_selection_points(v->ndim, retrivalSize, pointArray);
      //free(pointArray); // user has to free this
      return result;
      break;      
    }
  default:
    log_error("Error: Type of selection is not supported!\n\n");
    return NULL;
  }
}

ADIOS_QUERY* getFirstLeaf(ADIOS_QUERY* q) {
  if (q == NULL) {
    return NULL;
  }

  if ((q->left == NULL) && (q->right == NULL)) {
    return q;
  }
  return getFirstLeaf(q->left);
}

int  adios_query_fastbit_evaluate(ADIOS_QUERY* q,
				  int incomingTimestep,
				  uint64_t batchSize, 
				  ADIOS_SELECTION* outputBoundary, 
				  ADIOS_SELECTION** result)
{
  casestudyLogger_init();
  casestudyLogger_starts("queryArrived. initfastbit");

  adios_query_fastbit_init();

  /*
  if (q->_onTimeStep < 0) {
    log_error(":: Error: need to call evaluate first! Exit.\n");
    return -1;
  }
  */

  /*if (assertTimeStepValidWithQuery(q) != 0) {
    return -1;
    }*/

  if (batchSize == 0) {
    log_debug(":: ==> will not fetch. batchsize=0\n");
    return -1;
  }

  int timeStep = adios_get_actual_timestep(q, incomingTimestep);

  call_fastbit_evaluate(q, timeStep, 0);
  log_debug("::\t max=%llu  lastRead=%llu batchsize=%llu\n", q->maxResultsDesired, q->resultsReadSoFar, batchSize);

  uint64_t retrivalSize = q->maxResultsDesired - q->resultsReadSoFar;
  if (retrivalSize == 0) {
     log_debug(":: ==> no more results to fetch\n");
     return 0;
  }

  if (retrivalSize > batchSize) {
      retrivalSize = batchSize;
  }

  //uint64_t coordinates[retrivalSize];
  uint64_t* coordinates = (uint64_t*) calloc(retrivalSize, sizeof(uint64_t));
  if (coordinates == 0) {
    log_error("Unable to allocate for coordiantes.");
    return -1;
  }

  struct timespec startT;
  casestudyLogger_getRealtime(&startT);

  if (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle != 0) {
    fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinates, retrivalSize, q->resultsReadSoFar);
    casestudyLogger_idx_writeout(&startT, "getCoordinates");
  } else {
    bitarray_getHits(q->dataSlice, coordinates, BITNSLOTS(q->rawDataSize), q->resultsReadSoFar, retrivalSize);
  }

  q->resultsReadSoFar += retrivalSize;
  
  if (outputBoundary == 0) {
    ADIOS_QUERY* firstLeaf = getFirstLeaf(q);
    if ((firstLeaf == NULL) || (firstLeaf->varinfo == NULL)) {
	log_error(":: Error: unable to get a valid first leaf! Exit. \n");
	free(coordinates);
	return -1;
      }
    if (firstLeaf->sel == NULL) {
      *result = getSpatialCoordinatesDefault(firstLeaf->varinfo, coordinates, retrivalSize, timeStep);
    } else {
      *result = getSpatialCoordinates(firstLeaf->sel, coordinates, retrivalSize, firstLeaf->varinfo, timeStep);
    }
    free(coordinates);
  } else {
    //*result = getSpatialCoordinates(outputBoundary, coordinates, retrivalSize);
    // variable needs to be in place to handle the block information
    // not sure wheather this is well defined case of combined query?! but the first varibale will be used for block information calculation
    *result = getSpatialCoordinates(outputBoundary, coordinates, retrivalSize, getFirstLeaf(q)->varinfo, timeStep);

    free(coordinates);
    if (*result == 0) {
      return -1;
    }
  }
  casestudyLogger_frame_writeout(&startT, "bitarrayGetHits");

  // print results
  /*
  int i=0; 
  log_debug("\n:: coordinates: [\n");
  for (i=0; i<retrivalSize; i++) {  
    if (i<100) {
      log_debug("%lld ", coordinates[i]);
    } else {
      break;
    }
  }
  log_debug("]\n\n");
  */
  casestudyLogger_frame_print();
  casestudyLogger_ends("evaluation_total");

  if (q->resultsReadSoFar == q->maxResultsDesired) {
    return 0;
  } else {
    return 1;
  }
}

//void adios_fastbit_free_query(ADIOS_QUERY* query) 
void  adios_query_fastbit_free(ADIOS_QUERY* query) 
{  
  if (query == NULL) {
    return;
  }

  /*
  log_debug(":: free %s  has parent? %d\n", query->condition, query->hasParent);
  
  free(query->predicateValue);
  free(query->condition);
  
  //adios_selection_delete(query->_sel);
  common_read_free_varinfo(query->varinfo);

  // can free _queryInternal only once
  //if(query->hasParent == 0) {
    //fastbit_selection_free(query->queryInternal);
  clear_fastbit_internal(query);
  free(query->queryInternal);
    //}
  free(query->dataSlice);
  query->dataSlice = 0; 

  free(query);
  query = 0;
  */

  FASTBIT_INTERNAL* s = (FASTBIT_INTERNAL*)(query->queryInternal);  
  if (query->hasParent == 0) {
    if ((s!= NULL) && (s->_idxFile != NULL)) {
      common_read_close(s->_idxFile);
    }
  }

  clear_fastbit_internal(query);
  free(query->queryInternal);

  //fastbit_iapi_free_all();

}


/*

*/






void getVarName(const char* sliceStr, char** varName, char** dimDef)
{                                           
  //char str[strlen(sliceStr)];
  //strncpy(str, sliceStr, strlen(sliceStr));
  char* str = strdup(sliceStr);
  //str[strlen(sliceStr)]=0;                                                                                                                                                      
  char* rest = strtok(str, "[");

  *varName = strdup(rest);

  rest = strtok(NULL, "]"); 

  if (rest == NULL) {
    //strcpy(dimDef, "");
    *dimDef = NULL;
  } else {
    *dimDef = strdup(rest);
  }

  free(str);
}



/*
*/

//
// dimDef is asusmed to be: start0 count0, start1 count1, .. , startN countN
//
void createBox(ADIOS_VARINFO* v, char* dimDef, uint64_t* start, uint64_t* count)
{
  // assigns default:
  int i;

  for (i=0 ; i<v->ndim; i++) {
    start[i] = 0;
    count[i] = v->dims[i];
  }

 if ((dimDef == NULL) || (strlen(dimDef) == 0)) {
    log_debug("::\t There is no restriction on dimension.\n");
    //use default 
    return;
  }
  
  const char* comma = ",";

  // int column=-1, columnStarts=-1, columnEnds=0;

  char *dimSpecStart, *end;

  dimSpecStart = end = dimDef;
  i = 0;

  while (*end) {
    switch (*end) {
    case ',':
      *end = '\0';      
      //*column = atoi(dimSpecStart);
      
      count[i]=atol(dimSpecStart);
      dimSpecStart = end+1;
      *end = ',';
      i++;
      break;
    case ':':
      *end = '\0';
      //*columnStarts=atoi(dimSpecStart);
      start[i] = atol(dimSpecStart);      
      *end = ':';
      dimSpecStart = end+1;
      //break;
    }
    end++;
  }

  count[i] = atol(dimSpecStart);

  for (i=0; i<v->ndim; i++) {
    if (count[i] <= 0) {
      count[i] = v->dims[i] - start[i];
    }
  }
  /*
  if (*column ==-1) {
    *column = atoi(dimSpecStart);
  } else if (*columnStarts >= 0) {    
    *columnEnds = atoi(dimSpecStart);
  }
  */

}
/*
void createBox(ADIOS_VARINFO* v, char* dimDef, uint64_t* start, uint64_t* count)
{
  int column=-1, columnStarts=-1, columnEnds=0;
  parseSlice(dimDef, &column, &columnStarts, &columnEnds);
  log_debug(" column= %d, %d:%d \n\n", column, columnStarts, columnEnds);

  int j;
  for (j = 0; j < v->ndim; j++) {
    if (column == -1) {
      start[j] = 0;
      count[j] = v->dims[j];
    } else if (j == column) {
      
      if (columnStarts > 0) {
	start[j] = columnStarts;
      } else {
	start[j] = 0;
      }
      if (columnEnds > 0) {
	count[j] = columnEnds - start[j];
      } else {
	count[j] = v->dims[j] - start[j];
      }
      
      //start[j] = column;
      //count[j] = 1;
    } else {// not the column
      start[j] = 0;
      count[j] = v->dims[j];
    }
  }
  //start[v->ndim]=0;
  //count[v->n dim]=0;
}
*/

enum ADIOS_PREDICATE_MODE getOp(const char* opStr) 
{
  if ((strcmp(opStr, ">=") == 0) || (strcmp(opStr, "GE") == 0)) {
    return ADIOS_GTEQ;
  } else if ((strcmp(opStr, "<=") == 0) || (strcmp(opStr, "LE") == 0)) {
    return ADIOS_LTEQ;
  } else if ((strcmp(opStr, "<") == 0) || (strcmp(opStr, "LT") == 0)) {
    return ADIOS_LT;
  } else if ((strcmp(opStr, ">") == 0) || (strcmp(opStr, "GT") == 0)) {
    return ADIOS_GT;
  } else if ((strcmp(opStr, "=") == 0) || (strcmp(opStr, "EQ") == 0)) {
    return ADIOS_EQ;
  } else { // if (strcmp(opStr, "!=") == 0) {
    return ADIOS_NE;
  }
}

/*
ADIOS_QUERY* getQuery(const char* condition, ADIOS_FILE* f) 
{
    
    char* varStr;// = malloc(sizeof(char) * strlen(condition)); 

    //char opStr[5];
    char* opStr ;// = malloc(sizeof(char)*5);
    //char value[strlen(condition)];
    char* valueStr; // = malloc(sizeof(char)*strlen(condition));
    
    //parseVar(condition, &varStr, &opStr, &valueStr); // these values are ok in parseVar but just valueStr became null after this call!!!??

    //char str[strlen(condition)];
    //strncpy(str, condition, strlen(condition));
    //str[strlen(condition)] = 0;
    char* str = strdup(condition);
   

    char * pch;
    pch = strtok (str," ");  


    varStr = strdup(pch);
    
    
    pch = strtok (NULL, " ");
    opStr = strdup(pch);
    
    
    pch = strtok (NULL, " ");
    valueStr = strdup(pch);
    
    char* varName; //[strlen(varStr)];
    char* dimDef; //[strlen(varStr)];
    getVarName(varStr, &varName, &dimDef);


    ADIOS_VARINFO* v = common_read_inq_var(f, varName);
    //ADIOS_VARINFO* v = getAdiosVariable(f, varName);
    
    if (v == NULL) {
      free(valueStr);free(opStr);free(varStr);free(varName);free(str);
      if (dimDef != NULL) 
	{free(dimDef);}
      return NULL;
    }

    uint64_t* start = malloc(sizeof(uint64_t)*v->ndim);
    uint64_t* count = malloc(sizeof(uint64_t)*v->ndim);

    createBox(v, dimDef, start, count);

    ADIOS_SELECTION* sel =adios_selection_boundingbox (v->ndim, start, count);
    
    common_read_free_varinfo(v);

    ADIOS_QUERY* q = adios_query_create(f, varName, sel, getOp(opStr), valueStr);


    free(valueStr);free(opStr);free(varStr);free(varName);free(str);
    //free(start); free(count); // if deleted, then adios_sel values would be affected
    if (dimDef != NULL) 
      {free(dimDef);}

    log_debug("::\t query created for: %s\n", condition);
    return q;
}
*/






