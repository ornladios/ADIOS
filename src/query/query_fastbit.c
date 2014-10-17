#include "public/adios_query.h"
   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/common_read.h"
#include "core/adios_logger.h"
#include "fastbit_adios.h"
#include <iapi.h>
#include <math.h>

typedef struct {
  double* _keys; 
  int64_t* _offsets; 
  uint32_t* _bms;
  char* _arrayName;
  //void* _rawData;
  FastBitSelectionHandle _handle;
} FASTBIT_INTERNAL;
 

void create_fastbit_internal(ADIOS_QUERY* q) 
{
  if (q->queryInternal == NULL) {
     FASTBIT_INTERNAL* internal = malloc(sizeof(FASTBIT_INTERNAL));
     internal->_keys = NULL;
     internal->_offsets = NULL;
     internal->_bms = NULL;
     internal->_arrayName = NULL;
     internal->_handle = NULL;

     q->queryInternal = internal;
  }
  
  if (q->left != NULL) {
    create_fastbit_internal(q->left);
  } 
  if (q->right != NULL) {
    create_fastbit_internal(q->right);
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

  s->_keys = NULL; s->_bms = NULL; s->_offsets = NULL;
  
  fastbit_iapi_free_array_by_addr(query->dataSlice);

  if (s->_arrayName != NULL) {
    fastbit_iapi_free_array(s->_arrayName);
    free(s->_arrayName);
    s->_arrayName = NULL;
  }
  
  if (query->hasParent == 0) {
    fastbit_selection_free(s->_handle); 
  }   

  //s = NULL;
  //free(s);
}


void clear_fastbit_internal_recursive(ADIOS_QUERY* query) 
{
  clear_fastbit_internal(query);
  if (query->hasParent != 0) {
    fastbit_selection_free( ((FASTBIT_INTERNAL*)(query->queryInternal))->_handle); 
  }
 
  if (query->left != NULL) {
    clear_fastbit_internal(query->left);
  }
  if (query->right != NULL) {
    clear_fastbit_internal(query->right);
  }
}

ADIOS_QUERY* getFirstLeaf(ADIOS_QUERY* q);
void getHandle(int timeStep, int blockIdx, ADIOS_FILE* idxFile, ADIOS_QUERY* q);


void assertValue(char* input, char* endptr) {
  if (*endptr != '\0')  
    if ((errno == ERANGE) || (errno != 0) || (endptr == input) || (*endptr != '\0')) {
      //perror("strtol");
        printf("FastbitQueryEvaluation Exit due to :invalid type value: %s\n", input);
	exit(EXIT_FAILURE);
    }
}

long getMilliseconds() {
  time_t          s;  // Seconds
  struct timespec spec;

  clock_gettime(CLOCK_REALTIME, &spec);

  s  = spec.tv_sec;
  long ms = round(spec.tv_nsec / 1.0e6); // Convert nanoseconds to milliseconds

  return ms;
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
  ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = fastbit_selection_osr(arrayName, compareOp, vv);
  fastbit_adios_util_checkNotNull(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, arrayName);
}




void setCombinedQueryInternal(ADIOS_QUERY* q) {
    ADIOS_QUERY* left = (ADIOS_QUERY*)(q->left);
    ADIOS_QUERY* right = (ADIOS_QUERY*)(q->right);

    if (q->combineOp == ADIOS_QUERY_OP_AND) {
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle =  fastbit_selection_combine(((FASTBIT_INTERNAL*)left->queryInternal)->_handle, FastBitCombineAnd, 
										    ((FASTBIT_INTERNAL*)right->queryInternal)->_handle);
    } else {
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle =  fastbit_selection_combine(((FASTBIT_INTERNAL*)left->queryInternal)->_handle, FastBitCombineOr, 
										    ((FASTBIT_INTERNAL*)right->queryInternal)->_handle);
    }    
}

void getCoordinateFromPoints(uint64_t pos, const ADIOS_SELECTION_POINTS_STRUCT* sel, uint64_t* coordinates) 
{
  int i=0; 
  for(i=0; i<sel->ndim; i++) {
    (coordinates)[i] = sel->points[pos*(sel->ndim)+i];
  }
}


void getCoordinateFromBlock(uint64_t pos, const ADIOS_VARBLOCK* sel, int n, uint64_t* coordinates) 
{
  //log_debug("pos = %ld, n=%d \n", pos, n);
  if (n == 1) {
      coordinates[n-1] = pos + sel->start[n-1];
      return ;
  } 
 
  uint64_t lastDimSize= sel->count[n-1];     
  uint64_t res  = pos % lastDimSize;

  //log_debug("      lastDim = %ld, res=%d \n", lastDimSize, res);
  coordinates[n-1] = res + sel->start[n-1];
  uint64_t stepUp = (pos - res)/lastDimSize;

  //log_debug("      coordinate[%d]=%ld\n", n-1, coordinates[n-1]);
  getCoordinateFromBlock(stepUp, sel, n-1, coordinates);  
}

//check point coordinates
//offset in the bounding box needs to be taken account
void getCoordinateFromBox(uint64_t pos, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* sel, int n, uint64_t* coordinates) 
{
  //log_debug("pos = %ld, n=%d \n", pos, n);
  if (n == 1) {
      coordinates[n-1] = pos + sel->start[n-1];
      return ;
  } 
 
  uint64_t lastDimSize= sel->count[n-1];     
  uint64_t res  = pos % lastDimSize;

  //log_debug("      lastDim = %ld, res=%d \n", lastDimSize, res);
  coordinates[n-1] = res + sel->start[n-1];
  uint64_t stepUp = (pos - res)/lastDimSize;

  //log_debug("      coordinate[%d]=%ld\n", n-1, coordinates[n-1]);

  getCoordinateFromBox(stepUp, sel, n-1, coordinates);  
}

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

 

int64_t getPosInBox(const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* sel, int n, uint64_t* spatialCoordinates) 
{
  if (sel->ndim <= 0) {
    return -1;
  }

  int i=0;
  if (n == sel->ndim) { // check validation once
    for (i=0; i<sel->ndim; i++) {
      uint64_t min = sel->start[i];
      uint64_t max = sel->count[i] + min;
      if (spatialCoordinates[i] >= max) {
	return -1;
      }
      if (spatialCoordinates[i] < min) {
	return -1;
      } 
    }
  }

  // spatial Coordinate is valid
  if (n == 1) {
    return spatialCoordinates[0]-sel->start[0];
  }

  return (spatialCoordinates[n-1]-sel->start[n-1]) + sel->count[n-1]*getPosInBox(sel, n-1, spatialCoordinates);
}

int64_t getPosInVariable(const ADIOS_VARINFO* v, int n, uint64_t* spatialCoordinates) 
{
  if (v->ndim <= 0) {
    return -1;
  }
  
  int i=0;

  if (n == 1) {
    //log_debug("\n.....n=1 => %lld\n", spatialCoordinates[0]);
    return spatialCoordinates[0];
  }

  int64_t result =  spatialCoordinates[n-1] + v->dims[n-1]*getPosInVariable(v, n-1, spatialCoordinates); 
  //log_debug(".... n=%d, %lld + %ld * f[%d]\n", n, spatialCoordinates[n-1],v->dims[n-2], n-1);
  return result;
}


int64_t getRelativeIdxInBoundingBox(uint64_t currPosInBlock, const ADIOS_SELECTION_BOUNDINGBOX_STRUCT* bb, const ADIOS_VARBLOCK* blockSel)
{
    uint64_t spatialCoordinates[bb->ndim];
    getCoordinateFromBlock(currPosInBlock, blockSel, bb->ndim, spatialCoordinates);
    return getPosInBox(bb, bb->ndim, spatialCoordinates);
}

int64_t getRelativeIdxInVariable(uint64_t currPosInBlock, const ADIOS_VARINFO* v, const ADIOS_VARBLOCK* blockSel)
{
    uint64_t spatialCoordinates[v->ndim];
    getCoordinateFromBlock(currPosInBlock, blockSel, v->ndim, spatialCoordinates);
    //log_debug("   spatial=[%lld, %lld, %lld]    ", spatialCoordinates[0],spatialCoordinates[1],spatialCoordinates[2]);
    return getPosInVariable(v, v->ndim, spatialCoordinates);
}

int evaluateWithIdxOnBoundingBox(ADIOS_FILE* idxFile, ADIOS_QUERY* q, int timeStep)
{  
  ADIOS_SELECTION* sel = q->sel;
  ADIOS_VARINFO* v = q->varinfo;

  create_fastbit_internal(q);

  if (v == NULL) {
    ADIOS_QUERY* left = (ADIOS_QUERY*)(q->left);
    ADIOS_QUERY* right = (ADIOS_QUERY*)(q->right);

    if (evaluateWithIdxOnBoundingBox(idxFile, left, timeStep) < 0) {
      return -1;
    }
    if (evaluateWithIdxOnBoundingBox(idxFile, right, timeStep) < 0) {
      return -1;
    }

    setCombinedQueryInternal(q);
  } else {
    // is a leaf
    int blockStart=0; // relative
    int blockEnd = 0; // relative
    int i=0;

    const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = NULL;
    if (sel == NULL) {
      blockStart=0;
      blockEnd = v->nblocks[timeStep] -1;      
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
    }

    uint16_t* bitSlice = malloc((q->rawDataSize)* sizeof(uint16_t));

    for (i=0; i<q->rawDataSize; i++) {
      bitSlice[i] = 0;
    }

    uint64_t currBlockIdx = blockStart;

    uint64_t junk=0;
    for (currBlockIdx=blockStart; currBlockIdx <= blockEnd; currBlockIdx++) {
      getHandle(timeStep, currBlockIdx, idxFile, q);	      
      if (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle == 0) {
	log_error(">> Unable to construct fastbit query with NULL \n");
	return -1;
      }
      uint64_t count = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 
      
      junk += count;
      log_debug("block: %d hits = %lld, sum of hits so far: %lld\n", currBlockIdx, count, junk);
      //i = currBlockIdx-blockStart;
      uint64_t  coordinateArray[count];
      fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinateArray, count, 0);      

      //fastbit_selection_free(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);
      ////fastbit_iapi_free_array_by_addr(q->_dataSlice); // if attached index    
      //fastbit_iapi_free_array_by_addr(q->dataSlice); // if attached index       
      clear_fastbit_internal_recursive(q);

      int k=0;
      for (k=0; k<count; k++) {
	uint64_t currPosInBlock = coordinateArray[k];
	int64_t currPos = 0;
	if (bb != NULL) {
	  currPos = getRelativeIdxInBoundingBox(currPosInBlock, bb, &(v->blockinfo[currBlockIdx]));
	} else {
	  currPos = getRelativeIdxInVariable(currPosInBlock, v,  &(v->blockinfo[currBlockIdx]));
	}

	//log_warn("%lld th in block[%d],   =>  in actual %lld, limit: %ld \n", currPosInBlock, currBlockIdx, currPos, q->rawDataSize);
	//if ((currPos >= 0) && (currPos < q->rawDataSize)) {
	if (currPos >= 0) {
	  bitSlice[currPos] = 1;
	}
      }

      log_debug("----\n");
    }

    char bitsArrayName[50+strlen(q->condition)];
    sprintf(bitsArrayName, "%ld-%d-%s-%d", getMilliseconds(), v->varid, q->condition, timeStep);
    //return fastbit_selection_create(dataType, dataOfInterest, dataSize, compareOp, &vv);

    free(q->dataSlice);
    q->dataSlice = bitSlice;
    fastbit_iapi_register_array(bitsArrayName, FastBitDataTypeUShort, q->dataSlice, q->rawDataSize);
    //printData(q->_dataSlice, adios_unsigned_short, q->_rawDataSize); 
    //fastbit_selection_free(q->_queryInternal);
    ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = fastbit_selection_osr(bitsArrayName, FastBitCompareGreater, 0.5);
    fastbit_adios_util_checkNotNull(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, bitsArrayName);
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


void getHandleFromBlockAtLeafQuery(int timeStep, int blockIdx, ADIOS_FILE* idxFile, ADIOS_QUERY* q) 
{
  //double *keys = NULL; int64_t*offsets= NULL; uint32_t *bms = NULL;  
    uint64_t nk, no, nb;
    
    ADIOS_VARINFO* v = q->varinfo;
    
    ADIOS_FILE* dataFile = q->file;

    /*
    // read data from dataFile
    ADIOS_SELECTION* box = adios_selection_writeblock(blockIdx);
    common_read_inq_var_blockinfo(dataFile, v);
    */
    if (v->blockinfo == NULL) {
      common_read_inq_var_blockinfo(dataFile, v);
    }
    uint64_t blockSize = fastbit_adios_util_getBlockSize(v, blockIdx);

    /*
    free(q->_dataSlice);
    q->_dataSlice = malloc(adios_type_size(v->type, v->value)*blockSize);

    common_read_schedule_read_byid(dataFile, box, v->varid, timeStep, 1, NULL, q->_dataSlice);
    common_read_perform_reads(dataFile,1);
    
    //printData(q->_dataSlice, v->type, blockSize);
    */

    char blockDataName[40+strlen(q->condition)];
    sprintf(blockDataName, "%d-%s-%d-%d-%ld", v->varid, q->condition, timeStep, blockIdx, getMilliseconds());

    FASTBIT_INTERNAL* itn = (FASTBIT_INTERNAL*)(q->queryInternal);
    if (fastbit_adios_util_readFromIndexFile(idxFile, v, timeStep, blockIdx, &(itn->_keys), &nk, &(itn->_offsets), &no, &(itn->_bms), &nb) < 0) 
    {
      // no idx for this variable, read from file:
      free(q->dataSlice);
      q->dataSlice = malloc(adios_type_size(v->type, v->value)*blockSize);
      
      ADIOS_SELECTION* box = adios_selection_writeblock(blockIdx);   
      common_read_inq_var_blockinfo(dataFile, v);        
      common_read_schedule_read_byid(dataFile, box, v->varid, timeStep, 1, NULL, q->dataSlice);
      common_read_perform_reads(dataFile,1);

      FastBitDataType  dataType = fastbit_adios_util_getFastbitDataType(q->varinfo->type);
      FastBitCompareType compareOp = fastbit_adios_util_getFastbitCompareType(q->predicateOp);

      setQueryInternal(q, compareOp, dataType, blockSize, blockDataName);
      adios_selection_delete(box);
      return;
    }
    
    //int err = fastbit_iapi_register_array(blockDataName, fastbit_adios_util_getFastbitDataType(v->type), q->_dataSlice, blockSize);
    uint64_t nv = blockSize;
    itn->_arrayName = malloc(strlen(blockDataName)+2);
    sprintf(itn->_arrayName, "%s", blockDataName);
    
    int ierr = fastbit_iapi_register_array_index_only(blockDataName, fastbit_adios_util_getFastbitDataType(v->type), &nv, 1 , itn->_keys, nk, itn->_offsets, no, itn->_bms, mybmreader);

      /*
    if (ierr != 0) {
      log_error(" registering array failed. fastbit err code = %ld\n", ierr);
    }
    //printData(bms, adios_unsigned_integer, nb);

    ierr = fastbit_iapi_attach_index (blockDataName, keys, nk, offsets, no, bms, mybmreader);
      */
    if (ierr < 0) {
      log_error(" reattaching index failed. fastbit err code = %ld\n", ierr);
      //result = ierr;
    } else {
      ((FASTBIT_INTERNAL*)(q->queryInternal))->_handle = createHandle(q, blockDataName); //fastbit_selection_osr(blockDataName, getFastbitCompareType(q->_op), q->_value);
      fastbit_adios_util_checkNotNull(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, blockDataName);
    }

    //free(bms);
    //free(keys);
    //free(offsets);    

    //q->_queryInternal = result;
    //int64_t doubleC = fastbit_selection_evaluate(q->_queryInternal);
    //no need to return result here as we are using query internal to store
}



void printQueryData(ADIOS_QUERY* q, FastBitDataType dataType, int timeStep) {
  uint64_t dataSize = q->rawDataSize;
  int j;
  int batchSize = 31;
  log_debug ("::\t %s At timestep: %llu datasize=%llu \n\t\t   raw data:  [", q->condition, timeStep, dataSize);
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
	log_debug("%d ", ((uint64_t  *)(q->dataSlice))[j]);
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


int readWithTimeStepNoIdx(ADIOS_QUERY* q, int timeStep) {
  FastBitDataType  dataType = fastbit_adios_util_getFastbitDataType(q->varinfo->type);
  FastBitCompareType compareOp = fastbit_adios_util_getFastbitCompareType(q->predicateOp);

  int errorCode = common_read_schedule_read_byid (q->file, q->sel, q->varinfo->varid, 
                                                  timeStep, 1, NULL,  q->dataSlice);
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
  //adios_free_varinfo(q->_var);

  char datasetName[strlen(q->condition) + 40];
  sprintf(datasetName, "%s-%d-%ld", q->condition, timeStep, getMilliseconds());  
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
    getHandleFromBlockAtLeafQuery(timeStep, wb->index, idxFile, q);
  } else {
    blockSelectionFastbitHandle(idxFile, q->left, timeStep);
    blockSelectionFastbitHandle(idxFile, q->right, timeStep);

    setCombinedQueryInternal(q);
  }
}


void getHandle(int timeStep, int blockIdx, ADIOS_FILE* idxFile, ADIOS_QUERY* q) 
{
  ADIOS_FILE* dataFile = q->file;
  ADIOS_VARINFO* v = q->varinfo;

  FastBitSelectionHandle result = NULL;
  if (v == NULL) {
    getHandle(timeStep, blockIdx, idxFile, q->left);
    getHandle(timeStep, blockIdx, idxFile, q->right);

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
    getHandleFromBlockAtLeafQuery(timeStep, blockIdx, idxFile, q);
  }
      
}


int64_t  applyIndexIfExists (ADIOS_QUERY* q) 
{
  if (q->onTimeStep == gCurrentTimeStep) {
    //log_debug("::\t query index data has been read for timestep: %d\n", gCurrentTimeStep);
    return 0;
  } 

  // call fastbit_estimate_num_hits(selection)
  ADIOS_QUERY* leaf = getFirstLeaf(q);
  ADIOS_FILE* f = leaf->file;
  const char* basefileName = f->path;

  int64_t result = -1;
  MPI_Comm comm_dummy = MPI_COMM_SELF;
  ADIOS_FILE* idxFile = fastbit_adios_util_getFastbitIndexFileToRead(basefileName, comm_dummy);
    
  if (idxFile != NULL) {
      if ((leaf->sel == NULL) || (leaf->sel->type == ADIOS_SELECTION_BOUNDINGBOX)) {
	  evaluateWithIdxOnBoundingBox(idxFile,  q, gCurrentTimeStep);
	  result = fastbit_selection_estimate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);	
      } else if (leaf->sel->type == ADIOS_SELECTION_WRITEBLOCK) {
	  blockSelectionFastbitHandle(idxFile, q, gCurrentTimeStep);
	  result= fastbit_selection_estimate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle);       
      } 

      if (result > -1) {
	q->onTimeStep = gCurrentTimeStep;
	q->maxResultsDesired = 0;
	q->resultsReadSoFar = 0;

	//return result;
      } // otherwise, use no idx method
      common_read_close(idxFile);      
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
  /*
  int currentFileStep = leaf->file->current_step;
  if (currentFileStep > 0) {
    if (timestep != currentFileStep) {
      log_debug("timestep given %d is different from files: %d, correct to file time step.\n", timestep, currentFileStep);
      adios_query_set_timestep(currentFileStep);
      return 0;
    }
    return 0;
  } // == 0, can either be from read_open() or read_open_file(), cann't not distinguish
  */

  if (leaf->varinfo->nsteps <= timestep) {
    log_debug("timestep %d is more than variables limit: %d, can not evaluate.\n", timestep, leaf->varinfo->nsteps);
    return -1;
  }
  return 0;
}

int64_t adios_query_fastbit_estimate(ADIOS_QUERY* q) //, int timeStep) 
{
  if (q == NULL) {
    return -1;
  }

  int timeStep = gCurrentTimeStep;

  adios_query_fastbit_init();
  if (assertTimeStepValidWithQuery(q) != 0) {
    return -1;
  }

  create_fastbit_internal(q);

  int64_t estimate = applyIndexIfExists(q);
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
 

int64_t adios_query_fastbit_evaluate(ADIOS_QUERY* q, int timeStep, uint64_t _maxResult) 
{
  create_fastbit_internal(q);
  
  int64_t estimate = applyIndexIfExists(q);
  if (estimate < 0) {   // use no idx
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

  if ((q->queryInternal == 0) || (((FASTBIT_INTERNAL*)(q->queryInternal))->_handle == 0)) {
    log_error(">>  Unable to use fastbit to evaluate NULL query.\n"); 
    return -1;
  }
  int64_t numHits = fastbit_selection_evaluate(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle); 
  log_debug(":: ==> fastbit_evaluate() num of hits found for [%s] = %lld, at timestep %d \n", q->condition, numHits, timeStep);  

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
	log_debug("%d ", spatialCoordinates[k]);
      }
      log_debug("]\n");

}

void fillUp(int dimSize, uint64_t* spatialCoordinates, uint64_t i, uint64_t* pointArray) 
{
  int k=0;

  for (k = 0; k < dimSize; k++) {	  
    uint64_t idx = i * dimSize + k;
    //log_debug(" points[%d] = %lld ", idx, spatialCoordinates[k]);
    pointArray[idx] = spatialCoordinates[k];
  }	
  //log_debug("\n");
}

ADIOS_SELECTION* getSpatialCoordinatesDefault(ADIOS_VARINFO* var, uint64_t* coordinates, uint64_t retrivalSize)
{
  uint64_t arraySize = retrivalSize * (var->ndim);
  uint64_t* pointArray = (uint64_t*) (malloc(arraySize  * sizeof(uint64_t)));

  int i;
  for (i=0; i<retrivalSize; i++) {
    uint64_t spatialCoordinates[var->ndim];
    getCoordinateFromVariable(coordinates[i], var, var->ndim, spatialCoordinates);
    
    fillUp(var->ndim, spatialCoordinates, i, pointArray);
  }
  ADIOS_SELECTION* result =  common_read_selection_points(var->ndim, retrivalSize, pointArray);
  //free(pointArray); // user has to free this
  return result;
}

ADIOS_SELECTION* getSpatialCoordinates(ADIOS_SELECTION* outputBoundary, uint64_t* coordinates, uint64_t retrivalSize, ADIOS_VARINFO* v)
{
  int k = 0;
  uint64_t i=0;

  switch (outputBoundary->type) {
  case  ADIOS_SELECTION_BOUNDINGBOX:    
    {
      const ADIOS_SELECTION_BOUNDINGBOX_STRUCT *bb = &(outputBoundary->u.bb);

      uint64_t arraySize = retrivalSize * (bb->ndim);
      uint64_t* pointArray = (uint64_t*) (malloc(arraySize  * sizeof(uint64_t)));
      
      for (i=0; i<retrivalSize; i++) {
	   uint64_t spatialCoordinates[bb->ndim];
	   getCoordinateFromBox(coordinates[i], bb, bb->ndim, spatialCoordinates);

	   fillUp(bb->ndim, spatialCoordinates, i, pointArray);
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
	fillUp(points->ndim, spatialCoordinates, i, pointArray);
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
	   int absBlockCounter = query_utils_getGlobalWriteBlockId(wb->index, gCurrentTimeStep, v);
	   getCoordinateFromBlock(coordinates[i], &(v->blockinfo[absBlockCounter]), v->ndim, spatialCoordinates);

	   fillUp(v->ndim, spatialCoordinates, i, pointArray);
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

  if (q->varinfo != NULL) {
    return q;
  }
  return getFirstLeaf(q->left);
}

int  adios_query_fastbit_get_selection(ADIOS_QUERY* q, 
				      uint64_t batchSize, 
				      ADIOS_SELECTION* outputBoundary, 
				      ADIOS_SELECTION** result)
{
  /*
  if (q->_onTimeStep < 0) {
    log_error(":: Error: need to call evaluate first! Exit.\n");
    return -1;
  }
  */
  adios_query_fastbit_init();

  if (assertTimeStepValidWithQuery(q) != 0) {
    return -1;
  }

  adios_query_fastbit_evaluate(q, gCurrentTimeStep, 0);
  //log_debug("::\t max=%llu _lastRead=%llu\n", q->_maxResultDesired, q->_lastRead);
  if (batchSize == 0) {
    log_debug(":: ==> will not fetch. batchsize=0\n");
    return 0;
  }

  uint64_t retrivalSize = q->maxResultsDesired - q->resultsReadSoFar;
  if (retrivalSize == 0) {
     result = 0;
     log_debug(":: ==> no more results to fetch\n");
     return 0;
  }

  if (retrivalSize > batchSize) {
      retrivalSize = batchSize;
  }

  uint64_t coordinates[retrivalSize];

  fastbit_selection_get_coordinates(((FASTBIT_INTERNAL*)(q->queryInternal))->_handle, coordinates, retrivalSize, q->resultsReadSoFar);
    
  q->resultsReadSoFar += retrivalSize;
  
  if (outputBoundary == 0) {
    ADIOS_QUERY* firstLeaf = getFirstLeaf(q);
    if ((firstLeaf == NULL) || (firstLeaf->varinfo == NULL)) {
	log_error(":: Error: unable to get a valid first leaf! Exit. \n");
	return -1;
      }
    if (firstLeaf->sel == NULL) {
      *result = getSpatialCoordinatesDefault(firstLeaf->varinfo, coordinates, retrivalSize);
    } else {
      *result = getSpatialCoordinates(firstLeaf->sel, coordinates, retrivalSize, firstLeaf->varinfo);
    }
  } else {
    //*result = getSpatialCoordinates(outputBoundary, coordinates, retrivalSize);
    // variable needs to be in place to handle the block information
    // not sure wheather this is well defined case of combined query?! but the first varibale will be used for block information calculation
    *result = getSpatialCoordinates(outputBoundary, coordinates, retrivalSize, getFirstLeaf(q)->varinfo);

    if (*result == 0) {
      return -1;
    }
  }
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
  if (q->resultsReadSoFar == q->maxResultsDesired) {
    return 0;
  } else {
    return 1;
  }
}

//void adios_fastbit_free_query(ADIOS_QUERY* query) 
int  adios_query_fastbit_free(ADIOS_QUERY* query) 
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
  clear_fastbit_internal(query);
  free(query->queryInternal);

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
    
    adios_free_varinfo(v);

    ADIOS_QUERY* q = adios_query_create(f, varName, sel, getOp(opStr), valueStr);


    free(valueStr);free(opStr);free(varStr);free(varName);free(str);
    //free(start); free(count); // if deleted, then adios_sel values would be affected
    if (dimDef != NULL) 
      {free(dimDef);}

    log_debug("::\t query created for: %s\n", condition);
    return q;
}
*/
 /*
void queryDetail(ADIOS_QUERY* q, int timeStep) {
    int64_t estimated = adios_query_estimate(q);
    log_debug("::\t query estimated = %llu \n", estimated);

    int64_t numHits = adios_query_fastbit_evaluate(q, timeStep, 100000);
    log_debug("::\t query evaluated = %llu \n", numHits);

    uint64_t batchSize = 50;
    
    while (q->_maxResultDesired - q->_lastRead > 0) {      
      //log_debug("::\t max=%llu _lastRead=%llu\n", q->_maxResultDesired, q->_lastRead);
      ADIOS_SELECTION* t;
      ADIOS_SELECTION* bound;
      adios_query_get_selection(q, batchSize, bound, &t);
      free(t->u.points.points);
      adios_selection_delete(t);
      //log_debug("::\t      max=%llu _lastRead=%llu\n", q->_maxResultDesired, q->_lastRead);
    }
}
 */






