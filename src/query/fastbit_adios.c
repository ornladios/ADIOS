#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/common_read.h"
#include "core/adios_logger.h"
#include <iapi.h>

#include "fastbit_adios.h"

const char* _prefix="";

unsigned long _startMillis = 0;
unsigned long _milestoneMillis = 0;

//struct timespec _bmsVisits;
//struct timespec _idxVisits;

CollectionPoint _bmsPoint;
CollectionPoint _idxPoint;
CollectionPoint _proPoint; // fastbit_evaluate & getCoordinate
CollectionPoint _framePoint;

void casestudyLogger_initPoints(CollectionPoint* p) 
{
  p->_accumulator.tv_nsec = 0;
  p->_accumulator.tv_sec = 0;

  p->_counter = 0;
}


void casestudyLogger_init() 
{
  casestudyLogger_initPoints(&_bmsPoint);
  casestudyLogger_initPoints(&_idxPoint);
  casestudyLogger_initPoints(&_proPoint);
  casestudyLogger_initPoints(&_framePoint);
}

extern void casestudyLogger_print(CollectionPoint* p, const char* msg)
{
#ifdef TIMESTUDY
  printf("%s: %llu sec %llu millisec  accumulated, visited: %lu times \n",  msg,  (unsigned long long)p->_accumulator.tv_sec, (unsigned long long)p->_accumulator.tv_nsec/((unsigned long)1000000), p->_counter);
#else
  log_debug("%s: %llu sec %llu millisec  accumulated, visited: %lu times \n",  msg,  (unsigned long long)p->_accumulator.tv_sec, (unsigned long long)p->_accumulator.tv_nsec/((unsigned long)1000000), p->_counter);
#endif
}

extern void casestudyLogger_starts(const char* ref)
{
  _startMillis = fastbit_adios_getCurrentTimeMillis();
  _milestoneMillis = 0;

  /*_bmsVisits.tv_nsec=0;
  _bmsVisits.tv_sec=0;
  _idxVisits.tv_nsec=0;
  _idxVisits.tv_sec=0;
  */
#ifdef TIMESTUDY
  printf("==> casestudy %s starts at: %llu sec\n", ref, _startMillis/1000);
#endif
}


extern void casestudyLogger_bms_print()
{
  //printf("bitmap visit took: %llu sec %llu millisec  accumulated \n",  _bmsVisits.tv_sec, _bmsVisits.tv_nsec/((unsigned long)1000000));
  casestudyLogger_print(&_bmsPoint, "bitmap visit took");
}

extern void casestudyLogger_idx_print()
{
  //printf("idx access took: %llu sec %llu millisec  accumulated \n",  _idxVisits.tv_sec, _idxVisits.tv_nsec/((unsigned long)1000000));
  casestudyLogger_print(&_idxPoint, "idx access took:");
}

extern void casestudyLogger_pro_print()
{
  //printf("idx access took: %llu sec %llu millisec  accumulated \n",  _idxVisits.tv_sec, _idxVisits.tv_nsec/((unsigned long)1000000));
  casestudyLogger_print(&_proPoint, "fastbit query process took:");
}

extern void casestudyLogger_frame_print()
{
  //printf("idx access took: %llu sec %llu millisec  accumulated \n",  _idxVisits.tv_sec, _idxVisits.tv_nsec/((unsigned long)1000000));
  casestudyLogger_print(&_framePoint, "query frame took:");
}

extern void casestudyLogger_ends(const char* ref)
{
  unsigned long endMillis = fastbit_adios_getCurrentTimeMillis();
  unsigned long diffmillis = endMillis - _startMillis;

#ifdef TIMESTUDY
  unsigned long diffsec = diffmillis/1000;
  if (diffsec > 60) {
    unsigned long diffmin = diffsec/60;
    diffsec = diffsec - diffmin*60;
    printf("%s took: %lu millisecs = %lu min & %lu sec \n", ref, diffmillis, diffmin, diffsec);
  } else {
    printf("%s took: %lu millisecs = %lu sec\n", ref, diffmillis, diffsec);
  }
#endif
}

extern void casestudyLogger_getRealtime(struct timespec* spec)
{
#ifdef CLOCK_MONOTONIC
  clock_gettime(CLOCK_MONOTONIC, spec);
#else
  clock_gettime(CLOCK_REALTIME, spec);
#endif
}

extern void casestudyLogger_setPrefix(const char* prefix)
{
  _prefix = prefix;

  unsigned long currMillis = fastbit_adios_getCurrentTimeMillis();
  unsigned long diffMillis = currMillis-_startMillis;
  unsigned long diffsec = diffMillis/1000;

  unsigned long ministepMillis = currMillis - _milestoneMillis;
  if (_milestoneMillis == 0) {
    ministepMillis = diffMillis;
  } 
    

#ifdef TIMESTUDY
  if (diffsec > 60) {
    unsigned long diffmin = diffsec/60;
    diffsec = diffsec - diffmin*60;
    printf("on %s, time passed since start: %lu millisecs = %lu min & %lu sec, ministep=%lu milliseconds \n", prefix, diffMillis, diffmin, diffsec, ministepMillis);
  } else {
    printf("on %s, time passed since start: %lu millisecs = %lu sec, ministep=%lu milliseconds\n", prefix, diffMillis, diffsec, ministepMillis);
  }
#endif
  _milestoneMillis = currMillis;
}

/*
extern void casestudyLogger_idx_writeout(struct timespec* start, 
					 //struct timespec* end,
					 const char* desc)
{
  struct timespec end; 
  casestudyLogger_getRealtime(&end);

  struct timespec diff; 

  diff.tv_nsec = end.tv_nsec - start->tv_nsec; 
  diff.tv_sec = end.tv_sec - start->tv_sec; 

  if(diff.tv_nsec < 0 ){
    diff.tv_nsec += 1000000000; 
    diff.tv_sec = diff.tv_sec - 1;
  } 

  //printf("%s %s %llu sec %llu nanosec \n", _prefix, desc, diff.tv_sec, diff.tv_nsec);

  _idxVisits.tv_sec  += diff.tv_sec;
  _idxVisits.tv_nsec += diff.tv_nsec;

  unsigned long sec = _idxVisits.tv_nsec/(unsigned long)1000000000;
  if (sec > 0) {
    _idxVisits.tv_sec += sec;
    _idxVisits.tv_nsec -= sec*((unsigned long)1000000000);
  } 

}


extern void casestudyLogger_bms_writeout(struct timespec* start, 
                              //struct timespec* end,
			      const char* desc)
{
  struct timespec end; 
  casestudyLogger_getRealtime(&end);

  struct timespec diff; 

  diff.tv_nsec = end.tv_nsec - start->tv_nsec; 
  diff.tv_sec = end.tv_sec - start->tv_sec; 

  if(diff.tv_nsec < 0 ){
    diff.tv_nsec += 1000000000; 
    diff.tv_sec = diff.tv_sec - 1;
  } 

  //printf("%s %s %llu sec %llu nanosec \n", _prefix, desc, diff.tv_sec, diff.tv_nsec);

  _bmsVisits.tv_sec  += diff.tv_sec;
  _bmsVisits.tv_nsec += diff.tv_nsec;

  unsigned long sec = _bmsVisits.tv_nsec/(unsigned long)1000000000;
  if (sec > 0) {
    _bmsVisits.tv_sec += sec;
    _bmsVisits.tv_nsec -= sec*((unsigned long)1000000000);
  } 
  //printf("%s %s %llu sec %llu millisec  accumulated \n", _prefix, desc, _bmsVisits.tv_sec, _bmsVisits.tv_nsec/((unsigned long)1000000));
}
*/

extern void casestudyLogger_writeout(CollectionPoint* p,
				     struct timespec* start, 
				     const char* desc)
{
  struct timespec end; 
  casestudyLogger_getRealtime(&end);

  struct timespec diff; 

  diff.tv_nsec = end.tv_nsec - start->tv_nsec; 
  diff.tv_sec = end.tv_sec - start->tv_sec; 

  if(diff.tv_nsec < 0 ){
    diff.tv_nsec += 1000000000; 
    diff.tv_sec = diff.tv_sec - 1;
  } 

  //printf("%s %s %llu sec %llu nanosec \n", _prefix, desc, diff.tv_sec, diff.tv_nsec);

  p->_accumulator.tv_sec  += diff.tv_sec;
  p->_accumulator.tv_nsec += diff.tv_nsec;

  unsigned long sec = p->_accumulator.tv_nsec/(unsigned long)1000000000;
  if (sec > 0) {
    p->_accumulator.tv_sec += sec;
    p->_accumulator.tv_nsec -= sec*((unsigned long)1000000000);
  } 
  //printf("%s %s %llu sec %llu millisec  accumulated \n", _prefix, desc, _bmsVisits.tv_sec, _bmsVisits.tv_nsec/((unsigned long)1000000));
  
  p->_counter ++;
}

extern void casestudyLogger_idx_writeout(struct timespec* start, 
					 //struct timespec* end,
					 const char* desc)
{
  casestudyLogger_writeout(&_idxPoint, start, desc);
}

extern void casestudyLogger_pro_writeout(struct timespec* start, 
					 //struct timespec* end,
					 const char* desc)
{
  casestudyLogger_writeout(&_proPoint, start, desc);
}

extern void casestudyLogger_frame_writeout(struct timespec* start, 
					 //struct timespec* end,
					 const char* desc)
{
  casestudyLogger_writeout(&_framePoint, start, desc);
  //printf("%s %llu sec %llu millisec  accumulated \n", desc, _framePoint._accumulator.tv_sec, _framePoint._accumulator.tv_nsec/((unsigned long)1000000));
}

extern void casestudyLogger_bms_writeout(struct timespec* start, 
					 //struct timespec* end,
					 const char* desc)
{
  casestudyLogger_writeout(&_bmsPoint, start, desc);
}


//
//
//
long fastbit_adios_getCurrentTimeMillis() 
{
  time_t          s;  // Seconds
  struct timespec spec;

#ifdef CLOCK_MONOTONIC
  clock_gettime(CLOCK_MONOTONIC, &spec);
#else
  clock_gettime(CLOCK_REALTIME, &spec);
#endif
  s  = spec.tv_sec;
  long ms = round(spec.tv_nsec/1.0e6) + s*1000; // Convert nanoseconds to milliseconds
  return ms;
}


void fastbit_adios_util_checkNotNull(void* fastbitHandle, const char* arrayName) {
  if (fastbitHandle == NULL) {
    log_error(" >> Unable to create handle on fastbit, ref: %s\n", arrayName);
  } else {
    log_debug(" >> created handle on fastbit, ref: %s\n", arrayName);
  }
}

int fastbit_adios_util_getRelativeBlockNumForPoint(ADIOS_VARINFO* v,  uint64_t* point, int timestep) 
{
  int i=0;
  int j=0;

  if (v->nsteps == 1) {
    // 
    // if file is read through adios_read_open(), knows only about the current time step
    //
    int result = -1;    
    for (i=0; i<v->sum_nblocks; i++) {      
      if (result >=0 ) {	
	break;
      }
	
      ADIOS_VARBLOCK curr = v->blockinfo[i];
      for (j=0; j<v->ndim; j++) {
	uint64_t begin = curr.start[j];
	uint64_t end   = curr.start[j]+curr.count[j];

	if ((begin <= point[j]) && (point[j] < end)) {
	  result = i; // relative to the timestep                                                                                                                                      
	  //result = sum; //if want to return abs block number                                                                                                                         
	} else {
	  result = -1;
	  break;
	}
      }
    }
    return result;
  }

  //
  // if v has info of all time steps;
  //
  int totalBlocksInTimeStep = v->nblocks[timestep];
  int sum=0;

  for (i=0;i<timestep; i++) {
    sum += v->nblocks[i];
  }

  int result = -1;
  for (i=0; i<totalBlocksInTimeStep; i++) {
    if (result >= 0) {
      break;
    }
    ADIOS_VARBLOCK curr = v->blockinfo[sum];
    for (j=0; j<v->ndim; j++) {
      int begin = curr.start[j];
      int end   = curr.start[j]+curr.count[j];

      if ((begin <= point[j]) && (point[j] < end)) {
	result = i; // relative to the timestep
	//result = sum; //if want to return abs block number
      } else {
	result = -1;
	break;
      }
    }
    sum++;
  }

  return result;
}


char *fastbit_adios_util_getFastbitIndexFileName(const char* dataFileLoc) 
{
  int len = strlen(dataFileLoc);			   
  char  idxFileNamePad [len];
  char *idxFileName = malloc (len*sizeof(char)+2); // .idx is 1 more than .bp
  
  strncpy(idxFileNamePad, dataFileLoc, len-3);
  idxFileNamePad[len-3]=0;
  sprintf(idxFileName, "%s.idx", idxFileNamePad); 
  return idxFileName;
}

int fastbit_adios_util_FastbitIndexFileExists(const char* dataFileLoc)
{
    char *idxFileName = fastbit_adios_util_getFastbitIndexFileName (dataFileLoc);
    int retval = query_utils_file_exists (idxFileName);
    free (idxFileName);
    return retval;
}

ADIOS_FILE* fastbit_adios_util_getFastbitIndexFileToRead(const char* dataFileLoc, MPI_Comm comm) 
{
  char *idxFileName = fastbit_adios_util_getFastbitIndexFileName (dataFileLoc);
  ADIOS_FILE* f = NULL;

  if (query_utils_file_exists(idxFileName)) 
  {
      // turn off logging
      //int log_level = adios_verbose_level;
      //adios_verbose_level = 0;
      f = common_read_open_file (idxFileName, ADIOS_READ_METHOD_BP, comm);
      // turn back on logging
      //adios_verbose_level = log_level;
      // reset possible error
      //adios_clear_error();
      if (!f) {
          log_warn ("Could not open FastBit index file '%s'. "
                  "Will use FastBit evaluation on data directly\n",
                  idxFileName);
      }
  } else {
      log_warn ("No FastBit index file '%s' was found. "
              "Will use FastBit evaluation on data directly\n",
              idxFileName);
  }
  free (idxFileName);
  return f;
}




FastBitDataType fastbit_adios_util_getFastbitDataType(enum ADIOS_DATATYPES type) 
{  
  switch (type)
    {
    case adios_unsigned_byte:
      return FastBitDataTypeUByte;
      break;

    case adios_byte:
      return FastBitDataTypeByte;
      break;

    case adios_short:
      return FastBitDataTypeShort;
      break;

    case adios_unsigned_short:
      return FastBitDataTypeUShort;
      break;

    case adios_integer:
      return FastBitDataTypeInt;
      break;

    case adios_unsigned_integer:
      return FastBitDataTypeUInt;
      break;

    case adios_long:
      return FastBitDataTypeLong;
      break;

    case adios_unsigned_long:
      return FastBitDataTypeULong;
      break;

    case adios_string:
      return FastBitDataTypeUnknown;
      break;

    case adios_real:
      return FastBitDataTypeFloat;
      break;

    case adios_double:
      return FastBitDataTypeDouble;
      break;

    case adios_long_double:
    //sprintf (s, "%Lg", ((long double *) data)[idx]);
    case adios_complex:
    //sprintf (s, "(%g, %g)", ((float *) data)[2*idx], ((float *) data)[2*idx+1]);   	       
    case adios_double_complex:
    //sprintf (s, "(%lg, %lg)", ((double *) data)[2*idx], ((double *) data)[2*idx+1]);	       
    return FastBitDataTypeDouble;
    }

}

FastBitCompareType fastbit_adios_util_getFastbitCompareType(enum ADIOS_PREDICATE_MODE op) 
{
    switch (op) 
    {
    case ADIOS_LT:
      return FastBitCompareLess;
      break;
    case ADIOS_LTEQ:
      return FastBitCompareLessEqual;
      break;
    case ADIOS_GT:
      return FastBitCompareGreater;
      break;
    case ADIOS_GTEQ:
      return FastBitCompareGreaterEqual;
      break;
    case ADIOS_EQ:
      return FastBitCompareEqual;
      break;
    case ADIOS_NE:
      return FastBitCompareNotEqual;
      break;
    }
}

// k is numbered from 1 to sum_nblocks
//uint64_t getBlockDataSize(ADIOS_VARINFO* v, int k) // k = blockNumber 
uint64_t fastbit_adios_util_getBlockSize(ADIOS_VARINFO* v, int timestep, int relativeBlockIdx) // k = blockNumber 
{
  //uint64_t blockBytes = common_read_type_size (v->type, v->value);
  uint64_t blockSize = 1;
  int j=0;

  if (v->ndim <= 0) {
    return blockSize;
  }

  int k = 0; // absBlockIdx
  for (j=0; j<timestep; j++) {
    k+= v->nblocks[j];
  }

  k += relativeBlockIdx;
  
  if (k > v->sum_nblocks-1) {
    return 0;
  }

  log_debug("\n blockinfo[%d]: [ ", k);
  
  for (j=0; j<v->ndim; j++) 
    {  
      blockSize *= v->blockinfo[k].count[j];
      log_debug("%llu:%llu ", v->blockinfo[k].start[j], v->blockinfo[k].count[j]);
    }
  
  log_debug("]\n");
  
  //  log_debug("\t\t   block %d, bytes: %llu \n", k, blockBytes);      
  
  return blockSize;
}


static const char * value_to_string (enum ADIOS_DATATYPES type, void * data, int idx)
{
    static char s [100];
    s [0] = 0;


  switch (type)
    {
    case adios_unsigned_byte:
      sprintf (s, "%u", ((uint8_t *) data)[idx]);
      break;

    case adios_byte:
      sprintf (s, "%d", ((int8_t *) data)[idx]);
      break;

    case adios_short:
      sprintf (s, "%hd", ((int16_t *) data)[idx]);
      break;

    case adios_unsigned_short:
      sprintf (s, "%hu", ((uint16_t *) data)[idx]);
      break;

    case adios_integer:
      sprintf (s, "%d", ((int32_t *) data)[idx]);
      break;

    case adios_unsigned_integer:
      sprintf (s, "%u", ((uint32_t *) data)[idx]);
      break;
      
    case adios_long:
      sprintf (s, "%lld", ((int64_t *) data)[idx]);
      break;
      
    case adios_unsigned_long:
      sprintf (s, "%llu", ((uint64_t *) data)[idx]);
      break;
      
    case adios_real:
      sprintf (s, "%g", ((float *) data)[idx]);
      break;
      
    case adios_double:
      sprintf (s, "%lg", ((double *) data)[idx]);
      break;
      
    case adios_long_double:
      sprintf (s, "%Lg", ((long double *) data)[idx]);
      break;
      
    case adios_string:
      return (char*) ((char *)data+idx);
      break;

    case adios_complex:
      sprintf (s, "(%g, %g)",
	       ((float *) data)[2*idx], ((float *) data)[2*idx+1]);
      break;
      
    case adios_double_complex:
      sprintf (s, "(%lg, %lg)",
	       ((double *) data)[2*idx], ((double *) data)[2*idx+1]);
      break;
    }

  return s;
}

void getVarNames(char* bmsVarName, char* keyVarName, char* offsetName, ADIOS_VARINFO* v, int timestep, int blockNum) {

  /*
  sprintf(bmsVarName, "bms-%d-%d-block-%d", v->varid, timestep, blockNum);
  sprintf(keyVarName, "key-%d-%d-block-%d", v->varid, timestep, blockNum);
  sprintf(offsetName, "offset-%d-%d-block-%d", v->varid, timestep, blockNum);
  */

  sprintf(bmsVarName, "bms-%d-%d-box-%d", v->varid, timestep, blockNum);
  sprintf(keyVarName, "key-%d-%d-box-%d", v->varid, timestep, blockNum);
  sprintf(offsetName, "offset-%d-%d-box-%d", v->varid, timestep, blockNum);

  /*
  sprintf(bmsVarName, "bms-%d-%d-%d", v->varid, timestep, blockNum);
  sprintf(keyVarName, "key-%d-%d-%d", v->varid, timestep, blockNum);
  sprintf(offsetName, "offset-%d-%d-%d", v->varid, timestep, blockNum);
  */

}
//
//
// caller frees keys, offsets and bms.
//
//
int fastbit_adios_util_readNoBMSFromIndexFile(ADIOS_FILE* idxFile, ADIOS_VARINFO* v, int timestep, int blockNum, 
					      double** keys, uint64_t* nk, int64_t** offsets, uint64_t* no,
					      char** bmsVarName)
{
  char keyVarName[100];
  char offsetName[100];

  *bmsVarName = malloc(100);

  getVarNames(*bmsVarName, keyVarName, offsetName, v, timestep, blockNum);

  /*
  sprintf(*bmsVarName, "bms-%d-%d-%d", v->varid, timestep, blockNum);
  sprintf(keyVarName, "key-%d-%d-%d", v->varid, timestep, blockNum);
  sprintf(offsetName, "offset-%d-%d-%d", v->varid, timestep, blockNum);
  */
  log_debug("reading from index file: %s for variables: %s %s %s \n", idxFile->path, *bmsVarName, keyVarName, offsetName);

  ADIOS_VARINFO * keyV = common_read_inq_var (idxFile, keyVarName);
  ADIOS_VARINFO * offsetV = common_read_inq_var (idxFile, offsetName);

  if ((keyV == 0) || (offsetV == 0)) {
    log_warn("WARN: no index for this variable.\n");
    return -1;
  }

  uint64_t keySize    = (uint64_t)(keyV->dims[0] * common_read_type_size (keyV->type, keyV->value));
  uint64_t offsetSize = (uint64_t)(offsetV->dims[0] * common_read_type_size (offsetV->type, offsetV->value));
    
  *offsets = malloc(offsetSize);
  *keys    = malloc(keySize);

  uint64_t start[] = {0};
  uint64_t count_key[] = {keyV->dims[0]};
  uint64_t count_offset[] = {offsetV->dims[0]};

  ADIOS_SELECTION* keySel = common_read_selection_boundingbox(keyV->ndim, start, count_key);
  ADIOS_SELECTION* offsetSel = common_read_selection_boundingbox(offsetV->ndim, start, count_offset);

  // idx file has one timestep
  common_read_schedule_read(idxFile, keySel, keyVarName, 0, 1, NULL, *keys);
  common_read_schedule_read(idxFile, offsetSel, offsetName, 0, 1, NULL, *offsets);

  common_read_perform_reads(idxFile,1);

  *nk = keyV->dims[0];
  *no = offsetV->dims[0];

  log_debug(" /key/offset data: length=%lld/%lld\n", *nk, *no);
  
  //printData(*bms, bmsV->type, *nb);
  common_read_selection_delete(keySel);
  common_read_free_varinfo(keyV);

  common_read_selection_delete(offsetSel);
  common_read_free_varinfo(offsetV);

  return 0;
}

int fastbit_adios_util_readFromIndexFile(ADIOS_FILE* idxFile, ADIOS_VARINFO* v, int timestep, int blockNum, 
					 double** keys, uint64_t* nk, int64_t** offsets, uint64_t* no,
					 uint32_t** bms, uint64_t* nb)

{
  char bmsVarName[100];
  char keyVarName[100];
  char offsetName[100];

  getVarNames(bmsVarName, keyVarName,  offsetName, v, timestep, blockNum);
  /*
  sprintf(bmsVarName, "bms-%d-%d-%d", v->varid, timestep, blockNum);
  sprintf(keyVarName, "key-%d-%d-%d", v->varid, timestep, blockNum);
  sprintf(offsetName, "offset-%d-%d-%d", v->varid, timestep, blockNum);
  */
  log_debug("reading from index file: %s for variables: %s %s %s \n", idxFile->path, bmsVarName, keyVarName, offsetName);

  ADIOS_VARINFO * bmsV = common_read_inq_var (idxFile, bmsVarName);
  ADIOS_VARINFO * keyV = common_read_inq_var (idxFile, keyVarName);
  ADIOS_VARINFO * offsetV = common_read_inq_var (idxFile, offsetName);

  if ((bmsV == 0) || (keyV == 0) || (offsetV == 0)) {
    log_warn("WARN: no index for this variable.\n");
    return -1;
  }

  uint64_t bmsSize    = (uint64_t)(bmsV->dims[0] * common_read_type_size (bmsV->type, bmsV->value));
  uint64_t keySize    = (uint64_t)(keyV->dims[0] * common_read_type_size (keyV->type, keyV->value));
  uint64_t offsetSize = (uint64_t)(offsetV->dims[0] * common_read_type_size (offsetV->type, offsetV->value));
    
  *offsets = malloc(offsetSize);
  *keys    = malloc(keySize);
  *bms     = malloc(bmsSize);

  uint64_t start[] = {0};
  uint64_t count_bms[] = {bmsV->dims[0]};
  uint64_t count_key[] = {keyV->dims[0]};
  uint64_t count_offset[] = {offsetV->dims[0]};

  ADIOS_SELECTION* bmsSel = common_read_selection_boundingbox(bmsV->ndim, start, count_bms);
  ADIOS_SELECTION* keySel = common_read_selection_boundingbox(keyV->ndim, start, count_key);
  ADIOS_SELECTION* offsetSel = common_read_selection_boundingbox(offsetV->ndim, start, count_offset);

  // idx file has one timestep
  common_read_schedule_read(idxFile, bmsSel, bmsVarName, 0, 1, NULL, *bms);
  common_read_schedule_read(idxFile, keySel, keyVarName, 0, 1, NULL, *keys);
  common_read_schedule_read(idxFile, offsetSel, offsetName, 0, 1, NULL, *offsets);

  common_read_perform_reads(idxFile,1);

  *nk = keyV->dims[0];
  *no = offsetV->dims[0];
  *nb = bmsV->dims[0];

  log_debug(" bms/key/offset data: length=%lld/%lld/%lld\n", *nb, *nk, *no);
  
  //printData(*bms, bmsV->type, *nb);
  common_read_selection_delete(bmsSel);
  common_read_free_varinfo(bmsV);

  common_read_selection_delete(keySel);
  common_read_free_varinfo(keyV);

  common_read_selection_delete(offsetSel);
  common_read_free_varinfo(offsetV);

  return 0;
}
 
void fastbit_adios_util_printData(void* data, enum ADIOS_DATATYPES type, uint64_t size)
{
  /*
see if blocks are read by bounding boxes as in blockinfo
or is lined as 
		timestep n (1,...nblocks)
		do not see how i (in sum_nblocks) will reflect timestep...
*/

  int i=0;
  int max = 10;
  if (max > size) {
    max = size;
  }
  log_debug("  \tfirst %d data out of %lld:[", max, size);
  for (i=0; i<max; i++) {
    log_debug("%s ", value_to_string(type, data, i));
  }
  log_debug("]\n");
}
