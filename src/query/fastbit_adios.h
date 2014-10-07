#ifndef __FASTBIT_ADIOS_H__
#define __FASTBIT_ADIOS_H__

#ifdef __cplusplus
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "public/adios_read.h"
#include <iapi.h>

#include "public/adios_query.h"

/** A simple reader to be used by FastBit for index reconstruction.  In
    this simple case, the first argument is the whole array storing all the
    serialized bitmaps.  This first argument can be used to point to a data
    structure pointing to any complex object type necassary.
*/
//
// this static function is from fastbit example/tiapi.c
//
static int mybmreader(void *ctx, uint64_t start,uint64_t count, uint32_t *buf)
{
  const uint32_t *bms = (uint32_t*)ctx + start;
  unsigned j;
  for (j = 0; j < count; ++ j) {
    buf[j] = bms[j];
  }
  return 0;
}

int fastbit_adios_util_getRelativeBlockNumForPoint(ADIOS_VARINFO* v,  uint64_t* point, int timestep);
void fastbit_adios_util_checkNotNull(void* fastbitHandle, const char* arrayName);
  
char* fastbit_adios_util_getFastbitIndexFileName(const char* dataFileLoc);
int fastbit_adios_util_FastbitIndexFileExists(const char* dataFileLoc);
ADIOS_FILE* fastbit_adios_util_getFastbitIndexFileToRead(const char* dataFileLoc, MPI_Comm comm);

void fastbit_adios_util_printData(void* data, enum ADIOS_DATATYPES type, uint64_t size);

FastBitDataType fastbit_adios_util_getFastbitDataType(enum ADIOS_DATATYPES type);

FastBitCompareType fastbit_adios_util_getFastbitCompareType(enum ADIOS_PREDICATE_MODE op);

const char * fastbit_adios_util_value_to_string (enum ADIOS_DATATYPES type, void * data, int idx);

uint64_t fastbit_adios_util_getAdiosBlockSize(ADIOS_VARINFO* v, int k); // k = blockNumber;

  
int fastbit_adio_util_readFromFastbitIndexFile(ADIOS_FILE* idxFile, ADIOS_VARINFO* v, int timestep, int blockNum, 
			     double** keys, uint64_t* nkeys, int64_t** offsets, uint64_t* no,
			     uint32_t** bms, uint64_t* nb);
  
#ifdef __cplusplus
}
#endif

#endif /* __ADIOS_QUERY_H__ */
