#ifndef __FASTBIT_ADIOS_H__
#define __FASTBIT_ADIOS_H__

#ifdef __cplusplus
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adios_read.h"
#include <iapi.h>

#include "adios_query.h"

FastBitDataType getFastbitDataType(enum ADIOS_DATATYPES type);

FastBitCompareType getFastbitCompareType(enum ADIOS_PREDICATE_MODE op);

#ifdef __cplusplus
}
#endif

#endif /* __ADIOS_QUERY_H__ */
