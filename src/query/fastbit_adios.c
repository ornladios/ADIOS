#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adios_read.h"
#include <iapi.h>

#include "fastbit_adios.h"

FastBitDataType getFastbitDataType(enum ADIOS_DATATYPES type) 
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

FastBitCompareType getFastbitCompareType(enum ADIOS_PREDICATE_MODE op) 
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
