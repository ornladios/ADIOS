#ifndef _BINPACK_GENERAL_H
#define _BINPACK_GENERAL_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef BP_USEUNSIGNED

#  if defined(BP_USELONGLONG)
typedef unsigned long long bp_ulonglong_t;
#define BP_FMTSTR "%llu"
#  elif defined(BP_USELONG)
typedef unsigned long bp_ulonglong_t;
#define BP_FMTSTR "%lu"
#  else
typedef unsigned int bp_ulonglong_t;
#define BP_FMTSTR "%u"
#  endif

typedef unsigned int bp_uint_t;

#else // No BP_USEUNSIGNED

#  if defined(BP_USELONGLONG)
typedef long long bp_ulonglong_t;
#define BP_FMTSTR "%lld"
#  elif defined(BP_USELONG)
typedef long bp_ulonglong_t;
#define BP_FMTSTR "%ld"
#  else
typedef int bp_ulonglong_t;
#define BP_FMTSTR "%d"
#  endif

typedef int bp_uint_t;

#endif // BP_USEUNSIGNED

typedef unsigned char bp_uchar_t;

#define BP_ENDTAG 1000
#define BP_NUMTAGS 15
// BP_NUMTAGS must be updated if new tags are added!

// bit-wise OR with rank to signify a global array
#define ADIOS_GLOBAL_ARRAY_FLAG 0x80000000
#define ADIOS_GLOBAL_ARRAY_REMOVE 0x7FFFFFFF

enum TAG_t {NULL_TAG = 0
           ,GRP_TAG = 1
           ,DST_TAG = 2
           ,SCR_TAG = 3
           ,DIR_TAG = 4
           ,VAL_TAG = 5
           ,DSTVAL_TAG = 6
           ,DSTATRN_TAG = 7
           ,DSTATRS_TAG = 8
           ,GRPATRN_TAG = 9
           ,GRPATRS_TAG = 10
           };

enum vartype_t {bp_char = 0
               ,bp_short = 1
               ,bp_int = 2
               ,bp_long = 3
               ,bp_longlong = 4
               ,bp_float = 5
               ,bp_double = 6
               ,bp_longdouble = 7
               ,bp_pointer = 8
               ,bp_string = 9
               ,bp_complex = 10
               ,bp_uchar = 50
               ,bp_ushort = 51
               ,bp_uint = 52
               ,bp_ulong = 53
               ,bp_ulonglong = 54
               };

struct adios_bp_dimension_struct
{
    int local_bound;
    int global_bound;
    int global_offset;
};

struct adios_bp_element_struct
{
    enum TAG_t tag;
    int size;
    char * name;
    char * path;
    enum vartype_t type;
    int ranks;
    struct adios_bp_dimension_struct * dims;
    void * data;
};

// _BINPACK_GENERAL_H
#endif
