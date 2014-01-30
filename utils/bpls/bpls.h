/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */


#include <stdint.h>
#include "adios_read.h"

/* definitions for bpls.c */
#define myfree(p) if (p) { free(p); p=NULL; }

#ifndef HAVE_STRNDUP
#  define strndup(str,len) strdup(str)
#endif

typedef int bool;
#define false 0
#define true  1

#define CUT_TO_BYTE(x) (x < 0 ? 0 : (x > 255 ? 255 : x))

#define MAX_DIMS 16
#define MAX_MASKS 10
#define MAX_BUFFERSIZE (10*1024*1024)

// how to print one data item of an array
//enum PrintDataType {STRING, INT, FLOAT, DOUBLE, COMPLEX};

void init_globals(void);
void processDimSpecs(void);
void parseDimSpec(char *str, int *dims);
int compile_regexp_masks(void);
void printSettings(void);
int  doList(const char *path);
void mergeLists(int nV, char **listV, int nA, char **listA, char **mlist, bool *isVar);
int  readVar(ADIOS_FILE *fp, ADIOS_VARINFO *vi, const char *name, bool timed);
int cmpstringp(const void *p1, const void *p2);
bool grpMatchesMask(char *name);
bool matchesAMask(char *name);
int  print_start(const char *fname);
void print_slice_info(int ndim, uint64_t *dims, int timed, int nsteps, uint64_t *s, uint64_t *c);
int print_data(void *data, int item, enum ADIOS_DATATYPES adiosvartype, bool allowformat);
int print_dataset(void *data, enum ADIOS_DATATYPES adiosvartype, 
               uint64_t *s, uint64_t *c, int tdims, int* ndigits);  
void print_endline(void);
void print_stop(void);
int print_data_hist(ADIOS_VARINFO * vi, char * varname);
int print_data_characteristics(void * min, void * max, double * avg, double * std_dev, enum ADIOS_DATATYPES adiosvartype, bool allowformat);
void print_decomp(ADIOS_VARINFO *vi);
