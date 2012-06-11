#ifndef UTIL_H_
#define UTIL_H_

#include <stdint.h>
#include "public/adios_types.h"
#include "public/adios_selection.h"

struct PairStruct {
    char * name;
    char * value;
    struct PairStruct * next;
};
typedef struct PairStruct PairStruct;

typedef struct read_request
{
    int rank; // the rank of processor which issued this request. 
    ADIOS_SELECTION * sel;
    int varid;
    int from_steps;
    int nsteps;
    void * data;
    void * priv; // private structure for each read method
/*
    int ndims;
    uint64_t * start;
    uint64_t * count;
    void * data;
    uint64_t size;
    int file_idx;
    uint64_t offset;
    void * parent;
*/
    struct read_request * next;
} read_request;

/* Process a ;-separated and possibly multi-line text and 
   create a list of name=value pairs from each 
   item which has a "name=value" pattern. Whitespaces are removed. 
   Input is not modified. Space is allocated;
   Also, simple "name" or "name=" patterns are processed and 
   returned with value=NULL. 
*/
PairStruct * text_to_name_value_pairs (const char * text);
void free_name_value_pairs (PairStruct * pairs);
/* Reverse the order in an array in place.
   use swapping from Fortran/column-major order to ADIOS-read-api/C/row-major order and back
*/
void swap_order(int n, uint64_t *array, int *timedim);
void change_endianness( void *data, uint64_t slice_size, enum ADIOS_DATATYPES type);
void copy_data (void *dst, void *src,
                int idim,
                int ndim,
                uint64_t* size_in_dset,
                uint64_t* ldims,
                const uint64_t * readsize,
                uint64_t dst_stride,
                uint64_t src_stride,
                uint64_t dst_offset,
                uint64_t src_offset,
                uint64_t ele_num,
                int      size_of_type
                );
void alloc_namelist (char ***namelist, int length);
void free_namelist (char **namelist, int length);
void list_insert_reader (read_request ** h, read_request * q);
ADIOS_SELECTION * copy_selection (ADIOS_SELECTION * sel);
#endif
