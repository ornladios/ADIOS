#ifndef UTIL_H_
#define UTIL_H_

#include <stdint.h>
#include <inttypes.h>
#include "public/adios_types.h"
#include "public/adios_selection.h"
#include "core/a2sel.h"

typedef struct read_request
{
    ADIOS_SELECTION * sel;
    int varid;
    int from_steps;
    int nsteps;
    void * data;
    uint64_t datasize; // size of selection to hold data
// above is the common fields that all read method will use
    void * priv; // private structure for each read method
    struct read_request * next;
} read_request;

/* Reverse the order in an array in place.
   use swapping from Fortran/column-major order to ADIOS-read-api/C/row-major order and back
*/
void swap_order(int n, uint64_t *array, int *timedim);
void change_endianness( void *data, uint64_t slice_size, enum ADIOS_DATATYPES type);

/* Copy data from one n-dimensional block to another, where the two blocks logically somewhat overlap.
 */
void adios_util_copy_data (void *dst, void *src,
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
                           int      size_of_type,
                           enum ADIOS_FLAG change_endiness,
                           enum ADIOS_DATATYPES type
);

void list_insert_read_request_tail (read_request ** h, read_request * q);
void list_insert_read_request_next (read_request ** h, read_request * q);
void list_append_read_request_list (read_request ** h, read_request * q);
void list_free_read_request (read_request * h);
int list_get_length (read_request * h);

void * bufdup(const void *buf, uint64_t elem_size, uint64_t count);

#endif
