#ifndef __DECOMPOSE_H_
#define __DECOMPOSE_H_

#include <stdint.h>

/* Decompose a variable among several processors 

   count/start: array of size/offset in each dimension 
   writesize:   sum (count)

*/
void decompose (int numproc, int rank, int ndim, uint64_t *dims, int *decomp_values,
                /*OUT*/ uint64_t *count,
                /*OUT*/ uint64_t *start,
                /*OUT*/ uint64_t *writesize);

#endif
