#include "decompose.h"

/* Silly decompose: let process 0 do everything */
void decompose (int numproc, int rank, int ndim, uint64_t *dims, int *decomp_values,
                /*OUT*/ uint64_t *count,
                /*OUT*/ uint64_t *start,
                /*OUT*/ uint64_t *writesize)
{
    int i;
    if (rank == 0) 
        *writesize = 1;
    else 
        *writesize = 0;

    for (i=0; i<ndim; i++)
    {
        if (rank == 0) {
            count[i] = dims[i];
            *writesize *= dims[i];
        } else {
            count[i] = 0;
        }
        start[i] = 0;
    }
}


