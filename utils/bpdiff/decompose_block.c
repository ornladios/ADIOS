/* 
Decompose arrays in all dimensions (block), evenly among 
the processes. 

Does not care how small the blocks become, so small arrays
with many writers may suffer with this approach.

*/

#include "decompose.h"
#include "utils.h"

void decompose (int numproc, int rank, int ndim, uint64_t *dims, 
                int *np, // number of processes in each dimension
                /*OUT*/ uint64_t *count,
                /*OUT*/ uint64_t *start,
                /*OUT*/ uint64_t *writesize)
{
    int i;
    int pos[10]; // rank's position in each dimensions

    if (ndim == 0) {
        // scalars -> rank 0 writes them
        if (rank == 0) 
            *writesize = 1;
        else 
            *writesize = 0;
        return;
    }

    /* calculate this process' position in the n-dim space
        0 1 2
        3 4 5
        6 7 8

        for 1D: 
        posx = rank/1             ! 1st dim: 0, 1, 2...,rank-1 are in the same X position
    
        for 2D: 
        posx = mod(rank, npx)     ! 1st dim: 0, npx, 2npx... are in the same X position
        posy = rank/(npx)         ! 2nd dim: npx processes belong into one dim

        for 3D: 
        posx = mod(rank, npx)     ! 1st dim: 0, npx, 2npx... are in the same X position
        posy = mod(rank/npx, npy) ! 2nd dim: (0, npx-1) have the same dim (so divide with npx first)
        posz = rank/(npx*npy)     ! 3rd dim: npx*npy processes belong into one dim
    */
    int nps = 1;
    for (i=0; i<ndim-1; i++)
    {
        pos[i] = (rank / nps) % np[i];
        nps *= np[i];
    }
    pos[i] = rank / nps;

    char ints[256];
    ints_to_str(ndim, pos, ints);
    if (pos[ndim-1] >= np[ndim-1]) {
        print("rank %d: position in %d-D space = %s ---> Out of bound process\n", rank, ndim, ints);
    } else {
        //print("rank %d: position in %d-D space = %s\n", rank, ndim, ints);
    }

    /* Decompose each dimension according to the position */
    *writesize = 1;
    for (i=0; i<ndim; i++)
    {
        if (pos[ndim-1] >= np[ndim-1]) 
        {
            // this process gets nothing to read
            count[i] = 0;
            start[i] = 0;
        } 
        else 
        {
            count[i] = dims[i] / np[i];
            start[i] = count[i] * pos[i];
            if (pos[i] == np[i]-1) {
                // last one in the dimension may need to read more than the rest
                count[i] = dims[i] - count[i]*(np[i]-1);
            }
        }
        *writesize *= count[i];
    }
    int64s_to_str(ndim, count, ints);
    //print("rank %d: ldims   in %d-D space = %s\n", rank, ndim, ints);
    int64s_to_str(ndim, start, ints);
    //print("rank %d: offsets in %d-D space = %s\n", rank, ndim, ints);
}


