#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include "decompose.h"
#include "utils.h"

int ndim, nproc;
uint64_t *dims;
int *np;
int rank; // unused, required by print0 in utils.h

void printUsage(char *prgname)
{
    print0("Usage: %s N M dim1 ... dimN np1 ... npN\n"
           "    N           Number of dimensions\n"
           "    M           Number of processes\n"
           "    dim1..dimN  Dimensions\n"
           "    np1..npN    list of decomposition numbers e.g. 32 8 4\n"
           "            Decomposition values in each dimension of an array\n"
           "            The product of these number must be less then the number\n"
           "            of processes. Processes whose rank is higher than the\n"
           "            product, will not write anything.\n"
           "               Arrays with less dimensions than the number of values,\n"
           "            will be decomposed with using the appropriate number of\n"
           "            values.\n"
        ,prgname);
}

int processArgs(int argc, char ** argv)
{   
    int i, j, nd, prod; 
    char *end;

    if (argc < 3) {
        printUsage (argv[0]);
        return 1;
    }

    errno = 0;
    ndim = strtol(argv[1], &end, 10);
    if (errno || (end != 0 && *end != '\0')) {
        print0 ("ERROR: Invalid decomposition number in argument 1: '%s'\n", argv[j]);
        printUsage(argv[0]);
        return 1;
    }

    errno = 0;
    nproc = strtol(argv[2], &end, 10);
    if (errno || (end != 0 && *end != '\0')) {
        print0 ("ERROR: Invalid decomposition number in argument 2: '%s'\n", argv[j]);
        printUsage(argv[0]);
        return 1;
    }

    if (argc != 3+2*ndim) {
        print0 ("ERROR: Expected number of arguments is %d\n", 2*ndim+2);
        printUsage (argv[0]);
        return 1;
    }
    
    dims = (uint64_t*) malloc (sizeof(uint64_t) * ndim);
    np   = (int*) malloc (sizeof(int) * ndim);

    nd = 0;
    j = 3;
    for (nd=0; nd < ndim; nd++) { // get ndim dimensions
        errno = 0;
        dims[nd] = strtoll(argv[j], &end, 10);
        if (errno || (end != 0 && *end != '\0')) {
            print0 ("ERROR: Invalid decomposition number in argument %d: '%s'\n",
                    j, argv[j]);
            printUsage(argv[0]);
            return 1;
        }
        j++;
    }

    for (nd=0; nd < ndim; nd++) { // get ndim decomposition values
        errno = 0;
        np[nd] = strtoll(argv[j], &end, 10);
        if (errno || (end != 0 && *end != '\0')) {
            print0 ("ERROR: Invalid decomposition number in argument %d: '%s'\n",
                    j, argv[j]);
            printUsage(argv[0]);
            return 1;
        }
        j++;
    }

    prod = 1;
    for (i=0; i<ndim; i++) {
        prod *= np[i];
    }

    if (prod > nproc) {
        print0 ("ERROR: Product of decomposition numbers %d > number of processes %d\n",
                prod, nproc);
        printUsage(argv[0]);
        return 1;
    }

    return 0;
}

int main (int argc, char ** argv)
{
    char ints[256];
    uint64_t count[10], start[10], totalsize;
    int i, nd;

    if (processArgs(argc, argv)) 
        return 1;

    print("# of proc: %d\n", nproc);
    ints_to_str(ndim, np, ints);
    print("decomposition: %s\n", ints);

    for (nd = ndim; nd > 0; nd--)
    {
        int64s_to_str(nd, dims, ints);
        print("--------------------------------------------\n");
        print("%d  dimensions: %s\n", nd, ints);
        for (i=0; i<nproc; i++) {
            decompose (nproc, i, nd, dims, np, count, start, &totalsize);
        }
    }

    return 0;
}


