#ifndef UTIL_MPI_H_
#define UTIL_MPI_H_

#include <stdint.h>
#include "public/adios_mpi.h"

// This helper routine returns a vector of unique NID's
int get_unique_nids (MPI_Comm comm, uint32_t ** nids);

#endif
