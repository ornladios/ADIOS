#include "mpi.h"

int output_init(MPI_Comm comm, int bufsizeMB);

int output_define(int nx, int ny, int gnx, int gny, int offsx, int offsy);

int output_dump(char *filename, int step, void *data);

int output_finalize (int rank);

