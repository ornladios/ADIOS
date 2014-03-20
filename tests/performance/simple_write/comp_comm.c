#include <stdio.h>
#include <stdlib.h>
#include "comp_comm.h"

double *data;
// communication array
#define NELEMS 8192    // 32 KB of integers
static int    sendbuf[NELEMS];     // array to do communication with it
static int    recvbuf[NELEMS];     // array to do communication with it

void comp_comm_init(MPI_Comm comm)
{
    int rank, i;
    MPI_Comm_rank (comm, &rank);
    for (i = 0; i < NELEMS; i++)
        sendbuf[i] = rank;
}

void do_calc_unit(double * data, int nx, int ny)
{
    int i,j;
    double tmp;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny/2; j++) {
            tmp = data[i*ny + j];
            data[i*ny + j] = data[i*ny + (ny-j-1)];
            data[i*ny + (ny-j-1)] = tmp;
        }
    }
}


void do_comm_unit(MPI_Comm comm)
{
    // rank 0: reduce data from everybody
    MPI_Allreduce (sendbuf, recvbuf, NELEMS, MPI_INTEGER, MPI_MAX, comm);
    //MPI_Bcast (sendbuf, NELEMS, MPI_INTEGER, 0, comm);
}

