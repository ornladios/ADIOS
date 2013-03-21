/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write a variable along with a structured mesh. 
 * Note that the mesh dimensions depend on the rank.
*/
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char * argv[] ) 
{
    char        filename[]="structured.bp";
    int         rank, size, i, j;
    int         NX = 10;  
    double      t[NX], mean = 0;
    MPI_Comm    comm = MPI_COMM_WORLD;
    char      * str = "Jul, 2012";
    int         nspace=2;

    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    float       X[size*NX];
    float       Y[size][NX];
    float       XY[2*size*NX];
    int         size2 = 2*size*NX;
    int idx = 0;
    for (i = 0; i< size; i++)
    {
        for (j=0; j<NX; j++)
        {
            X[idx] = i + j;
            Y[i][j] = i - j;
            XY[2*idx] = X[idx];
            XY[2*idx+1]=Y[i][j];
            idx++;
        }
    }
   
    for (i = 0; i< NX; i++) {
        t[i] = i+rank;
    }
    mean /= NX;

    adios_init ("structured.xml", comm);

    adios_open (&adios_handle, "schema", filename, "w", &comm);

    adios_groupsize = 4 \
                      + 4 \
                + 4 \
                + 4 \
                + 4 \
                + 8 \
                + strlen(str) \
                + sizeof(float) * (size) * (NX) \
                + sizeof(float) * (size) * (NX) \
                + sizeof(float) * (2*NX*size) \
                + sizeof(double) * (NX);
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "size2", &size2);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "mean", &mean);
    adios_write (adios_handle, "nspace", &nspace);
    adios_write (adios_handle, "X", X);
    adios_write (adios_handle, "Y", Y);
    adios_write (adios_handle, "XY",  XY);
    adios_write (adios_handle, "date", str);
    adios_write (adios_handle, "temperature", t);

    adios_close (adios_handle);

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();

    return 0;
}
