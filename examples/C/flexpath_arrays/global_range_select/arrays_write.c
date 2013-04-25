/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"

/*************************************************************/
/*          Example of writing arrays in ADIOS               */
/*                                                           */
/*        Similar example is manual/2_adios_write.c          */
/*************************************************************/
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i;
    int         NX = 10; 
    int         NY = 1;
    double      t[NX];
    MPI_Comm    comm = MPI_COMM_WORLD;

    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
    
    fprintf(stderr, "starting writer.\n");
    strcpy (filename, "arrays");
    adios_init ("arrays.xml");
    
    
    int ii;
    for(ii = 0; ii<30; ii++){
      for (i = 0; i < NX; i++)
        t[i] = rank * NX + i*ii;
      fprintf(stderr, "open\n");
      adios_open (&adios_handle, "temperature", filename, "w", &comm);
      fprintf(stderr, "scalar write\n");
      adios_write (adios_handle, "NX", &NX);
      adios_write (adios_handle, "NY", &NY);
      fprintf(stderr, "array writes\n");
      adios_write (adios_handle, "size", &size);
      adios_write (adios_handle, "rank", &rank);
      adios_write (adios_handle, "var_2d_array", t);
      fprintf(stderr, "in app: rank: %d, size: %d, NX: %d\n", rank, size, NX);
      adios_close (adios_handle);
      fprintf(stderr, "commited write %d\n", ii);
    }
    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}
