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
    /* application data structures */
    char        filename [256];
    int         rank, i, j;
    int         NX = 10, NY = 100; 
    double      t[NX][NY];
    int         p[NX];

    /* MPI and ADIOS data structures */
    MPI_Comm    comm = MPI_COMM_WORLD;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    /* MPI and ADIOS setup */
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    strcpy (filename, "arrays");
    adios_init ("arrays.xml", comm);
    
    /* write data 200 times so that it forms a stream */
    int ii=0;
    for(ii=0; ii<30; ii++) {
        
        /* must open eveery time for stream api to work */
        adios_open (&adios_handle, "arrays", filename, "w", comm);
        
        /* initialize 2d array to:
            (rank*1000000) + (stream step * 1000) + offset */
        for (i = 0; i < NX; i++)
            for (j = 0; j< NY; j++)
                t[i][j] = (rank * 1000000) + (ii * 1000) + (i*NY) + j;

        /* initialize 1d array similarly */
        for (i = 0; i < NX; i++)
            p[i] = (rank*1000000) + (ii * 1000) + i;
        
        /* groupsize registration (necessary?) */
        adios_groupsize = 4 + 4 + 8 * (NX) * (NY) + 4 * (NX);
        adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);

        /* write a step */
        adios_write (adios_handle, "NX", &NX);
        adios_write (adios_handle, "NY", &NY);
        adios_write (adios_handle, "var_double_2Darray", t);
        adios_write (adios_handle, "var_int_1Darray", p);
        //fprintf(stderr, "2d arr %p 1d arr %p\n", t, p);
   
        /* commit the write */
        adios_close (adios_handle);
        printf("Committed Rank=%d Step=%d\n\n", rank, ii);
    }

    /* insure all writers are done before shutdown */
    MPI_Barrier (comm);

    /* shutdown adios and mpi */
    adios_finalize (rank);
    MPI_Finalize ();

    return 0;
}
