/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*************************************************************/
/*          Example of reading arrays in ADIOS               */
/*    which were written from the same number of processors  */
/*                                                           */
/*        Similar example is manual/2_adios_read.c           */
/*************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, j;
    int         NX, NY; 
    double      *t;
    int         *p;
    MPI_Comm    comm = MPI_COMM_WORLD;

    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle, adios_buf_size;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);


    strcpy (filename, "arrays.bp");
    adios_init ("arrays.xml");

    /* First read in the scalars to calculate the size of the arrays */
    adios_open (&adios_handle, "arrays", filename, "r", &comm);
    adios_groupsize = 0;
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_buf_size = 4;
    adios_read (adios_handle, "NX", &NX, adios_buf_size);
    adios_read (adios_handle, "NY", &NY, adios_buf_size);
    adios_close (adios_handle);
    /* Note that we have to close to perform the reading of the variables above */

    printf("rank=%d: NX=%d NY=%d\n", rank, NX, NY);

    /* Allocate space for the arrays */
    t = (double *) malloc (NX*NY*sizeof(double));
    p = (int *) malloc (NX*sizeof(int));

    /* Read the arrays */
    adios_open (&adios_handle, "arrays", filename, "r", &comm);
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_buf_size = 8 * (NX) * (NY);
    adios_read (adios_handle, "var_double_2Darray", t, adios_buf_size);
    adios_buf_size = 4 * (NX);
    adios_read (adios_handle, "var_int_1Darray", p, adios_buf_size);
    adios_close (adios_handle);

    printf("rank=%d: p = [%d", rank, p[0]);
    for (i = 1; i < NX; i++)
        printf(", %d", p[i]);
    printf("]\n");
    
    printf("rank=%d: t[5,*] = [%6.2f", rank, t[5*NY]);
    for (j = 1; j < NY; j++)
        printf(", %6.2f", t[5*NY+j]);
    printf("]\n");

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();

    return 0;
}
