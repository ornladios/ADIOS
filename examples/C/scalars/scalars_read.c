/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*********************************************************************/
/*   Example of reading various types of scalar variables in ADIOS   */
/*********************************************************************/
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"

typedef struct complex
{
   float r;
   float i; 
} complex;

typedef struct double_complex
{
   double r;
   double i; 
} double_complex;


int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i;
    int         NX = 10; 
    double      t[NX];
    MPI_Comm    comm = MPI_COMM_WORLD;

    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle, adios_buf_size;

    int8_t  v1 = 0;
    int16_t v2 = 0;
    int32_t v3 = 0;
    int64_t v4 = 0;

    uint8_t  v5 = 0;
    uint16_t v6 = 0;
    uint32_t v7 = 0;
    uint64_t v8 = 0;

    float v9 = 0.0;
    double v10 = 0.0;

    char v11[20];

    complex v12;
    v12.r = 0.0;
    v12.i = 0.0;

    double_complex v13;
    v13.r = 0.0;
    v13.i = 0.0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);

    strcpy (filename, "scalars.bp");

    adios_init ("scalars.xml");
    adios_open (&adios_handle, "scalars", filename, "r", &comm);
#include "gread_scalars.ch"
    adios_close (adios_handle);

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();

    if (rank == 0) {
        printf("byte        v1  = %d\n", v1);
        printf("short       v2  = %d\n", v2);
        printf("integer     v3  = %d\n", v3);
        printf("long        v4  = %lld\n", v4);

        printf("uns.byte    v5  = %u\n", v5);
        printf("uns.short   v6  = %u\n", v6);
        printf("uns.int     v7  = %u\n", v7);
        printf("uns.long    v8  = %llu\n", v8);

        printf("float       v9  = %g\n", v9);
        printf("double      v10 = %g\n", v10);

        printf("string      v11 = %s\n", v11);

        printf("complex     v12 = (%g, i%g)\n", v12.r, v12.i);
        printf("dbl-complex v13 = (%g, i%g)\n", v13.r, v13.i);
    }

    return 0;
}
