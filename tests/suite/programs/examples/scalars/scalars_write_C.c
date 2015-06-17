/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*************************************************************/
/*   Example of writing various types of variable in ADIOS   */
/*************************************************************/
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
    int         rank;
    MPI_Comm    comm = MPI_COMM_WORLD;

    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    int8_t v1 = -4;
    int16_t v2 = -3;
    int32_t v3 = -2;
    int64_t v4 = -1;

    uint8_t v5 = 1;
    uint16_t v6 = 2;
    uint32_t v7 = 3;
    uint64_t v8 = 4;

    float v9 = 5.0;
    double v10 = 6.0;

    char * v11 = "ADIOS example";

    complex v12;
    v12.r = 8.0;
    v12.i = 9.0;

    double_complex v13;
    v13.r = 10.0;
    v13.i = 11.0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);

    strcpy (filename, "scalars_C.bp");

    /* adios_open() opens a "group in a file", here the "scalars" group.   
       GWRITE is the convenient way to write all variables defined in the
       xml file but of course one can write the individual adios_write() 
       statements here too 
    */
    adios_init ("scalars_C.xml", comm);
    adios_open (&adios_handle, "scalars", filename, "w", comm);
#include "gwrite_scalars.ch"
    adios_close (adios_handle);

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}
