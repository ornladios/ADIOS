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
#include "adios_read.h"

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
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel1=NULL;

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

    char v11[256];

    complex v12;
    v12.r = 0.0;
    v12.i = 0.0;

    double_complex v13;
    v13.r = 0.0;
    v13.i = 0.0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);

    strcpy (filename, "scalars.bp");

    adios_read_init_method (method, comm, "verbose=3");
    ADIOS_FILE * f = adios_read_open (filename, method, comm, ADIOS_LOCKMODE_NONE, 0.0);

    adios_schedule_read (f, sel1, "var_byte",           0, 1, &v1);
    adios_schedule_read (f, sel1, "var_short",          0, 1, &v2);
    adios_schedule_read (f, sel1, "var_int",            0, 1, &v3);
    adios_schedule_read (f, sel1, "var_long",           0, 1, &v4);
    adios_schedule_read (f, sel1, "var_ubyte",          0, 1, &v5);
    adios_schedule_read (f, sel1, "var_ushort",         0, 1, &v6);
    adios_schedule_read (f, sel1, "var_uint",           0, 1, &v7);
    adios_schedule_read (f, sel1, "var_ulong",          0, 1, &v8);
    adios_schedule_read (f, sel1, "var_real",           0, 1, &v9);
    adios_schedule_read (f, sel1, "var_double",         0, 1, &v10);
    /* note that a string is an array and thus v11 a pointer already, 
       so we pass the v11 instead of &v11 here */
    adios_schedule_read (f, sel1, "var_string",         0, 1, v11);
    adios_schedule_read (f, sel1, "var_complex",        0, 1, &v12);
    adios_schedule_read (f, sel1, "var_double_complex", 0, 1, &v13);
    adios_perform_reads (f,1);

    if (rank == 0) {
        printf("byte        v1  = %d\n", v1);
        printf("short       v2  = %d\n", v2);
        printf("integer     v3  = %d\n", v3);
        printf("long        v4  = %" PRId64 "\n", v4);

        printf("uns.byte    v5  = %u\n", v5);
        printf("uns.short   v6  = %u\n", v6);
        printf("uns.int     v7  = %u\n", v7);
        printf("uns.long    v8  = %" PRIu64 "\n", v8);

        printf("float       v9  = %g\n", v9);
        printf("double      v10 = %g\n", v10);

        printf("string      v11 = %s\n", v11);

        printf("complex     v12 = (%g, i%g)\n", v12.r, v12.i);
        printf("dbl-complex v13 = (%g, i%g)\n", v13.r, v13.i);
    }

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    MPI_Finalize ();

    return 0;
}
