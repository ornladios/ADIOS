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


/*************************************************************/
/*   Example of writing various types of variable in ADIOS   */
/*************************************************************/
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i;
    int         NX = 10; 
    double      t[NX];
    MPI_Comm    comm = MPI_COMM_WORLD;

    int         adios_err;
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

    // the size of 'long double' is 8-bytes using PGI compiler
    //long double v11 = 7.0;

    const char * v12 = "ADIOS example";

    complex v13;
    v13.r = 8.0;
    v13.i = 9.0;

    double_complex v14;
    v14.r = 10.0;
    v14.i = 11.0;

    MPI_Init (&argc, &argv);

    strcpy (filename, "restart.bp");
    adios_init ("config.xml");
    adios_open (&adios_handle, "my_group", filename, "w");

    adios_groupsize = sizeof (int8_t)         \
                    + sizeof (int16_t)        \
                    + sizeof (int32_t)        \
                    + sizeof (int64_t)        \
                    + sizeof (uint8_t)        \
                    + sizeof (uint16_t)       \
                    + sizeof (uint32_t)       \
                    + sizeof (uint64_t)       \
                    + sizeof (float)          \
                    + sizeof (double)         \
                    + sizeof (long double)    \
                    + strlen (v12) + 1        \
                    + sizeof (complex)        \
                    + sizeof (double_complex);
 
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize, &comm);

    adios_write (adios_handle, "var_byte", &v1);
    adios_write (adios_handle, "var_short", &v2);
    adios_write (adios_handle, "var_int", &v3);
    adios_write (adios_handle, "var_long", &v4);

    adios_write (adios_handle, "var_ubyte", &v5);
    adios_write (adios_handle, "var_ushort", &v6);
    adios_write (adios_handle, "var_uint", &v7);
    adios_write (adios_handle, "var_ulong", &v8);

    adios_write (adios_handle, "var_real", &v9);
    adios_write (adios_handle, "var_double", &v10);
//  adios_write (adios_handle, "var_long_double", &v11);

    adios_write (adios_handle, "var_string", v12);
    adios_write (adios_handle, "var_complex", &v13);
    adios_write (adios_handle, "var_double_complex", &v14);

    adios_close (adios_handle);

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}
