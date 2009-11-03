/* ADIOS C Example: write a global array from N processors with gwrite
 *
 * How to run: mpirun -np <N> adios_global
 * Output: adios_global.bp
 * ADIOS config file: adios_global.xml
 *
*/
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i;
    int         NX = 10, level = 1; 
    double      t[NX], mean = 0;
    MPI_Comm    comm = MPI_COMM_WORLD;
    const char * str = "Nov, 2009";

    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    for (i = 0; i < NX; i++)
    {
        t[i] = rank * NX + i;
        mean += t[i];
    }

    mean /= NX;

    strcpy (filename, "restart.bp");

    adios_init ("adios_attributes.xml");

    adios_open (&adios_handle, "temperature", filename, "w");

    adios_groupsize = 4 \
                    + 4 \
                    + 4 \
                    + 8 \
                    + 8 * (1) * (NX) \
                    + strlen (str) + 1;

    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize, &comm);

    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "mean", &mean);
    adios_write (adios_handle, "temperature", t);
    adios_write (adios_handle, "date", str);
     
    adios_close (adios_handle);

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}
