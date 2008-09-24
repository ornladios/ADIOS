#include <stdio.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv) 
{
char        filename [256];
int         rank;
int         NX = 10; 
double      t[NX];

/* ADIOS variables declarations for matching gwrite_temperature.ch */
int         adios_err;
uint64_t    adios_groupsize, adios_totalsize;
int64_t     adios_handle;
MPI_Comm comm;
 
int         color, key;
int         size;
MPI_Init (&argc, &argv);
MPI_Comm_rank (MPI_COMM_WORLD, &rank);

/* MPI_Comm_split paritions the world group into disjointed 2 subgroups, 
  * the processes are ranked in terms of the argument key  
  *  a new communicator comm is returned for this specific grid configuration
  */
color = rank % 2;
key = rank / 2;
MPI_Comm_split (MPI_COMM_WORLD, color, key, &comm);
MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
    adios_init ("config_global.xml");
    sprintf (filename, "restart_%5.5d.bp", color);
    for (int it =0; it < 5; it++) {
        for (int i = 0; i < NX; i++)
            t [i] = 10 *(key*2+color)+i;
        /* every P/2 processes write into the same file 
         * there are 2 files generated. 
         */
        adios_open (&adios_handle, "temperature", filename, "a");
        #include "gwrite_temperature.ch"
        adios_close (adios_handle);
    }
    adios_finalize (rank);
    MPI_Finalize ();
return 0;
}
