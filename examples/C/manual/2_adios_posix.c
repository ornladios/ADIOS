/*
 * ADIOS example from the User's manual
 *
 * Write a separate file from each process by using 
 * the ADIOS POSIX method
 *
 */
#include <stdio.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank;
    int         NX = 10;
    double      t[NX];
    int         i;
    
    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;
    MPI_Comm  * comm =  MPI_COMM_SELF;
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    for (i=0; i<NX; i++)
        t[i] = rank*NX + i;

    sprintf (filename, "restart_%5.5d.bp", rank);
    adios_init ("config.xml");
    adios_open (&adios_handle, "temperature", filename, "w", &comm);
    #include "gwrite_temperature.ch"
    adios_close (adios_handle);
    adios_finalize (rank);
    MPI_Finalize ();
    return 0;
}

