/*
 * ADIOS example from the User's manual
 *
 * Read back data written by 2_adios_posix
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
    char        result[1024], s[32];
    int         i;
    
    /* ADIOS variables declarations for matching gread_temperature.ch */
    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize, adios_buf_size;
    int64_t     adios_handle;
    MPI_Comm    comm =  MPI_COMM_SELF;
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    sprintf (filename, "restart_%5.5d.bp", rank);
    adios_init ("config.xml");
    adios_open (&adios_handle, "temperature", filename, "r");
    #include "gread_temperature.ch"
    adios_close (adios_handle);
    adios_finalize (rank);
    MPI_Finalize ();

    sprintf(result, "rank=%d t=[%g", rank, t[0]);
    for (i=1; i<NX; i++) {
        sprintf (s, ",%g", t[i]);
        strcat (result, s);
    }
    printf("%s]\n", result);

    return 0;
}    

