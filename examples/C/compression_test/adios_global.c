/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

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

#define NAME_LEN	256

void usage()
{
    printf("Usage: mpirun -np 4 ./adios_global <IO_Compression>\n");
    printf("Example: mpirun -np 4 ./adios_global mpi_zlib\n");
}

int main (int argc, char ** argv)
{
    char        filename [256];
    int         rank, size, i;
    int         NX = 10;
    double      t[NX];
    MPI_Comm    comm = MPI_COMM_WORLD;

    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    if(argc < 2)
    {
        usage();
        return -1;
    }

    char* option = argv[1];
    // char option[NAME_LEN] = {0};
    char bp_file_name[NAME_LEN] = {0};
    char xml_file_name[NAME_LEN] = {0};

    snprintf(bp_file_name, NAME_LEN-1, "output/%s.bp", option);
    snprintf(xml_file_name, NAME_LEN-1, "conf/%s.xml", option);

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    for (i = 0; i < NX; i++)
        t[i] = rank * NX + i;

//	strcpy (filename, "adios_global.bp");
    strcpy (filename, bp_file_name);


    adios_init (xml_file_name, comm);

    adios_open (&adios_handle, "temperature", filename, "w", comm);
//	#include "gwrite_temperature.ch"
    adios_groupsize = 4 \
                        + 4 \
                        + 4 \
                        + 8 * (1) * (NX);

    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "temperature", t);
    adios_close (adios_handle);

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}
