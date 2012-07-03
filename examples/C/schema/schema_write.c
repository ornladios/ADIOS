/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write some attributes along with variables
*/
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"
int main (int argc, char ** argv[] ) 
{
    if (argc != 2) /* argc should be 2 for correct execution */
    {
        /* We print argv[0] assuming it is the program name */
        printf("usage: %s meshname\n", argv[0]);
    }  
    else
    {
        char        filename [256];
        char        meshname [256];
        char        xmlfilename[256];
        int         rank, size, i;
        int         NX = 10, level = 1; 
        double      t[NX], mean = 0;
        MPI_Comm    comm = MPI_COMM_WORLD;
        const char * str = "Jul, 2012";

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

        strcpy (meshname, argv[1]);

        if (!strcmp(meshname,"uniform") || !strcmp(meshname,"structured") || !strcmp(meshname,"rectilinear") || !strcmp(meshname,"unstructured"))
        {
            //strcpy (meshname, "uniform");
            //strcpy (filename, "schema");
            strcpy (filename, meshname);
            strcat (filename, ".bp");

            strcpy (xmlfilename,meshname);
            strcat (xmlfilename,".xml");
            adios_init (xmlfilename);

            adios_open (&adios_handle, "schema", filename, "w", &comm);
#include "gwrite_temperature.ch"
            adios_close (adios_handle);

            MPI_Barrier (comm);

            adios_finalize (rank);

            MPI_Finalize ();
        }
        else
        {
            /* We print that the mesh name is invalid */
            printf("meshname: %s is invalid. Please choose from the following options: uniform, structured, rectilinear or unstructured.\n", meshname);
        }
    }
    return 0;
}
