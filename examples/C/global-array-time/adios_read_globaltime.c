/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, j, k;
    MPI_Comm    comm = MPI_COMM_WORLD;
    void * data = NULL;
    uint64_t start[3], count[3], bytes_read = 0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    ADIOS_FILE * f = adios_fopen ("adios_globaltime.bp", comm);
    if (f == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_GROUP * g = adios_gopen (f, "restart");
    if (g == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * v = adios_inq_var (g, "temperature");

    // read in two timesteps
    data = malloc (2 * v->dims[1] * v->dims[2] * sizeof (double));
    if (data == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    // read in timestep 'rank'
    start[0] = rank % 13;
    count[0] = 1;

    start[1] = 0;
    count[1] = v->dims[1];

    start[2] = 0;
    count[2] = v->dims[2];
       
    bytes_read = adios_read_var (g, "temperature", start, count, data);

    printf("rank=%d: ", rank);
    for (i = 0; i < 1; i++) {
        printf ("[%lld,0:%lld,0:%lld] = [", start[0]+i, v->dims[1], v->dims[2]);   
        for (j = 0; j < v->dims[1]; j++) {
            printf (" [");
            for (k = 0; k < v->dims[2]; k++) {
                printf ("%g ", ((double *)data) [ i * v->dims[1] * v->dims[2] + j * v->dims[2] + k]);
            }
            printf ("]");
        }
        printf (" ]\t");
    }
    printf ("\n");

    free (data);
    adios_free_varinfo (v);

    adios_gclose (g);
    adios_fclose (f);

    MPI_Barrier (comm);

    MPI_Finalize ();
    return 0;
}


