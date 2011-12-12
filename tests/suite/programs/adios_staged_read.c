/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* This code is used to test staged-read method.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mpi.h"
#include "adios_read.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, j, NX;
    MPI_Comm    comm = MPI_COMM_WORLD;
    void * data = NULL, * data1 = NULL, * data2 = NULL;
    uint64_t start[2], count[2], bytes_read = 0, slice_size;

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    /* set the read method to staged read */
    adios_set_read_method (ADIOS_READ_METHOD_BP_STAGED1);

    ADIOS_FILE * f = adios_fopen ("adios_amr_write.bp", comm);

    if (!f)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_GROUP * g = adios_gopen (f, "temperature");
    if (!g)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * v = adios_inq_var (g, "temperature");

    start[0] = 0;
    count[0] = v->dims[0];

    slice_size = v->dims[1]/size;

    start[1] = slice_size * rank;
    if (rank == size - 1)
    {
        slice_size = slice_size + v->dims[1] % size;
    }
    count[1] = slice_size;

    data = malloc (slice_size * v->dims[0] * sizeof (double));
    assert (data);

    bytes_read = adios_read_var (g, "temperature", start, count, data);

    adios_gclose (g);

    adios_fclose (f);

    if (rank == 0)
    {
        for (i = 0; i < v->dims[0]; i++)
        {
            for (j = 0; j < slice_size; j++)
            {
                printf (" %7.5g", * ((double *)data + i * slice_size + j));
            }

            printf ("\n");
        }
    }

    free (data);
    adios_free_varinfo (v);

    MPI_Finalize ();

    return 0;
}
