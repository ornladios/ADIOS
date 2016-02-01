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
    int         rank, size, i, j, k, token;
    MPI_Comm    comm = MPI_COMM_WORLD;
    MPI_Status  status;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;
    void * data = NULL;
    uint64_t start[3], count[3], step = 0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=3");

    /* adios_read_open_file() allows for seeing all timesteps in the file */
    ADIOS_FILE * f = adios_read_open_file ("adios_globaltime.bp", method, comm);
    if (f == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * v = adios_inq_var (f, "temperature");

    // read in two timesteps
    data = malloc (2 * v->dims[0] * v->dims[1] * sizeof (double));
    if (data == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    // read in timestep 'rank' (up to 12)
    step = rank % 13;

    start[0] = 0;
    count[0] = v->dims[0];

    start[1] = 0;
    count[1] = v->dims[1];

    /* Read a subset of the temperature array */
    sel = adios_selection_boundingbox (v->ndim, start, count);
    /*    2 steps from 'step' */
    adios_schedule_read (f, sel, "temperature", step, 2, data);
    adios_perform_reads (f, 1);

    if (rank == 0) 
        printf ("Array size of temperature [0:%" PRIu64 ",0:%" PRIu64 "]\n", v->dims[0], v->dims[1]);

    if (rank > 0) {
        MPI_Recv (&token, 1, MPI_INT, rank-1, 0, comm, &status);
    }

    printf("------------------------------------------------\n");
    printf("rank=%d: \n", rank);
    for (i = 0; i < 2; i++) {
        printf ("step %" PRIu64 " = [\n", step+i);
        for (j = 0; j < v->dims[0]; j++) {
            printf (" [");
            for (k = 0; k < v->dims[1]; k++) {
                printf ("%g ", ((double *)data) [ i * v->dims[0] * v->dims[1] + j * v->dims[1] + k]);
            }
            printf ("]\n");
        }
        printf ("]\n");
    }
    printf ("\n");

    if (rank < size-1) {
        MPI_Send (&token, 1, MPI_INT, rank+1, 0, comm);
    }

    free (data);
    adios_free_varinfo (v);

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);
    MPI_Finalize ();
    return 0;
}


