/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 *
 * This code is using the generic read API, which can read in
 * arbitrary slices of an array and thus we can read in an array
 * on arbitrary number of processes (provided our code is smart 
 * enough to do the domain decomposition).
 *
 * Run this example after adios_global, which generates 
 * adios_global.bp. Run this example on equal or less 
 * number of processes since we decompose only on one 
 * dimension of the global array here. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"
#include "adios_read.h"

int main (int argc, char ** argv) 
{
    int         rank, size, i, j, npl, token;
    MPI_Comm    comm = MPI_COMM_WORLD;
    MPI_Status  status;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;
    void * data = NULL;
    uint64_t start[1], count[1];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=3");

    ADIOS_FILE * f = adios_read_open ("adios_global_no_xml.bp", method, 
                                      comm, ADIOS_LOCKMODE_NONE, 0);
    if (f == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * v = adios_inq_var (f, "temperature");

    /* Using less readers to read the global array back, i.e., non-uniform */
    uint64_t slice_size = v->dims[0]/size;
    start[0] = slice_size * rank;
    if (rank == size-1) /* last rank may read more lines */
        slice_size = slice_size + v->dims[0]%size;
    count[0] = slice_size;

    data = malloc (slice_size * sizeof (double));
    if (data == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    /* Read a subset of the temperature array */
    sel = adios_selection_boundingbox (v->ndim, start, count);
    adios_schedule_read (f, sel, "temperature", 0, 1, data);
    adios_perform_reads (f, 1);

    if (rank > 0) {
        MPI_Recv (&token, 1, MPI_INT, rank-1, 0, comm, &status);
    }

    printf (" ======== Rank %d ========== \n", rank);
    npl = 10;
    for (i = 0; i < slice_size; i+=npl) {
        printf ("[%4.4" PRIu64 "]  ", rank*slice_size+i);
        for (j= 0; j < npl; j++) {
            printf (" %6.6g", * ((double *)data + i + j));
        }
        printf ("\n");
    }
    fflush(stdout);
    sleep(1);

    if (rank < size-1) {
        MPI_Send (&token, 1, MPI_INT, rank+1, 0, comm);
    }

    free (data);

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);
    MPI_Finalize ();
    return 0;
}
