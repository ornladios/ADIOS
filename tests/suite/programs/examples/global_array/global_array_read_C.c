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
 * global_array_C.bp. Run this example on equal or less 
 * number of processes since we decompose only on one 
 * dimension of the global array here. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"
#include "core/adios_logger.h"

int main (int argc, char ** argv) 
{
    int         rank, size, i, j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;
    void * data = NULL;
    uint64_t start[2], count[2];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=4");
    adios_logger_open ("log_read_C", rank);

    ADIOS_FILE * f = adios_read_open ("global_array_C.bp", method, comm, ADIOS_LOCKMODE_NONE, 0);
    if (f == NULL)
    {
        log_error ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * v = adios_inq_var (f, "temperature");

    /* Using less readers to read the global array back, i.e., non-uniform */
    uint64_t slice_size = v->dims[0]/size;
    start[0] = slice_size * rank;
    if (rank == size-1) /* last rank may read more lines */
        slice_size = slice_size + v->dims[0]%size;
    count[0] = slice_size;

    start[1] = 0;
    count[1] = v->dims[1];
       

    data = malloc (slice_size * v->dims[1] * sizeof (double));
    if (data == NULL)
    {
        log_error (stderr, "malloc failed.\n");
        return -1;
    }

    /* Read a subset of the temperature array */
    sel = adios_selection_boundingbox (v->ndim, start, count);
    adios_schedule_read (f, sel, "temperature", 0, 1, data);
    adios_perform_reads (f, 1);

    for (i = 0; i < slice_size; i++) {
        log_test ("rank %d: [%lld,%d:%lld]", rank, start[0]+i, 0, slice_size);
        for (j = 0; j < v->dims[1]; j++)
            log_test (" %6.6g", * ((double *)data + i * v->dims[1] + j));
        log_test ("\n");
    }

    free (data);

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);
    adios_logger_close();
    MPI_Finalize ();
    return 0;
}
