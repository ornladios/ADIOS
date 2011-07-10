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
    uint64_t start[2], count[2], bytes_read = 0;
    struct timeval t0, t1;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_set_read_method (ADIOS_READ_METHOD_BP_SUBFILE);

    MPI_Barrier (MPI_COMM_WORLD);
    gettimeofday (&t0, NULL);

    MPI_Comm_rank (comm, &rank);

    ADIOS_FILE * f = adios_fopen ("adios_global.bp", comm);

    if (f == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_GROUP * g = adios_gopen (f, "temperature");
    if (g == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * v = adios_inq_var (g, "temperature");
v->dims[1] = 5;
    /* Using less readers to read the global array back, i.e., non-uniform */
    uint64_t slice_size = v->dims[0]/size;
//slice_size = 15;
    start[0] = slice_size * rank;
//start[0] = 1;

    if (rank == size-1)
        slice_size = slice_size + v->dims[0]%size;

    count[0] = slice_size;

    start[1] = 1;
    count[1] = v->dims[1];

    data = malloc (slice_size * v->dims[1] * sizeof (double));
    if (data == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    data1 = malloc (slice_size * v->dims[1] * sizeof (double));
    if (data1 == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    data2 = malloc (slice_size * v->dims[1] * sizeof (double));
    if (data2 == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    //FIXME: temperature and pressure has the same data
    bytes_read = adios_read_var (g, "NX", start, count, &NX);
    bytes_read = adios_read_var (g, "temperature", start, count, data);
    bytes_read = adios_read_var (g, "pressure", start, count, data1);
    bytes_read = adios_read_var (g, "speed", start, count, data2);

    adios_gclose (g);

    adios_fclose (f);

    for (i = 0; i < slice_size; i++) {
        printf ("rank %3d: temperature [%lld,%d:%lld]", rank, start[0]+i, 0, slice_size);
        for (j = 0; j < v->dims[1]; j++)
            printf (" %6.4g", * ((double *)data + i * v->dims[1] + j));
        printf ("\n");
    }

    for (i = 0; i < slice_size; i++) {
        printf ("rank %3d: pressure    [%lld,%d:%lld]", rank, start[0]+i, 0, slice_size);
        for (j = 0; j < v->dims[1]; j++)
            printf (" %6.4g", * ((double *)data1 + i * v->dims[1] + j));
        printf ("\n");
    }

    for (i = 0; i < slice_size; i++) {
        printf ("rank %3d: speed       [%lld,%d:%lld]", rank, start[0]+i, 0, slice_size);
        for (j = 0; j < v->dims[1]; j++)
            printf (" %6.4g", * ((double *)data2 + i * v->dims[1] + j));
        printf ("\n");
    }

    MPI_Barrier (MPI_COMM_WORLD);
    gettimeofday (&t1, NULL);

    free (data);
    free (data1);
    free (data2);

    adios_free_varinfo (v);

    if (rank == 0)
    {
        printf ("IO time = %f \n", t1.tv_sec - t0.tv_sec + (double)(t1.tv_usec - t0.tv_usec)/1000000 );
    }
    MPI_Finalize ();
    return 0;
}
