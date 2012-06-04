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
#include "public/adios_read.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, j, NX = 16;
    MPI_Comm    comm = MPI_COMM_WORLD;
    void * data = NULL, * data1 = NULL, * data2 = NULL;
    uint64_t start[2], count[2], bytes_read = 0;
    struct timeval t0, t1;

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
#if 0
    adios_set_read_method (ADIOS_READ_METHOD_BP_STAGED1);

    MPI_Barrier (MPI_COMM_WORLD);
    gettimeofday (&t0, NULL);

    MPI_Comm_rank (comm, &rank);

    ADIOS_FILE * f = adios_fopen ("adios_global.bp", comm);
//    ADIOS_FILE * f = adios_fopen ("adios_amr_write.bp", comm);

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

    /* Using less readers to read the global array back, i.e., non-uniform */
    uint64_t slice_size = v->dims[1]/size;
    start[0] = 0;
//printf ("slice_size = %llu\n", slice_size);
    if (rank == size-1)
        slice_size = slice_size + v->dims[1]%size;

    count[0] = v->dims[0];

    start[1] = rank * slice_size;
    count[1] = slice_size;

    data = malloc (slice_size * v->dims[0] * sizeof (double));
    if (data == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    //FIXME: temperature and pressure has the same data
//    bytes_read = adios_read_var (g, "NX", start, count, &NX);
    bytes_read = adios_read_var (g, "temperature", start, count, data);

    adios_gclose (g);

    adios_fclose (f);
//printf ("NX = %d\n", NX);


/*
if (rank == 0)
{
    for (i = 0; i < v->dims[0]; i++) {
//        printf ("rank %3d: temperature [%lld,%d:%lld]", rank, start[0]+i, 0, slice_size);
        for (j = 0; j < slice_size; j++)
            printf (" %7.5g", * ((double *)data + i * slice_size  + j));
        printf ("\n");
    }
}
*/
    MPI_Barrier (MPI_COMM_WORLD);
    gettimeofday (&t1, NULL);

    free (data);

    adios_free_varinfo (v);

    if (rank == 0)
    {
        printf ("IO time = %f \n", t1.tv_sec - t0.tv_sec + (double)(t1.tv_usec - t0.tv_usec)/1000000 );
    }
#endif
    MPI_Finalize ();
    return 0;
}
