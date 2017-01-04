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
#include "adios.h"

int main (int argc, char ** argv) 
{
    int         rank, size, i, j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;
    void * mesh = NULL;
    double * data = NULL, * grad = NULL, * R = NULL, * Z = NULL;
    uint64_t start[2], count[2];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=3");

    ADIOS_FILE * f = adios_read_open_file ("xgc.3d.08000.bp", method, comm);
    if (f == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * v = adios_inq_var (f, "dpot");

//    printf ("%lu, %lu\n", v->dims[0], v->dims[1]);

    start[0] = 0;
    count[0] = v->dims[0];

    start[1] = rank;
    count[1] = 1;
       

    data = malloc (count[0] * count[1] * sizeof (double));
    grad = malloc (count[0] * count[1] * sizeof (double));
    R = malloc (count[0] * count[1] * sizeof (double));
    Z = malloc (count[0] * count[1] * sizeof (double));
    if (data == NULL || grad == NULL || R == NULL || Z == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    /* Read a subset of the temperature array */
    sel = adios_selection_boundingbox (v->ndim, start, count);
    adios_schedule_read (f, sel, "dpot", 0, 1, data);
    adios_perform_reads (f, 1);
/*
    for (i = 0; i < count[0] * count[1]; i++) {
            printf (" %e", * ((double *)data + i ));
    }
*/
//    printf ("\n");

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);

    // read mesh
    adios_read_init_method (method, comm, "verbose=3");

    ADIOS_FILE * fmesh = adios_read_open_file ("xgc.mesh.bp", method, comm);
    if (fmesh == NULL)
    {
        printf ("%s\n", adios_errmsg());
        return -1;
    }

    ADIOS_VARINFO * conn = adios_inq_var (fmesh, "/cell_set[0]/node_connect_list");

//    printf ("conn: %lu, %lu\n", conn->dims[0], conn->dims[1]);


    start[0] = 0;
    count[0] = conn->dims[0];

    start[1] = 0;
    count[1] = conn->dims[1];


    mesh = malloc (count[0] * count[1] * sizeof (int));
    if (mesh == NULL)
    {
        fprintf (stderr, "malloc failed.\n");
        return -1;
    }

    /* Read a subset of the temperature array */
    sel = adios_selection_boundingbox (conn->ndim, start, count);
    adios_schedule_read (fmesh, sel, "/cell_set[0]/node_connect_list", 0, 1, mesh);
    start[0] = 0;
    start[1] = 0;

    count[0] = v->dims[0];
    count[1] = 1;

    sel = adios_selection_boundingbox (v->ndim, start, count);
    adios_schedule_read (fmesh, sel, "/coordinates/values", 0, 1, R);


    start[0] = 0;
    start[1] = 1;

    count[0] = v->dims[0];
    count[1] = 1;

    sel = adios_selection_boundingbox (v->ndim, start, count);
    adios_schedule_read (fmesh, sel, "/coordinates/values", 0, 1, Z);

    adios_perform_reads (fmesh, 1);

/*
    for (i = 0; i < count[0]; i++) {
        int n1 = * ((int *) mesh + i * 3);
        int n2 = * ((int *) mesh + i * 3 + 1);
        int n3 = * ((int *) mesh + i * 3 + 2);

        grad[n1] = grad[n2] = grad[n3] = 
              (double)(abs (data[n1] - data[n2]) + abs(data[n1] - data[n3]) + abs(data[n2] - data[n3])) / 3;
    }
*/
    adios_read_close (fmesh);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);


    int NX = v->dims[0], GX = v->dims[0], OX = 0;
    int MY = conn->dims[0], MX = conn->dims[1];

    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    adios_init ("test_xgc.xml", comm);

    adios_open (&adios_handle, "field", "dpot.bp", "w", comm);
    adios_groupsize = 6 * 4 + 3 * 8 * NX + 4 * MX * MY;
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);

    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "GX", &GX);
    adios_write (adios_handle, "OX", &OX);
    adios_write (adios_handle, "MX", &MX);
    adios_write (adios_handle, "MY", &MY);
    adios_write (adios_handle, "mesh", mesh);
    adios_write (adios_handle, "R", R);
    adios_write (adios_handle, "Z", Z);
    adios_write (adios_handle, "dpot", data);

    adios_close (adios_handle);

    free (data);
    free (grad);
    free (mesh);

    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}
