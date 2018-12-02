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
#include <assert.h>

#include "mpi.h"
#include "adios_read.h"
#include "adios.h"
#include "mgard_capi.h"


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

    if (argc != 2) printf ("Wrong arguments\n");
    double tolerance = atof(argv[1]);

    printf ("tolerance = %f\n", tolerance);

    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    int step;
    double iotime = 0.0;

    for (step = 0; step < 1; step++)
    {
        char fname[64];

        //sprintf (fname, "xgc.3d.080%02d.bp", step);
        sprintf (fname, "larger_data/xgc.3d.00480.bp");
        //sprintf (fname, "xgc-f0/xgc.dataspaces.f0analysis.00005.bp");
        adios_read_init_method (method, comm, "verbose=3");

        ADIOS_FILE * f = adios_read_open_file (fname, method, comm);
        if (f == NULL)
        {
            printf ("%s\n", adios_errmsg());
            return -1;
        }

        ADIOS_VARINFO * v = adios_inq_var (f, "dpot");
        //ADIOS_VARINFO * v = adios_inq_var (f, "f0_f");

        start[0] = 0;
        count[0] = v->dims[0]; 

        start[1] = rank;
        count[1] = 4;
       

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

        adios_read_close (f);
        MPI_Barrier (comm);
        adios_read_finalize_method (method);

        double * data_compressed = 0;
#if 0
        // Compress the data using ZFP
        int compressed_size = compressor (data, count[0]*count[1], tolerance, &data_compressed);
        printf ("compression ratio  = %f\n", ((double) count[0]*count[1]*8) / compressed_size);

        decompressor (data, count[0]*count[1], tolerance,
                      data_compressed, compressed_size);
#endif

#if 1
        int iflag = 1; //0 -> float, 1 -> double
        int out_size;
        unsigned char* mgard_comp_buff = mgard_compress (iflag, data, &out_size,
                                               count[1], count[0],  &tolerance );

        printf ("In size:  %10ld  Out size: %10d  Compression ratio: %f\n", v->dims[0]*8, out_size, (double) v->dims[0]*8/out_size);
#endif
    }

    //if (rank == 0) printf ("io time = %f\n", iotime);
    MPI_Finalize ();
    return 0;
}
