/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mpi.h"
#include "adios_read.h"

#ifdef WITH_NCSU_TIMER
#include <timer.h>
#endif

int main (int argc, char ** argv)
{
    int         i, datasize,  ndim;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel1;
    double * data = NULL;
    uint64_t npoints, * points;

    MPI_Init (&argc, &argv);
#ifdef WITH_NCSU_TIMER
    timer_init();
#endif
    adios_read_init_method (method, comm, NULL);

    ADIOS_FILE * f = adios_read_open_file ("adios_global.bp", method, comm);
    ADIOS_VARINFO * varinfo = adios_inq_var (f, "temperature");
    if (varinfo)
    {
        int nranks;

        ndim = varinfo->ndim;
        assert(ndim == 2);

        nranks = varinfo->dims[0];
        assert(varinfo->dims[1] == 10);

        datasize = (nranks / 2) * varinfo->dims[1] * sizeof(double);
        data = malloc (datasize);

        npoints = 2 * nranks;
        points = malloc(npoints * ndim * sizeof(uint64_t));
        for (i = 0; i < npoints; i += 2) {
            points[i * ndim + 0] = i/2; // rank/row
            points[i * ndim + 1] = 2; // col
            points[(i+1) * ndim + 0] = i/2; // rank/row
            points[(i+1) * ndim + 1] = 6; // col
        }

        sel1 = adios_selection_points(ndim, npoints, points);

        adios_schedule_read (f, sel1, "temperature", 0, 1, data);
        adios_perform_reads (f, 1);

        printf("Points read (columns 2 and 6 in all rows):\n");
        for (i = 0; i < npoints; i++) {
            printf("(%" PRIu64 ",%" PRIu64 ") = %.0lf\n", points[i*ndim+0], points[i*ndim+1], data[i]);
        }

        adios_selection_delete (sel1);
    }

    adios_free_varinfo (varinfo);
    adios_read_close (f);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
#ifdef WITH_NCSU_TIMER
    timer_finalize();
#endif
    MPI_Finalize ();
    return 0;
}
