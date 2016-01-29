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

#ifdef WITH_NCSU_TIMER
#include <timer.h>
#endif
#include "mpi.h"
#include "adios_read.h"

int main (int argc, char ** argv)
{
    int         i, j, datasize;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel1;
    double * data = NULL;
    uint64_t start[2], count[2];

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

        assert(varinfo->ndim == 2);

        nranks = varinfo->dims[0];
        assert(nranks % 4 == 0);
        assert(varinfo->dims[1] == 10);

        datasize = (nranks / 2) * varinfo->dims[1] * sizeof(double);
        data = malloc (datasize);

        start[0] = nranks / 4;
        start[1] = 2;
        count[0] = nranks / 2;
        count[1] = 6;

        sel1 = adios_selection_boundingbox (varinfo->ndim, start, count);

        adios_schedule_read (f, sel1, "temperature", 0, 1, data);
        adios_perform_reads (f, 1);

        printf("Subvolume at (%" PRIu64 ",%" PRIu64 ") of size (%" PRIu64 ",%" PRIu64 "):\n", start[0], start[1], count[0], count[1]);
        for (i = 0; i < count[0]; i++) {
            printf("[ ");
            for (j = 0; j < count[1]; j++) {
                printf("%.0lf ", data[i * count[1] + j]);
            }
            printf("]\n");
        }

        adios_selection_delete (sel1);
    }

    adios_free_varinfo (varinfo);
    adios_read_close (f);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);


#ifdef WITH_NCSU_TIMER
    printf("[TIMERS] ");
    timer_result_t *results = timer_get_results_sorted();
    for (i = 0; i < timer_get_num_timers(); i++) {
        printf("%s: %0.4lf ", results[i].name, results[i].time);
    }
    printf("\n");
    free(results);
#endif

#ifdef WITH_NCSU_TIMER
    timer_finalize();
#endif
    MPI_Finalize ();
    return 0;
}
