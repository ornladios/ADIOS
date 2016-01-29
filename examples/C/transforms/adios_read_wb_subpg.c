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
    int         i, datasize, ndim;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * wbsel;
    double * data = NULL;

    MPI_Init (&argc, &argv);
#ifdef WITH_NCSU_TIMER
    timer_init();
#endif

    adios_read_init_method (method, comm, NULL);

    ADIOS_FILE * f = adios_read_open_file ("adios_global.bp", method, comm);
    ADIOS_VARINFO * varinfo = adios_inq_var (f, "temperature");
    adios_inq_var_blockinfo(f, varinfo);
    if (varinfo)
    {
        ndim = varinfo->ndim;
        assert(ndim == 2);

        //int nranks = varinfo->dims[0];
        assert(varinfo->dims[1] == 10);

        assert(varinfo->nsteps >= 1);
        assert(varinfo->nblocks[0] >= 3);
        assert(varinfo->blockinfo[2].count[0] == 1);
        assert(varinfo->blockinfo[2].count[1] == 10);
        assert(varinfo->blockinfo[2].start[0] == 2);
        assert(varinfo->blockinfo[2].start[1] == 0);

        datasize = 1 * 6 * sizeof(double);
        data = (double *)malloc(datasize);

        wbsel = adios_selection_writeblock(2);
        wbsel->u.block.is_sub_pg_selection = 1;
        wbsel->u.block.element_offset = 2;
        wbsel->u.block.nelements = 6;

        adios_schedule_read (f, wbsel, "temperature", 0, 1, data);
        adios_perform_reads (f, 1);

        printf("Sub-PG writeblock for block %d reading elements in linear range [%" PRIu64 ", %" PRIu64 "):\n",
               wbsel->u.block.index, wbsel->u.block.element_offset, wbsel->u.block.element_offset + wbsel->u.block.nelements);
        printf("[ ");
        for (i = 0; i < wbsel->u.block.nelements; i++)
            printf("%.0lf ", data[i]);
        printf("]\n");

        adios_selection_delete (wbsel);
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
