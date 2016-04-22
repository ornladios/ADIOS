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

int main (int argc, char ** argv) 
{
    int         rank, size, i, j, datasize=0;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;
    void * data = NULL;
   
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=4");

    ADIOS_FILE * f = adios_read_open_file ("adios_global.bp", method, comm);
    ADIOS_VARINFO * varinfo = adios_inq_var (f, "temperature");
    if (varinfo)
    {
        printf ("nsteps = %d\n", varinfo->nsteps);
        for (i = 0; i < varinfo->nsteps; i++)
        {
            printf ("step: %d has %d", i, varinfo->nblocks[i]);
        }
        printf ("\n");

        adios_inq_var_blockinfo (f, varinfo);
        for (i = 0; i < varinfo->sum_nblocks; i++)
        {
            printf ("block[%d]: ", i);
            printf ("start=");
            for (j = 0; j < varinfo->ndim; j++)
            {
                printf ("%" PRIu64, varinfo->blockinfo[i].start[j]);
            }
            printf ("count=");
            for (j = 0; j < varinfo->ndim; j++)
            {
                printf ("%" PRIu64, varinfo->blockinfo[i].count[j]);
            }
            printf ("\n");
        }

        datasize = 8;
        for (i = 0; i < varinfo->ndim; i++)
        {
            datasize *= varinfo->blockinfo[rank].count[i];
        }

        data = malloc (datasize);

        sel = adios_selection_writeblock (rank);
        adios_schedule_read (f, sel, "temperature", 0, 1, data);
        adios_perform_reads (f, 1);
        adios_selection_delete (sel);
    }

    adios_free_varinfo (varinfo);
    adios_read_close (f);

    printf ("data:\n");
    for (i = 0; i < datasize/8; i ++)
    {
        printf ("%7.4f ", * ((double *)data + i));
    }
    printf ("\n");
    free (data);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    MPI_Finalize ();
    return 0;
}
