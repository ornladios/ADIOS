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
    int         rank, size, i, datasize;
    //int j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP_AGGREGATE;
    ADIOS_SELECTION * sel1;
    //ADIOS_SELECTION * sel2;
    //ADIOS_VARCHUNK * chunk = 0;
    void * data = NULL;
    uint64_t start[2], count[2]; 
    //uint64_t npoints, * points;
   
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "max_chunk_size=4;verbose=3;num_aggregators=2");

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
        printf ("ndim = %d\n", varinfo->ndim);
        printf ("dims = (");
        datasize = 8;
        for (i = 0; i < varinfo->ndim; i++)
        {
            datasize *= varinfo->dims[i];
            printf ("%" PRIu64, varinfo->dims[i]);
            if (i != varinfo->ndim - 1)
            {
                printf (",");
            }
        }
        printf (")\n");
        data = malloc (datasize);
/*
        adios_inq_var_blockinfo (f, varinfo);
        for (i = 0; i < varinfo->sum_nblocks; i++)
        {
            printf ("block[%d]: ", i);
            printf ("start=");
            for (j = 0; j < varinfo->ndim; j++)
            {
                printf ("%lu ", varinfo->blockinfo[i].start[j]);
            }
            printf ("count=");
            for (j = 0; j < varinfo->ndim; j++)
            {
                printf ("%lu ", varinfo->blockinfo[i].count[j]);
            }
            printf ("\n");
        }
*/
        for (i = 0; i < varinfo->ndim; i++)
        {
            start[i] = 0;
            count[i] = varinfo->dims[i];
        }

        sel1 = adios_selection_boundingbox (varinfo->ndim, start, count);
/*
        npoints = 1;
        for (i = 0; i < varinfo->ndim; i++)
        {
            npoints *= count[i];
        }

        points = (uint64_t *) malloc (npoints * varinfo->ndim * 8);
        for (i = 0; i < npoints; i++)
        {
            uint64_t temp = i;
            for (j = varinfo->ndim - 1; j > -1; j--)
            {
                points[i * varinfo->ndim + j] = temp % count[j];
                temp = temp/count[j];
            }
        }
*/
/*
        for (i = 0; i < npoints; i++)
        {
            printf ("(");
            for (j = 0; j < varinfo->ndim; j++)
            {
                printf ("%lu ",points[i * varinfo->ndim + j]);
            }
            printf (")");
        }

        printf ("\n");
*/
/*
        sel2 = adios_selection_points (varinfo->ndim, npoints, points);
*/
        adios_schedule_read (f, sel1, "temperature", 0, 1, data);
/*
        adios_schedule_read (f, sel2, "temperature", 0, 1, data);
*/
        adios_perform_reads (f, 1);
#if 0
        while (adios_check_reads (f, &chunk) > 0)
        {
            datasize = 1;
            for (i = 0; i < varinfo->ndim; i++)
            {
                datasize *= chunk->sel->u.bb.count[i];
            }
            printf ("data:\n");
/*
            for (i = 0; i < datasize; i ++)
            {
                printf ("%7.4f ", * ((double *)chunk->data + i));
            }
            printf ("\n");
*/
            for (i = 0; i < 10; i ++)
            {
                printf ("%7.4f ", * ((double *)chunk->data + i));
            }
            printf ("\n");

            adios_free_chunk (chunk);
        }
#endif
        adios_selection_delete (sel1);
/*
        adios_selection_delete (sel2);
*/
    }

    adios_free_varinfo (varinfo);
    adios_read_close (f);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    MPI_Finalize ();
    return 0;
}
