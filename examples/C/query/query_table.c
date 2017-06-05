/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: query over columns of a table (2D variables)
 *
 * How to run: ./query_table
 * Output: standard output
 *
*/

#ifndef _NOMPI
#define _NOMPI
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "adios.h"       /* includes dummy MPI */
#include "adios_query.h" /* includes the read API */

const char filename[] = "table.bp";

void read_rows (ADIOS_FILE *f, ADIOS_SELECTION *hits, int ncols) 
{
    /* Demonstrates the more involved work to read the whole rows of the table
       where the query was true. The point list returned by the evaluation
       (on any result column) can be used here, we just need the row numbers,
       i.e. the first dimensions of the points 
     */
    int n, i;
    int Npoints = hits->u.points.npoints;
    uint64_t *points = hits->u.points.points;

    /* Just read all data at once here, assuming it fits in memory */
    float *data = (float*) malloc (Npoints * ncols * sizeof(float));

    /* We use one bounding box selection of a single row repeatedly, 
       modifying the offset based on the actual point's row number.
     */
    uint64_t offs[] = {0,0};  
    uint64_t cnt[] = {1,ncols}; // select one row
    ADIOS_SELECTION *box = adios_selection_boundingbox(2, offs, cnt);
    
    printf ("\nRead all rows that match the query:\n");
    for (n=0; n<Npoints; n++) {
        box->u.bb.start[0] = points[2*n]; // set actual row 
        adios_schedule_read (f, box, "A", 0, 1, data+n*ncols);
    }
    adios_perform_reads (f, 1);

    /* Get the column names, just for fun */
    ADIOS_VARINFO *v = adios_inq_var (f, "Columns");
    char *Columns = (char*) malloc (v->dims[0]*v->dims[1]);
    //printf ("Allocate %" PRIu64 " bytes for column names\n", v->dims[0]*v->dims[1]);
    adios_schedule_read (f, NULL, "Columns", 0, 1, Columns);
    adios_perform_reads (f, 1);

    /* Print result */
    printf ("\n");
    for (i=0; i<ncols; i++) {
        printf ("%s ", Columns + i*v->dims[1]);
    }
    printf ("\n----------------------------------------------------------------------------\n");
    for (n=0; n<Npoints; n++) {
        for (i=0; i<ncols; i++) {
            printf ("%5g      ", data[n*ncols+i]);
        }
        printf ("\n");
    }
    printf ("\n");

    /* Free memory */
    free (data);
    adios_selection_delete(box);
    adios_free_varinfo (v);
    free (Columns);
}

void print_points (ADIOS_SELECTION *hits, float *KE, uint64_t *wboffs)
{
    int n;
    int Npoints = hits->u.points.npoints;
    uint64_t *points = hits->u.points.points;
    /* Note: hits->u.points.ndim is 2 in this example */

    printf ("\nHit           i       j    Kinetic E\n");
    printf ("----------------------------------------------\n");
    for (n=0; n<Npoints; n++) {
        printf ("  %3d      %4" PRIu64 "    %4" PRIu64 "      %g\n",
                n, points[2*n]+wboffs[0],points[2*n+1]+wboffs[1],KE[n]);
    }
    printf ("\n");
}

void query_columns(ADIOS_FILE* f, enum ADIOS_QUERY_METHOD method, ADIOS_VARINFO *vi)
{
    int nrows = vi->dims[0];
    int ncols = vi->dims[1];
    printf("\n====== querying over columns of a table  =======\n");
    uint64_t offs1[] = {0,1};  // element
    uint64_t offs2[] = {0,2};  // potential
    uint64_t offs3[] = {0,3};  // kinetic energy
    uint64_t cnt[] = {nrows, 1};  
    ADIOS_SELECTION* col1 = adios_selection_boundingbox(2, offs1, cnt);
    ADIOS_SELECTION* col2 = adios_selection_boundingbox(2, offs2, cnt);
    ADIOS_SELECTION* col3 = adios_selection_boundingbox(2, offs3, cnt);

    ADIOS_QUERY *q1, *q2, *q;
    q1 = adios_query_create(f, col1, "A", ADIOS_EQ, "0"); // select Carbon
    q2 = adios_query_create(f, col2, "A", ADIOS_LTEQ, "96"); // select Potential <= 96
    q  = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);
    printf("Query: %s\n",q->condition);
    adios_query_set_method (q, method);

    adios_inq_var_blockinfo(f, vi);

    if (q!= NULL) {
        int timestep = 0;
        int64_t batchSize = 20;
        int nBatches = 1;
        while (1) {
            /* Evaluate query, get the list of points (of a limited number at once) */
            ADIOS_QUERY_RESULT *result = adios_query_evaluate(q, col3, timestep, batchSize);

            if (result->status == ADIOS_QUERY_RESULT_ERROR) {
                printf ("Query evaluation failed with error: %s\n", adios_errmsg());
                break;
            }

            if (result->nselections == 0) {
                printf ("Zero results returned in batch %d\n", nBatches);
                break;
            }

            printf("Number of hits returned in batch %d = %" PRIu64 " points in %d containers\n",
                    nBatches, result->npoints, result->nselections);

            int n;
            for (n = 0; n < result->nselections; n++)
            {
                ADIOS_SELECTION* hits = &(result->selections[n]);
                const ADIOS_SELECTION_POINTS_STRUCT * pts = &(hits->u.points);
                uint64_t * wboffs = calloc (pts->ndim, sizeof(uint64_t));
                if (pts->container_selection &&
                        pts->container_selection->type == ADIOS_SELECTION_WRITEBLOCK)
                {
                    int i;
                    int blockidx = pts->container_selection->u.block.index;
                    // calculate actual block index if multiple timesteps are available
                    for (i = 0; i < timestep-1; i++)
                        blockidx += vi->nblocks[i];
                    // now record the offset of this block in global space
                    // point coordinates are relative to block
                    for (i = 0; i < pts->ndim; ++i) {
                        wboffs[i] = vi->blockinfo[blockidx].start[i];
                    }
                }

                /* Read the data of those points */
                float *KE = (float *) malloc (sizeof(double)*pts->npoints);
                adios_schedule_read (f, hits, "A", timestep, 1, KE);
                adios_perform_reads (f, 1);

                print_points (hits, KE, wboffs);
                free (KE);
                free (wboffs);

                read_rows (f, hits, ncols);
            } 

            /* free resources used only in this batch */
            for (n=0; n < result->nselections; n++) {
                adios_selection_delete (&(result->selections[n]));
            }

            if (result->status == ADIOS_QUERY_NO_MORE_RESULTS) {
                free (result);
                break;
            }
            free (result);
            nBatches++;
        }

    }
    /* free resources used for the query */
    adios_query_free(q1);
    adios_query_free(q2);
    adios_query_free(q);
    adios_selection_delete(col1);
    adios_selection_delete(col2);
    adios_selection_delete(col3);
}

void printUsage(char *prgname)
{
    printf ("Usage: %s [fastbit|alacrity]\n"
           "  Choose the query method to use.\n"
           "  For ALACRITY, you need to build write_table with and ADIOS which has ALACRITY transformation.\n"
           "  For FastBit, you need to run 'adios_index_fastbit table.bp' to generate the index 'table.idx'.\n"
           ,prgname);
}

int main (int argc, char ** argv) 
{
    ADIOS_FILE * f;
    MPI_Comm    comm_dummy = 0;  /* MPI_Comm is defined through adios.h/adios_read.h */
    adios_read_init_method(ADIOS_READ_METHOD_BP,0,"");
    enum ADIOS_QUERY_METHOD query_method;

    adios_read_init_method(ADIOS_READ_METHOD_BP,0,"");

    if (!adios_query_is_method_available(ADIOS_QUERY_METHOD_ALACRITY) && 
        !adios_query_is_method_available(ADIOS_QUERY_METHOD_FASTBIT)) 
    {
        printf ("This query test on tabular data is only supported by accurate "
                "point-based query methods like FASTBIT and ALACRITY. "
                "No such method is available in this ADIOS build.\n");
        return 1;
    }
    if (argc > 1) {
        if (!strncasecmp (argv[1], "alacrity", 8)) {
            if (adios_query_is_method_available(ADIOS_QUERY_METHOD_ALACRITY)) {
                query_method = ADIOS_QUERY_METHOD_ALACRITY;
                printf ("Set query method to ALACRITY\n");
            } else {
                printf ("ERROR: The ALACRITY method is not available in this ADIOS build.\n"
                        "Try FASTBIT but first run 'adios_index_fastbit table.bp'\n");
                return 1;
            }
        } else if (!strncasecmp (argv[1], "fastbit", 7)) {
            if (adios_query_is_method_available(ADIOS_QUERY_METHOD_FASTBIT)) {
                query_method = ADIOS_QUERY_METHOD_FASTBIT;
                printf ("Set query method to FASTBIT\n");
            } else {
                printf ("ERROR: The FASTBIT method is not available in this ADIOS build.\n"
                        "Try ALACRITY but first run the write_table code with alacrity transformation!\n");
                return 1;
            }
        } else {
            printUsage(argv[0]);
            return 1;
        }
    } else {
        printUsage(argv[0]);
        return 1;
    }

    printf("File : %s\n",filename);
    f = adios_read_open_file (filename, ADIOS_READ_METHOD_BP, comm_dummy);
    if (f == NULL) {
        printf ("::%s\n", adios_errmsg());
        return 1;
    }

    ADIOS_VARINFO *v = adios_inq_var (f, "A");
    if (v == NULL) {
        printf ("Error: did not find variable A in file %s::%s\n",
                filename, adios_errmsg());
        return 2;
    }
    if (v->ndim != 2) {
        printf ("Error: Variable A is expected to have 2 dimensions\n");
        return 2;
    }
    if (v->dims[1] < 4) {
        printf ("Error: Table A is expected to have at least 4 columns\n");
        return 2;
    }

    printf ("Variable A has %" PRIu64 " rows and %" PRIu64 " columns\n", v->dims[0], v->dims[1]);

    query_columns(f, query_method, v);
    adios_free_varinfo(v);
    adios_read_close(f);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    return 0;
}
