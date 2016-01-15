/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: query over multiple variables 
 *
 * How to run: ./query_vars
 * Output: standard output
 *
*/

#ifndef _NOMPI
#define _NOMPI
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "adios.h"       /* includes dummy MPI */
#include "adios_query.h" /* includes the read API */

const char filename[] = "vars.bp";

void print_points (ADIOS_SELECTION *hits, double *P, double *V, double *T)
{
    int n;
    int Npoints = hits->u.points.npoints;
    uint64_t *points = hits->u.points.points;
    /* Note: hits->u.points.ndim is 2 in this example */

    printf ("\nHit           i       j        P      V     T\n");
    printf ("----------------------------------------------\n");
    for (n=0; n<Npoints; n++) {
        printf ("  %3d      %4lld    %4lld      %g   %g   %g\n",
                n, points[2*n],points[2*n+1],P[n],V[n],T[n]);
    }
    printf ("\n");
}


void print_block (ADIOS_VARINFO *vi, int blockid, double *P, double *V, double *T)
{
    int n;
    printf ("    Block %d of dimensions { ", blockid);
    for (n=0; n < vi->ndim; n++) {
        printf ("%lld ", vi->blockinfo[blockid].count[n]);
    }
    printf ("}\n");
}

void query_OneBoundBoxForAllVars(ADIOS_FILE* f, enum ADIOS_QUERY_METHOD method)
{
    printf("\n====== querying with one bound box for all variables =======\n");
    uint64_t start[] = {0,0};
    uint64_t count[] = {5,6};

    ADIOS_SELECTION* box = adios_selection_boundingbox(2, start, count);
    ADIOS_QUERY *q1, *q2, *q3, *q4, *q;
    q1 = adios_query_create(f, box, "P", ADIOS_GT, "80.0");
    q2 = adios_query_create(f, box, "P", ADIOS_LT, "90.0");
    q3 = adios_query_combine(q1, ADIOS_QUERY_OP_AND, q2);
    q4 = adios_query_create(f, box, "V", ADIOS_LTEQ, "50.0");
    q  = adios_query_combine(q3, ADIOS_QUERY_OP_AND, q4);
    printf("File : %s\n",filename);
    printf("Query: %s\n",q->condition);
    adios_query_set_method (q, method);

    if (q!= NULL) {
        int timestep = 0;
        int64_t batchSize = 20;
        int nBatches = 1;
        while (1) {
            /* Evaluate query, get the list of points or write blocks (of a limited number at once) */
            ADIOS_QUERY_RESULT *result = adios_query_evaluate(q, box, timestep, batchSize);

            if (result->status == ADIOS_QUERY_RESULT_ERROR) {
                printf ("Query evaluation failed with error: %s\n", adios_errmsg());
                break;
            }

            if (result->nresults == 0) {
                printf ("Zero results returned in batch %d\n", nBatches);
                break;
            }

            if (result->selections->type == ADIOS_SELECTION_POINTS)
            {
                // we have one selection which contains the points
                ADIOS_SELECTION* hits = result->selections;
                printf("Number of hits returned in batch %d = %lld \n",nBatches,
                        hits->u.points.npoints);

                if (hits->u.points.npoints > 0) {
                    /* Read the data of those points */
                    double *P = (double *) malloc (sizeof(double)*hits->u.points.npoints);
                    double *V = (double *) malloc (sizeof(double)*hits->u.points.npoints);
                    double *T = (double *) malloc (sizeof(double)*hits->u.points.npoints);
                    adios_schedule_read (f, hits, "P", timestep, 1, P);
                    adios_schedule_read (f, hits, "V", timestep, 1, V);
                    adios_schedule_read (f, hits, "T", timestep, 1, T);
                    adios_perform_reads (f, 1);

                    print_points (hits, P, V, T);
                    free (P);
                    free (V);
                    free (T);
                }

                /* free resources used in case of point based queries */
                free(hits->u.points.points);
            }
            else if (result->selections->type == ADIOS_SELECTION_WRITEBLOCK)
            {
                // we have multiple selections, each one is a writeblock
                printf("Number of blocks returned in batch %d = %d \n",nBatches,
                        result->nresults);
                ADIOS_VARINFO *vP = adios_inq_var (f, "P");
                adios_inq_var_blockinfo(f, vP);

                int i;
                for (i=0; i < result->nresults; i++)
                {
                    uint64_t nelems = 1;
                    int j;
                    for (j=0; j < vP->ndim; j++) {
                        nelems *= vP->blockinfo[result->selections[i].u.block.index].count[j];
                    }
                    double *P = (double *) malloc (sizeof(double)*nelems);
                    double *V = (double *) malloc (sizeof(double)*nelems);
                    double *T = (double *) malloc (sizeof(double)*nelems);

                    adios_schedule_read (f, &result->selections[i], "P", timestep, 1, P);
                    adios_schedule_read (f, &result->selections[i], "V", timestep, 1, V);
                    adios_schedule_read (f, &result->selections[i], "T", timestep, 1, T);
                    adios_perform_reads (f, 1);

                    print_block (vP, result->selections[i].u.block.index, P, V, T);
                    free (P);
                    free (V);
                    free (T);

                }
            }

            if (result->status == ADIOS_QUERY_NO_MORE_RESULTS) {
                break;
            }
            free (result->selections);
            free (result);
            nBatches++;
        }

    }
    /* free resources used for the query */
    adios_query_free(q1);
    adios_query_free(q2);
    adios_query_free(q3);
    adios_query_free(q4);
    adios_query_free(q);
    adios_selection_delete(box);
}

int main (int argc, char ** argv) 
{
    ADIOS_FILE * f;
    MPI_Comm    comm_dummy = 0;  /* MPI_Comm is defined through adios.h/adios_read.h */
    enum ADIOS_QUERY_METHOD method;

    adios_read_init_method(ADIOS_READ_METHOD_BP,0,"");
    if (adios_query_is_method_available(ADIOS_QUERY_METHOD_ALACRITY)) {
        method = ADIOS_QUERY_METHOD_ALACRITY;
        printf ("Set query method to ALACRITY\n");
    } else if (adios_query_is_method_available(ADIOS_QUERY_METHOD_FASTBIT)) {
        method = ADIOS_QUERY_METHOD_FASTBIT;
        printf ("Set query method to FASTBIT\n");
    } else {
        method = ADIOS_QUERY_METHOD_UNKNOWN;
        printf ("Let the query engine select the query method\n");
    }

    f = adios_read_open_file (filename, ADIOS_READ_METHOD_BP, comm_dummy);
    if (f == NULL) {
        printf ("::%s\n", adios_errmsg());
        return 1;
    }

    query_OneBoundBoxForAllVars(f, method);
    //testNoBoxOnSelection(f);
    // testDefaultBoundBox(f);
    //testMultiBoundBox(f);
    //testAllDifferentBoundBoxes(f);
    //testUseOneWriteBlock(f, 0); 
    adios_read_close(f);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    return 0;
}
