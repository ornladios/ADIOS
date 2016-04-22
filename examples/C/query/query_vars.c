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
        printf ("  %3d      %4" PRIu64 "    %4" PRIu64 "      %g   %g   %g\n",
                n, points[2*n],points[2*n+1],P[n],V[n],T[n]);
    }
    printf ("\n");
}

void query_OneBoundBoxForAllVars(ADIOS_FILE* f)
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

    if (q!= NULL) {
        int timestep = 0;
        int64_t batchSize = 20;
        int nBatches = 1;
        while (1) {
            /* Evaluate query, get the list of points (of a limited number at once) */
            ADIOS_SELECTION* hits = NULL;
            int hasMore = adios_query_evaluate(q, box, timestep, batchSize, &hits);
            printf("Number of hits returned in batch %d = %" PRIu64 " \n",nBatches, hits->u.points.npoints);

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

            /* free resources used only in this batch */
            free(hits->u.points.points);
            adios_selection_delete(hits);

            if (hasMore <= 0) {
                break;
            }
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
    adios_read_init_method(ADIOS_READ_METHOD_BP,0,"");

    f = adios_read_open_file (filename, ADIOS_READ_METHOD_BP, comm_dummy);
    if (f == NULL) {
        printf ("::%s\n", adios_errmsg());
        return 1;
    }

    query_OneBoundBoxForAllVars(f); 
    //testNoBoxOnSelection(f);
    // testDefaultBoundBox(f);
    //testMultiBoundBox(f);
    //testAllDifferentBoundBoxes(f);
    //testUseOneWriteBlock(f, 0); 
    adios_read_close(f);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    return 0;
}
