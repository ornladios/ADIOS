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
#include <string.h>
#include "adios.h"       /* includes dummy MPI */
#include "adios_query.h" /* includes the read API */

const char filename[] = "vars.bp";

void print_points (ADIOS_SELECTION *hits, uint64_t *wboffs, double *P, double *V, double *T)
{
    int n;
    int Npoints = hits->u.points.npoints;
    uint64_t *points = hits->u.points.points;
    /* Note: hits->u.points.ndim is 2 in this example */

    printf ("\nHit           i       j        P      V     T\n");
    printf ("----------------------------------------------\n");
    for (n=0; n<Npoints; n++) {
        printf ("  %3d      %4" PRIu64 "    %4" PRIu64 "      %g   %g   %g\n",
                n, points[2*n]+wboffs[0],points[2*n+1]+wboffs[1],
                P[n],V[n],T[n]);
    }
    printf ("\n");
}


void print_block (ADIOS_VARINFO *vi, int blockid, double *P, double *V, double *T)
{
    int n;
    printf ("    Block %d of dimensions { ", blockid);
    for (n=0; n < vi->ndim; n++) {
        printf ("%" PRIu64 " ", vi->blockinfo[blockid].count[n]);
    }
    printf ("}\n");
}

void query_OneBoundBoxForAllVars(ADIOS_FILE* f, enum ADIOS_QUERY_METHOD method)
{
    printf("\n====== querying with one bound box for all variables =======\n");
    uint64_t start[] = {0,0};
    uint64_t count[] = {5,6};
    int i;
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
    ADIOS_VARINFO *vP = adios_inq_var (f, "P");
    adios_inq_var_blockinfo(f, vP);


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

            if (result->nselections == 0) {
                printf ("Zero results returned in batch %d\n", nBatches);
                break;
            }

            if (result->selections->type == ADIOS_SELECTION_POINTS)
            {
                // we have selection(s) and each one contains the points
                printf("Number of hits returned in batch %d = %" PRIu64 " points in %d containers\n",
                        nBatches, result->npoints, result->nselections);

                for (i=0; i < result->nselections; i++)
                {
                    ADIOS_SELECTION* hits = &(result->selections[i]);
                    const ADIOS_SELECTION_POINTS_STRUCT * pts = &(hits->u.points);
                    uint64_t * wboffs = calloc (pts->ndim, sizeof(uint64_t));
                    if (pts->container_selection &&
                            pts->container_selection->type == ADIOS_SELECTION_WRITEBLOCK)
                    {
                        int i;
                        int blockidx = pts->container_selection->u.block.index;
                        // calculate actual block index if multiple timesteps are available
                        for (i = 0; i < timestep-1; i++)
                            blockidx += vP->nblocks[i];
                        // now record the offset of this block in global space
                        // point coordinates are relative to block
                        for (i = 0; i < pts->ndim; ++i) {
                            wboffs[i] = vP->blockinfo[blockidx].start[i];
                        }
                    }

                    /* Read the data of those points */
                    double *P = (double *) malloc (sizeof(double)*hits->u.points.npoints);
                    double *V = (double *) malloc (sizeof(double)*hits->u.points.npoints);
                    double *T = (double *) malloc (sizeof(double)*hits->u.points.npoints);
                    adios_schedule_read (f, hits, "P", timestep, 1, P);
                    adios_schedule_read (f, hits, "V", timestep, 1, V);
                    adios_schedule_read (f, hits, "T", timestep, 1, T);
                    adios_perform_reads (f, 1);

                    print_points (hits, wboffs, P, V, T);
                    free (P);
                    free (V);
                    free (T);
                    free (wboffs);

                }
            }
            else if (result->selections->type == ADIOS_SELECTION_WRITEBLOCK)
            {
                // we have multiple selections, each one is a writeblock
                printf("Number of blocks returned in batch %d = %d \n",nBatches,
                        result->nselections);

                for (i=0; i < result->nselections; i++)
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

            /* free resources used only in this batch */
            for (i=0; i < result->nselections; i++) {
                adios_selection_delete (&(result->selections[i]));
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
    adios_query_free(q3);
    adios_query_free(q4);
    adios_query_free(q);
    adios_selection_delete(box);
    adios_free_varinfo(vP);
}

void printUsage(char *prgname)
{
    printf ("Usage: %s [fastbit|alacrity]\n"
           "  Choose the query method to use.\n"
           "  For ALACRITY, you need to build write_vars with and ADIOS which has ALACRITY transformation.\n"
           "  For FastBit, you need to run 'adios_index_fastbit vars.bp' to generate the index 'vars.idx'.\n"
           ,prgname);
}

int main (int argc, char ** argv) 
{
    ADIOS_FILE * f;
    MPI_Comm    comm_dummy = 0;  /* MPI_Comm is defined through adios.h/adios_read.h */
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
        } else if (!strncasecmp (argv[1], "minmax", 7)) {
            if (adios_query_is_method_available(ADIOS_QUERY_METHOD_MINMAX)) {
                query_method = ADIOS_QUERY_METHOD_MINMAX;
                printf ("Set query method to MINMAX\n");
            } else {
                printf ("ERROR: The MINMAX method is not available in this ADIOS build.\n"
                        "Try ALACRITY or FASTBIT query methods\n");
                return 1;
            }
        } else {
            printUsage(argv[0]);
            return 1;
        }
    } else {
        query_method = ADIOS_QUERY_METHOD_UNKNOWN;
        //printUsage(argv[0]);
        //return 1;
    }

    f = adios_read_open_file (filename, ADIOS_READ_METHOD_BP, comm_dummy);
    if (f == NULL) {
        printf ("::%s\n", adios_errmsg());
        return 1;
    }

    query_OneBoundBoxForAllVars(f, query_method);
    //testNoBoxOnSelection(f);
    // testDefaultBoundBox(f);
    //testMultiBoundBox(f);
    //testAllDifferentBoundBoxes(f);
    //testUseOneWriteBlock(f, 0); 
    adios_read_close(f);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    return 0;
}
