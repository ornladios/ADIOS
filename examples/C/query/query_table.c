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
    int32_t *data = (int32_t*) malloc (Npoints * ncols * sizeof(int32_t));

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
            printf ("%5d      ", data[n*ncols+i]);
        }
        printf ("\n");
    }
    printf ("\n");

    /* Free memory */
    free (data);
    adios_selection_delete(box);
    adios_free_varinfo (v);

}

void print_points (ADIOS_SELECTION *hits, int *KE)
{
    int n;
    int Npoints = hits->u.points.npoints;
    uint64_t *points = hits->u.points.points;
    /* Note: hits->u.points.ndim is 2 in this example */

    printf ("\nHit           i       j    Kinetic E\n");
    printf ("----------------------------------------------\n");
    for (n=0; n<Npoints; n++) {
        printf ("  %3d      %4" PRIu64 "    %4" PRIu64 "      %d\n",
                n, points[2*n],points[2*n+1],KE[n]);
    }
    printf ("\n");
}

void query_columns(ADIOS_FILE* f, int nrows, int ncols)
{

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

    if (q!= NULL) {
        int timestep = 0;
        int64_t batchSize = 20;
        int nBatches = 1;
        while (1) {
            /* Evaluate query, get the list of points (of a limited number at once) */
            ADIOS_SELECTION* hits = NULL;
            int hasMore = adios_query_evaluate(q, col3, timestep, batchSize, &hits);
            printf("Number of hits returned in batch %d = %" PRIu64 " \n",nBatches, hits->u.points.npoints);

            if (hits->u.points.npoints > 0) {
                /* Read the data of those points */
                int *KE = (int *) malloc (sizeof(double)*hits->u.points.npoints);
                adios_schedule_read (f, hits, "A", timestep, 1, KE);
                adios_perform_reads (f, 1);

                print_points (hits, KE);
                free (KE);

                read_rows (f, hits, ncols);
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
    adios_query_free(q);
    adios_selection_delete(col1);
    adios_selection_delete(col2);
    adios_selection_delete(col3);
}

int main (int argc, char ** argv) 
{
    ADIOS_FILE * f;
    MPI_Comm    comm_dummy = 0;  /* MPI_Comm is defined through adios.h/adios_read.h */
    adios_read_init_method(ADIOS_READ_METHOD_BP,0,"");

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

    query_columns(f, v->dims[0], v->dims[1]); 
    adios_read_close(f);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    return 0;
}
