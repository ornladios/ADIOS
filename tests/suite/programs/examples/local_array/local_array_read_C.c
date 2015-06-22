/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*************************************************************/
/*          Example of reading arrays in ADIOS               */
/*    which were written from the same number of processors  */
/*                                                           */
/*        Similar example is manual/2_adios_read.c           */
/*************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include "mpi.h"
#include "adios_read.h"
#include "core/adios_logger.h"

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, j;
    int         NX, NY; 
    double      *t;
    int         *p;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_read_init_method (method, comm, "verbose=3");
    adios_logger_open ("log_read_C", rank);

    strcpy (filename, "local_array_C.bp");
    ADIOS_FILE * f = adios_read_open (filename, method, comm, ADIOS_LOCKMODE_NONE, 0);

    /* Specify a selection that points to a specific writer's block */
    sel = adios_selection_writeblock (rank);

    /* First get the scalars to calculate the size of the arrays */
    adios_schedule_read (f, sel, "NX", 0, 1, &NX);
    adios_schedule_read (f, sel, "NY", 0, 1, &NY);
    adios_perform_reads (f, 1);

    log_test("rank=%d: NX=%d NY=%d\n", rank, NX, NY);

    /* Allocate space for the arrays */
    t = (double *) malloc (NX*NY*sizeof(double));
    p = (int *) malloc (NX*sizeof(int));

    /* Read the arrays */
    adios_schedule_read (f, sel, "var_double_2Darray", 0, 1, t);
    adios_schedule_read (f, sel, "var_int_1Darray", 0, 1, p);
    adios_perform_reads (f, 1);

    /* At this point, we have the data in memory */
    log_test("rank=%d: p = [%d", rank, p[0]);
    for (i = 1; i < NX; i++)
        log_test(", %d", p[i]);
    log_test("]\n");
    
    log_test("rank=%d: t[5,*] = [%6.2f", rank, t[5*NY]);
    for (j = 1; j < NY; j++)
        log_test(", %6.2f", t[5*NY+j]);
    log_test("]\n");

    free (t);
    free (p);

    adios_read_close (f);
    MPI_Barrier (comm);
    adios_read_finalize_method (method);
    adios_logger_close();
    MPI_Finalize ();

    return 0;
}
