/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write variables along with an unstructured mesh. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"




int main (int argc, char ** argv ) 
{
    int         rank, size, i, y;
    MPI_Comm    comm = MPI_COMM_WORLD;

    int         npoints, num_cells;
    float       *points;
    int         *cells;
    double      *N; // node centered variable
    double      *C; // cell centered variable
	

    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    int n = size; // number of points in y direction
    npoints = n*(size+1);       // global number of points
    num_cells = 2*(n-1)*size;   // global number of cells

    // local mesh + data
    points = malloc (sizeof(float) * 2*n*2);  // 2n points with 2 (x-y) coords
    cells = malloc (sizeof(int) * 2*(n-1)*3); // 2(n-1) triangles
    N = malloc (sizeof(double) * 2*n);           // 2n data points on the points
    C = malloc (sizeof(double) * 2*(n-1));       // 2(n-1) data points on the cells

    // generate local data 
    // points (rank,0), (rank,1),...,(rank,n-1), 
    //        (rank+1,0) ... (rank+1,n-1)
    int lp = 2*n; // number of points in this local array
    int op = n*rank; // offset in global points array (next process overwrites our second half)
    for (y=0; y<n; y++) {
        points[2*y] = rank;
        points[2*y+1] = y+rank*0.1;
    }
    // second half are ghost cells with next process
    for (y=0; y<n; y++) {
        points[2*n+2*y] = rank+1;
        points[2*n+2*y+1] = y+(rank+1)*0.1;
    }

    // cells (0--n--n+1) (0--n+1--1) (1--n+1--n+2) (1--n+2--2)... (n-2, 2n-1, n-1) in local terms
    // 0 point local is rank*n in global = op
    int lc = 2*(n-1);
    int oc = lc*rank;
    for (i=0; i<n-1; i++) {
        // two triangles define one rectangle 
        int p=6*i;
        cells[p+0] = op+i; cells[p+1] = op+n+i;   cells[p+2] = op+n+i+1;
        cells[p+3] = op+i; cells[p+4] = op+n+i+1; cells[p+5] = op+i+1;
    }

    for (i=0; i<lp; i++)
        N[i] = rank*n+i;

    for (i=0; i<lc; i++)
        C[i] = rank*lc+i;

    adios_init ("triangle2d.xml", comm);

    adios_open (&adios_handle, "tri", "triangle2d.bp", "w", comm);

    adios_groupsize = 7*sizeof(int) \
	+ sizeof(double) * (lp) \
	+ sizeof(double) * (lc) \
	+ sizeof(float) * (lp) * (2) \
	+ sizeof(int) * (lc) * (3);

    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "nproc", &size);
    adios_write (adios_handle, "npoints", &npoints);
    adios_write (adios_handle, "num_cells", &num_cells);
    adios_write (adios_handle, "lp", &lp);
    adios_write (adios_handle, "op", &op);
    adios_write (adios_handle, "lc", &lc);
    adios_write (adios_handle, "oc", &oc);
    adios_write (adios_handle, "cells", cells);
    adios_write (adios_handle, "points", points);
    adios_write (adios_handle, "N", N);
    adios_write (adios_handle, "C", C);

    adios_close (adios_handle);

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}
