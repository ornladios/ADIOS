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
#include <time.h>
#include "mpi.h"
#include "adios_read.h"
#include <stdint.h>

#define BLOCKING 1
#define NON_BLOCKING 0

#define Z 1
#define Y 2
#define X 4

#define NTIMES 1

#define RS_SCALE (1.0 / (1.0 + RAND_MAX))

double drand (void) 
{
    double d;
    do {
       d = ((((rand () * RS_SCALE) + rand ()) * RS_SCALE + rand ()) * RS_SCALE + rand ()) * RS_SCALE; 
    } while (d >= 1); /* Round off */ 
    return d;
}  

uint64_t irand (uint64_t x) 
{
    return ((x) * drand ());
}

void read_points (char filename [], char varname [], uint32_t npoints)
{
    printf ("===%s===\n", __FUNCTION__);
    // What are the first 2 dimensions
    enum ADIOS_READ_METHOD  method   = ADIOS_READ_METHOD_BP;
    MPI_Comm                comm     = MPI_COMM_WORLD;
    ADIOS_SELECTION         *sel1;
    ADIOS_FILE              *f       = adios_read_open_file (filename, method, comm);
    ADIOS_VARINFO           *varinfo = adios_inq_var (f, varname);

    assert (varinfo);

    //uint32_t                nblocks = varinfo->sum_nblocks;
    uint32_t                i       = 0;
    uint32_t                j       = 0;
    uint32_t                t       = 0;
    uint64_t                *points = 0;
    double                  *data   = 0;

    uint32_t                ndim = varinfo->ndim;

    // This code should only work for 3D data
    // assert (varinfo->ndim == 3);

    /*
    printf ("Dimensions for %s are %s: ", filename, varname);
    for (i = 0; i < varinfo->ndim; i ++) {
        printf ("%d ", varinfo->dims [i]);        
    }
    printf ("\n");
    printf ("Timesteps: %d\n", f->last_step + 1);
    assert (f->last_step == 0);
    */
 
    data    = (double *) malloc (npoints * sizeof (double));
    points  = (uint64_t *) malloc (npoints * varinfo->ndim * sizeof (uint64_t));

    assert (data);
    assert (points);

    for (t = 0; t < NTIMES; t ++) {
        // Generate a random point 
        for (i = 0; i < npoints; i ++) {
            for (j = 0; j < ndim; j ++) {
                points [i * ndim + j] = irand (varinfo->dims [j]);
            }
        }

        sel1 = adios_selection_points (ndim, npoints, points);
        
        adios_schedule_read (f, sel1, varname, 0, 1, data);
        adios_perform_reads (f, BLOCKING);
   
        // Print the points. For double checking 
        /*
        for (i = 0; i < npoints; i++) {
            printf("(%llu,%llu) = %.0lf\n", points [i * ndim + 0], points [i * ndim + 1], data [i]);
        }
        */
        adios_selection_delete (sel1);
    }

    free (points);
    free (data);

    adios_free_varinfo (varinfo);
    adios_read_close (f);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
}

void read_bounding_box (char filename [], char varname [], uint64_t counts [], uint8_t plane, uint8_t is_aligned)
{
    printf ("===%s===\n", __FUNCTION__);
    // What are the first 2 dimensions
    enum ADIOS_READ_METHOD  method   = ADIOS_READ_METHOD_BP;
    MPI_Comm                comm     = MPI_COMM_WORLD;
   
    ADIOS_SELECTION         *sel1;

    ADIOS_FILE              *f       = adios_read_open_file (filename, method, comm);
    ADIOS_VARINFO           *varinfo = adios_inq_var (f, varname);

    assert (varinfo);
    adios_inq_var_blockinfo (f, varinfo);

    //uint32_t                nblocks = varinfo->sum_nblocks;
   
    uint32_t                i       = 0;
    uint32_t                j       = 0;
    uint32_t                t       = 0;
    uint64_t                *starts = 0;
    double                  *data   = 0;

    uint32_t                ndim = varinfo->ndim;
    uint64_t                *pg_dims = varinfo->blockinfo [0].count;
    uint64_t                npoints = 1;

    // This code should only work for 3D data
    // assert (varinfo->ndim == 3);

    printf ("%u dimensions for %s are %s: ", ndim, filename, varname);
    for (i = 0; i < varinfo->ndim; i ++) {
        printf ("%lld ", varinfo->dims [i]);        
    }
    printf ("\n");
    printf ("Timesteps: %d\n", f->last_step + 1);

    // assert (f->last_step == 0);

    #if 1 
        counts [0] = 64;
        counts [1] = 64;
        counts [2] = 64;
    #endif

    for (j = 0; j < ndim; j ++) {
        npoints *= counts [j];
    }

    printf ("npoints: %lld\n", npoints);

    data    = (double *) malloc (npoints * sizeof (double));
    starts  = (uint64_t *) malloc (varinfo->ndim * sizeof (uint64_t));

    assert (data);
    assert (starts);

    assert (npoints);

    for (t = 0; t < NTIMES; t ++) {
        // Generate a random point 
        for (j = 0; j < ndim; j ++) {
            if (plane & (1 << j)) {
                starts [j] = irand (varinfo->dims [j] - counts [j]);
            } else {
                starts [j] = irand (varinfo->dims [j] - 1);
            }
        }

        if (is_aligned == 1) {
            // I need to align every dimension
            for (j = 0; j < ndim; j ++) {
                if (starts [j] % pg_dims [j] == 0) continue;
                starts [j] = (starts [j] / pg_dims [j]) * pg_dims [j];
            }
        } else {
            for (j = 0; j < ndim; j ++) {
                if ((pg_dims [j] > 1) && (starts [j] % pg_dims [j] == 0)) {
                    starts [j] = 1 + irand (pg_dims [j] - 1); 
                }
            }
        }

        #if 1 
            printf ("Starts: ");
            starts [0] = 0;
            starts [1] = 128;
            starts [2] = 28;
            
            for (j = 0; j < ndim; j ++) {
                printf ("%llu ", starts [j]);
            }
            printf ("\n");

            printf ("Counts: ");
            starts [0] = 0;
            for (j = 0; j < ndim; j ++) {
                printf ("%llu ", counts [j]);
            }
            printf ("\n");
        #endif

        sel1 = adios_selection_boundingbox (varinfo->ndim, starts, counts);
        
        adios_schedule_read (f, sel1, varname, 0, 1, data);
        adios_perform_reads (f, BLOCKING);
        
        adios_selection_delete (sel1);
    }

    free (starts);
    free (data);

    adios_free_varinfo (varinfo);
    adios_read_close (f);

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
}

/*
 * Main assumption is that there is going to be only a single timestep and this
 * is for a 3D data. 
 * (1) Read random points
 * (2) Read random XY plane
 * (3) Read random YZ plane
 * (4) Read random XZ plane
 * (5) Read random XYZ subvolume
 */
int main (int argc, char *argv[]) 
{
    MPI_Init (&argc, &argv);
    srand (time (NULL));
    uint64_t counts [] = {1, 10, 10};

    if (argc >= 2) {
        read_bounding_box (argv [1], argv [2], counts, Y + Z, 1);
        // read_points (argv [1], argv [2], 10);
    } else {
        printf ("Usage: <%s> <filename> <variable name>\n", argv [0]);
    }

    MPI_Finalize ();
    return 0;
}
