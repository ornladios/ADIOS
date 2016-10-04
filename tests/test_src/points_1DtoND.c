/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: 
 *  Create 1D points in N-dimensional bounding box containers and test the
 *  function adios_selection_points_1DtoND() to convert them to N-D points
 *
 * How to run: ./points_1DtoND
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include "public/adios.h"
#include "public/adios_read.h"
#include "public/adios_query.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define log(...) fprintf (stderr, "[line %3d]: ", __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);
#define printE(...) fprintf (stderr, "[line %3d]: ERROR: ", __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);


static const int MAX_POINTS = 1000000;


ADIOS_SELECTION *pts;
ADIOS_SELECTION *box;
uint64_t boxstart[3];
uint64_t boxcount[3];

uint64_t *p1D = NULL;      // 1D test point coordinates
uint64_t *pND = NULL;      // converted point coordinates
uint64_t *expected = NULL; // expected values of point coordinates

int test_convert (ADIOS_SELECTION *pts1D, int ndim, int global, uint64_t *expected, enum ADIOS_ERRCODES expected_err)
{
    int nerr = 0;
    int i, d;

    ADIOS_SELECTION *ptsND = adios_selection_points_1DtoND(pts1D, global);

    if (expected_err != adios_errno)
    {
        printE ("    Expected error code %d but got %d: %s\n", expected_err, adios_errno, adios_get_last_errmsg());
        return 1;
    }

    if (expected_err == err_no_error)
    {
        if (ptsND->u.points.ndim != ndim)
        {
            printE ("    Expected %d-dimensional points but got %d-dimensional points\n", ndim, ptsND->u.points.ndim);
            return 1;
        }

        uint64_t *p = ptsND->u.points.points;
        uint64_t *e = expected;

        for (i = 0; i < pts1D->u.points.npoints; i++)
        {
            for (d = 0; d < ndim; d++)
            {
                if (*p != *e)
                {
                    if (nerr < 10) {
                        printE ("    Point %d, dimension %d expected value %" PRIu64 " but got %" PRIu64 "\n", i, d, *e, *p);
                    }
                    nerr++;
                }
                p++;
                e++;
            }
        }
        if (nerr >= 10) {
            printE ("    Too many errors (%d), not printed all\n", nerr);
        }

        if (global)
        {
            // container should be NULL
            if (ptsND->u.points.container_selection != NULL)
            {
                printE ("    Expected NULL container selection when converting points to global offsets but received non-NULL\n");
                nerr++;
            }
        }
        else
        {
            // container should be the same as the input container
            if (ptsND->u.points.container_selection == NULL)
            {
                printE ("    Container selection is NULL when converting points to local offsets\n");
                nerr++;
            }
            else
            {
                for (d = 0; d < ndim; d++)
                {
                    if (pts1D->u.points.container_selection->u.bb.start[d] != ptsND->u.points.container_selection->u.bb.start[d])
                    {
                        printE ("    Result container selection in dimension %d expected value %" PRIu64 " but got %" PRIu64 "\n",
                                d, pts1D->u.points.container_selection->u.bb.start[d], ptsND->u.points.container_selection->u.bb.start[d]);
                        nerr++;
                        break;
                    }
                    if (pts1D->u.points.container_selection->u.bb.count[d] != ptsND->u.points.container_selection->u.bb.count[d])
                    {
                        printE ("    Result container selection in dimension %d expected value %" PRIu64 " but got %" PRIu64 "\n",
                                d, pts1D->u.points.container_selection->u.bb.count[d], ptsND->u.points.container_selection->u.bb.count[d]);
                        nerr++;
                        break;
                    }
                }
            }
        }
    }
    else
    {
        log ("     Expected error detected");
        if (ptsND != NULL)
        {
            printf (" but the results is not NULL");
            nerr++;
        }
        printf ("\n");
    }

    if (ptsND)
        adios_selection_delete (ptsND);
    return nerr;
}

int test ()
{
    int err=0;

    int maxpoints = 10*MAX_POINTS;
    int global;

    while (!p1D || !pND || !expected)
    {
        maxpoints = maxpoints / 10;
        if (!p1D)
            p1D = (uint64_t *) malloc (MAX_POINTS * sizeof(uint64_t));
        if (!pND)
            pND = (uint64_t *) malloc (3 * MAX_POINTS * sizeof(uint64_t));
        if (!expected)
            expected = (uint64_t *) malloc (3 * MAX_POINTS * sizeof(uint64_t));
    }


    /*
     * Tests that should succeed
     */

    // Convert zero 1D point in 2D box
    p1D[0] = 0;
    pts = adios_selection_points(1, 0, p1D);
    boxstart[0] = 3; boxstart[1] = 4;
    boxcount[0] = 5; boxcount[1] = 4;
    box = adios_selection_boundingbox (2, boxstart, boxcount);
    pts->u.points.container_selection = box;
    global = 0;
    log ("  Convert 1D point list of zero points in 2D box to local 2D point\n");
    err = test_convert (pts, 2, global, expected, err_no_error);
    if (err) {
        adios_selection_delete(pts);
        goto endtest;
    }
    global = 1;
    log ("  Convert 1D point list of zero points in 2D box to global 2D point\n");
    err = test_convert (pts, 2, global, expected, err_no_error);
    adios_selection_delete(pts);
    if (err)
        goto endtest;


    // Convert a single 1D point in 2D box
    // box 5x4 with start offset (3,4)
    // local conversion:  point 9 -> (2,1)
    // global conversion:  point 9 -> (5,5)
    p1D[0] = 9;
    pts = adios_selection_points(1, 1, p1D);
    boxstart[0] = 3; boxstart[1] = 4;
    boxcount[0] = 5; boxcount[1] = 4;
    box = adios_selection_boundingbox (2, boxstart, boxcount);
    pts->u.points.container_selection = box;
    global = 0;
    expected[0] = 2;     expected[1] = 1;
    log ("  Convert single 1D point in 2D box to local 2D point\n");
    err = test_convert (pts, 2, global, expected, err_no_error);
    if (err) {
        adios_selection_delete(pts);
        goto endtest;
    }
    global = 1;
    expected[0] = 5;     expected[1] = 5;
    log ("  Convert single 1D point in 2D box to global 2D point\n");
    err = test_convert (pts, 2, global, expected, err_no_error);
    adios_selection_delete(pts);
    if (err)
        goto endtest;


    // Convert a two 1D points in 2D box
    // box 5x4 with start offset (3,4)
    // local conversion:  point 7 -> (1,3), 8 -> (2,0)
    // global conversion: point 7 -> (4,7), 8 -> (5,4)
    p1D[0] = 7;
    p1D[1] = 8;
    pts = adios_selection_points(1, 1, p1D);
    boxstart[0] = 3; boxstart[1] = 4;
    boxcount[0] = 5; boxcount[1] = 4;
    box = adios_selection_boundingbox (2, boxstart, boxcount);
    pts->u.points.container_selection = box;
    global = 0;
    expected[0] = 1;     expected[1] = 3;
    expected[2] = 2;     expected[3] = 0;
    log ("  Convert two 1D points in 2D box to local 2D point\n");
    err = test_convert (pts, 2, global, expected, err_no_error);
    if (err) {
        adios_selection_delete(pts);
        goto endtest;
    }
    global = 1;
    expected[0] = 4;     expected[1] = 7;
    expected[2] = 5;     expected[3] = 4;
    log ("  Convert two 1D points in 2D box to global 2D point\n");
    err = test_convert (pts, 2, global, expected, err_no_error);
    adios_selection_delete(pts);
    if (err)
        goto endtest;



    // Convert a single 1D point in 3D box
    // box 5x4x2, with start offset (3,4,7)
    // local conversion:  point 9 -> (1,0,1)
    // global conversion:  point 9 -> (4,4,8)
    p1D[0] = 9;
    pts = adios_selection_points(1, 1, p1D);
    boxstart[0] = 3; boxstart[1] = 4; boxstart[2] = 7;
    boxcount[0] = 5; boxcount[1] = 4; boxcount[2] = 2;
    box = adios_selection_boundingbox (3, boxstart, boxcount);
    pts->u.points.container_selection = box;
    global = 0;
    expected[0] = 1;     expected[1] = 0;    expected[2] = 1;
    log ("  Convert single 1D point in 3D box to local 3D point\n");
    err = test_convert (pts, 3, global, expected, err_no_error);
    if (err) {
        adios_selection_delete(pts);
        goto endtest;
    }
    global = 1;
    expected[0] = 4;     expected[1] = 4;    expected[2] = 8;
    log ("  Convert single 1D point in 3D box to global 3D point\n");
    err = test_convert (pts, 3, global, expected, err_no_error);
    adios_selection_delete(pts);
    if (err)
        goto endtest;


    // Convert two 1D points in 3D box
    // box 5x4x2, with start offset (3,4,7)
    // local conversion:  point 9 -> (1,0,1), 20 -> (2,2,0)
    // global conversion:  point 9 -> (4,4,8), 20 -> (5,6,7)
    p1D[0] = 9;
    p1D[1] = 20;
    pts = adios_selection_points(1, 1, p1D);
    boxstart[0] = 3; boxstart[1] = 4; boxstart[2] = 7;
    boxcount[0] = 5; boxcount[1] = 4; boxcount[2] = 2;
    box = adios_selection_boundingbox (3, boxstart, boxcount);
    pts->u.points.container_selection = box;
    global = 0;
    expected[0] = 1;     expected[1] = 0;    expected[2] = 1;
    expected[3] = 2;     expected[4] = 2;    expected[5] = 0;
    log ("  Convert two 1D points in 3D box to local 3D point\n");
    err = test_convert (pts, 3, global, expected, err_no_error);
    if (err) {
        adios_selection_delete(pts);
        goto endtest;
    }
    global = 1;
    expected[0] = 4;     expected[1] = 4;    expected[2] = 8;
    expected[3] = 5;     expected[4] = 6;    expected[5] = 7;
    log ("  Convert two 1D points in 3D box to global 3D point\n");
    err = test_convert (pts, 3, global, expected, err_no_error);
    adios_selection_delete(pts);
    if (err)
        goto endtest;


    // Convert all 1D points in 3D box
    int d3 = maxpoints / (70*90);
    int d2 = maxpoints / d3 / 90;
    int d1 = maxpoints / (d2*d3);
    int n = d1*d2*d3;
    int i,j,k;

    for (i = 0; i < n; ++i) {
        p1D[i] = i;
    }

    pts = adios_selection_points(1, n, p1D);
    boxstart[0] = 100; boxstart[1] = 200; boxstart[2] = 300;
    boxcount[0] = d1; boxcount[1] = d2; boxcount[2] = d3;
    box = adios_selection_boundingbox (3, boxstart, boxcount);
    pts->u.points.container_selection = box;
    global = 0;
    uint64_t *e = expected;
    for (i = 0; i < d1; ++i) {
        for (j = 0; j < d2; ++j) {
            for (k = 0; k < d3; ++k) {
                e[0] = i;
                e[1] = j;
                e[2] = k;
                e = e + 3;
            }
        }
    }
    log ("  Convert all %d 1D points in a %dx%dx%d box to local 3D points\n", n, d1,d2,d3);
    err = test_convert (pts, 3, global, expected, err_no_error);
    if (err) {
        adios_selection_delete(pts);
        goto endtest;
    }
    global = 1;
    e = expected;
    for (i = 0; i < d1; ++i) {
        for (j = 0; j < d2; ++j) {
            for (k = 0; k < d3; ++k) {
                e[0] = e[0] + boxstart[0];
                e[1] = e[1] + boxstart[1];
                e[2] = e[2] + boxstart[2];
                e = e + 3;
            }
        }
    }
    log ("  Convert all %d 1D points in a %dx%dx%d box to global 3D points\n", n, d1,d2,d3);
    err = test_convert (pts, 3, global, expected, err_no_error);
    adios_selection_delete(pts);
    if (err)
        goto endtest;



    /*
     * Tests that should fail
     */

    // Convert a single point with NULL container selection
    // This function should return NULL and print an error
    p1D[0] = 0;
    pts = adios_selection_points(1, 1, p1D);
    pts->u.points.container_selection = NULL;
    global = 0;
    expected[0] = 0;
    log ("  Convert single 1D point without container selection, expect error\n");
    err = test_convert (pts, 2, global, expected, err_invalid_selection);
    adios_selection_delete(pts);
    if (err)
        goto endtest;


    // Convert a single 2D point with NULL container selection
    // This function should return NULL and print an error
    p1D[0] = 0; p1D[1] = 0;
    pts = adios_selection_points(2, 1, p1D);
    boxstart[0] = 3; boxstart[1] = 4;
    boxcount[0] = 5; boxcount[1] = 4;
    box = adios_selection_boundingbox (2, boxstart, boxcount);
    pts->u.points.container_selection = box;
    global = 0;
    expected[0] = 0;
    log ("  Convert single 2D point in 2D box, expect error\n");
    err = test_convert (pts, 2, global, expected, err_invalid_selection);
    adios_selection_delete(pts);
    if (err)
        goto endtest;

    // Convert a single point in writeblock 2
    // This function should return NULL and print an error
    ADIOS_SELECTION *wblock = adios_selection_writeblock(2);

    p1D[0] = 0;
    pts = adios_selection_points(1, 1, p1D);
    pts->u.points.container_selection = wblock;
    global = 0;
    expected[0] = 0;
    log ("  Convert single 1D point in a writeblock, expect error\n");
    err = test_convert (pts, 2, global, expected, err_invalid_selection);
    adios_selection_delete(pts);
    //adios_selection_delete (wblock); // deleted when pts is deleted
    if (err)
        goto endtest;


endtest:
    return err;
}

int main (int argc, char ** argv)
{
    int err = test();
    if (!err) {
        log ("All tests succeeded\n")
    } else {
        log ("A test failed, consecutive tests were not executed\n")
    }
    return err;
}
