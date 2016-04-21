/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: 
 *  Write a 2D array of 4 2D blocks, multiple steps. Each block is a 5x5 array.
 *  The whole array is 10x10 and is patterned like this: 10*step+<blockid>.<row><col>
 *
$ bpls -l read_points_2d.bp -d data -n 10 -f "%5.2f "
  real     data                       2*{10, 10} = 0 / 13.44 / 6.72 / 5.12545
    (0,0,0)     0.00   0.01   0.02   0.03   0.04  |  1.00   1.01   1.02   1.03   1.04
    (0,1,0)     0.10   0.11   0.12   0.13   0.14  |  1.10   1.11   1.12   1.13   1.14
    (0,2,0)     0.20   0.21   0.22   0.23   0.24  |  1.20   1.21   1.22   1.23   1.24
    (0,3,0)     0.30   0.31   0.32   0.33   0.34  |  1.30   1.31   1.32   1.33   1.34
    (0,4,0)     0.40   0.41   0.42   0.43   0.44  |  1.40   1.41   1.42   1.43   1.44
    (0,5,0)     2.00   2.01   2.02   2.03   2.04  |  3.00   3.01   3.02   3.03   3.04
    (0,6,0)     2.10   2.11   2.12   2.13   2.14  |  3.10   3.11   3.12   3.13   3.14
    (0,7,0)     2.20   2.21   2.22   2.23   2.24  |  3.20   3.21   3.22   3.23   3.24
    (0,8,0)     2.30   2.31   2.32   2.33   2.34  |  3.30   3.31   3.32   3.33   3.34
    (0,9,0)     2.40   2.41   2.42   2.43   2.44  |  3.40   3.41   3.42   3.43   3.44
    ----------------------------------------------+----------------------------------
    (1,0,0)    10.00  10.01  10.02  10.03  10.04  | 11.00  11.01  11.02  11.03  11.04
    (1,1,0)    10.10  10.11  10.12  10.13  10.14  | 11.10  11.11  11.12  11.13  11.14
    (1,2,0)    10.20  10.21  10.22  10.23  10.24  | 11.20  11.21  11.22  11.23  11.24
    (1,3,0)    10.30  10.31  10.32  10.33  10.34  | 11.30  11.31  11.32  11.33  11.34
    (1,4,0)    10.40  10.41  10.42  10.43  10.44  | 11.40  11.41  11.42  11.43  11.44
    (1,5,0)    12.00  12.01  12.02  12.03  12.04  | 13.00  13.01  13.02  13.03  13.04
    (1,6,0)    12.10  12.11  12.12  12.13  12.14  | 13.10  13.11  13.12  13.13  13.14
    (1,7,0)    12.20  12.21  12.22  12.23  12.24  | 13.20  13.21  13.22  13.23  13.24
    (1,8,0)    12.30  12.31  12.32  12.33  12.34  | 13.30  13.31  13.32  13.33  13.34
    (1,9,0)    12.40  12.41  12.42  12.43  12.44  | 13.40  13.41  13.42  13.43  13.44

 *
 *  Then test every possible reading of points
 *    with a bounding box container, and without
 *    with a writeblock container
 *    1D and 2D points (local offset in contiguous space vs N-dim points)
 *    multiple timesteps
 *
 * How to run: ./read_points <N> <steps>
 * It writes N*N 2D blocks organized into a NxN 2D array.
 * Output: read_points_2d.bp
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

#define log(...) fprintf (stderr, "[rank=%3.3d, line %d]: ", rank, __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);
#define printE(...) fprintf (stderr, "[rank=%3.3d, line %d]: ERROR: ", rank, __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);

/* user arguments */
int N = 1;       // organize blocks in NxN shape
int NSTEPS = 1;  // number of output steps

static const char FILENAME[] = "read_points_2d.bp";
#define VALUE(rank, step) (step * 1000 + rank + 1)

static const int ldim1 = 5;
static const int ldim2 = 5;
int gdim1, gdim2;
int offs1, offs2;

int64_t       m_adios_group;

/* Variables to write */
float  a2[ldim1*ldim2];

/* Variables to read */
float  r2[ldim1*ldim2];

MPI_Comm    comm = MPI_COMM_SELF; // dummy comm for sequential code
int rank;
int size;


void set_gdim()
{
    gdim1 = N*ldim1;
    gdim2 = N*ldim2;
}

void set_offsets (int row, int col)
{
	offs1 = row*ldim1;
	offs2 = col*ldim2;
}

void fill_block(int step, int row, int col)
{
	int n;
    float v_intpart = 10*step + row*N + col;
    float v;
    int i, j, k;

    n = ldim1 * ldim2;
    //log ("  Fill up array of %d elements starting from value %f...\n",n,v_start);
    k = 0;
    for (i=0; i<ldim1; i++) {
    	v = v_intpart + i*0.1;
        //log ("      row %d starts from value %f... (element %d)\n",i,v,k);
    	for (j=0; j<ldim2; j++) {
    		a2[k] = v;
    		k++;
    		v += 0.01;
    	}
    }
}


void Usage() 
{
    printf("Usage: read_points <N> <nsteps>\n"
            "    <N>:       Number of blocks in each of X and Y direction\n"
    		"    <nsteps>:  Number of write cycles (to same file)\n");
}

void define_vars ();
int write_file (int step);
int read_points ();

int main (int argc, char ** argv) 
{
    int err,i ; 

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    if (argc == 1)
    {
        // this case is for the test harness. otherwise this should be calling for Usage();
        N = 2;
        NSTEPS = 2;
        printf("Running read_points <N=%d> <nsteps=%d>\n", N, NSTEPS);
    }
    else
    {
        if (argc < 3) { Usage(); return 1; }

        errno = 0;
        i = strtol (argv[1], NULL, 10);
        if (errno || i < 1) { printf("Invalid 1st argument %s\n", argv[1]); Usage(); return 1;}
        N = i;

        errno = 0;
        i = strtol (argv[2], NULL, 10);
        if (errno || i < 1) { printf("Invalid 2nd argument %s\n", argv[2]); Usage(); return 1;}
        NSTEPS = i;
    }
    adios_init_noxml (comm);
    err = adios_read_init_method(ADIOS_READ_METHOD_BP, comm, "verbose=2");
    if (err) {
        printE ("%s\n", adios_errmsg());
    }

    adios_declare_group (&m_adios_group, "read_points", "", adios_flag_yes);
    adios_select_method (m_adios_group, "POSIX", "", "");


    define_vars();
    set_gdim();
    
    for (i=0; i<NSTEPS; i++) {
        if (!err) {
            err = write_file (i); 
        }
    }

    if (!err)
        err = read_points ();

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    adios_finalize (rank);
    MPI_Finalize ();
    return err;
}

void define_vars ()
{
    int i;

    adios_define_var (m_adios_group, "ldim1", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "ldim2", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "gdim1", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "gdim2", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "offs1", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "offs2", "", adios_integer, 0, 0, 0);

    for (i=0; i<N*N; i++) {
    	adios_define_var (m_adios_group, "data", "", adios_real,
    			"ldim1,ldim2",
				"gdim1,gdim2",
				"offs1,offs2");
    }
}

int write_file (int step) 
{
    int64_t       fh;
    int           i, j;
    double        tb, te;

    log ("Write step %d to %s\n", step, FILENAME);
    adios_open (&fh, "read_points", FILENAME, (step ? "a" : "w"), comm);
    
    /*
    groupsize  = (4 + nblocks*2) * sizeof(int);             // dimensions
    log ("  groupsize calculated = %llu\n", groupsize);
    groupsize += nblocks * ldim1 * ldim2 * sizeof(float);     // 2D  blocks
    log ("  groupsize calculated = %llu\n", groupsize);

    adios_group_size (fh, groupsize, &totalsize);
    log ("  groupsize %llu, totalsize %llu\n", groupsize, totalsize);
    */

    tb = MPI_Wtime();
    for (i=0; i<N; i++) {
    	for (j=0; j<N; j++) {
    		set_offsets (i, j);
    		fill_block (step, i, j);
    		adios_write (fh, "gdim1", &gdim1);
    		adios_write (fh, "gdim2", &gdim2);
    		adios_write (fh, "ldim1", &ldim1);
    		adios_write (fh, "ldim2", &ldim2);
    		adios_write (fh, "offs1", &offs1);
    		adios_write (fh, "offs2", &offs2);
    		/*
    		int k=0, l, m;
    		for (l=0; l<ldim1; l++) {
    	        printf ("      a[%d,*] = ",l);
    	    	for (m=0; m<ldim2; m++) {
        	        printf ("%2.0f ",a2[k]);
    	    		k++;
    	    	}
    	    	printf ("\n");
    	    }
	    	printf ("\n");
	    	*/
    		adios_write (fh, "data", a2);
     	}
    }
    adios_close (fh);
    te = MPI_Wtime();

    if (rank==0) {
        log ("  Write time for step %d was %6.3lf seconds\n", step, te-tb);
    }
    MPI_Barrier (comm);
    return 0;
}

#define CHECK_VARINFO(VARNAME, NDIM, NSTEPS) \
    vi = adios_inq_var (f, VARNAME); \
    if (vi == NULL) { \
        printE ("No such variable: %s\n", VARNAME); \
        err = 101; \
        goto endread; \
    } \
    if (vi->ndim != NDIM) { \
        printE ("Variable %s has %d dimensions, but expected %d\n", VARNAME, vi->ndim, NDIM); \
        err = 102; \
        goto endread; \
    } \
    if (vi->nsteps != NSTEPS) { \
        printE ("Variable %s has %d steps, but expected %d\n", VARNAME, vi->nsteps, NSTEPS); \
        err = 103; \
        /*goto endread; */\
    } \
    adios_free_varinfo (vi);


void reset_readvars()
{
    size_t n;
    n = (size_t)ldim1 * (size_t)ldim2;
    memset (r2,  0, n*sizeof(float));
}

int nearlyEqual(float a, float b, float epsilon) {
    float absA = fabsf(a);
    float absB = fabsf(b);
    float diff = fabsf(a - b);

    if (a == b) { // shortcut, handles infinities
        return 1;
    } else if (a == 0 || b == 0 || diff < FLT_MIN) {
        // a or b is zero or both are extremely close to it
        // relative error is less meaningful here
        return diff < (epsilon * FLT_MIN);
    } else { // use relative error
        return diff / fmin((absA + absB), FLT_MAX) < epsilon;
    }
}

int test_read (ADIOS_FILE *f, ADIOS_SELECTION *pts, int from_steps, int nsteps,
        float *expected)
{
    float *data = (float *) malloc (nsteps * pts->u.points.npoints * sizeof(float));
    adios_schedule_read(f, pts, "data", from_steps, nsteps, data);
    adios_perform_reads(f, 1);
    int s, i, idx;

    // print it out just for manual testing's sake
    for (s = 0; s < nsteps; ++s) {
        log ("    Step %2d: ", s);
        idx = s * pts->u.points.npoints;
        for (i = 0; i < pts->u.points.npoints; ++i) {
            if (nearlyEqual(data[idx], expected[idx], 0.000001)) {
                printf (" %6.2f  ", data[idx]);
            } else {
                printf ("**%5.2f* ", data[idx]);
              }
            idx++;
        }
        printf("\n");
    }

    for (s = 0; s < nsteps; ++s) {
        idx = s * pts->u.points.npoints;
        for (i = 0; i < pts->u.points.npoints; ++i) {
            if (!nearlyEqual(data[idx], expected[idx], 0.000001))
            {
//                printE ("Point %d in step %d value = %5.2f but was expected %5.2f\n",
                printE ("Point %d in step %d value = %16.12f but was expected %16.12f\n",
                        i, s, data[idx], expected[idx]);
                free (data);
                return 110;
            }
            idx++;
        }
    }


    free (data);
    return 0;
}

int read_points ()
{
    ADIOS_FILE * f;
    ADIOS_VARINFO * vi;
    int err=0;

    reset_readvars();

    log ("Open %s for reading\n", FILENAME);
    f = adios_read_open_file (FILENAME, ADIOS_READ_METHOD_BP, comm);
    if (f == NULL) {
        printE ("Error at opening file: %s\n", adios_errmsg());
        return 1;
    }

    log ("  Check variable definitions in %s\n", FILENAME);
    CHECK_VARINFO("data", 2, NSTEPS)
    MPI_Barrier (comm);

    ADIOS_SELECTION *pts;
    uint64_t start[100];
    uint64_t count[100];
    ADIOS_SELECTION *box;
    uint64_t boxstart[2];
    uint64_t boxcount[2];
    ADIOS_SELECTION *wblock;
    float expected[32]; // expected values of points

    /*
     * Points without containers based tests
     */

    // Test 1
    // Read a single point with 2D coordinates in global space
    // middle point of first 5x5 block, = 0.22
    start[0] = 2; start[1] = 2;
    pts = adios_selection_points(2, 1, start);
    expected[0] = 0.22;
    log ("  Read single 2D global point at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;

    /*
     * Points in BOUNDINGBOX based tests
     */
    // Test 2
    // Read a single point with 2D coordinates in an 5x5 bounding box
    // middle point of that 5x5 block, = 3.00
    // Limit the query to the bounding box of the middle of the global array
    boxstart[0] = 3; boxstart[1] = 3;
    boxcount[0] = 5; boxcount[1] = 5;
    box = adios_selection_boundingbox (2, boxstart, boxcount);

    start[0] = 2; start[1] = 2; // 2D point of middle of a 5x5 array
    pts = adios_selection_points(2, 1, start);
    pts->u.points.container_selection = box;
    expected[0] = 3.00;
    log ("  Read single 2D point in a 2D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    //adios_selection_delete (box); // deleted when pts is deleted
    if (err)
        goto endread;


    // Test 3
    // Read a single point with 1D offset in an 5x5 bounding box
    // middle point of that 5x5 block, = 3.00
    // Limit the query to the bounding box of the middle of the global array
    // This is the same actual point as in Test 2
    boxstart[0] = 3; boxstart[1] = 3;
    boxcount[0] = 5; boxcount[1] = 5;
    box = adios_selection_boundingbox (2, boxstart, boxcount);

    start[0] = 12; // offset of middle of a 5x5 array
    pts = adios_selection_points(1, 1, start);
    pts->u.points.container_selection = box;
    expected[0] = 3.00;
    log ("  Read single 1D point in a 2D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;

    // Test 4-5
    // Read several single points with 2D coordinates in an 5x5 bounding box
    // back diagonal of that 5x5 block, = 1.32  1.41  3.00  2.14  2.23
    // Limit the query to the bounding box of the middle of the global array
    boxstart[0] = 3; boxstart[1] = 3;
    boxcount[0] = 5; boxcount[1] = 5;
    box = adios_selection_boundingbox (2, boxstart, boxcount);

    // back diagonal of the 5x5 array
    start[0] = 0; start[1] = 4;
    start[2] = 1; start[3] = 3;
    start[4] = 2; start[5] = 2;
    start[6] = 3; start[7] = 1;
    start[8] = 4; start[9] = 0;
    pts = adios_selection_points(2, 5, start);
    pts->u.points.container_selection = box;
    expected[0] = 1.32;
    expected[1] = 1.41;
    expected[2] = 3.00;
    expected[3] = 2.14;
    expected[4] = 2.23;
    log ("  Read back diagonal of box with five 2D points in a 2D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    if (err) {
        adios_selection_delete(pts);
        goto endread;
    }
    expected[5] = 11.32;
    expected[6] = 11.41;
    expected[7] = 13.00;
    expected[8] = 12.14;
    expected[9] = 12.23;
    log ("  Read back diagonal of box with five 2D points in a 2D bounding box at two steps\n");
    err = test_read (f, pts, 0, 2, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;

    // Test 6-7
    // Read several single points with 1D offset in an 5x5 bounding box
    // back diagonal of that 5x5 block, = 1.32  1.41  3.00  2.14  2.23
    // Limit the query to the bounding box of the middle of the global array
    // This is the same actual point as in Test 4
    boxstart[0] = 3; boxstart[1] = 3;
    boxcount[0] = 5; boxcount[1] = 5;
    box = adios_selection_boundingbox (2, boxstart, boxcount);

    // back diagonal of the 5x5 array as 1D offsets
    start[0] = 4;
    start[1] = 8;
    start[2] = 12;
    start[3] = 16;
    start[4] = 20;
    pts = adios_selection_points(1, 5, start);
    pts->u.points.container_selection = box;
    expected[0] = 1.32;
    expected[1] = 1.41;
    expected[2] = 3.00;
    expected[3] = 2.14;
    expected[4] = 2.23;
    log ("  Read back diagonal of box with five 1D offsets in a 2D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    if (err) {
        adios_selection_delete(pts);
        goto endread;
    }
    expected[5] = 11.32;
    expected[6] = 11.41;
    expected[7] = 13.00;
    expected[8] = 12.14;
    expected[9] = 12.23;
    log ("  Read back diagonal of box with five 1D points in a 2D bounding box at two steps\n");
    err = test_read (f, pts, 0, 2, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;



    /*
     * Points in WRITEBLOCK based tests
     */


    // Test 8
    // Read a single point with 2D coordinates in writeblock 2 (third block)
    // middle point of that 5x5 block, = 2.22
    wblock = adios_selection_writeblock(2);

    start[0] = 2; start[1] = 2; // 2D point of middle of a 5x5 array
    pts = adios_selection_points(2, 1, start);
    pts->u.points.container_selection = wblock;
    expected[0] = 2.22;
    log ("  Read single 2D point in a 2D writeblock at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    //adios_selection_delete (wblock); // deleted when pts is deleted
    if (err)
        goto endread;


    // Test 9
    // Read a single point with 1D offset in an writeblock 2
    // middle point of that 5x5 block, = 2.22
    // This is the same actual point as in Test 8
    wblock = adios_selection_writeblock(2);

    start[0] = 12; // offset of middle of a 5x5 array
    pts = adios_selection_points(1, 1, start);
    pts->u.points.container_selection = wblock;
    expected[0] = 2.22;
    log ("  Read single 1D point in a 2D writeblock at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;

    // Test 10-11
    // Read several single points with 2D coordinates in writeblock 2
    // back diagonal of that 5x5 block, = 2.04 2.13 2.22 2.31 2.40
    wblock = adios_selection_writeblock(2);

    // back diagonal of the 5x5 array
    start[0] = 0; start[1] = 4;
    start[2] = 1; start[3] = 3;
    start[4] = 2; start[5] = 2;
    start[6] = 3; start[7] = 1;
    start[8] = 4; start[9] = 0;
    pts = adios_selection_points(2, 5, start);
    pts->u.points.container_selection = wblock;
    expected[0] = 2.04;
    expected[1] = 2.13;
    expected[2] = 2.22;
    expected[3] = 2.31;
    expected[4] = 2.40;
    log ("  Read back diagonal of box with five 2D points in a 2D writeblock at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    if (err) {
        adios_selection_delete(pts);
        goto endread;
    }

    expected[5] = 12.04;
    expected[6] = 12.13;
    expected[7] = 12.22;
    expected[8] = 12.31;
    expected[9] = 12.40;
    log ("  Read back diagonal of box with five 2D points in a 2D writeblock at two steps\n");
    err = test_read (f, pts, 0, 2, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;

    // Test 12-13
    // Read several single points with 1D offset in writeblock 2
    // back diagonal of that 5x5 block, = 2.04 2.13 2.22 2.31 2.40
    // This is the same actual point as in Test 10-11
    wblock = adios_selection_writeblock(2);

    // back diagonal of the 5x5 array as 1D offsets
    start[0] = 4;
    start[1] = 8;
    start[2] = 12;
    start[3] = 16;
    start[4] = 20;
    pts = adios_selection_points(1, 5, start);
    pts->u.points.container_selection = wblock;
    expected[0] = 2.04;
    expected[1] = 2.13;
    expected[2] = 2.22;
    expected[3] = 2.31;
    expected[4] = 2.40;
    log ("  Read back diagonal of box with five 1D offsets in a 2D writeblock at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    if (err) {
        adios_selection_delete(pts);
        goto endread;
    }
    expected[5] = 12.04;
    expected[6] = 12.13;
    expected[7] = 12.22;
    expected[8] = 12.31;
    expected[9] = 12.40;
    log ("  Read back diagonal of box with five 1D points in a 2D writeblock at two steps\n");
    err = test_read (f, pts, 0, 2, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;

endread:

    adios_read_close(f);
    MPI_Barrier (comm);
    return err;
}

