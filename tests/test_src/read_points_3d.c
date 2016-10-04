/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: 
 *  Write a 3D array of 8 3D blocks, multiple steps. Each block is a 5x5x5 array.
 *  The whole array is 10x10x10 and is patterned like this: 10*step+<blockid>.<row><col><z>
 *
 *  Front 2D slice of the 3D array:
 *  $ bpls -l read_points_3d.bp -d data -n 10 -f "%6.3f" -s "0,0,0,0" -c "1,10,10,1"
  real     data                       2*{10, 10, 10} = 0 / 17.444 / 8.722 / 5.50184
    slice (0:0, 0:9, 0:9, 0:0)
    (0,0,0,0)     0.000  0.010  0.020  0.030  0.040  1.000  1.010  1.020  1.030  1.040
    (0,1,0,0)     0.100  0.110  0.120  0.130  0.140  1.100  1.110  1.120  1.130  1.140
    (0,2,0,0)     0.200  0.210  0.220  0.230  0.240  1.200  1.210  1.220  1.230  1.240
    (0,3,0,0)     0.300  0.310  0.320  0.330  0.340  1.300  1.310  1.320  1.330  1.340
    (0,4,0,0)     0.400  0.410  0.420  0.430  0.440  1.400  1.410  1.420  1.430  1.440
    (0,5,0,0)     2.000  2.010  2.020  2.030  2.040  3.000  3.010  3.020  3.030  3.040
    (0,6,0,0)     2.100  2.110  2.120  2.130  2.140  3.100  3.110  3.120  3.130  3.140
    (0,7,0,0)     2.200  2.210  2.220  2.230  2.240  3.200  3.210  3.220  3.230  3.240
    (0,8,0,0)     2.300  2.310  2.320  2.330  2.340  3.300  3.310  3.320  3.330  3.340
    (0,9,0,0)     2.400  2.410  2.420  2.430  2.440  3.400  3.410  3.420  3.430  3.440
 *
 *  Back 2D slice of the 3D array
 *  $ bpls -l read_points_3d.bp -d data -n 10 -f "%6.3f" -s "0,0,0,-1" -c "1,10,10,1"
  real     data                       2*{10, 10, 10} = 0 / 17.444 / 8.722 / 5.50184
    slice (0:0, 0:9, 0:9, 9:9)
    (0,0,0,9)     4.004  4.014  4.024  4.034  4.044  5.004  5.014  5.024  5.034  5.044
    (0,1,0,9)     4.104  4.114  4.124  4.134  4.144  5.104  5.114  5.124  5.134  5.144
    (0,2,0,9)     4.204  4.214  4.224  4.234  4.244  5.204  5.214  5.224  5.234  5.244
    (0,3,0,9)     4.304  4.314  4.324  4.334  4.344  5.304  5.314  5.324  5.334  5.344
    (0,4,0,9)     4.404  4.414  4.424  4.434  4.444  5.404  5.414  5.424  5.434  5.444
    (0,5,0,9)     6.004  6.014  6.024  6.034  6.044  7.004  7.014  7.024  7.034  7.044
    (0,6,0,9)     6.104  6.114  6.124  6.134  6.144  7.104  7.114  7.124  7.134  7.144
    (0,7,0,9)     6.204  6.214  6.224  6.234  6.244  7.204  7.214  7.224  7.234  7.244
    (0,8,0,9)     6.304  6.314  6.324  6.334  6.344  7.304  7.314  7.324  7.334  7.344
    (0,9,0,9)     6.404  6.414  6.424  6.434  6.444  7.404  7.414  7.424  7.434  7.444

 * Next slice in Z-coord is similar with all values +0.001
 *  Then test every possible reading of points
 *    with a bounding box container, and without
 *    with a writeblock container
 *    1D and 3D points (local offset in contiguous space vs N-dim points)
 *    multiple timesteps
 *
 * How to run: ./read_points_3d
 * It writes 5x5x5 3D blocks organized into a 10x10x10 3D array.
 * Output: read_points.bp
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
int N = 2;       // organize blocks in NxNxN shape
int NSTEPS = 2;  // number of output steps

static const char FILENAME[] = "read_points_3d.bp";
#define VALUE(rank, step) (step * 1000 + rank + 1)

#define LDIM1 5
#define LDIM2 5
#define LDIM3 5

static const int ldim1 = LDIM1;
static const int ldim2 = LDIM2;
static const int ldim3 = LDIM3;

int gdim1, gdim2, gdim3;
int offs1, offs2, offs3;

int64_t       m_adios_group;

/* Variables to write */
float  a3[LDIM1*LDIM2*LDIM3];

/* Variables to read */
float  r3[LDIM1*LDIM2*LDIM3];

MPI_Comm    comm = MPI_COMM_SELF; // dummy comm for sequential code
int rank;
int size;


void set_gdim()
{
    gdim1 = N*ldim1;
    gdim2 = N*ldim2;
    gdim3 = N*ldim3;
}

void set_offsets (int row, int col, int z)
{
    offs1 = row*ldim1;
    offs2 = col*ldim2;
    offs3 = z*ldim3;
}

void fill_block(int step, int row, int col, int z)
{
    int n;
    float v_intpart = 10*step + row*N + col + z*N*N;
    float v;
    int i, j, k, idx;

    n = ldim1 * ldim2 * ldim3;
    //log ("  Fill up array of %d elements starting from value %f...\n",n,v_start);
    idx = 0;
    for (i=0; i<ldim1; i++) {
        //v = v_intpart + i*0.1;
        //log ("      row %d starts from value %f... (element %d)\n",i,v,k);
        for (j=0; j<ldim2; j++) {
            v = v_intpart + i*0.1 + j*0.01;
            for (k=0; k<ldim3; k++) {
                a3[idx] = v;
                idx++;
                v += 0.001;
            }
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

    adios_declare_group (&m_adios_group, "read_points", "", adios_stat_default);
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
    adios_define_var (m_adios_group, "ldim3", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "gdim1", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "gdim2", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "gdim3", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "offs1", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "offs2", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "offs3", "", adios_integer, 0, 0, 0);

    for (i=0; i<N*N*N; i++) {
        adios_define_var (m_adios_group, "data", "", adios_real,
                "ldim1,ldim2,ldim3",
                "gdim1,gdim2,gdim3",
                "offs1,offs2,offs3");
    }
}

int write_file (int step) 
{
    int64_t       fh;
    int           i, j, k;
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
            for (k=0; k<N; k++) {
                set_offsets (i, j, k);
                fill_block (step, i, j, k);
                adios_write (fh, "gdim1", &gdim1);
                adios_write (fh, "gdim2", &gdim2);
                adios_write (fh, "gdim3", &gdim3);
                adios_write (fh, "ldim1", &ldim1);
                adios_write (fh, "ldim2", &ldim2);
                adios_write (fh, "ldim3", &ldim3);
                adios_write (fh, "offs1", &offs1);
                adios_write (fh, "offs2", &offs2);
                adios_write (fh, "offs3", &offs3);
                adios_write (fh, "data", a3);
            }
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
    memset (r3,  0, n*sizeof(float));
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
                printf (" %7.3f  ", data[idx]);
            } else {
                printf ("**%6.3f* ", data[idx]);
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
    CHECK_VARINFO("data", 3, NSTEPS)
    MPI_Barrier (comm);

    ADIOS_SELECTION *pts;
    uint64_t start[100];
    uint64_t count[100];
    ADIOS_SELECTION *box;
    uint64_t boxstart[3];
    uint64_t boxcount[3];
    ADIOS_SELECTION *wblock;
    float expected[32]; // expected values of points

    /*
     * Points without containers based tests
     */

    // Test 1
    // Read a single point with 3D coordinates in global space
    // middle point of first 5x5x5 block, = 0.22
    start[0] = 2; start[1] = 2; start[2] = 2;
    pts = adios_selection_points(3, 1, start);
    expected[0] = 0.222;
    log ("  Read single 3D global point at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;

    /*
     * Points in BOUNDINGBOX based tests
     */
    log ("  --------------- Points in Boundinx Boxes -----------------------  \n");
    // Test 2
    // Read a single point with 3D coordinates in an 5x5x5 bounding box
    // middle point of the 5x5x5 block, = 7.00
    // Limit the query to the bounding box
    boxstart[0] = 3; boxstart[1] = 3; boxstart[2] = 3;
    boxcount[0] = 5; boxcount[1] = 5; boxcount[2] = 5;
    box = adios_selection_boundingbox (3, boxstart, boxcount);

    start[0] = 2; start[1] = 2; start[2] = 2; // 3D point of middle of a 5x5x5 array
    pts = adios_selection_points(3, 1, start);
    pts->u.points.container_selection = box;
    expected[0] = 7.000;
    log ("  Read single 3D point in a 3D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    //adios_selection_delete (box); // deleted when pts is deleted
    if (err)
        goto endread;


    // Test 3
    // Read a single point with 1D offset in an 5x5x5 bounding box
    // middle point of that 5x5x5 block, = 7.00
    // Limit the query to the bounding box of the middle of the global array
    // This is the same actual point as in Test 2
    boxstart[0] = 3; boxstart[1] = 3; boxstart[2] = 3;
    boxcount[0] = 5; boxcount[1] = 5; boxcount[2] = 5;
    box = adios_selection_boundingbox (3, boxstart, boxcount);

    start[0] = 62; // offset of middle of a 5x5x5 array: 2*25 + 2*5 + 2
    pts = adios_selection_points(1, 1, start);
    pts->u.points.container_selection = box;
    expected[0] = 7.000;
    log ("  Read single 1D point in a 3D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;


    // Test 4-5
    // Read several single points with 3D coordinates in an 5x5x5 bounding box
    // back diagonal of that 5x5x5 block, = 5.332  5.411  7.000  2.144  2.233
    // Limit the query to the bounding box of the middle of the global array
    /* See: $ bpls -l read_points_3d.bp -d data -n 5 -s "0,3,3,3" -c "1,5,5,5" -f "%6.3f "
      real     data     2*{10, 10, 10} = 0 / 17.444 / 8.722 / 5.50184
      slice (0:0, 3:7, 3:7, 3:7)
      (0,3,3,3)     0.333   0.334   4.330   4.331   4.332
      (0,3,4,3)     0.343   0.344   4.340   4.341   4.342
      (0,3,5,3)     1.303   1.304   5.300   5.301   5.302
      (0,3,6,3)     1.313   1.314   5.310   5.311   5.312
      (0,3,7,3)     1.323   1.324   5.320   5.321  [5.322] <-- first point, bottom corner of first 2d slice
      (0,4,3,3) ...
      ...             v--------------------------------------- last point, top corner
      (0,7,3,3)    [2.233]  2.234   6.230   6.231   6.232
      (0,7,4,3)     2.243   2.244   6.240   6.241   6.242
      (0,7,5,3)     3.203   3.204   7.200   7.201   7.202
      (0,7,6,3)     3.213   3.214   7.210   7.211   7.212
      (0,7,7,3)     3.223   3.224   7.220   7.221   7.222
    */
    boxstart[0] = 3; boxstart[1] = 3; boxstart[2] = 3;
    boxcount[0] = 5; boxcount[1] = 5; boxcount[2] = 5;
    box = adios_selection_boundingbox (3, boxstart, boxcount);

    // back diagonal of the 5x5x5 array
    start[0]  = 0; start[1]  = 4; start[2]  = 4;
    start[3]  = 1; start[4]  = 3; start[5]  = 3;
    start[6]  = 2; start[7]  = 2; start[8]  = 2;
    start[9]  = 3; start[10] = 1; start[11] = 1;
    start[12] = 4; start[13] = 0; start[14] = 0;

    pts = adios_selection_points(3, 5, start);
    pts->u.points.container_selection = box;
    expected[0] = 5.322;
    expected[1] = 5.411;
    expected[2] = 7.000;
    expected[3] = 2.144;
    expected[4] = 2.233;
    log ("  Read back diagonal of box with five 3D points in a 3D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    if (err) {
        adios_selection_delete(pts);
        goto endread;
    }
    expected[5] = 15.322;
    expected[6] = 15.411;
    expected[7] = 17.000;
    expected[8] = 12.144;
    expected[9] = 12.233;
    log ("  Read back diagonal of box with five 3D points in a 3D bounding box at two steps\n");
    err = test_read (f, pts, 0, 2, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;


    // Test 6-7
    // Read several single points with 1D offset in an 5x5x5 bounding box
    // back diagonal of that 5x5x5 block, = 5.332  5.411  7.000  2.144  2.233
    // Limit the query to the bounding box of the middle of the global array
    boxstart[0] = 3; boxstart[1] = 3; boxstart[2] = 3;
    boxcount[0] = 5; boxcount[1] = 5; boxcount[2] = 5;
    box = adios_selection_boundingbox (3, boxstart, boxcount);

    // back diagonal of the 5x5x5 array
    start[0] =                 4*ldim2 + 4;
    start[1] = 1*ldim2*ldim3 + 3*ldim2 + 3;
    start[2] = 2*ldim2*ldim3 + 2*ldim2 + 2;
    start[3] = 3*ldim2*ldim3 + 1*ldim2 + 1;
    start[4] = 4*ldim2*ldim3;
    pts = adios_selection_points(1, 5, start);
    pts->u.points.container_selection = box;
    expected[0] = 5.322;
    expected[1] = 5.411;
    expected[2] = 7.000;
    expected[3] = 2.144;
    expected[4] = 2.233;
    log ("  Read back diagonal of box with five 1D offsets in a 2D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    if (err) {
        adios_selection_delete(pts);
        goto endread;
    }
    expected[5] = 15.322;
    expected[6] = 15.411;
    expected[7] = 17.000;
    expected[8] = 12.144;
    expected[9] = 12.233;
    log ("  Read back diagonal of box with five 1D points in a 2D bounding box at two steps\n");
    err = test_read (f, pts, 0, 2, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;


    // Test 8
    // Read "center cross" with 3D offset in an 5x5x5 bounding box
    // points concentrated in 3x1x3 2D X-Z center plane to test the reduction of reading box size
    // 5-point center of that 3x1x3 block, =          5.400
    //                                          3.004 7.000 7.001
    //                                                7.100
    // The box-reduction should read a 3x3=9 element box instead of
    // the full 5x5x5=125 elements
    // (actually read_var_bb then will read more than this because it reads contiguous arrays)
    // Limit the query to the bounding box of the middle of the global array
    /*
     * See the cross in file:
    $ bpls -l read_points_3d.bp -d data -n 5 -s "0,4,5,5" -c "1,1,1,1" -f "%6.3f " -n 1
    real     data                       2*{10, 10, 10} = 0 / 17.444 / 8.722 / 5.50184
    slice (0:0, 4:4, 5:5, 5:5)
    (0,4,5,5)     5.400

    $ bpls -l read_points_3d.bp -d data -n 5 -s "0,5,5,4" -c "1,1,1,3" -f "%6.3f " -n 1
    real     data                       2*{10, 10, 10} = 0 / 17.444 / 8.722 / 5.50184
    slice (0:0, 5:5, 5:5, 4:6)
    (0,5,5,4)     3.004
    (0,5,5,5)     7.000
    (0,5,5,6)     7.001

    $ bpls -l read_points_3d.bp -d data -n 5 -s "0,6,5,5" -c "1,1,1,1" -f "%6.3f " -n 1
    real     data                       2*{10, 10, 10} = 0 / 17.444 / 8.722 / 5.50184
    slice (0:0, 6:6, 5:5, 5:5)
    (0,6,5,5)     7.100
     */
    boxstart[0] = 3; boxstart[1] = 3; boxstart[2] = 3;
    boxcount[0] = 5; boxcount[1] = 5; boxcount[2] = 5;
    box = adios_selection_boundingbox (3, boxstart, boxcount);

    // 5-point center cross in X-Z plane of the 5x5x5 array
    start[0]  = 1; start[1]  = 2; start[2]  = 2;
    start[3]  = 2; start[4]  = 2; start[5]  = 1;
    start[6]  = 2; start[7]  = 2; start[8]  = 2;
    start[9]  = 2; start[10] = 2; start[11] = 3;
    start[12] = 3; start[13] = 2; start[14] = 2;
    pts = adios_selection_points(3, 5, start);
    pts->u.points.container_selection = box;
    expected[0] = 5.400;
    expected[1] = 3.004;
    expected[2] = 7.000;
    expected[3] = 7.001;
    expected[4] = 7.100;
    log ("  Read back center cross in X-Z plane with five 3D points in a 3D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;


    // Test 9
    // Read "center cross" with 1D offset in an 5x5x5 bounding box
    // points concentrated in 3x1x3 2D X-Z center plane to test the reduction of reading box size
    // 5-point center of that 3x1x3 block, = 5.400 3.004 7.000 7.001 7.100
    // The box-reduction should read a 3x5x5=75 element box instead of
    // the full 5x5x5=125 elements
    // Limit the query to the bounding box of the middle of the global array
    // Same test as Test 8 but with 1D points

    boxstart[0] = 3; boxstart[1] = 3; boxstart[2] = 3;
    boxcount[0] = 5; boxcount[1] = 5; boxcount[2] = 5;
    box = adios_selection_boundingbox (3, boxstart, boxcount);

    // 5-point center cross in X-Z plane of the 5x5x5 array
    start[0] = 1*ldim2*ldim3 + 2*ldim2 + 2;
    start[1] = 2*ldim2*ldim3 + 2*ldim2 + 1;
    start[2] = 2*ldim2*ldim3 + 2*ldim2 + 2;
    start[3] = 2*ldim2*ldim3 + 2*ldim2 + 3;
    start[4] = 3*ldim2*ldim3 + 2*ldim2 + 2;

    pts = adios_selection_points(1, 5, start);
    pts->u.points.container_selection = box;
    expected[0] = 5.400;
    expected[1] = 3.004;
    expected[2] = 7.000;
    expected[3] = 7.001;
    expected[4] = 7.100;
    log ("  Read back center cross in X-Z plane with five 1D points in a 3D bounding box at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;


    /*
     * Points in WRITEBLOCK based tests
     */
    log ("  --------------- Points in WriteBlocks -----------------------  \n");

    // Test 8
    // Read a single point with 3D coordinates in writeblock 2 (third block)
    // block 2: [0:4, 5:9, 0:4] = 1 / 1.444/ 1.222/ 0.142134
    // middle point of that 5x5x5 block, = 1.222
    wblock = adios_selection_writeblock(2);

    start[0] = 2; start[1] = 2; start[2] = 2; // 3D point of middle of a 5x5x5 writeblock
    pts = adios_selection_points(3, 1, start);
    pts->u.points.container_selection = wblock;
    expected[0] = 1.222;
    log ("  Read single 3D point in a 3D writeblock at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    //adios_selection_delete (wblock); // deleted when pts is deleted
    if (err)
        goto endread;


    // Test 9
    // Read a single point with 1D offset in an writeblock 2
    // middle point of that 5x5x5 block, = 1.222
    // This is the same actual point as in Test 8
    wblock = adios_selection_writeblock(2);

    start[0] = 62; // offset of middle of a 5x5x5 array
    pts = adios_selection_points(1, 1, start);
    pts->u.points.container_selection = wblock;
    expected[0] = 1.222;
    log ("  Read single 1D point in a 3D writeblock at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;


    // Test 10-11
    // Read several single points with 3D coordinates in writeblock 2
    // back diagonal of that 5x5x5 block, = 1.044 1.133 1.222 1.311 1.400
    /* See
     * $ bpls -l read_points_3d.bp -d data -n 5 -s "0,0,5,0" -c "1,5,5,5" -f "%6.3f "
       real     data                       2*{10, 10, 10} = 0 / 17.444 / 8.722 / 5.50184
    slice (0:0, 0:4, 5:9, 0:4)
    (0,0,5,0)     1.000   1.001   1.002   1.003   1.004
    (0,0,6,0)     1.010   1.011   1.012   1.013   1.014
    (0,0,7,0)     1.020   1.021   1.022   1.023   1.024
    (0,0,8,0)     1.030   1.031   1.032   1.033   1.034
    (0,0,9,0)     1.040   1.041   1.042   1.043  [1.044] <-- first point
    ...
                    v--------------------------------------- last point
    (0,4,5,0)    [1.400]  1.401   1.402   1.403   1.404
    (0,4,6,0)     1.410   1.411   1.412   1.413   1.414
    (0,4,7,0)     1.420   1.421   1.422   1.423   1.424
    (0,4,8,0)     1.430   1.431   1.432   1.433   1.434
    (0,4,9,0)     1.440   1.441   1.442   1.443   1.444
     *
     *
     */

    wblock = adios_selection_writeblock(2);

    // back diagonal of the 5x5x5 array
    start[0]  = 0; start[1]  = 4; start[2]  = 4;
    start[3]  = 1; start[4]  = 3; start[5]  = 3;
    start[6]  = 2; start[7]  = 2; start[8]  = 2;
    start[9]  = 3; start[10] = 1; start[11] = 1;
    start[12] = 4; start[13] = 0; start[14] = 0;
    pts = adios_selection_points(3, 5, start);
    pts->u.points.container_selection = wblock;
    expected[0] = 1.044;
    expected[1] = 1.133;
    expected[2] = 1.222;
    expected[3] = 1.311;
    expected[4] = 1.400;
    log ("  Read back diagonal of box with five 3D points in a 3D writeblock at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    if (err) {
        adios_selection_delete(pts);
        goto endread;
    }

    expected[5] = 11.044;
    expected[6] = 11.133;
    expected[7] = 11.222;
    expected[8] = 11.311;
    expected[9] = 11.400;
    log ("  Read back diagonal of box with five 3D points in a 3D writeblock at two steps\n");
    err = test_read (f, pts, 0, 2, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;


    // Test 12-13
    // Read several single points with 1D offset in writeblock 2
    // back diagonal of that 5x5x5 block, = 1.044 1.133 1.222 1.311 1.400
    // This is the same actual point as in Test 10-11
    wblock = adios_selection_writeblock(2);

    // back diagonal of the 5x5x5 array as 1D offsets
    start[0] =                 4*ldim2 + 4;
    start[1] = 1*ldim2*ldim3 + 3*ldim2 + 3;
    start[2] = 2*ldim2*ldim3 + 2*ldim2 + 2;
    start[3] = 3*ldim2*ldim3 + 1*ldim2 + 1;
    start[4] = 4*ldim2*ldim3;
    pts = adios_selection_points(1, 5, start);
    pts->u.points.container_selection = wblock;
    expected[0] = 1.044;
    expected[1] = 1.133;
    expected[2] = 1.222;
    expected[3] = 1.311;
    expected[4] = 1.400;
    log ("  Read back diagonal of box with five 1D offsets in a 3D writeblock at one step\n");
    err = test_read (f, pts, 0, 1, expected);
    if (err) {
        adios_selection_delete(pts);
        goto endread;
    }
    expected[5] = 11.044;
    expected[6] = 11.133;
    expected[7] = 11.222;
    expected[8] = 11.311;
    expected[9] = 11.400;
    log ("  Read back diagonal of box with five 1D points in a 3D writeblock at two steps\n");
    err = test_read (f, pts, 0, 2, expected);
    adios_selection_delete(pts);
    if (err)
        goto endread;

endread:

    adios_read_close(f);
    MPI_Barrier (comm);
    return err;
}

