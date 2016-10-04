/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: 
 *  Write a 2D array of 2D blocks. Each block is a 5x5 array.
 *  The whole array is patterned like this:
 *
 *  1 2 3 4 ...
 *  2 3 4 ...
 *  3 4 ...
 *  4 ..
 *
 *  Then test the minmax query method if it returns the correct results
 *
 * How to run: ./query_minmax <N> <steps>
 * It writes N*N 2D blocks organized into a NxN 2D array.
 * Output: query_minmax.bp
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
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

static const char FILENAME[] = "query_minmax.bp";
#define VALUE(rank, step) (step * 1000 + rank + 1)

#define LDIM1 5
#define LDIM2 5
static const int ldim1 = LDIM1;
static const int ldim2 = LDIM2;
int gdim1, gdim2;
int offs1, offs2;

int64_t       m_adios_group;

/* Variables to write */
float  a2[LDIM1*LDIM2];

/* Variables to read */
float  r2[LDIM1*LDIM2];

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
    float v_start = 1.0*(step+1) + row*ldim1 + col*ldim2;
    float v;
    int i, j, k;

    n = ldim1 * ldim2;
    //log ("  Fill up array of %d elements starting from value %f...\n",n,v_start);
    k = 0;
    for (i=0; i<ldim1; i++) {
    	v = v_start + i*1.0;
        //log ("      row %d starts from value %f... (element %d)\n",i,v,k);
    	for (j=0; j<ldim2; j++) {
    		a2[k] = v;
    		k++;
    		v++;
    	}
    }
}


void Usage() 
{
    printf("Usage: query_minmax <N> <nsteps>\n"
            "    <N>:       Number of blocks in each of X and Y direction\n"
    		"    <nsteps>:  Number of write cycles (to same file)\n");
}

void define_vars ();
int write_file (int step);
int query_as_file ();

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
        printf("Running query_minmax <N=%d> <nsteps=%d>\n", N, NSTEPS);
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

    adios_declare_group (&m_adios_group, "query_minmax", "", adios_stat_default);
    adios_select_method (m_adios_group, "POSIX", "", "");


    define_vars();
    set_gdim();
    
    for (i=0; i<NSTEPS; i++) {
        if (!err) {
            err = write_file (i); 
        }
    }

    if (!err)
        err = query_as_file ();

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
    uint64_t       groupsize=0, totalsize;
    int           nblocks = N*N;
    int           i, j;
    double        tb, te;

    log ("Write step %d to %s\n", step, FILENAME);
    adios_open (&fh, "query_minmax", FILENAME, (step ? "a" : "w"), comm);
    
    groupsize  = (4 + nblocks*2) * sizeof(int);             // dimensions
    log ("  groupsize calculated = %llu\n", groupsize);
    groupsize += nblocks * ldim1 * ldim2 * sizeof(float);     // 2D  blocks
    log ("  groupsize calculated = %llu\n", groupsize);

    adios_group_size (fh, groupsize, &totalsize);
    log ("  groupsize %llu, totalsize %llu\n", groupsize, totalsize);

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

static const char minstr[] = "11.0";
static const char maxstr[] = "21.0";

int query_test (ADIOS_FILE *f, ADIOS_SELECTION *boxsel)
{
    int nerr = 0;
    int step;
    ADIOS_QUERY  *q1, *q2, *q;
    double        tb, te;
    double        teb, t_eval; // time for just evaluating query for one step

    q1 = adios_query_create (f, boxsel, "data", ADIOS_GTEQ, minstr);
    q2 = adios_query_create (f, boxsel, "data", ADIOS_LTEQ, maxstr);
    q  = adios_query_combine (q1, ADIOS_QUERY_OP_AND, q2);

    if (q == NULL)  {
        log ("ERROR: Query creation failed: %s\n", adios_errmsg());
        return 1;
    }

    // We can call this with unknown too, just testing the default behavior here
    adios_query_set_method (q, ADIOS_QUERY_METHOD_MINMAX);

    for (step=0; step<NSTEPS; step++) {
        tb = MPI_Wtime();
        t_eval = 0;

        // retrieve the whole query result at once
        int64_t batchSize = N*N;
        log ("    set upper limit to total number of blocks in array  = %lld\n", batchSize);

        teb = MPI_Wtime();
        ADIOS_QUERY_RESULT *result = adios_query_evaluate(q, boxsel, step, batchSize);
        t_eval += MPI_Wtime() - teb;

        if (result->status == ADIOS_QUERY_RESULT_ERROR) {
            log ("ERROR: Query evaluation failed with error: %s\n", adios_errmsg());
            nerr++;
        }
        else if (result->status == ADIOS_QUERY_HAS_MORE_RESULTS)
        {
            log ("ERROR: Query retrieval failure: "
                   "it says it has more results to retrieve, although we "
                   "tried to get all at once\n");
            nerr++;
        }

        log ("    Query returned %d blocks as result\n", result->nselections);
        free (result->selections);

        MPI_Barrier (comm);
        te = MPI_Wtime();
        log ("  Processing time for step %d was %6.3lfs, query evaluation was %6.3lfs\n",
                step, te-tb, t_eval);
    }

    adios_query_free(q);
    adios_query_free(q2);
    adios_query_free(q1);
    return nerr;
}


int query_as_file ()
{
    ADIOS_FILE * f;
    ADIOS_VARINFO * vi;
    int err=0;
    double        tb, te;

    uint64_t start[2] = {0,0};
    uint64_t count[2] = {gdim1,gdim2};
    
    reset_readvars();

    log ("Query data in %s\n", FILENAME);
    tb = MPI_Wtime();
    f = adios_read_open_file (FILENAME, ADIOS_READ_METHOD_BP, comm);
    if (f == NULL) {
        printE ("Error at opening file: %s\n", adios_errmsg());
        return 1;
    }
    te = MPI_Wtime();
    log ("  File opening time was %6.3lfs\n", te-tb);

    log ("  Check variable definitions... %s\n", FILENAME);
    CHECK_VARINFO("data", 2, NSTEPS)
    MPI_Barrier (comm);


    // limit the query to the bounding box of this process
    ADIOS_SELECTION *boxsel = adios_selection_boundingbox (2, start, count);
    log ("  Query variable content with a bounding box selection...\n");
    err += query_test (f, boxsel);
    adios_selection_delete (boxsel);

    // test with NULL bounding box selection, too
    if (!err) {
        log ("  Query variable content with NULL as bounding box selection...\n");
        err += query_test (f, NULL);
    }

endread:


    adios_read_close(f);
    MPI_Barrier (comm);
    return err;
}

