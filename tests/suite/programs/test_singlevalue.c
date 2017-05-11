/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: 
 *  Write the first variable only once: one process and one step only. 
 *  Then write some other variables over multiple steps.
 *  Then read them all and check if they are correct. 
 *
 *  This situation broke reading step-by-step because bp_utils.c:get_time() 
 *  uses the the first variable to find the next timestep for all.
 *
 * How to run: mpirun -np <N> test_singlevalue
 * Output: test_singlevalue.bp
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "adios.h"
#include "adios_read.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define log(...) fprintf (stderr, "[rank=%3.3d, line %d]: ", rank, __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);
#define printE(...) fprintf (stderr, "[rank=%3.3d, line %d]: ERROR: ", rank, __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);

static const int NSTEPS = 3;
static const char FILENAME[] = "test_singlevalue.bp";
#define VALUE(rank, step) (rank * 100 + step)
#define VALUE0(step) (step)

/* Variables to write */
int a0;
int  *a1;


/* Variables to read */
int r0, rt0;
int  *r1, *rt1;

static const int ldim1 = 7;
int gdim1;
int offs1;

MPI_Comm    comm = MPI_COMM_WORLD;
int rank;
int size;

void alloc_vars()
{
    int n;

    gdim1 = size*ldim1;
    offs1 = rank*ldim1;
    n = ldim1;
    a1  = (int*) malloc (n * sizeof(int));
    r1  = (int*) malloc (n * sizeof(int));
    rt1 = (int*) malloc (n * sizeof(int));
}

void set_vars(int step)
{
    int n, i;
    int v = VALUE(rank, step);

    a0 = v;

    n = ldim1;
    for (i=0; i<n; i++) a1[i] = v;
}

void fini_vars()
{
    free (a1);
    free (r1);
    free (rt1);
}


int write_file (int step);
int read_file ();
int read_stream ();

int main (int argc, char ** argv) 
{
    int err,i ; 

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    alloc_vars();
    adios_init_noxml(comm);
    err = adios_read_init_method(ADIOS_READ_METHOD_BP, comm, "verbose=2");
    if (err) {
        printE ("%s\n", adios_errmsg());
    }

    int64_t       m_adios_group;
    adios_declare_group (&m_adios_group, "once", "", adios_stat_default);
    adios_select_method (m_adios_group, "MPI", "", "");
    adios_define_var (m_adios_group, "once", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "onceperstep", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "scalar", "", adios_integer, 0, 0, 0);
    char gdimstr[32], ldimstr[32], offsstr[32];
    sprintf(gdimstr, "%d", gdim1);
    sprintf(ldimstr, "%d", ldim1);
    sprintf(offsstr, "%d", offs1);
    adios_define_var (m_adios_group, "a1", "", adios_integer,
            ldimstr, gdimstr, offsstr);
    
    for (i=0; i<NSTEPS; i++) {
        if (!err) {
            set_vars (i);
            err = write_file (i); 
        }
    }

    if (!err)
        err = read_file (); 

    if (!err)
        err = read_stream (); 

    adios_finalize (rank);
    fini_vars();
    MPI_Finalize ();
    return err;
}


int write_file (int step) 
{
    int64_t       fh;
    uint64_t       groupsize=0, totalsize;

    log ("Write step %d to %s\n", step, FILENAME);
    adios_open (&fh, "once", FILENAME, (step ? "a" : "w"), comm);
    
    if (rank == 0) {
        if (step == 0) {
            adios_write (fh, "once", &a0);
        }
        adios_write (fh, "onceperstep", &a0);
    }
    adios_write (fh, "scalar", &a0);
    adios_write (fh, "a1", a1);

    adios_close (fh);
    MPI_Barrier (comm);
    return 0;
}


#define CHECK_VARINFO(VARNAME, NDIM, NSTEPS) \
    vi = adios_inq_var (f, VARNAME); \
    if (vi == NULL) { \
        printE ("No such variable " VARNAME "\n"); \
        err = 101; \
        goto endread; \
    } \
    if (vi->ndim != NDIM) { \
        printE ("Variable " VARNAME " has %d dimensions, but expected %d\n", vi->ndim, NDIM); \
        err = 102; \
        goto endread; \
    } \
    if (vi->nsteps != NSTEPS) { \
        printE ("Variable " VARNAME " has %d steps, but expected %d\n", vi->nsteps, NSTEPS); \
        err = 103; \
        /*goto endread; */\
    } \
    adios_free_varinfo (vi);

#define CHECK_SCALAR(VARNAME, VAR, VALUE, STEP) \
    if (VAR != VALUE) { \
        printE (#VARNAME " step %d: wrote %d but read %d\n",STEP,VALUE,VAR);\
        err = 104; \
        /*goto endread;*/\
    }

#define CHECK_ARRAY(VARNAME,A,N,VALUE,STEP,i) \
    for (i=0;i<N;i++) \
        if (A[i] != VALUE) { \
            printE (#VARNAME "[%d] step %d: wrote %d but read %d\n",i,STEP,VALUE,A[i]);\
            err = 104; \
            /*goto endread;*/\
            break; \
        }

void reset_readvars()
{
    int n;
    
    r0  = -1;
    rt0 = -1;

    n = ldim1;
    memset (r1,  -1, n*sizeof(int));
    memset (rt1, -1, n*sizeof(int));

}

int read_file ()
{
    ADIOS_SELECTION *sel0,*sel1;
    ADIOS_FILE * f;
    ADIOS_VARINFO * vi;
    int err=0,v,v0,i,n;
    int iMacro; // loop variable in macros

    uint64_t start[1] = {offs1};
    uint64_t count[1] = {ldim1};
    uint64_t ndim;
    
    reset_readvars();

    log ("Read and check data in %s\n", FILENAME);
    f = adios_read_open_file (FILENAME, ADIOS_READ_METHOD_BP, comm);
    if (f == NULL) {
        printE ("Error at opening file: %s\n", adios_errmsg());
        return 1;
    }

    sel0 = adios_selection_boundingbox (0, start, count); 
    sel1 = adios_selection_boundingbox (1, start, count); 


    log ("  Check variable definitions... %s\n", FILENAME);

    CHECK_VARINFO("once", 0, 1)
    CHECK_VARINFO("onceperstep", 0, NSTEPS)
    CHECK_VARINFO("scalar", 0, NSTEPS)
    CHECK_VARINFO("a1", 1, NSTEPS)

    log ("  Check variable 'once'...\n");
    v0 = VALUE0(0);
    adios_schedule_read (f, sel0, "once",  0, 1, &r0);
    adios_perform_reads (f, 1);
    CHECK_SCALAR (once,  r0,  v0, i) // scalar is from writer rank 0, not this rank!

    log ("  Check variables 'onceperstep', 'scalar', 'a1', steps=%d...\n",NSTEPS);
    for (i=0; i<NSTEPS; i++) {
        v = VALUE(rank,i);
        v0 = VALUE0(i);
        log ("    Step %d value %d\n", i, v);

        adios_schedule_read (f, sel0, "onceperstep",  i, 1, &r0);
        adios_schedule_read (f, sel0, "scalar", i, 1, &rt0);
        adios_schedule_read (f, sel1, "a1",  i, 1, r1);
        adios_perform_reads (f, 1);

        CHECK_SCALAR (onceperstep,  r0,  v0, i) // scalar is from writer rank 0, not this rank!
        CHECK_SCALAR (scalar, rt0, v0, i) // so value is v0 at step 'i', not v
        CHECK_ARRAY (a1,  r1,  ldim1, v, i, iMacro)
    } 

endread:

    adios_selection_delete (sel0);
    adios_selection_delete (sel1);

    adios_read_close(f);
    MPI_Barrier (comm);
    return err;
}





int read_stream ()
{
    ADIOS_SELECTION *sel0,*sel1;
    ADIOS_FILE * f;
    ADIOS_VARINFO * vi;
    int err=0,v,v0,i,n;
    int iMacro; // loop variable in macros

    uint64_t start[3] = {offs1};
    uint64_t count[3] = {ldim1};
    uint64_t ndim;
    
    reset_readvars();

    log ("Read as stream and check data in %s\n", FILENAME);
    f = adios_read_open (FILENAME, ADIOS_READ_METHOD_BP, comm,
                         ADIOS_LOCKMODE_NONE, 0.0);
    if (f == NULL) {
        printE ("Error at opening file as stream: %s\n", adios_errmsg());
        return 1;
    }

    sel0 = adios_selection_boundingbox (0, start, count); 
    sel1 = adios_selection_boundingbox (1, start, count); 

    n = 0;
    while (n < NSTEPS && adios_errno != err_end_of_stream) {
        log ("  Step %d\n", n);

        log ("    Check variable definitions... %s\n", FILENAME);
        if (n == 0) {
            CHECK_VARINFO("once", 0, 1)
        }
        CHECK_VARINFO("onceperstep", 0, 1)
        CHECK_VARINFO("scalar", 0, 1)
        CHECK_VARINFO("a1", 1, 1)

        v = VALUE(rank,n);
        v0 = VALUE0(n);

        if (n == 0) {
            log ("  Check variable 'once'...\n");
            adios_schedule_read (f, sel0, "once",  n, 1, &r0);
            adios_perform_reads (f, 1);
            CHECK_SCALAR (once,  r0,  v0, n) // scalar is from writer rank 0, not this rank!
        }

        log ("  Check variables 'onceperstep', 'scalar', 'a1', Step %d value %d\n", n, v);

        adios_schedule_read (f, sel0, "onceperstep",  0, 1, &r0);
        adios_schedule_read (f, sel0, "scalar", 0, 1, &rt0);
        adios_schedule_read (f, sel1, "a1",  0, 1, r1);
        adios_perform_reads (f, 1);

        CHECK_SCALAR (onceperstep,  r0,  v0, n) // scalar is from writer rank 0, not this rank!
        CHECK_SCALAR (scalar, rt0, v0, n) // so value is v0 at 'n', not v
        CHECK_ARRAY (a1,  r1,  ldim1, v, n, iMacro)

        if (adios_advance_step (f, 0, 0.0) >= 0)
            n = f->current_step;
        else
            n++; //just to end the loop
    } 

endread:

    adios_selection_delete (sel0);
    adios_selection_delete (sel1);

    adios_read_close(f);
    MPI_Barrier (comm);
    return err;
}
