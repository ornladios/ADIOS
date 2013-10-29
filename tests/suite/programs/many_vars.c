/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: 
 *  Write a huge number of variables
 *  Then read them all and check if they are correct. 
 *
 * How to run: mpirun -np <N> many_vars <nvars> <blocks per process> <steps>
 * Output: many_vars.bp
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include "mpi.h"
#include "adios.h"
#include "adios_read.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define log(...) fprintf (stderr, "[rank=%3.3d, line %d]: ", rank, __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);
#define printE(...) fprintf (stderr, "[rank=%3.3d, line %d]: ERROR: ", rank, __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);

int NVARS = 1;
int NBLOCKS = 1;
int NSTEPS = 1;
static const char FILENAME[] = "many_vars.bp";
#define VALUE(rank, step, block) (step * 10000 + 10*rank + block)

/* Variables to write */
int  *a2;

static const int ldim1 = 5;
static const int ldim2 = 5;
int gdim1, gdim2;
int offs1, offs2;

int64_t       m_adios_group;

/* Variables to read */
int  *r2;

MPI_Comm    comm = MPI_COMM_WORLD;
int rank;
int size;
char ** varnames;

void alloc_vars()
{
    int n,i;

    n = ldim1 * ldim2;
    a2  = (int*) malloc (n * sizeof(int));
    r2  = (int*) malloc (n * sizeof(int));
    varnames = (char**) malloc (NVARS * sizeof(char*));
    for (i=0; i<NVARS; i++) {
        varnames[i] = (char*) malloc (16);
    }

    /* make varnames like v001,v002,.. */
    int digit=1, d=10;
    while (NVARS/d > 0) {
        d *= 10;
        digit++;
    }

    char fmt[16];
    sprintf (fmt, "v%%%d.%dd",digit,digit);
    printf ("ftm=[%s]\n", fmt);
    for (i=0; i<NVARS; i++) {
        //sprintf(varnames[i], "v%-*d", digit, i);
        sprintf(varnames[i], fmt, i);
    }
    printf ("varname[0]=%s\n", varnames[0]);
    printf ("varname[%d]=%s\n", NVARS-1, varnames[NVARS-1]);
}

void set_gdim()
{
    gdim1 = size*ldim1;
    gdim2 = NBLOCKS*ldim2;
}

void set_vars(int step, int block)
{
    int n, i;
    int v = VALUE(rank, step, block);

    offs1 = rank*ldim1;
    offs2 = block*ldim2;

    n = ldim1 * ldim2;
    for (i=0; i<n; i++) a2[i] = v;
}

void fini_vars()
{
    int i;
    free (a2);
    free (r2);
    for (i=0; i<NVARS; i++) {
        free(varnames[i]);
    }
    free(varnames);
}

void Usage() 
{
    printf("Usage: many_vars <nvars> <nblocks> <nsteps>\n" 
            "    <nvars>:   Number of variables to generate\n"
            "    <nblocks>: Number of blocks per process to write\n"
            "    <nsteps>:  Number of write cycles (to same file)\n");
}

void define_vars ();
int write_file (int step);
int read_file ();

int main (int argc, char ** argv) 
{
    int err,i ; 

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    if (argc < 2) { Usage(); return 1; }

    errno = 0;
    i = strtol (argv[1], NULL, 10);
    if (errno || i < 1) { printf("Invalid 1st argument %s\n", argv[1]); Usage(); return 1;}
    NVARS = i;

    if (argc > 2) {
        i = strtol (argv[2], NULL, 10);
        if (errno || i < 1) { printf("Invalid 2nd argument %s\n", argv[2]); Usage(); return 1;}
        NBLOCKS = i;
    }

    if (argc > 3) {
        i = strtol (argv[3], NULL, 10);
        if (errno || i < 1) { printf("Invalid 3rd argument %s\n", argv[3]); Usage(); return 1;}
        NSTEPS = i;
    }

    alloc_vars();
    adios_init_noxml (comm);
    adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10);
    err = adios_read_init_method(ADIOS_READ_METHOD_BP, comm, "verbose=2");
    if (err) {
        printE ("%s\n", adios_errmsg());
    }

    adios_declare_group (&m_adios_group, "multiblock", "iter", adios_flag_yes);
    adios_select_method (m_adios_group, "MPI", "", "");


    define_vars();
    set_gdim();
    
    for (i=0; i<NSTEPS; i++) {
        if (!err) {
            err = write_file (i); 
        }
    }

    //if (!err)
    //    err = read_file (); 

    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    adios_finalize (rank);
    fini_vars();
    MPI_Finalize ();
    return err;
}

void define_vars ()
{
    int i, block;

    adios_define_var (m_adios_group, "ldim1", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "ldim2", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "gdim1", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "gdim2", "", adios_integer, 0, 0, 0);

    for (block=0; block<NBLOCKS; block++) {
        adios_define_var (m_adios_group, "offs1", "", adios_integer, 0, 0, 0);
        adios_define_var (m_adios_group, "offs2", "", adios_integer, 0, 0, 0);

        for (i=0; i<NVARS; i++) {
            adios_define_var (m_adios_group, varnames[i], "", adios_integer, 
                    "iter,ldim1,ldim2",
                    "gdim1,gdim2",
                    "offs1,offs2");
        }
    }
}

int write_file (int step) 
{
    int64_t       fh;
    uint64_t       groupsize=0, totalsize;
    int           block, v, i;

    log ("Write step %d to %s\n", step, FILENAME);
    adios_open (&fh, "multiblock", FILENAME, (step ? "a" : "w"), comm);
    
    groupsize  = (4 + NBLOCKS*2) * sizeof(int);             // dimensions 
    groupsize += NVARS * NBLOCKS * ldim1 * ldim2 * sizeof(int);     // 2D  blocks
    //groupsize +=1024;

    adios_group_size (fh, groupsize, &totalsize);
    log ("  groupsize %lld, totalsize %lld\n", groupsize, totalsize);

    adios_write (fh, "gdim1", &gdim1);
    adios_write (fh, "gdim2", &gdim2);
    adios_write (fh, "ldim1", &ldim1);
    adios_write (fh, "ldim2", &ldim2);

    for (block=0; block<NBLOCKS; block++) {
        v = VALUE(rank, step, block);
        log ("  Write block %d, value %d to %s\n", block, v, FILENAME);
        set_vars (step, block);

        adios_write (fh, "offs1", &offs1);
        adios_write (fh, "offs2", &offs2);

        for (i=0; i<NVARS; i++) {
            adios_write (fh, varnames[i], a2);
        }
    }

    adios_close (fh);
    MPI_Barrier (comm);
    return 0;
}

#if 0
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
    
    n = ldim1 * ldim2;
    memset (r2,  -1, n*sizeof(int));
}

int read_file ()
{
    ADIOS_SELECTION *sel2;
    ADIOS_FILE * f;
    ADIOS_VARINFO * vi;
    int err=0,v,v0,i,n;
    int nsteps_a, nsteps_b, nsteps_c;
    int iMacro; // loop variable in macros

    uint64_t start[3] = {offs1,offs2,offs3};
    uint64_t count[3] = {ldim1,ldim2,ldim3};
    uint64_t ndim;
    
    reset_readvars();

    log ("Read and check data in %s\n", FILENAME);
    f = adios_read_open_file (FILENAME, ADIOS_READ_METHOD_BP, comm);
    if (f == NULL) {
        printE ("Error at opening file: %s\n", rank, adios_errmsg());
        return 1;
    }

    sel0 = adios_selection_boundingbox (0, start, count); 
    sel1 = adios_selection_boundingbox (1, start, count); 
    sel2 = adios_selection_boundingbox (2, start, count); 
    sel3 = adios_selection_boundingbox (3, start, count); 

    log ("  Check variable definitions... %s\n", FILENAME);
    nsteps_a = NSTEPS / 2 + NSTEPS % 2;
    nsteps_b = NSTEPS / 2;
    nsteps_c = NSTEPS;

    CHECK_VARINFO("a0", 0, nsteps_a)
    CHECK_VARINFO("b0", 0, nsteps_b)
    CHECK_VARINFO("c0", 0, nsteps_c)
    CHECK_VARINFO("at0", 0, nsteps_a)
    CHECK_VARINFO("bt0", 0, nsteps_b)
    CHECK_VARINFO("ct0", 0, nsteps_c)
    CHECK_VARINFO("a1", 1, nsteps_a)
    CHECK_VARINFO("b1", 1, nsteps_b)
    CHECK_VARINFO("c1", 1, nsteps_c)
    CHECK_VARINFO("at1", 1, nsteps_a)
    CHECK_VARINFO("bt1", 1, nsteps_b)
    CHECK_VARINFO("ct1", 1, nsteps_c)
    CHECK_VARINFO("a2", 2, nsteps_a)
    CHECK_VARINFO("b2", 2, nsteps_b)
    CHECK_VARINFO("c2", 2, nsteps_c)
    CHECK_VARINFO("a3", 3, nsteps_a)
    CHECK_VARINFO("b3", 3, nsteps_b)
    CHECK_VARINFO("c3", 3, nsteps_c)

    log ("  Check variables c0,ct0,c1,ct1,c2,c3, steps=%d...\n",nsteps_c);
    for (i=0; i<nsteps_c; i++) {
        v = VALUE(rank,i);
        v0 = VALUE0(i);
        log ("    Step %d value %d\n", i, v);

        adios_schedule_read (f, sel0, "c0",  i, 1, &r0);
        adios_schedule_read (f, sel0, "ct0", i, 1, &rt0);
        adios_schedule_read (f, sel1, "c1",  i, 1, r1);
        adios_schedule_read (f, sel1, "ct1", i, 1, rt1);
        adios_schedule_read (f, sel2, "c2",  i, 1, r2);
        adios_schedule_read (f, sel3, "c3",  i, 1, r3);
        adios_perform_reads (f, 1);

        CHECK_SCALAR (c0,  r0,  v0, i) // scalar is from writer rank 0, not this rank!
        CHECK_SCALAR (ct0, rt0, v0, i) // so value is v0 at step 'i', not v
        CHECK_ARRAY (c1,  r1,  ldim1, v, i, iMacro)
        CHECK_ARRAY (ct1, rt1, ldim1, v, i, iMacro)
        CHECK_ARRAY (c2,  r2,  ldim1*ldim2, v, i, iMacro)
        CHECK_ARRAY (c3,  r3,  ldim1*ldim2*ldim3, v, i, iMacro)
    } 

    log ("  Check variables a0,at0,a1,at1,a2,a3, steps=%d...\n",nsteps_a);
    for (i=0; i<nsteps_a; i++) {
        v = VALUE(rank,i*2);
        v0 = VALUE0(i*2);
        log ("    Step %d value %d\n", i, v);

        adios_schedule_read (f, sel0, "a0",  i, 1, &r0);
        adios_schedule_read (f, sel0, "at0", i, 1, &rt0);
        adios_schedule_read (f, sel1, "a1",  i, 1, r1);
        adios_schedule_read (f, sel1, "at1", i, 1, rt1);
        adios_schedule_read (f, sel2, "a2",  i, 1, r2);
        adios_schedule_read (f, sel3, "a3",  i, 1, r3);
        adios_perform_reads (f, 1);

        CHECK_SCALAR (a0,  r0,  v0, i)
        CHECK_SCALAR (at0, rt0, v0, i)
        CHECK_ARRAY (a1,  r1,  ldim1, v, i, iMacro)
        CHECK_ARRAY (at1, rt1, ldim1, v, i, iMacro)
        CHECK_ARRAY (a2,  r2,  ldim1*ldim2, v, i, iMacro)
        CHECK_ARRAY (a3,  r3,  ldim1*ldim2*ldim3, v, i, iMacro)
    } 

    log ("  Check variables b0,bt0,b1,bt1,b2,b3, steps=%d...\n",nsteps_b);
    for (i=0; i<nsteps_b; i++) {
        v = VALUE(rank,i*2)+1;
        v0 = VALUE0(i*2)+1;
        log ("    Step %d value %d\n", i, v);

        adios_schedule_read (f, sel0, "b0",  i, 1, &r0);
        adios_schedule_read (f, sel0, "bt0", i, 1, &rt0);
        adios_schedule_read (f, sel1, "b1",  i, 1, r1);
        adios_schedule_read (f, sel1, "bt1", i, 1, rt1);
        adios_schedule_read (f, sel2, "b2",  i, 1, r2);
        adios_schedule_read (f, sel3, "b3",  i, 1, r3);
        adios_perform_reads (f, 1);

        CHECK_SCALAR (b0,  r0,  v0, i)
        CHECK_SCALAR (bt0, rt0, v0, i)
        CHECK_ARRAY (b1,  r1,  ldim1, v, i, iMacro)
        CHECK_ARRAY (bt1, rt1, ldim1, v, i, iMacro)
        CHECK_ARRAY (b2,  r2,  ldim1*ldim2, v, i, iMacro)
        CHECK_ARRAY (b3,  r3,  ldim1*ldim2*ldim3, v, i, iMacro)
    } 

endread:

    adios_selection_delete (sel0);
    adios_selection_delete (sel1);
    adios_selection_delete (sel2);
    adios_selection_delete (sel3);

    adios_read_close(f);
    MPI_Barrier (comm);
    return err;
}
#endif
