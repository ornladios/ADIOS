/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: 
 *  Write some local arrays with global dimension and then
 *  read them using global array reading
 *
 * How to run: mpirun -np <N> selections
 * Output: selections.bp
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include "adios.h"
#include "adios_read.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

//#define DEBUG_PRINTS

#define log(...) fprintf (stderr, "[rank=%3.3d, line %d]: ", rank, __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);
#define printE(...) fprintf (stderr, "[rank=%3.3d, line %d]: ERROR: ", rank, __LINE__); fprintf (stderr, __VA_ARGS__); fflush(stderr);

static const int NSTEPS = 1;
static const char FILENAME[] = "joinedarray.bp";
#define VALUE1D(rank,step,i) (100*(rank)+10*(step)+i)
#define VALUE2D(rank,step,i,j) (1000*(rank)+100*(step)+10*(i)+j)
#define VALUE3D(rank,step,i,j,k) (10000*(rank)+1000*(step)+100*(i)+10*(j)+k)

/* Variables to write */
int a0;
int  *a1;
int  *a2;
int  *a3;

/* Variables to read */
int r0;
int  *r1, *r2, *r3;

static const int ldim1 = 3;
static const int ldim2 = 5;
static const int ldim3 = 7;
int gdim1, gdim2, gdim3;
int offs1, offs2, offs3;

MPI_Comm    comm = MPI_COMM_WORLD;
int         rank;
int         size;
int64_t     m_adios_group;
enum ADIOS_READ_METHOD read_method = ADIOS_READ_METHOD_BP;

void alloc_vars()
{
    int n;

    gdim1 = size*ldim1;
    gdim2 = ldim2;
    gdim3 = ldim3;

    offs1 = rank*ldim1;
    offs2 = 0;
    offs3 = 0;

    n = ldim1;
    a1  = (int*) malloc (n * sizeof(int));
    r1  = (int*) malloc (n * sizeof(int));

    n = ldim1 * ldim2;
    a2  = (int*) malloc (n * sizeof(int));
    r2  = (int*) malloc (n * sizeof(int));

    n = ldim1 * ldim2 * ldim3;
    a3  = (int*) malloc (n * sizeof(int));
    r3  = (int*) malloc (n * sizeof(int));
}

void set_vars(int step)
{
    int i, j, k;

    a0 = VALUE1D(rank,step,0);

    for (i=0; i<ldim1; i++) {
        a1[i] = VALUE1D(rank,step,i);
        for (j=0; j<ldim2; j++) {
            a2[i*ldim2+j] = VALUE2D(rank,step,i,j);
            for (k=0; k<ldim3; k++) {
                a3[i*ldim2*ldim3+j*ldim3+k] = VALUE3D(rank,step,i,j,k);
            }
        }
    }
}

void fini_vars()
{
    free (a1);
    free (r1);
    free (a2);
    free (r2);
    free (a3);
    free (r3);
}


void define_vars();
int write_file (int step);
int read_file ();

void Usage() 
{
    printf("Usage: joinedarray <method> [w|r]\n"
            "    <method>:  Method, 0: file, 1: staging with dataspaces\n"
            "    w: do write only\n"
            "    r: do read only\n");
}   

int main (int argc, char ** argv) 
{
    int err, step ; 
    int do_write = 1;
    int do_read = 1;
    int m = 0;
    char write_method[16];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    if (argc > 1) {
        errno=0;
        m = strtol (argv[1], NULL, 10);
        if (errno) { printf("Invalid 1st argument %s\n", argv[1]); Usage(); return 1;}
    }
    if (argc > 2) { 
        if (argv[2][0] == 'w' || argv[2][0] == 'W') {
            do_read = 0;
        } else if (argv[2][0] == 'r' || argv[2][0] == 'R') {
            do_write = 0;
        } else {
            printE ("Invalid command line argument %s. Allowed ones:\n"
                    " w: do write only\n"
                    " r: do read only\n", argv[2]);
            MPI_Finalize ();
            return 1;
        }
    }

    if (m==0) {
        read_method = ADIOS_READ_METHOD_BP;
        strcpy(write_method,"MPI");
    } else {
        read_method = ADIOS_READ_METHOD_DATASPACES;
        strcpy(write_method,"DATASPACES");
    }

    
    log ("Writing: %s method=%s\n"
         "Reading: %s method=%d\n", 
        (do_write ? "yes" : "no"), write_method,
        (do_read ? "yes" : "no"), read_method);

    alloc_vars();
    if (do_write) {
        adios_init_noxml (comm);
    }
    if (do_read) {
        err = adios_read_init_method(read_method, comm, "verbose=2");
        if (err) {
            printE ("%s\n", adios_errmsg());
        }
    }

    
    if (do_write) {
        adios_declare_group (&m_adios_group, "joineddim", "", adios_stat_default);
        adios_select_method (m_adios_group, write_method, "verbose=2", "");

        define_vars();

        for (step=0; step<NSTEPS; step++) {
            if (!err) {
                set_vars (step);
                err = write_file (step); 
                sleep(1);
            }
        }
        adios_free_group (m_adios_group);
    }

    if (!err && do_read)
        err = read_file ();


    if (do_read) {
        adios_read_finalize_method (read_method);
    }
    fini_vars();
    if (do_write) {
        adios_finalize (rank);
    }
    MPI_Finalize ();
    return err;
}


void define_vars ()
{
    int i;

    adios_define_var (m_adios_group, "ldim1", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "ldim2", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "ldim3", "", adios_integer, 0, 0, 0);

    adios_define_var (m_adios_group, "a1", "", adios_integer,
            "ldim1", "JoinedDim", "");

    adios_define_var (m_adios_group, "a2_1", "", adios_integer,
            "ldim1,ldim2", "JoinedDim,ldim2", "");

    adios_define_var (m_adios_group, "a2_2", "", adios_integer,
            "ldim1,ldim2", "ldim1,JoinedDim", "");

    adios_define_var (m_adios_group, "a3_1", "", adios_integer,
            "ldim1,ldim2,ldim3", "JoinedDim,ldim2,ldim3", "");

    adios_define_var (m_adios_group, "a3_2", "", adios_integer,
            "ldim1,ldim2,ldim3", "ldim1,JoinedDim,ldim3", "");

    adios_define_var (m_adios_group, "a3_3", "", adios_integer,
            "ldim1,ldim2,ldim3", "ldim1,ldim2,JoinedDim", "");
}

int write_file (int step) 
{
    int64_t       fh;
    uint64_t       groupsize=0, totalsize;

    log ("Write step %d to %s\n", step, FILENAME);
    adios_open (&fh, "joineddim", FILENAME, (step ? "a" : "w"), comm);

    adios_write (fh, "ldim1", &ldim1);
    adios_write (fh, "ldim2", &ldim2);
    adios_write (fh, "ldim3", &ldim3);

    adios_write (fh, "a1", a1);
    adios_write (fh, "a2_1", a2);
    adios_write (fh, "a2_2", a2);
    adios_write (fh, "a3_1", a3);
    adios_write (fh, "a3_2", a3);
    adios_write (fh, "a3_3", a3);

    adios_close (fh);
    MPI_Barrier (comm);
    return 0;
}

void reset_readvars()
{
    int n;
    
    r0  = -1;

    n = ldim1;
    memset (r1,  -1, n*sizeof(int));

    n = ldim1 * ldim2;
    memset (r2,  -1, n*sizeof(int));

    n = ldim1 * ldim2 * ldim3;
    memset (r3,  -1, n*sizeof(int));
}

int read_file()
{
    ADIOS_SELECTION *sel1,*sel2_1,*sel2_2,*sel3_1,*sel3_2,*sel3_3;
    ADIOS_FILE * f;
    ADIOS_VARINFO * vi;
    int err=0,n,n1, i,j,k;
    int nsteps_a, nsteps_b, nsteps_c;
    int v; 

    uint64_t start[3] = {0,0,0};
    uint64_t count[3] = {ldim1,ldim2,ldim3};
    
    reset_readvars();

    log ("Read and check data in %s using point selections\n", FILENAME);
    f = adios_read_open (FILENAME, read_method, comm,
                         ADIOS_LOCKMODE_CURRENT, 0.0);
    if (f == NULL) {
        printE ("Error at opening file: %s\n", adios_errmsg());
        return 1;
    }

    start[0] = rank*ldim1;
    sel1 = adios_selection_boundingbox(1, start, count);

    start[0] = rank*ldim1;
    start[1] = 0;
    sel2_1 = adios_selection_boundingbox(2, start, count);

    start[0] = 0;
    start[1] = rank*ldim2;
    sel2_2 = adios_selection_boundingbox(2, start, count);

    start[0] = rank*ldim1;
    start[1] = 0;
    start[2] = 0;
    sel3_1 = adios_selection_boundingbox(3, start, count);

    start[0] = 0;
    start[1] = rank*ldim2;
    start[2] = 0;
    sel3_2 = adios_selection_boundingbox(3, start, count);

    start[0] = 0;
    start[1] = 0;
    start[2] = rank*ldim3;
    sel3_3 = adios_selection_boundingbox(3, start, count);

    n1=0;
    while (n1 < NSTEPS && adios_errno != err_end_of_stream) {
        n1++;
        log ("  Step %d\n", f->current_step);


        log ("  Check 1D variable a1...\n");
        adios_schedule_read (f, sel1, "a1",  0, 1, r1);
        adios_perform_reads (f, 1);

#ifdef DEBUG_PRINTS
        fprintf(stderr, "1D result: {");
        for (i=0; i<ldim1; i++) {
            fprintf (stderr, "%d ", r1[i]);
        }
        fprintf(stderr, "}\n");
#endif

        for (i=0; i<ldim1; i++) {
            v = VALUE1D(rank,f->current_step,i);
            if (r1[i] != v) {
                printE ("Error: a1[%d]=%d  !=  read=%d\n", i, v, r1[i]); 
                //goto endread;
            }
        }



        log ("  Check 2D variable a2_1...\n");
        adios_schedule_read (f, sel2_1, "a2_1",  0, 1, r2);
        adios_perform_reads (f, 1);

#ifdef DEBUG_PRINTS
        fprintf(stderr, "2D result :\n");
        n = 0;
        for (i=0; i<ldim1; i++) {
            fprintf (stderr, "row=%2d- = {", i);
            for (j=0; j<ldim2; j++) {
                fprintf (stderr, "%d ", r2[n]); 
                n++;
            }
            fprintf(stderr, "}\n");
        }
#endif

        n = 0;
        for (i=0; i<ldim1; i++) {
            for (j=0; j<ldim2; j++) {
                v = VALUE2D(rank,f->current_step,i,j);
                if (v != r2[n]) {
                    printE ("Error: a2_2[%d,%d]=%d  !=  read=%d\n", i, j, v, r2[n]);
                    //goto endread;
                }
                n++;
            }
        }


        log ("  Check 2D variable a2_2...\n");
        adios_schedule_read (f, sel2_2, "a2_2",  0, 1, r2);
        adios_perform_reads (f, 1);

#ifdef DEBUG_PRINTS
        fprintf(stderr, "2D result :\n");
        n = 0;
        for (i=0; i<ldim1; i++) {
            fprintf (stderr, "row=%2d- = {", i);
            for (j=0; j<ldim2; j++) {
                fprintf (stderr, "%d ", r2[n]);
                n++;
            }
            fprintf(stderr, "}\n");
        }
#endif

        n = 0;
        for (i=0; i<ldim1; i++) {
            for (j=0; j<ldim2; j++) {
                v = VALUE2D(rank,f->current_step,i,j);
                if (v != r2[n]) {
                    printE ("Error: a2_2[%d,%d]=%d  !=  read=%d\n", i, j, v, r2[n]);
                    //goto endread;
                }
                n++;
            }
        }



        log ("  Check 3D variable a3_1...\n");
        adios_schedule_read (f, sel3_1, "a3_1",  0, 1, r3);
        adios_perform_reads (f, 1);

#ifdef DEBUG_PRINTS
        fprintf(stderr, "3D selection :\n");
        n=0;
        for (i=0; i<ldim1; i++) {
            for (j=0; j<ldim2; j++) {
                fprintf (stderr, "[%d,%d] = {", i, j);
                for (k=0; k<ldim3; k++) {
                    fprintf (stderr, "%d ", r3[n]);
                    n++;
                }
                fprintf(stderr, "}\n");
            }
        }
#endif
        n = 0;
        for (i=0; i<ldim1; i++) {
            for (j=0; j<ldim2; j++) {
                for (k=0; k<ldim3; k++) {
                    v = VALUE3D(rank,f->current_step,i,j,k);
                    if (v != r3[n]) {
                        printE ("Error: a3_2[%d,%d,%d]=%d  !=  read=%d\n", i, j, k, v, r3[n]);
                        //goto endread;
                    }
                    n++;
                }
            }
        }



        log ("  Check 3D variable a3_2...\n");
        adios_schedule_read (f, sel3_2, "a3_2",  0, 1, r3);
        adios_perform_reads (f, 1);

#ifdef DEBUG_PRINTS
        fprintf(stderr, "3D selection :\n");
        n=0;
        for (i=0; i<ldim1; i++) {
            for (j=0; j<ldim2; j++) {
                fprintf (stderr, "[%d,%d] = {", i, j);
                for (k=0; k<ldim3; k++) {
                    fprintf (stderr, "%d ", r3[n]);
                    n++;
                }
                fprintf(stderr, "}\n");
            }
        }
#endif
        n = 0;
        for (i=0; i<ldim1; i++) {
            for (j=0; j<ldim2; j++) {
                for (k=0; k<ldim3; k++) {
                    v = VALUE3D(rank,f->current_step,i,j,k);
                    if (v != r3[n]) {
                        printE ("Error: a3_2[%d,%d,%d]=%d  !=  read=%d\n", i, j, k, v, r3[n]);
                        //goto endread;
                    }
                    n++;
                }
            }
        }



        log ("  Check 3D variable a3_3...\n");
        adios_schedule_read (f, sel3_3, "a3_3",  0, 1, r3);
        adios_perform_reads (f, 1);

#ifdef DEBUG_PRINTS
        fprintf(stderr, "3D selection :\n");
        n=0;
        for (i=0; i<ldim1; i++) {
            for (j=0; j<ldim2; j++) {
                fprintf (stderr, "[%d,%d] = {", i, j); 
                for (k=0; k<ldim3; k++) {
                    fprintf (stderr, "%d ", r3[n]); 
                    n++;
                }
                fprintf(stderr, "}\n");
            }
        }
#endif
        n = 0;
        for (i=0; i<ldim1; i++) {
            for (j=0; j<ldim2; j++) {
                for (k=0; k<ldim3; k++) {
                    v = VALUE3D(rank,f->current_step,i,j,k);
                    if (v != r3[n]) {
                        printE ("Error: a3_3[%d,%d,%d]=%d  !=  read=%d\n", i, j, k, v, r3[n]);
                        //goto endread;
                    }
                    n++;
                }
            }
        }



        if (n1 < NSTEPS)
        {
            adios_advance_step (f, 0, -1.0);
        }
    }


endread:

    adios_selection_delete (sel1);
    adios_selection_delete (sel2_1);
    adios_selection_delete (sel2_2);

    adios_read_close(f);
    MPI_Barrier (comm);
    return err;
}

