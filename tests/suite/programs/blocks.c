/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: 
 *  Write a global array over time, with multiple blocks per process
 *    Similar to examples/C/global-array-time/adios_global_time_noxml.c
 *  Then open for reading and check if the blockinfo for each block is correct. 
 *  Do the reading twice, once with opening as file (all steps at once) and 
 *     once as streaming (step-by-step)
 *
 * How to run: mpirun -np <N> blocks
 * Output: blocks.bp
 * Exit code: the number of errors found (0=OK)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adios.h"
#include "adios_read.h"
#include "adios_error.h"

const static char fname[] = "blocks.bp";
const static MPI_Comm  comm = MPI_COMM_WORLD;
static int rank, size;
static int nerrors = 0;

int write_blocks ();
void print_written_info();
int read_all ();
int read_stepbystep ();

/* Remember (on rank 0) what was written (from all process) to check against it at reading */
static int nblocks_per_step;
static int nsteps;
static int * block_offset;  // block_offset[ step*nblocks_per_step + i ] is i-th block offset written in "step".
static int * block_count;   // block_count [ step*nblocks_per_step + i ] is i-th block size written in "step".


int main (int argc, char ** argv) 
{
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    nerrors = 0;
    write_blocks();
    if (!rank) {
        print_written_info(); // this is just for debug to check if rank 0 stores the correct values
        read_all();
        read_stepbystep();
    }

    MPI_Barrier (comm);
    MPI_Finalize ();
    free (block_offset);
    free (block_count);
    return nerrors;
}

int write_blocks () 
{
    int         NX, G, O; 
    double      *t;
    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    int         it, i, r;
    uint64_t    adios_groupsize, adios_totalsize;

    if (!rank) printf ("------- Write blocks -------\n");
    // We will have "3 steps * 2 blocks per process * number of processes" blocks
    nsteps = 3;
    nblocks_per_step = 2;
    block_offset = (int*) malloc (sizeof(int) * nsteps * nblocks_per_step * size);
    block_count  = (int*) malloc (sizeof(int) * nsteps * nblocks_per_step * size);

    adios_init_noxml (comm);
    adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10);

    int64_t       m_adios_group;
    int64_t       m_adios_file;

    adios_declare_group (&m_adios_group, "restart", "", adios_flag_yes);
    adios_select_method (m_adios_group, "MPI", "", "");

    adios_define_var (m_adios_group, "NX"
            ,"", adios_integer
            ,0, 0, 0);

    adios_define_var (m_adios_group, "G"
            ,"", adios_integer
            ,0, 0, 0);

    /* have to define O and temperature as many times as we 
       write them within one step (twice) */
    for (it=0; it < nblocks_per_step; it++) {
        adios_define_var (m_adios_group, "O"
                ,"", adios_integer
                ,0, 0, 0);

        adios_define_var (m_adios_group, "t"
                ,"", adios_double
                ,"NX", "G", "O");
    }

    for (it =0; it < nsteps; it++) {
        if (!rank) printf ("Step %d:\n", it);
        NX = 10+it;
        G = nblocks_per_step * NX * size;

        t = (double *) malloc (NX*sizeof(double));

        for (i = 0; i < NX; i++)
            t[i] = rank + it*0.1 + 0.01;

        MPI_Barrier (comm);
        adios_open (&m_adios_file, "restart", fname, "a", comm);
        adios_groupsize = 4 + 4 + 4 + NX * 8
            + 4 + 4 + 4 + NX * 8;
        adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);

        adios_write(m_adios_file, "NX", (void *) &NX);
        adios_write(m_adios_file, "G", (void *) &G);
        O = rank * nblocks_per_step * NX;
        adios_write(m_adios_file, "O", (void *) &O);
        adios_write(m_adios_file, "t", t);

        printf ("rank %d: block 1: size=%d, offset=%d\n", rank, NX, O);
        for (r = 0; r < size; r++) {
            block_count  [it*nblocks_per_step*size + nblocks_per_step*r] = NX; 
            block_offset [it*nblocks_per_step*size + nblocks_per_step*r] = r * nblocks_per_step * NX; 
        }

        for (i = 0; i < NX; i++)
            t[i] += 0.01;

        O = rank * nblocks_per_step * NX + NX;
        adios_write(m_adios_file, "O", (void *) &O);
        adios_write(m_adios_file, "t", t);

        printf ("rank %d: block 2: size=%d, offset=%d\n", rank, NX, O);
        for (r = 0; r < size; r++) {
            block_count  [it*nblocks_per_step*size + nblocks_per_step*r + 1] = NX; 
            block_offset [it*nblocks_per_step*size + nblocks_per_step*r + 1] = r * nblocks_per_step * NX + NX; 
        }

        adios_close (m_adios_file);
        MPI_Barrier (comm);

        free(t);
    }

    adios_finalize (rank);

    return 0;
}

void print_written_info()
{
    int s, r, b;
    printf ("------- Information on rank 0 --------\n");
    for (s = 0; s < nsteps; s++) {
        printf ("Step %d:\n", s);
        for (r = 0; r < size; r++) {
            for (b = 0; b < nblocks_per_step; b++) {
                printf ("rank %d: block %d: size=%d, offset=%d\n", r, b+1, 
                        block_count  [s*nblocks_per_step*size + nblocks_per_step*r + b],
                        block_offset [s*nblocks_per_step*size + nblocks_per_step*r + b]
                       );
            }
        }
    }
}

int print_varinfo (ADIOS_FILE *f, int start_step) 
{
    ADIOS_VARINFO * v;
    int i,j,k;

    v = adios_inq_var (f, "t");
    adios_inq_var_blockinfo (f, v);

    printf ("ndim = %d\n",  v->ndim);
    printf ("dims[%llu]\n",  v->dims[0]);
    printf ("nsteps = %d\n",  v->nsteps);
    printf ("sum_nblocks = %d\n",  v->sum_nblocks);
    for (i = 0; i < v->nsteps; i++) {
        printf ("  nblocks[%d] = %d\n", i, v->nblocks[i]);
        for (j = 0; j < v->nblocks[i]; j++) {
            printf("    block %2d: [%lld:%lld]", j,
                        v->blockinfo[j].start[0],
                        v->blockinfo[j].start[0] + v->blockinfo[j].count[0]-1);
            
            if (v->blockinfo[j].start[0] != block_offset [(start_step+i)*nblocks_per_step*size + j] ||
                v->blockinfo[j].count[0] != block_count  [(start_step+i)*nblocks_per_step*size + j] ) 
            {
                nerrors++;
                printf ("\tERROR: expected [%lld:%lld]",
                    block_offset [(start_step+i)*nblocks_per_step*size + j],
                    block_offset [(start_step+i)*nblocks_per_step*size + j] + 
                      block_count  [(start_step+i)*nblocks_per_step*size + j] -1
                );
            }
            printf("]\n");
        }
    }
    adios_free_varinfo (v);
}

int read_all ()
{
    ADIOS_FILE * f;
    float timeout_sec = 0.0; 
    int steps = 0;
    int retval = 0;
    MPI_Comm    comm = MPI_COMM_SELF;

    adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=3");
    printf ("--------- Read as file  ------------\n");
    f = adios_read_open_file (fname, ADIOS_READ_METHOD_BP, comm);
    if (f == NULL) {
        printf ("Error at opening file: %s\n", adios_errmsg());
        retval = adios_errno;
    }
    else
    {
        /* Processing all the steps at once */
        print_varinfo (f, 0);
        adios_read_close (f);
    }
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    return retval;
}

int read_stepbystep ()
{
    ADIOS_FILE * f;
    float timeout_sec = 0.0; 
    int steps = 0;
    int retval = 0;
    MPI_Comm    comm = MPI_COMM_SELF;

    adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=3");
    printf ("--------- Read as stream  ------------\n");
    f = adios_read_open (fname, ADIOS_READ_METHOD_BP,
                          comm, ADIOS_LOCKMODE_NONE, timeout_sec);
    if (adios_errno == err_file_not_found)
    {
        printf ("Stream not found after waiting %f seconds: %s\n",
                timeout_sec, adios_errmsg());
        retval = adios_errno;
    }
    else if (adios_errno == err_end_of_stream)
    {
        printf ("Stream terminated before open. %s\n", adios_errmsg());
        retval = adios_errno;
    }
    else if (f == NULL) {
        printf ("Error at opening stream: %s\n", adios_errmsg());
        retval = adios_errno;
    }
    else
    {
        /* Processing loop over the steps (we are already in the first one) */
        while (adios_errno != err_end_of_stream) {
            steps++; // steps start counting from 1
            printf ("Step: %d\n", f->current_step);
            print_varinfo (f, f->current_step);

            // advance to 1) next available step with 2) blocking wait
            adios_advance_step (f, 0, timeout_sec);
            if (adios_errno == err_step_notready)
            {
                //printf ("No new step arrived within the timeout. Quit. %s\n",
                //        adios_errmsg());
                break; // quit while loop
            }
        }
        adios_read_close (f);
    }
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    //printf ("We have processed %d steps\n", steps);
    return retval;
}
