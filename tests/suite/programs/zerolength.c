/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: Test writing/reading 0 length blocks with/without transforms
 *  Write a global array over time, process rank 0 writes 0 size in the global space.
 *  Do this with a transformation.
 *  Then open for reading and try to read it back. 
 *
 * How to run: mpirun -np <N> zerolength
 * Output: zerolength.bp
 * Exit code: the number of errors found (0=OK)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include "adios.h"
#include "adios_read.h"
#include "adios_error.h"
#include "adios_transform_methods.h"

static char fname[256];
const static MPI_Comm  comm = MPI_COMM_WORLD;
static int rank, size;
static int nerrors = 0;

int write_data (char * transformdef);
void print_written_info();
int read_all ();
int read_stepbystep ();
int read_scalar ();
int read_scalar_stepbystep ();

/* Remember (on rank 0) what was written (from all process) to check against it at reading */
static int nsteps = 2;
static uint64_t * block_offset;  // block_offset[ step*size + i ] is i-th block offset written in "step".
static uint64_t * block_count;   // block_count [ step*size + i ] is i-th block size written in "step".
static uint64_t * gdims;  // gdims[i] is the global dimension in i-th "step".

#define VALUE(step, rank, i) (step + rank*0.1 + 0.001*(i%100))

int main (int argc, char ** argv) 
{
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
    adios_init_noxml (comm);
    adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=3");
    nerrors = 0;

    ADIOS_AVAILABLE_TRANSFORM_METHODS * transforms = adios_available_transform_methods();
    if (transforms) {
        int i;
        for (i=0; i<transforms->ntransforms; i++)
        {
            if (!rank)
                printf("=======================================================\n"
                       "Test with transform \"%s\"\t: %s\n",  transforms->name[i],
                        transforms->description[i]);

            char transformstr[256];
            if (!strcmp(transforms->name[i], "zfp"))
            {
                snprintf(transformstr, sizeof(transformstr), "zfp:accuracy=0.0001");
            }
            else if (!strcmp(transforms->name[i], "sz"))
            {
                snprintf(transformstr, sizeof(transformstr), "sz:absolute=0.0001");
            }
            else
            {
                snprintf(transformstr, sizeof(transformstr), "%s", transforms->name[i]);
            }

            snprintf(fname, sizeof(fname), "zerolength_%s.bp", transforms->name[i]);
            MPI_Barrier (comm);
            write_data(transformstr);
            if (!rank) {
                print_written_info(); // this is just for debug to check if rank 0 stores the correct values
            }
            nerrors += read_all();
            MPI_Barrier (comm);
            sleep(1);
        }
        adios_available_transform_methods_free(transforms);
    }

    MPI_Barrier (comm);
    adios_read_finalize_method (ADIOS_READ_METHOD_BP);
    adios_finalize (rank);
    int total_errors;
    MPI_Allreduce(&nerrors, &total_errors, 1, MPI_INT, MPI_SUM, comm);
    MPI_Finalize ();
    free (block_offset);
    free (block_count);
    free (gdims);
    if (!rank) printf ("----------- Done. Found %d errors -------\n", total_errors);
    return total_errors;
}

int write_data (char * transformdef)
{
    int         NX, G, O; 
    double      *t;
    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    int         it, i, r;
    uint64_t    adios_groupsize, adios_totalsize;

    if (!rank) printf ("------- Write blocks to file %s -------\n", fname);
    // We will have "nsteps * number of processes" blocks
    block_offset = (uint64_t*) malloc (sizeof(uint64_t) * nsteps * size);
    block_count  = (uint64_t*) malloc (sizeof(uint64_t) * nsteps * size);
    gdims        = (uint64_t*) malloc (sizeof(uint64_t) * nsteps);


    int64_t       m_adios_group;
    int64_t       m_adios_file;
    int64_t       var_t;

    adios_declare_group (&m_adios_group, transformdef, "", adios_stat_default);
    adios_select_method (m_adios_group, "MPI_AGGREGATE", "num_aggregators=2;num_ost=2", "");

    adios_define_var (m_adios_group, "NX", "", adios_integer, 0, 0, 0); 
    adios_define_var (m_adios_group, "G", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "O", "", adios_integer, 0, 0, 0);
    var_t = adios_define_var (m_adios_group, "t", "", adios_double, "NX", "G", "O");

    adios_set_transform (var_t, transformdef);

    for (it =0; it < nsteps; it++) {

        if (!rank) printf ("Step %d:\n", it);

        NX = 100; // we will change this for rank 0 below
        t = (double *) malloc (NX*sizeof(double));

        for (i = 0; i < NX; i++)
            t[i] = VALUE(it, rank, i);

        G = NX * (size-1);

        // first block is 0 length
        block_count [it*size] = 0;
        block_offset [it*size] = 0;
        for (r = 1; r < size; r++) {
            block_count  [it*size + r] = NX;
            block_offset [it*size + r] = (r-1) * NX;
        }

        if (!rank) {
            NX = 0;
            O = 0;
        } else {
            O = (rank-1) * NX;
        }
        gdims [it] = G;


        printf ("rank %d: size=%d, offset=%d\n", rank, NX, O);
        MPI_Barrier (comm);
        if (it==0)
            adios_open (&m_adios_file, transformdef, fname, "w", comm);
        else
            adios_open (&m_adios_file, transformdef, fname, "a", comm);

        adios_write(m_adios_file, "NX", (void *) &NX);
        adios_write(m_adios_file, "G", (void *) &G);
        adios_write(m_adios_file, "O", (void *) &O);
        adios_write(m_adios_file, "t", t);

        adios_close (m_adios_file);
        MPI_Barrier (comm);

        free(t);
    }


    return 0;
}

void print_written_info()
{
    int s, r, b;
    printf ("\n------- Information recorded on rank 0 (read will compare to this info)  --------\n");
    for (s = 0; s < nsteps; s++) {
        printf ("Step %d:\n", s);
        printf ("  Global dim = %" PRIu64 "\n", gdims[s]);
        for (r = 0; r < size; r++) {
            printf ("  rank %d: size=%" PRIu64 ", offset=%" PRIu64 "\n",
                    r,
                    block_count  [s*size + r],
                    block_offset [s*size + r]
            );
        }
    }
}

int close_enough(double v1, double v2, double limit)
{
    if (v1-v2 < limit && v2-v1 < limit)
        return 1;
    return 0;
}

int read_data (ADIOS_FILE *f, int step) //uint64_t count, uint64_t offset)
{
    int i;
    ADIOS_SELECTION *sel;
    uint64_t count  = block_count [step*size+rank]; 
    uint64_t offset = block_offset [step*size+rank];
    printf ("rank %d bounding box = %" PRIu64 
            " elements from offset %" PRIu64 "\n", 
            rank, count, offset);

    double *t = (double *) malloc ((count+1)*sizeof(double)); 

    sel = adios_selection_boundingbox (1, &offset, &count);
    adios_schedule_read (f, sel, "t", step, 1, t);
    adios_perform_reads (f, 1);

    for (i = 0; i < count; i++)
    {
        if (!close_enough(t[i], VALUE(step, rank, i), 0.0001))
        {
           printf ("ERROR: step %d rank %d element %d expected value %9.6f but is %9.6f\n",
                   step, rank, i, VALUE(step, rank, i), t[i]);
           return 1;
        }
    }

    adios_selection_delete (sel);
    return 0;
}


int read_all ()
{
    ADIOS_FILE * f;
    int steps = 0;
    int retval = 0;


    if (!rank)
        printf ("\n--------- Read as file %s  ------------\n", fname);
    f = adios_read_open_file (fname, ADIOS_READ_METHOD_BP, comm);
    if (f == NULL) {
        printf ("Error at opening file: %s\n", adios_errmsg());
        retval = adios_errno;
    }
    else
    {
        int s;
        for (s = 0; s < nsteps; s++) {
            if (!rank) printf ("Step %d:\n", s);
            int err = read_data (f, s);
            MPI_Barrier(comm);
            MPI_Allreduce(&err, &retval, 1, MPI_INT, MPI_MAX, comm);
            if (!rank) printf ("\n");
            if (retval)
                break;
        }
        adios_read_close (f);
    }
    return retval;
}

