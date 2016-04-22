/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2015.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C test: append data to an existing file and check reading it back 
   Runs it multiple times to get a file with multiple timesteps over multiple code runs.
   Then checks correct reading of a scalar variable 
   (only the values from rank 0, and only it's block 0 writes.
   Read test is similar to execute "bpls -d append.bp value"
 */

#include <stdio.h>
#include <stdlib.h>   // malloc()
#include <string.h>
#include <sys/stat.h> // struct stat sb, stat()

#include "adios.h"
#include "adios_read.h"
#include "adios_types.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

static MPI_Comm    comm = MPI_COMM_WORLD;
static int  rank;
static int  size;
static char filename [256];



int write_data (int start_step, int nsteps, int nblocks) 
{
    int         i, j, step, block;

    float       value; // == 1.0*rank    + 0.1 * (start_step+step)   + 0.01 * block

    /* ADIOS variables declarations for matching gwrite_temperature.ch */
    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;


    strcpy (filename, "append.bp");

    adios_init_noxml (comm);
    adios_set_max_buffer_size (1);

    int64_t       m_adios_group;
    int64_t       m_adios_file;

    adios_declare_group (&m_adios_group, "append", "", adios_flag_yes);
    adios_select_method (m_adios_group, "MPI_LUSTRE", "", "");

    if (rank == 0) {
        adios_define_var (m_adios_group, "nproc" ,"", adios_integer ,0, 0, 0);
        adios_define_var (m_adios_group, "nblocks" ,"", adios_integer ,0, 0, 0);
        adios_define_var (m_adios_group, "nsteps" ,"", adios_integer ,0, 0, 0);
    }

    for (i=0;i<nblocks;i++) {
        adios_define_var (m_adios_group, "value" ,"", adios_real, "", "", "");
    }

    char mode[2] = "a"; 
    if (start_step == 0)
        mode[0] = 'w'; // very first time, create file

    printf ("Write data to %s\n", filename);
    for (step=0; step < nsteps; step++) {
        printf ("  Write data step %d\n", start_step+step);

        adios_open (&m_adios_file, "append", filename, mode, comm);
        mode[0] = 'a';

        adios_groupsize = nblocks * sizeof(float) + 4 + 4 + 4;

        adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);
        if (rank == 0) {
            adios_write(m_adios_file, "nproc", (void *) &size);
            adios_write(m_adios_file, "nblocks", (void *) &nblocks);
            adios_write(m_adios_file, "nsteps", (void *) &nsteps);
        }

        /* now we will write the data for each block */
        for (block=0;block<nblocks;block++) {
            value = 1.0*rank + 0.1*(start_step+step) + 0.01*block;
            adios_write(m_adios_file, "value", &value);
        }

        adios_close (m_adios_file);
    }

    MPI_Barrier (comm);

    adios_finalize (rank);

    return 0;
}

// call this function from a single process
int read_data (int nsteps, int nblocks)
{
    int retval = 0;
    ADIOS_FILE *f;
    printf ("Read data from %s\n", filename);
    adios_read_init_method(ADIOS_READ_METHOD_BP, comm, "verbose=2");
    f = adios_read_open_file (filename, ADIOS_READ_METHOD_BP, comm);
    if (f != NULL) 
    {
        // check number of steps first
        ADIOS_VARINFO *vi = adios_inq_var (f, "value");
        if (nsteps != vi->nsteps) 
        {
            printf ("ERROR: number of steps for variable 'value' in file %s is %d but expected %d\n",
                    filename, vi->nsteps, nsteps);
            retval = 1;
        } 
        else 
        {
            // read array of scalar variable 'value' over time
            // like what "bpls -d append.bp value" does
            // This only reads the first block from rank 0 over time...
            printf ("  Read 'value' over %d nsteps. Only the values written by rank 0, only block 0\n", nsteps);
            ADIOS_SELECTION *sel = NULL;
            float *v = (float *) malloc (nsteps * sizeof(float));
            adios_schedule_read (f, sel, "value", 0, nsteps, v);
            adios_perform_reads (f, 1);

            int j;
            float expv;
            for (j = 0; j < nsteps; j++) {
                // expected value = (rank) 0 + 0.1*step + (block) 0
                expv = j * 0.1;
                if (v[j] != expv) {
                    printf ("    ERROR: expected value for element %d of 'value' is %f but got %f\n", 
                            j, expv, v[j]);
                    retval = 1;
                }
            }

            adios_free_varinfo (vi);
            adios_selection_delete (sel);
            free (v);
        }
    } 
    else 
    {
        printf ("  File %s could not be opened: \n", filename, adios_errmsg());
        retval = 1;
    }
    adios_read_finalize_method(ADIOS_READ_METHOD_BP);
    return retval;
}

int main (int argc, char ** argv) 
{
    int nblocks = 1; // number of record-blocks written per process in one step
    int nsteps = 3;
    int retval = 0;

    MPI_Comm    comm = MPI_COMM_WORLD;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    strcpy (filename, "append.bp");

    write_data (0, nsteps, nblocks);
    write_data (nsteps, nsteps, nblocks);

    if (rank == 0)
    {
        retval = read_data(2*nsteps, nblocks);
    }


    MPI_Finalize ();
    return retval;
}
