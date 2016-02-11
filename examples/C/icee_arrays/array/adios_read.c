/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS ICEE Example
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> 
#include <libgen.h>
#include <locale.h>
#include "mpi.h"
#include "adios.h"
#include "adios_read.h"
#include "adios_error.h"

void usage(const char *argv0)
{
    fprintf(stderr, "usage: %s\n", argv0);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "-w [ MPI | ICEE ]\n");
    fprintf(stderr, "-r [ BP | ICEE ]\n");
    exit (1);
}

int main (int argc, char ** argv) 
{
    int c;
    opterr = 0;

    char*  cm_host = "localhost";
    int    cm_port = 59997;
    char*  cm_remote_host = "localhost";
    int    cm_remote_port = 59999;
    int    is_multi_writers = 0;
    char*  remote_list = "";
    char*  attr_list = "";
    char   initstring [512];
    int    verbose_level = 3;
    char*  cm_transport = "TCP";
    float  timeout_sec = 10.0; 
    int    use_native_contact = 0;
    int    is_passive = 0;

    char*  adios_write_method = "MPI";
    enum ADIOS_READ_METHOD adios_read_method = ADIOS_READ_METHOD_BP;

    while ((c = getopt (argc, argv, "h:p:s:t:u:a:w:r:v:T:o:nP")) != -1)
    {
        switch (c)
        {
        case 'h':
            cm_host = optarg;
            break;
        case 'p':
            cm_port = atoi(optarg);
            break;
        case 's':
            cm_remote_host = optarg;
            break;
        case 't':
            cm_remote_port = atoi(optarg);
            break;
        case 'u':
            is_multi_writers = 1;
            remote_list = optarg;
            break;
        case 'a':
            is_multi_writers = 1;
            attr_list = optarg;
            break;
        case 'w':
            adios_write_method = optarg;
            break;
        case 'r':
            if (strcmp(optarg, "BP") ==0) {
                adios_read_method = ADIOS_READ_METHOD_BP;
            } else if (strcmp(optarg, "ICEE") == 0) {
                adios_read_method = ADIOS_READ_METHOD_ICEE;
            } else {
                fprintf(stderr, "No read method: %s\n", optarg);
            }
            break;
        case 'v':
            verbose_level = atoi(optarg);
            break;
        case 'T':
            cm_transport = optarg;
            break;
        case 'o':
            timeout_sec = atoi(optarg);
            break;
        case 'n':
            use_native_contact = 1;
            break;
        case 'P':
            is_passive = 1;
            break;
        default:
            usage(basename(argv[0]));
            break;
        }
    }

    int         rank, size, i;
    MPI_Comm    comm = MPI_COMM_WORLD;
    ADIOS_FILE * f;
    ADIOS_VARINFO * v;
    ADIOS_SELECTION * sel;

    int steps = 0;
    int retval = 0;

    void * data = NULL;
    uint64_t start[2], count[2];

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);


    if (!is_multi_writers)
        sprintf(initstring, "verbose=%d;cm_host=%s;cm_port=%d;cm_remote_host=%s;cm_remote_port=%d;transport=%s;use_native_contact=%d;is_passive=%d", 
                verbose_level, cm_host, cm_port+rank, cm_remote_host, cm_remote_port, cm_transport, use_native_contact, is_passive);
    else
        sprintf(initstring, "verbose=%d;cm_host=%s;cm_port=%d;remote_list=%s;attr_list=%s;transport=%s;use_native_contact=%d;is_passive=%d", 
                verbose_level, cm_host, cm_port+rank, remote_list, attr_list, cm_transport, use_native_contact, is_passive);
        

    adios_read_init_method (adios_read_method, comm, initstring);

    f = adios_read_open ("adios_globaltime.bp", adios_read_method,
                         comm, ADIOS_LOCKMODE_NONE, timeout_sec);
    if (adios_errno == err_file_not_found)
    {
        printf ("rank %d: Stream not found after waiting %f seconds: %s\n",
                rank, timeout_sec, adios_errmsg());
        retval = adios_errno;
    }
    else if (adios_errno == err_end_of_stream)
    {
        printf ("rank %d: Stream terminated before open. %s\n", rank, adios_errmsg());
        retval = adios_errno;
    }
    else if (f == NULL) {
        printf ("rank %d: Error at opening stream: %s\n", rank, adios_errmsg());
        retval = adios_errno;
    }
    else
    {
        /* process file here... */
        v = adios_inq_var (f, "temperature");

        uint64_t slice_size = v->dims[0]/size;
        if (rank == size-1)
            slice_size = slice_size + v->dims[0]%size;

        start[0] = rank * slice_size;
        count[0] = slice_size;

        data = malloc (slice_size * sizeof(double));
	assert(data != NULL);

        /* Processing loop over the steps (we are already in the first one) */
        while (adios_errno != err_end_of_stream) {
            steps++; // steps start counting from 1

            sel = adios_selection_boundingbox (v->ndim, start, count);
            adios_schedule_read (f, sel, "temperature", 0, 1, data);
            adios_perform_reads (f, 1);

            double sum = 0.0;
            for (i = 0; i < slice_size; i++) 
            {
                sum += *((double *)data + i);
            }

            printf("Step:%d, rank=%d: sum(data[%" PRIu64 ":%" PRId64 "]) = %.01f\n",
                   f->current_step, rank, start[0], start[0]+count[0]-1, sum);

            // advance to 1) next available step with 2) blocking wait
            adios_advance_step (f, 0, timeout_sec);
            if (adios_errno == err_step_notready)
            {
                printf ("rank %d: No new step arrived within the timeout. Quit. %s\n",
                        rank, adios_errmsg());
                break; // quit while loop
            }
            
        }
        
        adios_read_close (f);
    }
    
    if (rank==0) 
        printf ("We have processed %d steps\n", steps);
    
    adios_read_finalize_method (adios_read_method);
    free (data);
    MPI_Finalize ();
    
    return retval;
}

