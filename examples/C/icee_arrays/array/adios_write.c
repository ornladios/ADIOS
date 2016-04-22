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

#ifdef DMALLOC
#include "dmalloc.h"
#endif

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
    int    cm_port = 59999;
    char*  cm_remote_host = "localhost";
    int    cm_remote_port = 59997;
    int    max_client = 1;
    char   initstring [512];
    int    verbose_level = 3;
    char*  cm_transport = "TCP";
    int    interval_sec = 5;
	int    NX = 1000;
    int    is_passive = 0;

    char*  adios_write_method = "MPI";
    enum ADIOS_READ_METHOD adios_read_method = ADIOS_READ_METHOD_BP;

    while ((c = getopt (argc, argv, "h:p:s:t:m:w:r:v:T:i:n:P")) != -1)
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
        case 'm':
            max_client = atoi(optarg);
            break;
        case 'w':
            adios_write_method = optarg;
            break;
        case 'r':
            if (strcmp(optarg, "BP") == 0) {
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
        case 'i':
            interval_sec = atoi(optarg);
            break;
        case 'n':
            NX = atoi(optarg);
            break;
        case 'P':
            is_passive = 1;
            break;
        default:
            usage(basename(argv[0]));
            break;
        }
    }

	char        filename [256];
	int         rank, size, i;
	int         G, O; 
	double      *t = (double *) malloc(NX * sizeof(double));
	assert(t != NULL);
	MPI_Comm    comm = MPI_COMM_WORLD;

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	int         it;
	uint64_t    adios_groupsize, adios_totalsize;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

    setlocale(LC_NUMERIC, "");

    G = NX * size;
    if (rank==0) printf("NX = %d\n", NX);

	strcpy (filename, "adios_globaltime.bp");

    sprintf(initstring, "verbose=%d;cm_host=%s;cm_port=%d;max_client=%d;transport=%s;is_passive=%d", 
            verbose_level, cm_host, cm_port+rank, max_client, cm_transport, is_passive);

	adios_init_noxml (comm);
    adios_set_max_buffer_size (10);

    int64_t       m_adios_group;
    int64_t       m_adios_file;

    adios_declare_group (&m_adios_group, "restart", "", adios_flag_yes);
    adios_select_method (m_adios_group, adios_write_method, initstring, "");

    adios_define_var (m_adios_group, "NX"
                      ,"", adios_integer
                      ,0, 0, 0);
  
    adios_define_var (m_adios_group, "G"
                      ,"", adios_integer
                      ,0, 0, 0);

    adios_define_var (m_adios_group, "O"
                      ,"", adios_integer
                      ,0, 0, 0);
    
    adios_define_var (m_adios_group, "temperature"
                      ,"", adios_double
                      ,"NX", "G", "O");

    for (it =0; it < 5; it++) 
    {

        for (i = 0; i < NX; i++)
            t[i] = rank + it + 1.0;

        MPI_Barrier(comm); 
        double t_start = MPI_Wtime();

        adios_open (&m_adios_file, "restart", filename, "a", comm);
        adios_groupsize = 4 + 4 + 4 + NX * 8;
        adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);

        adios_write(m_adios_file, "NX", (void *) &NX);
        adios_write(m_adios_file, "G", (void *) &G);
        O = rank * NX;
        adios_write(m_adios_file, "O", (void *) &O);
        adios_write(m_adios_file, "temperature", t);

        adios_close (m_adios_file);

        MPI_Barrier(comm);
        double t_end = MPI_Wtime();
        double t_elap = t_end - t_start;

        if (rank==0)
            printf("[%d] Elapsed %.03f seconds, throughput %'.03f KB/sec\n", 
                   it, t_elap, (double)adios_groupsize/t_elap/1024.0);
        
        sleep(interval_sec);
    }

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}
