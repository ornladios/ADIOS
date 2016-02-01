/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write a global array from N processors with gwrite
 *
 * How to run: mpirun -np <N> adios_global
 * Output: adios_global.bp
 * ADIOS config file: adios_global.xml
 *
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "adios.h"
//#include "adios_mpi.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

void print_usage (int argc, char ** argv)
{
    printf ("%s [nsteps]\n"
            "    nsteps: Number of steps\n",
            argv[0]
           );
}

int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank, size, i;
	int         NX = 100, G, O; 
        int         nsteps = 5;
	double      t[NX];
	MPI_Comm    comm = MPI_COMM_WORLD;

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	int         it;
	uint64_t    adios_groupsize;
	//uint64_t    adios_totalsize;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

        if (argc > 1) {
            char *ep;
            errno = 0;
            long int n = strtol (argv[1], &ep, 10); 
            if (!errno) {
                if (*ep != '\0') {
                    printf ("Argument is not an integer: %s. Use default nsteps = %d\n", argv[1], nsteps);
                } else {
                    nsteps = n;
                }
            } else {
                printf ("Error when converting argument '%s': %s\n"
                        "Use default nsteps = %d\n",
                        argv[1], strerror(errno), nsteps);
            }
        }

        G = 2 * NX * size;

	strcpy (filename, "adios_globaltime.bp");

	adios_init_noxml (comm);
        adios_set_max_buffer_size (10);

        int64_t       m_adios_group;
        int64_t       m_adios_file;

        adios_declare_group (&m_adios_group, "restart", "", adios_flag_yes);
        adios_select_method (m_adios_group, "MPI", "", "");
        //adios_select_method (m_adios_group, "MPI_AGGREGATE", "num_ost=2;num_aggregators=2;aggregation_type=2;verbose=3", "");

        adios_define_var (m_adios_group, "NX"
                     ,"", adios_integer
                     ,0, 0, 0);
  
        adios_define_var (m_adios_group, "G"
                     ,"", adios_integer
                     ,0, 0, 0);

        /* have to define O and temperature as many times as we 
           write them within one step (twice) */
        for (it=0; it < 2; it++) {
            adios_define_var (m_adios_group, "O"
                         ,"", adios_integer
                         ,0, 0, 0);
    
            int64_t var_id = adios_define_var (m_adios_group, "temperature"
                         ,"", adios_double
                         ,"NX", "G", "O");
            adios_set_transform (var_id, "none");
        }

        for (it = 0; it < nsteps; it++) {

            for (i = 0; i < NX; i++)
                t[i] = rank + it*0.1 + 0.01;

        	adios_open (&m_adios_file, "restart", filename, (it==0 ? "w" : "a"), comm);
        	adios_groupsize = 4 + 4 + 4 + NX * 8
                        	+ 4 + 4 + 4 + NX * 8;
        	//adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);

        	adios_write(m_adios_file, "NX", (void *) &NX);
        	adios_write(m_adios_file, "G", (void *) &G);
        	O = rank * 2 * NX;
        	adios_write(m_adios_file, "O", (void *) &O);
        	adios_write(m_adios_file, "temperature", t);

        	for (i = 0; i < NX; i++)
            		t[i] += 0.01;

        	O = rank * 2 * NX + NX;
        	adios_write(m_adios_file, "O", (void *) &O);
        	adios_write(m_adios_file, "temperature", t);

        	adios_close (m_adios_file);
                MPI_Barrier (comm);
        }

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}
