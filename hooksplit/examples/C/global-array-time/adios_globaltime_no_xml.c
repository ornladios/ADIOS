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
#include <string.h>
#include "mpi.h"
#include "adios.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank, size, i;
	int         NX = 10, G, O; 
	double      t[NX];
	MPI_Comm    comm = MPI_COMM_WORLD;

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	int         adios_err, it;
	uint64_t    adios_groupsize, adios_totalsize;
	int64_t     adios_handle;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

        G = 2 * NX * size;

	for (i = 0; i < NX; i++)
	    t[i] = rank * NX + i;

	strcpy (filename, "adios_globaltime.bp");

	adios_init_noxml ();
        adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10);

        int64_t       m_adios_group;
        int64_t       m_adios_file;

        adios_declare_group (&m_adios_group, "restart", "iter", adios_flag_yes);
        adios_select_method (m_adios_group, "MPI", "", "");

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
                     ,"iter,NX", "G", "O");

        for (it =0; it < 5; it++) {

        	adios_open (&m_adios_file, "restart", filename, "a", &comm);
        	adios_groupsize = 4 + 4 + 4 + NX * 8
                        	+ 4 + 4 + 4 + NX * 8;
        	adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);

        	adios_write(m_adios_file, "NX", (void *) &NX);
        	adios_write(m_adios_file, "G", (void *) &G);
        	O = rank * 2 * NX;
        	adios_write(m_adios_file, "O", (void *) &O);
        	adios_write(m_adios_file, "temperature", t);

        	for (i = 0; i < NX; i++)
            		t[i]++;

        	adios_write(m_adios_file, "NX", (void *) &NX);
        	adios_write(m_adios_file, "G", (void *) &G);
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
