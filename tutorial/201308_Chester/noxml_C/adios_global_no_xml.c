/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write a global array from N processors with gwrite
 *
 * How to run: mpirun -np <N> adios_global_no_xml
 * Output: adios_global_no_xml.bp
 * ADIOS config file: None
 *
*/

/* This example will write out 2 sub blocks of the variable temperature
   and place these in the global array.
   This example illustrates both the use of sub blocks in writing, and
   the usage of the ADIOS non-xml API's
*/

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"
#include "adios_types.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank, size, i, block;
	int         NX = 100, Global_bounds, Offsets; 
	double      t[NX], b[NX];
        int         sub_blocks = 3;
	MPI_Comm    comm = MPI_COMM_WORLD;

	/* ADIOS variables */
	uint64_t    adios_groupsize, adios_totalsize;
        int64_t     m_adios_group;
        int64_t     m_adios_file;


	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

        Global_bounds = sub_blocks * NX * size;

	strcpy (filename, "adios_global_no_xml.bp");

	adios_init_noxml (comm);
        adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10);

        adios_declare_group (&m_adios_group, "restart", "iter", adios_flag_yes);
        adios_select_method (m_adios_group, "MPI", "", "");


        adios_define_var (m_adios_group, "NX"
			,"", adios_integer
			,0, 0, 0);
   
	adios_define_var (m_adios_group, "Global_bounds"
			,"", adios_integer
			,0, 0, 0);

        for (i=0;i<sub_blocks;i++) {
   
           adios_define_var (m_adios_group, "Offsets"
                        ,"", adios_integer
                        ,0, 0, 0);
   
           adios_define_var (m_adios_group, "temperature"
                        ,"", adios_double
                        ,"NX", "Global_bounds", "Offsets");

           adios_define_var (m_adios_group, "blocks"
                        ,"", adios_double
                        ,"NX", "Global_bounds", "Offsets");
        }
   
        adios_open (&m_adios_file, "restart", filename, "w", comm);

        adios_groupsize = sub_blocks * (4 + 4 + 4 + 2 * NX * 8);

        adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);
	adios_write(m_adios_file, "NX", (void *) &NX);
	adios_write(m_adios_file, "Global_bounds", (void *) &Global_bounds);
/* now we will write the data for each sub block */
        for (block=0;block<sub_blocks;block++) {

           Offsets = rank * sub_blocks * NX + block*NX;
           adios_write(m_adios_file, "Offsets", (void *) &Offsets);

           for (i = 0; i < NX; i++) {
               t[i] = Offsets + i;
               b[i] = rank * NX * sub_blocks + block*50+ 1;
           }

           adios_write(m_adios_file, "temperature", t);
           adios_write(m_adios_file, "blocks", b);
        }

        adios_close (m_adios_file);

        MPI_Barrier (comm);

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}
