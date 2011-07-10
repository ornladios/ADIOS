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
	int         L[] = {1,2,2,2};
	int			G[] = {4,2,2,2};
	int 		O[] = {0,0,0,0};
	MPI_Comm    comm = MPI_COMM_WORLD;

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	int         adios_err;
	uint64_t    adios_groupsize, adios_totalsize;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);
	int 		t[][][][] ={{{rank,rank},{rank,rank}},{{rank,rank},{rank,rank}}};

	strcpy (filename, "chaji_no_xml.bp");

	adios_init_noxml ();
	adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10);

	int64_t       m_adios_group;
	int64_t       m_adios_file;

	adios_declare_group (&m_adios_group, "restart", "", adios_flag_yes);
	adios_select_method (m_adios_group, "MPI", "", "");

	//int i=0;
	for(int i=0; i<2; i++){
		adios_define_var (m_adios_group, "L0","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "L1","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "L2","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "L3","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "G0","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "G1","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "G2","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "G3","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "O0","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "O1","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "O2","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "O3","", 2, 0, 0, 0);
		adios_define_var (m_adios_group, "temperature","", 2, "L0,L1,L2,L3", "G0,G1,G2,G3", "O0,O1,O2,O3");
		if(rank!=0)
			break;
	}



	adios_open (&m_adios_file, "restart", filename, "w", &comm);

	adios_groupsize = 12*4 + 32 * 4;
	if(rank==0)
		adios_groupsize *= 2;

	adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);

	adios_write(m_adios_file, "L0", (void *) &L[0]);
	adios_write(m_adios_file, "L1", (void *) &L[1]);
	adios_write(m_adios_file, "L2", (void *) &L[2]);
	adios_write(m_adios_file, "L3", (void *) &L[3]);
	adios_write(m_adios_file, "G0", (void *) &G[0]);
	adios_write(m_adios_file, "G1", (void *) &G[1]);
	adios_write(m_adios_file, "G2", (void *) &G[2]);
	adios_write(m_adios_file, "G3", (void *) &G[3]);
	adios_write(m_adios_file, "O1", (void *) &O[1]);
	adios_write(m_adios_file, "O2", (void *) &O[2]);
	adios_write(m_adios_file, "O3", (void *) &O[3]);


	if(rank==1){
		O[0]=1;
		adios_write(m_adios_file, "O0", (void *) &O[0]);
		adios_write(m_adios_file, "temperature", &t[0]);
	}else if(rank==2){
		O[0]=2;
		adios_write(m_adios_file, "O0", (void *) &O[0]);
		adios_write(m_adios_file, "temperature", &t[0]);
	}else if(rank==0){
		O[0]=0;
		adios_write(m_adios_file, "O0", (void *) &O[0]);
		adios_write(m_adios_file, "temperature", &t[0]);

		adios_write(m_adios_file, "L0", (void *) &L[0]);
		adios_write(m_adios_file, "L1", (void *) &L[1]);
		adios_write(m_adios_file, "L2", (void *) &L[2]);
		adios_write(m_adios_file, "L3", (void *) &L[3]);
		adios_write(m_adios_file, "G0", (void *) &G[0]);
		adios_write(m_adios_file, "G1", (void *) &G[1]);
		adios_write(m_adios_file, "G2", (void *) &G[2]);
		adios_write(m_adios_file, "G3", (void *) &G[3]);
		adios_write(m_adios_file, "O1", (void *) &O[1]);
		adios_write(m_adios_file, "O2", (void *) &O[2]);
		adios_write(m_adios_file, "O3", (void *) &O[3]);
		O[0]=3;
		adios_write(m_adios_file, "O0", (void *) &O[0]);
		adios_write(m_adios_file, "temperature", &t[0]);
	}

	adios_close (m_adios_file);

	MPI_Barrier (comm);

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}

