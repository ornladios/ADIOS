/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Query Example write part: 
 *   write a few global array from a single processor
 *   to be queried using query_vars.c
 *
 * How to run: write_vars
 * Output: vars.bp
 * ADIOS config file: None
 *
*/

/* This example will write out three 2D variables, P, V and T.
*/
#ifndef _NOMPI
#define _NOMPI
#endif

#include <stdio.h>
#include <string.h>
#include "adios.h"
#include "adios_types.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

const int NX = 5;
const int NY = 6;

double T[5][6] = { 
    {1.1, 1.2, 1.3, 1.4, 1.5, 1.6},
    {2.1, 2.2, 2.3, 2.4, 2.5, 2.6},
    {3.1, 3.2, 3.3, 3.4, 3.5, 3.6},
    {4.1, 4.2, 4.3, 4.4, 4.5, 4.6},
    {5.1, 5.2, 5.3, 5.4, 5.5, 5.6},
};

double P[5][6] = { 
    {41.1, 61.2, 81.3, 81.4, 91.5, 31.6},
    {42.1, 62.2, 82.3, 82.4, 92.5, 32.6},
    {43.1, 63.2, 83.3, 83.4, 93.5, 33.6},
    {44.1, 64.2, 84.3, 84.4, 94.5, 34.6},
    {45.1, 65.2, 85.3, 85.4, 95.5, 35.6},
};

double V[5][6] = { 
    {41.1, 41.2, 41.3, 41.4, 41.5, 41.6},
    {45.1, 45.2, 45.3, 45.4, 45.5, 45.6},
    {49.1, 49.2, 49.3, 49.4, 49.5, 49.6},
    {54.1, 54.2, 54.3, 54.4, 54.5, 54.6},
    {55.1, 55.2, 55.3, 55.4, 55.5, 55.6},
};

int main (int argc, char ** argv) 
{
	MPI_Comm    comm = 0; // dummy mpi 

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	uint64_t  adios_groupsize, adios_totalsize;
	int64_t   g;
	int64_t   f;
	int64_t   Tid, Pid, Vid; // variable IDs
	char dimstr[32];

	sprintf (dimstr, "%d,%d", NX, NY);

	adios_init_noxml (comm);
	adios_set_max_buffer_size (1);

	adios_declare_group (&g, "vars", "", adios_flag_yes);
	adios_select_method (g, "POSIX", "", "");

	Tid = adios_define_var (g, "T" ,"", adios_double, dimstr, dimstr, "0,0");
	adios_set_transform (Tid, "none");
	Pid = adios_define_var (g, "P" ,"", adios_double, dimstr, dimstr, "0,0");
	adios_set_transform (Pid, "none");
	Vid = adios_define_var (g, "V" ,"", adios_double, dimstr, dimstr, "0,0");
	adios_set_transform (Vid, "none");


	adios_open (&f, "vars", "vars.bp", "w", comm);
	adios_groupsize = 3*NX*NY*sizeof(double);
	adios_group_size (f, adios_groupsize, &adios_totalsize);
	adios_write (f, "T", T);
	adios_write (f, "P", P);
	adios_write (f, "V", V);
	adios_close (f);

	adios_finalize (0);
	return 0;
}
