/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Query Example write part: 
 *   write a few global array from a single processor
 *   to be queried using query_table.c
 *
 * How to run: write_table
 * Output: table.bp
 * ADIOS config file: None
 *
*/

/* This example will write out a 2D table like particle data,
   where each row is data of one particle, and 
   each column is a property of the particles.
   It's completely fake data with integer values for demonstration.
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

/* table for Atoms, named A.
            unique ID
            element type (C, O, N etc)
            potential   
            kinetic energy  
            position X,Y and Z 

   Note: This is not some data with any physical meaning...
*/

const int NX = 10;
const int NY = 7;

int32_t A[10][7] = { 
/*   id, Elem,   P,   KE,   x,  y,  z */
    {0,   0,     97,  15,   8,  0,  7},
    {1,   0,     96,  16,   5, -1,  7},
    {2,   1,     32,   1,   7,  1,  7},
    {3,   2,    200,   8,   6,  0,  7},

    {4,   0,     94,  14,   1,  3,  5},
    {5,   0,     96,  13,   2,  4,  5},
    {6,   2,    201,   2,   3,  3,  5},
    {7,   1,     37,   9,   4,  2,  5},
    {8,   0,     98,  15,   5,  1,  5},
    {9,   0,     99,  14,   6,  2,  5},
};

const int Columns_length = 11; /* length of string + 1 for the terminating 0 */
char Columns[7][11] = {
    "ID        ",
    "Element   ",
    "Potential ",
    "Kinetic E ",
    "Position X",
    "Position Y",
    "Position Z",
};

const int n_of_elements = 3;
const int Elements_length = 9; /* length of string + 1 for the terminating 0 */
char Elements[3][9] = {
    "Carbon  ",
    "Nitrogen",
    "Oxygen  "
};


int main (int argc, char ** argv) 
{
	MPI_Comm    comm = 0; // dummy mpi 

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	uint64_t  adios_groupsize, adios_totalsize;
	int64_t   g;
	int64_t   f;
	char dimstr[32];

	adios_init_noxml (comm);
	adios_set_max_buffer_size (1);

	adios_declare_group (&g, "table", "", adios_flag_yes);
	adios_select_method (g, "POSIX", "", "");

	sprintf (dimstr, "%d,%d", NX, NY);
	adios_define_var (g, "A" ,"", adios_integer, dimstr, dimstr, "0,0");
	sprintf (dimstr, "%d,%d", n_of_elements, Elements_length);
	adios_define_var (g, "Elements" ,"", adios_byte, dimstr, dimstr, "0,0");
	sprintf (dimstr, "%d,%d", NY, Columns_length);
	adios_define_var (g, "Columns" ,"", adios_byte, dimstr, dimstr, "0,0");


	adios_open (&f, "table", "table.bp", "w", comm);

	adios_groupsize = NX*NY*sizeof(int32_t)           /* size of A */
	+ n_of_elements * Elements_length /* size of Elements */
	+ NY * Columns_length;            /* size of Columns */

	adios_group_size (f, adios_groupsize, &adios_totalsize);
	adios_write (f, "A", A);
	adios_write (f, "Elements", Elements);
	adios_write (f, "Columns", Columns);
	adios_close (f);

	adios_finalize (0);
	return 0;
}
