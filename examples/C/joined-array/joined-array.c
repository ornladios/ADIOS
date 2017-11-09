/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write a 2D table as local array from N processors and
 * let ADIOS present it at reading as a joined global 2D array
 *
 * How to run: mpirun -np <N> joined-array
 * Output: joined-array.bp
 *
*/

/* This example will write out a 2D variable called 'table'. Processes don't
   define the offsets of their pieces in a global 2D array, instead they let
   ADIOS merge the tables at the slow dimension to present one global 2D table 
   at read time.

   Running on 4 processes the example output is the following:

$ mpirun -np 4 ./joined-array
$ bpls -la joined-array.bp -d table
  double   table  {18, 6} = 0 / 3.5 / null  / null 
    ( 0,0)    0 0.1 0.2 0.3 0.4 0.5
    ( 1,0)    0 0.1 0.2 0.3 0.4 0.5
    ( 2,0)    0 0.1 0.2 0.3 0.4 0.5
    ( 3,0)    1 1.1 1.2 1.3 1.4 1.5
    ( 4,0)    1 1.1 1.2 1.3 1.4 1.5
    ( 5,0)    1 1.1 1.2 1.3 1.4 1.5
    ( 6,0)    1 1.1 1.2 1.3 1.4 1.5
    ( 7,0)    2 2.1 2.2 2.3 2.4 2.5
    ( 8,0)    2 2.1 2.2 2.3 2.4 2.5
    ( 9,0)    2 2.1 2.2 2.3 2.4 2.5
    (10,0)    2 2.1 2.2 2.3 2.4 2.5
    (11,0)    2 2.1 2.2 2.3 2.4 2.5
    (12,0)    3 3.1 3.2 3.3 3.4 3.5
    (13,0)    3 3.1 3.2 3.3 3.4 3.5
    (14,0)    3 3.1 3.2 3.3 3.4 3.5
    (15,0)    3 3.1 3.2 3.3 3.4 3.5
    (16,0)    3 3.1 3.2 3.3 3.4 3.5
    (17,0)    3 3.1 3.2 3.3 3.4 3.5


   The details (per process data) are the following:
$ bpls -la joined-array.bp -dD cols
  integer  cols   scalar = 6
        step 0: 4 instances available
               6 6 6 6
$ bpls -la joined-array.bp -dD rows
  integer  rows   scalar = 3
        step 0: 4 instances available
               3 4 5 6
$ bpls -la joined-array.bp -dD table
  double   table  {18, 6} = 0 / 3.5 / null  / null 
        step 0: 
          block 0: [ 0: 2, 0:5] = 0 / 0.5/ N/A / N/A 
    (0,0)    0 0.1 0.2 0.3 0.4 0.5
    (1,0)    0 0.1 0.2 0.3 0.4 0.5
    (2,0)    0 0.1 0.2 0.3 0.4 0.5
          block 1: [ 3: 6, 0:5] = 1 / 1.5/ N/A / N/A 
    (0,0)    1 1.1 1.2 1.3 1.4 1.5
    (1,0)    1 1.1 1.2 1.3 1.4 1.5
    (2,0)    1 1.1 1.2 1.3 1.4 1.5
    (3,0)    1 1.1 1.2 1.3 1.4 1.5
          block 2: [ 7:11, 0:5] = 2 / 2.5/ N/A / N/A 
    (0,0)    2 2.1 2.2 2.3 2.4 2.5
    (1,0)    2 2.1 2.2 2.3 2.4 2.5
    (2,0)    2 2.1 2.2 2.3 2.4 2.5
    (3,0)    2 2.1 2.2 2.3 2.4 2.5
    (4,0)    2 2.1 2.2 2.3 2.4 2.5
          block 3: [12:17, 0:5] = 3 / 3.5/ N/A / N/A 
    (0,0)    3 3.1 3.2 3.3 3.4 3.5
    (1,0)    3 3.1 3.2 3.3 3.4 3.5
    (2,0)    3 3.1 3.2 3.3 3.4 3.5
    (3,0)    3 3.1 3.2 3.3 3.4 3.5
    (4,0)    3 3.1 3.2 3.3 3.4 3.5
    (5,0)    3 3.1 3.2 3.3 3.4 3.5


*/

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"
#include "adios_types.h"


int main (int argc, char ** argv) 
{
	char        filename[] = "joined-array.bp";
	int         rank, size, i,j,n;
	int         ROWS; // different number of rows per process
	int         COLS = 6;
	double     *table;
	int64_t     varid;
	MPI_Comm    comm = MPI_COMM_WORLD;


	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

    ROWS = rank + 3;
	table = (double *) malloc (ROWS * COLS * sizeof(double));

	adios_init_noxml (comm);

	int64_t       m_adios_group;
	int64_t       m_adios_file;

	adios_declare_group (&m_adios_group, "joined", "", adios_stat_default);
	adios_select_method (m_adios_group, "MPI", "verbose=3", "");

	adios_define_var (m_adios_group, "rows" ,"", adios_integer ,0, 0, 0);
	adios_define_var (m_adios_group, "cols" ,"", adios_integer ,0, 0, 0);
    varid = adios_define_var (m_adios_group, "table" ,"", adios_double
				,"rows,cols", "JoinedDim,cols", "0,0");

	adios_open (&m_adios_file, "joined", filename, "w", comm);

    n = 0;
    for (i = 0; i < ROWS; i++)
    {
        for (j = 0; j < COLS; j++)
        {
            table[n] = (double)rank + (double)j/10.0;
            n++;
        }
    }

    adios_write(m_adios_file, "rows", &ROWS);
    adios_write(m_adios_file, "cols", &COLS);
    adios_write(m_adios_file, "table", table);
	adios_close (m_adios_file);

	MPI_Barrier (comm);

	adios_finalize (rank);

	free (table);
	MPI_Finalize ();
	return 0;
}
