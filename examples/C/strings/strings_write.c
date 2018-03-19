/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: Write strings as a table
 *
 * How to run: mpirun -np <N> strings
 * Output: strings.bp
 * ADIOS config file: None
 *
*/

/* This example will write out many strings into an adios array.
   The data is expected to be stored on every process in a 
   1D array of M char* pointers. The pointers themselves contain the strings.
   The local dimension of the adios array should be declared as M.

   The global dimension and offset could be calculated accordingly but
   since there are different numbers of strings per process, we use the
   joined array trick and let the reader present the per-process
   string tables as one big global table.
*/

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include "mpi.h"
#include "public/adios.h"
#include "public/adios_types.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

char * texts[] = {
    "first_in_rank_", "second_in_rank_", "third_in_rank_", 
    "fourth_in_rank_", "fifth_in_rank_", "sixth_in_rank_"
};


char ** create_table (int n, int rank)
{
    char** table = (char **) malloc (n * sizeof(char*));
    int i;
    printf(" Table on rank %d: {", rank);
    for (i=0; i<n; ++i)
    {
       int idx = i%sizeof(texts);
       table[i] = (char*) malloc ((strlen(texts[idx]) + 5) * sizeof(char));
       sprintf (table[i], "%s%4.4d", texts[idx], rank);
       printf(" \"%s\"", table[i]);
    }
    printf(" }\n");
    return table;
}

void free_table(int n, char ** table)
{
    int i;
    for (i=0; i<n; ++i)
        free(table[i]);
    free(table);
}

int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, block;
    MPI_Comm    comm = MPI_COMM_WORLD;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    /* Create the strings */
    int nstrings = (rank % 3) + 2;
    char** table = create_table(nstrings, rank);

    strcpy (filename, "strings.bp");

    /* Define ADIOS output */
    adios_init_noxml (comm);
    int64_t m_adios_group;
    int64_t m_adios_file;

    char sizeStr[32];
    snprintf (sizeStr, sizeof(sizeStr), "%d", size);


    adios_declare_group (&m_adios_group, "strings", "", adios_stat_default);
    adios_select_method (m_adios_group, "MPI", "verbose=3", "");

    adios_define_var (m_adios_group, "N","", adios_integer,
                      "1", "JoinedDim", "0");

    adios_define_var (m_adios_group, "Table", "", adios_string,
                      sizeStr, "JoinedDim", "0");

    adios_open (&m_adios_file, "strings", filename, "w", comm);
    adios_write(m_adios_file, "N", (void *) &nstrings);
    adios_write(m_adios_file, "Table", table);

    adios_close (m_adios_file);
    MPI_Barrier (comm);
    adios_finalize (rank);
    free (table);
    MPI_Finalize ();
    return 0;
}
