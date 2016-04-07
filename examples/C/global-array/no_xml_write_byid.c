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
#include "public/adios.h"
#include "public/adios_types.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

int main (int argc, char ** argv) 
{
	char        filename [256];
	int         rank, size, i, j;
	int         NX = 100, gb, offset;  //local/global/offset
	double      t[NX];
        int         nblocks = 3;
	MPI_Comm    comm = MPI_COMM_WORLD;
        char g_str[100], o_str[100], l_str[100];
        // attributes (from C variables)
        int someints[5] = {5,4,3,2,1};
        double somedoubles[5] = {5.55555, 4.4444, 3.333, 2.22, 1.1};

	/* ADIOS variables declarations for matching gwrite_temperature.ch */
	uint64_t    adios_groupsize, adios_totalsize;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (comm, &rank);
	MPI_Comm_size (comm, &size);

        gb = nblocks * NX * size;
        sprintf (g_str, "%d", gb);
        sprintf (l_str, "%d", NX);

	strcpy (filename, "no_xml_write_byid.bp");

	adios_init_noxml (comm);
        adios_set_max_buffer_size (10);

        int64_t       m_adios_group;
        int64_t       m_adios_file;
        int64_t       var_ids[nblocks];

        adios_declare_group (&m_adios_group, "restart", "iter", adios_flag_yes);
        adios_select_method (m_adios_group, "MPI", "", "");

        for (i = 0; i < nblocks; i++) 
        {
            offset = rank * nblocks * NX + i * NX;
            sprintf (o_str, "%d", offset);
            var_ids[i] = adios_define_var (m_adios_group, "temperature"
                                          ,"", adios_double
                                          ,l_str, g_str, o_str
                                          );
            adios_set_transform (var_ids[i], "none");
        }

        // add some attributes
        adios_define_attribute_byvalue (m_adios_group, 
                "single_string","", adios_string,  1, "A single string attribute");
        char *strings[] = {"X","Yy","ZzZ"};
        adios_define_attribute_byvalue (m_adios_group, 
                "three_strings","", adios_string_array,  3, strings);
        adios_define_attribute_byvalue (m_adios_group, 
                "single_int",   "", adios_integer, 1, &someints);
        adios_define_attribute_byvalue (m_adios_group, 
                "single_double","", adios_double,  1, &somedoubles);
        adios_define_attribute_byvalue (m_adios_group, 
                "five_ints",    "", adios_integer, 5, &someints);
        adios_define_attribute_byvalue (m_adios_group, 
                "five_double",  "", adios_double,  5, &somedoubles);



        adios_open (&m_adios_file, "restart", filename, "w", comm);

        adios_groupsize = nblocks * (4 + 4 + 4 + NX * 8);

        adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);
/* now we will write the data for each sub block */
        for (i = 0; i < nblocks; i++)
        {
           offset = rank * nblocks * NX + i * NX;
           for (j = 0; j < NX; j++)
               t[j] = offset + j;

           adios_write_byid(m_adios_file, var_ids[i], t);
        }

        adios_close (m_adios_file);

        MPI_Barrier (comm);

	adios_finalize (rank);

	MPI_Finalize ();
	return 0;
}
