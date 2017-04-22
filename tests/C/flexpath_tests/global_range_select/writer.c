/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"
#include "adios.h"

#include "misc.h"
#include "utils.h"
#include "test_common.h"
/*************************************************************/
/*          Example of writing arrays in ADIOS               */
/*                                                           */
/*        Similar example is manual/2_adios_write.c          */
/*************************************************************/
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i;
    int         NX = 40; 
    int         NY = 1;
    double      t[NX];
    MPI_Comm    comm = MPI_COMM_WORLD;

    int64_t     adios_handle;
    struct adios_tsprt_opts adios_opts;
    int err_count = 0;

    GET_ENTRY_OPTIONS(adios_opts, "Runs writers. It is recommended to run as many writers as readers.");

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    SET_ERROR_IF_NOT_ZERO(adios_init(adios_opts.xml_adios_init_filename, comm), err_count);
    RET_IF_ERROR(err_count, rank);
    
    int test_scalar = rank * 1000;
    int group_num;
    for(group_num = 0; group_num<20; group_num++){
      for (i = 0; i < NX; i++)
        t[i] = rank * NX + i + 100*group_num;
    
      adios_open (&adios_handle, "temperature", FILE_NAME, "w", comm);
    
      test_scalar++;
      adios_write (adios_handle, "NX", &NX);
      adios_write (adios_handle, "NY", &NY);
      adios_write (adios_handle, "test_scalar", &test_scalar);
      adios_write (adios_handle, "size", &size);
      adios_write (adios_handle, "rank", &rank);
      adios_write (adios_handle, "var_2d_array", t);

      fprintf(stderr, "Writer side Rank=%d: test_scalar: %d step: %d, t[0] = %6.2f\n",rank, test_scalar, group_num, t[0]);
		//int j; // CFLAGS may not be C99 (sigh)
		/*for(j=1; j < NX; j++) {
		    fprintf(stderr, ", %6.2f", t[j]);
		}
		fprintf(stderr, "]\n");
               */
    
      adios_close (adios_handle);
      fprintf(stderr, "Rank=%d commited write %d\n", rank, group_num);
      //sleep(11);
    }
    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}

