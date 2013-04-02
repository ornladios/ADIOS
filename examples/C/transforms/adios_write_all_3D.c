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
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#define MAX_DIMS 5

struct dimensions {
    uint8_t ndims;
    uint32_t dims [MAX_DIMS];
    uint8_t element_size;
};

typedef struct dimensions dim_t;

// Given the input file, you want to divide the data into different PG sizes,
// with different transforms

void adios_write_pg (char transform [], char input_dir [], uint8_t nvars, char **vars,
                       dim_t data_dim, dim_t pg_dim)
{
    int         rank, size;
    int         i = 0;
    int         t = 0;
    uint32_t    pg_var_size = pg_dim.element_size;
    uint32_t    data_var_size = data_dim.element_size;

    uint32_t    ntimesteps = 1;
    char        varfile [nvars][256];
    FILE        *fp [nvars];

    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    char input_xml [256];
    char output_bp_file [256];

    MPI_Comm    comm = MPI_COMM_WORLD;
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    printf ("pg_var_size = %u, data_var_size = %u\n", pg_var_size, data_var_size);
    // Get the var size and the pg size
    for (i = 0; i < pg_dim.ndims; i ++) {
        pg_var_size *= pg_dim.dims [i];
    }

    for (i = 0; i < data_dim.ndims; i ++) {
        data_var_size *= data_dim.dims [i];
    }
    printf ("pg_var_size = %u, data_var_size = %u\n", pg_var_size, data_var_size);

    // Read the XML file specific to this transform and PG size
    sprintf (input_xml, "%s/%s.xml", input_dir, transform);
    adios_init (input_xml);

    // Name the output bp file based on the name of the transform
    sprintf (output_bp_file, "%s/%s_%d.bp", input_dir, transform, pg_var_size);

    ntimesteps = data_var_size / pg_var_size;

    // Open the input raw data file for each variable
    for (i = 0; i < nvars; i ++) {
        sprintf (varfile [i], "%s/%s.bin", input_dir, vars [i]);
        fp [i] = fopen (varfile [i], "rb");

        assert (fp [i] != 0);
    }

    // Calculate the length, and allocate memory
    uint32_t NX = pg_var_size / pg_dim.element_size;
    char *pg_var_data = (char *) malloc (pg_var_size);

    printf ("ntimesteps = %d, NX = %u\n", ntimesteps, NX);

    for (t = 0; t < ntimesteps; t ++) {

        adios_groupsize = 4 \
                        + 4 \
                        + 4 \
                        + 8 * (1) * (NX) \
                        + 8 * (1) * (NX) \
                        + 8 * (1) * (NX) \
                        + 8 * (1) * (NX) ;

        if (t == 0) {
            adios_open (&adios_handle, "S3D", output_bp_file, "w", &comm);
        } else {
            adios_open (&adios_handle, "S3D", output_bp_file, "a", &comm);
        }

        adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);

        adios_write (adios_handle, "NX", &NX);
        adios_write (adios_handle, "size", &size);
        adios_write (adios_handle, "rank", &rank);

        // fread and adios_write for each variable
        for (i = 0; i < nvars; i ++) {
            fread (pg_var_data, sizeof (char), pg_var_size, fp [i]);
            adios_write (adios_handle, vars [i], pg_var_data);
        }

        adios_close (adios_handle);
    }

    free (pg_var_data);
    for (i = 0; i < nvars; i ++) {
        fclose (fp [i]);
    }

    adios_finalize (rank);

    return ;
}

int main (int argc, char ** argv)
{
    MPI_Init (&argc, &argv);

    dim_t data_dim;
    dim_t pg_dim;

    data_dim.ndims = 4;
    data_dim.dims [0] = 32;
    data_dim.dims [1] = 256;
    data_dim.dims [2] = 256;
    data_dim.dims [3] = 128;
    data_dim.element_size = 8;

    pg_dim.ndims = 3;
    pg_dim.dims [0] = 128;
    pg_dim.dims [1] = 128;
    pg_dim.dims [2] = 128;
    pg_dim.element_size = 8;

    char *vars [4] = {"temp", "uvel", "vvel", "wvel"};
    if (argc >= 2) {
        adios_write_pg (argv [2], argv [1], 4, vars, data_dim, pg_dim);
    } else {
        printf ("Usage: %s <base directory> <transform> <pg size (bytes)>\n", argv [0]);
    }

    MPI_Finalize ();
    return 0;
}

