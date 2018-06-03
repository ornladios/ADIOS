/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*************************************************************/
/*          Example of reading arrays in ADIOS               */
/*    which were written from the same number of processors  */
/*                                                           */
/*        Similar example is manual/2_adios_read.c           */
/*************************************************************/
#include "adios.h"
#include "mpi.h"
#include "public/adios_read.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "misc.h"
#include "test_common.h"
#include "utils.h"

int main(int argc, char **argv)
{
    int rank;
    int NX, NY;
    int writer_size;
    double *t;
    MPI_Comm comm = MPI_COMM_WORLD;
    diag_t return_val;

    struct err_counts err = {0, 0};
    struct adios_tsprt_opts adios_opts;
    struct test_info test_result = {TEST_PASSED, "global_range_select"};

    GET_ENTRY_OPTIONS(
        adios_opts,
        "Runs readers.");
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);

    SET_ERROR_IF_NOT_ZERO(adios_read_init_method(adios_opts.method, comm,
                                                 adios_opts.adios_options),
                          err.adios);
    RET_IF_ERROR(err.adios, rank);

    /* schedule_read of a scalar. */
    int test_scalar = -1;
    ADIOS_FILE *afile = adios_read_open(FILE_NAME, adios_opts.method, comm,
                                        ADIOS_LOCKMODE_NONE, 0.0);

    int ii = 0;
    while (adios_errno != err_end_of_stream) {

        /* get a bounding box - rank 0 for now*/
        ADIOS_VARINFO *nx_info = adios_inq_var(afile, "NX");
        ADIOS_VARINFO *ny_info = adios_inq_var(afile, "NY");
        ADIOS_VARINFO *writer_size_info = adios_inq_var(afile, "size");

        if (nx_info->value) {
            NX = *((int *)nx_info->value);
        } else {
            test_failed(test_result.name, rank);
            return_val = DIAG_ERR;
        }

        if (ny_info->value) {
            NY = *((int *)ny_info->value);
        } else {
            test_failed(test_result.name, rank);
            return_val = DIAG_ERR;
        }

        if (writer_size_info->value) {
            writer_size = *((int *)writer_size_info->value);
        } else {
            test_failed(test_result.name, rank);
            return_val = DIAG_ERR;
        }

        int sel_index = rank % writer_size;
        int k;

        /* Allocate space for the arrays */
        ADIOS_SELECTION *global_range_select;
        ADIOS_SELECTION *scalar_block_select;

        int nelem = writer_size * NX * NY;
        size_t arr_size = sizeof(double) * nelem;
        t = (double *)malloc(arr_size);
        memset(t, 0, arr_size);

        uint64_t *start_array = malloc(sizeof(uint64_t) * (NY + 1));
        uint64_t *count_array = malloc(sizeof(uint64_t) * (NY + 1));
        start_array[0] = 0;
        count_array[0] = writer_size;
        for (k = 1; k < (NY + 1); ++k) {
            start_array[k] = 0;
            count_array[k] = NX;
        }

        global_range_select =
            adios_selection_boundingbox((NY + 1), start_array, count_array);
        scalar_block_select = adios_selection_writeblock(sel_index);

        /* Read the arrays */
        adios_schedule_read(afile, global_range_select, "var_2d_array", 0, 1,
                            t);
        adios_schedule_read(afile, scalar_block_select, "test_scalar", 0, 1,
                            &test_scalar);

        adios_perform_reads(afile, 1);

        adios_release_step(afile);
        adios_advance_step(afile, 0, 30);
        ii++;
    }

just_clean:
    // clean everything
    free(t);
    t = NULL;

close_adios:
    CLOSE_ADIOS_READER(afile, adios_opts.method);

    return 0;
}
