#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "public/adios.h"
#include "public/adios_read.h"

#define ADIOS_USE_READ_API_V1

#define BP_FILENAME "transforms.bp"
#define BP_GROUP "transforms"

#define RAW_VAR "rawdouble"
#define XFORM_VAR "transformdouble"
#define TIMESTEP "iter"

#define N 100

int main(int argc, char **argv) {
    int64_t fd;
    uint64_t total_size;

    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    adios_init("transforms.xml");

    double *arr = malloc(N * sizeof(double));
    double *readarr = malloc(N * sizeof(double));

    // Write the data to file
    adios_open(&fd, BP_GROUP, BP_FILENAME, "w", &comm);
    adios_group_size(fd, 6 * N * sizeof(double), &total_size);
    adios_write(fd, RAW_VAR, arr);
    adios_write(fd, XFORM_VAR, arr);
    adios_close(fd);

    adios_finalize(0);

    // Now verify the data
    adios_set_read_method(ADIOS_READ_METHOD_BP);
    ADIOS_FILE *af = adios_fopen(BP_FILENAME, comm);
    ADIOS_GROUP *ag = adios_gopen(af, BP_GROUP);

    uint64_t zeros[] = {0, 0};
    int64_t readlen;

    // Check the raw data
    ADIOS_VARINFO *vi_raw = adios_inq_var(ag, RAW_VAR);
    assert(vi_raw);
    assert(vi_raw->type == adios_double);
    assert(vi_raw->ndim == 2); // One "dummy" dimension that always == 1, and the real dimension
    assert(vi_raw->dims[0] == 1);
    assert(vi_raw->dims[1] == N);

    readlen = adios_read_var_byid(ag, vi_raw->varid, &zeros, vi_raw->dims, readarr);

    assert(readlen == N * sizeof(double));
    assert(memcmp(readarr, arr, N * sizeof(double)) == 0);

    adios_free_varinfo(vi_raw);
    // Done checking the raw data

    // Check the transformed data
    ADIOS_VARINFO *vi_xform = adios_inq_var(ag, XFORM_VAR);
    assert(vi_xform);
    assert(vi_xform->type == adios_double);
    assert(vi_xform->ndim == 2); // One "dummy" dimension that always == 1, and the real dimension
    assert(vi_xform->dims[0] == 1);
    assert(vi_xform->dims[1] == N);

    readlen = adios_read_var_byid(ag, vi_xform->varid, &zeros, vi_xform->dims, readarr);

    assert(readlen == N * sizeof(double));
    assert(memcmp(readarr, arr, N * sizeof(double)) == 0);

    adios_free_varinfo(vi_xform);
    // Done checking the transformed data

    // Cleanup
    adios_gclose(ag);
    adios_fclose(af);

    adios_finalize(0);
    MPI_Finalize();
}
