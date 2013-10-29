#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "adios.h"
#include "adios_read.h"

#define ADIOS_USE_READ_API_V1

#define BP_FILENAME "transforms.bp"
#define BP_GROUP "transforms"

#define RAW_VAR "rawdouble"
#define XFORM_VAR "transformdouble"
#define TIMESTEP "iter"

#define N 1048576

const MPI_Comm comm = MPI_COMM_WORLD;

static void write_test_file(double *arr) {
    int64_t fd;
    uint64_t total_size;

    // Write the data to file
    adios_open(&fd, BP_GROUP, BP_FILENAME, "w", &comm);
    adios_group_size(fd, 2 * N * sizeof(double), &total_size);
    adios_write(fd, RAW_VAR, arr);
    adios_write(fd, XFORM_VAR, arr);
    adios_close(fd);
}

static void check_raw_var(ADIOS_GROUP *ag, double *arr) {
    int64_t readlen;
    double *readarr = malloc(N * sizeof(double));
    memset(readarr, 0, N*sizeof(double));

    // Check the raw data
    ADIOS_VARINFO *vi_raw = adios_inq_var(ag, RAW_VAR);
    assert(vi_raw);
    assert(vi_raw->type == adios_double);
    assert(vi_raw->ndim == 1); // Apparently the "dummy" dimension has gone away for timed arrays that are written only once
    assert(vi_raw->dims[0] == N);

    readlen = adios_read_var_byid(ag, vi_raw->varid, (uint64_t[]){0}, vi_raw->dims, readarr);

    assert(readlen == N * sizeof(double));
    assert(memcmp(readarr, arr, N * sizeof(double)) == 0);

    // Cleanup
    adios_free_varinfo(vi_raw);
    free(readarr);
}

static void check_xform_var(ADIOS_GROUP *ag, double *arr) {
    int64_t readlen;
    double *readarr = malloc(N * sizeof(double));
    memset(readarr, 0, N*sizeof(double));

    // Check the transformed data
    ADIOS_VARINFO *vi_xform = adios_inq_var(ag, XFORM_VAR);
    assert(vi_xform);
    assert(vi_xform->type == adios_double);
    assert(vi_xform->ndim == 1); // Apparently the "dummy" dimension has gone away for timed arrays that are written only once
    assert(vi_xform->dims[0] == N);

    readlen = adios_read_var_byid(ag, vi_xform->varid, (uint64_t[]){0}, vi_xform->dims, readarr);

    assert(readlen == N * sizeof(double));
    assert(memcmp(readarr, arr, N * sizeof(double)) == 0);

    // Cleanup
    adios_free_varinfo(vi_xform);
    free(readarr);
}

static void read_test_file(double *arr) {
    adios_set_read_method(ADIOS_READ_METHOD_BP);
    ADIOS_FILE *af = adios_fopen(BP_FILENAME, comm);
    ADIOS_GROUP *ag = adios_gopen(af, BP_GROUP);

    check_raw_var(ag, arr);
    check_xform_var(ag, arr);

    // Cleanup
    adios_gclose(ag);
    adios_fclose(af);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    adios_init("transforms.xml", comm);

    double *arr = malloc(N * sizeof(double));
    memset(arr, 123, N * sizeof(double));

    write_test_file(arr);
    read_test_file(arr);

    adios_finalize(0);
    MPI_Finalize();
}
