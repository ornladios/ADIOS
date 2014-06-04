/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS output functions for write.c */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include "output.h"
#include "timing.h"

#include "mpi.h"
#include "hdf5.h"  // USE sequential HDF5 here

#define RANK            2
#define DATASETNAME     "xy" 

hsize_t dimsf[RANK];              /* global dataset dimensions */
hsize_t count[RANK];              /* hyperslab selection parameters */
hsize_t offset[RANK];
int rank;


int output_init(MPI_Comm comm, int bufsizeMB)
{
    MPI_Comm_rank (comm, &rank);
    return 0;
}

int output_define(int nx, int ny, int gnx, int gny, int offsx, int offsy)
{
    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = gnx;
    dimsf[1] = gny;

    count[0] = nx;
    count[1] = ny;

    offset[0] = offsx;
    offset[1] = offsy;

    return 0;
}

int output_dump(char *filename, int step, void *data)
{
    hid_t   file_id, dset_id;         /* file and dataset identifiers */
    hid_t   space;      /* file and memory dataspace identifiers */

    herr_t      status;
    double t1, t2;

    char fname[256];
    snprintf (fname, sizeof(fname), "%s.%d", filename, rank);

    /*
     * Create a new file collectively and release property list identifier.
     */
    t1 = MPI_Wtime();
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    t2 = MPI_Wtime();
    Tio_open[step] = t2-t1;
   
    /*
     * Create the dataspace for the dataset.
     */
    space = H5Screate_simple(RANK, count, NULL); 

    /*
     * Create the dataset with default properties 
     */
    dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_DOUBLE, space,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    t1 = MPI_Wtime();
    Tio_group[step] = t1-t2;
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, data);
    t2 = MPI_Wtime();
    Tio_write[step] = t2-t1;

    H5Dclose(dset_id);
    H5Sclose(space);
    H5Fclose(file_id);

    t1 = MPI_Wtime();
    Tio_close[step] = t1-t2;

    if (status < 0)
        return -1;
    else
        return 0;
}

int output_finalize (int rank)
{
    return 0;
}

