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
#include "hdf5.h"

#define RANK            2
#define DATASETNAME     "xy" 

hsize_t dimsf[RANK];              /* dataset dimensions */
hsize_t count[RANK];              /* hyperslab selection parameters */
hsize_t offset[RANK];

hid_t   plist_id1;                /* property list identifier */
hid_t   plist_id2;                /* property list identifier */

int output_init(MPI_Comm comm, int bufsizeMB)
{
    /* 
     * Set up file access property list with parallel I/O access
     */
    plist_id1 = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id1, comm, MPI_INFO_NULL);

    /*
     * Create property list for collective dataset write.
     */
    plist_id2 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id2, H5FD_MPIO_COLLECTIVE);

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
    hid_t   filespace, memspace;      /* file and memory dataspace identifiers */

    herr_t      status;
    double t1, t2;

    /*
     * Create a new file collectively and release property list identifier.
     */
    t1 = MPI_Wtime();
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id1);
    t2 = MPI_Wtime();
    Tio_open[step] = t2-t1;
   
    /*
     * Create the dataspace for the dataset.
     */
    filespace = H5Screate_simple(RANK, dimsf, NULL); 
    memspace = H5Screate_simple(RANK, count, NULL);

    /*
     * Create the dataset with default properties and close filespace.
     */
    dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_DOUBLE, filespace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    /*
     * Select hyperslab in the file.
     */
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    t1 = MPI_Wtime();
    Tio_group[step] = t1-t2;
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace,
                      plist_id2, data);
    t2 = MPI_Wtime();
    Tio_write[step] = t2-t1;

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
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

    H5Pclose(plist_id1);
    H5Pclose(plist_id2);

    return 0;
}

