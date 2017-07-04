/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: read global arrays from a BP file
 *
 * This code is using the generic read API, which can read in
 * arbitrary slices of an array and thus we can read in an array
 * on arbitrary number of processes (provided our code is smart 
 * enough to do the domain decomposition).
 *
 * Run this example after adios_global, which generates 
 * adios_global.bp. Run this example on equal or less 
 * number of processes since we decompose only on one 
 * dimension of the global array here. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mpi.h"
#include "adios_read.h"
#include "adios.h"
#include "zfp.h"

int
compressor (double * array, int nx, double tolerance,
          double ** array_compressed)
{
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream * stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_double;
    field = zfp_field_1d (array, type, nx);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open (NULL);

    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision, type); */
    zfp_stream_set_accuracy (zfp, tolerance, type);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size (zfp, field);
    buffer = malloc (bufsize);
    printf ("buffer size = %lu\n", bufsize);
    assert (buffer);

    /* associate bit stream with allocated buffer */
    stream = stream_open (buffer, bufsize);
    zfp_stream_set_bit_stream (zfp, stream);
    zfp_stream_rewind (zfp);

    /* compress array and output compressed stream */
    zfpsize = zfp_compress (zfp, field);
    assert (zfpsize);

    /* clean up */
    zfp_field_free (field);
    zfp_stream_close (zfp);
    stream_close (stream);

    * array_compressed = (double *) buffer;

    return zfpsize;
}

int
decompressor (double * array, int nx, double tolerance,
              double * array_compressed, size_t array_size_compressed)
{
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    bitstream * stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_double;
    field = zfp_field_1d (array, type, nx);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open (NULL);
    zfp_stream_set_accuracy (zfp, tolerance, type);

    /* associate bit stream with allocated buffer */
    stream = stream_open (array_compressed, array_size_compressed);
    zfp_stream_set_bit_stream (zfp, stream);
    zfp_stream_rewind (zfp);

    assert (zfp_decompress(zfp, field));

    /* clean up */
    zfp_field_free (field);
    zfp_stream_close (zfp);
    stream_close (stream);

    return 0;
}

int main (int argc, char ** argv) 
{
    int         rank, size, i, j;
    MPI_Comm    comm = MPI_COMM_WORLD;
    enum ADIOS_READ_METHOD method = ADIOS_READ_METHOD_BP;
    ADIOS_SELECTION * sel;
    void * mesh = NULL;
    double * data = NULL, * R = NULL, * Z = NULL;
    uint64_t start[5], count[5], nelements = 1;

    MPI_Init (&argc, &argv);

    if (argc != 2) printf ("Wrong arguments\n");
    double tolerance = atof(argv[1]);

    printf ("tolerance = %f\n", tolerance);

    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    if (size != 32)
    {
        printf ("This needs to run with 32 processors.");
        //return -1;
    }

    int step;
    double iotime = 0.0;

    for (step = 0; step < 1; step++)
    {
        char fname[64];

        //sprintf (fname, "xgc.3d.080%02d.bp", step);
        //sprintf (fname, "larger_data/xgc.3d.00480.bp");
        sprintf (fname, "xgc-f0/xgc.dataspaces.f0analysis.00005.bp");
        adios_read_init_method (method, comm, "verbose=3");

        ADIOS_FILE * f = adios_read_open_file (fname, method, comm);
        if (f == NULL)
        {
            printf ("%s\n", adios_errmsg());
            return -1;
        }

        //ADIOS_VARINFO * v = adios_inq_var (f, "dpot");
        ADIOS_VARINFO * v = adios_inq_var (f, "f0_f");

        start[0] = 0;
        start[1] = 0;
        start[2] = rank;
        start[3] = 0;
        start[4] = 0;

        count[0] = v->dims[0];
        count[1] = v->dims[1];
        count[2] = 1;
        count[3] = v->dims[3];
        count[4] = v->dims[4];

        nelements = count[0] * count[1] * count[2] * count[3] * count[4];
        data = malloc (nelements * sizeof (double));
        if (data == NULL)
        {
            fprintf (stderr, "malloc failed.\n");
            return -1;
        }

        /* Read a subset of the temperature array */
        sel = adios_selection_boundingbox (v->ndim, start, count);
        adios_schedule_read (f, sel, "f0_f", 0, 1, data);
        adios_perform_reads (f, 1);

        adios_read_close (f);
        MPI_Barrier (comm);
        adios_read_finalize_method (method);

        double * data_compressed = 0;
        // Compress the data using ZFP
        int compressed_size = compressor (data, nelements, tolerance, &data_compressed);
        printf ("compression ratio  = %f\n", ((double) nelements*8) / compressed_size);

        decompressor (data, nelements, tolerance,
                      data_compressed, compressed_size);

        free (data_compressed);
#if 0
        // read mesh
        adios_read_init_method (method, comm, "verbose=3");

        ADIOS_FILE * fmesh = adios_read_open_file ("larger_data/xgc.mesh.bp", method, comm);
        if (fmesh == NULL)
        {
            printf ("%s\n", adios_errmsg());
            return -1;
        }

        ADIOS_VARINFO * conn = adios_inq_var (fmesh, "/cell_set[0]/node_connect_list");

        start[0] = 0;
        count[0] = conn->dims[0];

        start[1] = 0;
        count[1] = conn->dims[1];


        mesh = malloc (count[0] * count[1] * sizeof (int));
        if (mesh == NULL)
        {
            fprintf (stderr, "malloc failed.\n");
            return -1;
        }

        /* Read a subset of the temperature array */
        sel = adios_selection_boundingbox (conn->ndim, start, count);
        adios_schedule_read (fmesh, sel, "/cell_set[0]/node_connect_list", 0, 1, mesh);
        start[0] = 0;
        start[1] = 0;

        count[0] = v->dims[0];
        count[1] = 1;

        sel = adios_selection_boundingbox (v->ndim, start, count);
        adios_schedule_read (fmesh, sel, "/coordinates/values", 0, 1, R);


        start[0] = 0;
        start[1] = 1;

        count[0] = v->dims[0];
        count[1] = 1;

        sel = adios_selection_boundingbox (v->ndim, start, count);
        adios_schedule_read (fmesh, sel, "/coordinates/values", 0, 1, Z);

        adios_perform_reads (fmesh, 1);

        adios_read_close (fmesh);
        MPI_Barrier (comm);
        adios_read_finalize_method (method);
#endif
#if 1
        int NX = v->dims[3], GX = v->dims[3], OX = 0;
//        int MY = conn->dims[0], MX = conn->dims[1];

        uint64_t    adios_groupsize, adios_totalsize;
        int64_t     adios_handle;

        MPI_Barrier (comm);
        double start_io_time = MPI_Wtime ();

        adios_init_noxml (comm);
        adios_set_max_buffer_size (100);

        int64_t       m_adios_group;
        int64_t       m_adios_file;

        adios_declare_group (&m_adios_group, "field", "iter", adios_flag_yes);
        adios_select_method (m_adios_group, "MPI", "", "");


        adios_define_var (m_adios_group, "NX"
                        ,"", adios_integer
                        ,0, 0, 0);

        adios_define_var (m_adios_group, "OX"
                        ,"", adios_integer
                        ,0, 0, 0);

        adios_define_var (m_adios_group, "GX"
                        ,"", adios_integer
                        ,0, 0, 0);



        int64_t varid;
        char local_dims[128], global_dims[128], offsets[128];
        sprintf (local_dims, "%lu,%lu,%d,%lu,%lu", v->dims[0],
                                               v->dims[1],
                                               1,
                                               v->dims[3],
                                               v->dims[4]);
        sprintf (global_dims, "%lu,%lu,%lu,%lu,%lu", v->dims[0],
                                               v->dims[1],
                                               v->dims[2], 
                                               v->dims[3],
                                               v->dims[4]);
        sprintf (offsets, "%d,%d,%d,%d,%d", 0,
                                            0,
                                            rank,
                                            0,
                                            0);

        varid = adios_define_var (m_adios_group, "f0_f"
                        ,"", adios_double
                        ,local_dims, global_dims, offsets);
        adios_set_transform (varid, "none");

        sprintf (fname, "f0_f.zfp.%d.bp", step);

        adios_open (&adios_handle, "field", fname, "w", comm);
        adios_groupsize = 6 * 4 + 8 * v->dims[0] * v->dims[1] * v->dims[3] * v->dims[4];
        adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);

        adios_write (adios_handle, "NX", &NX);
        adios_write (adios_handle, "OX", &OX);
        adios_write (adios_handle, "GX", &GX);
        adios_write (adios_handle, "f0_f", data);
 
        adios_close (adios_handle);

        free (data);

        adios_finalize (rank);

        MPI_Barrier (MPI_COMM_WORLD);
        double end_io_time = MPI_Wtime ();
    
        iotime += end_io_time - start_io_time;
#endif
    }

    //if (rank == 0) printf ("io time = %f\n", iotime);
    MPI_Finalize ();
    return 0;
}
