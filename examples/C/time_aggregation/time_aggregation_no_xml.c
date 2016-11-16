/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"

const int  NX = 10;
const int  NY = 5;
const char diagfilename[] = "diag.bp";
const char diag2filename[] = "diag2.bp";
const char ckptfilename[] = "ckpt.bp";

const MPI_Comm comm = MPI_COMM_WORLD;
int rank, size;

void write_diag (int step, double * p)
{
    int64_t     adios_handle;
    if (rank==0) printf("    write diagnostics\n");

    if (step==1)
        adios_open (&adios_handle, "diagnostics", diagfilename, "w", comm);
    else
        adios_open (&adios_handle, "diagnostics", diagfilename, "a", comm);

    adios_write (adios_handle, "NY", &NY);
    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "pressure", p);

    adios_close (adios_handle);
}

void write_diag2 (int step, double * t)
{
    int64_t     adios_handle;
    if (rank==0) printf("    write diag2\n");

    if (step==1)
        adios_open (&adios_handle, "diag2", diag2filename, "w", comm);
    else
        adios_open (&adios_handle, "diag2", diag2filename, "a", comm);

    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "t0", t);

    adios_close (adios_handle);
}

void write_checkpoint (int step, double * p, double *t)
{
    int64_t     adios_handle;
    if (rank==0) printf("    Checkpointing at step %d\n", step);
    adios_open (&adios_handle, "checkpoint", ckptfilename, "w", comm);
    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "NY", &NY);
    adios_write (adios_handle, "size", &size);
    adios_write (adios_handle, "rank", &rank);
    adios_write (adios_handle, "step", &step);
    adios_write (adios_handle, "temperature", t);
    adios_write (adios_handle, "pressure", p);
    adios_close (adios_handle);
}

void define_groups ()
{
    int64_t  g_diag, g_diag2, g_ckpt;

    // Group diagnosis
    adios_declare_group (&g_diag, "diagnostics", "", adios_stat_default);
    adios_define_var (g_diag, "NY",  "", adios_integer, 0, 0, 0);
    adios_define_var (g_diag, "size","", adios_integer, 0, 0, 0);
    adios_define_var (g_diag, "rank","", adios_integer, 0, 0, 0);
    int64_t var_p = adios_define_var (g_diag, "pressure", "", adios_double,
                                       "1,NY", "size,NY", "rank,0");
    adios_set_transform (var_p, "none");


    // Group diag2
    adios_declare_group (&g_diag2, "diag2", "", adios_stat_default);
    adios_define_var (g_diag2, "size","", adios_integer, 0, 0, 0);
    adios_define_var (g_diag2, "rank","", adios_integer, 0, 0, 0);
    int64_t var_t0 = adios_define_var (g_diag2, "t0", "", adios_double,
                                       "1,1", "size,1", "rank,0");
    adios_set_transform (var_t0, "none");


    // Group checkpoint
    adios_declare_group (&g_ckpt, "checkpoint", "", adios_stat_default);
    adios_define_var (g_ckpt, "NX",  "", adios_integer, 0, 0, 0);
    adios_define_var (g_ckpt, "NY",  "", adios_integer, 0, 0, 0);
    adios_define_var (g_ckpt, "size","", adios_integer, 0, 0, 0);
    adios_define_var (g_ckpt, "rank","", adios_integer, 0, 0, 0);
    adios_define_var (g_ckpt, "step","", adios_integer, 0, 0, 0);
    adios_define_var (g_ckpt, "temperature", "", adios_double,
                      "1,NX", "size,NX", "rank,0");
    adios_define_var (g_ckpt, "pressure", "", adios_double,
                      "1,NY", "size,NY", "rank,0");


    adios_select_method (g_diag,  "MPI", "verbose=3", "");
    adios_select_method (g_diag2, "MPI", "verbose=3", "");
    adios_select_method (g_ckpt,  "MPI", "verbose=3", "");
    //adios_select_method (m_ckpt, "MPI_AGGREGATE", "num_ost=2;num_aggregators=2;aggregation_type=2;verbose=3", "");

    adios_set_time_aggregation (g_diag, 12000, g_ckpt);
    adios_set_time_aggregation (g_diag2, 32000, g_ckpt);
}

int main (int argc, char ** argv) 
{

    int         i, it;
    double      t[NX];
    double      p[NY];

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_init_noxml (comm);
    adios_set_max_buffer_size (10);
    define_groups();


    for (it = 1; it <= 100; it++)
    {
        if (rank==0) printf("Timestep %d...\n", it);

        for (i = 0; i < NX; i++)
            t[i] = it*1000.0 + rank*NX + i;

        for (i = 0; i < NY; i++)
            p[i] = it*1000.0 + rank*NY + i;

        write_diag(it, p);
        write_diag2(it, t);

        if ( it%30 == 0) {
            write_checkpoint(it, p, t);
        }

        MPI_Barrier (comm);
        if (rank==0) printf("    step completed\n");
    }

    MPI_Barrier (comm);
    adios_finalize (rank);
    MPI_Finalize ();
    return 0;
}
