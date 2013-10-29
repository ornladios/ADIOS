#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"

#define NDIMS	3
// #define DIM_GLOBAL	256
// #define DIM_LOCAL		128

int DIM_GLOBAL = 256;
int DIM_LOCAL = 128;

#define NAME_LEN		256
#define GROUP_NAME	"genarray3D"

void usage()
{
    printf("Usage: mpirun -np 8 ./genarray3D <dim_global> <dim_local> <IO_Compression>\n");
    printf("Example: mpirun -np 8 ./genarray3D 32 16 mpi_zlib\n");
}

double dclock ()
{
    struct timeval tv;
    gettimeofday(&tv,0);
    return (double) tv.tv_sec + (double) tv.tv_usec * 1e-6;
}

int main (int argc, char ** argv)
{
    int i = 0;

    if(argc < 4)
    {
        printf("wrong args\n");
        usage();
        return -1;
    }

    DIM_GLOBAL = atoi (argv[1]);
    DIM_LOCAL = atoi (argv[2]);
    char* option = argv[3];

    char bp_file_name[NAME_LEN] = {0};
    char xml_file_name[NAME_LEN] = {0};

    snprintf(bp_file_name, NAME_LEN-1, "output/%s.bp", option);
    snprintf(xml_file_name, NAME_LEN-1, "conf/%s.xml", option);

    // MPI related intialization
    int rank, nproc;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &nproc);

    double t1 = 0.0;
    double t2 = 0.0;
    double t3 = 0.0;
    double t4 = 0.0;

    // variable dimensions
    int gndx = DIM_GLOBAL;
    int gndy = DIM_GLOBAL;
    int gndz = DIM_GLOBAL;

    int ndx = DIM_LOCAL;
    int ndy = DIM_LOCAL;
    int ndz = DIM_LOCAL;

    int npx = gndx / ndx;
    int npy = gndy / ndy;
    int npz = gndz / ndz;

    if(nproc != npx * npy * npz)
    {
        printf("process num error! nproc != npx * npy * npz\n");
        MPI_Finalize();
        return -1;
    }

    int posx = rank / (npx * npy);
    int posy = rank % (npx * npy) / npy;
    int posz = rank % (npx * npy) % npy;

    // posx = mod(rank, npx)     // 1st dim easy: 0, npx, 2npx... are in the same X position
    // posy = mod(rank/npx, npy) // 2nd dim: (0, npx-1) have the same dim (so divide with npx first)
    // posz = rank/(npx*npy)     // 3rd dim: npx*npy processes belong into one dim
    int offx = posx * ndx;
    int offy = posy * ndy;
    int offz = posz * ndz;

    int timesteps = 0;

    srand(0); // all procs generate the same random datasets

    double* double_xyz = (double*) malloc (sizeof(double) * ndx * ndy * ndz);
    for(i = 0; i < ndx * ndy * ndz; i++)
    {
        double_xyz[i] = (double) rand () / RAND_MAX;
    }

    int adios_err;
    uint64_t adios_groupsize, adios_totalsize;
    int64_t adios_handle;

    if(rank == 0)
        t3 = dclock();

    MPI_Barrier(comm);

    t1 = dclock();

    adios_init (xml_file_name, comm);
    adios_open (&adios_handle, GROUP_NAME, bp_file_name, "w", comm);

    //////////////////////////////////////////////////////////////////////////////////////
    adios_groupsize = 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 4
                    + 8 * (ndx) * (ndy) * (ndz)
                    + 8 * (ndx) * (ndy) * (ndz);
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write (adios_handle, "gndx", &gndx);
    adios_write (adios_handle, "gndy", &gndy);
    adios_write (adios_handle, "gndz", &gndz);
    adios_write (adios_handle, "nproc", &nproc);
    adios_write (adios_handle, "npx", &npx);
    adios_write (adios_handle, "npy", &npy);
    adios_write (adios_handle, "npz", &npz);
    adios_write (adios_handle, "offx", &offx);
    adios_write (adios_handle, "offy", &offy);
    adios_write (adios_handle, "offz", &offz);
    adios_write (adios_handle, "ndx", &ndx);
    adios_write (adios_handle, "ndy", &ndy);
    adios_write (adios_handle, "ndz", &ndz);
    adios_write (adios_handle, "temperature", double_xyz);
    adios_write (adios_handle, "preasure", double_xyz);

    //////////////////////////////////////////////////////////////////////////////////////

    adios_close (adios_handle);

    /*
    t2 = dclock();

    double tt = t2 - t1;

    MPI_Barrier (comm);

    if(rank == 0)
    {
        t4 = dclock();
    }
    */

    adios_finalize (rank);

    /*
    double* all_tt = (double*) malloc (sizeof(double) * nproc);

    // calling MPI_Gather
    int rtn = MPI_Gather (&tt, 1, MPI_DOUBLE, all_tt, 1, MPI_DOUBLE, 0, comm);
    MPI_Barrier (comm);
    if(rank == 0)
    {
        int k = 0;
        double sum = 0.0;
        for(k = 0; k < nproc; k++)
        {
            // printf("proc %d time %f\n", k, all_tt[k]);
            sum += all_tt[k];
        }

        printf("%s average_write_time %f\n", xml_file_name, sum / nproc);
        printf("%s total_write_time %f\n", xml_file_name, t4 - t3);
    }

    if(all_tt)
    {
        free(all_tt);
    }
    */

    MPI_Finalize ();

    if(double_xyz)
    {
        free(double_xyz);
    }


    return 0;
}
