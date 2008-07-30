#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

// mpi
#include "mpi.h"

#include "adios.h"

int main (int argc, char ** argv)
{
    char * type_name = "restart";
    char * filename = "restart.bp";
    long long io_handle;  /* io handle */
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;

    int var_x1 = 101;
    int var_x2 = 102;
    int z_dim_size = 10;
    float z_dim [z_dim_size];

    int r_var_x1;
    int r_var_x2;
    int r_zsize;
    float r_z [z_dim_size];

    int node = 0;

    uint64_t total;

    z_dim [0] = 10.0;
    z_dim [1] = 11.0;
    z_dim [2] = 12.0;
    z_dim [3] = 13.0;
    z_dim [4] = 14.0;
    z_dim [5] = 15.0;
    z_dim [6] = 16.0;
    z_dim [7] = 17.0;
    z_dim [8] = 18.0;
    z_dim [9] = 19.0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    adios_init ("config_c.xml");

    adios_open (&io_handle, type_name, filename, "w");
    adios_group_size (io_handle, 4 + 4 + 4 + 4 + 4 * 10 + 4, &total, &comm);
    adios_write (io_handle, "comm", &comm);
    adios_write (io_handle, "/mype", &var_x1);
    adios_write (io_handle, "/test/mype", &var_x2);
    adios_write (io_handle, "zionsize", &z_dim_size);
    adios_write (io_handle, "zion", z_dim);
    adios_write (io_handle, "zion1", z_dim);
    adios_write (io_handle, "node-attr", &node);
    adios_close (io_handle);

    printf ("rank: %d write completed\n", rank);

    MPI_Barrier (MPI_COMM_WORLD);

    adios_open (&io_handle, type_name, filename, "r");
    adios_write (io_handle, "comm", &comm);
    //adios_write (io_handle, "zionsize", &z_dim_size);
    adios_read (io_handle, "/mype", &r_var_x1);
    adios_read (io_handle, "/test/mype", &r_var_x2);
    adios_read (io_handle, "zionsize", &r_zsize);
    adios_read (io_handle, "zion", &r_z);
    adios_close (io_handle);

    MPI_Barrier (MPI_COMM_WORLD);

    if (   var_x1 != r_var_x1
        || var_x2 != r_var_x2
        || r_zsize != r_zsize
        || r_z [0] != z_dim [0]
        || r_z [1] != z_dim [1]
        || r_z [2] != z_dim [2]
        || r_z [3] != z_dim [3]
        || r_z [4] != z_dim [4]
        || r_z [5] != z_dim [5]
        || r_z [6] != z_dim [6]
        || r_z [7] != z_dim [7]
        || r_z [8] != z_dim [8]
        || r_z [9] != z_dim [9]
       )
    {
        printf ("rank: %d mismatch in reading\n", rank);
    }
    else
    {
        printf ("rank: %d read matches write\n", rank);
    }

#if 0
    var_x = 11;
    adios_open (&io_handle, type_name, filename, "a");
    adios_group_size (io_handle, 4 + 4 + 4 + 4 + 4 * 10 + 4, &total, &comm);
    adios_write (io_handle, "comm", &comm);
    adios_write (io_handle, "mype", &var_x);
    adios_write (io_handle, "zionsize", &z_dim_size);
    adios_write (io_handle, "zion", z_dim);
    adios_close (io_handle);

    MPI_Barrier (MPI_COMM_WORLD);
#endif

    adios_finalize (node);
    MPI_Finalize ();

    return 0;
}
