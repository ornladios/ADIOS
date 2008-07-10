#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// mpi
#include "mpi.h"

#include "adios.h"

int main (int argc, char ** argv)
{
    char * type_name = "restart";
    char * filename = "restart.bp";
    long long io_handle;  /* io handle */

    int var_x = 10;
    int z_dim_size = 10;
    float z_dim [z_dim_size];

    int r_var_x;
    int r_zsize;
    float r_z [z_dim_size];

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

    int node = 0;

    MPI_Init (&argc, &argv);
    adios_init ("config_c.xml"); //, MPI_COMM_WORLD, MPI_COMM_SELF, MPI_INFO_NULL);

    adios_open (&io_handle, type_name, filename, "w");
    adios_write (io_handle, "mype", &var_x);
    adios_write (io_handle, "zionsize", &z_dim_size);
    adios_write (io_handle, "zion", z_dim);
    adios_write (io_handle, "node-attr", &node);
    adios_close (io_handle);

    MPI_Barrier (MPI_COMM_WORLD);

    adios_open (&io_handle, type_name, filename, "r");
    adios_read (io_handle, "mype", &r_var_x);
    adios_read (io_handle, "zionsize", &r_zsize);
    adios_read (io_handle, "zion", &r_z);
    adios_close (io_handle);

    if (   var_x != r_var_x
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
        printf ("mismatch in reading\n");
    }
    else
    {
        printf ("read matches write\n");
    }

    MPI_Barrier (MPI_COMM_WORLD);

    var_x = 11;
    adios_open (&io_handle, type_name, filename, "a");
    adios_write (io_handle, "mype", &var_x);
    adios_write (io_handle, "zionsize", &z_dim_size);
    adios_write (io_handle, "zion", z_dim);
    adios_close (io_handle);

    MPI_Barrier (MPI_COMM_WORLD);

    adios_finalize (node);
    MPI_Finalize ();

    return 0;
}
