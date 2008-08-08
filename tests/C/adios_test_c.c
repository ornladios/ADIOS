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
    long long io_handle;  // io handle
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;

    int var_x1 = 101;
    int var_x2 = 102;

    int z_dim_size1 = 10;
    int z_dim_size2 = 2;
    int z_dim_size3 = 5;

    float * z_dim1;
    float * z_dim2;
    float * z_dim3;

    z_dim1 = malloc (z_dim_size1 * sizeof (float));
    z_dim2 = malloc (z_dim_size1 * z_dim_size2 * sizeof (float));
    z_dim3 = malloc (z_dim_size2 * z_dim_size3 * sizeof (float));
    int r_var_x1;
    int r_var_x2;
    int r_zsize;
    float r_z [z_dim_size1];

    int node = 0;

    uint64_t total;

    z_dim1 [0] = 10.0;
    z_dim1 [1] = 11.0;
    z_dim1 [2] = 12.0;
    z_dim1 [3] = 13.0;
    z_dim1 [4] = 14.0;
    z_dim1 [5] = 15.0;
    z_dim1 [6] = 16.0;
    z_dim1 [7] = 17.0;
    z_dim1 [8] = 18.0;
    z_dim1 [9] = 19.0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    adios_init ("config_c.xml");

#if 1
printf ("XXXXXXXXXXXXXXXX do a write XXXXXXXXXXXXXXXXX\n");
    adios_open (&io_handle, type_name, filename, "w");
    adios_group_size (io_handle,  4 + 4
                                + 4 * z_dim_size1
                                + 4 + z_dim_size1 * z_dim_size2
                                + 4 + 4 * z_dim_size2 * z_dim_size3
                                + 4
                     ,&total, &comm
                     );
    adios_write (io_handle, "/mype", &var_x1);
    adios_write (io_handle, "/test/mype", &var_x2);

    adios_write (io_handle, "zion1", z_dim1);

    adios_write (io_handle, "zionsize2", &z_dim_size1);
    adios_write (io_handle, "zion2", z_dim2);

    adios_write (io_handle, "zionsize3", &z_dim_size3);
    adios_write (io_handle, "zion3", z_dim3);

    adios_write (io_handle, "node-attr", &node);

    adios_close (io_handle);

    printf ("rank: %d write completed\n", rank);

    MPI_Barrier (MPI_COMM_WORLD);
#endif

#if 0
printf ("XXXXXXXXXXXXXXXX do a read XXXXXXXXXXXXXXXXX\n");

    adios_open (&io_handle, type_name, filename, "r");
    adios_group_size (io_handle, 0, &total, &comm);
    adios_read (io_handle, "/mype", &r_var_x1, 4);
    adios_read (io_handle, "/test/mype", &r_var_x2, 4);
    adios_read (io_handle, "zionsize2", &r_zsize, 4);
    adios_read (io_handle, "zion1", r_z, 4 * 10);
    adios_close (io_handle);

    MPI_Barrier (MPI_COMM_WORLD);

    if (   var_x1 != r_var_x1
        || var_x2 != r_var_x2
        || r_zsize != z_dim_size
        || r_z [0] != z_dim1 [0]
        || r_z [1] != z_dim1 [1]
        || r_z [2] != z_dim1 [2]
        || r_z [3] != z_dim1 [3]
        || r_z [4] != z_dim1 [4]
        || r_z [5] != z_dim1 [5]
        || r_z [6] != z_dim1 [6]
        || r_z [7] != z_dim1 [7]
        || r_z [8] != z_dim1 [8]
        || r_z [9] != z_dim1 [9]
       )
    {
        int i;
        printf ("rank: %d mismatch in reading\n", rank);
        printf ("r_var_x1: %d (%d)\n", r_var_x1, var_x1);
        printf ("r_var_x2: %d (%d)\n", r_var_x2, var_x2);
        printf ("r_zsize: %d (%d)\n", r_zsize, z_dim_size);
        for (i = 0; i < 10; i++)
            printf ("z_dim [%d]: %f (%f)\n", i, r_z [i], z_dim [i]);
    }
    else
    {
        printf ("rank: %d read matches write\n", rank);
    }
#endif

#if 0
printf ("XXXXXXXXXXXXXXXX do an append XXXXXXXXXXXXXXXXX\n");
    var_x1 = 11;
    adios_open (&io_handle, type_name, filename, "a");
    adios_group_size (io_handle,  4 + 4
                                + 4 * z_dim_size1
                                + 4 + z_dim_size1 * z_dim_size2
                                + 4
                     ,&total, &comm
                     );
    adios_write (io_handle, "/mype", &var_x1);
    adios_write (io_handle, "/test/mype", &var_x2);

    adios_write (io_handle, "zion1", z_dim1);

    adios_write (io_handle, "zionsize2", &z_dim_size);
    adios_write (io_handle, "zion2", z_dim2);

    adios_write (io_handle, "node-attr", &node);

    adios_close (io_handle);

    printf ("rank: %d append completed\n", rank);

#endif

    MPI_Barrier (MPI_COMM_WORLD);

    adios_finalize (node);
    MPI_Finalize ();

    return 0;
}
