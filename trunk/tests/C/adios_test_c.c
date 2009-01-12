#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

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

    int byte_test_length = 2600;
    char byte_test [byte_test_length + 1];
    char r_byte_test [byte_test_length + 1];
    byte_test [byte_test_length] = 0;
    r_byte_test [byte_test_length] = 0;

    int var_x1 = 101;
    int var_x2 = 102;

    int zionsize1 = 10; // fixed size
    int zionsize2 = 2;  // attr_dim
    int zionsize3 = 5;  // attr_dim2

    float * zion1;
    float * zion2;
    float * zion3;

    zion1 = malloc (zionsize1 * sizeof (float));
    zion2 = malloc (zionsize2 * zionsize2 * sizeof (float));
    zion3 = malloc (zionsize2 * zionsize3 * sizeof (float));
    assert (zion1);
    assert (zion2);
    assert (zion3);
    memset (zion1, 0, zionsize1 * sizeof (float));
    memset (zion2, 0, zionsize2 * zionsize2 * sizeof (float));
    memset (zion3, 0, zionsize2 * zionsize3 * sizeof (float));
    int r_var_x1;
    int r_var_x2;
    int r_zsize;
    float r_z [zionsize1];

    int node = 0;

    uint64_t total;
    int i;
    int j;

    zion1 [0] = 10.0;
    zion1 [1] = 11.0;
    zion1 [2] = 12.0;
    zion1 [3] = 13.0;
    zion1 [4] = 14.0;
    zion1 [5] = 15.0;
    zion1 [6] = 16.0;
    zion1 [7] = 17.0;
    zion1 [8] = 18.0;
    zion1 [9] = 19.0;

    for (i = 0; i < 100; i++)
        for (j = 0; j < 26; j++)
            byte_test [i * 26 + j] = 'a' + j;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (!adios_init ("config_c.xml"))
        return -1;

#if 1
printf ("XXXXXXXXXXXXXXXX do a write XXXXXXXXXXXXXXXXX\n");
    adios_open (&io_handle, type_name, filename, "w");
    adios_group_size (io_handle,  4 + 4
                                + 4 * zionsize1
                                + 4 + 4 * zionsize2 * zionsize2
                                + 4 + 4 * zionsize2 * zionsize3
                                + 4
                                + 4 + 2600
                     ,&total, &comm
                     );
    adios_write (io_handle, "/mype", &var_x1);
    adios_write (io_handle, "/test/mype", &var_x2);

    adios_write (io_handle, "zion1", zion1);

    adios_write (io_handle, "zionsize2", &zionsize2);
    adios_write (io_handle, "zion2", zion2);

    adios_write (io_handle, "zionsize3", &zionsize3);
    adios_write (io_handle, "zion3", zion3);

    adios_write (io_handle, "node-attr", &node);

    adios_write (io_handle, "byte_test_length", &byte_test_length);
    adios_write (io_handle, "byte_test", byte_test);

    adios_close (io_handle);

    printf ("rank: %d write completed\n", rank);

    MPI_Barrier (MPI_COMM_WORLD);
#endif

#if 1
printf ("XXXXXXXXXXXXXXXX do a read XXXXXXXXXXXXXXXXX\n");

    adios_open (&io_handle, type_name, filename, "r");
    adios_group_size (io_handle, 0, &total, &comm);
    adios_read (io_handle, "/mype", &r_var_x1, 4);
    adios_read (io_handle, "/test/mype", &r_var_x2, 4);
    adios_read (io_handle, "zionsize2", &r_zsize, 4);
    adios_read (io_handle, "zion1", r_z, 4 * 10);
    adios_read (io_handle, "byte_test", r_byte_test, 26 * 100);
    adios_close (io_handle);

    MPI_Barrier (MPI_COMM_WORLD);

    for (i = 0; i < byte_test_length; i++)
            if (r_byte_test [i] != byte_test [i])
            {
                printf ("byte_test doesn't match %d\n", i);
                printf ("byte_test:\n%s\n", byte_test);
                printf ("r_byte_test:\n%s\n", r_byte_test);
                break;
            }

    if (   var_x1 != r_var_x1
        || var_x2 != r_var_x2
        || r_zsize != zionsize2
        || r_z [0] != zion1 [0]
        || r_z [1] != zion1 [1]
        || r_z [2] != zion1 [2]
        || r_z [3] != zion1 [3]
        || r_z [4] != zion1 [4]
        || r_z [5] != zion1 [5]
        || r_z [6] != zion1 [6]
        || r_z [7] != zion1 [7]
        || r_z [8] != zion1 [8]
        || r_z [9] != zion1 [9]
       )
    {
        int i;
        printf ("rank: %d mismatch in reading\n", rank);
        printf ("r_var_x1: %d (%d)\n", r_var_x1, var_x1);
        printf ("r_var_x2: %d (%d)\n", r_var_x2, var_x2);
        printf ("r_zsize: %d (%d)\n", r_zsize, zionsize1);
        for (i = 0; i < 10; i++)
            printf ("z_dim [%d]: %f (%f)\n", i, r_z [i], zion1 [i]);
    }
    else
    {
        printf ("rank: %d read matches write\n", rank);
    }
#endif

#if 1
printf ("XXXXXXXXXXXXXXXX do an append XXXXXXXXXXXXXXXXX\n");
    var_x1 = 11;
    adios_open (&io_handle, type_name, filename, "a");
    adios_group_size (io_handle,  4 + 4
                                + 4 * zionsize1
                                + 4 + 4 * zionsize2 * zionsize2
                                + 4
                     ,&total, &comm
                     );
    adios_write (io_handle, "/mype", &var_x1);
    adios_write (io_handle, "/test/mype", &var_x2);

    adios_write (io_handle, "zion1", zion1);

    adios_write (io_handle, "zionsize2", &zionsize2);
    adios_write (io_handle, "zion2", zion2);

    adios_write (io_handle, "node-attr", &node);

    adios_close (io_handle);

    printf ("rank: %d append completed\n", rank);

#endif

    MPI_Barrier (MPI_COMM_WORLD);

    adios_finalize (node);
    MPI_Finalize ();

    free (zion1);
    free (zion2);
    free (zion3);

    return 0;
}
