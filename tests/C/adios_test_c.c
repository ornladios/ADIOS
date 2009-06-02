#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>

// mpi
#include "mpi.h"

#include "adios.h"

#include <sys/time.h>
/* Subtract the `struct timeval' values X and Y,
   storing the result in RESULT.
   Return 1 if the difference is negative, otherwise 0.  */

int
timeval_subtract (result, x, y)
     struct timeval *result, *x, *y;
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

#define DO_WRITE 1
#define DO_READ 0
#define DO_APPEND 0
#define MEMORY_THIEF 0

int main (int argc, char ** argv)
{
    char * type_name = "restart";
    char * filename = "restart.bp";
    int64_t io_handle;  // io handle
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    int size;
    struct timeval time_start;
    struct timeval time_end;
    struct timeval time_diff;
    struct timeval * time_diff_all;

    uint64_t byte_test_length = 2LL * 1024 * 1024 * 1024;
    uint64_t memory_thief_length =  1024 * 1024 * 1024  // 1 GB
                             + 128 * 1024 * 1024   // 128 MB
                             - byte_test_length;
printf ("byte_test_length: %llu\n", byte_test_length);
#if DO_WRITE
    char * byte_test = 0;
    byte_test = malloc (byte_test_length);
    if (byte_test == 0)
    {
        fprintf (stderr, "Error allocating memory for write byte_test: %llu\n"
                ,byte_test_length
                );
        exit (-1);
    }
    byte_test [byte_test_length] = 0;
#endif
#if DO_READ
    char * r_byte_test = 0;
    r_byte_test = malloc (byte_test_length);
    if (byte_test == 0)
    {
        fprintf (stderr, "Error allocating memory for read byte_test: %llu\n"
                ,byte_test_length
                );
        exit (-1);
    }
    r_byte_test = malloc (byte_test_length);
    r_byte_test [byte_test_length] = 0;
#endif

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

    // allocate a big block of memory to stymie unwanted local caching
    // that would make the numbers suspect
#if MEMORY_THIEF
    char * memory_thief = malloc (1 * 1024 * 1024 * 1024);
#endif

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    if (!adios_init ("config_c.xml"))
        return -1;

    if (rank == 0)
    time_diff_all = (struct timeval *) malloc (size * sizeof (struct timeval));
    else
        time_diff_all = 0;

#if DO_WRITE
//printf ("XXXXXXXXXXXXXXXX do a write XXXXXXXXXXXXXXXXX\n");
    gettimeofday (&time_start, NULL);
    adios_open (&io_handle, type_name, filename, "w");
    adios_group_size (io_handle, 4 + byte_test_length, &total, &comm);
#if 0
    adios_group_size (io_handle,  4 + 4
                                + 4 * zionsize1
                                + 4 + 4 * zionsize2 * zionsize2
                                + 4 + 4 * zionsize2 * zionsize3
                                + 4
                                + 4 + byte_test_length
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

#endif
printf ("1\n");
    adios_write (io_handle, "byte_test_length", &byte_test_length);
printf ("1\n");
    adios_write (io_handle, "byte_test", byte_test);
printf ("1\n");

    //printf ("A rank: %d write completed\n", rank);
    adios_close (io_handle);
printf ("1\n");
    //printf ("B rank: %d write completed\n", rank);
    gettimeofday (&time_end, NULL);

    //printf ("C rank: %d write completed\n", rank);

    MPI_Barrier (MPI_COMM_WORLD);

    timeval_subtract (&time_diff, &time_end, &time_start);
    MPI_Gather (&time_diff, sizeof (struct timeval), MPI_BYTE
               ,time_diff_all, sizeof (struct timeval), MPI_BYTE
               ,0, MPI_COMM_WORLD
               );
    if (rank == 0)
    {
        memcpy (&(time_diff_all [0]), &time_diff, sizeof (struct timeval));
        int i;
        int max_sec = 0;
        int max_usec = 0;
        int min_sec = INT_MAX;
        int min_usec = INT_MAX;
        printf ("proc\tsec\tusec\n");
        for (i = 0; i < size; i++)
        {
            printf ("%06d\t%02d\t%06d\n", i, time_diff_all [i].tv_sec
                   ,time_diff_all [i].tv_usec
                   );

            if (time_diff_all [i].tv_sec >= max_sec)
            {
                if (   time_diff_all [i].tv_usec > max_usec
                    || time_diff_all [i].tv_sec > max_sec
                   )
                {
                    max_sec = time_diff_all [i].tv_sec;
                    max_usec = time_diff_all [i].tv_usec;
                }
            }

            if (time_diff_all [i].tv_sec <=  min_sec)
            {
                if (   time_diff_all [i].tv_usec < min_usec
                    || time_diff_all [i].tv_sec < min_sec
                   )
                {
                    min_sec = time_diff_all [i].tv_sec;
                    min_usec = time_diff_all [i].tv_usec;
                }
            }
        }

        printf ("max time: %02d.%06d\n", max_sec, max_usec);
        printf ("Aggregate GB/sec: %f\n", (1.0 * byte_test_length * size / (1024 * 1024 * 1024)) / (max_sec + (max_usec / 1000000.0)));

        printf ("\n\n");
    }

    MPI_Barrier (MPI_COMM_WORLD);
#endif

#if DO_READ
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

#if DO_APPEND
for (int i = 0; i < 3; i++)
{
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
}

#endif

    MPI_Barrier (MPI_COMM_WORLD);

    adios_finalize (node);
    MPI_Finalize ();

    free (zion1);
    free (zion2);
    free (zion3);

#if DO_WRITE
    free (byte_test);
#endif
#if DO_READ
    free (r_byte_test);
#endif

#if MEMORY_THIEF
    free (memory_thief);
#endif

    return 0;
}
