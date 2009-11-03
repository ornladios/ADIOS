#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/vfs.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>
#include <signal.h>
#include <assert.h>

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "buffer.h"

static int adios_adaptive_initialized = 0;

#define DO_INDEX_COLLECTION 0
#define BLAST_OPENS 1

#define PRINT_MESSAGES 0

#define COLLECT_METRICS 1

struct adios_adaptive_data_struct
{
    int f;
    MPI_File fh;
    MPI_Request req;
    MPI_Status status;
    MPI_Comm group_comm;
    int rank;
    int size;

    struct adios_bp_buffer_struct_v1 b;

    struct adios_index_process_group_struct_v1 * old_pg_root;
    struct adios_index_var_struct_v1 * old_vars_root;
    struct adios_index_attribute_struct_v1 * old_attrs_root;

    uint64_t vars_start;
    uint64_t vars_header_size;
    uint64_t biggest_size;     // largest writer's data size (round up)
    uint16_t storage_targets;  // number of storage targets being used
    uint16_t split_groups;     // number of files to split output into
    uint16_t split_group_size; // how big the groups are at max (smaller at end)

    // these correspond to the series of parameters available to this transport
    // and control how it performs
    int16_t files_number;       // # files user is creating as part of output
    int16_t max_storage_targets;// number of OSTs in parallel filesystem
    int16_t max_stripe_count;   // max OSTs per file (system max)
    int16_t min_stripe_count;   // min OSTs per file (user desire)
    int16_t overlap_factor;     // percentage of OSTs to overlap for efficiency
    enum {split_files_unknown
         ,split_files_min
         ,split_files_max
         } split_files_count;   // how many 'split' files to generate
    int16_t split_target_count; // number of OSTs for each 'split' file

    pthread_t coordinator;     // only used if in this proc is coord
    pthread_t sub_coordinator; // only used if in this proc is sub coord
    pthread_t writer;
    pthread_t probe;

    pthread_mutex_t mpi_mutex; // need to make sure only one thread does comm
    pthread_mutex_t writer_mutex; // need to make sure only one thread does comm
    pthread_mutex_t sub_coordinator_mutex;
    pthread_mutex_t coordinator_mutex;

    pthread_cond_t mpi_cond;
    pthread_cond_t writer_cond;
    pthread_cond_t sub_coordinator_cond;
    pthread_cond_t coordinator_cond;

    int groups;          // how many groups we split into
    int group;           // which group are we a member of
    int group_size;      // how large our group is
    int sub_coord_rank;  // what is the rank of our sub coordinator
    int coord_rank;      // what is the rank of the coordinator
    int stripe_size;     // how big each stripe piece is

    struct adios_file_struct * fd; // link to what was passed in
    struct adios_method_struct * method; // link to main struct

    // adaptive support stuff start
    // messaging between threads ([0] is msg, others are parameters)
#define PARAMETER_COUNT 6
    Queue mpi_flag;
    Queue writer_flag;
    Queue sub_coordinator_flag;
    Queue coordinator_flag;

    // we need all three so that we only shutdown the messaging when the
    // most significant process (coord vs. sub coord vs. writer) is finished
    int coord_shutdown_flag;
    int sub_coord_shutdown_flag;
    int writer_shutdown_flag;
    int probe_shutdown_flag; // tell the probe thread to die
};

#if COLLECT_METRICS
// see adios_adaptive_finalize for what each represents
struct timeval_writer
{
    struct timeval t;  // time value
    int pid;            // process id
};

int timeval_writer_compare (const void * left, const void * right)
{
    struct timeval_writer * l = (struct timeval_writer *) left;
    struct timeval_writer * r = (struct timeval_writer *) right;

    if (l->pid < r->pid) return -1;
    if (l->pid == r->pid) return 0;
    if (l->pid > r->pid) return 1;
}

struct timing_metrics
{
    struct timeval t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;
    struct timeval t14, t16, t18, t19, t20, t21, t22, t23, t25, t26;
    struct timeval t27, t28;

    uint64_t write_count; // number used
    uint64_t write_size;  // number allocated
    struct timeval * t24;

    uint64_t do_write_count; // number used
    uint64_t do_write_size;  // number allocated
    struct timeval_writer * t15;  // track in sub coord do_write msgs

    uint64_t write_complete_count; // number used
    uint64_t write_complete_size;  // number allocated
    struct timeval_writer * t17;  // track in sub coord write_complete msgs

    uint64_t mpi_queue_count; // number used
    uint64_t mpi_queue_size;  // number allocated
    struct timeval_writer * t29;  // track in sub coord write_complete msgs

    uint64_t iprobe_count; // number used (t25)
    uint64_t send_count;   // number used (t26)
    uint64_t recv_count;   // number used (t27)
};

static struct timing_metrics timing;

// Subtract the `struct timeval' values X and Y,
// storing the result in RESULT.
// Return 1 if the difference is negative, otherwise 0.
static int timeval_subtract (struct timeval * result
                            ,struct timeval * x, struct timeval * y
                            )
{
  // Perform the carry for the later subtraction by updating y.
  if (x->tv_usec < y->tv_usec)
  {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000)
  {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  // Compute the time remaining to wait.
  // tv_usec is certainly positive.
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  // Return 1 if result is negative.
  return x->tv_sec < y->tv_sec;
}

static void timeval_add (struct timeval * result
                        ,struct timeval * x, struct timeval * y
                        )
{
    result->tv_usec = x->tv_usec + y->tv_usec;
    result->tv_sec = x->tv_sec + y->tv_sec;
    if (result->tv_usec > 1000000)
    {
        result->tv_usec -= 1000000;
        result->tv_sec++;
    }
}

static void print_metric (FILE * f, struct timing_metrics * t, int iteration, int rank, int size, int sub_coord_rank);

static
void print_metrics (struct adios_adaptive_data_struct * md, int iteration)
{
    MPI_Barrier (md->group_comm);
    if (md->rank == 0)
    {
        int i;
        struct timing_metrics * t;
        int * sub_coord_ranks;

        t = malloc (sizeof (struct timing_metrics) * md->size);
        sub_coord_ranks = malloc (sizeof (int) * md->size);

        memcpy (&t [0], &timing, sizeof (struct timing_metrics));

        // get the bulk data
        MPI_Gather (&timing, sizeof (struct timing_metrics), MPI_BYTE
                   ,t, sizeof (struct timing_metrics), MPI_BYTE
                   ,0, md->group_comm
                   );

        // get the write timing
        int * index_sizes = malloc (4 * md->size);
        int * index_offsets = malloc (4 * md->size);
        uint32_t total_size = 0;
        char * recv_buffer = 0;
        char * recv_buffer1 = 0;
        char * recv_buffer2 = 0;
        char * recv_buffer3 = 0;
        for (i = 0; i < md->size; i++)
        {
            index_sizes [i] = t [i].write_count * sizeof (struct timeval);
            index_offsets [i] = total_size;
            total_size += t [i].write_count * sizeof (struct timeval);
        }

        recv_buffer = malloc (total_size + 1);

        MPI_Gatherv (t [0].t24, 0, MPI_BYTE
                    ,recv_buffer, index_sizes, index_offsets, MPI_BYTE
                    ,0, md->group_comm
                    );

        t [0].t24 = timing.t24;
        for (i = 1; i < md->size; i++)
        {
            t [i].t24 = (struct timeval *) (recv_buffer + index_offsets [i]);
        }

        // get the DO_WRITE start times
        total_size = 0;
        for (i = 0; i < md->size; i++)
        {
            index_sizes [i] = t [i].do_write_count
                                        * sizeof (struct timeval_writer);
            index_offsets [i] = total_size;
            total_size += t [i].do_write_count * sizeof (struct timeval_writer);
        }

        recv_buffer1 = malloc (total_size + 1);

        MPI_Gatherv (t [0].t15, 0, MPI_BYTE
                    ,recv_buffer1, index_sizes, index_offsets, MPI_BYTE
                    ,0, md->group_comm
                    );
        t [0].t15 = timing.t15;
        for (i = 1; i < md->size; i++)
        {
            t [i].t15 = (struct timeval_writer *)
                                            (recv_buffer1 + index_offsets [i]);
        }

        // get the WRITE_COMPLETE start times
        total_size = 0;
        for (i = 0; i < md->size; i++)
        {
            index_sizes [i] = t [i].write_complete_count
                                             * sizeof (struct timeval_writer);
            index_offsets [i] = total_size;
            total_size += t [i].write_complete_count
                                             * sizeof (struct timeval_writer);
        }

        recv_buffer2 = malloc (total_size + 1);

        MPI_Gatherv (t [0].t17, 0, MPI_BYTE
                    ,recv_buffer2, index_sizes, index_offsets, MPI_BYTE
                    ,0, md->group_comm
                    );
        t [0].t17 = timing.t17;
        for (i = 1; i < md->size; i++)
        {
            t [i].t17 = (struct timeval_writer *)
                                           (recv_buffer2 + index_offsets [i]);
        }

        // get the mpi queue dequeue times
        total_size = 0;
        for (i = 0; i < md->size; i++)
        {
            index_sizes [i] = t [i].mpi_queue_count
                                            * sizeof (struct timeval_writer);
            index_offsets [i] = total_size;
            total_size += t [i].mpi_queue_count
                                           * sizeof (struct timeval_writer);
        }

        recv_buffer3 = malloc (total_size + 1);

        MPI_Gatherv (t [0].t29, 0, MPI_BYTE
                    ,recv_buffer3, index_sizes, index_offsets, MPI_BYTE
                    ,0, md->group_comm
                    );
        t [0].t29 = timing.t29;
        for (i = 1; i < md->size; i++)
        {
            t [i].t29 = (struct timeval_writer *)
                                          (recv_buffer3 + index_offsets [i]);
        }

        // get the sub coordinator ranks
        sub_coord_ranks [0] = md->sub_coord_rank;
        MPI_Gather (&md->sub_coord_rank, 1, MPI_INTEGER
                   ,sub_coord_ranks, 1, MPI_INTEGER
                   ,0, md->group_comm
                   );

        // print the detailed metrics
        FILE * f = fopen ("adios_metrics", "a");
        for (i = 0; i < md->size; i++)
        {
            print_metric (f, &t [i], iteration, i, md->size, sub_coord_ranks [i]);
        }
        fclose (f);

        free (sub_coord_ranks);
        free (t);
        free (recv_buffer);
        free (recv_buffer1);
        free (recv_buffer2);
        free (recv_buffer3);
        free (index_sizes);
        free (index_offsets);
    }
    else
    {
        // send the bulk data
        MPI_Gather (&timing, sizeof (struct timing_metrics), MPI_BYTE
                   ,&timing, sizeof (struct timing_metrics), MPI_BYTE
                   ,0, md->group_comm
                   );

        // send the write timing
        MPI_Gatherv (timing.t24, timing.write_count * sizeof (struct timeval)
                    ,MPI_BYTE
                    ,0, 0, 0, 0
                    ,0, md->group_comm
                    );

        // send the write start times
        qsort (timing.t15, timing.do_write_count, sizeof (struct timeval_writer)
              ,timeval_writer_compare
              );
        MPI_Gatherv (timing.t15, timing.do_write_count
                                              * sizeof (struct timeval_writer)
                    ,MPI_BYTE
                    ,0, 0, 0, 0
                    ,0, md->group_comm
                    );

        // send the write complete times
        qsort (timing.t17, timing.write_complete_count
              ,sizeof (struct timeval_writer)
              ,timeval_writer_compare
              );
        MPI_Gatherv (timing.t17, timing.write_complete_count
                                              * sizeof (struct timeval_writer)
                    ,MPI_BYTE
                    ,0, 0, 0, 0
                    ,0, md->group_comm
                    );

        // send the write complete times
        qsort (timing.t29, timing.mpi_queue_count
              ,sizeof (struct timeval_writer)
              ,timeval_writer_compare
              );
        MPI_Gatherv (timing.t29, timing.mpi_queue_count
                                             * sizeof (struct timeval_writer)
                    ,MPI_BYTE
                    ,0, 0, 0, 0
                    ,0, md->group_comm
                    );

        // send the sub coordinator rank
        MPI_Gather (&md->sub_coord_rank, 1, MPI_INTEGER
                   ,NULL, 1, MPI_INTEGER
                   ,0, md->group_comm
                   );
    }
    MPI_Barrier (md->group_comm);
}

static void print_metric (FILE * f, struct timing_metrics * t, int iteration, int rank, int size, int sub_coord_rank)
{
    struct timeval diff;
    if (rank == 0)
    {
        timeval_subtract (&diff, &t->t2, &t->t1);
        fprintf (f, "cc\t%2d\tFile create (stripe setup):\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t->t6, &t->t5);
        fprintf (f, "dd\t%2d\tMass file open:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t->t5, &t->t16);
        fprintf (f, "ee\t%2d\tBuild file offsets:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
    }
    if (rank == size - 1)
    {
        timeval_subtract (&diff, &t->t4, &t->t3);
        fprintf (f, "bb\t%2d\tCreate threads:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t->t10, &t->t9);
        fprintf (f, "ff\t%2d\tGlobal index creation:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t->t10, &t->t7);
        fprintf (f, "gg\t%2d\tAll writes complete (w/ local index):\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t->t11, &t->t0);
        fprintf (f, "hh\t%2d\tTotal time:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t->t10, &t->t0);
        fprintf (f, "xx\t%2d\tCoord total time:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
    }
    if (rank == sub_coord_rank)
    {
        timeval_subtract (&diff, &t->t13, &t->t0);
        fprintf (f, "yy\t%2d\tSub coord total time:\t%6d\t%02d.%06d\n"
               ,iteration, rank, diff.tv_sec, diff.tv_usec);
    }

    timeval_subtract (&diff, &t->t13, &t->t12);
    fprintf (f, "ii\t%2d\tLocal index creation:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t11, &t->t14);
    fprintf (f, "jj\t%2d\tThread shutdown:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t22, &t->t21);
    fprintf (f, "kk\t%2d\tshould buffer time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t19, &t->t23);
    fprintf (f, "ll\t%2d\tclose startup time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t19, &t->t0);
    fprintf (f, "mm\t%2d\tsetup time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t14, &t->t20);
    fprintf (f, "nn\t%2d\tcleanup time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t21, &t->t0);
    fprintf (f, "oo\t%2d\topen->should_buffer time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t24 [0], &t->t22);
    fprintf (f, "pp\t%2d\tshould_buffer->write1 time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    int i;
    for (i = 0; i < t->write_count - 1; i++)
    {
        timeval_subtract (&diff, &t->t24 [i + 1], &t->t24 [i]);
        fprintf (f, "qq[%i]\t%2d\twrite1->write2 time:\t%6d\t%02d.%06d\n"
               ,i, iteration, rank, diff.tv_sec, diff.tv_usec);
    }

    int w = 0;
    int c = 0;
    while (w < t->do_write_count && c < t->write_complete_count)
    {
        if (t->t15 [w].pid != t->t17 [c].pid)
        {
            if (t->t15 [w].pid < t->t17 [c].pid)
            {
                fprintf (f, "do_write < complete\n");
                w++;
            }
            else
            {
                fprintf (f, "do_write > complete\n");
                c++;
            }
        }
        else
        {
            timeval_subtract (&diff, &t->t17 [c].t, &t->t15 [w].t);
            fprintf (f, "AAA\t%2d\twrite->complete %d time:\t%6d\t%02d.%06d\n"
                   ,iteration, t->t15 [w].pid, rank, diff.tv_sec, diff.tv_usec);
            w++;
            c++;
        }
    }

    w = 0;
    c = 0;
    while (w < t->do_write_count && c < t->mpi_queue_count)
    {
        if (t->t15 [w].pid != t->t29 [c].pid)
        {
            if (t->t15 [w].pid < t->t29 [c].pid)
            {
                fprintf (f, "do_write (pid: %d) < mpi queue (pid: %d)\n", t->t15 [w].pid, t->t29 [c].pid);
                w++;
            }
            else
            {
                fprintf (f, "do_write (pid: %d) > mpi queue (pid: %d) \n", t->t15 [w].pid, t->t29 [c].pid);
                c++;
            }
        }
        else
        {
            timeval_subtract (&diff, &t->t29 [c].t, &t->t15 [w].t);
            fprintf (f, "BBB\t%2d\tMPI outbound queue %d time:\t%6d\t%02d.%06d\n"
                   ,iteration, t->t15 [w].pid, rank, diff.tv_sec, diff.tv_usec);
            w++;
            c++;
        }
    }

    timeval_subtract (&diff, &t->t23, &t->t24 [t->write_count - 1]);
    fprintf (f, "rr\t%2d\twrite[n]->close start time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    fprintf (f, "ss\t%2d\tMPI_Iprobe time:\t%6d\tcount: %d\t%02d.%06d\n"
           ,iteration, rank, t->iprobe_count, t->t25.tv_sec, t->t25.tv_usec);

    fprintf (f, "tt\t%2d\tMPI_Send time:\t%6d\tcount: %d\t%02d.%06d\n"
           ,iteration, rank, t->send_count, t->t26.tv_sec, t->t26.tv_usec);

    fprintf (f, "uu\t%2d\tMPI_Recv time:\t%6d\tcount: %d\t%02d.%06d\n"
           ,iteration, rank, t->recv_count, t->t27.tv_sec, t->t27.tv_usec);

    fprintf (f, "aa\t%2d\tprocess write time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, t->t8.tv_sec, t->t8.tv_usec);

    timeval_subtract (&diff, &t->t11, &t->t23);
    fprintf (f, "vv\t%2d\tclose start to shutdown time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t28, &t->t23);
    fprintf (f, "ww\t%2d\tclose total time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t->t18, &t->t0);
    fprintf (f, "zz\t%2d\twriter total time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, diff.tv_sec, diff.tv_usec);
}
#endif

struct thread_struct
{
    struct adios_method_struct * method; // link to main struct
    struct adios_file_struct * fd;       // link to what was passed in
    struct adios_adaptive_data_struct * md;
};

struct flag_struct
{
    int source;
    int tag;
    char * msg;
};

struct mpi_flag_struct
{
    int target;
    int tag;
    uint64_t size;
    char * msg;
};

#define COPY_ALL_PARAMS(dst,src) \
{ \
int i; \
for (i = 0; i < PARAMETER_COUNT; i++) \
dst [i] = src [i]; \
}

#define INIT_PARAMS(x) \
{ \
int i; \
for (i = 0; i < PARAMETER_COUNT; i++) \
x [i] = NO_FLAG; \
}

// used to message between threads without a mutex or MPI message
enum MESSAGE_FLAGS
{
     NO_FLAG                 = -1  // no msg
    ,SHUTDOWN_FLAG           = -2  // exit threads
    ,DO_WRITE_FLAG           = -3  // start a write
    ,REGISTER_COMPLETE       = -4  // registration suceed, start write
    ,WRITE_COMPLETE          = -5  // write has finished
    ,ADAPTIVE_WRITE_START    = -6  // start an adaptive write to this group
    ,WRITERS_BUSY            = -7  // all writers are either done or busy
    ,OVERALL_WRITE_COMPLETE  = -8  // all writers are done
    ,SEND_INDEX              = -9  // tell the writers to send their index
    ,REGISTER_FLAG           = -10 // register with the parent
    ,INDEX_SIZE              = -11 // size of index from sub to coord
    ,START_WRITES            = -12 // start the writing process
    ,INDEX_BODY              = -13 // the index contents
};

static
const char * message_to_string (enum MESSAGE_FLAGS m)
{
    static char x [20];
    switch (m)
    {
        case NO_FLAG: return "NO_FLAG";
        case SHUTDOWN_FLAG: return "SHUTDOWN_FLAG";
        case DO_WRITE_FLAG: return "DO_WRITE_FLAG";
        case REGISTER_COMPLETE: return "REGISTER_COMPLETE";
        case WRITE_COMPLETE: return "WRITE_COMPLETE";
        case ADAPTIVE_WRITE_START: return "ADAPTIVE_WRITE_START";
        case WRITERS_BUSY: return "WRITERS_BUSY";
        case OVERALL_WRITE_COMPLETE: return "OVERALL_WRITE_COMPLETE";
        case SEND_INDEX: return "SEND_INDEX";
        case REGISTER_FLAG: return "REGISTER_FLAG";
        case INDEX_SIZE: return "INDEX_SIZE";
        case START_WRITES: return "START_WRITES";
        case INDEX_BODY: return "INDEX_BODY";
        default: sprintf (x, "unknown (%d)", m); return x;
    }
}

// s == 'c' for coordinator
// s == 's' for sub coordinator
static
const char * message_to_string_full (uint64_t * msg, char s)
{
    static char x [100];
    switch (msg [0])
    {
        case NO_FLAG:
            return "NO_FLAG";

        case SHUTDOWN_FLAG:
            return "SHUTDOWN_FLAG";

        case DO_WRITE_FLAG:
            sprintf (x, "DO_WRITE_FLAG (tg: %lld tgr: %lld g: %lld gr: %lld "
                        "offset: %lld)"
                    ,msg [1], msg [2], msg [3], msg [4], msg [5]
                    );
            return x;

        case REGISTER_COMPLETE:
            return "REGISTER_COMPLETE";

        case WRITE_COMPLETE:
            sprintf (x, "WRITE_COMPLETE (tg: %lld wg: %lld eo: %lld "
                        "is: %lld writer: %d)"
                    ,msg [1], msg [2], msg [3], msg [4], msg [5]
                    );
            return x;

        case ADAPTIVE_WRITE_START:
            sprintf (x, "ADAPTIVE_WRITE_START (tg: %lld tscr: %lld o: %lld)"
                    ,msg [1], msg [2], msg [3]
                    );
            return x;

        case WRITERS_BUSY:
            sprintf (x, "WRITERS_BUSY (tg: %lld wg: %lld)", msg [1], msg [2]);
            return x;

        case OVERALL_WRITE_COMPLETE:
            sprintf (x, "OVERALL_WRITE_COMPLETE (eo: %lld)", msg [1]);
            return x;

        case SEND_INDEX:
            return "SEND_INDEX";

        case REGISTER_FLAG:
            sprintf (x, "REGISTER_FLAG (1: %lld 2: %lld)", msg [1], msg [2]);
            return x;

        case INDEX_SIZE:
            sprintf (x, "INDEX_SIZE (g: %lld is: %lld)", msg [1], msg [2]);
            return x;

        case START_WRITES:
            return "START_WRITES";

        case INDEX_BODY:
            return "INDEX_BODY";

        default:
            sprintf (x, "unknown (%d)", msg [0]);
            return x;
    }
}

// used to specify which thread is the target for the MPI messages
enum MPI_TAG
{
     TAG_WRITER          = 0
    ,TAG_SUB_COORDINATOR = 1
    ,TAG_COORDINATOR     = 2
    ,TAG_FILE_OPEN       = 3
    ,TAG_SUB_COORDINATOR_INDEX_BODY = 4
    ,TAG_COORDINATOR_INDEX_BODY = 5
};

static
void calc_groups (int rank, int size, int groups
                 ,int * group, int * group_size, int * sub_coord_rank
                 ,int * coord_rank
                 )
{
    *group_size = size / groups;
    int larger_groups = size % groups;
    if (rank < larger_groups * (*group_size + 1) || !larger_groups)
    {
        if (larger_groups)
            (*group_size)++;
        *group = rank / (*group_size);
        *sub_coord_rank =   *group_size * *group
                          + *group_size - 1;
    }
    else
    {
        *group =     (larger_groups)
                   + (rank - (larger_groups * (*group_size + 1)))
                 / *group_size;
        *sub_coord_rank =   (larger_groups * (*group_size + 1) - 1)
                          + (   (*group - larger_groups + 1)
                              * *group_size
                             );
    }
    if (*group == groups - 1)
        *sub_coord_rank = size - 1;
    *coord_rank = size - 1;
}

static void buffer_write (char ** buffer, uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,const void * data, uint64_t size
                         )
{
    if (*buffer_offset + size > *buffer_size || *buffer == 0)
    {
        char * b = realloc (*buffer, *buffer_offset + size + 1000);
        if (b)
        {
            *buffer = b;
            *buffer_size = (*buffer_offset + size + 1000);
        }
        else
        {
            fprintf (stderr, "Cannot allocate memory in buffer_write.  "
                             "Requested: %llu\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}
// adaptive support stuff end

static void set_stripe_size (struct adios_adaptive_data_struct * md
                            ,struct adios_file_struct * fd
                            ,const char * filename
                            );

static void adios_var_to_comm (const char * comm_name
                              ,enum ADIOS_FLAG host_language_fortran
                              ,void * data
                              ,MPI_Comm * comm
                              )
{
    if (data)
    {
        int t = *(int *) data;

        if (!comm_name)
        {
            if (!t)
            {
                fprintf (stderr, "communicator not provided and none "
                                 "listed in XML.  Defaulting to "
                                 "MPI_COMM_SELF\n"
                        );

                *comm = MPI_COMM_SELF;
            }
            else
            {
                if (host_language_fortran == adios_flag_yes)
                {
                    *comm = MPI_Comm_f2c (t);
                }
                else
                {
                    *comm = *(MPI_Comm *) data;
                }
            }
        }
        else
        {
            if (!strcmp (comm_name, ""))
            {
                if (!t)
                {
                    fprintf (stderr, "communicator not provided and none "
                                     "listed in XML.  Defaulting to "
                                     "MPI_COMM_SELF\n"
                            );

                    *comm = MPI_COMM_SELF;
                }
                else
                {
                    if (host_language_fortran == adios_flag_yes)
                    {
                        *comm = MPI_Comm_f2c (t);
                    }
                    else
                    {
                        *comm = *(MPI_Comm *) data;
                    }
                }
            }
            else
            {
                if (!t)
                {
                    fprintf (stderr, "communicator not provided but one "
                                     "listed in XML.  Defaulting to "
                                     "MPI_COMM_WORLD\n"
                            );

                    *comm = MPI_COMM_WORLD;
                }
                else
                {
                    if (host_language_fortran == adios_flag_yes)
                    {
                        *comm = MPI_Comm_f2c (t);
                    }
                    else
                    {
                        *comm = *(MPI_Comm *) data;
                    }
                }
            }
        }
    }
    else
    {
        fprintf (stderr, "coordination-communication not provided. "
                         "Using MPI_COMM_WORLD instead\n"
                );

        *comm = MPI_COMM_WORLD;
    }
}

void adios_adaptive_init (const char * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_adaptive_data_struct * md =
                  (struct adios_adaptive_data_struct *) method->method_data;
    if (!adios_adaptive_initialized)
    {
        adios_adaptive_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_adaptive_data_struct));
    md = (struct adios_adaptive_data_struct *) method->method_data;
    md->f = -1;
    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->rank = 0;
    md->size = 0;
    md->group_comm = MPI_COMM_NULL;
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
    md->vars_start = 0;
    md->vars_header_size = 0;
    md->biggest_size = 0;
    md->storage_targets = 0;
    md->split_groups = 1;
    md->split_group_size = 0;

    md->files_number = 1;          // always can be 1 by default
    md->max_storage_targets = -1;
    md->max_stripe_count = -1;
    md->min_stripe_count = 1;      // always can do only 1 by default
    md->overlap_factor = 0;        // always can not overlap by default
    md->split_target_count = -1;
    md->split_files_count = split_files_unknown;

    md->group = -1;
    md->sub_coord_rank = -1;
    md->fd = 0;
    md->method = 0;

    queue_init (&md->mpi_flag, free);
    queue_init (&md->writer_flag, free);
    queue_init (&md->sub_coordinator_flag, free);
    queue_init (&md->coordinator_flag, free);

    pthread_mutex_init (&md->mpi_mutex, NULL);
    pthread_mutex_init (&md->writer_mutex, NULL);
    pthread_mutex_init (&md->sub_coordinator_mutex, NULL);
    pthread_mutex_init (&md->coordinator_mutex, NULL);

    pthread_cond_init (&md->mpi_cond, NULL);
    pthread_cond_init (&md->writer_cond, NULL);
    pthread_cond_init (&md->sub_coordinator_cond, NULL);
    pthread_cond_init (&md->coordinator_cond, NULL);

    md->coord_shutdown_flag = NO_FLAG;
    md->sub_coord_shutdown_flag = NO_FLAG;
    md->writer_shutdown_flag = NO_FLAG;
    md->probe_shutdown_flag = NO_FLAG;

    // parse the parameters into key=value segments for optional settings
    if (parameters)
    {
        int param_len = strlen (parameters);
        if (param_len > 0)
        {
            char * p = strdup (parameters);
            char * token = strtok (p, ";");

            while (token)
            {
                char * equal_pos = strchr (token, '=');
                if (!equal_pos)
                {
                    continue;
                }
                int equal = equal_pos - token + 1;
                int len = strlen (token);
                char * key = malloc (len + 1);
                char * value = malloc (len + 1);
                strncpy (key, token, equal);
                key [equal - 1] = 0;
                strncpy (value, equal_pos + 1, len - equal);
                value [len - equal] = 0;

                if (key && value)
                {
                    if (strcasecmp (key, "max_storage_targets") == 0)
                    {
                        int v = atoi (value);
                        md->max_storage_targets = v;
                    }
                    else if (strcasecmp (key, "max_stripe_count") == 0)
                    {
                        int v = atoi (value);
                        md->max_stripe_count = v;
                    }
                    else if (strcasecmp (key, "min_stripe_count") == 0)
                    {
                        int v = atoi (value);
                        md->min_stripe_count = v;
                    }
                    else if (strcasecmp (key, "files_number") == 0)
                    {
                        int v = atoi (value);
                        if (v < 1)
                        {
                            fprintf (stderr, "ADAPTIVE: files_number %d "
                                             "too small. defaulting to 1.\n"
                                    ,v
                                    );

                            v = 1;
                        }
                        md->files_number = v;
                    }
                    else if (strcasecmp (key, "overlap_factor") == 0)
                    {
                        int v = atoi (value);
                        if (v < 0)
                        {
                            fprintf (stderr, "ADAPTIVE: overlap_factor %d "
                                             "too small. defaulting to 0.\n"
                                    ,v
                                    );

                            v = 0;
                        } else
                        if (v > 99)
                        {
                            fprintf (stderr, "ADAPTIVE: overlap_factor %d "
                                             "too large. defaulting to 99.\n"
                                    ,v
                                    );

                            v = 99;
                        }
                        md->overlap_factor = v;
                    }
                    else if (strcasecmp (key, "split_target_count") == 0)
                    {
                        int v = atoi (value);
                        md->split_target_count = v;
                    }
                    else if (strcasecmp (key, "split_files_count") == 0)
                    {
                        if (strcasecmp (value, "min") == 0)
                        {
                            md->split_files_count = split_files_min;
                        }
                        else if (strcasecmp (value, "max") == 0)
                        {
                            md->split_files_count = split_files_max;
                        }
                        else
                        {
                            fprintf (stderr, "unkown split_files_count: %s "
                                             "(parameter ignored)\n"
                                    ,value
                                    );
                            md->split_files_count = split_files_unknown;
                        }
                    }
                    else
                    {
                        fprintf (stderr, "ADAPTIVE parameter: key: {%s} "
                                         "value: {%s} not recognized. Ignored\n"
                               ,key, value
                               );
                    }
                }

                free (key);
                free (value);

                token = strtok (NULL, ";");
            }

            free (p);
        }
    }

    adios_buffer_struct_init (&md->b);
}

static void * sub_coordinator_main (void * param);
static void * coordinator_main (void * param);
static void * writer_main (void * param);
static void * probe_main (void * param);

int adios_adaptive_open (struct adios_file_struct * fd
                   ,struct adios_method_struct * method
                   )
{
    struct adios_adaptive_data_struct * md =
                   (struct adios_adaptive_data_struct *) method->method_data;

#if COLLECT_METRICS
    gettimeofday (&timing.t0, NULL);
#endif
    // we have to wait for the group_size (should_buffer) to get the comm
    // before we can do an open for any of the modes
    md->fd = fd;
    md->method = method;

    // reset our shutdown state to not shut down yet
    md->coord_shutdown_flag = NO_FLAG;
    md->sub_coord_shutdown_flag = NO_FLAG;
    md->writer_shutdown_flag = NO_FLAG;
    md->probe_shutdown_flag = NO_FLAG;

#if COLLECT_METRICS
    timing.write_count = 0;
    timing.write_size = 0;
    if (timing.t24) free (timing.t24);
    timing.t24 = 0;

    timing.do_write_count = 0;
    timing.do_write_size = 0;
    if (timing.t15) free (timing.t15);
    timing.t15 = 0;

    timing.write_complete_count = 0;
    timing.write_complete_size = 0;
    if (timing.t17) free (timing.t17);
    timing.t17 = 0;

    timing.mpi_queue_count = 0;
    timing.mpi_queue_size = 0;
    if (timing.t29) free (timing.t29);
    timing.t29 = 0;

    timing.iprobe_count = 0;
    memset (&timing.t25, 0, sizeof (struct timeval));

    timing.send_count = 0;
    memset (&timing.t26, 0, sizeof (struct timeval));

    timing.recv_count = 0;
    memset (&timing.t27, 0, sizeof (struct timeval));
#endif

    return 1;
}

static
void build_offsets (struct adios_bp_buffer_struct_v1 * b
                   ,MPI_Offset * offsets, int size, char * group_name
                   ,struct adios_index_process_group_struct_v1 * pg_root
                   )
{
    while (pg_root)
    {
        if (!strcasecmp (pg_root->group_name, group_name))
        {
            MPI_Offset size = 0;

            if (pg_root->next)
            {
                size = pg_root->next->offset_in_file - pg_root->offset_in_file;
            }
            else
            {
                size = b->pg_index_offset - pg_root->offset_in_file;
            }

            offsets [pg_root->process_id * 3] = pg_root->offset_in_file;
            offsets [pg_root->process_id * 3 + 1] = size;
            offsets [pg_root->process_id * 3 + 2] = b->version;
        }

        pg_root = pg_root->next;
    }
}

static void
adios_build_file_offset (struct adios_adaptive_data_struct *md
                        ,struct adios_file_struct *fd, char *name)
{
#define SCATTER_PARAMS 5
    if (md->group_comm != MPI_COMM_NULL)
    {
        int err;
        if (md->rank == 0)
        {
            // make one space for offset and one for size
            uint64_t * offsets = malloc(sizeof (uint64_t)
                                           * md->size * SCATTER_PARAMS);
            int i;

            offsets [0] = fd->write_size_bytes;
            MPI_Gather (offsets, 1, MPI_LONG_LONG
                       ,offsets, 1, MPI_LONG_LONG
                       ,0, md->group_comm);

            // find the largest and use that as a basis for stripe
            // size for each process writing
            uint64_t biggest_size = 0;
            for (i = 0; i < md->size; i++)
            {
                if (offsets [i] > biggest_size)
                    biggest_size = offsets [i];
            }
            // now round up to the next stripe size increment (Lustre: 64 KB)
#define STRIPE_INCREMENT (64 * 1024)
            if (biggest_size % (STRIPE_INCREMENT))
            {
                biggest_size = (  ((biggest_size / STRIPE_INCREMENT) + 1)
                                * STRIPE_INCREMENT
                               );
            }
            // also round up the base_offset, just in case
            if (fd->base_offset % (STRIPE_INCREMENT))
            {
                fd->base_offset = (  ((fd->base_offset / STRIPE_INCREMENT) + 1)
                                   * STRIPE_INCREMENT
                                  );
            }
            md->biggest_size = biggest_size;
#undef STRIPE_INCREMENT
            md->groups = (md->size > md->max_storage_targets)
                                         ? md->max_storage_targets : md->size;
            md->split_groups = md->groups;
            calc_groups (md->rank, md->size, md->groups
                        ,&md->group, &md->group_size, &md->sub_coord_rank
                        ,&md->coord_rank
                        );
            set_stripe_size (md, fd, name);
            offsets [1] = biggest_size;
            offsets [2] = md->storage_targets;
            offsets [3] = md->split_groups;
            // need to use this calc since last grouping may be incomplete
            // with a previous group having a piece later in the file.
            //md->b.pg_index_offset =   last_offset + biggest_size;

            md->stripe_size = biggest_size;
            err = MPI_Bcast (offsets, SCATTER_PARAMS, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            assert (err == MPI_SUCCESS);
            fd->base_offset = offsets [0];
            fd->pg_start_in_file = fd->base_offset;
            free (offsets);
        }
        else
        {
            uint64_t offset [SCATTER_PARAMS];
            offset [0] = fd->write_size_bytes;
            int i;
            for (i = 0; i < SCATTER_PARAMS; i++) offset [i] = 0;

            MPI_Gather (offset, 1, MPI_LONG_LONG
                       ,offset, 1, MPI_LONG_LONG
                       ,0, md->group_comm
                       );

            err = MPI_Bcast (offset, SCATTER_PARAMS, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            assert (err == MPI_SUCCESS);

            md->biggest_size = offset [0];
            md->stripe_size = offset [1];
            md->storage_targets = offset [2];
            md->split_groups = offset [3];
            md->groups = offset [3];
            md->groups = (md->size > md->max_storage_targets) ?
                                           md->max_storage_targets : md->size;
            calc_groups (md->rank, md->size, md->groups
                        ,&md->group, &md->group_size, &md->sub_coord_rank
                        ,&md->coord_rank
                        );

            fd->pg_start_in_file = fd->base_offset;
        }
    }
    else
    {
        md->b.pg_index_offset = fd->write_size_bytes;
    }
}

// calc_stripe_info figures out how many OSTs to use based on either an
// assumption of one file for all processes or a different count if the
// parameters are supplied in the XML file. These parameters are defined as
// follows:

// files_number - the number of files being written simultaneously. This also
//                implies that all MPI processes are involved in the write.
// max_storage_targets - the number of OSTs in the system.
//                       On ewok, this is 12. On jaguarpf, this is 672.
// max_stripe_count - the maximum number of OSTs available for a single file.
//                    This is 160 on jaguarpf and 12 on ewok.
// min_stripe_count - the fewest OSTs to use per file. This will override the
//                    overlap_factor, if necessary.
// overlap_factor - what percentage of the allocated portion of the OSTs should
//                  be allowed to overlap with the next file set. For example,
//                  a value of 50 means that half of the next set of OSTs will
//                  also be used for this set (e.g., set 0= 0-15, set 1= 10-25,
//                  set 2=20-35, ...).

// The way to put it into the XML, which is currently pretty unforgiving, is
// like this:
// <transport method="ADAPTIVE" group="restart">max_storage_targets=672;max_stripe_count=160;files_number=3;overlap_factor=50;min_stripe_count=10</transport>

// This says to divide the 672 OSTs so that there are 16 groups. Each group
// will consist of 1/3 portion of the whole plus 1/6 as an overlap factor.
// This will not exceed the max_stripe_count. (672/3 = 224 + 50% = 336, but
// max is 160 so limited to 160 [0-159, 112-272, 225-385, etc.]).

// If this were set for 600 files, then it would be 10 OSTs per file at an
// offset of 2 [0-9, 2-11, 4-13, etc.]

// Initial assumption is that based on the rank, we can guess which set of OSTs
// to use. If this won't work, then we need to add a communication to exchange
// information so that the processes can make that decision (and MPI_All_to_all
// of the rank is sufficient so that each main process can determine the
// ordering and then calculate which set to use).

// Once we figure out the stripe info, we may determine that we are not using
// all possible OSTs and therefore limiting parallelism in writing. To address
// this, two other options may be supplied
// split_target_count - the number of OSTs to use if maximizing the number
//                      of files we create.
// split_files_count - {min|max} Either maximize the number of files (using
//                     the split_target_count to guide the number (no shared
//                     OSTs) or minimize while still using all of the OSTs.
static void calc_stripe_info (struct adios_adaptive_data_struct * md
                             ,int * stripe_offset
                             ,int * stripe_count
                             )
{
    // these 2 are the required parameters to decide to do stripe manipulation
    if (   md->max_storage_targets > 0
        && md->max_stripe_count > 0
       )
    {
        int targets_per_file = md->max_storage_targets / md->files_number;
        if (md->max_storage_targets % md->files_number)
            targets_per_file++;

        int overlap_count = (int) (  targets_per_file
                                   * (md->overlap_factor/100.0)
                                  );

        int net_targets_per_file = targets_per_file + overlap_count;

        if (net_targets_per_file < md->min_stripe_count)
            net_targets_per_file = md->min_stripe_count;
        if (net_targets_per_file > md->max_stripe_count)
            net_targets_per_file = md->max_stripe_count;

        int range = 0;
        int number = 0;

        range = md->size / md->files_number;
        if (md->size % md->files_number)
            range++;
        // leave number the same so that we don't falsely overlap
        number = md->rank / range;

        // offset should start based on the split position or number of files
        // if we split, then this is the base and size of this range of files
        *stripe_offset = targets_per_file * number;
        *stripe_count = net_targets_per_file;

        // split, if we must
        if (md->split_files_count != split_files_unknown)
        {
            switch (md->split_files_count)
            {
                case split_files_min:
                {
                     // if we can split to our advantage, do so
                     if (  net_targets_per_file * md->files_number
                         < md->max_storage_targets
                        )
                     {
                         int groups = md->max_storage_targets / *stripe_count;
                         if (md->max_storage_targets % *stripe_count)
                             groups++;
                         if (md->max_storage_targets > md->size)
                         {
                             groups = md->size / *stripe_count;
                             if (md->size % *stripe_count)
                                 groups++;
                         }
                         int group = md->rank / (md->size / groups);
                         int group_size = md->size / groups;
                         int group_rank = md->rank % group_size;

                         md->split_groups = groups;
                         *stripe_count = group_size;
                     }
                     break;
                }

                case split_files_max:
                {
                    // if we can split to our advantage, do so
                    if (net_targets_per_file > md->split_target_count)
                    {
                        int groups =   md->max_storage_targets
                                     / md->split_target_count;
                        if (md->max_storage_targets % md->split_target_count)
                            groups++;
                        if (md->max_storage_targets > md->size)
                        {
                            groups = md->size / md->split_target_count;
                            if (md->size % md->split_target_count)
                                groups++;
                        }
                        int group = md->rank / (md->size / groups);

                        md->split_groups = groups;
                        *stripe_count = md->split_target_count;
                    }
                    break;
                }
            }
        }
    }
    else
    {
        *stripe_count = UINT16_MAX;
    }
}

// LUSTRE Structure
// from /usr/include/lustre/lustre_user.h
#define LUSTRE_SUPER_MAGIC 0x0BD00BD0
#  define LOV_USER_MAGIC 0x0BD10BD0
#  define LL_IOC_LOV_SETSTRIPE  _IOW ('f', 154, long)
#  define LL_IOC_LOV_GETSTRIPE  _IOW ('f', 155, long)
#define O_LOV_DELAY_CREATE 0100000000

struct lov_user_ost_data {           // per-stripe data structure
        uint64_t l_object_id;        // OST object ID
        uint64_t l_object_gr;        // OST object group (creating MDS number)
        uint32_t l_ost_gen;          // generation of this OST index
        uint32_t l_ost_idx;          // OST index in LOV
} __attribute__((packed));
struct lov_user_md {                 // LOV EA user data (host-endian)
        uint32_t lmm_magic;          // magic number = LOV_USER_MAGIC_V1
        uint32_t lmm_pattern;        // LOV_PATTERN_RAID0, LOV_PATTERN_RAID1
        uint64_t lmm_object_id;      // LOV object ID
        uint64_t lmm_object_gr;      // LOV object group
        uint32_t lmm_stripe_size;    // size of stripe in bytes
        uint16_t lmm_stripe_count;   // num stripes in use for this object
        uint16_t lmm_stripe_offset;  // starting stripe offset in lmm_objects
        struct lov_user_ost_data  lmm_objects[0]; // per-stripe data
} __attribute__((packed));

// do the magic ioctl calls to set Lustre's stripe size
// generate the number of groups, group size, and number of OSTs per file.
// These are returned in md->split_groups, md->split_size, & md->storage_targets
static void set_stripe_size (struct adios_adaptive_data_struct * md
                            ,struct adios_file_struct * fd
                            ,const char * filename
                            )
{
    struct statfs fsbuf;
    int err;

    int f;
    int old_mask;
    int perm;

#if COLLECT_METRICS
    gettimeofday (&timing.t1, NULL);
#endif
    old_mask = umask (0);
    umask (old_mask);
    perm = old_mask ^ 0666;

    if (fd->mode == adios_mode_write)
        unlink (filename); // cleanup old stuff
    f = open (filename, O_RDONLY | O_CREAT | O_LOV_DELAY_CREATE, perm);
    // Note: Since each file might have different write_buffer,
    // So we will reset write_buffer even buffer_size != 0
    err = statfs (filename, &fsbuf);
    if (!err && fsbuf.f_type == LUSTRE_SUPER_MAGIC)
    {
        if (f != -1)
        {
            int i;
            int stripe_count;
            int stripe_offset;
            struct lov_user_md lum;
            lum.lmm_magic = LOV_USER_MAGIC;
            // get what Lustre assigns by default
            err = ioctl (f, LL_IOC_LOV_GETSTRIPE, (void *) &lum);
            stripe_count = lum.lmm_stripe_count;
            stripe_offset = lum.lmm_stripe_offset;
            calc_stripe_info (md, &stripe_offset, &stripe_count);

            // fixup for our desires
            lum.lmm_magic = LOV_USER_MAGIC;
            lum.lmm_pattern = 0;
            lum.lmm_stripe_size = md->biggest_size;
            lum.lmm_stripe_count = stripe_count; // from calc_stripe_info
            lum.lmm_stripe_offset = stripe_offset;
            err = ioctl (f, LL_IOC_LOV_SETSTRIPE, (void *) &lum);
            lum.lmm_stripe_count = 0;
            err = ioctl (f, LL_IOC_LOV_GETSTRIPE, (void *) &lum);
            // if err != 0, the must not be Lustre
            if (err == 0)
            {
                md->storage_targets = lum.lmm_stripe_count;
            }
            close (f);

            md->storage_targets = 1;

            int * f_split = malloc (sizeof (int) * md->split_groups);
            char split_format [7] = ".%d";
            char split_name [7];
            char * new_name = malloc (strlen (filename) + 7 + 1);
            for (i = 0; i < md->split_groups; i++)
            {
                sprintf (split_name, split_format, i);
                sprintf (new_name, "%s%s", filename, split_name);
                if (fd->mode == adios_mode_write)
                    unlink (new_name);  // clean up old stuff
                f_split [i] = open (new_name
                                   ,O_RDONLY | O_CREAT | O_LOV_DELAY_CREATE
                                   ,perm
                                   );
            }

            for (i = 0; i < md->split_groups; i++)
            {
                lum.lmm_magic = LOV_USER_MAGIC;
                lum.lmm_pattern = 0;
                lum.lmm_stripe_size = md->biggest_size;
                lum.lmm_stripe_count = md->storage_targets;
                lum.lmm_stripe_offset =   stripe_offset
                                        + i * md->storage_targets;
                err = ioctl (f_split [i], LL_IOC_LOV_SETSTRIPE
                            ,(void *) &lum
                            );
                close (f_split [i]);
            }
            free (f_split);
            unlink (filename);
        }
    }
#if COLLECT_METRICS
    gettimeofday (&timing.t2, NULL);
#endif
}

static void setup_threads (struct thread_struct * t)
{
    struct adios_adaptive_data_struct * md = t->md;

    uint64_t msg [PARAMETER_COUNT];
    MPI_Request req;
    MPI_Status status;
#if COLLECT_METRICS
    gettimeofday (&timing.t3, NULL);
#endif

    // we can get the threads setup and get everyone registered though.
    if (md->rank == md->coord_rank)
    {
        pthread_create (&md->coordinator, NULL, coordinator_main, (void *) t);
    }
    if (md->rank == md->sub_coord_rank)
    {
        pthread_create (&md->sub_coordinator, NULL, sub_coordinator_main
                       ,(void *) t
                       );
    }
    pthread_create (&md->writer, NULL, writer_main, (void *) t);
    pthread_create (&md->probe, NULL, probe_main, (void *) t);
#if COLLECT_METRICS
    gettimeofday (&timing.t4, NULL);
#endif
}

enum ADIOS_FLAG adios_adaptive_should_buffer (struct adios_file_struct * fd
                                        ,struct adios_method_struct * method
                                        ,void * comm
                                        )
{
    int i;
    struct adios_adaptive_data_struct * md =
                   (struct adios_adaptive_data_struct *) method->method_data;
    char * name;
    int err;
    int flag;    // used for coordinating the MPI_File_open

    int previous;
    int current;
    int next;
#if COLLECT_METRICS
    gettimeofday (&timing.t21, NULL);
#endif

    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);

    adios_var_to_comm (fd->group->group_comm
                      ,fd->group->adios_host_language_fortran
                      ,comm
                      ,&md->group_comm
                      );
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }
    fd->group->process_id = md->rank;

    if (md->rank == md->size - 1)
        next = -1;
    else
        next = md->rank + 1;
    previous = md->rank - 1;
    current = md->rank;

    fd->base_offset = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL, &md->fh
                                    );
                if (err != MPI_SUCCESS)
                {
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    fprintf (stderr, "MPI open read failed for %s: '%s'\n"
                            ,name, e
                            );
                    free (name);

                    return adios_flag_no;
                }

                MPI_Offset file_size;
                MPI_File_get_size (md->fh, &file_size);
                md->b.file_size = file_size;

                adios_init_buffer_read_version (&md->b);
                MPI_File_seek (md->fh, md->b.file_size - md->b.length
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.length, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_version (&md->b, &md->b.version);

                adios_init_buffer_read_index_offsets (&md->b);
                // already in the buffer
                adios_parse_index_offsets_v1 (&md->b);

                adios_init_buffer_read_process_group_index (&md->b);
                MPI_File_seek (md->fh, md->b.pg_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.pg_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_process_group_index_v1 (&md->b
                                                   ,&md->old_pg_root
                                                   );

#if 1
                adios_init_buffer_read_vars_index (&md->b);
                MPI_File_seek (md->fh, md->b.vars_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

                adios_init_buffer_read_attributes_index (&md->b);
                MPI_File_seek (md->fh, md->b.attrs_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.attrs_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_attributes_index_v1 (&md->b, &md->old_attrs_root);
#endif

                fd->base_offset = md->b.end_of_pgs;
            }

            if (   md->group_comm != MPI_COMM_NULL
                && md->group_comm != MPI_COMM_SELF
               )
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size * 3
                                                  );
                    memset (offsets, 0, sizeof (MPI_Offset) * md->size * 3);

                    // go through the pg index to build the offsets array
                    build_offsets (&md->b, offsets, md->size
                                  ,fd->group->name, md->old_pg_root
                                  );
                    MPI_Scatter (offsets, 3, MPI_LONG_LONG
                                ,offsets, 3, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
                    md->b.read_pg_offset = offsets [0];
                    md->b.read_pg_size = offsets [1];
                    free (offsets);
                }
                else
                {
                    MPI_Offset offset [3];
                    offset [0] = offset [1] = offset [2] = 0;

                    MPI_Scatter (offset, 3, MPI_LONG_LONG
                                ,offset, 3, MPI_LONG_LONG
                                ,0, md->group_comm
                                );

                    md->b.read_pg_offset = offset [0];
                    md->b.read_pg_size = offset [1];
                    md->b.version = offset [2];
                }
            }

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                // note rank 0 is already open
                // don't open it again here

                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mpi_mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mpi_mutex);
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mpi_mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mpi_mutex);
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }
            if (next != -1)
                MPI_Wait (&md->req, &md->status);

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                        ,name, e
                        );
                free (name);

                return adios_flag_no;
            }

            break;
        }

        case adios_mode_write:
        {
            fd->base_offset = 0;
            fd->pg_start_in_file = 0;

            // this needs to be before our build_file_offsets because
            // that process will attempt to create the file from scratch
            // and set the stripe parameters
            if (previous == -1)
            {
                unlink (name); // make sure clean
            }

            // figure out the offsets and create the file with proper striping
            // before the MPI_File_open is called
#if COLLECT_METRICS
            gettimeofday (&timing.t16, NULL);
#endif
            adios_build_file_offset (md, fd, name);

            if (md->split_groups != 1)
            {
                // if we need to do a file split, we need to fixup the name
                name = realloc (name, (  strlen (method->base_path)
                                       + strlen (fd->name) + 1 + 6
                                      )
                               ); // 6 extra for '.XXXXX' file number
                char split_format [7] = ".%d";
                char split_name [7];
                sprintf (split_name, split_format, md->group);
                strcat (name, split_name);
            }
#if COLLECT_METRICS
            gettimeofday (&timing.t5, NULL);
#endif
#if BLAST_OPENS
            md->f = open (name, O_WRONLY | O_LARGEFILE);
#else
            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                md->f = open (name, O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC);
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mpi_mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, TAG_FILE_OPEN
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mpi_mutex);
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INTEGER, previous, TAG_FILE_OPEN
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mpi_mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, TAG_FILE_OPEN
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mpi_mutex);
                }
                md->f = open (name, O_WRONLY | O_LARGEFILE | O_CREAT | O_TRUNC);
            }
            if (next != -1)
                MPI_Wait (&md->req, &md->status);
#endif
            if (md->f == -1)
            {
                printf ("File open error for %s: %s\n", name, strerror (errno));
                free (name);

                return adios_flag_no;
            }
#if COLLECT_METRICS
            gettimeofday (&timing.t6, NULL);
#endif

            break;
        }

        case adios_mode_append:
        {
            int old_file = 1;
            adios_buffer_struct_clear (&md->b);

            if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL, &md->fh
                                    );

                if (err != MPI_SUCCESS)
                {
                    old_file = 0;
                    err = MPI_File_open (MPI_COMM_SELF, name
                                        ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                        ,MPI_INFO_NULL, &md->fh
                                        );
                    if (err != MPI_SUCCESS)
                    {
                        char e [MPI_MAX_ERROR_STRING];
                        int len = 0;
                        memset (e, 0, MPI_MAX_ERROR_STRING);
                        MPI_Error_string (err, e, &len);
                        fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                                ,name, e
                                );
                        free (name);

                        return adios_flag_no;
                    }
                }
                MPI_Bcast (&old_file, 1, MPI_INTEGER, 0, md->group_comm);
            }
            else
            {
                if (md->group_comm != MPI_COMM_NULL)
                    MPI_Bcast (&old_file, 1, MPI_INTEGER, 0, md->group_comm);
            }

            if (old_file)
            {
                if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
                {
                    if (err != MPI_SUCCESS)
                    {
                        md->b.file_size = 0;
                    }
                    else
                    {
                        MPI_Offset file_size;
                        MPI_File_get_size (md->fh, &file_size);
                        md->b.file_size = file_size;
                    }

                    adios_init_buffer_read_version (&md->b);
                    MPI_File_seek (md->fh, md->b.file_size - md->b.length
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.length, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_version (&md->b, &md->b.version);

                    adios_init_buffer_read_index_offsets (&md->b);
                    // already in the buffer
                    adios_parse_index_offsets_v1 (&md->b);

                    adios_init_buffer_read_process_group_index (&md->b);
                    MPI_File_seek (md->fh, md->b.pg_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.pg_size, MPI_BYTE
                                  ,&md->status
                                  );

                    adios_parse_process_group_index_v1 (&md->b
                                                       ,&md->old_pg_root
                                                       );

                    // find the largest time index so we can append properly
                    struct adios_index_process_group_struct_v1 * p;
                    uint32_t max_time_index = 0;
                    p = md->old_pg_root;
                    while (p)
                    {
                        if (p->time_index > max_time_index)
                            max_time_index = p->time_index;
                        p = p->next;
                    }
                    fd->group->time_index = ++max_time_index;
                    MPI_Bcast (&fd->group->time_index, 1, MPI_INTEGER, 0
                              ,md->group_comm
                              );

                    adios_init_buffer_read_vars_index (&md->b);
                    MPI_File_seek (md->fh, md->b.vars_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

                    adios_init_buffer_read_attributes_index (&md->b);
                    MPI_File_seek (md->fh, md->b.attrs_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.attrs_size
                                  ,MPI_BYTE, &md->status
                                  );
                    adios_parse_attributes_index_v1 (&md->b
                                                    ,&md->old_attrs_root
                                                    );

                    fd->base_offset = md->b.end_of_pgs;
                    fd->pg_start_in_file = fd->base_offset;
                }
                else
                {
                    fd->base_offset = 0;
                    fd->pg_start_in_file = 0;
                    MPI_Bcast (&fd->group->time_index, 1, MPI_INTEGER, 0
                              ,md->group_comm
                              );

                }

                MPI_File_close (&md->fh);
            }
            else
            {
                fd->base_offset = 0;
                fd->pg_start_in_file = 0;
            }

            // figure out the offsets and create the file with proper striping
            // before the MPI_File_open is called
            adios_build_file_offset (md, fd, name);

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                // we know it exists, because we created it if it didn't
                // when reading the old file so can just open wronly
                // but adding the create for consistency with write mode
                // so it is easier to merge write/append later
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mpi_mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mpi_mutex);
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    //pthread_mutex_lock (&md->mpi_mutex);
                    MPI_Isend (&flag, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                    //pthread_mutex_unlock (&md->mpi_mutex);
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }
            if (next != -1)
                MPI_Wait (&md->req, &md->status);

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                        ,name, e
                        );
                free (name);

                return adios_flag_no;
            }

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);

            free (name);

            return adios_flag_no;
        }
    }

    free (name);

    if (fd->shared_buffer == adios_flag_no && fd->mode != adios_mode_read)
    {
        // write the process group header
        adios_write_process_group_header_v1 (fd, fd->write_size_bytes);

        MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
        MPI_File_write (md->fh, fd->buffer, fd->bytes_written, MPI_BYTE
                       ,&md->status
                       );
        int count;
        MPI_Get_count (&md->status, MPI_BYTE, &count);
        if (count != fd->bytes_written)
        {
            fprintf (stderr, "a:MPI method tried to write %llu, "
                             "only wrote %d\n"
                    ,fd->bytes_written
                    ,count
                    );
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);

        // setup for writing vars
        adios_write_open_vars_v1 (fd);
        md->vars_start = fd->base_offset;
        md->vars_header_size = fd->offset;
        fd->base_offset += fd->offset;
        MPI_File_seek (md->fh, md->vars_header_size, MPI_SEEK_CUR);
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);
    }

#if COLLECT_METRICS
    gettimeofday (&timing.t22, NULL);
#endif
    return fd->shared_buffer;
}

void adios_adaptive_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,void * data
                     ,struct adios_method_struct * method
                     )
{
    struct adios_adaptive_data_struct * md =
                     (struct adios_adaptive_data_struct *) method->method_data;

    if (v->got_buffer == adios_flag_yes)
    {
        if (data != v->data)  // if the user didn't give back the same thing
        {
            if (v->free_data == adios_flag_yes)
            {
                free (v->data);
                adios_method_buffer_free (v->data_size);
            }
        }
        else
        {
            // we already saved all of the info, so we're ok.
            return;
        }
    }

    if (fd->shared_buffer == adios_flag_no)
    {
        // var payload sent for sizing information
        adios_write_var_header_v1 (fd, v);

        MPI_File_write (md->fh, fd->buffer, fd->bytes_written
                       ,MPI_BYTE, &md->status
                       );
        int count;
        MPI_Get_count (&md->status, MPI_BYTE, &count);
        if (count != fd->bytes_written)
        {
            fprintf (stderr, "b:MPI method tried to write %llu, "
                             "only wrote %d\n"
                    ,fd->bytes_written
                    ,count
                    );
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);

        // write payload
        // adios_write_var_payload_v1 (fd, v);
        uint64_t var_size = adios_get_var_size (v, fd->group, v->data);
        MPI_File_write (md->fh, v->data, var_size, MPI_BYTE, &md->status);
        MPI_Get_count (&md->status, MPI_BYTE, &count);
        if (count != var_size)
        {
            fprintf (stderr, "c:MPI method tried to write %llu, "
                             "only wrote %d\n"
                    ,var_size
                    ,count
                    );
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);
    }
#if COLLECT_METRICS
    if (timing.write_count == timing.write_size)
    {
        timing.write_size += 10;
        timing.t24 = realloc (timing.t24, sizeof (struct timeval)
                                          * timing.write_size
                             );
        assert (timing.t24 != 0);
    }
    gettimeofday (&(timing.t24 [timing.write_count++]), NULL);
#endif
}

void adios_adaptive_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
    uint64_t mem_allowed;

    if (*size == 0)
    {
        *buffer = 0;

        return;
    }

    if (v->data && v->free_data)
    {
        adios_method_buffer_free (v->data_size);
        free (v->data);
    }

    mem_allowed = adios_method_buffer_alloc (*size);
    if (mem_allowed == *size)
    {
        *buffer = malloc (*size);
        if (!*buffer)
        {
            adios_method_buffer_free (mem_allowed);
            fprintf (stderr, "Out of memory allocating %llu bytes for %s\n"
                    ,*size, v->name
                    );
            v->got_buffer = adios_flag_no;
            v->free_data = adios_flag_no;
            v->data_size = 0;
            v->data = 0;
            *size = 0;
            *buffer = 0;
        }
        else
        {
            v->got_buffer = adios_flag_yes;
            v->free_data = adios_flag_yes;
            v->data_size = mem_allowed;
            v->data = *buffer;
        }
    }
    else
    {
        adios_method_buffer_free (mem_allowed);
        fprintf (stderr, "OVERFLOW: Cannot allocate requested buffer of %llu "
                         "bytes for %s\n"
                ,*size
                ,v->name
                );
        *size = 0;
        *buffer = 0;
    }
}

void adios_adaptive_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )
{
    v->data = buffer;
    v->data_size = buffer_size;
}

static void adios_mpi_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
    struct adios_adaptive_data_struct * md =
                    (struct adios_adaptive_data_struct *) method->method_data;
    struct adios_var_struct * v = fd->group->vars;

    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer = 0;
    data.buffer_len = 0;

    switch (md->b.version)
    {
        case 1:
        {
            // the three section headers
            struct adios_process_group_header_struct_v1 pg_header;
            struct adios_vars_header_struct_v1 vars_header;
            struct adios_attributes_header_struct_v1 attrs_header;

            struct adios_var_header_struct_v1 var_header;
            struct adios_var_payload_struct_v1 var_payload;
            struct adios_attribute_struct_v1 attribute;

            int i;

            adios_init_buffer_read_process_group (&md->b);
            MPI_File_seek (md->fh, md->b.read_pg_offset
                          ,MPI_SEEK_SET
                          );
            MPI_File_read (md->fh, md->b.buff, md->b.read_pg_size, MPI_BYTE
                          ,&md->status
                          );
            adios_parse_process_group_header_v1 (&md->b, &pg_header);

            adios_parse_vars_header_v1 (&md->b, &vars_header);

            for (i = 0; i < vars_header.count; i++)
            {
                memset (&var_payload, 0
                       ,sizeof (struct adios_var_payload_struct_v1)
                       );
                adios_parse_var_data_header_v1 (&md->b, &var_header);

                struct adios_var_struct * v1 = v;
                while (v1)
                {
                    if (   strcasecmp (var_header.name, v1->name)
                        || strcasecmp (var_header.path, v1->path)
                       )
                    {
                        v1 = v1->next;
                    }
                    else
                    {
                        break;
                    }
                }

                if (v1)
                {
                    var_payload.payload = v1->data;
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,&var_payload
                                                    ,v1->data_size
                                                    );
                }
                else
                {
                    printf ("MPI read: skipping name: %s path: %s\n"
                           ,var_header.name, var_header.path
                           );
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,NULL, 0
                                                    );
                }

                adios_clear_var_header_v1 (&var_header);
            }

#if 1
            adios_parse_attributes_header_v1 (&md->b, &attrs_header);

            for (i = 0; i < attrs_header.count; i++)
            {
                adios_parse_attribute_v1 (&md->b, &attribute);
                adios_clear_attribute_v1 (&attribute);
            }
#endif
            adios_clear_process_group_header_v1 (&pg_header);

            break;
        }

        default:
            fprintf (stderr, "MPI read: file version unknown: %u\n"
                    ,md->b.version
                    );
            return;
    }

    adios_buffer_struct_clear (&md->b);
}

void adios_adaptive_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_adaptive_data_struct * md =
                (struct adios_adaptive_data_struct *) method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;
#if COLLECT_METRICS
    gettimeofday (&timing.t23, NULL);
    static int iteration = 0;
#endif

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            // read the index to find the place to start reading
            adios_mpi_do_read (fd, method);
            struct adios_var_struct * v = fd->group->vars;
            while (v)
            {
                v->data = 0;
                v = v->next;
            }

            break;
        }

        case adios_mode_write:
        {
            MPI_Request req;
            MPI_Status status;
            int count;
            uint64_t msg [PARAMETER_COUNT];
            int err;
            int message_available = 0;
            struct thread_struct t;

            t.method = method;
            t.fd = fd;
            t.md = md;

            setup_threads (&t);

            if (md->rank == md->coord_rank)
            {
                struct flag_struct * f = malloc (sizeof (struct flag_struct));
                f->source = md->rank;
                f->tag = TAG_COORDINATOR;
                f->msg = malloc (8 * PARAMETER_COUNT);
                ((uint64_t *) f->msg) [0] = START_WRITES;
                pthread_mutex_lock (&md->coordinator_mutex);
                queue_enqueue (&md->coordinator_flag, f);
                pthread_mutex_unlock (&md->coordinator_mutex);
                pthread_cond_signal (&md->coordinator_cond);
            }
#if COLLECT_METRICS
            gettimeofday (&timing.t19, NULL);
#endif

            struct timespec ts;
            struct timeval tp;
            struct timeval interval;
            interval.tv_sec = 0;
            interval.tv_usec = 500000; // half a second wait
            int shutdown_flag = 0;
#if 0
char xx_b [10];
sprintf (xx_b, "%d.mpi", md->rank);
FILE * fp_mpi = fopen (xx_b, "w");
#endif
            do
            {
                int msg;
                pthread_mutex_lock (&md->mpi_mutex);
                gettimeofday (&tp, NULL);
                timeval_add (&tp, &tp, &interval);
                ts.tv_sec = tp.tv_sec;
                ts.tv_nsec = tp.tv_usec * 1000;

#if PRINT_MESSAGES
                printf ("cond wait MAIN %d\n", md->rank);
#endif
                pthread_cond_timedwait (&md->mpi_cond, &md->mpi_mutex, &ts);
#if PRINT_MESSAGES
                printf ("signaled MAIN %d\n", md->rank);
#endif
                MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, md->group_comm
                           ,&msg, &status
                           );
                pthread_mutex_unlock (&md->mpi_mutex);
                if (msg)
                {
#if 0
fprintf (fp_mpi, "get lock to receive\n");
fflush (fp_mpi);
#endif
                    pthread_mutex_lock (&md->mpi_mutex);
                    struct flag_struct * fs = malloc
                                               (sizeof (struct flag_struct));
                    MPI_Get_count (&status, MPI_BYTE, &count);
                    fs->source = status.MPI_SOURCE;
                    fs->tag = status.MPI_TAG;
                    fs->msg = malloc (count);
#if COLLECT_METRICS
                    struct timeval a, b, diff;
                    gettimeofday (&a, NULL);
#endif
                    MPI_Recv (fs->msg, count, MPI_BYTE, fs->source, fs->tag
                             ,md->group_comm, &status
                             );
#if 0
fprintf (fp_mpi, "received from %d %s\n", fs->source, message_to_string (((uint64_t *) fs->msg) [0]));
fflush (fp_mpi);
#endif
                    pthread_mutex_unlock (&md->mpi_mutex);
#if COLLECT_METRICS
                    gettimeofday (&b, NULL);
                    timeval_subtract (&diff, &b, &a);
                    timeval_add (&timing.t27, &timing.t27, &diff);
                    timing.recv_count++;
#endif

                    switch (fs->tag)
                    {
                        case TAG_WRITER:
                            pthread_mutex_lock (&md->writer_mutex);
                            queue_enqueue (&md->writer_flag, fs);
                            pthread_mutex_unlock (&md->writer_mutex);
                            pthread_cond_signal (&md->writer_cond);
                            break;

                        case TAG_SUB_COORDINATOR:
                        case TAG_SUB_COORDINATOR_INDEX_BODY:
                            pthread_mutex_lock (&md->sub_coordinator_mutex);
                            queue_enqueue (&md->sub_coordinator_flag, fs);
                            pthread_mutex_unlock (&md->sub_coordinator_mutex);
                            pthread_cond_signal (&md->sub_coordinator_cond);
                            break;

                        case TAG_COORDINATOR:
                        case TAG_COORDINATOR_INDEX_BODY:
                            pthread_mutex_lock (&md->coordinator_mutex);
                            queue_enqueue (&md->coordinator_flag, fs);
                            pthread_mutex_unlock (&md->coordinator_mutex);
                            pthread_cond_signal (&md->coordinator_cond);
                            break;

                        default:
                            printf ("Invalid tag: %d source: %d message "
                                    "size: %d\n"
                                   ,fs->tag, fs->source, count
                                   );
                            free (fs);
                            break;
                    }
                }

                while (queue_size (&md->mpi_flag) != 0)
                {
#if COLLECT_METRICS
                    struct timeval a, b, diff;
                    gettimeofday (&a, NULL);
#endif
                    pthread_mutex_lock (&md->mpi_mutex);
                    struct mpi_flag_struct * f;
                    queue_dequeue (&md->mpi_flag, &f);
#if 0
fprintf (fp_mpi, "sending msg to %d %s\n", f->target, message_to_string_full (((uint64_t *) f->msg), 'w'));
fflush (fp_mpi);
#endif
                    MPI_Send (f->msg, f->size, MPI_BYTE, f->target, f->tag
                             ,md->group_comm
                             );
                    pthread_mutex_unlock (&md->mpi_mutex);
#if COLLECT_METRICS
                    gettimeofday (&b, NULL);
                    timeval_subtract (&diff, &b, &a);
                    timing.send_count++;
                    timeval_add (&timing.t26, &timing.t26, &diff);

                    if (timing.mpi_queue_count == timing.mpi_queue_size)
                    {
                        timing.mpi_queue_size += 10;
                        timing.t29 = realloc (timing.t29
                                             ,sizeof (struct timeval_writer)
                                              * timing.mpi_queue_size
                                             );
                        assert (timing.t29 != 0);
                    }
                    timing.t29 [timing.mpi_queue_count].pid = f->target;
                    gettimeofday (&(timing.t29 [timing.mpi_queue_count++].t)
                                 ,NULL
                                 );
#endif
                    free (f->msg);
                    free (f);
                }

                if (md->rank == md->coord_rank)
                {
                    if (   md->coord_shutdown_flag == SHUTDOWN_FLAG
                        && md->sub_coord_shutdown_flag == SHUTDOWN_FLAG
                        && md->writer_shutdown_flag == SHUTDOWN_FLAG
                       )
                        shutdown_flag = SHUTDOWN_FLAG;
#if PRINT_MESSAGES
                    if (shutdown_flag == SHUTDOWN_FLAG)
                        printf ("xx: %2d coord finished\n", md->rank);
#endif
                }
                else
                {
                    if (md->rank == md->sub_coord_rank)
                    {
                        if (   md->sub_coord_shutdown_flag == SHUTDOWN_FLAG
                            && md->writer_shutdown_flag == SHUTDOWN_FLAG
                           )
                            shutdown_flag = SHUTDOWN_FLAG;
#if PRINT_MESSAGES
                        if (shutdown_flag == SHUTDOWN_FLAG)
                            printf ("xx: %2d(%d) sub coord finished\n"
                                   ,md->rank, md->group);
#endif
                    }
                    else
                    {
                        if (md->writer_shutdown_flag == SHUTDOWN_FLAG)
                        {
                            shutdown_flag = SHUTDOWN_FLAG;
#if PRINT_MESSAGES
                            if (shutdown_flag == SHUTDOWN_FLAG)
                                printf ("xx: %2d writer finished\n", md->rank);
#endif
                        }
                    }
                }
#if COLLECT_METRICS
                static int got_shutdown = 0;
                if (   !got_shutdown
                    && (   md->coord_shutdown_flag == SHUTDOWN_FLAG
                        || md->sub_coord_shutdown_flag == SHUTDOWN_FLAG
                        || md->writer_shutdown_flag == SHUTDOWN_FLAG
                       )
                   )
                {
                    gettimeofday (&timing.t20, NULL);
                    got_shutdown = 1;
                }
#endif
            } while (shutdown_flag != SHUTDOWN_FLAG);
#if 0
fclose (fp_mpi);
#endif

            md->probe_shutdown_flag = SHUTDOWN_FLAG;
            // finished with output. Shutdown the system
#if COLLECT_METRICS
            gettimeofday (&timing.t14, NULL);
#endif
            //err = pthread_cancel (md->probe);
            //err = pthread_kill (md->probe, SIGEND);
            err = pthread_join (md->probe, NULL);
            if (err != 0) printf ("thread probe kill error: %d\n", errno);
            err = pthread_join (md->writer, NULL);
            if (err != 0) printf ("thread join writer error: %d\n", err);

            if (md->rank == md->sub_coord_rank)
            {
                err = pthread_join (md->sub_coordinator, NULL);
                if (err != 0) printf ("thread join sub coord error: %d\n", err);
            }
            if (md->rank == md->coord_rank)
            {
                err = pthread_join (md->coordinator, NULL);
                if (err != 0) printf ("thread join sub coord error: %d\n", err);
            }
#if COLLECT_METRICS
            gettimeofday (&timing.t11, NULL);
#endif

            break;
        }

        case adios_mode_append:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = md->b.pg_index_offset;

            if (fd->shared_buffer == adios_flag_no)
            {
                MPI_Offset new_off;
                // set it up so that it will start at 0, but have correct sizes
                MPI_File_get_position (md->fh, &new_off);
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_vars_v1 (fd);
                // fd->vars_start gets updated with the size written
                MPI_File_seek (md->fh, md->vars_start, MPI_SEEK_SET);
                MPI_File_write (md->fh, fd->buffer, md->vars_header_size
                               ,MPI_BYTE, &md->status
                               );
                int count;
                MPI_Get_count (&md->status, MPI_BYTE, &count);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "d:MPI method tried to write %llu, "
                                     "only wrote %d\n"
                            ,md->vars_header_size
                            ,count
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
                adios_shared_buffer_free (&md->b);

                adios_write_open_attributes_v1 (fd);
                md->vars_start = new_off;
                md->vars_header_size = fd->offset;
                MPI_File_seek (md->fh, new_off + md->vars_header_size
                              ,MPI_SEEK_SET
                              ); // go back to end, but after attr header
                fd->base_offset += fd->offset;  // add size of header
                fd->offset = 0;
                fd->bytes_written = 0;

                while (a)
                {
                    adios_write_attribute_v1 (fd, a);
                    MPI_File_write (md->fh, fd->buffer, fd->bytes_written
                                   ,MPI_BYTE, &md->status
                                   );
                    MPI_Get_count (&md->status, MPI_BYTE, &count);
                    if (count != fd->bytes_written)
                    {
                        fprintf (stderr, "e:MPI method tried to write %llu, "
                                         "only wrote %d\n"
                                ,fd->bytes_written
                                ,count
                                );
                    }
                    fd->base_offset += count;
                    fd->offset = 0;
                    fd->bytes_written = 0;
                    adios_shared_buffer_free (&md->b);

                    a = a->next;
                }

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);
                MPI_File_seek (md->fh, md->vars_start, MPI_SEEK_SET);
                // fd->vars_start gets updated with the size written
                MPI_File_write (md->fh, fd->buffer, md->vars_header_size
                               ,MPI_BYTE, &md->status
                               );
                MPI_Get_count (&md->status, MPI_BYTE, &count);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "f:MPI method tried to write %llu, "
                                     "only wrote %d\n"
                            ,md->vars_header_size
                            ,count
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
            }

            // build index appending to any existing index
            adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root
                                 ,&md->old_attrs_root
                                 );
            // if collective, gather the indexes from the rest and call
            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    int * index_sizes = malloc (4 * md->size);
                    int * index_offsets = malloc (4 * md->size);
                    char * recv_buffer = 0;
                    uint32_t size = 0;
                    uint32_t total_size = 0;
                    int i;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, md->group_comm
                               );

                    for (i = 0; i < md->size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size + 1);

                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, md->group_comm
                                );

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    for (i = 1; i < md->size; i++)
                    {
                        md->b.buff = recv_buffer + index_offsets [i];
                        md->b.length = index_sizes [i];
                        md->b.offset = 0;

                        adios_parse_process_group_index_v1 (&md->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );
                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }
                    md->b.buff = buffer_save;
                    md->b.length = buffer_size_save;
                    md->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);
                }
                else
                {
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
                                         );

                    MPI_Gather (&buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->group_comm
                               );
                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->group_comm
                                );
                }
            }

            if (fd->shared_buffer == adios_flag_yes)
            {
                // everyone writes their data
                MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
                MPI_File_write (md->fh, fd->buffer, fd->bytes_written, MPI_BYTE
                               ,&md->status
                               );
            }

            if (md->rank == 0)
            {
                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,index_start, md->old_pg_root
                                     ,md->old_vars_root
                                     ,md->old_attrs_root
                                     );
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                MPI_File_seek (md->fh, md->b.pg_index_offset, MPI_SEEK_SET);
                MPI_File_write (md->fh, buffer, buffer_offset, MPI_BYTE
                               ,&md->status
                               );
            }

            free (buffer);

            adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);
            adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                                 ,md->old_attrs_root
                                 );
            new_pg_root = 0;
            new_vars_root = 0;
            new_attrs_root = 0;
            md->old_pg_root = 0;
            md->old_vars_root = 0;
            md->old_attrs_root = 0;

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);
        }
    }

    if (md && md->f)
    {
        close (md->f);
        md->f = -1;
    }

    if (md && md->fh)
        MPI_File_close (&md->fh);

    adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                         ,md->old_attrs_root
                         );
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;

#if COLLECT_METRICS
    gettimeofday (&timing.t28, NULL);
    print_metrics (md, iteration++);
#endif

    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
    {
        md->group_comm = MPI_COMM_NULL;
    }

    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->group_comm = MPI_COMM_NULL;
}

void adios_adaptive_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_adaptive_data_struct * md =
                  (struct adios_adaptive_data_struct *) method->method_data;

    if (adios_adaptive_initialized)
    {
        adios_adaptive_initialized = 0;

        queue_destroy (&md->mpi_flag);
        queue_destroy (&md->sub_coordinator_flag);
        queue_destroy (&md->coordinator_flag);
        queue_destroy (&md->writer_flag);
        pthread_mutex_destroy (&md->mpi_mutex);
        pthread_mutex_destroy (&md->writer_mutex);
        pthread_mutex_destroy (&md->sub_coordinator_mutex);
        pthread_mutex_destroy (&md->coordinator_mutex);
        pthread_cond_destroy (&md->mpi_cond);
        pthread_cond_destroy (&md->writer_cond);
        pthread_cond_destroy (&md->sub_coordinator_cond);
        pthread_cond_destroy (&md->coordinator_cond);
    }
}

void adios_adaptive_end_iteration (struct adios_method_struct * method)
{
}

void adios_adaptive_start_calculation (struct adios_method_struct * method)
{
}

void adios_adaptive_stop_calculation (struct adios_method_struct * method)
{
}

// each sub coordinator will do the following:
// 1. wait for a registration indication from all writers it manages
// 2. send each writer, one at a time, a note to write and wait for
//    a completion notice.
// 3. Possibly process coordinator messages to tell a writer to go to a
//    different file. Wait for response for writer and then respond to
//    coordinator
// 4. When done tell coordinator and wait for shutdown message
// 5. shutdown
static void * sub_coordinator_main (void * param)
{
    struct thread_struct * t = (struct thread_struct *) param;

    struct adios_adaptive_data_struct * md = t->md;

    MPI_Request req;
    MPI_Status status;
    int message_available;
    struct flag_struct * f = 0;
    uint64_t * msg = 0;   // shortcut for the message payload

    ////////////////////////////////////////////////////////////////////////
    // wait for START_WRITES
    ////////////////////////////////////////////////////////////////////////
    //while (queue_size (&md->sub_coordinator_flag) == 0)
    //    ;
    pthread_mutex_lock (&md->sub_coordinator_mutex);
#if PRINT_MESSAGES
    printf ("cond wait SUB COORDINATOR %d START_WRITES\n", md->group);
#endif
    pthread_cond_wait (&md->sub_coordinator_cond, &md->sub_coordinator_mutex);
#if PRINT_MESSAGES
    printf ("signaled SUB_COORDINATOR %d START_WRITES\n", md->group);
#endif
    queue_dequeue (&md->sub_coordinator_flag, &f);
    pthread_mutex_unlock (&md->sub_coordinator_mutex);
    msg = (uint64_t *) f->msg;
    assert (msg [0] == START_WRITES);
    free (f->msg);
    free (f);
    f = 0;
    msg = 0;

    ////////////////////////////////////////////////////////////////////////
    // do writing
    ////////////////////////////////////////////////////////////////////////
    int completed_writing = 0;  // track when we have notified coordinator we
                                // have completed (avoid multiple notifies)
    int writers_size = (int) (md->group_size * 1.20); // add 20% for adaptation
    int * writers = (int *) malloc (writers_size * sizeof (int));
    // start at first so that we can always write using next_writer
    int current_writer = md->sub_coord_rank - md->group_size + 1;
    int next_writer = current_writer + 1; // track the next one so the adaptive
                                      // writer knows who to start next
    int writers_served = 0;
    int indices_collected = 0;
    int active_writers = 0; // track how many of our procs are writing

    int * adaptive_writers = 0;       // keep track of adaptive ranks pending
    int adaptive_writers_size = 0;    // track size of array
    int adaptive_writers_being_served = 0;  // how many pending

    uint64_t current_offset = 0;

    int * index_sizes = (int *) malloc (sizeof (int) * md->size);
    int largest_index = 0;

    // for encoding the local index for writing and sending to the coord
    char * buffer = 0;
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;

    message_available = 0;
    uint64_t shutdown_flag = NO_FLAG;
#if PRINT_MESSAGES
    int xxx = 0;
#else
    int xxx = 0;
#endif
#if 0
    char buf_x [10];
    sprintf (buf_x, "%d.w", md->group);
    FILE * fp_w = fopen (buf_x, "w");
#endif
    int currently_writing = 0;
    int local_writer = 0;
    enum MESSAGE_FLAGS last_message = NO_FLAG;
    do
    {
        xxx++;
#if PRINT_MESSAGES
        //if (!(xxx++ % 10000000))
            printf ("AAAA %d\n", md->group);
#endif
        if (queue_size (&md->sub_coordinator_flag) != 0)
        {
            message_available = 1;
            pthread_mutex_lock (&md->sub_coordinator_mutex);
            queue_dequeue (&md->sub_coordinator_flag, &f);
            pthread_mutex_unlock (&md->sub_coordinator_mutex);
            msg = (uint64_t *) f->msg;
        }

        if (message_available)
        {
            last_message = (((uint64_t *) f->msg) [0]);
#if PRINT_MESSAGES
            printf ("sc: %d source: %2d %s A\n", md->group, f->source
                   ,message_to_string_full (msg, 's')
                   );
#endif
            message_available = 0;
            switch (((uint64_t *) f->msg) [0])
            {
                case WRITE_COMPLETE:
                {
#if 0
fprintf (fp_w, "%d\treceived %s\n", f->source, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
#if COLLECT_METRICS
                    if (   timing.write_complete_count
                        == timing.write_complete_size
                       )
                    {
                        timing.write_complete_size += 10;
                        timing.t17 = realloc (timing.t17
                                             ,sizeof (struct timeval_writer)
                                              * timing.write_complete_size
                                             );
                        assert (timing.t17 != 0);
                    }
                    timing.t17 [timing.write_complete_count].pid = msg [5];
                    gettimeofday (&(timing.t17
                                      [timing.write_complete_count++].t), NULL);
#endif
                    // if this group was writing it, we were tracking it
                    if (msg [2] == md->group)
                    {
                        active_writers--;
                    }

                    // if it was a purely local write
                    local_writer--;

                    // if the target group is our group (we wrote to our file)
                    if (msg [1] == md->group)
                    {
#if 0
fprintf (fp_w, "%d\tsave for index %s\n", msg [5], message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                        // save the writer in our list for indexing
                        if (writers_size <= writers_served)
                        {
                            writers_size = (int) (writers_size * 1.20);
                            writers = (int *) realloc (writers,   writers_size
                                                                * sizeof (int)
                                                      );
                        }
                        writers [writers_served] = msg [5];
                        index_sizes [writers_served] = msg [4];
                        if (index_sizes [writers_served] > largest_index)
                            largest_index = index_sizes [writers_served];

                        writers_served++;
                    }

                    // if it is adaptive, tell the coordinator it is done
                    if (msg [2] == md->group && msg [1] != msg [2])
                    {
#if 0
fprintf (fp_w, "%d\ttell coord adaptive done %s\n", md->rank, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                        // tell coordinator done with this one
                        // and if we have more capacity
                        // only if we were writing to our group
                        //if (msg [1] != msg [2])
                        {
                            if (md->rank != md->coord_rank)
                            {
                                struct mpi_flag_struct * f = malloc
                                           (sizeof (struct mpi_flag_struct));
                                f->target = md->coord_rank;
                                f->tag = TAG_COORDINATOR;
                                f->size = 8 * PARAMETER_COUNT;
                                f->msg = malloc (8 * PARAMETER_COUNT);
                                ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                                ((uint64_t *) f->msg) [1] = msg [1];
                                ((uint64_t *) f->msg) [2] = msg [2];
                                ((uint64_t *) f->msg) [3] = msg [3];
                                ((uint64_t *) f->msg) [4] = msg [4];
                                ((uint64_t *) f->msg) [5] = msg [5];
                                pthread_mutex_lock (&md->mpi_mutex);
                                queue_enqueue (&md->mpi_flag, f);
                                pthread_mutex_unlock (&md->mpi_mutex);
                                pthread_cond_signal (&md->mpi_cond);
                            }
                            else
                            {
                                struct flag_struct * f = malloc
                                               (sizeof (struct flag_struct));
                                f->source = md->rank;
                                f->tag = TAG_COORDINATOR;
                                f->msg = malloc (8 * PARAMETER_COUNT);
                                ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                                ((uint64_t *) f->msg) [1] = msg [1];
                                ((uint64_t *) f->msg) [2] = msg [2];
                                ((uint64_t *) f->msg) [3] = msg [3];
                                ((uint64_t *) f->msg) [4] = msg [4];
                                ((uint64_t *) f->msg) [5] = msg [5];
                                pthread_mutex_lock (&md->coordinator_mutex);
                                queue_enqueue (&md->coordinator_flag, f);
                                pthread_mutex_unlock (&md->coordinator_mutex);
                                pthread_cond_signal (&md->coordinator_cond);
                            }
                        }
                    }

                    // if it was local, start the next write or send
                    // the write complete for this group (at bottom).
                    // if we finished others and we were adaptive and there
                    // isn't a purely local writer currently in process,
                    // fall through that code as well to end our write.
                    if (   msg [1] == msg [2]
                        || (   next_writer > md->rank
                            && !local_writer
                           )
                       )
                    {
                        currently_writing = 0;
                        current_writer = next_writer++;
                    }
                    break;
                }

                case ADAPTIVE_WRITE_START:
                {
#if 0
fprintf (fp_w, "%d\treceived %s\n", f->source, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                    if (next_writer > md->rank)
                    {
                        // tell coordinator we are done and can't do it
                        if (md->rank != md->coord_rank)
                        {
                            struct mpi_flag_struct * f = malloc
                                           (sizeof (struct mpi_flag_struct));
                            f->target = md->coord_rank;
                            f->tag = TAG_COORDINATOR;
                            f->size = 8 * PARAMETER_COUNT;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = WRITERS_BUSY;
                            ((uint64_t *) f->msg) [1] = msg [1];
                            ((uint64_t *) f->msg) [2] = md->group;
                            pthread_mutex_lock (&md->mpi_mutex);
                            queue_enqueue (&md->mpi_flag, f);
                            pthread_mutex_unlock (&md->mpi_mutex);
                            pthread_cond_signal (&md->mpi_cond);
                        }
                        else
                        {
                            struct flag_struct * f = malloc
                                                (sizeof (struct flag_struct));
                            f->source = md->rank;
                            f->tag = TAG_COORDINATOR;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = WRITERS_BUSY;
                            ((uint64_t *) f->msg) [1] = msg [1];
                            ((uint64_t *) f->msg) [2] = md->group;
                            pthread_mutex_lock (&md->coordinator_mutex);
                            queue_enqueue (&md->coordinator_flag, f);
                            pthread_mutex_unlock (&md->coordinator_mutex);
                            pthread_cond_signal (&md->coordinator_cond);
                        }
                    }
                    else
                    {
                        if (next_writer < md->rank)
                        {
#if COLLECT_METRICS
                            if (timing.do_write_count == timing.do_write_size)
                            {
                                timing.do_write_size += 10;
                                timing.t15 = realloc (timing.t15
                                               ,sizeof (struct timeval_writer)
                                                * timing.do_write_size
                                               );
                                assert (timing.t15 != 0);
                            }
                            timing.t15 [timing.do_write_count].pid
                                                                 = next_writer;
                            gettimeofday (&(timing.t15
                                           [timing.do_write_count++].t), NULL);
#endif
                            struct mpi_flag_struct * f = malloc
                                           (sizeof (struct mpi_flag_struct));
                            f->target = next_writer;
                            f->tag = TAG_WRITER;
                            f->size = 8 * PARAMETER_COUNT;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = DO_WRITE_FLAG;
                            ((uint64_t *) f->msg) [1] = msg [1];
                            ((uint64_t *) f->msg) [2] = msg [2];
                            ((uint64_t *) f->msg) [3] = md->group;
                            ((uint64_t *) f->msg) [4] = md->rank;
                            ((uint64_t *) f->msg) [5] = msg [3];
#if 0
fprintf (fp_w, "%d\tsending %s\n", f->target, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                            pthread_mutex_lock (&md->mpi_mutex);
                            queue_enqueue (&md->mpi_flag, f);
                            pthread_mutex_unlock (&md->mpi_mutex);
                            pthread_cond_signal (&md->mpi_cond);
                            active_writers++;
                            next_writer++;
                        }
                        else
                        {
                            if (next_writer == md->rank)
                            {
#if COLLECT_METRICS
                                if (timing.do_write_count
                                                      == timing.do_write_size)
                                {
                                    timing.do_write_size += 10;
                                    timing.t15 = realloc (timing.t15
                                              ,sizeof (struct timeval_writer)
                                                * timing.do_write_size
                                              );
                                    assert (timing.t15 != 0);
                                }
                                timing.t15 [timing.do_write_count].pid
                                                                = next_writer;
                                gettimeofday (&(timing.t15
                                           [timing.do_write_count++].t), NULL);
#endif
                                active_writers++;
                                struct flag_struct * f = malloc
                                                (sizeof (struct flag_struct));
                                f->source = md->rank;
                                f->tag = TAG_WRITER;
                                f->msg = malloc (8 * PARAMETER_COUNT);
                                ((uint64_t *) f->msg) [0] = DO_WRITE_FLAG;
                                ((uint64_t *) f->msg) [1] = msg [1];
                                ((uint64_t *) f->msg) [2] = msg [2];
                                ((uint64_t *) f->msg) [3] = md->group;
                                ((uint64_t *) f->msg) [4] = md->rank;
                                ((uint64_t *) f->msg) [5] = msg [3];
#if 0
fprintf (fp_w, "%d\tsending %s\n", md->rank, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                                pthread_mutex_lock (&md->writer_mutex);
                                queue_enqueue (&md->writer_flag, f);
                                pthread_mutex_unlock (&md->writer_mutex);
                                pthread_cond_signal (&md->writer_cond);
                                next_writer++;
                            }
                        }
                    }
                    break;
                }

                case OVERALL_WRITE_COMPLETE:
                {
#if COLLECT_METRICS
                    gettimeofday (&timing.t12, NULL);
#endif
                    // tell writers to enter send index mode
                    int i;
                    char remote_writer = 0;
                    char local_writer = 0;
#if 0
char b1 [10];
sprintf (b1, "%d", md->group);
FILE * ff = fopen (b1, "w");
#endif
                    for (i = 0; i < writers_served; i++)
                    {
#if 0
fprintf (ff, "%d\n", writers [i]);
#endif
                        if (writers [i] != md->rank)
                        {
                            struct mpi_flag_struct * f = malloc
                                            (sizeof (struct mpi_flag_struct));
                            f->target = writers [i];
                            f->tag = TAG_WRITER;
                            f->size = 8 * PARAMETER_COUNT;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = SEND_INDEX;
                            pthread_mutex_lock (&md->mpi_mutex);
                            queue_enqueue (&md->mpi_flag, f);
                            pthread_mutex_unlock (&md->mpi_mutex);
                            remote_writer = 1;
                        }
                        else
                        {
                            struct flag_struct * f = malloc
                                                (sizeof (struct flag_struct));
                            f->source = md->rank;
                            f->tag = TAG_WRITER;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = SEND_INDEX;
                            pthread_mutex_lock (&md->writer_mutex);
                            queue_enqueue (&md->writer_flag, f);
                            pthread_mutex_unlock (&md->writer_mutex);
                            local_writer = 1;
                        }
                    }
#if 0
fclose (ff);
#endif
                    if (remote_writer)
                        pthread_cond_signal (&md->mpi_cond);
                    if (local_writer)
                        pthread_cond_signal (&md->writer_cond);

#if DO_INDEX_COLLECTION == 0
#if COLLECT_METRICS
                    gettimeofday (&timing.t13, NULL);
#endif
                    // once we are done sending the index, shutdown the thread
#if 0
fprintf (fp_w, "%d\tshutdown 2\n", md->rank);
fflush (fp_w);
#endif
                    shutdown_flag = SHUTDOWN_FLAG;
#endif

                    break;
                }

                case INDEX_BODY:
                {
                    indices_collected++;
#if DO_INDEX_COLLECTION
                    char * buf = malloc (largest_index + 1);
                    buf [largest_index] = 0;
                    struct adios_bp_buffer_struct_v1 b;
                    struct adios_index_process_group_struct_v1 * new_pg_root;
                    struct adios_index_var_struct_v1 * new_vars_root;
                    struct adios_index_attribute_struct_v1 * new_attrs_root;
                    new_pg_root = 0;
                    new_vars_root = 0;
                    new_attrs_root = 0;
                    for (i = 0; i < writers_served; i++)
                    {
                        buf [index_sizes [i]] = 0;
                        if (writers [i] != md->rank)
                        {
                            int message_available = 0;
                            while (!message_available)
                            {
#if COLLECT_METRICS
                                struct timeval a, b, diff;
                                gettimeofday (&a, NULL);
#endif
                                pthread_mutex_lock (&md->mpi_mutex);
                                MPI_Iprobe (MPI_ANY_SOURCE, TAG_SUB_COORDINATOR
                                           ,md->group_comm
                                           ,&message_available, &status
                                           );
                                pthread_mutex_unlock (&md->mpi_mutex);
#if COLLECT_METRICS
                                gettimeofday (&b, NULL);
                                iprobe_count++;
                                timeval_subtract (&diff, &b, &a);
                                timeval_add (&timing.t25, &timing.t25, &diff);
#endif
                            }

#if COLLECT_METRICS
                            struct timeval a, b, diff;
                            gettimeofday (&a, NULL);
#endif
                            pthread_mutex_lock (&md->mpi_mutex);
                            MPI_Recv (buf, index_sizes [i], MPI_BYTE
                                      ,writers [i]
                                      ,TAG_SUB_COORDINATOR, md->group_comm
                                      ,&status
                                      );

                            pthread_mutex_unlock (&md->mpi_mutex);
#if COLLECT_METRICS
                            gettimeofday (&b, NULL);
                            timeval_subtract (&diff, &b, &a);
                            timeval_add (&timing.t27, &timing.t27, &diff);
                            timing.recv_count++;
#endif
                            b.buff = buf;
                        }
                        else
                        {
                            uint64_t * flag;
                            while (queue_size (&md->sub_coordinator_flag) == 0)
                                ;
                            pthread_mutex_lock (&md->sub_coordinator_mutex);
                            queue_dequeue (&md->sub_coordinator_flag, &flag);
                            pthread_mutex_unlock (&md->sub_coordinator_mutex);
                            b.buff = (char *) (flag [1]);
                            //free (flag); // freed on writer side
                        }

                        // merge buf into the index
                        b.length = index_sizes [i];
                        b.offset = 0;

                        adios_parse_process_group_index_v1 (&b ,&new_pg_root);
                        adios_parse_vars_index_v1 (&b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&b ,&new_attrs_root);
                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }
                    free (buf);

                    uint64_t only_index_buffer_offset;
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
                                         );
                    only_index_buffer_offset = buffer_offset;
                    adios_write_version_v1 (&buffer, &buffer_size
                                           ,&buffer_offset
                                           );

                    lseek (md->f, md->stripe_size * msg [1], SEEK_SET);
                    ssize_t s = write (md->f, buffer, buffer_offset);
                    if (s != buffer_offset)
                    {
                        fprintf (stderr, "Need to do multi-write 3\n");
                    }

                    uint64_t index_size = buffer_offset;

                    // send index to the coordinator for global use
                    if (md->rank != md->coord_rank)
                    {
                        struct mpi_flag_struct * f = malloc
                                           (sizeof (struct mpi_flag_struct));
                        f->target = md->coord_rank;
                        f->tag = TAG_COORDINATOR;
                        f->size = 8 * PARAMETER_COUNT;
                        f->msg = malloc (8 * PARAMETER_COUNT);
                        ((uint64_t *) f->msg) [0] = INDEX_SIZE;
                        ((uint64_t *) f->msg) [1] = md->group;
                        ((uint64_t *) f->msg) [2] = only_index_buffer_offset;
                        pthread_mutex_lock (&md->mpi_mutex);
                        queue_enqueue (&md->mpi_flag, f);
                        pthread_mutex_unlock (&md->mpi_mutex);
                        pthread_cond_signal (&md->mpi_cond);

                        pthread_mutex_lock (&md->mpi_mutex);
                        MPI_Isend (buffer, only_index_buffer_offset, MPI_BYTE
                                  ,md->coord_rank, TAG_COORDINATOR
                                  ,md->group_comm, &req
                                  );

                        MPI_Wait (&req, &status);
                        pthread_mutex_unlock (&md->mpi_mutex);
                    }
                    else
                    {
printf ("FIX THIS\n");
                        uint64_t * flag = malloc (8 * PARAMETER_COUNT);
                        INIT_PARAMS(flag);
                        flag [0] = INDEX_SIZE;
                        flag [1] = md->group;
                        flag [2] = only_index_buffer_offset;
                        flag [3] = (uint64_t) buffer;
                        pthread_mutex_lock (&md->coordinator_mutex);
                        queue_enqueue (&md->coordinator_flag, flag);
                        pthread_mutex_unlock (&md->coordinator_mutex);
                        pthread_cond_signal (&md->coordinator_cond);
                    }
                    free (buffer);
#endif

#if 0
fprintf (fp_w, "%d\tshutdown 1\n", md->rank);
fflush (fp_w);
#endif
                    // once we are done sending the index, shutdown the thread
                    shutdown_flag = SHUTDOWN_FLAG;
#if DO_INDEX_COLLECTION
                    indices_collected = writers_served;
#endif

                    if (indices_collected == writers_served)
                    {
#if COLLECT_METRICS
                        gettimeofday (&timing.t13, NULL);
#endif
                    }
                    break;
                }

                case SHUTDOWN_FLAG:
                {
                    break;
                }
            }
#if PRINT_MESSAGES
            printf ("sc: %d source: %2d %s B\n", md->group, f->source
                   ,message_to_string_full (msg, 's')
                   );
#endif
        }
        if (!currently_writing)
        {
            if (current_writer < md->rank)
            {
#if COLLECT_METRICS
                if (timing.do_write_count == timing.do_write_size)
                {
                    timing.do_write_size += 10;
                    timing.t15 = realloc (timing.t15
                                         ,sizeof (struct timeval_writer)
                                          * timing.do_write_size
                                         );
                    assert (timing.t15 != 0);
                }
                timing.t15 [timing.do_write_count].pid = current_writer;
                gettimeofday (&(timing.t15 [timing.do_write_count++].t), NULL);
#endif
                local_writer++;
                struct mpi_flag_struct * f = malloc
                                           (sizeof (struct mpi_flag_struct));
                f->target = current_writer;
                f->tag = TAG_WRITER;
                f->size = 8 * PARAMETER_COUNT;
                f->msg = malloc (8 * PARAMETER_COUNT);
                ((uint64_t *) f->msg) [0] = DO_WRITE_FLAG;
                ((uint64_t *) f->msg) [1] = md->group;
                ((uint64_t *) f->msg) [2] = md->rank;
                ((uint64_t *) f->msg) [3] = md->group;
                ((uint64_t *) f->msg) [4] = md->rank;
                ((uint64_t *) f->msg) [5] = current_offset++;
#if 0
fprintf (fp_w, "%d\tsending %s\n", f->target, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                pthread_mutex_lock (&md->mpi_mutex);
                queue_enqueue (&md->mpi_flag, f);
                pthread_mutex_unlock (&md->mpi_mutex);
                pthread_cond_signal (&md->mpi_cond);
                active_writers++;
                currently_writing = 1;
            }
            else
            {
                if (current_writer == md->rank)
                {
#if COLLECT_METRICS
                    if (timing.do_write_count == timing.do_write_size)
                    {
                        timing.do_write_size += 10;
                        timing.t15 = realloc (timing.t15
                                             ,sizeof (struct timeval_writer)
                                              * timing.do_write_size
                                             );
                        assert (timing.t15 != 0);
                    }
                    timing.t15 [timing.do_write_count].pid = current_writer;
                    gettimeofday (&(timing.t15
                                          [timing.do_write_count++].t), NULL);
#endif
                    active_writers++;
                    local_writer++;
                    currently_writing = 1;
                    struct flag_struct * f = malloc
                                              (sizeof (struct flag_struct));
                    f->source = md->rank;
                    f->tag = TAG_WRITER;
                    f->msg = malloc (8 * PARAMETER_COUNT);
                    ((uint64_t *) f->msg) [0] = DO_WRITE_FLAG;
                    ((uint64_t *) f->msg) [1] = md->group;
                    ((uint64_t *) f->msg) [2] = md->rank;
                    ((uint64_t *) f->msg) [3] = md->group;
                    ((uint64_t *) f->msg) [4] = md->rank;
                    ((uint64_t *) f->msg) [5] = current_offset++;
#if 0
fprintf (fp_w, "%d\tsending %s\n", md->rank, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                    pthread_mutex_lock (&md->writer_mutex);
                    queue_enqueue (&md->writer_flag, f);
                    pthread_mutex_unlock (&md->writer_mutex);
                    pthread_cond_signal (&md->writer_cond);
                }
                else
                {
                    if (!completed_writing) // && !active_writers)
                    {
                        completed_writing = 1;
                        if (md->rank != md->coord_rank)
                        {
                            struct mpi_flag_struct * f = malloc
                                           (sizeof (struct mpi_flag_struct));
                            f->target = md->coord_rank;
                            f->tag = TAG_COORDINATOR;
                            f->size = 8 * PARAMETER_COUNT;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                            ((uint64_t *) f->msg) [1] = md->group;
                            ((uint64_t *) f->msg) [2] = md->group;
                            ((uint64_t *) f->msg) [3] = msg [3];
                            ((uint64_t *) f->msg) [4] = msg [4];
                            ((uint64_t *) f->msg) [5] = md->rank;
#if 0
fprintf (fp_w, "%d\tsending %s\n", md->rank, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                            pthread_mutex_lock (&md->mpi_mutex);
                            queue_enqueue (&md->mpi_flag, f);
                            pthread_mutex_unlock (&md->mpi_mutex);
                            pthread_cond_signal (&md->mpi_cond);
                        }
                        else
                        {
                            struct flag_struct * f = malloc
                                               (sizeof (struct flag_struct));
                            f->source = md->rank;
                            f->tag = TAG_COORDINATOR;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                            ((uint64_t *) f->msg) [1] = md->group; //msg [1];
                            ((uint64_t *) f->msg) [2] = md->group; //msg [2];
                            ((uint64_t *) f->msg) [3] = msg [3];
                            ((uint64_t *) f->msg) [4] = msg [4];
                            ((uint64_t *) f->msg) [5] = md->rank;
#if 0
fprintf (fp_w, "%d\tsending %s\n", md->rank, message_to_string_full ((uint64_t *) f->msg, 's'));
fflush (fp_w);
#endif
                            pthread_mutex_lock (&md->coordinator_mutex);
                            queue_enqueue (&md->coordinator_flag, f);
                            pthread_mutex_unlock (&md->coordinator_mutex);
                            pthread_cond_signal (&md->coordinator_cond);
                        }
                    }
                }
            }
        }

        if (f)
        {
            free (f->msg);
            free (f);
            f = 0;
            msg = 0;
        }
        if (shutdown_flag != SHUTDOWN_FLAG)
        {
            if (queue_size (&md->sub_coordinator_flag) == 0)
            {
                struct timespec ts;
                struct timeval tp;
                struct timeval interval;
                interval.tv_sec = 0;
                interval.tv_usec = 500000; // half a second wait

                pthread_mutex_lock (&md->sub_coordinator_mutex);
                gettimeofday (&tp, NULL);
                timeval_add (&tp, &tp, &interval);
                ts.tv_sec = tp.tv_sec;
                ts.tv_nsec = tp.tv_usec * 1000;

#if PRINT_MESSAGES
                printf ("cond wait SUB COORDINATOR %d\n", md->group);
#endif
                pthread_cond_timedwait (&md->sub_coordinator_cond
                                       ,&md->sub_coordinator_mutex, &ts
                                       );
#if PRINT_MESSAGES
                printf ("signaled SUB COORDINATOR %d\n", md->group);
#endif
                pthread_mutex_unlock (&md->sub_coordinator_mutex);
            }
//            if (xxx >= 300) {shutdown_flag = SHUTDOWN_FLAG; fprintf (stderr, "rank:\t%d\taborting sub coord timed out last message: %s\n", md->rank, message_to_string (last_message));}
        }
    } while (shutdown_flag != SHUTDOWN_FLAG);
#if 0
fclose (fp_w);
#endif
    free (writers);
    //free (index_sizes);

    while (   queue_size (&md->mpi_flag) != 0
           || queue_size (&md->coordinator_flag) != 0
          )
        ;
    md->sub_coord_shutdown_flag = SHUTDOWN_FLAG;
    pthread_cond_signal (&md->mpi_cond);
    pthread_exit (NULL);
    return NULL;
}

// one coordinator overall that only talks with the sub_coordinators
// 1. wait for registration from each of the sub_coordinators.
// 2. wait for messages from sub_coordinators indicating completion
// 3. Whenever a subcoordinator indicates it is done, remove it from the
//    pending set and add it to the adaptation set.
// 4. Go through the pending set telling one to tell it to use on of the
//    adaptation set. Wait for a response from that adaptive write.
// 5. When all sub_coordinators report completed, tell sub_coordinators to
//    build index, close, and send a copy of the index to the coordinator.
// 6. collect the sub indices and make the overall index file.
// 7. wait for shutdown message.
// 8. shutdown.
static void * coordinator_main (void * param)
{
    struct thread_struct * t = (struct thread_struct *) param;

    struct adios_adaptive_data_struct * md = t->md;

    struct flag_struct * f = 0;
    uint64_t  * msg = 0;
    int source;
    int message_available;
    MPI_Status status;
    MPI_Request req;
    int i;
    enum group_state
    {
         STATE_WRITING = 0
        ,STATE_COMPLETE = 1
        ,STATE_WRITERS_ALL_OCCUPIED = 2
    };

    int sub_coord_ranks [md->groups];

    char group_state [md->groups];
    int groups_complete = 0;
    uint64_t group_offset [md->groups];

    memset (group_offset, 0, md->groups * sizeof (uint64_t));
    memset (group_state, STATE_WRITING, md->groups);

    ////////////////////////////////////////////////////////////////////////
    // sub_coordinator rank setup
    ////////////////////////////////////////////////////////////////////////
    i = 0;
    int rank = 0;
    while (rank < md->size)
    {
        int group;
        int group_size;
        int sub_coord_rank;
        int coord_rank;

        calc_groups (rank, md->size, md->groups
                    ,&group, &group_size, &sub_coord_rank
                    ,&coord_rank
                    );
        sub_coord_ranks [i++] = sub_coord_rank;
        rank += group_size;
#if PRINT_MESSAGES
        printf ("calc sc: %d rank: %d\n", i - 1, sub_coord_ranks [i - 1]);
#endif
    }

    ////////////////////////////////////////////////////////////////////////
    // wait for START_WRITE to complete before continuing
    ////////////////////////////////////////////////////////////////////////
#if PRINT_MESSAGES
    printf ("wait for START_WRITES\n");
#endif
#if 1
    while (queue_size (&md->coordinator_flag) == 0)
        ;
    pthread_mutex_lock (&md->coordinator_mutex);
#else
    pthread_mutex_lock (&md->coordinator_mutex);
#if PRINT_MESSAGES
    printf ("cond wait COORDINATOR\n");
#endif
    pthread_cond_wait (&md->coordinator_cond, &md->coordinator_mutex);
#if PRINT_MESSAGES
    printf ("signaled COORDINATOR\n");
#endif
#endif
    queue_dequeue (&md->coordinator_flag, &f);
    pthread_mutex_unlock (&md->coordinator_mutex);
    msg = (uint64_t *) f->msg;
    assert (msg [0] == START_WRITES);
    free (f->msg);
    free (f);
#if COLLECT_METRICS
    gettimeofday (&timing.t7, NULL);
#endif

    for (i = 0; i < md->groups; i++)
    {
#if PRINT_MESSAGES
        printf ("START_WRITES to %d BEGIN\n", sub_coord_ranks [i]);
#endif
        if (sub_coord_ranks [i] != md->rank)
        {
            struct mpi_flag_struct * f = malloc
                                           (sizeof (struct mpi_flag_struct));
            f->target = sub_coord_ranks [i];
            f->tag = TAG_SUB_COORDINATOR;
            f->size = 8;
            f->msg = malloc (8);
            ((uint64_t *) f->msg) [0] = START_WRITES;
            pthread_mutex_lock (&md->mpi_mutex);
            queue_enqueue (&md->mpi_flag, f);
            pthread_mutex_unlock (&md->mpi_mutex);
        }
        else
        {
            f = malloc (sizeof (struct flag_struct));
            f->source = md->rank;
            f->tag = TAG_SUB_COORDINATOR;
            f->msg = malloc (8 * PARAMETER_COUNT);
            ((uint64_t *) f->msg) [0] = START_WRITES;
            pthread_mutex_lock (&md->sub_coordinator_mutex);
            queue_enqueue (&md->sub_coordinator_flag, f);
            pthread_mutex_unlock (&md->sub_coordinator_mutex);
            f = 0;
        }
        pthread_cond_signal (&md->mpi_cond);
        pthread_cond_signal (&md->sub_coordinator_cond);
#if PRINT_MESSAGES
        printf ("START_WRITES to %d END\n", sub_coord_ranks [i]);
#endif
    }

    ////////////////////////////////////////////////////////////////////////
    // Main message loop
    ////////////////////////////////////////////////////////////////////////
    uint64_t largest_index_size = 0;
    int index_sizes_received = 0;
    char * index_buf = 0;     // for receiving from remote
    ssize_t index_size = 0;

    char * buffer = 0;          // for building overall
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;

    struct adios_file_index_format_v2
    {
        uint32_t file_number;
        uint16_t name_len;
        char * name;
        uint64_t offset;
        uint64_t length;
    };

    uint64_t file_index_count = 0;
    char start_index_collection = 0;
    char index_collection_started = 0;
    uint64_t adaptive_writes_outstanding = 0;

    struct adios_file_index_format_v2 * file_index =
       (struct adios_file_index_format_v2 *)
              malloc (sizeof (struct adios_file_index_format_v2) * md->groups);

    message_available = 0;
    uint64_t shutdown_flag = NO_FLAG;
    struct timespec ts;
    struct timeval tp;
    struct timeval interval;
    interval.tv_sec = 0;
    interval.tv_usec = 500000; // half a second wait
#if PRINT_MESSAGES
    int xxx = 0;
#else
    int xxx = 0;
#endif
    enum MESSAGE_FLAGS last_message = NO_FLAG;
#if 0
FILE * fp_c = fopen ("c", "w");
#endif
    do
    {
        xxx++;
#if PRINT_MESSAGES
        //if (!(xxx++ % 10000000))
            printf ("BBBB adaptive_writes_outstanding: %lld\n"
                   ,adaptive_writes_outstanding
                   );
#endif
        if (queue_size (&md->coordinator_flag) == 0)
        {
            pthread_mutex_lock (&md->coordinator_mutex);
            gettimeofday (&tp, NULL);
            timeval_add (&tp, &tp, &interval);
            ts.tv_sec = tp.tv_sec;
            ts.tv_nsec = tp.tv_usec * 1000;

#if PRINT_MESSAGES
            printf ("cond wait COORDINATOR\n");
#endif
            pthread_cond_timedwait (&md->coordinator_cond
                                   ,&md->coordinator_mutex
                                   ,&ts
                                   );
#if PRINT_MESSAGES
            printf ("signaled COORDINATOR\n");
#endif
            pthread_mutex_unlock (&md->coordinator_mutex);
        }

        if (queue_size (&md->coordinator_flag) != 0)
        {
            message_available = 1;
            pthread_mutex_lock (&md->coordinator_mutex);
            queue_dequeue (&md->coordinator_flag, &f);
            pthread_mutex_unlock (&md->coordinator_mutex);
            msg = (uint64_t *) f->msg;
        }

        if (message_available)
        {
            last_message = msg [0];
#if PRINT_MESSAGES
            printf ("c: source: %2d msg: %s A\n"
                   ,f->source, message_to_string_full (msg, 'c')
                   );
#endif
            message_available = 0;
            switch (msg [0])
            {
                case WRITE_COMPLETE:
                {
#if 0
fprintf (fp_c, "adaptive_writes_outstanding: %lld msg: %s\n"
        ,adaptive_writes_outstanding, message_to_string_full (msg, 'c')
        );
fflush (fp_c);
#endif
                    // what part of the file finished writing last is unknown
                    // so only save the largest final offset
                    if (msg [3] > group_offset [msg [1]])
                        group_offset [msg [1]] = msg [3];

                    // we just finished this group, so
                    // start an adaptive write for this file
                    if (msg [1] == msg [2])
                    {
                        groups_complete++;
                        group_state [msg [1]] = STATE_COMPLETE;
                        if (groups_complete != md->groups)
                        {
                            int i = (msg [1] + 1) % md->groups;
                            int groups_tried = 0;
                            while (groups_tried != md->groups)
                            {
                                groups_tried++;
                                if (group_state [i] == STATE_WRITING)
                                {
                                    // tell the subcoordinator to write
                                    // to this file
                                    adaptive_writes_outstanding++;
#if 0
fprintf (fp_c, "adding adaptive write 1 (%d)\n", adaptive_writes_outstanding);
fflush (fp_c);
#endif
                                    if (sub_coord_ranks [i] != md->rank)
                                    {
                                        struct mpi_flag_struct * f = malloc
                                             (sizeof (struct mpi_flag_struct));
                                        f->target = sub_coord_ranks [i];
                                        f->tag = TAG_SUB_COORDINATOR;
                                        f->size = 8 * PARAMETER_COUNT;
                                        f->msg = malloc (8 * PARAMETER_COUNT);
                                        ((uint64_t *) f->msg) [0] =
                                                    ADAPTIVE_WRITE_START;
                                        ((uint64_t *) f->msg) [1] = msg [1];
                                        ((uint64_t *) f->msg) [2] =
                                                    sub_coord_ranks [msg [1]];
                                        ((uint64_t *) f->msg) [3] =
                                                    group_offset [msg [1]];
                                        pthread_mutex_lock (&md->mpi_mutex);
                                        queue_enqueue (&md->mpi_flag, f);
                                        pthread_mutex_unlock (&md->mpi_mutex);
                                        pthread_cond_signal (&md->mpi_cond);
                                    }
                                    else
                                    {
                                        struct flag_struct * f = malloc
                                              (sizeof (struct flag_struct));
                                        f->source = md->rank;
                                        f->tag = TAG_SUB_COORDINATOR;
                                        f->msg = malloc (8 * PARAMETER_COUNT);
                                        ((uint64_t *) f->msg) [0] =
                                                    ADAPTIVE_WRITE_START;
                                        ((uint64_t *) f->msg) [1] = msg [1];
                                        ((uint64_t *) f->msg) [2] =
                                                    sub_coord_ranks [msg [1]];
                                        ((uint64_t *) f->msg) [3] =
                                                    group_offset [msg [1]];
                                        pthread_mutex_lock
                                             (&md->sub_coordinator_mutex);
                                        queue_enqueue
                                               (&md->sub_coordinator_flag, f);
                                        pthread_mutex_unlock
                                              (&md->sub_coordinator_mutex);
                                        pthread_cond_signal
                                              (&md->sub_coordinator_cond);
                                    }
                                    break;
                                }
                                i = (i + 1) % md->groups;
                            }
                        }
                        else // start index collection if no adaptive left
                        {
                            if (!adaptive_writes_outstanding)
                            {
#if 0
fprintf (fp_c, "start_index_collection 1\n");
fflush (fp_c);
#endif
                                start_index_collection = 1;
                            }
                        }
                    }
                    else  // move to the next adaptive writer
                    {
                        // tell the target that we wrote to it and such
                        if (md->group != msg [1])
                        {
                            struct mpi_flag_struct * f = malloc
                                          (sizeof (struct mpi_flag_struct));
                            f->target = sub_coord_ranks [msg [1]];
                            f->tag = TAG_SUB_COORDINATOR;
                            f->size = 8 * PARAMETER_COUNT;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                            ((uint64_t *) f->msg) [1] = msg [1];
                            ((uint64_t *) f->msg) [2] = msg [2];
                            ((uint64_t *) f->msg) [3] = msg [3];
                            ((uint64_t *) f->msg) [4] = msg [4];
                            ((uint64_t *) f->msg) [5] = msg [5];
                            pthread_mutex_lock (&md->mpi_mutex);
                            queue_enqueue (&md->mpi_flag, f);
                            pthread_mutex_unlock (&md->mpi_mutex);
                            pthread_cond_signal (&md->mpi_cond);
                        }
                        else
                        {
                            struct flag_struct * f = malloc
                                             (sizeof (struct flag_struct));
                            f->source = msg [5];
                            f->tag = TAG_SUB_COORDINATOR;
                            f->msg = malloc (8 * PARAMETER_COUNT);
                            ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                            ((uint64_t *) f->msg) [1] = msg [1];
                            ((uint64_t *) f->msg) [2] = msg [2];
                            ((uint64_t *) f->msg) [3] = msg [3];
                            ((uint64_t *) f->msg) [4] = msg [4];
                            ((uint64_t *) f->msg) [5] = msg [5];
                            pthread_mutex_lock (&md->sub_coordinator_mutex);
                            queue_enqueue (&md->sub_coordinator_flag, f);
                            pthread_mutex_unlock (&md->sub_coordinator_mutex);
                            pthread_cond_signal (&md->sub_coordinator_cond);
                        }
                        adaptive_writes_outstanding--;
#if 0
fprintf (fp_c, "removing adaptive write 2 (%d)\n", adaptive_writes_outstanding);
fflush (fp_c);
#endif
                        int i = (msg [2] + 1) % md->groups;
                        int groups_tried = 0;
                        // if we have completed, don't even try to adapt
                        if (groups_complete == md->groups)
                            groups_tried = md->groups;
                        // make sure we can keep offloading
                        // by allowing using the same group over and over
                        while (groups_tried != md->groups)
                        {
                            groups_tried++;
                            if (group_state [i] == STATE_WRITING)
                            {
                                // tell the subcoordinator to write
                                // to this file
                                adaptive_writes_outstanding++;
#if 0
fprintf (fp_c, "adding adaptive write 3 (%d)\n", adaptive_writes_outstanding);
fflush (fp_c);
#endif
                                if (sub_coord_ranks [i] != md->rank)
                                {
                                    struct mpi_flag_struct * f = malloc
                                           (sizeof (struct mpi_flag_struct));
                                    f->target = sub_coord_ranks [i];
                                    f->tag = TAG_SUB_COORDINATOR;
                                    f->size = 8 * PARAMETER_COUNT;
                                    f->msg = malloc (8 * PARAMETER_COUNT);
                                    ((uint64_t *) f->msg) [0] =
                                                ADAPTIVE_WRITE_START;
                                    ((uint64_t *) f->msg) [1] = msg [1];
                                    ((uint64_t *) f->msg) [2] =
                                                sub_coord_ranks [msg [1]];
                                    ((uint64_t *) f->msg) [3] =
                                                group_offset [msg [1]];
                                    pthread_mutex_lock (&md->mpi_mutex);
                                    queue_enqueue (&md->mpi_flag, f);
                                    pthread_mutex_unlock (&md->mpi_mutex);
                                    pthread_cond_signal (&md->mpi_cond);
                                }
                                else
                                {
                                    struct flag_struct * f = malloc
                                                (sizeof (struct flag_struct));
                                    f->source = md->rank;
                                    f->tag = TAG_SUB_COORDINATOR;
                                    f->msg = malloc (8 * PARAMETER_COUNT);
                                    ((uint64_t *) f->msg) [0] =
                                                ADAPTIVE_WRITE_START;
                                    ((uint64_t *) f->msg) [1] = msg [1];
                                    ((uint64_t *) f->msg) [2] =
                                                sub_coord_ranks [msg [1]];
                                    ((uint64_t *) f->msg) [3] =
                                                group_offset [msg [1]];
                                    pthread_mutex_lock
                                               (&md->sub_coordinator_mutex);
                                    queue_enqueue
                                               (&md->sub_coordinator_flag, f);
                                    pthread_mutex_unlock
                                               (&md->sub_coordinator_mutex);
                                    pthread_cond_signal
                                              (&md->sub_coordinator_cond);
                                }

                                break;
                            }
                            i = (i + 1) % md->groups;
                        }
                        if (   !adaptive_writes_outstanding
                            && groups_complete == md->groups
                           )
                        {
#if 0
fprintf (fp_c, "start_index_collection 2\n");
fflush (fp_c);
#endif
                            start_index_collection = 1;
                        }
                    }
                    break;
                }

                case WRITERS_BUSY: // we told it to adaptive write, but
                                   // all are writing already for that group
                {
                    adaptive_writes_outstanding--;
#if 0
fprintf (fp_c, "removing adaptive write 4 (%d)\n", adaptive_writes_outstanding);
fflush (fp_c);
#endif
                    // the writer that told us this can't take any more
                    // adaptive requests so mark as occupied (only if still
                    // marked as writing)
                    if (group_state [msg [2]] == STATE_WRITING)
                        group_state [msg [2]] = STATE_WRITERS_ALL_OCCUPIED;

                    // look for the next group we can ask to write
                    int i = (msg [2] + 1) % md->groups;
                    while (i != msg [2])  // look at all
                    {
                        if (   i != msg [1]  // don't write to the source group
                            && group_state [i] == STATE_WRITING
                           )
                        {
                            // tell the subcoordinator to write
                            // to this file
                            adaptive_writes_outstanding++;
#if 0
fprintf (fp_c, "adding adaptive write 5 (%d)\n", adaptive_writes_outstanding);
fflush (fp_c);
#endif
                            if (sub_coord_ranks [i] != md->rank)
                            {
                                struct mpi_flag_struct * f = malloc
                                            (sizeof (struct mpi_flag_struct));
                                f->target = sub_coord_ranks [i];
                                f->tag = TAG_SUB_COORDINATOR;
                                f->size = 8 * PARAMETER_COUNT;
                                f->msg = malloc (8 * PARAMETER_COUNT);
                                ((uint64_t *) f->msg) [0] =
                                            ADAPTIVE_WRITE_START;
                                ((uint64_t *) f->msg) [1] = msg [1];
                                ((uint64_t *) f->msg) [2] =
                                            sub_coord_ranks [msg [1]];
                                ((uint64_t *) f->msg) [3] =
                                            group_offset [msg [1]];
                                pthread_mutex_lock (&md->mpi_mutex);
                                queue_enqueue (&md->mpi_flag, f);
                                pthread_mutex_unlock (&md->mpi_mutex);
                                pthread_cond_signal (&md->mpi_cond);

                                break;
                            }
                            else
                            {
                                struct flag_struct * f = malloc
                                               (sizeof (struct flag_struct));
                                f->source = md->rank;
                                f->tag = TAG_SUB_COORDINATOR;
                                f->msg = malloc (8 * PARAMETER_COUNT);
                                ((uint64_t *) f->msg) [0] =
                                            ADAPTIVE_WRITE_START;
                                ((uint64_t *) f->msg) [1] = msg [1];
                                ((uint64_t *) f->msg) [2] =
                                            sub_coord_ranks [msg [1]];
                                ((uint64_t *) f->msg) [3] =
                                            group_offset [msg [1]];
                                pthread_mutex_lock
                                              (&md->sub_coordinator_mutex);
                                queue_enqueue (&md->sub_coordinator_flag, f);
                                pthread_mutex_unlock
                                              (&md->sub_coordinator_mutex);
                                pthread_cond_signal (&md->sub_coordinator_cond);

                                break;
                            }
                        }
                        i = (i + 1) % md->groups;
                    }
                    if (   !adaptive_writes_outstanding
                        && groups_complete == md->groups
                       )
                    {
#if 0
fprintf (fp_c, "start_index_collection 3\n");
fflush (fp_c);
#endif
                        start_index_collection = 1;
                    }
                    break;
                }

                case INDEX_SIZE:
                case INDEX_BODY:
                {
#if COLLECT_METRICS
                    if (index_sizes_received == 0)
                        gettimeofday (&timing.t9, NULL);
#endif
                    int source_group = msg [1];
                    uint64_t proc_index_size = msg [2];

                    index_sizes_received++;

                    if (largest_index_size < proc_index_size)
                        largest_index_size = proc_index_size;

                    if (index_size < largest_index_size)
                    {
                        index_size = largest_index_size;
                        if (index_buf)
                            ;//free (index_buf);
                        index_buf = malloc (largest_index_size + 1);
                    }

                    char * new_name;
                    int new_name_len;

                    new_name = malloc (  strlen (md->method->base_path)
                                       + strlen (md->fd->name) + 1 + 6
                                      );
                    char split_format [10] = "%s%s.%lld";
                    sprintf (new_name, split_format, md->method->base_path
                            ,md->fd->name, source_group
                            );
                    new_name_len = strlen (new_name);
                    uint64_t buffer_offset_tmp = buffer_offset;
                    buffer_offset += 8 + 4 + 2 + new_name_len;

                    if (source_group != md->group)
                    {
#if COLLECT_METRICS
                        struct timeval a, b, diff;
                        gettimeofday (&a, NULL);
#endif
                        pthread_mutex_lock (&md->mpi_mutex);
                        MPI_Recv (index_buf, proc_index_size
                                  ,MPI_BYTE
                                  ,sub_coord_ranks [source_group]
                                  ,TAG_COORDINATOR
                                  ,md->group_comm, &status
                                  );
                        pthread_mutex_unlock (&md->mpi_mutex);
#if COLLECT_METRICS
                        gettimeofday (&b, NULL);
                        timeval_subtract (&diff, &b, &a);
                        timeval_add (&timing.t27, &timing.t27, &diff);
                        timing.recv_count++;
#endif
                        buffer_write (&buffer, &buffer_size, &buffer_offset
                                     ,index_buf, proc_index_size
                                     );
                    }
                    else
                    {
                        buffer_write (&buffer, &buffer_size, &buffer_offset
                                     ,(void *) msg [3]
                                     ,proc_index_size
                                     );
                    }

                    file_index [file_index_count].file_number =
                                                       file_index_count + 1;
                    file_index [file_index_count].name_len = new_name_len;
                    file_index [file_index_count].name = new_name;
                    file_index [file_index_count].offset = buffer_offset_tmp;
                    file_index [file_index_count].length = proc_index_size
                                                         + new_name_len
                                                         + 2   // name len
                                                         + 4   // number
                                                         + 8;  // length
                    buffer_write (&buffer, &buffer_size, &buffer_offset_tmp
                                 ,&file_index [file_index_count].offset, 8
                                 );
                    buffer_write (&buffer, &buffer_size, &buffer_offset_tmp
                                 ,&file_index [file_index_count].file_number, 4
                                 );
                    buffer_write (&buffer, &buffer_size, &buffer_offset_tmp
                                 ,&file_index [file_index_count].name_len, 2
                                 );
                    buffer_write (&buffer, &buffer_size, &buffer_offset_tmp
                                 ,&file_index [file_index_count].name
                                 ,new_name_len
                                 );
                    file_index_count++;

                    // if we have recieved all, write to the file
                    if (index_sizes_received == md->groups)
                    {
                        adios_write_version_v2 (&buffer, &buffer_size
                                               ,&buffer_offset
                                               );

                        char * new_name;
                        int old_mask;
                        int perm;

                        old_mask = umask (0);
                        umask (old_mask);
                        perm = old_mask ^ 0666;

                        new_name = malloc (  strlen (md->method->base_path)
                                           + strlen (md->fd->name) + 1
                                          );
                        char split_format [10] = "%s%s";
                        sprintf (new_name, split_format, md->method->base_path
                                ,md->fd->name
                                );
                        int f = open (new_name, O_WRONLY
                                              | O_LARGEFILE
                                              | O_CREAT
                                              | O_TRUNC
                                     ,perm
                                     );
                        if (f == -1) printf ("oops! %s\n", new_name);
                        ssize_t s = write (f, buffer, buffer_offset);
                        if (s != buffer_offset)
                        {
                            fprintf (stderr, "Need to do multi-write 4\n");
                        }
                        close (f);
                        //free (new_name);

                        //for (i = 0; i < md->groups; i++)
                        //    free (file_index [i].name);
                        //free (file_index);

                        // we are done so exit the thread
                        shutdown_flag = SHUTDOWN_FLAG;
#if COLLECT_METRICS
                        gettimeofday (&timing.t10, NULL);
#endif
                    }
                    break;
                }

                case SHUTDOWN_FLAG:
                {
                    break;
                }

                default:
                {
                    printf ("Unknown coordinator message: %lld\n", msg [0]);
                    break;
                }
            }

            ///////////////////////////////////////////////////////////////
            // we have finally gotten all of the writes completed, so
            // start the index collection
            ///////////////////////////////////////////////////////////////
            if (   start_index_collection
                && !index_collection_started
                && !adaptive_writes_outstanding
               )
            {
#if DO_INDEX_COLLECTION == 0
#if COLLECT_METRICS
                gettimeofday (&timing.t9, NULL);
#endif
#endif
                index_collection_started = 1;
                int i;
                for (i = 0; i < md->groups; i++)
                {
                    if (i != md->group)
                    {
                        struct mpi_flag_struct * f = malloc
                                            (sizeof (struct mpi_flag_struct));
                        f->target = sub_coord_ranks [i];
                        f->tag = TAG_SUB_COORDINATOR;
                        f->size = 8 * PARAMETER_COUNT;
                        f->msg = malloc (8 * PARAMETER_COUNT);
                        ((uint64_t *) f->msg) [0] = OVERALL_WRITE_COMPLETE;
                        ((uint64_t *) f->msg) [1] = group_offset [i];
                        pthread_mutex_lock (&md->mpi_mutex);
                        queue_enqueue (&md->mpi_flag, f);
                        pthread_mutex_unlock (&md->mpi_mutex);
                        pthread_cond_signal (&md->mpi_cond);
                    }
                    else
                    {
                        struct flag_struct * f = malloc
                                               (sizeof (struct flag_struct));
                        f->source = md->rank;
                        f->tag = TAG_SUB_COORDINATOR;
                        f->msg = malloc (8 * PARAMETER_COUNT);
                        ((uint64_t *) f->msg) [0] = OVERALL_WRITE_COMPLETE;
                        ((uint64_t *) f->msg) [1] = group_offset [i];
                        pthread_mutex_lock (&md->sub_coordinator_mutex);
                        queue_enqueue (&md->sub_coordinator_flag, f);
                        pthread_mutex_unlock (&md->sub_coordinator_mutex);
                        pthread_cond_signal (&md->sub_coordinator_cond);
                    }
                }
#if DO_INDEX_COLLECTION == 0
#if COLLECT_METRICS
                // while redundant with above, this will get called if
                // we are not doing the index. The other will get call if we do
                gettimeofday (&timing.t10, NULL);
#endif
                shutdown_flag = SHUTDOWN_FLAG;
#endif
            }
#if PRINT_MESSAGES
            printf ("c: source: %2d msg: %s B\n"
                   ,f->source, message_to_string_full (msg, 'c')
                   );
#endif
        }

        if (f)
        {
            free (f->msg);
            free (f);
            f = 0;
            msg = 0;
        }
//        if (xxx >= 300) {shutdown_flag = SHUTDOWN_FLAG; fprintf (stderr, "rank:\t%d\taborting coord timed out last message: %s\n", md->rank, message_to_string (last_message));}
    } while (shutdown_flag != SHUTDOWN_FLAG);
#if 0
fclose (fp_c);
#endif

    //if (index_buf)
    //    free (index_buf);
    if (buffer)
        ;//free (buffer);

    // if we are not doing index collection, exit the thread here
    // but wait for the overall_write_complete messages to be
    // delivered so that those procs can shut down too.
    while (   queue_size (&md->mpi_flag) != 0
           || queue_size (&md->sub_coordinator_flag) != 0
          )
        ;
    md->coord_shutdown_flag = SHUTDOWN_FLAG;
    pthread_cond_signal (&md->mpi_cond);
    pthread_exit (NULL);
    return NULL;
}

static void * writer_main (void * param)
{
    struct thread_struct * t = (struct thread_struct *) param;

    struct adios_method_struct * method = t->method;
    struct adios_adaptive_data_struct * md = t->md;
    struct adios_file_struct * fd = t->fd;

    struct flag_struct * f = 0;
    uint64_t * msg = 0;
    int message_available;
    MPI_Status status;
    MPI_Request req;

    char * buffer = 0;
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;
    uint64_t index_start = md->b.pg_index_offset;

    char * index_buffer = 0;
    uint64_t index_buffer_size = 0;
    uint64_t index_buffer_offset = 0;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    int rc;
    int shutdown_flag = NO_FLAG;
    struct timespec ts;
    struct timeval tp;
    struct timeval interval;
    interval.tv_sec = 0;
    interval.tv_usec = 500000; // half a second wait

    int xxx = 0;
    enum MESSAGE_FLAGS last_message = NO_FLAG;
    do
    {
        xxx++;
//        if (xxx >= 300) {shutdown_flag = SHUTDOWN_FLAG; fprintf (stderr, "rank:\t%d\taborting writer timed out last message: %s\n", md->rank, message_to_string (last_message));}
        //while (queue_size (&md->writer_flag) == 0)
        //    ;
        pthread_mutex_lock (&md->writer_mutex);
        gettimeofday (&tp, NULL);
        timeval_add (&tp, &tp, &interval);
        ts.tv_sec = tp.tv_sec;
        ts.tv_nsec = tp.tv_usec * 1000;
#if PRINT_MESSAGES
        printf ("cond wait WRITER %d\n", md->rank);
#endif
        rc = pthread_cond_timedwait (&md->writer_cond, &md->writer_mutex, &ts);
#if PRINT_MESSAGES
        printf ("signaled WRITER %d\n", md->rank);
#endif
        pthread_mutex_unlock (&md->writer_mutex);

        if (rc != ETIMEDOUT || queue_size (&md->writer_flag) != 0)
        {
            message_available = 1;
            pthread_mutex_lock (&md->writer_mutex);
            queue_dequeue (&md->writer_flag, &f);
            pthread_mutex_unlock (&md->writer_mutex);
            msg = (uint64_t *) f->msg;
#if PRINT_MESSAGES
            printf ("w: source: %2d rank: %2d msg: %s B\n"
                   ,f->source, md->rank, message_to_string_full (msg, 'w')
                   );
#endif
        }
        else
            continue;   // loop back to try again

        if (message_available) {
            last_message = msg [0];
        switch (msg [0])
        {
            case DO_WRITE_FLAG:
            {
#if COLLECT_METRICS
                struct timeval a, b;
                gettimeofday (&a, NULL);
#endif
                // set our base offset for building the index
                fd->base_offset = msg [5] * md->stripe_size;
                // build index appending to any existing index
                adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root
                                     ,&md->old_attrs_root
                                     );

                // we need the size of the buffer for responding to the write
                adios_write_index_v1 (&index_buffer, &index_buffer_size
                                     ,&index_buffer_offset
                                     ,0, md->old_pg_root
                                     ,md->old_vars_root
                                     ,md->old_attrs_root
                                     );

                if (msg [1] == msg [3]) // same file
                {
                    if (md->f == -1) printf ("we got a bad file handle\n");
                    lseek (md->f, md->stripe_size * msg [5], SEEK_SET);
                    ssize_t s = write (md->f, fd->buffer, fd->bytes_written);
                    if (s != fd->bytes_written)
                    {
                        fprintf (stderr, "Need to do multi-write 1 (tried: "
                                         "%llu wrote: %lld) errno %d\n"
                                ,fd->bytes_written, s, errno
                                );
                    }
                }
                else // we are adaptive writing and need to use the other file
                {
                    char * new_name;

                    new_name = malloc (  strlen (method->base_path)
                                       + strlen (fd->name) + 1 + 6
                                      ); // 6 extra for '.XXXXX' file number
                    char split_format [10] = "%s%s.%d";
                    sprintf (new_name, split_format, method->base_path
                            ,fd->name, msg [1]
                            );

                    int f = open (new_name, O_WRONLY | O_LARGEFILE);
                    if (f != -1)
                    {
                        lseek (f, md->stripe_size * msg [5], SEEK_SET);
                        ssize_t s = write (f, fd->buffer, fd->bytes_written);
                        if (s != fd->bytes_written)
                        {
                            fprintf (stderr, "Need to do multi-write 2\n");
                        }
                        close (f);
                    }
                    else
                    {
                        fprintf (stderr, "ADAPTIVE WRITE FAILURE. File: %s\n"
                                ,new_name
                                );
                    }
                }

                // respond to sub coord(s) we are done
                // cannot delay the outbound messages or we will end up in a
                // condition where the adaptive write is completed causing
                // the overall_write_complete to be sent before we send our
                // message to the owner of the file we wrote to.
                uint64_t new_offset = msg [5] + 1;
                f->source = msg [2]; // who to tell our index to
                if (md->rank != msg [4])
                {
                    struct mpi_flag_struct * f = malloc
                                            (sizeof (struct mpi_flag_struct));
                    f->target = msg [4];
                    f->tag = TAG_SUB_COORDINATOR;
                    f->size = 8 * PARAMETER_COUNT;
                    f->msg = malloc (8 * PARAMETER_COUNT);
                    ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                    ((uint64_t *) f->msg) [1] = msg [1];
                    ((uint64_t *) f->msg) [2] = msg [3];
                    assert (msg [3] == md->group);
                    ((uint64_t *) f->msg) [3] = new_offset;
                    ((uint64_t *) f->msg) [4] = index_buffer_offset;
                    ((uint64_t *) f->msg) [5] = md->rank;
                    pthread_mutex_lock (&md->mpi_mutex);
                    queue_enqueue (&md->mpi_flag, f);
                    pthread_mutex_unlock (&md->mpi_mutex);
                    pthread_cond_signal (&md->mpi_cond);
                }
                else
                {
                    struct flag_struct * f = malloc
                                              (sizeof (struct flag_struct));
                    f->source = md->rank;
                    f->tag = TAG_SUB_COORDINATOR;
                    f->msg = malloc (8 * PARAMETER_COUNT);
                    ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                    ((uint64_t *) f->msg) [1] = msg [1];
                    ((uint64_t *) f->msg) [2] = msg [3];
                    assert (msg [3] == md->group);
                    ((uint64_t *) f->msg) [3] = new_offset;
                    ((uint64_t *) f->msg) [4] = index_buffer_offset;
                    ((uint64_t *) f->msg) [5] = md->rank;
                    pthread_mutex_lock (&md->sub_coordinator_mutex);
                    queue_enqueue (&md->sub_coordinator_flag, f);
                    pthread_mutex_unlock (&md->sub_coordinator_mutex);
                    pthread_cond_signal (&md->sub_coordinator_cond);
                }

#if 0
                // tell our sub_coord if we are adaptive writing
                if (msg [1] != msg [3])
                {
                    if (md->rank != msg [4])
                    {
                        struct mpi_flag_struct * f = malloc
                                            (sizeof (struct mpi_flag_struct));
                        f->target = msg [4];
                        f->tag = TAG_SUB_COORDINATOR;
                        f->size = 8 * PARAMETER_COUNT;
                        f->msg = malloc (8 * PARAMETER_COUNT);
                        ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                        ((uint64_t *) f->msg) [1] = msg [1];
                        ((uint64_t *) f->msg) [2] = msg [3];
                        assert (msg [3] == md->group);
                        ((uint64_t *) f->msg) [3] = new_offset;
                        ((uint64_t *) f->msg) [4] = index_buffer_offset;
                        ((uint64_t *) f->msg) [5] = md->rank;
                        pthread_mutex_lock (&md->mpi_mutex);
                        queue_enqueue (&md->mpi_flag, f);
                        pthread_mutex_unlock (&md->mpi_mutex);
                        pthread_cond_signal (&md->mpi_cond);
                    }
                    else
                    {
                        struct flag_struct * f = malloc
                                               (sizeof (struct flag_struct));
                        f->source = md->rank;
                        f->tag = TAG_SUB_COORDINATOR;
                        f->msg = malloc (8 * PARAMETER_COUNT);
                        ((uint64_t *) f->msg) [0] = WRITE_COMPLETE;
                        ((uint64_t *) f->msg) [1] = msg [1];
                        ((uint64_t *) f->msg) [2] = msg [3];
                        assert (msg [3] == md->group);
                        ((uint64_t *) f->msg) [3] = new_offset;
                        ((uint64_t *) f->msg) [4] = index_buffer_offset;
                        ((uint64_t *) f->msg) [5] = md->rank;
                        pthread_mutex_lock (&md->sub_coordinator_mutex);
                        queue_enqueue (&md->sub_coordinator_flag, f);
                        pthread_mutex_unlock (&md->sub_coordinator_mutex);
                        pthread_cond_signal (&md->sub_coordinator_cond);
                    }
                }
#endif
                    
#if COLLECT_METRICS
                gettimeofday (&b, NULL);
                timeval_subtract (&timing.t8, &b, &a);
#endif
                break;
            }

            case SEND_INDEX:
            {
#if DO_INDEX_COLLECTION
                if (md->rank != msg [2])
                {
                    struct mpi_flag_struct * f = malloc
                                          (sizeof (struct mpi_flag_struct));
                    f->target = msg [2];
                    f->tag = TAG_SUB_COORDINATOR_INDEX_BODY;
                    f->size = index_buffer_offset;
                    f->msg = index_buffer;
                    pthread_mutex_lock (&md->mpi_mutex);
                    queue_enqueue (&md->mpi_flag);
                    pthread_mutex_unlock (&md->mpi_mutex);
                    pthread_cond_signal (&md->mpi_cond);
                }
                else
                {
                    struct flag_struct * f = malloc
                                             (sizeof (struct flag_struct)));
                    f->source = md->rank;
                    f->tag = TAG_SUB_COORDINATOR_INDEX_BODY;
                    f->msg = malloc (8 * PARAMETER_COUNT);
                    ((uint64_t *) f->msg) [0] = (uint64_t) index_buffer;
                    pthread_mutex_lock (&md->sub_coordinator_mutex);
                    queue_enqueue (&md->sub_coordinator_flag, f);
                    pthread_mutex_unlock (&md->sub_coordinator_mutex);
                    pthread_cond_signal (&md->sub_coordinator_cond);
                }
#endif

                if (buffer)
                {
                    free (buffer);
                    buffer = 0;
                    buffer_size = 0;
                    buffer_offset = 0;
                }

                adios_clear_index_v1 (new_pg_root, new_vars_root
                                     ,new_attrs_root
                                     );
                adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                                     ,md->old_attrs_root
                                     );
                new_pg_root = 0;
                new_vars_root = 0;
                new_attrs_root = 0;
                md->old_pg_root = 0;
                md->old_vars_root = 0;
                md->old_attrs_root = 0;

                shutdown_flag = SHUTDOWN_FLAG;
#if COLLECT_METRICS
                gettimeofday (&timing.t18, NULL);
#endif
                break;
            }
        }
        }
    }
    while (shutdown_flag != SHUTDOWN_FLAG);

    // make sure our message gets processed before we allow the thread to exit
    while (   queue_size (&md->mpi_flag) != 0
           || queue_size (&md->sub_coordinator_flag) != 0
          )
        ;
    md->writer_shutdown_flag = SHUTDOWN_FLAG;
    pthread_cond_signal (&md->mpi_cond);
    pthread_exit (NULL);
    return NULL;
}

static void * probe_main (void * param)
{
    struct thread_struct * t = (struct thread_struct *) param;

    struct adios_method_struct * method = t->method;
    struct adios_adaptive_data_struct * md = t->md;
    struct adios_file_struct * fd = t->fd;

    while (md->probe_shutdown_flag == NO_FLAG)
    {
        int msg;
        MPI_Status status;

#if COLLECT_METRICS
                struct timeval a, b, diff;
                gettimeofday (&a, NULL);
#endif
                pthread_mutex_lock (&md->mpi_mutex);
        MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, md->group_comm
                   ,&msg, &status
                   );
                pthread_mutex_unlock (&md->mpi_mutex);
#if COLLECT_METRICS
                gettimeofday (&b, NULL);
                timing.iprobe_count++;
                timeval_subtract (&diff, &b, &a);
                timeval_add (&timing.t25, &timing.t25, &diff);
#endif
        if (msg)
        {
#if PRINT_MESSAGES
            printf ("MPI msg available on %d\n", md->rank);
#endif
            pthread_cond_signal (&md->mpi_cond);
#if PRINT_MESSAGES
            printf ("MPI msg signaled on %d\n", md->rank);
#endif
        }
    }

#if PRINT_MESSAGES
    printf ("probe exit: %d\n", md->rank);
#endif
    pthread_exit (NULL);
    return NULL;
}
