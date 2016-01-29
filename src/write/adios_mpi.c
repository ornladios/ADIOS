/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#if defined(__APPLE__)
#    include <sys/param.h>
#    include <sys/mount.h>
#else
#    include <sys/vfs.h>
#endif
#define __USE_LINUX_IOCTL_DEFS
#include <sys/ioctl.h>
#include <assert.h>

// xml parser
#include <mxml.h>

#include "public/adios_mpi.h"
#include "public/adios_error.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"
#include "core/adios_logger.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

static int adios_mpi_initialized = 0;

#define COLLECT_METRICS 0


struct adios_MPI_data_struct
{
    MPI_File fh;
    MPI_Request req;
    MPI_Status status;
    MPI_Comm group_comm;
    MPI_Info info;      // set with base hints for Lustre
    int rank;
    int size;

    struct adios_bp_buffer_struct_v1 b;
    struct adios_index_struct_v1 * index;
};

#if COLLECT_METRICS
// see adios_adaptive_finalize for what each represents
struct timeval_writer
{
    struct timeval t;  // time value
    int pid;            // process id
};

static
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
    struct timeval t0, t5, t6, t7, t8, t11, t12, t13;
    struct timeval t14, t16, t19, t20, t21, t22, t23;
    struct timeval t27, t28;

    uint64_t write_count; // number used
    uint64_t write_size;  // number allocated
    struct timeval * t24;
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

static void print_metric (FILE * f, struct timing_metrics * t, int iteration, int rank, int size, int sub_coord_rank);

static
void print_metrics (struct adios_MPI_data_struct * md, int iteration)
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

        // print the detailed metrics
        FILE * f = fopen ("adios_metrics", "a");
        for (i = 0; i < md->size; i++)
        {
            print_metric (f, &t [i], iteration, i, md->size
                         ,sub_coord_ranks [i]
                         );
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
    }
    MPI_Barrier (md->group_comm);
}

static void print_metric (FILE * f, struct timing_metrics * t, int iteration, int rank, int size, int sub_coord_rank)
{
    struct timeval diff;
    if (rank == 0)
    {
        timeval_subtract (&diff, &t->t6, &t->t5);
        fprintf (f, "dd\t%2d\tMass file open:\t%02d.%06d\n"
                ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t->t5, &t->t16);
        fprintf (f, "ee\t%2d\tBuild file offsets:\t%02d.%06d\n"
                ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t->t11, &t->t0);
        fprintf (f, "hh\t%2d\tTotal time:\t%02d.%06d\n"
                ,iteration, diff.tv_sec, diff.tv_usec);
    }
    
    timeval_subtract (&diff, &t->t13, &t->t12);
    fprintf (f, "ii\t%2d\tLocal index creation:\t%6d\t%02d.%06d\n"
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
    
    timeval_subtract (&diff, &(t->t24 [0]), &t->t21);
    fprintf (f, "pp\t%2d\tshould_buffer->write1 time:\t%6d\t%02d.%06d\n"
            ,iteration, rank, diff.tv_sec, diff.tv_usec);

    fprintf (f, "aa\t%2d\tprocess write time:\t%6d\t%02d.%06d\n"
            ,iteration, rank, t->t8.tv_sec, t->t8.tv_usec);

    timeval_subtract (&diff, &t->t11, &t->t23);
    fprintf (f, "vv\t%2d\tclose start to shutdown time:\t%6d\t%02d.%06d\n"
           ,iteration, rank, (int)diff.tv_sec, (int)diff.tv_usec);

    timeval_subtract (&diff, &t->t28, &t->t23);
    fprintf (f, "ww\t%2d\tclose total time:\t%6d\t%02d.%06d\n"
            ,iteration, rank, (int)diff.tv_sec, (int)diff.tv_usec);
}
#endif


void adios_mpi_init (const PairStruct * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    if (!adios_mpi_initialized)
    {
        adios_mpi_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_MPI_data_struct));
    md = (struct adios_MPI_data_struct *) method->method_data;
    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    MPI_Info_create (&md->info);
    MPI_Info_set (md->info, "romio_ds_read", "disable");
    MPI_Info_set (md->info, "romio_ds_write", "disable");
    MPI_Info_set (md->info, "ind_wr_buffer_size", "16777216");
    md->rank = 0;
    md->size = 0;
    md->group_comm = method->init_comm; // unused here, adios_open will set the current comm
    md->index = adios_alloc_index_v1(1); // with hashtables

    adios_buffer_struct_init (&md->b);
#if COLLECT_METRICS
    // init the pointer for the first go around avoiding the bad free in open
    timing.t24 = 0;
#endif
}

static
void build_read_offsets (struct adios_bp_buffer_struct_v1 * b
                        ,MPI_Offset * offsets, int size, char * group_name
                        ,struct adios_index_struct_v1 * index
                        )
{
    struct adios_index_process_group_struct_v1 * pg_root = index->pg_root;
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

int adios_mpi_open (struct adios_file_struct * fd
                   ,struct adios_method_struct * method, MPI_Comm comm
                   )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

#if COLLECT_METRICS
    gettimeofday (&timing.t0, NULL); // only used on rank == size - 1, but we don't
                              // have the comm yet to get the rank/size
#endif
    adios_buffer_struct_clear (&md->b);

    md->group_comm = comm;
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }
    fd->group->process_id = md->rank;

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

    if (md->rank == md->size - 1)
        next = -1;
    else
        next = md->rank + 1;
    previous = md->rank - 1;
    current = md->rank;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,md->info, &md->fh
                                    );
                if (err != MPI_SUCCESS)
                {
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    adios_error (err_file_open_error, 
                                 "MPI method: open read failed for %s: '%s'\n",
                                 name, e);
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
                adios_parse_process_group_index_v1 (&md->b, &md->index->pg_root, &md->index->pg_tail);

#if 1
                adios_init_buffer_read_vars_index (&md->b);
                MPI_File_seek (md->fh, md->b.vars_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_vars_index_v1 (&md->b, &md->index->vars_root, 
                                           md->index->hashtbl_vars,
                                           &md->index->vars_tail);

                adios_init_buffer_read_attributes_index (&md->b);
                MPI_File_seek (md->fh, md->b.attrs_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.attrs_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_attributes_index_v1 (&md->b, &md->index->attrs_root);
#endif
                // md->b.end_of_pgs points to the end of the last PG in file
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
                    build_read_offsets (&md->b, offsets, md->size
                                        ,fd->group->name, md->index
                                        );
                    MPI_Scatter (offsets, 3, MPI_LONG_LONG
                                ,MPI_IN_PLACE, 3, MPI_LONG_LONG
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

                    MPI_Scatter (0, 3, MPI_LONG_LONG
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
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INT, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_RDONLY
                                    ,md->info
                                    ,&md->fh
                                    );
            }

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                adios_error (err_file_open_error, 
                        "MPI method, rank %d: open read failed for %s: '%s'\n", 
                        md->rank, name, e);
                free (name);

                return adios_flag_no;
            }

            break;
        }

        case adios_mode_write:
        {
#if COLLECT_METRICS                     
            gettimeofday (&timing.t16, NULL);
#endif

#if COLLECT_METRICS
            gettimeofday (&timing.t5, NULL);
#endif

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                MPI_File_delete (name, MPI_INFO_NULL);  // make sure clean

                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,md->info
                                    ,&md->fh
                                    );
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INT, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
                                    ,md->info
                                    ,&md->fh
                                    );
            }

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                adios_error (err_file_open_error, 
                        "MPI method, rank %d: open write failed for %s: '%s'\n", 
                        md->rank, name, e);
                free (name);

                return 0;
            }
#if COLLECT_METRICS
            gettimeofday (&timing.t6, NULL);
#endif

            break;
        }

        case adios_mode_append:
        case adios_mode_update:
        {
            int old_file = 1;

            if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,md->info, &md->fh
                                    );

                if (err != MPI_SUCCESS)
                {
                    old_file = 0;
                    MPI_File_close (&md->fh);
                    err = MPI_File_open (MPI_COMM_SELF, name
                                        ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                        ,md->info, &md->fh
                                        );
                    if (err != MPI_SUCCESS)
                    {
                        char e [MPI_MAX_ERROR_STRING];
                        int len = 0;
                        memset (e, 0, MPI_MAX_ERROR_STRING);
                        MPI_Error_string (err, e, &len);
                        adios_error (err_file_open_error, 
                                "MPI method, rank %d: open for append failed for %s: '%s'\n", 
                                md->rank, name, e);
                        free (name);

                        return 0;
                    }
                }
                MPI_Bcast (&old_file, 1, MPI_INT, 0, md->group_comm);
            }
            else
            {
                if (md->group_comm != MPI_COMM_NULL)
                    MPI_Bcast (&old_file, 1, MPI_INT, 0, md->group_comm);
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

                    adios_parse_process_group_index_v1 (&md->b, &md->index->pg_root, &md->index->pg_tail);

                    // find the largest time index so we can append properly
                    struct adios_index_process_group_struct_v1 * p;
                    uint32_t max_time_index = 0;
                    p = md->index->pg_root;
                    while (p)
                    {
                        if (p->time_index > max_time_index)
                            max_time_index = p->time_index;
                        p = p->next;
                    }
                    if (fd->mode == adios_mode_append) {
                        ++max_time_index;
                    }
                    fd->group->time_index = max_time_index;
                    MPI_Bcast (&fd->group->time_index, 1, MPI_INT, 0
                              ,md->group_comm
                              );

                    adios_init_buffer_read_vars_index (&md->b);
                    MPI_File_seek (md->fh, md->b.vars_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_vars_index_v1 (&md->b, &md->index->vars_root, 
                                               md->index->hashtbl_vars,
                                               &md->index->vars_tail);

                    adios_init_buffer_read_attributes_index (&md->b);
                    MPI_File_seek (md->fh, md->b.attrs_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.attrs_size
                                  ,MPI_BYTE, &md->status
                                  );
                    adios_parse_attributes_index_v1 (&md->b, &md->index->attrs_root);

                    // remember the end of the last PG in file. The new PG written now
                    // will have an offset updated from here. 
                    // md->b.end_of_pgs points to this position

                    MPI_File_close (&md->fh);
                }
                else
                {
                    MPI_Bcast (&fd->group->time_index, 1, MPI_INT, 0
                              ,md->group_comm
                              );
                }
            }
            else
            {
                if (md->rank == 0)
                    MPI_File_close (&md->fh);
            }


            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                // we know it exists, because we created it if it didn't
                // when reading the old file so can just open wronly
                // but adding the create for consistency with write mode
                // so it is easier to merge write/append later
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,md->info
                                    ,&md->fh
                                    );
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INT, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
                                    ,md->info
                                    ,&md->fh
                                    );
            }

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                adios_error (err_file_open_error, 
                        "MPI method, rank %d: open for append failed for %s: '%s'\n", 
                        md->rank, name, e);
                free (name);

                return 0;
            }

            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode, 
                         "MPI method: Unknown file mode requested: %d\n", 
                         fd->mode);

            free (name);

            return 0;
        }
    }

    free (name);

#if COLLECT_METRICS
    gettimeofday (&timing.t22, NULL);
#endif

#if COLLECT_METRICS
    timing.write_count = 0;
    timing.write_size = 0;
    if (timing.t24) free (timing.t24);
    timing.t24 = 0;
#endif

    return 1;
}


static void build_file_offsets (struct adios_MPI_data_struct *md,
                                struct adios_file_struct *fd)
{
    if (md->group_comm != MPI_COMM_NULL)
    {
        if (md->rank == 0)
        {
            // make one space for offset and one for size
            MPI_Offset * offsets = malloc(sizeof (MPI_Offset)
                                           * md->size);
            int i;

            offsets [0] = fd->bytes_written; // = the size of data in buffer on this processor
// mpixlc_r on Eugene doesn't support 64 bit mode. Therefore the following may have problem
// on Eugene for large data size since MPI_LONG_LONG is 32bit 
            MPI_Gather (&(fd->bytes_written), 1, MPI_LONG_LONG
                       ,offsets, 1, MPI_LONG_LONG
                       ,0, md->group_comm);


            uint64_t last_pgsize = offsets [0];
            offsets [0] = md->b.end_of_pgs; // 0 or where to append to existing data
            //printf ("offsets[%d] = {%llu", md->size, offsets[0]);
            for (i = 1; i < md->size; i++)
            {
                uint64_t this_offset = offsets [i];
                offsets [i] = offsets [i - 1] + last_pgsize;
                last_pgsize = this_offset;
                //printf ("  %llu", offsets[i]);
            }
            //printf (" }\n");
            md->b.pg_index_offset =   offsets [md->size - 1]
                                    + last_pgsize;
            //printf (" last_pgsize = %llu, pg_index_offset = %llu\n", last_pgsize, md->b.pg_index_offset);

            MPI_Scatter (offsets, 1, MPI_LONG_LONG
                        ,MPI_IN_PLACE, 1, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            fd->current_pg->pg_start_in_file = offsets[0];
            free (offsets);
        }
        else
        {
            MPI_Offset offset[1];
            offset[0] = fd->bytes_written;

            MPI_Gather (offset, 1, MPI_LONG_LONG
                       ,0, 1, MPI_LONG_LONG
                       ,0, md->group_comm
                       );

            MPI_Scatter (0, 1, MPI_LONG_LONG
                        ,offset, 1, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            fd->current_pg->pg_start_in_file = offset[0];
        }
        //printf ("rank %d: pg_start_in_file = %llu\n", md->rank, fd->current_pg->pg_start_in_file);
    }
    else
    {
        md->b.pg_index_offset = fd->bytes_written;
        fd->current_pg->pg_start_in_file = md->b.end_of_pgs; // 0 or where to append to existing data
    }
}


enum BUFFERING_STRATEGY adios_mpi_should_buffer (struct adios_file_struct * fd
                                                ,struct adios_method_struct * method
                                                )
{
    return stop_on_overflow;
}

void adios_mpi_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,const void * data
                     ,struct adios_method_struct * method
                     )
{
    if (v->got_buffer == adios_flag_yes)
    {
        if (data != v->data)  // if the user didn't give back the same thing
        {
            if (v->free_data == adios_flag_yes)
            {
                free (v->adata);
                adios_method_buffer_free (v->data_size);
            }
        }
        else
        {
            // we already saved all of the info, so we're ok.
            return;
        }
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

void adios_mpi_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
    uint64_t mem_allowed;
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;

    if (*size == 0)
    {
        *buffer = 0;

        return;
    }

    if (v->adata && v->free_data)
    {
        adios_method_buffer_free (v->data_size);
        free (v->adata);
    }

    mem_allowed = adios_method_buffer_alloc (*size);
    if (mem_allowed == *size)
    {
        *buffer = malloc (*size);
        if (!*buffer)
        {
            adios_method_buffer_free (mem_allowed);
            adios_error (err_no_memory,
                         "MPI method, rank %d: cannot allocate %llu bytes for variable %s\n",
                         md->rank, *size, v->name);
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
 
        adios_error (err_buffer_overflow,
                "MPI method, rank %d: OVERFLOW: Cannot get requested ADIOS buffer of %llu "
                "bytes for variable %s. Remaining buffer size was %llu\n",
                md->rank, *size, v->name, mem_allowed);
        *size = 0;
        *buffer = 0;
    }
}

void adios_mpi_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )
{
    v->data = v->adata = buffer;
    v->data_size = buffer_size;
}

static void adios_mpi_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;
    struct adios_var_struct * v = fd->group->vars;

    uint32_t version = md->b.version & ADIOS_VERSION_NUM_MASK;
    switch (version)
    {
        case 1:
        case 2:
        case 3:
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
                    var_payload.payload = v1->adata;
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,&var_payload
                                                    ,v1->data_size
                                                    );
                }
                else
                {
                    log_warn ("MPI method, rank %d: read: skipping name: %s path: %s\n",
                           md->rank, var_header.name, var_header.path);
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
            adios_error (err_invalid_file_version, 
                         "MPI method read: ADIOS file version unknown: %u\n",
                         version);
            return;
    }

    adios_buffer_struct_clear (&md->b);
}

void adios_mpi_buffer_overflow (struct adios_file_struct * fd, 
                                struct adios_method_struct * method)
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    adios_error (err_buffer_overflow, 
            "rank %d: MPI method only works with complete buffering of data between adios_open() "
            "and adios_close(). Portions of global arrays, that do not fit into the "
            "buffer on some processors will not be written by this method to %s\n", 
            md->rank, fd->name);
}


void adios_mpi_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    //struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;
#if COLLECT_METRICS
    gettimeofday (&timing.t23, NULL);
    timing.t19.tv_sec = timing.t23.tv_sec;
    timing.t19.tv_usec = timing.t23.tv_usec;
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
                v->data = v->adata = 0;
                v = v->next;
            }

            break;
        }

        case adios_mode_write:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            int err;

#if COLLECT_METRICS
            gettimeofday (&timing.t19, NULL);
            gettimeofday (&timing.t12, NULL);
#endif
            // figure out the offsets
            // before writing out the buffer and build the index based on target offsets
            build_file_offsets (md, fd);

            // build index appending to any existing index
            adios_build_index_v1 (fd, md->index);

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

                    recv_buffer = malloc (total_size);

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

                        adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                        // do not merge attributes from other processes from 1.4
                        /*
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );
                        */
                        adios_merge_index_v1 (md->index, new_pg_root, 
                                              new_vars_root, new_attrs_root, 0);
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
                                         ,0, md->index);

                    uint32_t tmp_buffer_size = (uint32_t) buffer_offset;
                    MPI_Gather (&tmp_buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->group_comm
                               );
                    MPI_Gatherv (buffer, buffer_offset, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->group_comm
                                );
                }
            }

#if COLLECT_METRICS
            gettimeofday (&timing.t13, NULL);
#endif
            // everyone writes their data
            MPI_File_seek (md->fh, fd->current_pg->pg_start_in_file, MPI_SEEK_SET);
            // if we need to write > 2 GB, need to do it in parts
            // since count is limited to MAX_MPIWRITE_SIZE (signed 32-bit max).
            uint64_t bytes_written = 0;
            int32_t to_write = 0;
            if (fd->bytes_written > MAX_MPIWRITE_SIZE)
            {
                to_write = MAX_MPIWRITE_SIZE;
            }
            else
            {
                to_write = (int32_t) fd->bytes_written;
            }

            while (bytes_written < fd->bytes_written)
            {
                err = MPI_File_write (md->fh, fd->buffer + bytes_written
                        ,to_write, MPI_BYTE, &md->status
                        );
                if (err != MPI_SUCCESS) 
                {              
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    adios_error (err_write_error, 
                            "MPI method, rank %d: adios_close(): writing of buffered data "
                            "[%llu..%llu] to file %s failed: '%s'\n",
                            md->rank, bytes_written, bytes_written+to_write-1, 
                            fd->name, e);       
                }
                bytes_written += to_write;
                if (fd->bytes_written > bytes_written)
                {
                    if (fd->bytes_written - bytes_written > MAX_MPIWRITE_SIZE)
                    {
                        to_write = MAX_MPIWRITE_SIZE;
                    }
                    else
                    {
                        to_write = fd->bytes_written - bytes_written;
                    }
                }
            }

            /* Rank 0 writes the index */
            if (md->rank == 0)
            {
                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,md->b.pg_index_offset, md->index);
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                MPI_File_seek (md->fh, md->b.pg_index_offset, MPI_SEEK_SET);
#if 0
                err = MPI_File_write (md->fh, buffer, buffer_offset, MPI_BYTE
                                     ,&md->status
                                     );
#endif
                {              
                    uint64_t total_written = 0;
                    uint64_t to_write = buffer_offset;
                    int write_len = 0;
                    int count;
                    char * buf_ptr = buffer;
                    while (total_written < buffer_offset)
                    {
                        write_len = (to_write > MAX_MPIWRITE_SIZE) ? MAX_MPIWRITE_SIZE : to_write;
#if COLLECT_METRICS
struct timeval a, b;
gettimeofday (&a, NULL);
#endif
                        err = MPI_File_write (md->fh, buf_ptr, write_len, MPI_BYTE, &md->status);
#if COLLECT_METRICS
gettimeofday (&b, NULL);
timeval_subtract (&timing.t8, &b, &a);
#endif
                        MPI_Get_count(&md->status, MPI_BYTE, &count);
                        if (count != write_len)
                        {
                            log_error ("MPI method, rank %d: Need to do multi-write 1 (tried: "
                                    "%d wrote: %d) errno %d\n",
                                    md->rank, write_len, count, errno);
                            err = count;
                            break;
                        }
                        total_written += count;
                        buf_ptr += count;
                        to_write -= count;
                        //err = total_written;
                    }
                }              
                if (err != MPI_SUCCESS) 
                {              
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    adios_error (err_write_error, 
                            "MPI method, rank %d: adios_close(): writing of index data "
                            "of %llu bytes to file %s failed: '%s'\n",
                            md->rank, buffer_offset, fd->name, e);       
                }
            }
#if COLLECT_METRICS
            gettimeofday (&timing.t20, NULL);
            gettimeofday (&timing.t14, NULL);
#endif

            if (buffer)
            {
                free (buffer);
                buffer = 0;
                buffer_size = 0;
                buffer_offset = 0;
            }

#if COLLECT_METRICS
            gettimeofday (&timing.t11, NULL);
#endif

            break;
        }

        case adios_mode_append:
        case adios_mode_update:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            int err;

            // figure out the offsets
            // before writing out the buffer and build the index based on target offsets
            build_file_offsets (md, fd);

            // build index appending to any existing index
            adios_build_index_v1 (fd, md->index);
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

                    recv_buffer = malloc (total_size);

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

                        adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                        // do not merge attributes from other processes from 1.4
                        /*
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );
                        */
                        adios_merge_index_v1 (md->index, new_pg_root, 
                                              new_vars_root, new_attrs_root, 0);
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
                                         ,0, md->index);

                    int _buffer_size = buffer_offset;

                    MPI_Gather (&_buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->group_comm
                               );
                    MPI_Gatherv (buffer, buffer_offset, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->group_comm
                                );
                }
            }

            /* Write data buffer */
            MPI_File_seek (md->fh, fd->current_pg->pg_start_in_file, MPI_SEEK_SET);
#if 0
            err = MPI_File_write (md->fh, fd->buffer, fd->bytes_written
                    ,MPI_BYTE, &md->status
                    );
#endif
            {              
                uint64_t total_written = 0;
                uint64_t to_write = fd->bytes_written;
                int write_len = 0;
                int count;
                char * buf_ptr = fd->buffer;
                while (total_written < fd->bytes_written)
                {
                    write_len = (to_write > MAX_MPIWRITE_SIZE) ? MAX_MPIWRITE_SIZE : to_write;
                    err = MPI_File_write (md->fh, buf_ptr, write_len, MPI_BYTE, &md->status);
                    MPI_Get_count(&md->status, MPI_BYTE, &count);
                    if (count != write_len)
                    {
                        err = count;
                        break;
                    }
                    total_written += count;
                    buf_ptr += count;
                    to_write -= count;
                    //err = total_written;
                }
            }              
            if (err != MPI_SUCCESS) 
            {              
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                adios_error (err_write_error, 
                        "MPI method, rank %d: adios_close(): writing of buffered data "
                        "of %llu bytes to file %s failed: '%s'\n",
                        md->rank, fd->bytes_written, fd->name, e);       
            }

            /* Rank 0 writes the index */
            if (md->rank == 0)
            {
                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,md->b.pg_index_offset, md->index);
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                MPI_File_seek (md->fh, md->b.pg_index_offset, MPI_SEEK_SET);
#if 0
                err = MPI_File_write (md->fh, buffer, buffer_offset, MPI_BYTE
                                     ,&md->status
                                     );
#endif
                {              
                    uint64_t total_written = 0;
                    uint64_t to_write = buffer_offset;
                    int write_len = 0;
                    int count;
                    char * buf_ptr = buffer;
                    while (total_written < buffer_offset)
                    {
                        write_len = (to_write > MAX_MPIWRITE_SIZE) ? MAX_MPIWRITE_SIZE : to_write;
                        err = MPI_File_write (md->fh, buf_ptr, write_len, MPI_BYTE, &md->status);
                        MPI_Get_count(&md->status, MPI_BYTE, &count);
                        if (count != write_len)
                        {
                            log_error ("MPI method, rank %d: Need to do multi-write 2 (tried: "
                                    "%d wrote: %d) errno %d\n",
                                    md->rank, write_len, count, errno);
                            err = count;
                            break;
                        }
                        total_written += count;
                        buf_ptr += count;
                        to_write -= count;
                        //err = total_written;
                    }
                }              
                if (err != MPI_SUCCESS) 
                {              
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    adios_error (err_write_error, 
                            "MPI method, rank %d: adios_close(): writing of index data "
                            "of %llu bytes to file %s failed: '%s'\n",
                            md->rank, buffer_offset, fd->name, e);       
                }
            }

            free (buffer);

            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode,
                    "MPI method: Unknown file mode: %d in adios_close()\n", 
                    fd->mode);
        }
    }

    if (md && md->fh)
    {
#if COLLECT_METRICS
        MPI_File_sync (md->fh);
#endif
        MPI_File_close (&md->fh);
    }

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

    adios_clear_index_v1 (md->index);
}

void adios_mpi_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    if (adios_mpi_initialized)
    {
        adios_mpi_initialized = 0;
        MPI_Info_free (&md->info);
    }
    adios_free_index_v1 (md->index);
}

void adios_mpi_end_iteration (struct adios_method_struct * method)
{
}

void adios_mpi_start_calculation (struct adios_method_struct * method)
{
}

void adios_mpi_stop_calculation (struct adios_method_struct * method)
{
}
