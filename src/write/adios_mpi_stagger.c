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
#ifndef __APPLE__
#include <sys/vfs.h>
#else
#include <sys/param.h>
#include <sys/mount.h>
#endif
#include <sys/ioctl.h>

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"

static int adios_mpi_stagger_initialized = 0;

#define COLLECT_METRICS 0

struct adios_MPI_data_struct
{
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
    MPI_Comm split_comm;
    int split_size;
    int split_rank;

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
};

#if COLLECT_METRICS
// see adios_adaptive_finalize for what each represents
struct timeval t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25;

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

static
void print_metrics (struct adios_MPI_data_struct * md, int iteration)
{
    struct timeval diff;
    if (md->rank == 0)
    {
        timeval_subtract (&diff, &t2, &t1);
        printf ("cc\t%2d\tFile create (stripe setup):\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t6, &t5);
        printf ("dd\t%2d\tMass file open:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t17, &t16);
        printf ("ee\t%2d\tBuild file offsets:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
    }
    if (md->rank == md->size - 1)
    {
        timeval_subtract (&diff, &t10, &t9);
        printf ("ff\t%2d\tGlobal index creation:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t8, &t7);
        printf ("gg\t%2d\tAll writes complete (w/ local index):\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);

        timeval_subtract (&diff, &t11, &t0);
        printf ("hh\t%2d\tTotal time:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
    }

    timeval_subtract (&diff, &t13, &t12);
    printf ("ii\t%2d\tLocal index creation:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t22, &t21);
    printf ("kk\t%2d\tshould buffer time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t19, &t23);
    printf ("ll\t%2d\tclose startup time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t19, &t0);
    printf ("mm\t%2d\tsetup time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t14, &t20);
    printf ("nn\t%2d\tcleanup time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t21, &t0);
    printf ("oo\t%2d\topen->should_buffer time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t24, &t21);
    printf ("pp\t%2d\tshould_buffer->write1 time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t25, &t24);
    printf ("qq1\t%2d\twrite1->write2 time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);

    timeval_subtract (&diff, &t23, &t25);
    printf ("qq2\t%2d\twrite2->close start time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
}
#endif

static void set_stripe_size (struct adios_file_struct * fd
                            ,struct adios_MPI_data_struct * md
                            ,const char * filename
                            );


void adios_mpi_stagger_init (const PairStruct * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    if (!adios_mpi_stagger_initialized)
    {
        adios_mpi_stagger_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_MPI_data_struct));
    md = (struct adios_MPI_data_struct *) method->method_data;
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
    md->split_comm = MPI_COMM_NULL;
    md->split_size = -1;
    md->split_rank = -1;

    md->files_number = 1;          // always can be 1 by default
    md->max_storage_targets = -1;
    md->max_stripe_count = -1;
    md->min_stripe_count = 1;      // always can do only 1 by default
    md->overlap_factor = 0;        // always can not overlap by default
    md->split_target_count = -1;
    md->split_files_count = split_files_unknown;

    // parse the parameters into key=value segments for optional settings
    if (method->parameters)
    {
        int param_len = strlen (method->parameters);
        if (param_len > 0)
        {
            char * p = strdup (method->parameters);
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
                char * key = malloc (len);
                char * value = malloc (len);
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
                            fprintf (stderr, "MPI_STAGGER: files_number %d "
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
                            fprintf (stderr, "MPI_STAGGER: overlap_factor %d "
                                             "too small. defaulting to 0.\n"
                                    ,v
                                    );

                            v = 0;
                        } else
                        if (v > 99)
                        {
                            fprintf (stderr, "MPI_STAGGER: overlap_factor %d "
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
                        fprintf (stderr, "MPI_STAGGER parameter: key: {%s} "
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

int adios_mpi_stagger_open (struct adios_file_struct * fd
                   ,struct adios_method_struct * method, MPI_Comm comm
                   )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

#if COLLECT_METRICS
    gettimeofday (&t0, NULL); // only used on rank == size - 1, but we don't
                              // have the comm yet to get the rank/size
#endif
    md->group_comm = comm;
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }

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
adios_mpi_build_file_offset(struct adios_MPI_data_struct *md,
                            struct adios_file_struct *fd, char *name)
{
#define SCATTER_PARAMS 5
    if (md->group_comm != MPI_COMM_NULL)
    {
        if (md->rank == 0)
        {
            // make one space for offset and one for size
            MPI_Offset * offsets = malloc(sizeof (MPI_Offset)
                                           * md->size * SCATTER_PARAMS);
            int i;

            offsets [0] = fd->write_size_bytes;
            MPI_Gather (offsets, 1, MPI_LONG_LONG
                       ,offsets, 1, MPI_LONG_LONG
                       ,0, md->group_comm);

// top section: make things a consistent stripe size
// bottom section: just pack the file
#if 1
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
            set_stripe_size (fd, md, name);
#undef STRIPE_INCREMENT
// linear offsets is top, storage target interleaved is bottom
#if 0
            offsets [0 + 0] = fd->base_offset;
            offsets [0 + 1] = biggest_size;
            offsets [0 + 2] = md->storage_targets;
            for (i = 1; i < md->size; i++)
            {
                offsets [i * SCATTER_PARAMS + 0] = offsets [(i - 1) * SCATTER_PARAMS + 0] + biggest_size;
                offsets [i * SCATTER_PARAMS + 1] = biggest_size;
                offsets [i * SCATTER_PARAMS + 2] = md->storage_targets;
            }
            md->b.pg_index_offset =   offsets [(md->size - 1) * SCATTER_PARAMS + 0]
                                    + biggest_size;
#else
            int groups = md->split_groups;
            // since the size of the last group may vary, we can't rely on
            // it for any decisions internally. Instead we will track what the
            // next group should be (or -1 for none) and treat that case
            // special. Hopefully, this won't be more than one group...
            int group = 0;
            int group_size = md->split_size;
                        int sub_groups = md->storage_targets;
                        // how many procs write to each OST
                        int sub_group_size = md->split_size / sub_groups;
                        if (md->split_size % sub_groups > 0)
                            sub_group_size++;
            uint64_t last_offset = 0;   // how far into file is last group?
            int fixed_already = 0;

            int r; // rank within the subgroup for offset purposes
            int sg; // subgroup within the offset
            int * oo = 0;

            for (i = 0; i < md->size; i++)
            {
                group = i / group_size;
                if (group == groups - 1 && !fixed_already)
                {
                    fixed_already = 1;
                    //last_offset = 0;
                    int remainder = md->size - (groups - 1) * group_size;
                    group_size = remainder;
                    sub_group_size = remainder / md->storage_targets;
                    if (remainder % md->storage_targets)
                        sub_group_size++;
                }
                int rank_in_group;
                if (group_size < md->split_size)
                {
                    rank_in_group = i - (md->size - group_size);
                }
                else
                {
                    rank_in_group = i % group_size;
                }
                // which OST for this group we write to
                int ost;
                if (group_size >= md->storage_targets)
                    ost = rank_in_group / (group_size / md->storage_targets);
                else
                    ost = rank_in_group;
                // how many procs write to each OST
                // (since this is updated by the MPI_Comm_split, it may be
                // as small as 1 so we need to take that into account)
                // which split file we are part of
                int sub_group = ost;
                // write order for this OST (duplicate # in each OST)
                int order = rank_in_group % sub_group_size;
                // overall file offset to force each process to write to the
                // desired OST
                int offset;
                if (   (   i % group_size == 0
                        && group_size == md->split_size
                       )
                    || (   group_size < md->split_size
                        && rank_in_group == 0
                       )
                   )
                {
                    oo = (int *) realloc (oo, sizeof (int) * group_size);
                    int j = 0;
                    r = 0;
                    sg = 0;
                    oo [0] = 0;
                    for (j = 1; j < group_size; j++)
                    {
                        sg += sub_group_size;
                        if (sg >= group_size)
                        {
                            r++;
                            sg = r;
                        }
                        oo [sg] = j;
                    }
                    offset = oo [rank_in_group];
                }
                else
                {
                    offset = oo [rank_in_group];
                }
                offsets [i * SCATTER_PARAMS + 0] =   fd->base_offset
                                      + (offset * biggest_size);
                if (last_offset < offsets [i * SCATTER_PARAMS + 0])
                    last_offset = offsets [i * SCATTER_PARAMS + 0];
                offsets [i * SCATTER_PARAMS + 1] = biggest_size;
                offsets [i * SCATTER_PARAMS + 2] = md->storage_targets;
                offsets [i * SCATTER_PARAMS + 3] = md->split_groups;
                offsets [i * SCATTER_PARAMS + 4] = md->split_size;
            }
            if (oo) free (oo);
            // need to use this calc since last grouping may be incomplete
            // with a previous group having a piece later in the file.
            md->b.pg_index_offset =   last_offset
                                    + biggest_size;
#endif
#else
            uint64_t biggest_size = 0;
            uint64_t last_offset = offsets [0];
            offsets [0] = fd->base_offset;
            for (i = 1; i < md->size; i++)
            {
                uint64_t this_offset = offsets [i];
                offsets [i] = offsets [i - 1] + last_offset;
                last_offset = this_offset;
            }
            md->b.pg_index_offset =   offsets [md->size - 1]
                                    + last_offset;
            md->biggest_size = biggest_size;
#endif
            MPI_Scatter (offsets, SCATTER_PARAMS, MPI_LONG_LONG
                        ,offsets, SCATTER_PARAMS, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            fd->base_offset = offsets [0];
            fd->pg_start_in_file = fd->base_offset;
            free (offsets);
        }
        else
        {
            MPI_Offset offset [SCATTER_PARAMS];
            offset [0] = fd->write_size_bytes;

            MPI_Gather (offset, 1, MPI_LONG_LONG
                       ,offset, 1, MPI_LONG_LONG
                       ,0, md->group_comm
                       );

            MPI_Scatter (offset, SCATTER_PARAMS, MPI_LONG_LONG
                        ,offset, SCATTER_PARAMS, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            fd->base_offset = offset [0];
            md->biggest_size = offset [1];
            md->storage_targets = offset [2];
            md->split_groups = offset [3];
            md->split_size = offset [4];
            fd->pg_start_in_file = fd->base_offset;
        }

        // if we should split, make a new comm
        if (md->split_groups > 1)
        {
            // how many groups we split into
            int groups = md->split_groups;

            // how big each group is max (round up)
            int group_size = md->split_size;
            // which group we are part of
            //int group = md->rank / (md->size / groups);
            int group = md->rank / group_size;
            // our rank within our group
            int group_rank = md->rank % group_size;
            // how many procs write to each OST for this group
            int procs_per_ost = group_size / md->storage_targets;

            md->split_group_size = md->split_size; // save old split size for
                                                   // calculating groups
            MPI_Comm_split (md->group_comm, group, group_rank, &md->split_comm);
            MPI_Comm_size (md->split_comm, &md->split_size);
            MPI_Comm_rank (md->split_comm, &md->split_rank);

            // which OST for this group we write to
            // need to take into account that the split size may be < storage
            // targets
            int ost;
            if (md->split_size >= md->storage_targets)
                ost = md->split_rank / (md->split_size / md->storage_targets);
            else 
                ost = md->split_rank / md->split_size;
            // how many files are we splitting into
            int sub_groups = md->storage_targets;
            // how many procs write to each OST
            int sub_group_size = md->split_size / sub_groups;
            if (md->split_size % sub_groups > 0)
                sub_group_size++;
            // which split file we are part of
            int sub_group = ost;
            // write order for this OST (duplicate # in each OST)
            int order = md->split_rank % sub_group_size;
            // overall file offset to force each process to write to the
            // desired OST

            // for the split files, put the index at the end. this way we don't
            // have to send it around to everyone or a separate send call.
            // This should also take into account the case where the last file
            // is smaller (uneven division of procs to files).
            if (group_rank == 0 && md->b.pg_index_offset == 0)
            {
                md->b.pg_index_offset = md->split_size * md->biggest_size;
            }
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
// <transport method="MPI_STAGGER" group="restart">max_storage_targets=672;max_stripe_count=160;files_number=3;overlap_factor=50;min_stripe_count=10</transport>

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
static void calc_stripe_info (struct adios_MPI_data_struct * md
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
static void set_stripe_size (struct adios_file_struct * fd
                            ,struct adios_MPI_data_struct * md
                            ,const char * filename
                            )
{
    struct statfs fsbuf;
    int err;

    int f;
    int old_mask;
    int perm;

#if COLLECT_METRICS
    gettimeofday (&t1, NULL);
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

            // figure out if we need to split to use all of the OSTs (round up)
            int targets_per_file;
            if (md->storage_targets > md->max_stripe_count)
                targets_per_file = md->max_stripe_count;
            else
                targets_per_file = md->storage_targets;
            // assuming our max targets per file, see if we have to split
            md->split_groups =   (md->max_storage_targets / md->files_number)
                               / targets_per_file;
            if ((md->max_storage_targets / md->files_number) % targets_per_file)
                md->split_groups++;

            // make sure we don't make more splits than we can handle with
            // our min stripe count
            if (md->size / md->split_groups < md->min_stripe_count)
                md->split_groups = md->size / md->min_stripe_count;
            if (md->size % md->min_stripe_count)
                md->split_groups++;

            // if we are splitting files, fixup now (need to do above to get
            // the real max storage targets from the filesystem).
            if (   md->split_groups > 1
                || md->split_files_count != split_files_unknown
               )
            {
                // adjust the split_groups based on how we want to split
                // if not specified (or min), just use all of the OSTs
                // in as few files as possible
                if (md->split_files_count == split_files_max)
                {
                    if (md->size / md->split_groups < md->split_target_count)
                    {
                        md->split_groups = md->size / md->split_target_count;
                        if (md->size % md->split_target_count)
                            md->split_groups++;
                    }
                }

                // adjust the storage targets to the number of groups
                md->storage_targets =   md->max_storage_targets
                                      / md->split_groups;
                if (md->storage_targets < md->min_stripe_count)
                    md->storage_targets = md->min_stripe_count;
                if (md->storage_targets > md->max_stripe_count)
                    md->storage_targets = md->max_stripe_count;

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

            // #groups employed for this file is md->split_groups
            // group_size for each part of this file is md->split_size
            md->split_size = md->size / md->split_groups;
            if (md->size % md->split_groups)
                md->split_size++;
        }
    }
#if COLLECT_METRICS         
    gettimeofday (&t2, NULL);
#endif
}

enum ADIOS_FLAG adios_mpi_stagger_should_buffer (struct adios_file_struct * fd
                                        ,struct adios_method_struct * method
                                        )
{
    int i;
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;
    char * name;
    int err;
    int flag;    // used for coordinating the MPI_File_open

    int previous;
    int current;
    int next;

#if COLLECT_METRICS
    gettimeofday (&t21, NULL);
#endif

    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);

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
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }

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
                MPI_File_delete (name, MPI_INFO_NULL);  // make sure clean
            }

#if COLLECT_METRICS                     
            gettimeofday (&t16, NULL);
#endif
            // figure out the offsets and create the file with proper striping
            // before the MPI_File_open is called
            adios_mpi_build_file_offset (md, fd, name);

            if (md->split_groups != 1)
            {
                // if we need to do a file split, we need to fixup the name
                name = realloc (name, (  strlen (method->base_path)
                                       + strlen (fd->name) + 1 + 6
                                      )
                               ); // 6 extra for '.XXXXX' file number
                // which group is rank / group_size
                // group_size is total_size / num_files
                // num_files is max_targets / storage_targets
                // use calculated storage_targets since min_targets might win
                int group;
                group = md->rank / md->split_group_size;
                char split_format [7] = ".%d";
                char split_name [7];
                sprintf (split_name, split_format, group);
                strcat (name, split_name);
            }
#if COLLECT_METRICS
            gettimeofday (&t17, NULL);
#endif
                      
#if COLLECT_METRICS   
            gettimeofday (&t5, NULL);
#endif 

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,MPI_INFO_NULL
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
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }

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
#if COLLECT_METRICS
            gettimeofday (&t6, NULL);
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
                    MPI_Bcast (&fd->group->time_index, 1, MPI_INT, 0
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
            adios_mpi_build_file_offset (md, fd, name);

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
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }

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
    gettimeofday (&t22, NULL);
#endif
    return fd->shared_buffer;
}

void adios_mpi_stagger_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,void * data
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;

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
        uint64_t var_size = adios_get_var_size (v, v->data);
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
    static int writes_seen = 0;

    if (writes_seen == 0) gettimeofday (&t24, NULL);
    else if (writes_seen == 1) gettimeofday (&t25, NULL);
    writes_seen++;
#endif
}

void adios_mpi_stagger_get_write_buffer (struct adios_file_struct * fd
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

void adios_mpi_stagger_read (struct adios_file_struct * fd
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
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;
    struct adios_var_struct * v = fd->group->vars;

    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer = 0;
    data.buffer_len = 0;

    switch (md->b.version & ADIOS_VERSION_NUM_MASK)
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

void adios_mpi_stagger_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;
#if COLLECT_METRICS
    gettimeofday (&t23, NULL);
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
            uint64_t index_start = md->b.pg_index_offset;

            // setup so we can use split files or not
            int index_rank = md->rank;
            int index_size = md->size;
            MPI_Comm index_comm = md->group_comm;

            // if we split, change to each new file for index
            if (md->split_comm != MPI_COMM_NULL)
            {
                index_comm = md->split_comm;
                MPI_Comm_size (index_comm, &index_size);
                MPI_Comm_rank (index_comm, &index_rank);
            }

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

#if COLLECT_METRICS
            gettimeofday (&t19, NULL);
#endif
#if COLLECT_METRICS
            gettimeofday (&t7, NULL);
#endif
#if COLLECT_METRICS
            gettimeofday (&t12, NULL);
#endif
            // build index appending to any existing index
            adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root
                                 ,&md->old_attrs_root
                                 );
            // if collective, gather the indexes from the rest and call
            if (md->group_comm != MPI_COMM_NULL)
            {
                if (index_rank == 0)
                {
                    int * index_sizes = malloc (4 * index_size);
                    int * index_offsets = malloc (4 * index_size);
                    char * recv_buffer = 0;
                    uint32_t size = 0;
                    uint32_t total_size = 0;
                    int i;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, index_comm
                               );

                    for (i = 0; i < index_size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    } 

                    recv_buffer = malloc (total_size);

                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, index_comm
                                );

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    int r = 0; // rank within the subgroup for offset purposes
                    int sg = 0; // subgroup within the offset
                    int * oo = 0;

                    for (i = 1; i < index_size; i++)
                    {
                        // if we are doing a split setup, we need to
                        // take into account our interleave for which
                        // index we place when so we keep them in order
                        int next_item = i;
                        if (md->split_comm != MPI_COMM_NULL)
                        {
                            // how many procs write to each OST
                            int sub_group_size;
                            if (md->split_size >= md->storage_targets)
                            {
                                sub_group_size =   md->split_size
                                                 / md->storage_targets;
                                if (md->split_size % md->storage_targets)
                                    sub_group_size++;
                            }
                            else
                                sub_group_size = 1; // rely on storage_targets

                            int x = i % sub_group_size;
                            int y = i / sub_group_size;
                            next_item = x * sub_group_size + y;

                            if (i == 1)
                            {
                                oo = (int *) realloc (oo,   sizeof (int)
                                                          * index_size
                                                     );
                                int j = 0;
                                r = 0;
                                sg = 0;
                                oo [0] = 0;
                                for (j = 1; j < index_size; j++)
                                {
                                    sg += sub_group_size;
                                    if (sg >= index_size)
                                    {
                                        r++;
                                        sg = r;
                                    }
                                    oo [j] = sg;
                                }
                                next_item = oo [i];
                            }
                            else
                            {
                                next_item = oo [i];
                            }
                        }
                        md->b.buff = recv_buffer + index_offsets [next_item];
                        md->b.length = index_sizes [next_item];
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
                                             ,new_attrs_root, 0
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }
#if COLLECT_METRICS
            gettimeofday (&t9, NULL);
#endif
// do global index creation
#if COLLECT_METRICS
            gettimeofday (&t10, NULL);
#endif
                    if (oo) free (oo);
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
                               ,0, index_comm
                               );
                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, index_comm
                                );
                }
            }
#if COLLECT_METRICS
            gettimeofday (&t13, NULL);
#endif

            if (fd->shared_buffer == adios_flag_yes)
            {
                // if we have a comm and an OST count, stagger writing
                if (md->group_comm != MPI_COMM_NULL && md->storage_targets)
                {
                    MPI_Comm write_comm = md->group_comm;
                    if (md->split_comm != MPI_COMM_NULL)
                        write_comm = md->split_comm;
                    // ordering (4 OSTs and 10 ranks):
                    // 0: A 1: B 2: C 3: D
                    // 4: A 5: B 6: C 7: D
                    // 8: A 9: B
                    // write all As, then Bs, then Cs, then Ds.
                    // each A will send to B, etc. to do the ordering
                    int next_rank = md->rank + 1;
                    int prev_rank = md->rank - 1;
                    int current_rank = md->rank;
                    int flag = 0;

                    if (md->split_groups > 1)
                    {
                        // how many groups we split into
                        int groups = md->split_groups;
                        // which group we are part of
                        int group = md->rank / groups;
                        // how big each group is
                        int group_size = md->split_group_size;
                        // our rank within our group
                        int group_rank = md->rank % group_size;
                        // how many procs write to each OST for this group
                        int procs_per_ost = group_size / md->storage_targets;
                        // how many files are we splitting into
                        int sub_groups = md->storage_targets;
                        // how many procs write to each OST
                        int sub_group_size = md->split_size / sub_groups;
                        if (md->split_size % sub_groups > 0)
                            sub_group_size++;

                        // change the prev/cur/next to be for each OST
                        current_rank = group_rank;
                        next_rank = current_rank + 1;
                        prev_rank = current_rank - 1;
                        if (group_rank == 0 || group_rank % sub_group_size == 0)
                            prev_rank = -1;
                        if (   group_rank % sub_group_size == sub_group_size - 1
                            || md->rank == md->size - 1
                            || md->split_rank == md->split_size - 1
                           )
                            next_rank = -1;
                    }
                    else
                    {
                        if (   next_rank >= md->size
                            ||    current_rank % md->storage_targets
                               == md->storage_targets - 1
                           )
                            next_rank = -1;
                        if (current_rank % md->storage_targets == 0)
                            prev_rank = -1;
                     }

                    if (prev_rank != -1)
                    {
                        MPI_Recv (&flag, 1, MPI_INT, prev_rank, prev_rank
                                 ,write_comm, &md->status
                                 );
                    }
                    //write
                    MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
                    MPI_File_write (md->fh, fd->buffer, fd->bytes_written
                                   ,MPI_BYTE, &md->status
                                   );
                    if (next_rank != -1)
                    {
                        //MPI_Isend (&flag, 1, MPI_INT, next_rank
                        //          ,current_rank, write_comm, &md->req
                        //          );
                        MPI_Send (&flag, 1, MPI_INT, next_rank
                                 ,current_rank, write_comm
                                 );
                    }
                }
                else
                {
                    // everyone writes their data
                    MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
                    MPI_File_write (md->fh, fd->buffer, fd->bytes_written
                                   ,MPI_BYTE, &md->status
                                   );
                }
            }

            if (index_rank == 0)
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
#if COLLECT_METRICS
            gettimeofday (&t8, NULL);
#endif
#if COLLECT_METRICS
            gettimeofday (&t20, NULL);
#endif
#if COLLECT_METRICS
            gettimeofday (&t14, NULL);
#endif
            if (buffer)
            {
                free (buffer);
                buffer = 0;
                buffer_size = 0;
                buffer_offset = 0;
            }

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
#if COLLECT_METRICS
            gettimeofday (&t11, NULL);
            t15.tv_sec = t11.tv_sec;
            t15.tv_usec = t11.tv_usec;
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
                                             ,new_attrs_root, 0
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

    if (md && md->fh)
        MPI_File_close (&md->fh);

    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
    {
        md->group_comm = MPI_COMM_NULL;
    }
    if (md && md->split_comm != MPI_COMM_NULL)
    {
        MPI_Comm_free (&md->split_comm);
        md->split_comm = MPI_COMM_NULL;
        md->split_size = -1;
        md->split_rank = -1;
    }

    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->group_comm = MPI_COMM_NULL;

    adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                         ,md->old_attrs_root
                         );
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
#if COLLECT_METRICS
    print_metrics (md, iteration++);
#endif
}

void adios_mpi_stagger_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
    if (adios_mpi_stagger_initialized)
        adios_mpi_stagger_initialized = 0;
}

void adios_mpi_stagger_end_iteration (struct adios_method_struct * method)
{
}

void adios_mpi_stagger_start_calculation (struct adios_method_struct * method)
{
}

void adios_mpi_stagger_stop_calculation (struct adios_method_struct * method)
{
}
