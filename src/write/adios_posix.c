           /* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

// see if we have MPI or other tools
#include "config.h"


// xml parser
#include <mxml.h>

#include "public/adios_mpi.h" // MPI or dummy MPI for seq. build
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"
#include "core/adios_logger.h"

#if defined(__APPLE__)  || defined(__CYGWIN__)
#    define O_LARGEFILE 0
#endif

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#define START_TIMER(t) adios_timing_go (fd->group->timing_obj, (t) ) 
#else
#define START_TIMER(t) ; 
#endif

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#define STOP_TIMER(t) adios_timing_stop (fd->group->timing_obj, (t) )
#else
#define STOP_TIMER(t) ;
#endif


static int adios_posix_initialized = 0;

struct adios_POSIX_data_struct
{
    // our file bits
    struct adios_bp_buffer_struct_v1 b;

    // old index structs we read in and have to be merged in
    struct adios_index_struct_v1 * index;

#ifdef HAVE_MPI
    // Metadata file handle
    int mf;
    MPI_Comm group_comm;
    int rank;
    int size;
#endif
    int g_have_mdf;
    int file_is_open; // = 1 if in append mode we leave the file open (close at finalize)
    int index_is_in_memory; // = 1 when index is kept in memory, no need to read from file. =1 after first 'append/update' is completed but not after first 'write'.
    uint64_t pg_start_next; // remember end of PG data for future append steps

    uint64_t total_bytes_written;  /* bytes including all PGs written during one open()..close() with overflows
          index position will be fd->current_pg->pg_start_in_file + total_bytes_written
          it is calculated but not used; fd->current_pg->pg_start_in_file will point to index beginning
          */
};


/* For each group and each method, init is called, 'method' is unique for each call */
void adios_posix_init (const PairStruct * parameters
                      ,struct adios_method_struct * method
                      )
{
    struct adios_POSIX_data_struct * p = 0;

    if (!adios_posix_initialized)
    {
        adios_posix_initialized = 1;
    }

    method->method_data = malloc (sizeof (struct adios_POSIX_data_struct));
    p = (struct adios_POSIX_data_struct *) method->method_data;
    adios_buffer_struct_init (&p->b);
    p->index = adios_alloc_index_v1(1); // with hashtables
#ifdef HAVE_MPI
    p->mf = 0;
    p->group_comm = MPI_COMM_NULL;
    p->rank = 0;
    p->size = 0;
#endif
    p->g_have_mdf = 1;
    p->file_is_open = 0;  // = 1 when posix file is open (used in append mode only)
    p->index_is_in_memory = 0; 
    p->pg_start_next = 0;
    p->total_bytes_written = 0;
}


// Indices for the timer object
#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
static int ADIOS_TIMER_COMM         = ADIOS_TIMING_MAX_USER_TIMERS + 0;
static int ADIOS_TIMER_IO           = ADIOS_TIMING_MAX_USER_TIMERS + 1;
static int ADIOS_TIMER_LOCALMD      = ADIOS_TIMING_MAX_USER_TIMERS + 2;
static int ADIOS_TIMER_GLOBALMD     = ADIOS_TIMING_MAX_USER_TIMERS + 3;
static int ADIOS_TIMER_AD_OPEN      = ADIOS_TIMING_MAX_USER_TIMERS + 4;
static int ADIOS_TIMER_AD_WRITE     = ADIOS_TIMING_MAX_USER_TIMERS + 5;
static int ADIOS_TIMER_AD_OVERFLOW  = ADIOS_TIMING_MAX_USER_TIMERS + 6;
static int ADIOS_TIMER_AD_CLOSE     = ADIOS_TIMING_MAX_USER_TIMERS + 7;
#endif

int adios_posix_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method, MPI_Comm comm
                     )
{
    char * subfile_name = 0;
    char * mdfile_name = 0;
    char * name_with_rank, rank_string[16];
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;

    char *temp_string, *m_size;

    temp_string = (char *) malloc (strlen (method->parameters) + 1);
    strcpy (temp_string, method->parameters);
    trim_spaces (temp_string);

    if ( (m_size = strstr (temp_string, "have_metadata_file")) )
    {
        char * m = strchr (m_size, '=');
        char * n = strtok (m, ";");

        if (!n)
            p->g_have_mdf = atoi (n + 1);
        else
            p->g_have_mdf = atoi (m + 1);
    }
    else
    {
        // by default, write metadata file. 
        p->g_have_mdf = 1;
    }

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
    int timer_count = 8;
    char ** timer_names = (char**) malloc (timer_count * sizeof (char*) );
    timer_names [0] = "Communication";
    timer_names [1] = "I/O";
    timer_names [2] = "Local metadata";
    timer_names [3] = "Global metadata";
    timer_names [4] = "ad_open";
    timer_names [5] = "ad_write";
    timer_names [6] = "ad_overflow";
    timer_names [7] = "ad_close";


    // Ensure both timing objects exist
    // timing_obj should get created at every open
    // prev_timing_obj should only be created at the first open
    if (fd->group)
    {
        if (!fd->group->timing_obj)
            fd->group->timing_obj = adios_timing_create (timer_count, timer_names);

        if (!fd->group->prev_timing_obj)
            fd->group->prev_timing_obj = adios_timing_create (timer_count, timer_names);
    }


#endif

START_TIMER (ADIOS_TIMER_AD_OPEN);

#ifdef HAVE_MPI
    // Need to figure out new the new fd->name, such as restart.bp.0, restart.bp.1....
    p->group_comm = comm;
    if (p->group_comm == MPI_COMM_NULL)
    {
        p->group_comm = MPI_COMM_SELF;
    }

    // if communicator is not MPI_COMM_NULL/MPI_COMM_SELF, subfiles will be generated in a dir.
    if (p->group_comm != MPI_COMM_SELF)
    {
        char * n = strrchr (fd->name, '/');
        if (!n)
        {
            n = fd->name;
        }
        else
        {
            n++;
        }

        MPI_Comm_rank (p->group_comm, &p->rank);
        MPI_Comm_size (p->group_comm, &p->size);
        fd->group->process_id = p->rank;

        sprintf (rank_string, "%d", p->rank);
        // fd->name + '.' + MPI rank + '\0'
        name_with_rank = malloc (strlen (n) + strlen (rank_string) + 2);
        sprintf (name_with_rank, "%s.%s",  n, rank_string);

        // e.g., subfile_name is restart.bp.dir/restart.bp.0
        subfile_name = malloc (strlen (fd->name)
                              + 5
                              + strlen (method->base_path)
                              + strlen (name_with_rank)
                              + 1
                              );
        sprintf (subfile_name, "%s%s%s%s"
                             , fd->name
                             , ".dir/"
                             , method->base_path
                             , name_with_rank
                             );

        mdfile_name = malloc (strlen (method->base_path)
                             + strlen (fd->name)
                             + 1
                             );
        sprintf (mdfile_name, "%s%s"
                            , method->base_path
                            , fd->name
                            );

        free (name_with_rank);
    }
    else
#endif
    {
        // if the communicator is MPI_COMM_SELF, there won't be metadata file generated.
        // The actually subfile name is the one supplied by the user
        subfile_name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
        sprintf (subfile_name, "%s%s", method->base_path, fd->name);
        mdfile_name = 0;
    }

#ifdef HAVE_MPI
    fd->subfile_index = p->rank; // Only if HAVE_MPI
#endif

    p->total_bytes_written = 0; // counts bytes written only in this open()..close() step

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            p->b.f = open (subfile_name, O_RDONLY | O_LARGEFILE);
            if (p->b.f == -1)
            {
                fprintf (stderr, "ADIOS POSIX: file not found: %s\n", fd->name);

                free (subfile_name);

                return 0;
            }
            struct stat s;
            if (fstat (p->b.f, &s) == 0)
                p->b.file_size = s.st_size;
            p->file_is_open = 1;

            break;
        }

        case adios_mode_write:
        {
#ifdef HAVE_MPI
            // create dir to keep all the subfiles
            if (p->group_comm != MPI_COMM_SELF)
            {
                if (p->rank == 0)
                {
                    char * dir_name = malloc (strlen (fd->name) + 4 + 1);
                    sprintf (dir_name, "%s%s"
                                     , fd->name
                                     , ".dir"
                                     ) ;

                    mkdir (dir_name, S_IRWXU | S_IRWXG);
                    free (dir_name);
                }

                MPI_Barrier (p->group_comm);
            }
#endif
            p->b.f = open (subfile_name, O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE
                                       ,  S_IRUSR | S_IWUSR
                                        | S_IRGRP | S_IWGRP
                                        | S_IROTH | S_IWOTH
                            );
            if (p->b.f == -1)
            {
                fprintf (stderr, "adios_posix_open failed for "
                                 "base_path %s, subfile name %s\n"
                        ,method->base_path, subfile_name
                        );

                free (subfile_name);
                free (mdfile_name);

                return 0;
            }

#ifdef HAVE_MPI
            // open metadata file
            START_TIMER (ADIOS_TIMER_GLOBALMD);
            if (p->group_comm != MPI_COMM_SELF && p->g_have_mdf)
            {
                if (p->rank == 0)
                {
                    p->mf = open (mdfile_name, O_WRONLY | O_CREAT | O_TRUNC | O_LARGEFILE
                                  ,  S_IRUSR | S_IWUSR
                                   | S_IRGRP | S_IWGRP
                                   | S_IROTH | S_IWOTH
                             );
                    if (p->mf == -1)
                    {
                        fprintf (stderr, "adios_posix_open failed for "
                                         "base_path %s, metadata file name %s\n"
                                ,method->base_path, mdfile_name
                                );

                        free (subfile_name);
                        free (mdfile_name);

                        return 0;
                    }
                }
            }
            STOP_TIMER (ADIOS_TIMER_GLOBALMD);
#endif
            p->file_is_open = 1;
            p->pg_start_next = 0;

            break;
        }

        case adios_mode_append:
        case adios_mode_update:
        {
            int old_file = 1;
            
#ifdef HAVE_MPI
            /* FIXME: we need to do this synchronisation only if we don't have
               the subdirectory already created. p->index_is_in_memory will be
               true from the second append call, at which point this subdir
               must exist */
            if (p->group_comm != MPI_COMM_SELF && !p->index_is_in_memory)
            {
                if (p->rank == 0)
                {
                    char * dir_name = malloc (strlen (fd->name) + 4 + 1);
                    sprintf (dir_name, "%s%s"
                                     , fd->name
                                     , ".dir"
                                     ) ;

                    mkdir (dir_name, S_IRWXU | S_IRWXG);
                    free (dir_name);
                }

                MPI_Barrier (p->group_comm);
            }
#endif
           

            // if file was not closed in previous append steps, then we are good, otherwise open file
            if (!p->file_is_open) 
            {
                p->b.f = open (subfile_name, O_RDWR | O_LARGEFILE);
                if (p->b.f == -1)
                {
                    old_file = 0;
                    p->b.f = open (subfile_name,  O_WRONLY | O_CREAT | O_LARGEFILE,
                                  S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
                    if (p->b.f == -1)
                    {
                        fprintf (stderr, "adios_posix_open failed to create  file %s\n" ,subfile_name);

                        free (subfile_name);
                        free (mdfile_name);

                        return 0;
                    }
                }
                p->file_is_open = 1;
            }

#ifdef HAVE_MPI
            // open metadata file by rank 0, if user wants to write it
            START_TIMER (ADIOS_TIMER_GLOBALMD);
            if (p->group_comm != MPI_COMM_SELF && 
                p->g_have_mdf &&
                p->rank == 0)
            {
                p->mf = open (mdfile_name, O_WRONLY | O_TRUNC | O_LARGEFILE
                        , S_IRUSR | S_IWUSR
                        | S_IRGRP | S_IWGRP
                        | S_IROTH | S_IWOTH
                        );
                if (p->mf == -1)
                {
                    p->mf = open (mdfile_name, O_WRONLY| O_CREAT | O_LARGEFILE
                            , S_IRUSR | S_IWUSR
                            | S_IRGRP | S_IWGRP
                            | S_IROTH | S_IWOTH
                            );
                    if (p->mf == -1)
                    {
                        fprintf (stderr, "adios_posix_open failed to create  file %s\n" ,mdfile_name);
                        free (subfile_name);
                        free (mdfile_name);

                        return 0;
                    }
                }
            }
            STOP_TIMER (ADIOS_TIMER_GLOBALMD);
#endif
            START_TIMER (ADIOS_TIMER_LOCALMD);
            if (old_file)
            {
                // There is previous data in file. Metadata is in memory 
                // from the second append step. Metadata is only in the 
                // existing file at first append (even if the file was
                // created by a 'w' mode write step in this run

                if (!p->index_is_in_memory)
                {
                    // now we have to read the old stuff so we can merge it
                    // in at the end and set the base_offset for the old index
                    // start
                    uint32_t version;
                    struct stat s;
                    if (fstat (p->b.f, &s) == 0)
                        p->b.file_size = s.st_size;
                    adios_posix_read_version (&p->b);
                    adios_parse_version (&p->b, &version);

                    switch (version & ADIOS_VERSION_NUM_MASK)
                    {
                        case 1:
                        case 2:
                        case 3:
                            // read the old stuff and set the base offset
                            adios_posix_read_index_offsets (&p->b);
                            adios_parse_index_offsets_v1 (&p->b);

                            adios_posix_read_process_group_index (&p->b);
                            adios_parse_process_group_index_v1 (&p->b, &p->index->pg_root, &p->index->pg_tail);

                            // find the largest time index so we can append properly
                            struct adios_index_process_group_struct_v1 * pg;
                            uint32_t max_time_index = 0;
                            pg = p->index->pg_root;
                            while (pg)
                            {
                                if (pg->time_index > max_time_index)
                                    max_time_index = pg->time_index;
                                pg = pg->next;
                            }
                            if (fd->mode == adios_mode_append) {
                                ++max_time_index;
                            }
                            fd->group->time_index = max_time_index;

                            adios_posix_read_vars_index (&p->b);
                            adios_parse_vars_index_v1 (&p->b, &p->index->vars_root, 
                                    p->index->hashtbl_vars,
                                    &p->index->vars_tail);

                            adios_posix_read_attributes_index (&p->b);
                            adios_parse_attributes_index_v1 (&p->b, &p->index->attrs_root);

                            // write data where previous steps finished (minus the index)
                            p->pg_start_next = p->b.end_of_pgs;

                            break;

                        default:
                            fprintf (stderr, "Unknown bp version: %d.  "
                                    "Cannot append\n"
                                    ,version
                                    );

                            free (subfile_name);
                            free (mdfile_name);

                            return 0;
                    }
                }
                else 
                {
                    // index is in memory, update time index and offsets
                    if (fd->mode == adios_mode_append) {
                        fd->group->time_index++;
                    }
                }
            }
            else 
            {
                // there is no previous data, start from offset 0
                p->pg_start_next = 0;
            }
            STOP_TIMER (ADIOS_TIMER_LOCALMD);
            //printf ("adios_posix_open append/update, old_file=%d, index_is_in_memory=%d, "
            //        "pg_start=%" PRId64 "\n",
            //        old_file, p->index_is_in_memory, p>pg_start_next);

            p->index_is_in_memory = 1; // to notify future append steps about the good news
            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);

            free (subfile_name);
            free (mdfile_name);

            return 0;
        }
    }

    if (subfile_name) 
    {
        free (subfile_name);
    }

    if (mdfile_name)
    {
        free (mdfile_name);
    }

    STOP_TIMER (ADIOS_TIMER_AD_OPEN);

    return 1;
}

enum BUFFERING_STRATEGY adios_posix_should_buffer (struct adios_file_struct * fd
                                                  ,struct adios_method_struct * method
                                                  )
{
    return continue_with_new_pg;
}

void adios_posix_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,const void * data
                       ,struct adios_method_struct * method
                       )
{
    START_TIMER (ADIOS_TIMER_AD_WRITE);
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
    STOP_TIMER (ADIOS_TIMER_AD_WRITE);
}

void adios_posix_get_write_buffer (struct adios_file_struct * fd
                                  ,struct adios_var_struct * v
                                  ,uint64_t * size
                                  ,void ** buffer
                                  ,struct adios_method_struct * method
                                  )
{
    // uint64_t mem_allowed;

    // if (*size == 0)
    // {
    //     *buffer = 0;

    //     return;
    // }

    // if (v->adata && v->free_data)
    // {
    //     adios_method_buffer_free (v->data_size);
    //     free (v->adata);
    // }

    // mem_allowed = adios_method_buffer_alloc (*size);
    // if (mem_allowed == *size)
    // {
    //     *buffer = malloc (*size);
    //     if (!*buffer)
    //     {
    //         adios_method_buffer_free (mem_allowed);
    //         fprintf (stderr, "Out of memory allocating %llu bytes for %s\n"
    //                 ,*size, v->name
    //                 );
    //         v->got_buffer = adios_flag_no;
    //         v->free_data = adios_flag_no;
    //         v->data_size = 0;
    //         v->data = 0;
    //         *size = 0;
    //         *buffer = 0;
    //     }
    //     else
    //     {
    //         v->got_buffer = adios_flag_yes;
    //         v->free_data = adios_flag_yes;
    //         v->data_size = mem_allowed;
    //         v->data = *buffer;
    //     }
    // }
    // else
    // {
    //     adios_method_buffer_free (mem_allowed);
    //     fprintf (stderr, "OVERFLOW: Cannot allocate requested buffer of %llu "
    //                      "bytes for %s\n"
    //             ,*size
    //             ,v->name
    //             );
    //     *size = 0;
    //     *buffer = 0;
    // }
}

void adios_posix_read (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * buffer
                      ,uint64_t buffer_size
                      ,struct adios_method_struct * method
                      )
{
    v->data = v->adata = buffer;
    v->data_size = buffer_size;
}

static void adios_posix_write_pg (struct adios_file_struct * fd
                                 ,struct adios_method_struct * method
                                 )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    int32_t to_write;
    uint64_t bytes_written = 0;

    // use offsets set at the end of previous append step
    // fd->current_pg->pg_start_in_file needs to be correctly set before
    // calling adios_build_index_v1()
    fd->current_pg->pg_start_in_file = p->pg_start_next;

    off_t offset = fd->current_pg->pg_start_in_file;

    /* FIXME: what is this check? Should never happen ? */
    assert (p->b.end_of_pgs <= fd->current_pg->pg_start_in_file);
    if (p->b.end_of_pgs > fd->current_pg->pg_start_in_file)
        offset = p->b.end_of_pgs;

    /*printf ("Write PG: pg_start_in_file = %" PRIu64 "  pg_start_next = %" PRIu64 "  "
            "buffer offset = %" PRIu64 "  total_bytes_written = %" PRIu64 "\n",
            fd->current_pg->pg_start_in_file, p->pg_start_next, fd->offset, p->total_bytes_written);*/

    lseek (p->b.f, offset, SEEK_SET);

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
        write (p->b.f, fd->buffer, to_write);
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

    p->total_bytes_written += bytes_written;
    // remember end of PG data (before index) in case overflow or close, or future append needs it
    // also, next PG will be written to this offset in overflow or close()
    p->pg_start_next += bytes_written; 

    /*printf ("Completed PG: pg_start_in_file = %" PRIu64 "  pg_start_next = %" PRIu64 "  "
            "buffer offset = %" PRIu64 "  total_bytes_written = %" PRIu64 "  bytes_written = %" PRIu64 "\n",
            fd->current_pg->pg_start_in_file, p->pg_start_next, fd->offset, p->total_bytes_written, bytes_written);*/

}

static void adios_posix_write_index (struct adios_file_struct * fd
        ,struct adios_method_struct * method
        ,char * buffer
        ,uint64_t buffer_size
        )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    // index location: p->pg_start_next points to the end of the last written PG
    /*printf ("Write Index: pg_start_next = %" PRIu64 "  "
            "buffer size = %" PRIu64 " total_bytes_written = %" PRIu64 "\n",
            p->pg_start_next, buffer_size, p->total_bytes_written);*/

    lseek (p->b.f, p->pg_start_next, SEEK_SET);
    write (p->b.f, buffer, buffer_size);
}

static void adios_posix_do_read (struct adios_file_struct * fd
        ,struct adios_method_struct * method
        )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
        method->method_data;
    struct adios_var_struct * v = fd->group->vars;

    uint32_t version = 0;

    adios_posix_read_version (&p->b);
    adios_parse_version (&p->b, &version);
    version &= ADIOS_VERSION_NUM_MASK;

    switch (version)
    {
        case 1:
        case 2:
        case 3:
        {
            struct adios_index_struct_v1 * index = adios_alloc_index_v1(0); // no hashtables
            struct adios_index_process_group_struct_v1 * pg_root = index->pg_root;
            struct adios_index_process_group_struct_v1 * pg_root_temp = 0;

            adios_posix_read_index_offsets (&p->b);
            adios_parse_index_offsets_v1 (&p->b);

            adios_posix_read_process_group_index (&p->b);
            adios_parse_process_group_index_v1 (&p->b, &pg_root, NULL);
#if 1
            adios_posix_read_vars_index (&p->b);
            adios_parse_vars_index_v1 (&p->b, &index->vars_root, NULL, NULL);

            adios_posix_read_attributes_index (&p->b);
            adios_parse_attributes_index_v1 (&p->b, &index->attrs_root);
#endif

            // the three section headers
            struct adios_process_group_header_struct_v1 pg_header;
            struct adios_vars_header_struct_v1 vars_header;
            struct adios_attributes_header_struct_v1 attrs_header;

            struct adios_var_header_struct_v1 var_header;
            struct adios_var_payload_struct_v1 var_payload;
            struct adios_attribute_struct_v1 attribute;

            int i;

            pg_root_temp = pg_root;
            while (pg_root_temp && pg_root_temp->next)
                pg_root_temp = pg_root_temp->next;

            p->b.read_pg_offset = pg_root_temp->offset_in_file;
            if (pg_root_temp->next)
            {
                p->b.read_pg_size =   pg_root_temp->next->offset_in_file
                                    - pg_root_temp->offset_in_file;
            }
            else
            {
                p->b.read_pg_size =   p->b.pg_index_offset
                                    - pg_root_temp->offset_in_file;
            }

            adios_posix_read_process_group (&p->b);
            adios_parse_process_group_header_v1 (&p->b, &pg_header);

            adios_parse_vars_header_v1 (&p->b, &vars_header);

            for (i = 0; i < vars_header.count; i++)
            {
                memset (&var_payload, 0
                       ,sizeof (struct adios_var_payload_struct_v1)
                       );
                adios_parse_var_data_header_v1 (&p->b, &var_header);

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
                        break;
                }

                if (v1)
                {
                    var_payload.payload = v1->adata;
                    adios_parse_var_data_payload_v1 (&p->b, &var_header
                                                    ,&var_payload
                                                    ,v1->data_size
                                                    );
                }
                else
                {
                    adios_parse_var_data_payload_v1 (&p->b, &var_header
                                                    ,NULL, 0
                                                    );
                }

                adios_clear_var_header_v1 (&var_header);
            }

#if 1
            adios_parse_attributes_header_v1 (&p->b, &attrs_header);

            for (i = 0; i < attrs_header.count; i++)
            {
                adios_parse_attribute_v1 (&p->b, &attribute);
                adios_clear_attribute_v1 (&attribute);
            }
#endif
            adios_clear_process_group_header_v1 (&pg_header);
            adios_clear_index_v1 (index);
            break;
        }

        default:
            fprintf (stderr, "POSIX read: file version unknown: %u\n", version);
            return;
    }

    adios_buffer_struct_clear (&p->b);
}

void adios_posix_buffer_overflow (struct adios_file_struct * fd, 
                                  struct adios_method_struct * method)
{
    START_TIMER (ADIOS_TIMER_AD_OVERFLOW);
    /* We get a full PG at this point that we can write out now. 
       The index should be built in close() after writing the rest of a PG there.
    */
    START_TIMER (ADIOS_TIMER_IO);
    adios_posix_write_pg (fd, method); // Buffered vars written here
    STOP_TIMER (ADIOS_TIMER_IO);

    STOP_TIMER (ADIOS_TIMER_AD_OVERFLOW);
}

void adios_posix_close (struct adios_file_struct * fd
                       ,struct adios_method_struct * method
                       )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    START_TIMER (ADIOS_TIMER_AD_CLOSE);

    switch (fd->mode)
    {
        case adios_mode_write:
        {
            // buffering or not, write the index
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;

            // write buffered data now
            START_TIMER (ADIOS_TIMER_IO);
            adios_posix_write_pg (fd, method); 
            STOP_TIMER (ADIOS_TIMER_IO);

            // Note: adios_posix_write_pg() sets the fd->current_pg->pg_start_in_file
            // to the correct offset, which is used in the build index. If you want 
            // to move the write_pg after building the index, set the offset here in advance:
            //    fd->current_pg->pg_start_in_file = p->pg_start_next;
            //    index_start = p->pg_start_next + fd->offset

            START_TIMER (ADIOS_TIMER_LOCALMD);
            uint64_t index_start = p->pg_start_next;
            // build new index for this step
            adios_build_index_v1 (fd, p->index);
            // if collective, gather the indexes from the rest and call
            // adios_merge_index_v1 (&new_pg_root, &new_vars_root, pg, vars);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, p->index);
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            STOP_TIMER (ADIOS_TIMER_LOCALMD);

#ifdef HAVE_MPI
            START_TIMER (ADIOS_TIMER_GLOBALMD);
            if (p->group_comm != MPI_COMM_SELF && p->g_have_mdf)
            {
                if (p->rank == 0)
                {
                    int * index_sizes = malloc (4 * p->size);
                    int * index_offsets = malloc (4 * p->size);
                    char * recv_buffer = 0;
                    int i;
                    uint32_t size = 0, total_size = 0;

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, p->group_comm
                               );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    for (i = 0; i < p->size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);
                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, p->group_comm
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    char * buffer_save = p->b.buff;
                    uint64_t buffer_size_save = p->b.length;
                    uint64_t offset_save = p->b.offset;

                    for (i = 1; i < p->size; i++)
                    {
                        p->b.buff = recv_buffer + index_offsets [i];
                        p->b.length = index_sizes [i];
                        p->b.offset = 0;

                        adios_parse_process_group_index_v1 (&p->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&p->b, &new_vars_root, NULL, NULL);
                        // do not merge attributes from other processes from 1.4
                        /*
                        adios_parse_attributes_index_v1 (&p->b
                                                        ,&new_attrs_root
                                                        );
                        */

                        adios_merge_index_v1 (p->index, new_pg_root, 
                                              new_vars_root, new_attrs_root, 0);
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }

                    p->b.buff = buffer_save;
                    p->b.length = buffer_size_save;
                    p->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);

                    char * global_index_buffer = 0;
                    uint64_t global_index_buffer_size = 0;
                    uint64_t global_index_buffer_offset = 0;
                    uint64_t global_index_start = 0;
                    uint16_t flag = 0;

                    adios_write_index_v1 (&global_index_buffer, &global_index_buffer_size
                                         ,&global_index_buffer_offset, global_index_start
                                         ,p->index);

                    flag |= ADIOS_VERSION_HAVE_SUBFILE;

                    adios_write_version_flag_v1 (&global_index_buffer
                                                ,&global_index_buffer_size
                                                ,&global_index_buffer_offset
                                                ,flag
                                                );
                    START_TIMER (ADIOS_TIMER_IO);
                    ssize_t s = write (p->mf, global_index_buffer, global_index_buffer_offset);
                    STOP_TIMER (ADIOS_TIMER_IO);
                    if (s != global_index_buffer_offset)
                    {
                        fprintf (stderr, "POSIX method tried to write %" PRIu64 ", "
                                         "only wrote %" PRId64 ". %s:%d\n"
                                         ,fd->bytes_written
                                         ,(int64_t)s
                                         ,__func__, __LINE__
                                );
                    }

                    close (p->mf);
                    free (global_index_buffer);
                }
                else
                {
                    START_TIMER (ADIOS_TIMER_COMM);
                    // Added this explicit cast to avoid truncation of low-order bytes on BGP
                    int i_buffer_size = (int) buffer_size;
                    MPI_Gather (&i_buffer_size, 1, MPI_INT
                               ,0, 0, MPI_INT
                               ,0, p->group_comm
                               );

                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, p->group_comm
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);
                }
            }
            STOP_TIMER (ADIOS_TIMER_GLOBALMD);
#endif

            // write buffered index now
            START_TIMER (ADIOS_TIMER_IO);
            adios_posix_write_index (fd, method, buffer, buffer_offset); 
            STOP_TIMER (ADIOS_TIMER_IO);

            // close the file assuming we are done in 'w' mode
            adios_posix_close_internal (&p->b);
            p->file_is_open = 0;
            // in 'w' mode we forget about index, first append needs to read it from file
            adios_clear_index_v1 (p->index); 
            // notify future append steps that write mode does not keep the index in memory
            p->index_is_in_memory = 0; 

            free (buffer);

            break;
        }

        case adios_mode_append:
        case adios_mode_update:
        {

            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;

            // write buffered data now
            START_TIMER (ADIOS_TIMER_IO);
            adios_posix_write_pg (fd, method); 
            STOP_TIMER (ADIOS_TIMER_IO);

            // Note: adios_posix_write_pg() sets the fd->current_pg->pg_start_in_file
            // to the correct offset, which is used in the build index. If you want 
            // to move the write_pg after building the index, set the offset here in advance:
            //    fd->current_pg->pg_start_in_file = p->pg_start_next;
            //    index_start = p->pg_start_next + fd->offset

            START_TIMER (ADIOS_TIMER_LOCALMD);
            uint64_t index_start = p->pg_start_next;
            // build index for current step and merge it into the 
            // existing index for previous steps
            struct adios_index_struct_v1 * current_index;
            current_index = adios_alloc_index_v1(1); // no hashtables
            adios_build_index_v1 (fd, current_index);
            // new timestep so sorting is not needed during merge
            adios_merge_index_v1 (p->index, current_index->pg_root, 
                    current_index->vars_root, current_index->attrs_root, 0);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, p->index);
            // free current_index structure but do not clear it's content, which is merged
            // into p->index
            adios_free_index_v1 (current_index);
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            STOP_TIMER (ADIOS_TIMER_LOCALMD);

#ifdef HAVE_MPI
            START_TIMER (ADIOS_TIMER_GLOBALMD);
            if (p->group_comm != MPI_COMM_SELF && p->g_have_mdf)
            {
                if (p->rank == 0)
                {
                    int * index_sizes = malloc (4 * p->size);
                    int * index_offsets = malloc (4 * p->size);
                    char * recv_buffer = 0;
                    int i;
                    uint32_t size = buffer_size;
                    uint32_t total_size = 0;

                    // Need to make a temporary copy of p->index and merge
                    // into that. p->index (or rank 0)  must be kept intact for 
                    // future append steps.
                    // Therefore rank 0 will parse it's own buffer too
                    // and merge it into a new index
                    struct adios_index_struct_v1 * gindex;
                    gindex = adios_alloc_index_v1(1); // no hashtables

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, p->group_comm
                               );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    for (i = 0; i < p->size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                        //printf ("index %d offset %d  size %d  total size %u\n",
                        //         i, index_offsets[i], index_sizes[i], total_size);
                    }

                    recv_buffer = malloc (total_size);
                    memcpy (recv_buffer, buffer, index_sizes[0]);

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, p->group_comm
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    char * buffer_save = p->b.buff;
                    uint64_t buffer_size_save = p->b.length;
                    uint64_t offset_save = p->b.offset;

                    for (i = 0; i < p->size; i++)
                    {
                        p->b.buff = recv_buffer + index_offsets [i];
                        p->b.length = index_sizes [i];
                        p->b.offset = 0;

                        //printf ("parse index %d offset %d  size %d\n",
                        //        i, index_offsets[i], index_sizes[i]);
                        adios_parse_process_group_index_v1 (&p->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&p->b, &new_vars_root, NULL, NULL);
                        // do not merge attributes from other processes from 1.4
                        if (i==0) {
                            adios_parse_attributes_index_v1 (&p->b, &new_attrs_root);
                        }

                        // global index would become unsorted on main aggregator during merging 
                        // so sort timesteps in this case (appending)
                        adios_merge_index_v1 (gindex, new_pg_root, 
                                              new_vars_root, new_attrs_root, 1);
                    
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }

                    p->b.buff = buffer_save;
                    p->b.length = buffer_size_save;
                    p->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);

                    char * global_index_buffer = 0;
                    uint64_t global_index_buffer_size = 0;
                    uint64_t global_index_buffer_offset = 0;
                    uint64_t global_index_start = 0;
                    uint16_t flag = 0;

                    adios_write_index_v1 (&global_index_buffer, &global_index_buffer_size
                                         ,&global_index_buffer_offset, global_index_start
                                         ,gindex);

                    flag |= ADIOS_VERSION_HAVE_SUBFILE;

                    adios_write_version_flag_v1 (&global_index_buffer
                                                ,&global_index_buffer_size
                                                ,&global_index_buffer_offset
                                                ,flag
                                                );

                    START_TIMER (ADIOS_TIMER_IO);
                    ssize_t s = write (p->mf, global_index_buffer, global_index_buffer_offset);
                    STOP_TIMER (ADIOS_TIMER_IO);
                    if (s != global_index_buffer_offset)
                    {
                        fprintf (stderr, "POSIX method tried to write %" PRIu64 ", "
                                         "only wrote %" PRId64 ", Mode: a. %s:%d\n"
                                         ,global_index_buffer_offset
                                         ,(int64_t)s
                                         ,__func__, __LINE__
                                );
                    }

                    close (p->mf);

                    free (global_index_buffer);
                    adios_clear_index_v1 (gindex);
                    adios_free_index_v1 (gindex);
                }
                else
                {
                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&buffer_size, 1, MPI_INT
                               ,0, 0, MPI_INT
                               ,0, p->group_comm
                               );

                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, p->group_comm
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);
                }
            }
            STOP_TIMER (ADIOS_TIMER_GLOBALMD);
#endif

            // write buffered index now
            START_TIMER (ADIOS_TIMER_IO);
            adios_posix_write_index (fd, method, buffer, buffer_offset); 
            STOP_TIMER (ADIOS_TIMER_IO);

            free (buffer);

            break;
        }

        case adios_mode_read:
        {
            // read the index to find the place to start reading
            adios_posix_do_read (fd, method);
            struct adios_var_struct * v = fd->group->vars;
            while (v)
            {
                v->data = v->adata = 0;
                v = v->next;
            }

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);

            return;
        }
    }

    STOP_TIMER (ADIOS_TIMER_AD_CLOSE);

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS

    //Finished timing this cycle, swap the timing buffers
    adios_timing_destroy(fd->group->prev_timing_obj);
    fd->group->prev_timing_obj = fd->group->timing_obj;
    fd->group->timing_obj = 0;

    // prev_timing_obj points to unwritten timing info, timing_obj is
    // ready to allocate at the next open

#endif

}

/* For each group's each method, a finalize function is called */ 
void adios_posix_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    if (p->file_is_open) {
        adios_clear_index_v1 (p->index); // append and update methods never cleared the index
        adios_posix_close_internal (&p->b);
        p->file_is_open = 0;
    }
    p->index_is_in_memory = 0; 

    adios_free_index_v1 (p->index);

    if (adios_posix_initialized)
        adios_posix_initialized = 0;
}

void adios_posix_end_iteration (struct adios_method_struct * method)
{
}

void adios_posix_start_calculation (struct adios_method_struct * method)
{
}

void adios_posix_stop_calculation (struct adios_method_struct * method)
{
}
