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

// see if we have MPI or other tools
#include "config.h"


// xml parser
#include <mxml.h>

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "buffer.h"

#if defined(__APPLE__) 
#    define O_LARGEFILE 0
#endif

static int adios_posix_initialized = 0;

struct adios_POSIX_data_struct
{
    // our file bits
    struct adios_bp_buffer_struct_v1 b;

    // old index structs we read in and have to be merged in
    struct adios_index_process_group_struct_v1 * old_pg_root;
    struct adios_index_var_struct_v1 * old_vars_root;
    struct adios_index_attribute_struct_v1 * old_attrs_root;

    uint64_t vars_start;
    uint64_t vars_header_size;
#ifdef HAVE_MPI
    // Metadata file handle
    int mf;
    MPI_Comm group_comm;
    int rank;
    int size;
#endif
};

#ifdef HAVE_MPI
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
                         "Using MPI_COMM_SELF instead\n"
                );

        *comm = MPI_COMM_SELF;
    }
}
#endif

void adios_posix_init (const char * parameters
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
    p->old_pg_root = 0;
    p->old_vars_root = 0;
    p->old_attrs_root = 0;
    p->vars_start = 0;
    p->vars_header_size = 0;
#ifdef HAVE_MPI
    p->mf = 0;
    p->group_comm = MPI_COMM_NULL;
    p->rank = 0;
    p->size = 0;
#endif
}

int adios_posix_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method, void * comm
                     )
{
    char * subfile_name;
    char * mdfile_name;
    int rank;
    char * name_with_rank, rank_string[16];
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
#ifdef HAVE_MPI
    // Need to figure out new the new fd->name, such as restart.bp.0, restart.bp.1....
    adios_var_to_comm (fd->group->group_comm
                      ,fd->group->adios_host_language_fortran
                      ,comm
                      ,&p->group_comm
                      );

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

    fd->subfile_index = p->rank;

    struct stat s;
    if (stat (subfile_name, &s) == 0)
        p->b.file_size = s.st_size;

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
            fd->base_offset = 0;
            fd->pg_start_in_file = 0;

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
            if (p->group_comm != MPI_COMM_SELF)
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
#endif
            fd->base_offset = 0;
            fd->pg_start_in_file = 0;

            break;
        }

        case adios_mode_append:
        {
            int old_file = 1;
#ifdef HAVE_MPI
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
            p->b.f = open (subfile_name, O_RDWR | O_LARGEFILE);
            if (p->b.f == -1)
            {
                old_file = 0;
                p->b.f = open (subfile_name,  O_WRONLY | O_CREAT | O_LARGEFILE
                                ,  S_IRUSR | S_IWUSR
                                 | S_IRGRP | S_IWGRP
                                 | S_IROTH | S_IWOTH
                                );
                if (p->b.f == -1)
                {
                    fprintf (stderr, "adios_posix_open failed for "
                                     "base_path %s, name %s\n"
                            ,method->base_path, fd->name
                            );

                    free (subfile_name);
                    free (mdfile_name);

                    return 0;
                }
            }
#ifdef HAVE_MPI
            // open metadata file
            if (p->group_comm != MPI_COMM_SELF)
            {
                if (p->rank == 0)
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
                            fprintf (stderr, "adios_posix_open failed for "
                                             "base_path %s, name %s\n"
                                    ,method->base_path, fd->name
                                    );

                            free (subfile_name);
                            free (mdfile_name);

                            return 0;
                        }
                    }
                }
            }
#endif
            if (old_file)
            {
                // now we have to read the old stuff so we can merge it
                // in at the end and set the base_offset for the old index
                // start
                uint32_t version;
                adios_posix_read_version (&p->b);
                adios_parse_version (&p->b, &version);

                switch (version & ADIOS_VERSION_NUM_MASK)
                {
                    case 1:
                        // read the old stuff and set the base offset
                        adios_posix_read_index_offsets (&p->b);
                        adios_parse_index_offsets_v1 (&p->b);

                        adios_posix_read_process_group_index (&p->b);
                        adios_parse_process_group_index_v1 (&p->b
                                                           ,&p->old_pg_root
                                                           );

                        // find the largest time index so we can append properly
                        struct adios_index_process_group_struct_v1 * pg;
                        uint32_t max_time_index = 0;
                        pg = p->old_pg_root;
                        while (pg)
                        {
                            if (pg->time_index > max_time_index)
                                max_time_index = pg->time_index;
                            pg = pg->next;
                        }
                        fd->group->time_index = ++max_time_index;

                        adios_posix_read_vars_index (&p->b);
                        adios_parse_vars_index_v1 (&p->b, &p->old_vars_root);

                        adios_posix_read_attributes_index (&p->b);
                        adios_parse_attributes_index_v1 (&p->b
                                                        ,&p->old_attrs_root
                                                        );

                        fd->base_offset = p->b.end_of_pgs;
                        fd->pg_start_in_file = p->b.end_of_pgs;
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

    free (subfile_name);
    free (mdfile_name);

    return 1;
}

enum ADIOS_FLAG adios_posix_should_buffer (struct adios_file_struct * fd
                                          ,struct adios_method_struct * method
                                          )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;

    if (fd->shared_buffer == adios_flag_no && fd->mode != adios_mode_read)
    {
        // write the process group header
        adios_write_process_group_header_v1 (fd, fd->write_size_bytes);

        lseek (p->b.f, fd->base_offset, SEEK_SET);
        ssize_t s = write (p->b.f, fd->buffer, fd->bytes_written);
        if (s != fd->bytes_written)
        {
            fprintf (stderr, "POSIX method tried to write %llu, "
                             "only wrote %llu\n"
                    ,fd->bytes_written
                    ,s
                    );
        }
        fd->base_offset += s;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&p->b);

        // setup for writing vars
        adios_write_open_vars_v1 (fd);
        p->vars_start = lseek (p->b.f, fd->offset, SEEK_CUR);  // save loc
        p->vars_header_size = p->vars_start - fd->base_offset;  // the size
        p->vars_start -= fd->offset; // adjust to start of header
        fd->base_offset += fd->offset;  // add the size of the vars header
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&p->b);
    }

    return fd->shared_buffer;   // buffer if there is space
}

void adios_posix_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,void * data
                       ,struct adios_method_struct * method
                       )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
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
        ssize_t s = write (p->b.f, fd->buffer, fd->bytes_written);
        if (s != fd->bytes_written)
        {
            fprintf (stderr, "POSIX method tried to write %llu, "
                             "only wrote %llu\n"
                    ,fd->bytes_written
                    ,s
                    );
        }
        fd->base_offset += s;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&p->b);

        // write payload
        // adios_write_var_payload_v1 (fd, v);
        uint64_t var_size = adios_get_var_size (v, fd->group, v->data);
        if (fd->base_offset + var_size > fd->pg_start_in_file + fd->write_size_bytes)
            fprintf (stderr, "adios_posix_write exceeds pg bound. File is corrupted. "
                             "Need to enlarge group size. \n"); 

        int32_t to_write;
        uint64_t bytes_written = 0;
        if (var_size > INT32_MAX)
        {
            to_write = INT32_MAX;
        }
        else
        {
            to_write = (int32_t) fd->bytes_written;
        }

        while (bytes_written < var_size)
        {
            bytes_written += write (p->b.f, v->data + bytes_written, to_write);
            if (var_size > bytes_written)
            {
                if (var_size - bytes_written > INT32_MAX)
                {
                    to_write = INT32_MAX;
                }
                else
                {
                    to_write = var_size - bytes_written;
                }
            }
        }

//        s = write (p->b.f, v->data, var_size);
        s = bytes_written;
        if (s != var_size)
        {
            fprintf (stderr, "POSIX method tried to write %llu, "
                             "only wrote %llu\n"
                    ,var_size
                    ,s
                    );
        }
        fd->base_offset += s;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&p->b);
    }
}

void adios_posix_get_write_buffer (struct adios_file_struct * fd
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

void adios_posix_read (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * buffer
                      ,uint64_t buffer_size
                      ,struct adios_method_struct * method
                      )
{
    v->data = buffer;
    v->data_size = buffer_size;
}

static void adios_posix_do_write (struct adios_file_struct * fd
                                 ,struct adios_method_struct * method
                                 ,char * buffer
                                 ,uint64_t buffer_size
                                 )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    int32_t to_write;
    uint64_t bytes_written = 0;

    if (fd->shared_buffer == adios_flag_yes)
    {
        lseek (p->b.f, p->b.end_of_pgs, SEEK_SET);
        if (p->b.end_of_pgs + fd->bytes_written > fd->pg_start_in_file + fd->write_size_bytes)
            fprintf (stderr, "adios_posix_write exceeds pg bound. File is corrupted. "
                             "Need to enlarge group size. \n");

        if (fd->bytes_written > INT32_MAX)
        {
            to_write = INT32_MAX;
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
                if (fd->bytes_written - bytes_written > INT32_MAX)
                {
                    to_write = INT32_MAX;
                }
                else
                {
                    to_write = fd->bytes_written - bytes_written;
                }
            }
        }
    }

    // index location calculation:
    // for buffered, base_offset = 0, fd->offset = write loc
    // for unbuffered, base_offset = write loc, fd->offset = 0
    // for append buffered, base_offset = start, fd->offset = size
    lseek (p->b.f, fd->base_offset + fd->offset, SEEK_SET);
    write (p->b.f, buffer, buffer_size);
}

static void adios_posix_do_read (struct adios_file_struct * fd
                                ,struct adios_method_struct * method
                                )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    struct adios_var_struct * v = fd->group->vars;

    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer = 0;
    data.buffer_len = 0;

    uint32_t version = 0;

    adios_posix_read_version (&p->b);
    adios_parse_version (&p->b, &version);

    switch (version & ADIOS_VERSION_NUM_MASK)
    {
        case 1:
        {
            struct adios_index_process_group_struct_v1 * pg_root = 0;
            struct adios_index_process_group_struct_v1 * pg_root_temp = 0;
            struct adios_index_var_struct_v1 * vars_root = 0;
            struct adios_index_attribute_struct_v1 * attrs_root = 0;

            adios_posix_read_index_offsets (&p->b);
            adios_parse_index_offsets_v1 (&p->b);

            adios_posix_read_process_group_index (&p->b);
            adios_parse_process_group_index_v1 (&p->b, &pg_root);
#if 1
            adios_posix_read_vars_index (&p->b);
            adios_parse_vars_index_v1 (&p->b, &vars_root);

            adios_posix_read_attributes_index (&p->b);
            adios_parse_attributes_index_v1 (&p->b, &attrs_root);
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
                    var_payload.payload = v1->data;
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
            adios_clear_index_v1 (pg_root, vars_root, attrs_root);
            break;
        }

        default:
            fprintf (stderr, "POSIX read: file version unknown: %u\n"
                    ,version
                    );
            return;
    }

    adios_buffer_struct_clear (&p->b);
}

void adios_posix_close (struct adios_file_struct * fd
                       ,struct adios_method_struct * method
                       )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    switch (fd->mode)
    {
        case adios_mode_write:
        {
            if (fd->shared_buffer == adios_flag_no)
            {
                off_t new_off;
                // set it up so that it will start at 0, but have correct sizes
                new_off = lseek (p->b.f, 0, SEEK_CUR);
                fd->offset = fd->base_offset - p->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_vars_v1 (fd);
                // fd->vars_start gets updated with the size written
                fd->offset = lseek (p->b.f, p->vars_start, SEEK_SET);
                ssize_t s = write (p->b.f, fd->buffer, p->vars_header_size);
                if (s != fd->vars_start)
                {
                    fprintf (stderr, "POSIX method tried to write %llu, "
                                     "only wrote %llu\n"
                            ,fd->vars_start
                            ,s
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
                adios_shared_buffer_free (&p->b);

                new_off = lseek (p->b.f, new_off, SEEK_SET);  // go back to end
                adios_write_open_attributes_v1 (fd);
                p->vars_start = lseek (p->b.f, fd->offset, SEEK_CUR); // save loc
                p->vars_header_size = p->vars_start - fd->base_offset;
                p->vars_start -= fd->offset; // adjust to start of header
                fd->base_offset += fd->offset;  // add size of header
                fd->offset = 0;
                fd->bytes_written = 0;

                while (a)
                {
                    adios_write_attribute_v1 (fd, a);
                    if (fd->base_offset + fd->bytes_written > fd->pg_start_in_file + fd->write_size_bytes)
                        fprintf (stderr, "adios_posix_write exceeds pg bound. File is corrupted. "
                                         "Need to enlarge group size. \n");
                    ssize_t s = write (p->b.f, fd->buffer, fd->bytes_written);
                    if (s != fd->bytes_written)
                    {
                        fprintf (stderr, "POSIX method tried to write %llu, "
                                         "only wrote %llu\n"
                                ,fd->bytes_written
                                ,s
                                );
                    }
                    fd->base_offset += s;
                    fd->offset = 0;
                    fd->bytes_written = 0;
                    adios_shared_buffer_free (&p->b);

                    a = a->next;
                }

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - p->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);
                fd->offset = lseek (p->b.f, p->vars_start, SEEK_SET);
                // fd->vars_start gets updated with the size written
                s = write (p->b.f, fd->buffer, p->vars_header_size);
                if (s != p->vars_header_size)
                {
                    fprintf (stderr, "POSIX method tried to write %llu, "
                                     "only wrote %llu\n"
                            ,p->vars_header_size
                            ,s
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
            }

            // buffering or not, write the index
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = fd->base_offset + fd->offset;

            // build index
            adios_build_index_v1 (fd, &p->old_pg_root, &p->old_vars_root
                                 ,&p->old_attrs_root
                                 );
            // if collective, gather the indexes from the rest and call
            // adios_merge_index_v1 (&new_pg_root, &new_vars_root, pg, vars);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, p->old_pg_root, p->old_vars_root
                                 ,p->old_attrs_root
                                 );
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_posix_do_write (fd, method, buffer, buffer_offset);
#ifdef HAVE_MPI
            if (p->group_comm != MPI_COMM_SELF)
            {
                if (p->rank == 0)
                {
                    int * index_sizes = malloc (4 * p->size);
                    int * index_offsets = malloc (4 * p->size);
                    char * recv_buffer = 0;
                    int i;
                    uint32_t size = 0, total_size = 0;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, p->group_comm
                               );

                    for (i = 0; i < p->size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);
                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, p->group_comm
                                );

                    char * buffer_save = p->b.buff;
                    uint64_t buffer_size_save = p->b.length;
                    uint64_t offset_save = p->b.offset;

                    for (i = 1; i < p->size; i++)
                    {
                        p->b.buff = recv_buffer + index_offsets [i];
                        p->b.length = index_sizes [i];
                        p->b.offset = 0;

                        adios_parse_process_group_index_v1 (&p->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&p->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&p->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (&p->old_pg_root
                                             ,&p->old_vars_root
                                             ,&p->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
                                             );
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
                                         ,p->old_pg_root, p->old_vars_root, p->old_attrs_root
                                         );

                    flag |= ADIOS_VERSION_HAVE_SUBFILE;

                    adios_write_version_flag_v1 (&global_index_buffer
                                                ,&global_index_buffer_size
                                                ,&global_index_buffer_offset
                                                ,flag
                                                );
                    ssize_t s = write (p->mf, global_index_buffer, global_index_buffer_offset);
                    if (s != global_index_buffer_offset)
                    {
                        fprintf (stderr, "POSIX method tried to write %llu, "
                                         "only wrote %llu\n"
                                         ,fd->bytes_written
                                         ,s
                                );
                    }

                    close (p->mf);
                }
                else
                {
                    MPI_Gather (&buffer_size, 1, MPI_INT
                               ,0, 0, MPI_INT
                               ,0, p->group_comm
                               );

                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, p->group_comm
                                );
                }
            }
#endif
            free (buffer);

            adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);

            break;
        }

        case adios_mode_append:
        {
            if (fd->shared_buffer == adios_flag_no)
            {
                off_t new_off;
                // set it up so that it will start at 0, but have correct sizes
                new_off = lseek (p->b.f, 0, SEEK_CUR);
                fd->offset = fd->base_offset - p->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_vars_v1 (fd);
                // fd->vars_start gets updated with the size written
                fd->offset = lseek (p->b.f, p->vars_start, SEEK_SET);
                ssize_t s = write (p->b.f, fd->buffer, p->vars_header_size);
                if (s != fd->vars_start)
                {
                    fprintf (stderr, "POSIX method tried to write %llu, "
                                     "only wrote %llu\n"
                            ,fd->vars_start
                            ,s
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
                adios_shared_buffer_free (&p->b);

                new_off = lseek (p->b.f, new_off, SEEK_SET);  // go back to end
                adios_write_open_attributes_v1 (fd);
                p->vars_start = lseek (p->b.f, fd->offset, SEEK_CUR); // save loc
                p->vars_header_size = p->vars_start - fd->base_offset;
                p->vars_start -= fd->offset; // adjust to start of header
                fd->base_offset += fd->offset;  // add size of header
                fd->offset = 0;
                fd->bytes_written = 0;

                while (a)
                {
                    adios_write_attribute_v1 (fd, a);
                    ssize_t s = write (p->b.f, fd->buffer, fd->bytes_written);
                    if (s != fd->bytes_written)
                    {
                        fprintf (stderr, "POSIX method tried to write %llu, "
                                         "only wrote %llu\n"
                                ,fd->bytes_written
                                ,s
                                );
                    }
                    fd->base_offset += s;
                    fd->offset = 0;
                    fd->bytes_written = 0;
                    adios_shared_buffer_free (&p->b);

                    a = a->next;
                }

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - p->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);
                fd->offset = lseek (p->b.f, p->vars_start, SEEK_SET);
                // fd->vars_start gets updated with the size written
                s = write (p->b.f, fd->buffer, p->vars_header_size);
                if (s != p->vars_header_size)
                {
                    fprintf (stderr, "POSIX method tried to write %llu, "
                                     "only wrote %llu\n"
                            ,p->vars_header_size
                            ,s
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
            }

            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = fd->base_offset + fd->offset;

            // build index
            adios_build_index_v1 (fd, &p->old_pg_root, &p->old_vars_root
                                 ,&p->old_attrs_root
                                 );
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, p->old_pg_root, p->old_vars_root
                                 ,p->old_attrs_root
                                 );
#ifdef HAVE_MPI
            if (p->group_comm != MPI_COMM_SELF)
            {
                if (p->rank == 0)
                {
                    int * index_sizes = malloc (4 * p->size);
                    int * index_offsets = malloc (4 * p->size);
                    char * recv_buffer = 0;
                    int i;
                    uint32_t size = 0, total_size = 0;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, p->group_comm
                               );

                    for (i = 0; i < p->size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);
                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, p->group_comm
                                );

                    char * buffer_save = p->b.buff;
                    uint64_t buffer_size_save = p->b.length;
                    uint64_t offset_save = p->b.offset;

                    for (i = 1; i < p->size; i++)
                    {
                        p->b.buff = recv_buffer + index_offsets [i];
                        p->b.length = index_sizes [i];
                        p->b.offset = 0;

                        adios_parse_process_group_index_v1 (&p->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&p->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&p->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (&p->old_pg_root
                                             ,&p->old_vars_root
                                             ,&p->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
                                             );
                    
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }

                    adios_sort_index_v1 (&p->old_pg_root
                                        ,&p->old_vars_root
                                        ,&p->old_attrs_root
                                        );

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
                                         ,p->old_pg_root, p->old_vars_root, p->old_attrs_root
                                         );

                    flag |= ADIOS_VERSION_HAVE_SUBFILE;

                    adios_write_version_flag_v1 (&global_index_buffer
                                                ,&global_index_buffer_size
                                                ,&global_index_buffer_offset
                                                ,flag
                                                );

                    ssize_t s = write (p->mf, global_index_buffer, global_index_buffer_offset);
                    if (s != global_index_buffer_offset)
                    {
                        fprintf (stderr, "POSIX method tried to write %llu, "
                                         "only wrote %llu, Mode: a\n"
                                         ,global_index_buffer_offset
                                         ,s
                                );
                    }

                    close (p->mf);

                    free (global_index_buffer);
                }
                else
                {
                    MPI_Gather (&buffer_size, 1, MPI_INT
                               ,0, 0, MPI_INT
                               ,0, p->group_comm
                               );

                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, p->group_comm
                                );
                }
            }
#endif
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_posix_do_write (fd, method, buffer, buffer_offset);

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
                v->data = 0;
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

    adios_posix_close_internal (&p->b);
    adios_clear_index_v1 (p->old_pg_root, p->old_vars_root, p->old_attrs_root);
    p->old_pg_root = 0;
    p->old_vars_root = 0;
    p->old_attrs_root = 0;
}

void adios_posix_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
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
