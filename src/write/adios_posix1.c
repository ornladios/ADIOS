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

static int adios_posix1_initialized = 0;

struct adios_POSIX1_data_struct
{
    // our file bits
    struct adios_bp_buffer_struct_v1 b;

    // old index structs we read in and have to be merged in
    struct adios_index_process_group_struct_v1 * old_pg_root;
    struct adios_index_var_struct_v1 * old_vars_root;
    struct adios_index_attribute_struct_v1 * old_attrs_root;

    uint64_t vars_start;
    uint64_t vars_header_size;
};

void adios_posix1_init (const char * parameters
                      ,struct adios_method_struct * method
                      )
{
    struct adios_POSIX1_data_struct * p = 0;

    if (!adios_posix1_initialized)
    {
        adios_posix1_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_POSIX1_data_struct));
    p = (struct adios_POSIX1_data_struct *) method->method_data;
    adios_buffer_struct_init (&p->b);
    p->old_pg_root = 0;
    p->old_vars_root = 0;
    p->old_attrs_root = 0;
    p->vars_start = 0;
    p->vars_header_size = 0;
}

int adios_posix1_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method, void * comm
                     )
{
    char * name;
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
                                                          method->method_data;

    // figure out the actual name of the file.
    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);
    struct stat s;
    if (stat (name, &s) == 0)
        p->b.file_size = s.st_size;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            p->b.f = open (name, O_RDONLY
#ifndef __APPLE__
| O_LARGEFILE
#endif
);
            if (p->b.f == -1)
            {
                fprintf (stderr, "ADIOS POSIX1: file not found: %s\n", fd->name);

                free (name);

                return 0;
            }
            fd->base_offset = 0;
            fd->pg_start_in_file = 0;

            break;
        }

        case adios_mode_write:
        {
            p->b.f = open (name, O_WRONLY | O_CREAT | O_TRUNC
#ifndef __APPLE__
| O_LARGEFILE
#endif
                            ,  S_IRUSR | S_IWUSR
                             | S_IRGRP | S_IWGRP
                             | S_IROTH | S_IWOTH
                            );
            if (p->b.f == -1)
            {
                fprintf (stderr, "adios_posix1_open failed for "
                                 "base_path %s, name %s\n"
                        ,method->base_path, fd->name
                        );

                free (name);

                return 0;
            }
            fd->base_offset = 0;
            fd->pg_start_in_file = 0;

            break;
        }

        case adios_mode_append:
        {
            int old_file = 1;
            p->b.f = open (name, O_RDWR
#ifndef __APPLE__
| O_LARGEFILE
#endif
);
            if (p->b.f == -1)
            {
                old_file = 0;
                p->b.f = open (name,  O_WRONLY | O_CREAT
#ifndef __APPLE__
| O_LARGEFILE
#endif
                                ,  S_IRUSR | S_IWUSR
                                 | S_IRGRP | S_IWGRP
                                 | S_IROTH | S_IWOTH
                                );
                if (p->b.f == -1)
                {
                    fprintf (stderr, "adios_posix1_open failed for "
                                     "base_path %s, name %s\n"
                            ,method->base_path, fd->name
                            );

                    free (name);

                    return 0;
                }
            }

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

                        free (name);

                        return 0;
                }
            }

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);

            free (name);

            return 0;
        }
    }

    free (name);

    return 1;
}

enum ADIOS_FLAG adios_posix1_should_buffer (struct adios_file_struct * fd
                                          ,struct adios_method_struct * method
                                          )
{
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
                                                          method->method_data;

    if (fd->shared_buffer == adios_flag_no && fd->mode != adios_mode_read)
    {
        // write the process group header
        adios_write_process_group_header_v1 (fd, fd->write_size_bytes);

        lseek (p->b.f, fd->base_offset, SEEK_SET);
        ssize_t s = write (p->b.f, fd->buffer, fd->bytes_written);
        if (s != fd->bytes_written)
        {
            fprintf (stderr, "POSIX1 method tried to write %llu, "
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

void adios_posix1_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,void * data
                       ,struct adios_method_struct * method
                       )
{
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
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
            fprintf (stderr, "POSIX1 method tried to write %llu, "
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
            fprintf (stderr, "adios_posix1_write exceeds pg bound. File is corrupted. "
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
            fprintf (stderr, "POSIX1 method tried to write %llu, "
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

void adios_posix1_get_write_buffer (struct adios_file_struct * fd
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

void adios_posix1_read (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * buffer
                      ,uint64_t buffer_size
                      ,struct adios_method_struct * method
                      )
{
    v->data = buffer;
    v->data_size = buffer_size;
}

static void adios_posix1_do_write (struct adios_file_struct * fd
                                 ,struct adios_method_struct * method
                                 ,char * buffer
                                 ,uint64_t buffer_size
                                 )
{
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
                                                          method->method_data;
    int32_t to_write;
    uint64_t bytes_written = 0;

    if (fd->shared_buffer == adios_flag_yes)
    {
        lseek (p->b.f, p->b.end_of_pgs, SEEK_SET);
        if (p->b.end_of_pgs + fd->bytes_written > fd->pg_start_in_file + fd->write_size_bytes)
            fprintf (stderr, "adios_posix1_write exceeds pg bound. File is corrupted. "
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

static void adios_posix1_do_read (struct adios_file_struct * fd
                                ,struct adios_method_struct * method
                                )
{
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
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
            fprintf (stderr, "POSIX1 read: file version unknown: %u\n"
                    ,version
                    );
            return;
    }

    adios_buffer_struct_clear (&p->b);
}

void adios_posix1_close (struct adios_file_struct * fd
                       ,struct adios_method_struct * method
                       )
{
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
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
                    fprintf (stderr, "POSIX1 method tried to write %llu, "
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
                        fprintf (stderr, "adios_posix1_write exceeds pg bound. File is corrupted. "
                                         "Need to enlarge group size. \n");
                    ssize_t s = write (p->b.f, fd->buffer, fd->bytes_written);
                    if (s != fd->bytes_written)
                    {
                        fprintf (stderr, "POSIX1 method tried to write %llu, "
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
                    fprintf (stderr, "POSIX1 method tried to write %llu, "
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
            adios_build_index_v1 (fd, &new_pg_root, &new_vars_root
                                 ,&new_attrs_root
                                 );
            // if collective, gather the indexes from the rest and call
            // adios_merge_index_v1 (&new_pg_root, &new_vars_root, pg, vars);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, new_pg_root, new_vars_root
                                 ,new_attrs_root
                                 );
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_posix1_do_write (fd, method, buffer, buffer_offset);

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
                    fprintf (stderr, "POSIX1 method tried to write %llu, "
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
                        fprintf (stderr, "POSIX1 method tried to write %llu, "
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
                    fprintf (stderr, "POSIX1 method tried to write %llu, "
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
            adios_build_index_v1 (fd, &new_pg_root, &new_vars_root
                                 ,&new_attrs_root
                                 );
            // merge in old indicies
            adios_merge_index_v1 (&p->old_pg_root, &p->old_vars_root
                                 ,&p->old_attrs_root
                                 ,new_pg_root, new_vars_root, new_attrs_root
                                 );
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, p->old_pg_root, p->old_vars_root
                                 ,p->old_attrs_root
                                 );
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_posix1_do_write (fd, method, buffer, buffer_offset);

            free (buffer);

            break;
        }

        case adios_mode_read:
        {
            // read the index to find the place to start reading
            adios_posix1_do_read (fd, method);
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

void adios_posix1_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
    if (adios_posix1_initialized)
        adios_posix1_initialized = 0;
}

void adios_posix1_end_iteration (struct adios_method_struct * method)
{
}

void adios_posix1_start_calculation (struct adios_method_struct * method)
{
}

void adios_posix1_stop_calculation (struct adios_method_struct * method)
{
}
