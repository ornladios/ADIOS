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

#include "public/adios_mpi.h" // MPI or dummy MPI
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"
#include "core/adios_logger.h"

static int adios_posix1_initialized = 0;

struct adios_POSIX1_data_struct
{
    // our file bits
    struct adios_bp_buffer_struct_v1 b;

    // old index structs we read in and have to be merged in
    struct adios_index_struct_v1 * index;
};

void adios_posix1_init (const PairStruct * parameters
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
    p->index = adios_alloc_index_v1(1); // with hashtables
}

int adios_posix1_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method, MPI_Comm comm
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
#ifndef __CYGWIN__
| O_LARGEFILE
#endif
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
#ifndef __CYGWIN__
| O_LARGEFILE
#endif
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
#ifndef __CYGWIN__
| O_LARGEFILE
#endif
#endif			   
);
            if (p->b.f == -1)
            {
                old_file = 0;
                p->b.f = open (name,  O_WRONLY | O_CREAT
#ifndef __APPLE__
#ifndef __CYGWIN__
| O_LARGEFILE
#endif
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
                        adios_parse_attributes_index_v1 (&p->b
                                                        ,&p->index->attrs_root
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

enum BUFFERING_STRATEGY adios_posix1_should_buffer (struct adios_file_struct * fd
                                                   ,struct adios_method_struct * method
                                                   )
{
    return continue_with_new_pg;
}

void adios_posix1_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,const void * data
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
    v->data = v->adata = buffer;
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
        if (p->b.end_of_pgs + fd->bytes_written > fd->pg_start_in_file + fd->offset)
            fprintf (stderr, "adios_posix1_write exceeds pg bound. File is corrupted. "
                             "Need to enlarge group size. \n");

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
            fprintf (stderr, "POSIX1 read: file version unknown: %u\n", version);
            return;
    }

    adios_buffer_struct_clear (&p->b);
}


void adios_posix1_buffer_overflow (struct adios_file_struct * fd, 
                                   struct adios_method_struct * method)
{
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
                                                 method->method_data;
    log_error ("POSIX1 method only works with complete buffering of data between adios_open() "
               "and adios_close(). Variables that do not fit into the buffer will not be "
               "written by this method to file %s\n", fd->name);
}

void adios_posix1_close (struct adios_file_struct * fd
                       ,struct adios_method_struct * method
                       )
{
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
                                                          method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    switch (fd->mode)
    {
        case adios_mode_write:
        {

            // buffering or not, write the index
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = fd->base_offset + fd->offset;

            // build index
            adios_build_index_v1 (fd, p->index);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, p->index);
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_posix1_do_write (fd, method, buffer, buffer_offset);

            free (buffer);

            break;
        }

        case adios_mode_append:
        case adios_mode_update:
        {

            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = fd->base_offset + fd->offset;

            // build index. Note: It merges with old indices already stored 
            //                    in p->index in adios_posix_open's append case
            adios_build_index_v1 (fd, p->index);
            // merge in old indicies
            //adios_merge_index_v1 (p->index, new_pg_root, 
            //                      new_vars_root, new_attrs_root
            //                     );
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, p->index);
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

    adios_posix_close_internal (&p->b);
    adios_clear_index_v1 (p->index);
}

void adios_posix1_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_POSIX1_data_struct * p = (struct adios_POSIX1_data_struct *)
        method->method_data;
    adios_free_index_v1 (p->index);
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
