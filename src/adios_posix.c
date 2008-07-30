#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

// xml parser
#include <mxml.h>

#include "binpack-general.h"
#include "binpack-utils.h"
#include "br-utils.h"
#include "bw-utils.h"
#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

static int adios_posix_initialized = 0;

struct adios_POSIX_data_struct
{
    struct adios_bp_buffer_struct_v1 b;
    uint64_t base_offset;   // where appends should start
};

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
    p->base_offset = 0;
}

int adios_posix_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    char name [STR_LEN];
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;

    sprintf (name, "%s%s", method->base_path, fd->name);
    struct stat s;
    if (stat (name, &s) == 0)
        p->b.file_size = s.st_size;
    if (fd->mode == adios_mode_read)
    {
        p->b.f = open64 (name, O_RDONLY);
        if (p->b.f == -1)
        {
            fprintf (stderr, "ADIOS POSIX: file not found: %s\n", fd->name);

            return 0;
        }
    }
    else
    {
        if (fd->mode == adios_mode_append)
        {
            p->base_offset = p->b.file_size;
            p->b.f = open (name, O_WRONLY);
            if (p->b.f == -1)
            {
                p->b.f = open64 (name,  O_WRONLY | O_CREAT
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

                    return 0;
                }
            }
        }
        else  // make sure we overwrite
        {
            p->b.f = open64 (name, O_WRONLY | O_CREAT | O_TRUNC
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

                return 0;
            }
        }
    }

    return 1;
}

int adios_posix_should_buffer (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              ,void * comm
                              )
{
    return 1;   // as far as we care, buffer
}

void adios_posix_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,void * data
                       ,struct adios_method_struct * method
                       )
{
    if (v->got_buffer)
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

    // nothing to do since we will use the shared buffer
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
fprintf (stderr, "Need to figure out the size of the var\n");
        //*size = adios_size_of_var (v, (void *) "");
        *size = 10000;
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
                      ,struct adios_method_struct * method
                      )
{
    v->data = buffer;
}

static void adios_posix_do_write (struct adios_file_struct * fd
                                 ,char * buffer
                                 ,uint64_t buffer_size
                                 ,struct adios_method_struct * method
                                 )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    write (p->b.f, fd->buffer, fd->bytes_written);
    write (p->b.f, buffer, buffer_size);
}

static void adios_posix_do_read (struct adios_file_struct * fd
                                ,struct adios_method_struct * method
                                )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;
    struct adios_var_struct * v = fd->group->vars;
    
    uint64_t element_size = 0;
    struct adios_bp_element_struct * element = NULL;
    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer = 0;
    data.buffer_len = 0;

    uint32_t version = 0;

    adios_posix_read_version (&p->b);
    adios_parse_version (&p->b, &version);

    switch (version)
    {
        case 1:
        {
            struct adios_index_process_group_struct_v1 * pg_root = 0;
            struct adios_index_var_struct_v1 * vars_root = 0;

            adios_posix_read_index_offsets (&p->b);
            adios_parse_index_offsets_v1 (&p->b);

            adios_posix_read_process_group_index (&p->b);
            adios_parse_process_group_index_v1 (&p->b, &pg_root);
#if 1
            adios_posix_read_vars_index (&p->b);
            adios_parse_vars_index_v1 (&p->b, &vars_root);
#endif
fprintf (stderr, "need to read from the last one?\n");

            // the three section headers
            struct adios_process_group_header_struct_v1 pg_header;
            struct adios_vars_header_struct_v1 vars_header;
            struct adios_attributes_header_struct_v1 attrs_header;

            struct adios_var_header_struct_v1 var_header;
            struct adios_var_payload_struct_v1 var_payload;
            struct adios_attribute_struct_v1 attribute;

            int i;

            adios_posix_read_process_group (&p->b, 1, pg_root);
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
                                                    );
                }
                else
                {
                    adios_parse_var_data_payload_v1 (&p->b, &var_header, NULL);
                }
            }

#if 1
            adios_parse_attributes_header_v1 (&p->b, &attrs_header);

            for (i = 0; i < attrs_header.count; i++)
            {
                adios_parse_attribute_v1 (&p->b, &attribute);
            }
#endif
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

    struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_var_struct_v1 * vars_root = 0;

    if (fd->mode == adios_mode_write)
    {
        char * buffer = 0;
        uint64_t buffer_size = 0;
        uint64_t buffer_offset = 0;
        uint64_t index_start = fd->offset;

        // build index
        adios_build_index_v1 (fd, &pg_root, &vars_root);
        // if collective, gather the indexes from the rest and call
        // adios_merge_index_v1 (&pg_root, &vars_root, new_pg, new_vars);
        adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                             ,index_start, pg_root, vars_root
                             );
        adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
        adios_posix_do_write (fd, buffer, buffer_offset, method);

        free (buffer);

        adios_clear_index_v1 (pg_root, vars_root);
    }

    if (fd->mode == adios_mode_append)
    {
        char * buffer = 0;
        uint64_t buffer_size = 0;
        uint64_t buffer_offset = 0;
        uint64_t index_start = fd->offset;

        // build index
        adios_build_index_v1 (fd, &pg_root, &vars_root);
        // merge in with old struct
        // merge in peers with old struct from initial read
        adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                             ,index_start, pg_root, vars_root
                             );
        adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
        adios_posix_do_write (fd, buffer, buffer_offset, method);

        free (buffer);

        adios_clear_index_v1 (pg_root, vars_root);
    }

    if (fd->mode == adios_mode_read)
    {
        // read the index to find the place to start reading
        adios_posix_do_read (fd, method);
        struct adios_var_struct * v = fd->group->vars;
        while (v)
        {
            v->data = 0;
            v = v->next;
        }
    }

    adios_posix_close_internal (&p->b);
    p->base_offset = 0;
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
