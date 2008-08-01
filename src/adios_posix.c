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
    // our file bits
    struct adios_bp_buffer_struct_v1 b;

    // old index structs we read in and have to be merged in
    struct adios_index_process_group_struct_v1 * old_pg_root;
    struct adios_index_var_struct_v1 * old_vars_root;
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
    p->old_pg_root = 0;
    p->old_vars_root = 0;
}

int adios_posix_open (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    char * name;
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;

    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);
    struct stat s;
    if (stat (name, &s) == 0)
        p->b.file_size = s.st_size;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            p->b.f = open64 (name, O_RDONLY);
            if (p->b.f == -1)
            {
                fprintf (stderr, "ADIOS POSIX: file not found: %s\n", fd->name);

                free (name);

                return 0;
            }
            fd->base_offset = 0;

            break;
        }

        case adios_mode_write:
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

                free (name);

                return 0;
            }
            fd->base_offset = 0;

            break;
        }

        case adios_mode_append:
        {
            int old_file = 1;
            p->b.f = open (name, O_RDWR);
            if (p->b.f == -1)
            {
                old_file = 0;
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

                switch (version)
                {
                    case 1:
                        // read the old stuff and set the base offset
                        adios_posix_read_index_offsets (&p->b);
                        adios_parse_index_offsets_v1 (&p->b);

                        adios_posix_read_process_group_index (&p->b);
                        adios_parse_process_group_index_v1 (&p->b
                                                           ,&p->old_pg_root
                                                           );

                        adios_posix_read_vars_index (&p->b);
                        adios_parse_vars_index_v1 (&p->b, &p->old_vars_root);

                        fd->base_offset = p->b.end_of_pgs;
                        break;

                    default:
                        fprintf (stderr, "Unkown bp version.  Cannot append\n");

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
                      ,struct adios_method_struct * method
                      )
{
    v->data = buffer;
}

static void adios_posix_do_write (struct adios_file_struct * fd
                                 ,struct adios_method_struct * method
                                 ,char * buffer
                                 ,uint64_t buffer_size
                                 )
{
    struct adios_POSIX_data_struct * p = (struct adios_POSIX_data_struct *)
                                                          method->method_data;

    lseek64 (p->b.f, fd->base_offset, SEEK_SET);
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
            struct adios_index_process_group_struct_v1 * pg_root_temp = 0;
            struct adios_index_var_struct_v1 * vars_root = 0;

            adios_posix_read_index_offsets (&p->b);
            adios_parse_index_offsets_v1 (&p->b);

            adios_posix_read_process_group_index (&p->b);
            adios_parse_process_group_index_v1 (&p->b, &pg_root);
#if 1
            adios_posix_read_vars_index (&p->b);
            adios_parse_vars_index_v1 (&p->b, &vars_root);
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
            while (pg_root && pg_root_temp->next)
                pg_root = pg_root->next;

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
            adios_clear_index_v1 (pg_root, vars_root);
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

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;

    switch (fd->mode)
    {
        case adios_mode_write:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = fd->base_offset + fd->offset;

            // build index
            adios_build_index_v1 (fd, &new_pg_root, &new_vars_root);
            // if collective, gather the indexes from the rest and call
            // adios_merge_index_v1 (&new_pg_root, &new_vars_root, pg, vars);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, new_pg_root, new_vars_root
                                 );
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_posix_do_write (fd, method, buffer, buffer_offset);

            free (buffer);

            adios_clear_index_v1 (new_pg_root, new_vars_root);

            break;
        }

        case adios_mode_append:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = fd->base_offset + fd->offset;

            // build index
            adios_build_index_v1 (fd, &new_pg_root, &new_vars_root);
            // merge in old indicies
            adios_merge_index_v1 (&p->old_pg_root, &p->old_vars_root
                                 ,new_pg_root, new_vars_root);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, p->old_pg_root, p->old_vars_root
                                 );
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
    adios_clear_index_v1 (p->old_pg_root, p->old_vars_root);
    p->old_pg_root = 0;
    p->old_vars_root = 0;
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
