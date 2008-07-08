#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// xml parser
#include <mxml.h>

#include "binpack-general.h"
#include "binpack-utils.h"
#include "br-utils.h"
#include "bw-utils.h"
#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"

static int adios_vtk_initialized = 0;

struct adios_File_data_struct
{
    long long handle;   // used for read
    int f;              // used for write
    void * buffer;
    unsigned long long buffer_size;
    unsigned long long start;
    int last_var_write_yes; // was the last item asked to write a write="yes"?
};

void adios_vtk_init (const char * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_File_data_struct * f;
    if (!adios_vtk_initialized)
    {
        adios_vtk_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_File_data_struct));
    f = (struct adios_File_data_struct *) method->method_data;
    f->f = -1;
    f->handle = 0;
    f->buffer = 0;
    f->buffer_size = 0;
    f->start = 0;
    f->last_var_write_yes = 0;
}

void adios_vtk_open (struct adios_file_struct * fd
                    ,struct adios_method_struct * method
                    )
{
    char name [STR_LEN];
    struct adios_File_data_struct * f = (struct adios_File_data_struct *)
                                                          method->method_data;

    sprintf (name, "%s%s", method->base_path, fd->name);
    if (fd->mode == adios_mode_read)
    {
        f->handle = br_fopen (fd->name);
        if (!f->handle)
        {
            fprintf (stderr, "file not found: %s\n", fd->name);

            return;
        }
    }
    else
    {
        if (fd->mode == adios_mode_append)
        {
            struct stat s;
            if (stat (name, &s) == 0)
                fd->base_offset += s.st_size;
            f->f = open (name, O_WRONLY);
            if (f->f == -1)
            {
                f->f = open (name,  O_WRONLY | O_CREAT
                            ,  S_IRUSR | S_IWUSR
                             | S_IRGRP | S_IWGRP
                             | S_IROTH | S_IWOTH
                            );
                if (f->f == -1)
                {
                    fprintf (stderr, "adios_vtk_open failed for "
                                     "base_path %s, name %s\n"
                            ,method->base_path, fd->name
                            );
                }
            }
        }
        else  // make sure we overwrite
        {
            f->f = open (name, O_WRONLY | O_CREAT | O_TRUNC
                        ,  S_IRUSR | S_IWUSR
                         | S_IRGRP | S_IWGRP
                         | S_IROTH | S_IWOTH
                        );
            if (f->f == -1)
            {
                fprintf (stderr, "adios_vtk_open failed for "
                                 "base_path %s, name %s\n"
                        ,method->base_path, fd->name
                        );
            }
        }
    }
}

void adios_vtk_write (struct adios_file_struct * fd
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

    if (v->copy_on_write == adios_flag_yes)
    {
        unsigned long long var_size;
        unsigned long long mem_allowed;

        var_size = adios_size_of_var (v, data);
        mem_allowed = adios_method_buffer_alloc (var_size);
        if (mem_allowed == var_size)
        {
            v->free_data = adios_flag_yes;
            v->data = adios_dupe_data (v, data);
            v->data_size = var_size;
        }
        else
        {
            adios_method_buffer_free (mem_allowed);
// need to handle overflow of allocatable space
            fprintf (stderr, "OVERFLOW in vtk when writing %s\n", v->name);
            v->free_data = adios_flag_no;
            v->data = 0;
            v->data_size = 0;
        }
    }
    else // otherwise, should be safe to hold pointer
    {
        v->data = data;
        v->free_data = adios_flag_no;
        v->data_size = 0;
    }
}

void adios_vtk_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,unsigned long long * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
    unsigned long long mem_allowed;
    
    if (*size == 0)
    {
        *size = adios_size_of_var (v, (void *) "");
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

// for now, just use the ones used for parse buffer.  Should they
// need to change, we're ready
static void vtk_read_pre_fetch (struct adios_bp_element_struct * element
                               ,void ** buffer, long long * buffer_size
                               ,void * private_data
                               )
{
    adios_pre_element_fetch (element, buffer, buffer_size, private_data);
}

// for now, just use the ones used for parse buffer.  Should they
// need to change, we're ready
static void vtk_read_post_fetch (struct adios_bp_element_struct * element
                                ,void * buffer, long long buffer_size
                                ,void * private_data
                                )
{
    adios_post_element_fetch (element, buffer, buffer_size, private_data);
}

void adios_vtk_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v
                    ,void * buffer
                    ,struct adios_method_struct * method
                    )
{
    v->data = buffer;
}

static void adios_vtk_do_write (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               )
{
    struct adios_File_data_struct * md = (struct adios_File_data_struct *)
                                         method->method_data;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_attribute_struct * a = fd->group->attributes;
    int overflow = 0;
    unsigned long long end = 0;

    unsigned long long size = adios_data_size (fd->group);
    unsigned long long mem_needed = 0;
    unsigned long long mem_allowed = 0;;

    if (size > md->buffer_size)
    {
        mem_needed = size - md->buffer_size;
        mem_allowed = adios_method_buffer_alloc (mem_needed);
        if (mem_needed != mem_allowed)
        {
            adios_method_buffer_free (mem_allowed);
            fprintf (stderr, "OVERFLOW!!\n");
            overflow = 1;
        }

        md->buffer_size = md->buffer_size + mem_allowed;
        free (md->buffer);
        md->buffer = malloc (md->buffer_size);
        if (!md->buffer)
            fprintf (stderr, "Cannot allocate buffer space\n");
    }

    while (v && !overflow)
    {
        overflow = adios_do_write_var (v
                                      ,md->buffer
                                      ,md->buffer_size
                                      ,md->start
                                      ,&end
                                      );

        md->start = end;
        v = v->next;
    }

    while (a && !overflow)
    {
        overflow = adios_do_write_attribute (a
                                            ,md->buffer
                                            ,md->buffer_size
                                            ,md->start
                                            ,&end
                                            );

        md->start = end;
        a = a->next;
    }

    if (overflow)
    {
        overflow = 0;
        long long file;
        char name [STR_LEN];

        sprintf (name, "%sOVERFLOW_%s", method->base_path, fd->name);
        bw_set_write (0); // set write to a file directly
        bw_fopen_ (name, &file);
        v = fd->group->vars;
        a = fd->group->attributes;

        while (v && !overflow)
        {
            overflow = adios_do_write_var (v
                                          ,&file
                                          ,9223372036854775807LL // LLONG_MAX
                                          ,0
                                          ,&end
                                          );

            v = v->next;
        }
        while (a && !overflow)
        {
            overflow = adios_do_write_attribute (a
                                                ,md->buffer
                                                ,md->buffer_size
                                                ,0
                                                ,&end
                                                );

            a = a->next;
        }

        bw_fclose_ (&file);
        bw_set_write (1); // set write to the buffer
    }
    else
    {
        lseek (md->f, fd->base_offset + fd->offset, SEEK_SET);
        write (md->f, md->buffer, end);
    }

    // clear out the cached data for the caller to be able to continue to work
    v = fd->group->vars;
    while (v)
    {
        if (v->free_data == adios_flag_yes)
        {
            if (v->data)
            {
                free (v->data);
                adios_method_buffer_free (v->data_size);
            }
        }

        v->data = 0;
        v->data_size = 0;
        v->free_data = adios_flag_no;
        v->got_buffer = adios_flag_no;
        v = v->next;
    }
}

static void adios_vtk_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
    struct adios_File_data_struct * md = (struct adios_File_data_struct *)
                                         method->method_data;
    struct adios_var_struct * v = fd->group->vars;
    
    unsigned long long element_size = 0;
    struct adios_bp_element_struct * element = NULL;
    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer = 0;
    data.buffer_len = 0;
    
    while ((element_size = br_get_next_element_specific (md->handle
                                                        ,vtk_read_pre_fetch
                                                        ,vtk_read_post_fetch
                                                        ,&data
                                                        ,&element
                                                        )
          ) != 0)
    {
        //printf ("element size: %d\n", element_size);
        //printf ("%s %s\n", adios_tag_to_string (element->tag), element->name);
        //printf ("\tPath: %s\n", element->path);

        br_free_element (element);
    }
    
    if (data.buffer)
        free (data.buffer);
}

void adios_vtk_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_File_data_struct * md = (struct adios_File_data_struct *)
                                                       method->method_data;

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
        adios_vtk_do_write (fd, method);

    if (fd->mode == adios_mode_read)
    {
        adios_vtk_do_read (fd, method);
        struct adios_var_struct * v = fd->group->vars;
        while (v)
        {
            v->data = 0;
            v = v->next;
        }
    }

    if (md->f != -1)
    {
        close (md->f);
    }
    else
    {
        if (md->handle)
        {
            br_fclose (md->handle);
        }
    }

    md->f = -1;
    md->handle = 0;
    md->start = 0;
    md->last_var_write_yes = 0;
}

void adios_vtk_finalize (int mype, struct adios_method_struct * method)
{
/* nothing to do here */
    if (adios_vtk_initialized)
        adios_vtk_initialized = 0;
}

void adios_vtk_end_iteration (struct adios_method_struct * method)
{
}

void adios_vtk_start_calculation (struct adios_method_struct * method)
{
}

void adios_vtk_stop_calculation (struct adios_method_struct * method)
{
}
