/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

// xml parser
#include <mxml.h>

#include "common_adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "adios_internals_mxml.h"
#include "buffer.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

extern struct adios_transport_struct * adios_transports;

///////////////////////////////////////////////////////////////////////////////
int common_adios_init (const char * config)
{
    // parse the config file
    return adios_parse_config (config);
}

///////////////////////////////////////////////////////////////////////////////
// all XML file pieces will be provided by another series of calls
int common_adios_init_noxml ()
{
    return adios_local_config ();
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_finalize (int mype)
{
    struct adios_method_list_struct * m;

    for (m = adios_get_methods (); m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_finalize_fn
           )
        {
            adios_transports [m->method->m].adios_finalize_fn (mype, m->method);
        }
    }

    adios_cleanup ();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_allocate_buffer (enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when
                                 ,uint64_t buffer_size)
{
    adios_buffer_size_requested_set (buffer_size * 1024 * 1024);
    adios_buffer_alloc_when_set (adios_buffer_alloc_when);

    return adios_set_buffer_size ();
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_open (int64_t * fd, const char * group_name
                ,const char * name, const char * file_mode, void * comm
               )
{
    int64_t group_id = 0;
    struct adios_file_struct * fd_p = (struct adios_file_struct *)
                                  malloc (sizeof (struct adios_file_struct));
    struct adios_group_struct * g = 0;
    struct adios_method_list_struct * methods = 0;
    enum ADIOS_METHOD_MODE mode;

    adios_common_get_group (&group_id, group_name);
    g = (struct adios_group_struct *) group_id;
    methods = g->methods;

    if (!strcasecmp (file_mode, "r"))
        mode = adios_mode_read;
    else
        if (!strcasecmp (file_mode, "w"))
            mode = adios_mode_write;
        else
            if (!strcasecmp (file_mode, "a"))
                mode = adios_mode_append;
            else
                if (!strcasecmp (file_mode, "u"))
                    mode = adios_mode_update;
                else
                {
                    fprintf (stderr, "adios_open: unknown file mode: %s\n"
                            ,file_mode
                            );

                    *fd = 0;

                    return 1;
                }

    fd_p->name = strdup (name);
    fd_p->subfile_index = -1; // subfile index is by default -1
    fd_p->group = g;
    fd_p->mode = mode;
    fd_p->data_size = 0;
    fd_p->buffer = 0;
    fd_p->offset = 0;
    fd_p->bytes_written = 0;
    fd_p->buffer_size = 0;
    fd_p->vars_start = 0;
    fd_p->vars_written = 0;
    fd_p->write_size_bytes = 0;
    fd_p->base_offset = 0;
    fd_p->pg_start_in_file = 0;

    if (mode != adios_mode_read)
        g->time_index++;

    while (methods)
    {
        if (   methods->method->m != ADIOS_METHOD_UNKNOWN
            && methods->method->m != ADIOS_METHOD_NULL
            && adios_transports [methods->method->m].adios_open_fn
           )
        {
            adios_transports [methods->method->m].adios_open_fn
                                                 (fd_p, methods->method, comm);
        }

        methods = methods->next;
    }

    *fd = (int64_t) fd_p;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_group_size (int64_t fd_p
                     ,uint64_t data_size
                     ,uint64_t * total_size
                     )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        fprintf (stderr, "Invalid handle passed to adios_group_size\n");

        return 1;
    }
    struct adios_method_list_struct * m = fd->group->methods;
    if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL)
    {
        // nothing to do so just return
        fd->shared_buffer = adios_flag_no;
        fd->write_size_bytes = 0;
        fd->buffer = 0;
        *total_size = 0;
        return 0;
    }

    fd->write_size_bytes = data_size;

    uint64_t overhead = adios_calc_overhead_v1 (fd);

    *total_size = data_size + overhead;

    // try to reserve a buffer using the adios_method_buffer_alloc
    // if it does not give the correct amount, overflow.  Make sure
    // the amount given is big enough for the first part of the file.

    fd->write_size_bytes += overhead;

    uint64_t allocated = adios_method_buffer_alloc (fd->write_size_bytes);
    if (allocated != fd->write_size_bytes)
    {
        fd->shared_buffer = adios_flag_no;

        fprintf (stderr, "adios_group_size (%s): Not buffering. "
                         "needs: %llu available: %llu.\n"
                ,fd->group->name, fd->write_size_bytes, allocated
                );
    }
    else
    {
        fd->shared_buffer = adios_flag_yes;
    }

    // call each transport method to coordinate the write and handle
    // if an overflow is detected.
    // now tell each transport attached that it is being written
    while (m)
    {
        enum ADIOS_FLAG should_buffer = adios_flag_yes;
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_should_buffer_fn
           )
        {
            should_buffer = adios_transports [m->method->m].
                                            adios_should_buffer_fn (fd
                                                                   ,m->method
                                                                   );
        }

        if (should_buffer == adios_flag_no)     // can't write directly since
            fd->shared_buffer = adios_flag_no;  // some might want to share

        m = m->next;
    }

    if (fd->shared_buffer == adios_flag_no)
    {
        adios_method_buffer_free (allocated);
        fd->buffer = 0;
        fd->offset = 0;
        fd->bytes_written = 0;
    }
    else
    {
        fd->buffer = malloc (fd->write_size_bytes);
        fd->buffer_size = fd->write_size_bytes;
        fd->offset = 0;
        fd->bytes_written = 0;
        if (!fd->buffer)
        {
            fprintf (stderr, "Cannot allocate %llu bytes for buffered output.\n",
                    fd->write_size_bytes);

            return 1;
        }
        else
        {
            // write the process group header
            adios_write_process_group_header_v1 (fd, *total_size);

            // setup for writing vars
            adios_write_open_vars_v1 (fd);
        }
    }

    // each var will be added to the buffer by the adios_write calls
    // attributes will be added by adios_close

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
/* common_adios_write is just a partial implementation. It expects filled out
 * structures. This is because C and Fortran implementations of adios_write are 
 * different for some part and this is the common part.
 */
int common_adios_write (struct adios_file_struct * fd, struct adios_var_struct * v, void * var)
{
    struct adios_method_list_struct * m = fd->group->methods;

    if (fd->shared_buffer == adios_flag_yes)
    {
        // var payload sent for sizing information
        adios_write_var_header_v1 (fd, v);

        // write payload
        adios_write_var_payload_v1 (fd, v);
    }

    // now tell each transport attached that it is being written
    while (m)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_write_fn
           )
        {
            adios_transports [m->method->m].adios_write_fn
                                   (fd, v, var, m->method);
        }

        m = m->next;
    }

    if (v->dimensions)
    {
        v->data = 0;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_get_write_buffer (int64_t fd_p, const char * name
                           ,uint64_t * size
                           ,void ** buffer
                           )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        fprintf (stderr, "Invalid handle passed to adios_get_write_buffer\n");

        return 1;
    }
    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    v = adios_find_var_by_name (v, name, fd->group->all_unique_var_names);

    if (!v)
    {
        fprintf (stderr
                ,"Bad var name (ignored): '%s' (%c%c%c)\n"
                ,name, name[0], name[1], name[2]
                );

        return 1;
    }

    if (fd->mode == adios_mode_read)
    {
        fprintf (stderr, "write attempted on %s in %s.  This was opened for"
                         " read\n"
                ,name , fd->name
                );

        return 1;
    }

    // since we are only getting one buffer, get it from the first
    // transport method that can provide it.
    while (m)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_get_write_buffer_fn
           )
        {
            adios_transports [m->method->m].adios_get_write_buffer_fn
                                (fd, v, size, buffer, m->method);
            m = 0;
        }
        else
            m = m->next;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_read (int64_t fd_p, const char * name, void * buffer
               ,uint64_t buffer_size
               )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        fprintf (stderr, "Invalid handle passed to adios_read\n");

        return 1;
    }
    struct adios_var_struct * v;
    struct adios_method_list_struct * m = fd->group->methods;

    if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL)
    {
        // nothing to do so just return
        return 0;
    }

    if (!(fd->mode == adios_mode_read))
    {
        fprintf (stderr, "read attempted on %s which was opened for write\n"
                ,fd->name
                );

        return 1;
    }

    v = adios_find_var_by_name (fd->group->vars, name
                               ,fd->group->all_unique_var_names
                               );
    if (v)
    {
        // since can only read from one place into the buffer,
        // read from the first transport method that can
        while (m)
        {
            if (   m->method->m != ADIOS_METHOD_UNKNOWN
                && m->method->m != ADIOS_METHOD_NULL
                && adios_transports [m->method->m].adios_read_fn
               )
            {
                adios_transports [m->method->m].adios_read_fn
                                     (fd, v, buffer, buffer_size, m->method);
                m = 0;
            }
            else
                m = m->next;
	}
    }
    else
    {
        fprintf (stderr, "var %s in file %s not found on read\n"
                ,name, fd->name
                );

        return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_set_path (int64_t fd_p, const char * path)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        fprintf (stderr, "Invalid handle passed to adios_set_path\n");

        return 1;
    }
    struct adios_group_struct * t = fd->group;
    struct adios_var_struct * v = t->vars;
    struct adios_attribute_struct * a = t->attributes;

    while (v)
    {
        if (v->path)
        {
            free (v->path);
        }

        v->path = strdup (path);

        v = v->next;
    }

    while (a)
    {
        if (a->path)
        {
            free (a->path);
        }

        a->path = strdup (path);

        a = a->next;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_set_path_var (int64_t fd_p, const char * path
                       ,const char * name
                       )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        fprintf (stderr, "Invalid handle passed to adios_set_path_var\n");

        return 1;
    }
    struct adios_group_struct * t = fd->group;
    struct adios_var_struct * v = t->vars;

    // check for vars and then attributes
    v = adios_find_var_by_name (t->vars, name, fd->group->all_unique_var_names);

    if (v)
    {
        if (v->path)
        {
            free (v->path);
        }

        v->path = strdup (path);
    }
    else
    {
        fprintf (stderr, "adios_set_path_var (path=%s, var=%s): var not found\n"
                ,path, name
                );

        return 1;
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// hint that we reached the end of an iteration (for asynchronous pacing)
int common_adios_end_iteration ()
{
    struct adios_method_list_struct * m;

    for (m = adios_get_methods (); m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_end_iteration_fn
           )
        {
            adios_transports [m->method->m].adios_end_iteration_fn
                                                (m->method);
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// hint to start communicating
int common_adios_start_calculation ()
{
    struct adios_method_list_struct * m;

    for (m = adios_get_methods (); m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_start_calculation_fn
           )
        {
            adios_transports [m->method->m].adios_start_calculation_fn
                                                  (m->method);
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// hint to stop communicating
int common_adios_stop_calculation ()
{
    struct adios_method_list_struct * m;

    for (m = adios_get_methods (); m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_stop_calculation_fn
           )
        {
            adios_transports [m->method->m].adios_stop_calculation_fn
                                                   (m->method);
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_close (int64_t fd_p)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        fprintf (stderr, "Invalid handle passed to adios_close\n");

        return 1;
    }
    struct adios_method_list_struct * m = fd->group->methods;
    if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL)
    {
        // nothing to do so just return
        return 0;
    }

    struct adios_attribute_struct * a = fd->group->attributes;
    struct adios_var_struct * v = fd->group->vars;

    if (fd->shared_buffer == adios_flag_yes)
    {
        adios_write_close_vars_v1 (fd);

        adios_write_open_attributes_v1 (fd);

        while (a)
        {
            adios_write_attribute_v1 (fd, a);

            a = a->next;
        }

        adios_write_close_attributes_v1 (fd);
    }

    // in order to get the index assembled, we need to do it in the
    // transport once we have collected all of the pieces

    // now tell all of the transports to write the buffer during close
    for (;m; m = m->next)
    {
        if (   m->method->m != ADIOS_METHOD_UNKNOWN
            && m->method->m != ADIOS_METHOD_NULL
            && adios_transports [m->method->m].adios_close_fn
           )
        {
            adios_transports [m->method->m].adios_close_fn
                                 (fd, m->method);
        }
    }

    if (fd->shared_buffer == adios_flag_yes)
    {
        adios_method_buffer_free (fd->write_size_bytes);
        free (fd->buffer);
        fd->buffer_size = 0;
        fd->buffer = 0;
        fd->offset = 0;
    }

    while (v)
    {
        v->write_offset = 0;
        if (v->data)
        {
            free (v->data);
            v->data = 0;
        }

        v = v->next;
    }

    while (fd->group->vars_written)
    {
        if (fd->group->vars_written->name)
            free (fd->group->vars_written->name);
        if (fd->group->vars_written->path)
            free (fd->group->vars_written->path);

        while (fd->group->vars_written->dimensions)
        {
            struct adios_dimension_struct * dimensions
                            = fd->group->vars_written->dimensions->next;

            free (fd->group->vars_written->dimensions);
            fd->group->vars_written->dimensions = dimensions;
        }

		// NCSU - Clear stat
        if (fd->group->vars_written->stats)
		{
            uint8_t j = 0, idx = 0;
            uint8_t c = 0, count = adios_get_stat_set_count(fd->group->vars_written->type);

            for (c = 0; c < count; c ++)
            {
                while (fd->group->vars_written->bitmap >> j)
                {
                    if ((fd->group->vars_written->bitmap >> j) & 1)
                    {
                        if (j == adios_statistic_hist)
                        {
                            struct adios_hist_struct * hist = (struct adios_hist_struct *) (fd->group->vars_written->stats[c][idx].data);
                            free (hist->breaks);
                            free (hist->frequencies);
                            free (hist);
                        }
                        else
                            free (fd->group->vars_written->stats[c][idx].data);

                        idx ++;
                    }
                    j ++;
                }
                free (fd->group->vars_written->stats[c]);
            }
            free (fd->group->vars_written->stats);
		}
        if (fd->group->vars_written->data)
            free (fd->group->vars_written->data);

        v = fd->group->vars_written->next;
        free (fd->group->vars_written);
        fd->group->vars_written = v;
    }

    if (fd->name)
    {
        free (fd->name);
        fd->name = 0;
    }

    free ((void *) fd_p);

    return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Methods normally only called by the XML parser
//////////////////////////////////////////////////////////////////////////////

// adios_common_declare_group is in adios_internals.c
// adios_common_define_var is in adios_internals.c
// adios_common_define_attribute is in adios_internals.c
// adios_common_select_method is in adios_internals.c
