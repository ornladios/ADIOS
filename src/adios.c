/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include "../config.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

// xml parser
#include <mxml.h>

#include "adios.h"
#include "common_adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "adios_internals_mxml.h"
#include "globals.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

extern struct adios_transport_struct * adios_transports;

int adios_set_application_id (int id)
{
    globals_adios_set_application_id (id);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int adios_init (const char * config)
{
    return common_adios_init (config);
}


///////////////////////////////////////////////////////////////////////////////
// all XML file pieces will be provided by another series of calls
int adios_init_noxml ()
{
    return common_adios_init_noxml ();
}

///////////////////////////////////////////////////////////////////////////////
int adios_finalize (int mype)
{
    return common_adios_finalize (mype);
}

///////////////////////////////////////////////////////////////////////////////
int adios_allocate_buffer (enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when
                          ,uint64_t buffer_size)
{
    return common_adios_allocate_buffer (adios_buffer_alloc_when, buffer_size);
}

///////////////////////////////////////////////////////////////////////////////
int adios_open (int64_t * fd, const char * group_name, const char * name
               ,const char * mode, void * comm
               )
{
    return common_adios_open (fd, group_name, name, mode, comm);
}

///////////////////////////////////////////////////////////////////////////////
int adios_group_size (int64_t fd_p, uint64_t data_size
                     ,uint64_t * total_size
                     )
{
    return common_adios_group_size (fd_p, data_size, total_size);
}

///////////////////////////////////////////////////////////////////////////////
/* This C api function is a bit different from the Fortran api funcion, but
 * they call the same common_adios_write()
 */
int adios_write (int64_t fd_p, const char * name, void * var)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;
    if (!fd)
    {
        fprintf (stderr, "Invalid handle passed to adios_write\n");
        return 1;
    }

    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL)
    {
        // nothing to do so just return
        return 0;
    }

    v = adios_find_var_by_name (v, name, fd->group->all_unique_var_names);

    if (!v)
    {
        fprintf (stderr, "Bad var name (ignored): '%s'\n", name);

        return 1;
    }

    if (fd->mode == adios_mode_read)
    {
        if (   strcasecmp (name, fd->group->group_comm)
            && v->is_dim != adios_flag_yes
           )
        {
            fprintf (stderr, "write attempted on %s in %s.  This was opened for read\n" ,name , fd->name);
            return 1;
        }
    }

    if (!var)
    {
        fprintf (stderr, "Invalid data: %s\n", name);

        return 1;
    }

    if (v->data)
    {
        free (v->data);
        v->data = 0;
    }

    // Q.L. 10-2010. To fix a memory leak problem.
    if (v->stats)
    {   
        int j, idx;
        int c, count = 1;

        if (v->type == adios_complex || v->type == adios_double_complex)
            count = 3;

        for (c = 0; c < count; c ++)
        {   
            j = idx = 0;
            while (v->bitmap >> j)
            {   
                if (v->bitmap >> j & 1)
                {   
                    if (j == adios_statistic_hist)
                    {   
                        struct adios_index_characteristics_hist_struct * hist =
                            (struct adios_index_characteristics_hist_struct *) v->stats[c][idx].data;
                        if (hist)
                        {   
                            free (hist->breaks);
                            free (hist->frequencies);
                            free (hist);
                            v->stats[c][idx].data = 0;
                        }
                    }
                    else
                    {
                        if (v->stats[c][idx].data)
                        {
                            free (v->stats[c][idx].data);
                            v->stats[c][idx].data = 0;
                        }
                    }
                    idx ++;
                }
                j ++;
            }
        }
    }

    if (v->dimensions)
    {
        v->data = var;
    }
    else
    {
        uint64_t element_size = adios_get_type_size (v->type, var);

        switch (v->type)
        {
            case adios_byte:
            case adios_short:
            case adios_integer:
            case adios_long:
            case adios_unsigned_byte:
            case adios_unsigned_short:
            case adios_unsigned_integer:
            case adios_unsigned_long:
            case adios_real:
            case adios_double:
            case adios_long_double:
            case adios_complex:
            case adios_double_complex:
                v->data = malloc (element_size);
                if (!v->data)
                {
                    fprintf (stderr, "cannot allocate %lld bytes to copy "
                                     "scalar %s\n"
                            ,element_size
                            ,v->name
                            );

                    return 0;
                }

                memcpy ((char *) v->data, var, element_size);
                break;

            case adios_string:
                v->data = malloc (element_size + 1);
                if (!v->data)
                {
                    fprintf (stderr, "cannot allocate %lld bytes to copy "
                                     "scalar %s\n"
                            ,element_size
                            ,v->name
                            );

                    return 0;
                }
                ((char *) v->data) [element_size] = 0;
                memcpy ((char *) v->data, var, element_size);
                break;

            default:
                v->data = 0;
                break;
        }
    }

    common_adios_write (fd, v, var);
    // v->data is set to NULL in the above call

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append) 
    {
        adios_copy_var_written (&fd->group->vars_written, v, fd);
    }

    return 0;
}


///////////////////////////////////////////////////////////////////////////////
int adios_get_write_buffer (int64_t fd_p, const char * name
                           ,uint64_t * size
                           ,void ** buffer
                           )
{
    return common_adios_get_write_buffer (fd_p, name, size, buffer);
}

///////////////////////////////////////////////////////////////////////////////
int adios_read (int64_t fd_p, const char * name, void * buffer
               ,uint64_t buffer_size
               )
{
    return common_adios_read (fd_p, name, buffer, buffer_size);
}

///////////////////////////////////////////////////////////////////////////////
int adios_set_path (int64_t fd_p, const char * path)
{
    return common_adios_set_path (fd_p, path);
}

///////////////////////////////////////////////////////////////////////////////
int adios_set_path_var (int64_t fd_p, const char * path, const char * name)
{
    return common_adios_set_path_var (fd_p, path, name);
}

///////////////////////////////////////////////////////////////////////////////
// hint that we reached the end of an iteration (for asynchronous pacing)
int adios_end_iteration ()
{
    return common_adios_end_iteration ();
}

///////////////////////////////////////////////////////////////////////////////
// hint to start communicating
int adios_start_calculation ()
{
    return common_adios_start_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
// hint to stop communicating
int adios_stop_calculation ()
{
    return common_adios_stop_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
int adios_close (int64_t fd_p)
{
    return common_adios_close (fd_p);
}

//////////////////////////////////////////////////////////////////////////////
// Methods normally only called by the XML parser
//////////////////////////////////////////////////////////////////////////////

// adios_common_declare_group is in adios_internals.c

///////////////////////////////////////////////////////////////////////////////
// group a list of vars into a composite group
int adios_declare_group (int64_t * id, const char * name
                        ,const char * time_index
                        ,enum ADIOS_FLAG stats
                        )
{
    int ret;
    ret = adios_common_declare_group (id, name, adios_flag_no
                                      ,""
                                      ,""
                                      ,time_index
                                      ,stats
                                      );
    if (ret == 1) {
        struct adios_group_struct * g = (struct adios_group_struct *) *id;
        g->all_unique_var_names = adios_flag_no;
    }
    return ret;
}


int adios_free_group (int64_t id)
{
    return adios_common_free_group (id);
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_var is in adios_internals.c

// declare a single var as an entry in a group
int adios_define_var (int64_t group_id, const char * name
                     ,const char * path, int type
                     ,const char * dimensions
                     ,const char * global_dimensions
                     ,const char * local_offsets
                     )
{
    return adios_common_define_var (group_id, name, path
                                   ,(enum ADIOS_DATATYPES) type
                                   ,dimensions
                                   ,global_dimensions, local_offsets
                                   );
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_define_attribute is in adios_internals.c

int adios_define_attribute (int64_t group, const char * name
                           ,const char * path, enum ADIOS_DATATYPES type
                           ,const char * value, const char * var
                           )
{
    return adios_common_define_attribute (group, name, path, type, value, var);
}

///////////////////////////////////////////////////////////////////////////////

// adios_common_select_method is in adios_internals_mxml.c
int adios_select_method (int64_t group, const char * method
                        ,const char * parameters
                        ,const char * base_path
                        )
{
    return adios_common_select_method_by_group_id (0, method, parameters, group
                                                  ,base_path, 0
                                                  );
}

