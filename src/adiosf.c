/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include "../config.h"
#include <string.h>
#include <unistd.h>
#include <stdint.h>

#include "common_adios.h"
#include "adios_internals.h"
#include "adios_internals_mxml.h"
#include "adios_transport_hooks.h"
#include "futils.h"
#include "globals.h"

#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif

extern int adios_errno;

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_set_application_id, adios_SET_APPLICATION_ID) (int *id, int * err)
{
    globals_adios_set_application_id (*id);
    if (err != 0) *err = 0;
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_init, adios_INIT) (const char * config, int * err, int config_size)
{
    char * buf1 = 0;

    buf1 = futils_fstr_to_cstr (config, config_size);
    if (buf1 != 0) {
        *err = common_adios_init (buf1);
        free (buf1);
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_init_noxml, adios_INIT_LOCAL) (int * err)
{
    *err = common_adios_init_noxml ();
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_finalize, adios_FINALIZE) (int * mype, int * err)
{
    *err = common_adios_finalize (*mype);
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_allocate_buffer, adios_ALLOCATE_BUFFER) (int *sizeMB, int * err)
{
    //FIX
    *err = common_adios_allocate_buffer (1 /*NOW*/, *sizeMB);
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_open, adios_OPEN) 
    (int64_t * fd, const char * group_name, const char * name
    ,const char * mode, void * comm, int * err
    ,int group_name_size, int name_size, int mode_size
    )
{
    char * buf1 = 0;
    char * buf2 = 0;
    char * buf3 = 0;

    buf1 = futils_fstr_to_cstr (group_name, group_name_size);
    buf2 = futils_fstr_to_cstr (name, name_size);
    buf3 = futils_fstr_to_cstr (mode, mode_size);

    if (buf1 != 0 && buf2 != 0 && buf3 != 0) {
        *err = common_adios_open (fd, buf1, buf2, buf3, comm);
        free (buf1);
        free (buf2);
        free (buf3);
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_group_size, adios_GROUP_SIZE) 
    (int64_t * fd_p, int64_t * data_size
    ,int64_t * total_size, int * err
    )
{
    *err = common_adios_group_size (*fd_p, (uint64_t) *data_size
                            ,(uint64_t *) total_size
                            );
}

///////////////////////////////////////////////////////////////////////////////
#include "stdio.h"
/* This Fortran api function is a bit different from the C api funcion, but
 * they call the same common_adios_write().
 * Difference: if the variable is string type then we need to convert
 * the void * var to a C string (add \0 to the end)
 * We rely on the assumption/fact that Fortran compilers pass on the 
 * length of a character array in an extra integer argument, even if
 * the C function declares a void* array in the argument list. 
 */
void FC_FUNC_(adios_write, adios_WRITE) 
    (int64_t * fd_p, const char * name, void * var, int * err
    ,int name_size, int var_size
    )
{
    struct adios_file_struct * fd = (struct adios_file_struct *) *fd_p;
    if (!fd)
    {
        fprintf (stderr, "Invalid handle passed to adios_write\n");
        *err = 1;
        return;
    }

    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL)
    {
        // nothing to do so just return
        *err = 0;
        return;
    }

    char * buf1 = 0;
    buf1 = futils_fstr_to_cstr (name, name_size);

    //printf("  -- adios_write: name=[%s] var size = %d\n", buf1, var_size);

    if (!buf1) {
        *err = -adios_errno;
        return;
    }

    v = adios_find_var_by_name (v, buf1, fd->group->all_unique_var_names);

    if (!v)
    {
        fprintf (stderr, "Bad var name (ignored): '%s'\n", buf1);
        *err = 1;
        free (buf1);
        return;
    }

    if (fd->mode == adios_mode_read)
    {
        if (   strcasecmp (buf1, fd->group->group_comm)
            && v->is_dim != adios_flag_yes
           )
        {
            fprintf (stderr, "write attempted on %s in %s.  This was opened for read\n" ,buf1 , fd->name);
            *err = 1;
            free (buf1);
            return;
        }
    }

    if (!var)
    {
        fprintf (stderr, "Invalid data: %s\n", buf1);
        *err = 1;
        free (buf1);
        return;
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
                        free (hist->breaks);
                        free (hist->frequencies);
                        free (hist);
                        v->stats[c][idx].data = 0;
                    }
                    else
                    {   
                        free (v->stats[c][idx].data);
                        v->stats[c][idx].data = 0;
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
                    fprintf (stderr, "cannot allocate %llu bytes to copy "
                                     "scalar %s\n"
                            ,element_size
                            ,v->name
                            );

                    *err = 1;
                    free (buf1);
                    return;
                }

                memcpy ((char *) v->data, var, element_size);
                break;
            case adios_string:
                v->data = futils_fstr_to_cstr (var, var_size);
                if (!v->data)
                {
                    fprintf (stderr, "cannot allocate %llu bytes to copy "
                                     "scalar %s\n"
                            ,element_size
                            ,v->name
                            );

                    *err = 1;
                    free (buf1);
                    return;
                }
                break;

            default:
                v->data = 0;
                break;
        }
    }

    *err = common_adios_write (fd, v, var);
    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        adios_copy_var_written (&fd->group->vars_written, v, fd);
    }

    free (buf1);
}


///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_get_write_buffer, adios_GET_WRITE_BUFFER) 
    (int64_t * fd_p, const char * name
    ,int64_t * size
    ,void ** buffer, int * err, int name_size
    )
{
    char * buf1 = 0;

    buf1 = futils_fstr_to_cstr (name, name_size);

    if (buf1 != 0) {
        *err = common_adios_get_write_buffer (*fd_p, buf1, (uint64_t *) size, buffer);
        free (buf1);
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_read, adios_READ) 
    (int64_t * fd_p, const char * name, void * buffer
    ,int64_t * buffer_size, int * err, int name_size
    )
{
    char * buf1 = 0;

    buf1 = futils_fstr_to_cstr (name, name_size);

    if (buf1 != 0) {
        *err = common_adios_read (*fd_p, buf1, buffer, *buffer_size);
        free (buf1);
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_set_path, adios_SET_PATH) 
    (int64_t * fd_p, const char * path, int * err, int path_size)
{
    char * buf1 = 0;

    buf1 = futils_fstr_to_cstr (path, path_size);

    if (buf1 != 0) {
        *err = common_adios_set_path (*fd_p, buf1);
        free (buf1);
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_set_path_var, adios_SET_PATH_VAR) 
    (int64_t * fd_p, const char * path, const char * name, int * err, int path_size, int name_size)
{
    char * buf1 = 0;
    char * buf2 = 0;

    buf1 = futils_fstr_to_cstr (path, path_size);
    buf2 = futils_fstr_to_cstr (name, name_size);

    if (buf1 != 0 && buf2 != 0) {
        *err = common_adios_set_path_var (*fd_p, buf1, buf2);
        free (buf1);
        free (buf2);
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
// hint that we reached the end of an iteration (for asynchronous pacing)
void FC_FUNC_(adios_end_iteration, adios_END_ITERATION) (int * err)
{
    *err = common_adios_end_iteration ();
}

///////////////////////////////////////////////////////////////////////////////
// hint to start communicating
void FC_FUNC_(adios_start_calculation, adios_START_CALCULATION) (int * err)
{
    *err = common_adios_start_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
// hint to stop communicating
void FC_FUNC_(adios_stop_calculation, adios_STOP_CALCULATION) (int * err)
{
    *err = common_adios_stop_calculation ();
}

///////////////////////////////////////////////////////////////////////////////
void FC_FUNC_(adios_close, adios_CLOSE) (int64_t * fd_p, int * err)
{
    *err = common_adios_close (*fd_p);
}

//////////////////////////////////////////////////////////////////////////////
// Methods normally only called by the XML parser
//////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// adios_common_declare_group is in adios_internals.c
// group a list of vars into a composite group
void FC_FUNC_(adios_declare_group, adios_DECLARE_GROUP) 
    (int64_t * id, const char * name
    ,const char * time_index, enum ADIOS_FLAG stats
    ,int * err, int name_size, int time_index_size
    )
{
    char * buf1 = 0;
    char * buf2 = 0;

    buf1 = futils_fstr_to_cstr (name, name_size);
    buf2 = futils_fstr_to_cstr (time_index, time_index_size);

    if (buf1 != 0 && buf2 != 0) {
        *err = adios_common_declare_group (id, buf1, adios_flag_yes, "", "", buf2, stats);
        free (buf1);
        free (buf2);
        if (*err == 1) {
            struct adios_group_struct * g = (struct adios_group_struct *) *id;
            g->all_unique_var_names = adios_flag_no;
        }
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
// adios_common_define_var is in adios_internals.c
// declare a single var as an entry in a group
void FC_FUNC_(adios_define_var, adios_DEFINE_VAR) 
    (int64_t * group_id, const char * name
    ,const char * path, int * type
    ,const char * dimensions
    ,const char * global_dimensions
    ,const char * local_offsets, int * err
    ,int name_size, int path_size, int dimensions_size
    ,int global_dimensions_size, int local_offsets_size
    )
{
    char * buf1 = 0;
    char * buf2 = 0;
    char * buf3 = 0;
    char * buf4 = 0;
    char * buf5 = 0;

    buf1 = futils_fstr_to_cstr (name, name_size);
    buf2 = futils_fstr_to_cstr (path, path_size);
    buf3 = futils_fstr_to_cstr (dimensions, dimensions_size);
    buf4 = futils_fstr_to_cstr (global_dimensions, global_dimensions_size);
    buf5 = futils_fstr_to_cstr (local_offsets, local_offsets_size);

    if (buf1 != 0 && buf2 != 0) {
        *err = adios_common_define_var (*group_id, buf1, buf2
                                       ,(enum ADIOS_DATATYPES) *type
                                       ,buf3, buf4, buf5
                                       );

        free (buf1);
        free (buf2);
        free (buf3);
        free (buf4);
        free (buf5);
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
// adios_common_define_attribute is in adios_internals.c
void FC_FUNC_(adios_define_attribute, adios_DEFINE_ATTRIBUTE) 
    (int64_t * group, const char * name
    ,const char * path, int * type, const char * value
    ,const char * var, int * err
    ,int name_size, int path_size, int value_size
    ,int var_size
    )
{
    char * buf1 = 0;
    char * buf2 = 0;
    char * buf3 = 0;
    char * buf4 = 0;

    buf1 = futils_fstr_to_cstr (name, name_size);
    buf2 = futils_fstr_to_cstr (path, path_size);
    buf3 = futils_fstr_to_cstr (value, value_size);
    buf4 = futils_fstr_to_cstr (var, var_size);

    if (buf1 != 0 && buf2 != 0 && buf3 != 0 && buf4 != 0) {
        *err = adios_common_define_attribute (*group, buf1, buf2
                                             ,(enum ADIOS_DATATYPES) *type, buf3
                                             ,buf4
                                             );

        free (buf1);
        free (buf2);
        free (buf3);
        free (buf4);
    } else {
        *err = -adios_errno;
    }
}

///////////////////////////////////////////////////////////////////////////////
// adios_common_select_method is in adios_internals_mxml.c
void FC_FUNC_(adios_select_method, adios_SELECT_METHOD) 
    (int64_t * group, const char * method
    ,const char * parameters, const char * base_path
    ,int * err, int method_size, int parameters_size
    ,int base_path_size
    )
{
    char * buf1 = 0;
    char * buf2 = 0;
    char * buf3 = 0;
    buf1 = futils_fstr_to_cstr (method, method_size);
    buf2 = futils_fstr_to_cstr (parameters, parameters_size);
    buf3 = futils_fstr_to_cstr (base_path, base_path_size);

    if (buf1 != 0 && buf2 != 0 && buf3 != 0) {
        struct adios_group_struct * g = (struct adios_group_struct *) (* group);
        *err = adios_common_select_method (0, buf1, buf2, g->name, buf3, 0);

        free (buf1);
        free (buf2);
        free (buf3);
    } else {
        *err = -adios_errno;
    }
}
