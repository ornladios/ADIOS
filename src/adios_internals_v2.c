/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <math.h>
#include <string.h>
#include <ctype.h>  /* isdigit() */
#include <unistd.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <arpa/inet.h>
#include <stdint.h>
#include <sys/stat.h>

#include "adios.h"
#include "adios_internals.h"
#include "adios_internals_v2.h"
#include "adios_bp_v2.h"

static struct adios_method_list_struct * adios_methods = 0;
static struct adios_group_list_struct * adios_groups = 0;

uint16_t adios_calc_var_overhead_v2 (struct adios_var_struct * v)
{
    uint16_t overhead = 0;

    struct adios_dimension_struct * d = v->dimensions;

    overhead += 8; // length of var entry
    overhead += 2; // member id
    overhead += 2; // length of name
    overhead += strlen (v->name); // name
    overhead += 2; // length of path
    overhead += strlen (v->path); // path
    overhead += 1; // datatype
    overhead += 1; // used as a dimension flag

    overhead += 1; // ranks
    overhead += 2; // dimensions length
    while (d)
    {
        overhead += 1; // var flag
        if (   d->dimension.id == 0
            && d->dimension.time_index == adios_flag_no
           )
        {
            overhead += 8; // value
        }
        else
        {
            overhead += 2; // member id
        }

        overhead += 1; // var flag
        if (   d->global_dimension.id == 0
            && d->global_dimension.time_index == adios_flag_no
           )
        {
            overhead += 8; // value
        }
        else
        {
            overhead += 2; // member id
        }

        overhead += 1; // var flag
        if (   d->local_offset.id == 0
            && d->local_offset.time_index == adios_flag_no
           )
        {
            overhead += 8; // value
        }
        else
        {
            overhead += 2; // member id
        }

        d = d->next;
    }
    overhead += adios_calc_var_characteristics_overhead (v);

    return overhead;
}

uint32_t adios_calc_attribute_overhead_v2 (struct adios_attribute_struct * a)
{
    uint32_t overhead = 0;

    overhead += 4; // attribute length
    overhead += 2; // member id
    overhead += 2; // length of name
    overhead += strlen (a->name); // name
    overhead += 2; // length of path
    overhead += strlen (a->path); // path
    overhead += 1; // var flag
    if (a->var)
        overhead += 2; // var member id
    else
    {
        overhead += 1; // datatype
        overhead += 4; // length of value
        overhead += adios_get_type_size (a->type, a->value); // value
    }

    return overhead;
}

uint64_t adios_calc_overhead_v2 (struct adios_file_struct * fd)
{
    uint64_t overhead = 0;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_attribute_struct * a = fd->group->attributes;
    struct adios_method_list_struct * m = fd->group->methods;

    overhead += 8; // process group length
    overhead += 1; // host language flag
    overhead += 2; // length of group name
    overhead += strlen (fd->group->name); // group name
    overhead += 2; // coordination var id
    overhead += 2; // length of time index name
    overhead += ((fd->group->time_index_name)
                    ? strlen (fd->group->time_index_name)
                    : 0
                );  // time index name
    overhead += 4; // time index

    overhead += 1; // count of methods employed
    overhead += 2; // length of methods section

    while (m)
    {
        overhead += 1; // method ID
        overhead += 2; // method params length
        overhead += strlen (m->method->parameters);
        m = m->next;
    }

    overhead += 2; // count of vars
    overhead += 8; // length of vars section

    while (v)
    {
        overhead += adios_calc_var_overhead_v2 (v);

        v = v->next;
    }

    overhead += 2; // attributes count
    overhead += 8; // attributes length

    while (a)
    {
        overhead += adios_calc_attribute_overhead_v2 (a);

        a = a->next;
    }

    return overhead;
}

int adios_write_process_group_header_v2 (struct adios_file_struct * fd
                                        ,uint64_t total_size
                                        )
{
    struct adios_group_struct * g = fd->group;

    uint8_t flag;
    struct adios_var_struct * var;
    uint16_t len;
    char * temp;

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &total_size, 8);

    flag = (g->adios_host_language_fortran == adios_flag_yes ? 'y' : 'n');
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &flag, 1);

    len = strlen (g->name);
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &len, 2);

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, g->name, len);

    var = adios_find_var_by_name (g->vars, g->group_by
                                 ,g->all_unique_var_names
                                 );
    if (var)
    {
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var->id, 2);
    }
    else
    {
        uint16_t i = 0;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &i, 2);
    }

    len = ((g->time_index_name) ? strlen (g->time_index_name) : 0);
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &len, 2);

    if (g->time_index_name)
    {
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,g->time_index_name, len
                     );
    }
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                 ,&g->time_index, 4
                 );

    struct adios_method_list_struct * m = fd->group->methods;
    uint8_t methods_count = 0;
    uint16_t methods_length = 0;
    while (m)
    {
        methods_count++;
        methods_length += 1 + 2 + strlen (m->method->parameters);

        m = m->next;
    }
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                 ,&methods_count, 1
                 );
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                 ,&methods_length, 2
                 );

    m = fd->group->methods;
    while (m)
    {
        uint16_t len = strlen (m->method->parameters);

        flag = (uint8_t) m->method->m;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &flag, 1);

        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &len, 2);

        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,m->method->parameters, len
                     );

        m = m->next;
    }

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return 0;
}

static void index_append_process_group_v2 (
                          struct adios_index_process_group_struct_v2 ** root
                         ,struct adios_index_process_group_struct_v2 * item
                         )
{
    while (root)
    {
        if (!*root)
        {
            *root = item;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}

static void index_append_var_v2 (struct adios_index_var_struct_v2 ** root
                                ,struct adios_index_var_struct_v2 * item
                                )
{
    while (root)
    {
        if (!*root)
        {
            *root = item;
            root = 0;
        }
        else
        {
            if (   !strcasecmp (item->group_name, (*root)->group_name)
                && !strcasecmp (item->var_name, (*root)->var_name)
                && !strcasecmp (item->var_path, (*root)->var_path)
                && item->type == (*root)->type
               )
            {
                if (    (*root)->characteristics_count
                      + item->characteristics_count
                    > (*root)->characteristics_allocated
                   )
                {
                    int new_items = (item->characteristics_count == 1)
                                         ? 100 : item->characteristics_count;
                    (*root)->characteristics_allocated = 
                            (*root)->characteristics_count + new_items;
                    void * ptr;
                    ptr = realloc ((*root)->characteristics
                            ,  (*root)->characteristics_allocated
                        * sizeof (struct adios_index_characteristic_struct_v2)
                            );

                    if (ptr)
                    {
                        (*root)->characteristics = ptr;
                    }
                    else
                    {
                        fprintf (stderr, "error allocating memory to build "
                                         "var index.  Index aborted\n"
                                );

                        return;
                    }
                }
                memcpy (&(*root)->characteristics
                                             [(*root)->characteristics_count]
                       ,item->characteristics
                       ,  item->characteristics_count
                        * sizeof (struct adios_index_characteristic_struct_v2)
                       );

                (*root)->characteristics_count += item->characteristics_count;

                free (item->characteristics);
                free (item->group_name);
                free (item->var_name);
                free (item->var_path);
                free (item);

                root = 0;  // exit the loop
            }
            else
            {
                root = &(*root)->next;
            }
        }
    }
}

static void index_append_attribute_v2
                                (struct adios_index_attribute_struct_v2 ** root
                                ,struct adios_index_attribute_struct_v2 * item
                                )
{
    while (root)
    {
        if (!*root)
        {
            *root = item;
            root = 0;
        }
        else
        {
            if (   !strcasecmp (item->group_name, (*root)->group_name)
                && !strcasecmp (item->attr_name, (*root)->attr_name)
                && !strcasecmp (item->attr_path, (*root)->attr_path)
               )
            {
                if (    (*root)->characteristics_count
                      + item->characteristics_count
                    > (*root)->characteristics_allocated
                   )
                {
                    int new_items = (item->characteristics_count == 1)
                                         ? 100 : item->characteristics_count;
                    (*root)->characteristics_allocated =
                                   (*root)->characteristics_count + new_items;
                    void * ptr;
                    ptr = realloc ((*root)->characteristics
                            ,  (*root)->characteristics_allocated
                       * sizeof (struct adios_index_characteristic_struct_v2)
                            );

                    if (ptr)
                    {
                        (*root)->characteristics = ptr;
                    }
                    else
                    {
                        fprintf (stderr, "error allocating memory to build "
                                         "attribute index.  Index aborted\n"
                                );

                        return;
                    }
                }
                memcpy (&(*root)->characteristics
                                              [(*root)->characteristics_count]
                       ,item->characteristics
                       ,  item->characteristics_count
                        * sizeof (struct adios_index_characteristic_struct_v2)
                       );

                (*root)->characteristics_count += item->characteristics_count;

                free (item->characteristics);
                free (item->group_name);
                free (item->attr_name);
                free (item->attr_path);
                free (item);

                root = 0;  // exit the loop
            }
            else
            {
                root = &(*root)->next;
            }
        }
    }
}

// p2 and v2 will be destroyed as part of the merge operation...
void adios_merge_index_v2 (struct adios_index_process_group_struct_v2 ** p1
                          ,struct adios_index_var_struct_v2 ** v1
                          ,struct adios_index_attribute_struct_v2 ** a1
                          ,struct adios_index_process_group_struct_v2 * p2
                          ,struct adios_index_var_struct_v2 * v2
                          ,struct adios_index_attribute_struct_v2 * a2
                          )
{
    // this will just add it on to the end and all should work fine
    index_append_process_group_v2 (p1, p2);

    // need to do vars attrs one at a time to merge them properly
    struct adios_index_var_struct_v2 * v_temp;
    struct adios_index_attribute_struct_v2 * a_temp;

    while (v2)
    {
        v_temp = v2->next;
        v2->next = 0;
        index_append_var_v2 (v1, v2);
        v2 = v_temp;
    }

    while (a2)
    {
        a_temp = a2->next;
        a2->next = 0;
        index_append_attribute_v2 (a1, a2);
        a2 = a_temp;
    }
}

static void adios_clear_process_groups_index_v2 (
                            struct adios_index_process_group_struct_v2 * root
                           )
{
    while (root)
    {
        struct adios_index_process_group_struct_v2 * temp = root->next;
        if (root->group_name)
            free (root->group_name);
        if (root->time_index_name)
            free (root->time_index_name);
        free (root);
        root = temp;
    }
}

static void adios_clear_vars_index_v2 (struct adios_index_var_struct_v2 * root)
{
    while (root)
    {
        int i;
        struct adios_index_var_struct_v2 * temp = root->next;

        if (root->group_name)
            free (root->group_name);
        if (root->var_name)
            free (root->var_name);
        if (root->var_path)
            free (root->var_path);
        for (i = 0; i < root->characteristics_count; i++)
        {
            if (root->characteristics [i].dims.count != 0)
                free (root->characteristics [i].dims.dims);
            if (root->characteristics [i].min)
                free (root->characteristics [i].min);
            if (root->characteristics [i].max)
                free (root->characteristics [i].max);
            if (root->characteristics [i].value)
                free (root->characteristics [i].value);
        }
        if (root->characteristics)
            free (root->characteristics);

        free (root);
        root = temp;
    }
}

static void adios_clear_attributes_index_v2
                                (struct adios_index_attribute_struct_v2 * root)
{
    while (root)
    {
        int i;
        struct adios_index_attribute_struct_v2 * temp = root->next;

        if (root->group_name)
            free (root->group_name);
        if (root->attr_name)
            free (root->attr_name);
        if (root->attr_path)
            free (root->attr_path);
        for (i = 0; i < root->characteristics_count; i++)
        {
            if (root->characteristics [i].dims.count != 0)
                free (root->characteristics [i].dims.dims);
            if (root->characteristics [i].min)
                free (root->characteristics [i].min);
            if (root->characteristics [i].max)
                free (root->characteristics [i].max);
            if (root->characteristics [i].value)
                free (root->characteristics [i].value);
        }
        if (root->characteristics)
            free (root->characteristics);

        free (root);
        root = temp;
    }
}

void adios_clear_index_v2 (struct adios_index_process_group_struct_v2 * pg_root
                          ,struct adios_index_var_struct_v2 * vars_root
                          ,struct adios_index_attribute_struct_v2 * attrs_root
                          )
{
    adios_clear_process_groups_index_v2 (pg_root);
    adios_clear_vars_index_v2 (vars_root);
    adios_clear_attributes_index_v2 (attrs_root);
}

static uint8_t count_dimensions (struct adios_dimension_struct * dimensions)
{
    uint8_t count = 0;

    while (dimensions)
    {
        count++;
        dimensions = dimensions->next;
    }

    return count;
}

static uint64_t cast_var_data_as_uint64 (const char * parent_name
                                        ,enum ADIOS_DATATYPES type
                                        ,void * data
                                        )
{
    if (!data)
    {
        fprintf (stderr, "cannot write var since dim %s not provided\n"
                ,parent_name
                );
        return 0;
    }

    switch (type)
    {
        case adios_byte:
            return (uint64_t) *(int8_t *) data;

        case adios_short:
            return (uint64_t) *(int16_t *) data;

        case adios_integer:
            return (uint64_t) *(int32_t *) data;

        case adios_long:
            return (uint64_t) *(int64_t *) data;

        case adios_unsigned_byte:
            return (uint64_t) *(uint8_t *) data;

        case adios_unsigned_short:
            return (uint64_t) *(uint16_t *) data;

        case adios_unsigned_integer:
            return (uint64_t) *(uint32_t *) data;

        case adios_unsigned_long:
            return (uint64_t) *(uint64_t *) data;

        case adios_real:
            return (uint64_t) *(float *) data;

        case adios_double:
            return (uint64_t) *(double *) data;

        case adios_long_double:
            return (uint64_t) *(long double *) data;

        case adios_string:
        case adios_complex:
        case adios_double_complex:
            fprintf (stderr, "Cannot convert type %s to integer for var %s\n"
                    ,adios_type_to_string_int (type), parent_name
                    );

            return 0;
    }
    return 0;
}

static uint64_t get_value_for_dim (struct adios_file_struct * fd
                                ,struct adios_dimension_item_struct * dimension
                                )
{
    uint64_t dim = 0;

    if (dimension->id != 0)
    {
        struct adios_var_struct * var = adios_find_var_by_id (fd->group->vars
                                                             ,dimension->id
                                                             );
        if (var)
        {
            if (var->data)
            {
                dim = cast_var_data_as_uint64 (var->name, var->type, var->data);
            }
            else
            {
                fprintf (stderr, "array dimension data missing\n");
            }
        }
        else
        {
            struct adios_attribute_struct * attr = adios_find_attribute_by_id
                                                        (fd->group->attributes
                                                        ,dimension->id
                                                        );
            if (attr)
            {
                if (attr->var)
                {
                    if (attr->var->data)
                    {
                        dim = cast_var_data_as_uint64 (attr->var->name
                                                      ,attr->var->type
                                                      ,attr->var->data
                                                      );
                    }
                    else
                    {
                        fprintf (stderr, "array dimension data missing\n");
                    }
                }
                else
                {
                    dim = cast_var_data_as_uint64 (attr->name, attr->type
                                                  ,attr->value
                                                  );
                }
            }
            else
            {
                fprintf (stderr, "invalid dimension member id: %d\n"
                        ,dimension->id
                        );
            }
        }
    }
    else
    {
        if (dimension->time_index == adios_flag_yes)
            dim = 1;
        else
            dim = dimension->rank;
    }

    return dim;
}

void adios_build_index_v2 (struct adios_file_struct * fd
                       ,struct adios_index_process_group_struct_v2 ** pg_root
                       ,struct adios_index_var_struct_v2 ** vars_root
                       ,struct adios_index_attribute_struct_v2 ** attrs_root
                       )
{
    struct adios_group_struct * g = fd->group;
    struct adios_var_struct * v = g->vars;
    struct adios_attribute_struct * a = g->attributes;
    struct adios_index_process_group_struct_v2 * g_item;

    uint64_t process_group_count = 0;
    uint16_t var_count = 0;

    g_item = (struct adios_index_process_group_struct_v2 *)
                malloc (sizeof (struct adios_index_process_group_struct_v2));
    g_item->group_name = (g->name ? strdup (g->name) : 0L);
    g_item->adios_host_language_fortran = g->adios_host_language_fortran;
    g_item->process_id = g->process_id;
    g_item->time_index_name = (g->time_index_name ? strdup (g->time_index_name) : 0L);
    g_item->time_index = g->time_index;
    g_item->offset_in_file = fd->pg_start_in_file;
    g_item->next = 0;

    // build the groups and vars index
    index_append_process_group_v2 (pg_root, g_item);

    while (v)
    {
        // only add items that were written to the index
        if (v->write_offset != 0)
        {
            struct adios_index_var_struct_v2 * v_index;
            v_index = malloc (sizeof (struct adios_index_var_struct_v2));
            v_index->characteristics = malloc (
                           sizeof (struct adios_index_characteristic_struct_v2)
                          );

            v_index->id = v->id;
            v_index->group_name = (g->name ? strdup (g->name) : 0L);
            v_index->var_name = (v->name ? strdup (v->name) : 0L);
            v_index->var_path = (v->path ? strdup (v->path) : 0L);
            v_index->type = v->type;
            v_index->characteristics_count = 1;
            v_index->characteristics_allocated = 1;
            v_index->characteristics [0].offset = v->write_offset;
            v_index->characteristics [0].payload_offset = v->write_offset + adios_calc_var_overhead_v2 (v);
            v_index->characteristics [0].min = 0;
            v_index->characteristics [0].max = 0;
            v_index->characteristics [0].value = 0;
            v_index->characteristics [0].dims.count = 0;

            uint64_t size = adios_get_type_size (v->type, v->data);
            switch (v->type)
            {
                case adios_byte:
                case adios_unsigned_byte:
                case adios_short:
                case adios_unsigned_short:
                case adios_integer:
                case adios_unsigned_integer:
                case adios_long:
                case adios_unsigned_long:
                case adios_real:
                case adios_double:
                case adios_long_double:
                case adios_complex:
                case adios_double_complex:
                    if (v->dimensions)
                    {
                        uint8_t c;
                        uint8_t j;
                        struct adios_dimension_struct * d = v->dimensions;
                        v_index->characteristics [0].min = malloc (size);
                        v_index->characteristics [0].max = malloc (size);
                        memcpy (v_index->characteristics [0].min, v->min, size);
                        memcpy (v_index->characteristics [0].max, v->max, size);
                        c = count_dimensions (v->dimensions);
                        v_index->characteristics [0].dims.count = c;
                        // (local, global, local offset)
                        v_index->characteristics [0].dims.dims = malloc
                            (3 * 8 * v_index->characteristics [0].dims.count);
                        for (j = 0; j < c; j++)
                        {
                            v_index->characteristics [0].dims.dims [j * 3 + 0] =
                                   get_value_for_dim (fd, &d->dimension);
                            v_index->characteristics [0].dims.dims [j * 3 + 1] =
                                   get_value_for_dim (fd, &d->global_dimension);
                            v_index->characteristics [0].dims.dims [j * 3 + 2] =
                                   get_value_for_dim (fd, &d->local_offset);

                            d = d->next;
                        }
                        v_index->characteristics [0].value = 0;
                    }
                    else
                    {
                        v_index->characteristics [0].min = 0;
                        v_index->characteristics [0].max = 0;
                        v_index->characteristics [0].value = malloc (size);
                        memcpy (v_index->characteristics [0].value, v->data
                               ,size
                               );
                        v_index->characteristics [0].dims.count = 0;
                        v_index->characteristics [0].dims.dims = 0;
                    }

                    break;

                case adios_string:
                {
                    v_index->characteristics [0].value = malloc (size + 1);
                    memcpy (v_index->characteristics [0].value, v->data, size);
                    ((char *) (v_index->characteristics [0].value)) [size] = 0;

                    break;
                }
            }
            v_index->next = 0;

            // this fn will either take ownership for free
            index_append_var_v2 (vars_root, v_index);
        }

        v = v->next;
    }

    while (a)
    {
        // only add items that were written to the index
        if (a->write_offset != 0)
        {
            struct adios_index_attribute_struct_v2 * a_index;
            a_index = malloc (sizeof (struct adios_index_attribute_struct_v2));
            a_index->characteristics = malloc (
                           sizeof (struct adios_index_characteristic_struct_v2)
                          );

            a_index->id = a->id;
            a_index->group_name = (g->name ? strdup (g->name) : 0L);
            a_index->attr_name = (a->name ? strdup (a->name) : 0L);
            a_index->attr_path = (a->path ? strdup (a->path) : 0L);
            a_index->type = a->type;
            a_index->characteristics_count = 1;
            a_index->characteristics_allocated = 1;
            uint64_t size = adios_get_type_size (a->type, a->value);

            a_index->characteristics [0].offset = a->write_offset;
            a_index->characteristics [0].payload_offset = a->write_offset + adios_calc_attribute_overhead_v2 (a);
            a_index->characteristics [0].min = 0;
            a_index->characteristics [0].max = 0;
            if (a->value)
            {
                a_index->characteristics [0].value = malloc (size + 1);
                ((char *) (a_index->characteristics [0].value)) [size] = 0;
                memcpy (a_index->characteristics [0].value, a->value, size);
            }
            else
            {
                a_index->characteristics [0].value = 0;
            }
            a_index->characteristics [0].dims.count = 0;
            a_index->characteristics [0].dims.dims = 0;
            if (a->var)
                a_index->characteristics [0].var_id = a->var->id;
            else
                a_index->characteristics [0].var_id = 0;

            a_index->next = 0;

            // this fn will either take ownership for free
            index_append_attribute_v2 (attrs_root, a_index);
        }

        a = a->next;
    }
}

int adios_write_index_v2 (char ** buffer
                         ,uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,uint64_t index_start
                         ,struct adios_index_process_group_struct_v2 * pg_root
                         ,struct adios_index_var_struct_v2 * vars_root
                         ,struct adios_index_attribute_struct_v2 * attrs_root
                         )
{
    uint64_t groups_count = 0;
    uint16_t vars_count = 0;
    uint16_t attrs_count = 0;

    uint64_t index_size = 0;
    uint64_t pg_index_start = index_start;
    uint64_t vars_index_start = 0;
    uint64_t attrs_index_start = 0;

    // we need to save the offset we will write the count and size
    uint64_t buffer_offset_start = 0; // since we realloc, we can't save a ptr

    // save for the process group index
    buffer_offset_start = *buffer_offset;

    *buffer_offset += (8 + 8); // save space for groups count and index size

    while (pg_root)
    {
        uint8_t flag;
        uint16_t len;
        uint16_t group_size = 0;
        uint64_t group_start = *buffer_offset;

        groups_count++;

        *buffer_offset += 2; // save space for the size

        len = strlen (pg_root->group_name);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        group_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset
                     ,pg_root->group_name, len
                     );
        index_size += len;
        group_size += len;

        flag = (pg_root->adios_host_language_fortran == adios_flag_yes ? 'y'
                                                                       : 'n'
               );
        buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
        index_size += 1;
        group_size += 1;

        buffer_write (buffer, buffer_size, buffer_offset
                     ,&pg_root->process_id, 4
                     );
        index_size += 4;
        group_size += 4;

        if (pg_root->time_index_name)
        {
            len = strlen (pg_root->time_index_name);
        }
        else
        {
            len = 0;
        }
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        group_size += 2;
        if (len)
        {
            buffer_write (buffer, buffer_size, buffer_offset
                         ,pg_root->time_index_name, len
                         );
        }
        index_size += len;
        group_size += len;

        buffer_write (buffer, buffer_size, buffer_offset
                     ,&pg_root->time_index, 4
                     );
        index_size += 4;
        group_size += 4;
        buffer_write (buffer, buffer_size, buffer_offset
                     ,&pg_root->offset_in_file, 8
                     );
        index_size += 8;
        group_size += 8;

        buffer_write (buffer, buffer_size, &group_start, &group_size, 2);

        pg_root = pg_root->next;
    }

    buffer_write (buffer, buffer_size, &buffer_offset_start, &groups_count, 8);
    buffer_write (buffer, buffer_size, &buffer_offset_start, &index_size, 8);

    buffer_offset_start = *buffer_offset; // save to write the vars_count/size
    vars_index_start = buffer_offset_start + index_start;
    index_size = 0;

    *buffer_offset += (2 + 8); // save space for count and size

    while (vars_root)
    {
        uint8_t flag;
        uint16_t len;
        uint32_t var_size = 0;
        uint64_t var_start = *buffer_offset;
        int i;

        vars_count++;

        *buffer_offset += 4; // save space for var length

        buffer_write (buffer, buffer_size, buffer_offset, &vars_root->id, 2);
        index_size += 2;
        var_size += 2;

        len = strlen (vars_root->group_name);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        var_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset
                     ,vars_root->group_name, len
                     );
        index_size += len;
        var_size += len;

        len = strlen (vars_root->var_name);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        var_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset
                     ,vars_root->var_name, len
                     );
        index_size += len;
        var_size += len;

        len = strlen (vars_root->var_path);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        var_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset
                     ,vars_root->var_path, len
                     );
        index_size += len;
        var_size += len;

        flag = vars_root->type;
        buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
        index_size += 1;
        var_size += 1;

        buffer_write (buffer, buffer_size, buffer_offset
                     ,&vars_root->characteristics_count, 8
                     );
        index_size += 8;
        var_size += 8;

        for (i = 0; i < vars_root->characteristics_count; i++)
        {
            uint64_t size;
            uint8_t characteristic_set_count = 0;
            uint32_t characteristic_set_length = 0;

            uint64_t characteristic_set_start = *buffer_offset;
            *buffer_offset += 1 + 4; // save space for characteristic count/len
            index_size += 1 + 4;
            var_size += 1 + 4;
            
            // add an offset characteristic for all vars
            characteristic_set_count++;
            flag = (uint8_t) adios_characteristic_offset;
            buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
            index_size += 1;
            var_size += 1;
            characteristic_set_length += 1;

            buffer_write (buffer, buffer_size, buffer_offset
                         ,&vars_root->characteristics [i].offset, 8
                         );
            index_size += 8;
            var_size += 8;
            characteristic_set_length += 8;

            // add a payload offset characteristic for all vars
            characteristic_set_count++;
            flag = (uint8_t) adios_characteristic_payload_offset;
            buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
            index_size += 1;
            var_size += 1;
            characteristic_set_length += 1;

            buffer_write (buffer, buffer_size, buffer_offset
                         ,&vars_root->characteristics [i].payload_offset, 8
                         );
            index_size += 8;
            var_size += 8;
            characteristic_set_length += 8;

            // depending on if it is an array or not, generate a different
            // additional set of characteristics
            size = adios_get_type_size (vars_root->type
                                       ,vars_root->characteristics [i].value
                                       );

            switch (vars_root->type)
            {
                case adios_byte:
                case adios_unsigned_byte:
                case adios_short:
                case adios_unsigned_short:
                case adios_integer:
                case adios_unsigned_integer:
                case adios_long:
                case adios_unsigned_long:
                case adios_real:
                case adios_double:
                case adios_long_double:
                case adios_complex:
                case adios_double_complex:
                    if (vars_root->characteristics [i].dims.count)
                    {
                        // add a dimensions characteristic
                        characteristic_set_count++;
                        flag = (uint8_t) adios_characteristic_dimensions;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&flag, 1
                                     );
                        index_size += 1;
                        var_size += 1;
                        characteristic_set_length += 1;

                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&vars_root->characteristics [i].dims.count
                                     ,1
                                     );
                        index_size += 1;
                        var_size += 1;
                        characteristic_set_length += 1;

                        len = 3 * 8 * vars_root->characteristics [i].dims.count;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&len, 2
                                     );
                        index_size += 2;
                        var_size += 2;
                        characteristic_set_length += 2;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,vars_root->characteristics [i].dims.dims
                                     ,len
                                     );
                        index_size += len;
                        var_size += len;
                        characteristic_set_length += len;

                        // add a min value characteristic
                        characteristic_set_count++;
                        flag = (uint8_t) adios_characteristic_min;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&flag, 1
                                     );
                        index_size += 1;
                        var_size += 1;
                        characteristic_set_length += 1;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,vars_root->characteristics [i].min, size
                                     );
                        index_size += size;
                        var_size += size;
                        characteristic_set_length += size;

                        // add a max value characteristic
                        characteristic_set_count++;
                        flag = (uint8_t) adios_characteristic_max;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&flag, 1
                                     );
                        index_size += 1;
                        var_size += 1;
                        characteristic_set_length += 1;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,vars_root->characteristics [i].max, size
                                     );
                        index_size += size;
                        var_size += size;
                        characteristic_set_length += size;
                    }
                    else
                    {
                        // add a value characteristic
                        characteristic_set_count++;
                        flag = (uint8_t) adios_characteristic_value;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&flag, 1
                                     );
                        index_size += 1;
                        var_size += 1;
                        characteristic_set_length += 1;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,vars_root->characteristics [i].value, size
                                     );
                        index_size += size;
                        var_size += size;
                        characteristic_set_length += size;
                    }
                    break;

                case adios_string:
                    {
                        // add a value characteristic
                        characteristic_set_count++;
                        flag = (uint8_t) adios_characteristic_value;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&flag, 1
                                     );
                        index_size += 1;
                        var_size += 1;
                        characteristic_set_length += 1;
                        if (vars_root->type == adios_string)
                        {
                            uint16_t len = (uint16_t) size;
                            buffer_write (buffer, buffer_size, buffer_offset
                                         ,&len, 2
                                         );
                            index_size += 2;
                            var_size += 2;
                            characteristic_set_length += 2;
                        }
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,vars_root->characteristics [i].value, size
                                     );
                        index_size += size;
                        var_size += size;
                        characteristic_set_length += size;
                    }
                    break;
            }
            // characteristics count/size prefix
            buffer_write (buffer, buffer_size, &characteristic_set_start
                         ,&characteristic_set_count, 1
                         );
            buffer_write (buffer, buffer_size, &characteristic_set_start
                         ,&characteristic_set_length, 4
                         );
        }

        buffer_write (buffer, buffer_size, &var_start, &var_size, 4);

        vars_root = vars_root->next;
    }

    // vars index count/size prefix
    buffer_write (buffer, buffer_size, &buffer_offset_start, &vars_count, 2);
    buffer_write (buffer, buffer_size, &buffer_offset_start, &index_size, 8);

    buffer_offset_start = *buffer_offset; // save to write the attrs_count/size
    attrs_index_start = buffer_offset_start + index_start;
    index_size = 0;

    *buffer_offset += (2 + 8); // save space for count and size

    while (attrs_root)
    {
        uint8_t flag;
        uint16_t len;
        uint32_t attr_size = 0;
        uint64_t attr_start = *buffer_offset;
        int i;

        attrs_count++;

        *buffer_offset += 4; // save space for attr length

        buffer_write (buffer, buffer_size, buffer_offset, &attrs_root->id, 2);
        index_size += 2;
        attr_size += 2;

        len = strlen (attrs_root->group_name);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        attr_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset
                     ,attrs_root->group_name, len
                     );
        index_size += len;
        attr_size += len;

        len = strlen (attrs_root->attr_name);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        attr_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset
                     ,attrs_root->attr_name, len
                     );
        index_size += len;
        attr_size += len;

        len = strlen (attrs_root->attr_path);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        attr_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset
                     ,attrs_root->attr_path, len
                     );
        index_size += len;
        attr_size += len;

        flag = attrs_root->type;
        buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
        index_size += 1;
        attr_size += 1;

        buffer_write (buffer, buffer_size, buffer_offset
                     ,&attrs_root->characteristics_count, 8
                     );
        index_size += 8;
        attr_size += 8;

        for (i = 0; i < attrs_root->characteristics_count; i++)
        {
            uint64_t size;
            uint8_t characteristic_set_count = 0;
            uint32_t characteristic_set_length = 0;

            uint64_t characteristic_set_start = *buffer_offset;
            *buffer_offset += 1 + 4; // save space for characteristic count/len
            index_size += 1 + 4;
            attr_size += 1 + 4;
            
            // add an offset characteristic for all attrs
            characteristic_set_count++;
            flag = (uint8_t) adios_characteristic_offset;
            buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
            index_size += 1;
            attr_size += 1;
            characteristic_set_length += 1;

            buffer_write (buffer, buffer_size, buffer_offset
                         ,&attrs_root->characteristics [i].offset, 8
                         );
            index_size += 8;
            attr_size += 8;
            characteristic_set_length += 8;

            // add a payload offset characteristic for all attrs
            characteristic_set_count++;
            flag = (uint8_t) adios_characteristic_payload_offset;
            buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
            index_size += 1;
            attr_size += 1;
            characteristic_set_length += 1;

            buffer_write (buffer, buffer_size, buffer_offset
                         ,&attrs_root->characteristics [i].payload_offset, 8
                         );
            index_size += 8;
            attr_size += 8;
            characteristic_set_length += 8;

            size = adios_get_type_size (attrs_root->type
                                       ,attrs_root->characteristics [i].value
                                       );

            if (attrs_root->characteristics [i].value != 0)
            {
                // add a value characteristic
                characteristic_set_count++;
                flag = (uint8_t) adios_characteristic_value;
                buffer_write (buffer, buffer_size, buffer_offset
                             ,&flag, 1
                             );
                index_size += 1;
                attr_size += 1;
                characteristic_set_length += 1;
                if (attrs_root->type == adios_string)
                {
                    uint16_t len = (uint16_t) size;
                    buffer_write (buffer, buffer_size, buffer_offset
                                 ,&len, 2
                                 );
                    index_size += 2;
                    attr_size += 2;
                    characteristic_set_length += 2;
                }
                buffer_write (buffer, buffer_size, buffer_offset
                             ,attrs_root->characteristics [i].value, size
                             );
                index_size += size;
                attr_size += size;
                characteristic_set_length += size;
            }
            if (attrs_root->characteristics [i].var_id != 0)
            {
                // add a var id characteristic
                characteristic_set_count++;
                flag = (uint8_t) adios_characteristic_var_id;
                buffer_write (buffer, buffer_size, buffer_offset
                             ,&flag, 1
                             );
                index_size += 1;
                attr_size += 1;
                characteristic_set_length += 1;
                buffer_write (buffer, buffer_size, buffer_offset
                             ,&attrs_root->characteristics [i].var_id, 2
                             );
                index_size += 2;
                attr_size += 2;
                characteristic_set_length += 2;
            }

            // characteristics count/size prefix
            buffer_write (buffer, buffer_size, &characteristic_set_start
                         ,&characteristic_set_count, 1
                         );
            buffer_write (buffer, buffer_size, &characteristic_set_start
                         ,&characteristic_set_length, 4
                         );
        }

        buffer_write (buffer, buffer_size, &attr_start, &attr_size, 4);

        attrs_root = attrs_root->next;
    }

    // attrs index count/size prefix
    buffer_write (buffer, buffer_size, &buffer_offset_start, &attrs_count, 2);
    buffer_write (buffer, buffer_size, &buffer_offset_start, &index_size, 8);

    // location of the beginning of the indexes (first proc groups then vars)
    buffer_write (buffer, buffer_size, buffer_offset, &pg_index_start, 8);
    buffer_write (buffer, buffer_size, buffer_offset, &vars_index_start, 8);
    buffer_write (buffer, buffer_size, buffer_offset, &attrs_index_start, 8);

    return 0;
}

int adios_write_version_v2 (char ** buffer
                           ,uint64_t * buffer_size
                           ,uint64_t * buffer_offset
                           )
{
    uint32_t test = 1;

    if (!*(char *) &test)
        test = 0x80000000;
    else
        test = 0;

    test += 2;   // master index file version

    test = htonl (test);

    buffer_write (buffer, buffer_size, buffer_offset, &test, 4);

    return 0;
}

static uint16_t calc_dimension_size (struct adios_dimension_struct * dimension)
{
    uint16_t size = 0;

    size += 1; // var (y or n)

    if (   dimension->dimension.id == 0
        && dimension->dimension.time_index == adios_flag_no
       )  // it is a number
    {
        size += 8;  // size of value
    }
    else   // it is a var
    {
        size += 2;  // size of var ID
    }

    size += 1; // var (y or n)

    if (   dimension->global_dimension.id == 0
        && dimension->global_dimension.time_index == adios_flag_no
       )  // it is a number
    {
        size += 8; // default to a rank
    }
    else
    {
        size += 2;  // size of var ID
    }

    size += 1; // var (y or n)

    if (   dimension->local_offset.id == 0
        && dimension->local_offset.time_index == adios_flag_no
       )  // it is a number
    {
        size += 8;  // default to a rank
    }
    else
    {
        size += 2;  // size of var ID
    }

    return size;
}

static uint16_t calc_dimensions_size (struct adios_dimension_struct * dimension)
{
    uint16_t size = 0;

    while (dimension)
    {
        size += calc_dimension_size (dimension);

        dimension = dimension->next;
    }

    return size;
}

static
uint64_t adios_write_dimension_v2 (struct adios_file_struct * fd
                                  ,struct adios_dimension_struct * dimension
                                  )
{
    uint64_t size = 0;
    uint8_t var;

    if (   dimension->dimension.id == 0
        && dimension->dimension.time_index == adios_flag_no
       )
    {
        var = 'n';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,&dimension->dimension.rank, 8
                     );
        size += 8;
    }
    else
    {
        var = 'y';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,&dimension->dimension.id, 2
                     );
        size += 2;
    }

    if (   dimension->global_dimension.id == 0
        && dimension->global_dimension.time_index == adios_flag_no
       )
    {
        var = 'n';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,&dimension->global_dimension.rank, 8
                     );
        size += 8;
    }
    else
    {
        var = 'y';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,&dimension->global_dimension.id, 2
                     );
        size += 2;
    }

    if (   dimension->local_offset.id == 0
        && dimension->local_offset.time_index == adios_flag_no
       )
    {
        var = 'n';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,&dimension->local_offset.rank, 8
                     );
        size += 8;
    }
    else
    {
        var = 'y';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,&dimension->local_offset.id, 2
                     );
        size += 2;
    }

    return size;
}

static
uint16_t adios_write_dimensions_v2 (struct adios_file_struct * fd
                                   ,struct adios_dimension_struct * dimensions
                                   )
{
    uint16_t size = 0;
    uint16_t dimensions_size = calc_dimensions_size (dimensions);
    uint8_t ranks = count_dimensions (dimensions);

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &ranks, 1);
    size += 1;
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dimensions_size, 2);
    size += 2;

    while (dimensions)
    {
        size += adios_write_dimension_v2 (fd, dimensions);

        dimensions = dimensions->next;
    }

    return size;
}

uint16_t adios_write_var_characteristics_dims_v2 (struct adios_file_struct * fd
                                                 ,struct adios_var_struct * v
                                                 )
{
    uint16_t total_size = 0;
    uint8_t dims_count = 0;
    uint16_t dims_length = 0;
    struct adios_dimension_struct * d = v->dimensions;
    uint64_t count_offset = fd->offset;

    fd->offset += 1;
    total_size += 1; // count

    fd->offset += 2;
    total_size += 2; // length

    while (d)
    {
        uint64_t dim = 0;

        dims_count++;

        dim = get_value_for_dim (fd, &d->dimension);
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dim, 8);
        total_size += 8;
        dims_length += 8;

        dim = get_value_for_dim (fd, &d->global_dimension);
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dim, 8);
        total_size += 8;
        dims_length += 8;

        dim = get_value_for_dim (fd, &d->local_offset);
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dim, 8);
        total_size += 8;
        dims_length += 8;

        d = d->next;
    }

    buffer_write (&fd->buffer, &fd->buffer_size, &count_offset, &dims_count, 1);
    buffer_write (&fd->buffer, &fd->buffer_size, &count_offset, &dims_length, 2);

    return total_size;
}

uint16_t adios_write_var_characteristics_v2 (struct adios_file_struct * fd
                                            ,struct adios_var_struct * v
                                            )
{
    uint8_t flag;
    uint64_t size;
    uint16_t len;
    uint8_t characteristic_set_count = 0;
    uint32_t characteristic_set_length = 0;
    uint64_t index_size = 0;

    uint64_t characteristic_set_start = fd->offset;
    fd->offset += 1 + 4; // save space for characteristic count/len
    index_size += 1 + 4;
    
    // depending on if it is an array or not, generate a different
    // additional set of characteristics
    size = adios_get_type_size (v->type, v->data);

    switch (v->type)
    {
        case adios_byte:
        case adios_unsigned_byte:
        case adios_short:
        case adios_unsigned_short:
        case adios_integer:
        case adios_unsigned_integer:
        case adios_long:
        case adios_unsigned_long:
        case adios_real:
        case adios_double:
        case adios_long_double:
        case adios_complex:
        case adios_double_complex:
            if (v->dimensions)
            {
                // add a dimensions characteristic
                characteristic_set_count++;
                flag = (uint8_t) adios_characteristic_dimensions;
                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                             ,&flag, 1
                             );
                index_size += 1;
                characteristic_set_length += 1;

                len = adios_write_var_characteristics_dims_v2 (fd, v);
                index_size += len;
                characteristic_set_length += len;

                // add a min value characteristic
                characteristic_set_count++;
                flag = (uint8_t) adios_characteristic_min;
                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                             ,&flag, 1
                             );
                index_size += 1;
                characteristic_set_length += 1;
                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                             ,v->min, size
                             );
                index_size += size;
                characteristic_set_length += size;

                // add a max value characteristic
                characteristic_set_count++;
                flag = (uint8_t) adios_characteristic_max;
                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                             ,&flag, 1
                             );
                index_size += 1;
                characteristic_set_length += 1;
                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                             ,v->max, size
                             );
                index_size += size;
                characteristic_set_length += size;
            }
            break;

        case adios_string:
            break;
    }
    // characteristics count/size prefix
    buffer_write (&fd->buffer, &fd->buffer_size, &characteristic_set_start
                 ,&characteristic_set_count, 1
                 );
    buffer_write (&fd->buffer, &fd->buffer_size, &characteristic_set_start
                 ,&characteristic_set_length, 4
                 );


    return index_size;
}

int adios_generate_var_characteristics_v2 (struct adios_file_struct * fd
                                    ,struct adios_var_struct * var
                                    )
{
    uint64_t total_size = adios_get_var_size (var, fd->group, var->data);
    uint64_t size = 0;

#if 1
#define MIN_MAX(a,b) \
{\
a * data = (a *) var->data; \
var->min = malloc (b); \
var->max = malloc (b); \
a * min = (a *) var->min; \
a * max = (a *) var->max; \
*min = data [0]; \
*max = data [0]; \
size++; \
while ((size * b) < total_size) \
{ \
    if (data [size] < *min) \
        *min = data [size]; \
    if (data [size] > *max) \
        *max = data [size]; \
    size++; \
} \
return 0; \
}
#else
#define MIN_MAX(a,b)\
{\
a * data = (a *) var->data; \
var->min = malloc (b); \
var->max = malloc (b); \
a * min = (a *) var->min; \
a * max = (a *) var->max; \
*min = data [0]; \
*max = data [0]; \
return 0; \
}
#endif

    switch (var->type)
    {
        case adios_byte:
            MIN_MAX(int8_t,1)

        case adios_unsigned_byte:
            MIN_MAX(uint8_t,1)

        case adios_short:
            MIN_MAX(int16_t,2)

        case adios_unsigned_short:
            MIN_MAX(uint16_t,2)

        case adios_integer:
            MIN_MAX(int32_t,4)

        case adios_unsigned_integer:
            MIN_MAX(uint32_t,4)

        case adios_long:
            MIN_MAX(int64_t,8)

        case adios_unsigned_long:
            MIN_MAX(uint64_t,8)

        case adios_real:
            MIN_MAX(float,4)

        case adios_double:
            MIN_MAX(double,8)

        case adios_long_double:
            MIN_MAX(long double,16)

        case adios_complex:
	{
            // complex = (float, float) = double in size
            float   minr, mini, maxr, maxi, ar, ai; 
            double  min, max, a;
            float    *data = (float *) var->data;   // view double array as float array 
            uint64_t step = 0, steps = total_size / 4;  // steps in float steps over double array!
            var->min = malloc ( 2*sizeof(float));
            var->max = malloc ( 2*sizeof(float));
            minr = data[0];
            mini = data[1];
            min  = (double)minr*minr+mini*mini;
            max  = (double)maxr*maxr+maxi*maxi;
            maxr = minr;
            maxi = maxi;
            step += 2; 
            // loop over the elements of the complex array (step is wrt float size!)
            while (step  < steps) {
                ar = data[step];
                ai = data[step+1];
                a  = (double)ar*ar+ai*ai;
                 
                if ( a < min) {
                    minr = ar; mini = ai; min  = a;
                }
                if ( a > max) {
                    maxr = ar; maxi = ai; max  = a;
                }
                step += 2; 
            } 
            ((float *) var->min)[0] = minr;
            ((float *) var->min)[1] = mini;
            ((float *) var->max)[0] = maxr;
            ((float *) var->max)[1] = maxi;
            return 0;
	}

        case adios_double_complex:
	{
            // complex = (double, double) = 16 bytes in size
            float   minr, mini, maxr, maxi, ar, ai; 
            long double  min, max, a; 
            /* FIXME: long double may not prevent overflow 
               PGI compiler: long double is 8 bytes like double
            */
            double   *data = (double *) var->data;   // view 16-byte array as double array 
            uint64_t step = 0, steps = total_size / 8;  // steps in double steps over 16-byte array!
            var->min = malloc ( 2*sizeof(double));
            var->max = malloc ( 2*sizeof(double));
            minr = data[0];
            mini = data[1];
            min  = (long double)minr*minr+mini*mini;
            max  = (long double)maxr*maxr+maxi*maxi;
            maxr = minr;
            maxi = maxi;
            step += 2; 
            // loop over the elements of the double complex array (step is wrt double size!)
            while (step  < steps) {
                ar = data[step];
                ai = data[step+1];
                a  = (long double)ar*ar+ai*ai;
                 
                if ( a < min) {
                    minr = ar; mini = ai; min  = a;
                }
                if ( a > max) {
                    maxr = ar; maxi = ai; max  = a;
                }
                step += 2; 
            } 
            ((double *) var->min)[0] = minr;
            ((double *) var->min)[1] = mini;
            ((double *) var->max)[0] = maxr;
            ((double *) var->max)[1] = maxi;
            return 0;

            return 0;
	}

        case adios_string:
        {
            var->min = 0;
            var->max = 0;

            return 0;
        }

        default:
	{
            uint64_t data = adios_unknown;
            var->min = 0;
            var->max = 0;
            return 0;
	}
    }
}

// data is only there for sizing
uint64_t adios_write_var_header_v2 (struct adios_file_struct * fd
                                   ,struct adios_var_struct * v
                                   )
{
    uint64_t total_size = 0;
    uint8_t flag;
    uint16_t len;

    uint64_t start = fd->offset;  // save to write the size
    v->write_offset = fd->offset + fd->base_offset; // save offset in file
    fd->offset += 8;              // save space for the size
    total_size += 8;              // makes final parsing easier

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &v->id, 2);
    total_size += 2;

    len = strlen (v->name);
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &len, 2);
    total_size += 2;

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, v->name, len);
    total_size += len;

    len = strlen (v->path);
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &len, 2);
    total_size += 2;

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, v->path, len);
    total_size += len;

    flag = v->type;
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &flag, 1);
    total_size += 1;

    flag = (v->is_dim == adios_flag_yes ? 'y' : 'n');
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &flag, 1);
    total_size += 1;

    total_size += adios_write_dimensions_v2 (fd, v->dimensions);

    adios_generate_var_characteristics_v2 (fd, v);
    total_size += adios_write_var_characteristics_v2 (fd, v);

    total_size += adios_get_var_size (v, fd->group, v->data); // payload

    buffer_write (&fd->buffer, &fd->buffer_size, &start, &total_size, 8);

    fd->vars_written++;

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return total_size;
}

int adios_write_var_payload_v2 (struct adios_file_struct * fd
                               ,struct adios_var_struct * var
                               )
{
    uint64_t size;

    // write payload
    size = adios_get_var_size (var, fd->group, var->data);
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, var->data, size);

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return 0;
}

int adios_write_attribute_v2 (struct adios_file_struct * fd
                             ,struct adios_attribute_struct * a
                             )
{
    uint64_t start;        // save the start to write the size
    uint32_t size = 0;
    uint16_t len = 0;
    uint8_t flag = 0;

    // save space for attr length
    start = fd->offset;
    a->write_offset = fd->offset + fd->base_offset; // save offset in file
    fd->offset += 4;

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &a->id, 2);
    size += 2;

    len = strlen (a->name);
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &len, 2);
    size += 2;

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, a->name, len);
    size += len;

    len = strlen (a->path);
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &len, 2);
    size += 2;

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, a->path, len);
    size += len;

    flag = (a->var ? 'y' : 'n');
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &flag, 1);
    size += 1;

    if (a->var)
    {
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,&a->var->id, 2
                     );
        size += 2;
    }
    else
    {
        flag = a->type;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &flag, 1);
        size += 1;

        uint32_t t = adios_get_type_size (a->type, a->value);
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &t, 4);
        size += 4;

        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                     ,a->value, t
                     );
        size += t;
    }

    // put in the size we have put in for this attribute
    buffer_write (&fd->buffer, &fd->buffer_size, &start, &size, 4);

    fd->vars_written++;

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return 0;
}

int adios_write_open_vars_v2 (struct adios_file_struct * fd)
{
    fd->vars_written = 0;

    // it is now setup to write the vars and then the attrs on close
    fd->vars_start = fd->offset;

    fd->offset += (2 + 8); // (count + size)

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return 0;
}

int adios_write_close_vars_v2 (struct adios_file_struct * fd)
{
    // close the var area (count and total size) and write the attributes
    uint64_t size = fd->offset - fd->vars_start;
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->vars_start, &fd->vars_written, 2);

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->vars_start, &size, 8);

    return 0;
}

int adios_write_open_attributes_v2 (struct adios_file_struct * fd)
{
    fd->vars_start = fd->offset;   // save the start of attr area for size
    fd->offset += (2 + 8);         // space to write the count and size
    fd->vars_written = 0;

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return 0;
}

int adios_write_close_attributes_v2 (struct adios_file_struct * fd)
{
    // write attribute count and total size
    uint64_t size = fd->offset - fd->vars_start;
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->vars_start, &fd->vars_written, 2);

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->vars_start, &size, 8);

    return 0;
}

// *****************************************************************************

static int adios_multiply_dimensions (uint64_t * size
                                     ,struct adios_var_struct * var
                                     ,enum ADIOS_DATATYPES type
                                     ,void * data
                                     )
{
    switch (type)
    {
        case adios_unsigned_byte:
            *size *= (*(uint8_t *) data);
            return 1;

        case adios_byte:
            *size *= (*(int8_t *) data);
            return 1;

        case adios_unsigned_short:
            *size *= (*(uint16_t *) data);
            return 1;

        case adios_short:
            *size *= (*(int16_t *) data);
            return 1;

        case adios_unsigned_integer:
            *size *= (*(uint32_t *) data);
            return 1;

        case adios_integer:
            *size *= (*(int32_t *) data);
            return 1;

        case adios_unsigned_long:
            *size *= (*(uint64_t *) data);
            return 1;

        case adios_long:
            *size *= (*(int64_t *) data);
            return 1;

        default:
            fprintf (stderr, "Invalid datatype for array dimension on "
                             "var %s: %s\n"
                    ,var->name
                    ,adios_type_to_string_int (type)
                    );

            return 0;
    }
}
