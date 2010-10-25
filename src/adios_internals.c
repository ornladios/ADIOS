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
#include "adios_bp_v1.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

extern struct adios_method_list_struct * adios_methods;
extern struct adios_group_list_struct * adios_groups;

int adios_int_is_var (const char * temp) // 1 == yes, 0 == no
{
    if (!temp)
        return 1;

    if (*temp == '-' || isdigit (*temp))
    {
        while (*temp)
        {
            if (isdigit (*temp))
                temp++;
            else
                return 1;
        }
    }
    else
        return 1;

    return 0;
}

int adios_int_is_num (char * temp) // 1 == yes, 0 == no
{
    char * extra = 0;

    strtod (temp, &extra);

    if (extra)
        return 0;
    else
        return 1;

    return 0;
}

struct adios_var_struct * adios_find_var_by_name (struct adios_var_struct * root
                                                 ,const char * name
                                                 ,enum ADIOS_FLAG unique_names
                                                 )
{
    int done = 0;
    struct adios_var_struct * var = 0;

    if (!name)
    {
        done = 1;
        root = 0;
    }

    while (!done && root)
    {
        char * compare_name = root->name;
        char * compare_name_path = root->name;
        if (unique_names == adios_flag_no)
        {
            compare_name_path = malloc (  strlen (root->name)
                                        + strlen (root->path)
                                        + 2 // null term and '/'
                                       );
            if (!strcmp (root->path, "/"))
                sprintf (compare_name_path, "/%s", root->name);
            else
                sprintf (compare_name_path, "%s/%s", root->path, root->name);
        }

        if (   !strcasecmp (name, compare_name)
            || (   unique_names == adios_flag_no
                && !strcasecmp (name, compare_name_path)
               )
           )
        {
            done = 1;
            var = root;
        }
        else
        {
            root = root->next;
        }

        if (unique_names == adios_flag_no)
        {
            free (compare_name_path);
        }
    }

    return var;
}

struct adios_attribute_struct * adios_find_attribute_by_name
                                        (struct adios_attribute_struct * root
                                        ,const char * name
                                        ,enum ADIOS_FLAG unique_names
                                        )
{
    int done = 0;
    struct adios_attribute_struct * attr = 0;

    if (!name)
    {
        done = 1;
        root = 0;
    }

    while (!done && root)
    {
        char * compare_name = root->name;
        char * compare_name_path = root->name;
        if (unique_names == adios_flag_no)
        {
            compare_name_path = malloc (  strlen (root->name)
                                        + strlen (root->path)
                                        + 2 // null term and '/'
                                       );
            if (!strcmp (root->path, "/"))
                sprintf (compare_name_path, "/%s", root->name);
            else
                sprintf (compare_name_path, "%s/%s", root->path, root->name);
        }

        if (   !strcasecmp (name, compare_name)
            || (   unique_names == adios_flag_no
                && !strcasecmp (name, compare_name_path)
               )
           )
        {
            done = 1;
            attr = root;
        }
        else
        {
            root = root->next;
        }

        if (unique_names == adios_flag_no)
        {
            free (compare_name_path);
        }
    }

    return attr;
}

struct adios_var_struct * adios_find_var_by_id (struct adios_var_struct * root
                                               ,uint16_t id
                                               )
{
    while (root)
    {
        if (root->id == id)
            return root;
        else
            root = root->next;
    }

    return NULL;
}

struct adios_attribute_struct * adios_find_attribute_by_id
                                         (struct adios_attribute_struct * root
                                         ,uint16_t id
                                         )
{
    while (root)
    {
        if (root->id == id)
            return root;
        else
            root = root->next;
    }

    return NULL;
}

int adios_parse_dimension (const char * dimension
                          ,const char * global_dimension
                          ,const char * local_offset
                          ,struct adios_group_struct * g
                          ,struct adios_dimension_struct * dim
                          )
{
    if (!dimension)
    {
        fprintf (stderr, "adios_parse_dimension: dimension not provided\n");

        return 0;
    }

    dim->dimension.rank = 0;
    dim->dimension.id = 0;
    dim->dimension.time_index = adios_flag_no;
    if (adios_int_is_var (dimension))
    {
        struct adios_var_struct * var = 0;
        dim->dimension.rank = 0;
        var = adios_find_var_by_name (g->vars, dimension
                                     ,g->all_unique_var_names
                                     );
        if (!var)
        {
            struct adios_attribute_struct * attr = 0;
            attr = adios_find_attribute_by_name (g->attributes, dimension
                                                ,g->all_unique_var_names
                                                );

            if (!attr)
            {
                if (   g->time_index_name
                    && !strcasecmp (g->time_index_name, dimension)
                   )
                {
                    dim->dimension.time_index = adios_flag_yes;
                }
                else
                {
                    fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                            ,dimension
                            );

                    return 0;
                }
            }
            else
            {
                if (attr->var)
                {
                    switch (attr->var->type)
                    {
                        case adios_string:
                        case adios_real:
                        case adios_double:
                        case adios_long_double:
                        case adios_complex:
                        case adios_double_complex:
                            fprintf (stderr, "config.xml: var dimension %s "
                                             "has an invalid type: %s\n"
                                    ,attr->name
                                    ,adios_type_to_string_int (attr->var->type)
                                    );
                            return 0;

                        default: // the integral numeric types are all fine
                            break;
                    }
                    attr->var->is_dim = adios_flag_yes;
                }
                else
                {
                    switch (attr->type)
                    {
                        case adios_string:
                        case adios_real:
                        case adios_double:
                        case adios_long_double:
                        case adios_complex:
                        case adios_double_complex:
                            fprintf (stderr, "config.xml: var dimension %s "
                                             "has an invalid type: %s\n"
                                    ,attr->name
                                    ,adios_type_to_string_int (attr->type)
                                    );
                            return 0;

                        default: // the integral numeric types are all fine
                            break;
                    }
                }
                dim->dimension.id = attr->id;
            }
        }
        else
        {
            switch (var->type)
            {
                case adios_string:
                case adios_real:
                case adios_double:
                case adios_long_double:
                case adios_complex:
                case adios_double_complex:
                    fprintf (stderr, "config.xml: var dimension %s "
                                     "has an invalid type: %s\n"
                            ,var->name
                            ,adios_type_to_string_int (var->type)
                            );
                    return 0;

                default: // the integral numeric types are all fine
                    break;
            }

            dim->dimension.id = var->id;
            var->is_dim = adios_flag_yes;
        }
    }
    else
    {
        dim->dimension.id = 0;
        dim->dimension.rank = atoi (dimension);
    }

    if (!global_dimension)
    {
        fprintf (stderr, "adios_parse_dimension: global_dimension not "
                         "provided\n"
                );

        return 0;
    }

    if (adios_int_is_var (global_dimension))
    {
        struct adios_var_struct * var = 0;
        dim->global_dimension.rank = 0;
        var = adios_find_var_by_name (g->vars, global_dimension
                                     ,g->all_unique_var_names
                                     );
        if (!var)
        {
            struct adios_attribute_struct * attr = 0;
            attr = adios_find_attribute_by_name (g->attributes, global_dimension
                                                ,g->all_unique_var_names
                                                );

            if (!attr)
            {
                if (   g->time_index_name
                    && !strcasecmp (g->time_index_name, global_dimension)
                   )
                {
                    dim->global_dimension.time_index = adios_flag_yes;
                }
                else
                {
                    fprintf (stderr, "config.xml: invalid global-bounds "
                                     "dimension: %s\n"
                            ,global_dimension
                            );

                    return 0;
                }
            }
            else
            {
                if (attr->var)
                {
                    switch (attr->var->type)
                    {
                        case adios_string:
                        case adios_real:
                        case adios_double:
                        case adios_long_double:
                        case adios_complex:
                        case adios_double_complex:
                            fprintf (stderr, "config.xml: var dimension %s "
                                             "has an invalid type: %s\n"
                                    ,attr->name
                                    ,adios_type_to_string_int (attr->var->type)
                                    );
                            return 0;

                        default: // the integral numeric types are all fine
                            break;
                    }
                    attr->var->is_dim = adios_flag_yes;
                }
                else
                {
                    switch (attr->type)
                    {
                        case adios_string:
                        case adios_real:
                        case adios_double:
                        case adios_long_double:
                        case adios_complex:
                        case adios_double_complex:
                            fprintf (stderr, "config.xml: var dimension %s "
                                             "has an invalid type: %s\n"
                                    ,attr->name
                                    ,adios_type_to_string_int (attr->type)
                                    );
                            return 0;

                        default: // the integral numeric types are all fine
                            break;
                    }
                }
                dim->global_dimension.id = attr->id;
            }
        }
        else
        {
            switch (var->type)
            {
                case adios_string:
                case adios_real:
                case adios_double:
                case adios_long_double:
                case adios_complex:
                case adios_double_complex:
                    fprintf (stderr, "config.xml: var dimension %s "
                                     "has an invalid type: %s\n"
                            ,var->name
                            ,adios_type_to_string_int (var->type)
                            );
                    return 0;

                default: // the integral numeric types are all fine
                    break;
            }
            var->is_dim = adios_flag_yes;
            dim->global_dimension.id = var->id;
        }
    }
    else
    {
        dim->global_dimension.id = 0;
        dim->global_dimension.rank = strtol (global_dimension, NULL, 10);
    }

    if (!local_offset)
    {
        fprintf (stderr, "adios_parse_dimension: local-offset not provided\n");

        return 0;
    }

    if (adios_int_is_var (local_offset))
    {
        struct adios_var_struct * var = 0;
        dim->local_offset.rank = 0;
        var = adios_find_var_by_name (g->vars, local_offset
                                     ,g->all_unique_var_names
                                     );
        if (!var)
        {
            struct adios_attribute_struct * attr = 0;
            attr = adios_find_attribute_by_name (g->attributes, local_offset
                                                ,g->all_unique_var_names
                                                );

            if (!attr)
            {
                if (   g->time_index_name
                    && !strcasecmp (g->time_index_name, local_offset)
                   )
                {
                    dim->local_offset.time_index = adios_flag_yes;
                }
                else
                {
                    fprintf (stderr, "config.xml: invalid var local_offset: "
                                     "%s\n"
                            ,local_offset
                            );

                    return 0;
                }
            }
            else
            {
                if (attr->var)
                {
                    switch (attr->var->type)
                    {
                        case adios_string:
                        case adios_real:
                        case adios_double:
                        case adios_long_double:
                        case adios_complex:
                        case adios_double_complex:
                            fprintf (stderr, "config.xml: var dimension %s "
                                             "has an invalid type: %s\n"
                                    ,attr->name
                                    ,adios_type_to_string_int (attr->var->type)
                                    );
                            return 0;

                        default: // the integral numeric types are all fine
                            break;
                    }
                    attr->var->is_dim = adios_flag_yes;
                }
                else
                {
                    switch (attr->type)
                    {
                        case adios_string:
                        case adios_real:
                        case adios_double:
                        case adios_long_double:
                        case adios_complex:
                        case adios_double_complex:
                            fprintf (stderr, "config.xml: var dimension %s "
                                             "has an invalid type: %s\n"
                                    ,attr->name
                                    ,adios_type_to_string_int (attr->type)
                                    );
                            return 0;

                        default: // the integral numeric types are all fine
                            break;
                    }
                }
                dim->local_offset.id = attr->id;
            }
        }
        else
        {
            switch (var->type)
            {
                case adios_string:
                case adios_real:
                case adios_double:
                case adios_long_double:
                case adios_complex:
                case adios_double_complex:
                    fprintf (stderr, "config.xml: var dimension %s "
                                     "has an invalid type: %s\n"
                            ,var->name
                            ,adios_type_to_string_int (var->type)
                            );
                    return 0;

                default: // the integral numeric types are all fine
                    break;
            }
            var->is_dim = adios_flag_yes;
            dim->local_offset.id = var->id;
        }
    }
    else
    {
        dim->local_offset.id = 0;
        dim->local_offset.rank = strtol (local_offset, NULL, 10);
    }

    return 1;
}

struct adios_method_list_struct * adios_get_methods ()
{
    return adios_methods;
}

struct adios_group_list_struct * adios_get_groups ()
{
    return adios_groups;
}


int adios_parse_scalar_string (enum ADIOS_DATATYPES type, char * value
                              ,void ** out
                              )
{
    char * end;

    switch (type)
    {
        case adios_byte:
        case adios_short:
        case adios_integer:
        {
            int errno_save = errno;
            long t = strtol (value, &end, 10);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "value: '%s' not valid integer\n"
                        ,value
                        );

                return 0;
            }
            else
            {
                switch (type)
                {
                    case adios_byte:
                        if (t < SCHAR_MIN || t > SCHAR_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string_int (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = malloc (1);
                            *((int8_t *) *out) = t;

                            return 1;
                        }
                    case adios_short:
                        if (t < SHRT_MIN || t > SHRT_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string_int (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = malloc (2);
                            *((int16_t *) *out) = t;

                            return 1;
                        }
                    case adios_integer:
                        if (t < INT_MIN || t > INT_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string_int (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = malloc (4);
                            *((int32_t *) *out) = t;

                            return 1;
                        }
                }
            }
        }
        case adios_long:
        {
            int errno_save = errno;
            int64_t t = strtoll (value, &end, 10);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "type is %s, value is out of range: '%s'\n"
                        ,adios_type_to_string_int (type), value
                        );

                return 0;
            }
            else
            {
                *out = malloc (8);
                *((int64_t *) *out) = t;

                return 1;
            }
        }
        case adios_unsigned_byte:
        case adios_unsigned_short:
        case adios_unsigned_integer:
        {
            int errno_save = errno;
            unsigned long t = strtoul (value, &end, 10);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "value: '%s' not valid integer\n"
                        ,value
                        );

                return 0;
            }
            else
            {
                switch (type)
                {
                    case adios_unsigned_byte:
                        if (t > UCHAR_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string_int (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = malloc (1);
                            *((uint8_t *) *out) = t;

                            return 1;
                        }
                    case adios_unsigned_short:
                        if (t > USHRT_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string_int (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = malloc (2);
                            *((uint16_t *) *out) = t;

                            return 1;
                        }
                    case adios_unsigned_integer:
                        if (t > UINT_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string_int (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = malloc (4);
                            *((uint32_t *) *out) = t;

                            return 1;
                        }
                }
            }
        }
        case adios_unsigned_long:
        {
            int errno_save = errno;
            uint64_t t = strtoull (value, &end, 10);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "type is %s, value is out of range: '%s'\n"
                        ,adios_type_to_string_int (type), value
                        );

                return 0;
            }
            else
            {
                *out = malloc (8);
                *((uint64_t *) *out) = t;

                return 1;
            }
        }
        case adios_real:
        {
            int errno_save = errno;
            float t = strtof (value, &end);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "type is %s, value is out of range: '%s'\n"
                        ,adios_type_to_string_int (type), value
                        );

                return 0;
            }
            else
            {
                *out = malloc (4);
                *((float *) *out) = t;

                return 1;
            }
        }
        case adios_double:
        {
            int errno_save = errno;
            double t = strtod (value, &end);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "type is %s, value is out of range: '%s'\n"
                        ,adios_type_to_string_int (type), value
                        );

                return 0;
            }
            else
            {
                *out = malloc (8);
                *((double *) *out) = t;

                return 1;
            }
        }
        case adios_long_double:
        {
            int errno_save = errno;
            long double t = strtold (value, &end);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "type is %s, value is out of range: '%s'\n"
                        ,adios_type_to_string_int (type), value
                        );

                return 0;
            }
            else
            {
                *out = malloc (16);
                *((long double *) *out) = t;
            }
        }
        case adios_string:
        {
            *out = (void *) strdup (value);

            return 1;
        }
        case adios_complex:
        {
            fprintf (stderr, "adios_complex type validation needs to be "
                             "implemented\n");
            return 1;
        }
        case adios_double_complex:
        {
            fprintf (stderr, "adios_double_complex type validation needs to "
                             "be implemented\n");
            return 1;
        }

        case adios_unknown:
        default:
            fprintf (stderr, "unknown type cannot be validated\n");

            return 0;
    }

    return 1;
}

int adios_common_define_attribute (int64_t group, const char * name
                                  ,const char * path
                                  ,enum ADIOS_DATATYPES type
                                  ,const char * value
                                  ,const char * var
                                  )
{
    struct adios_group_struct * g = (struct adios_group_struct *) group;
    struct adios_attribute_struct * attr = (struct adios_attribute_struct *)
                              malloc (sizeof (struct adios_attribute_struct));

    attr->name = strdup (name);
    attr->path = strdup (path);
    if (value)
    {
        if (type == adios_unknown)
        {
            fprintf (stderr, "config.xml: attribute element %s has invalid "
                             "type attribute\n"
                    ,name
                    );

            free (attr->name);
            free (attr->path);
            free (attr);

            return 0;
        }
        attr->type = type;
        if (adios_parse_scalar_string (type, (void *) value, &attr->value))
        {
            attr->var = 0;
        }
        else
        {
            fprintf (stderr, "config.xml: attribute element %s has invalid "
                             "value attribute: '%s'\n"
                    ,name, value
                    );

            free (attr->value);
            free (attr->name);
            free (attr->path);
            free (attr);

            return 0;
        }
    }
    else
    {
        attr->value = 0;
        attr->type = adios_unknown;
        attr->var = adios_find_var_by_name (g->vars, var
                                           ,g->all_unique_var_names
                                           );

        if (attr->var == 0)
        {
            fprintf (stderr, "config.xml: attribute element %s references "
                             "var %s that has not been defined.\n"
                    ,name, var
                    );

            free (attr->name);
            free (attr->path);
            free (attr);

            return 0;
        }
    }

    attr->next = 0;

    adios_append_attribute (&g->attributes, attr, ++g->member_count);

    return 1;
}

/*void adios_extract_string (char ** out, const char * in, int size)
{
    if (in && out)
    {
        *out = malloc (strlen (in) + 1);
        strcpy (*out, in);
    }
// for some Fortran implementations, we get a size for a string.
// for others (like PGI), we don't and it isn't null terminated
// unless we do it explicitly.  Assume that it is null terminated
// for now.
//
#if 0
    int i = 0;
    memcpy (out, in, size);
    while (i < size)
    {
        if (out [i] == ' ')
        {
            out [i] = 0;
            return;
        }
        else
            i++;
    }
    out [i] = 0;
#endif
}*/

void adios_append_method (struct adios_method_struct * method)
{
    struct adios_method_list_struct ** root = &adios_methods;

    while (root)
    {
        if (!*root)
        {
            struct adios_method_list_struct * new_node =
                 (struct adios_method_list_struct *)
                   malloc (sizeof (struct adios_method_list_struct));

            if (!new_node)
            {
                fprintf (stderr, "out of memory in adios_append_method\n");
            }
            new_node->method = method;
            new_node->next = 0;

            *root = new_node;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}

void adios_add_method_to_group (struct adios_method_list_struct ** root
                               ,struct adios_method_struct * method
                               )
{
    while (root)
    {
        if (!*root)
        {
            struct adios_method_list_struct * new_node =
                 (struct adios_method_list_struct *)
                   malloc (sizeof (struct adios_method_list_struct));

            if (!new_node)
            {
                fprintf (stderr, "out of memory in adios_append_method\n");
            }
            new_node->method = method;
            new_node->next = 0;

            *root = new_node;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}

void adios_append_group (struct adios_group_struct * group)
{
    struct adios_group_list_struct ** root = &adios_groups;
    int id = 1;

    while (root)
    {
        if (!*root)
        {
            struct adios_group_list_struct * new_node =
                 (struct adios_group_list_struct *)
                   malloc (sizeof (struct adios_group_list_struct));

            if (!new_node)
            {
                fprintf (stderr, "out of memory in adios_append_group\n");
            }
            group->id = id;
            new_node->group = group;
            new_node->next = 0;

            *root = new_node;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
            id++;
        }
    }
}

// return is whether or not the name is unique
enum ADIOS_FLAG adios_append_var (struct adios_var_struct ** root
                                 ,struct adios_var_struct * var
                                 ,uint16_t id
                                 )
{
    enum ADIOS_FLAG unique_names = adios_flag_yes;

    while (root)
    {
        if (   unique_names == adios_flag_yes
            && *root
            && !strcasecmp ((*root)->name, var->name)
           )
        {
            unique_names = adios_flag_no;
        }
        if (!*root)
        {
            var->id = id;
            *root = var;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }

    return unique_names;
}

void adios_append_dimension (struct adios_dimension_struct ** root
                            ,struct adios_dimension_struct * dimension
                            )
{
    while (root)
    {
        if (!*root)
        {
            *root = dimension;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}

void adios_append_attribute (struct adios_attribute_struct ** root
                            ,struct adios_attribute_struct * attribute
                            ,uint16_t id
                            )
{
    while (root)
    {
        if (!*root)
        {
            attribute->id = id;
            *root = attribute;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
// functions to support C & Fortran interface
///////////////////////////////////////////////////////////////////////////////
int adios_common_declare_group (int64_t * id, const char * name
                               ,enum ADIOS_FLAG host_language_fortran
                               ,const char * coordination_comm
                               ,const char * coordination_var
                               ,const char * time_index_name
                               ,enum ADIOS_FLAG stats
                               )
{
    struct adios_group_struct * g = (struct adios_group_struct *)
                             malloc (sizeof (struct adios_group_struct));

    g->name = strdup (name);
    g->adios_host_language_fortran = host_language_fortran;
    g->all_unique_var_names = adios_flag_yes;
    g->id = 0; // will be set in adios_append_group
    g->member_count = 0; // will be set in adios_append_group
    g->var_count = 0;
    g->vars = 0;
    g->vars_written = 0;
    g->attributes = 0;
    g->group_by = (coordination_var ? strdup (coordination_var) : 0L);
    g->group_comm = (coordination_comm ? strdup (coordination_comm) : 0L);
    g->time_index_name = (time_index_name ? strdup (time_index_name) : 0L);
    g->time_index = 0;
    g->stats_on = stats;
    g->process_id = 0;
    g->methods = 0;
    g->mesh = 0;

    *id = (int64_t) g;

    adios_append_group (g);

    return 1;
}

int adios_common_free_group (int64_t id)
{
    struct adios_group_list_struct * root = adios_groups;
    struct adios_group_list_struct * old_root = adios_groups;
    struct adios_group_struct * g = (struct adios_group_struct *) id;

    if (!root)
    {
        fprintf (stderr, "Err in adios_common_free_group()\n");
        return -1;
    }
    while (root && root->group->id != g->id)
    {
        old_root = root;
        root = root->next;
    };

    if (!root)
    {
        // Didn't find the group
        fprintf (stderr, "Err in adios_common_free_group()\n");
        return -1;
    }

    // old_root->root->root next
    if (old_root == adios_groups)
    {
        adios_groups =  root->next;
    }
    else
    {
        old_root->next = root->next;
    }

    if (g->name)
        free (g->name);

    while (g->vars)
    {
        struct adios_var_struct * vars = g->vars->next;

        if (g->vars->name)
            free (g->vars->name);
        if (g->vars->path)
            free (g->vars->path);

        while (g->vars->dimensions)
        {
            struct adios_dimension_struct * dimensions
                            = g->vars->dimensions->next;

            free (g->vars->dimensions);
            g->vars->dimensions = dimensions;
        }

        // NCSU - Clear Stat
        if (g->vars->stats)
		{
            uint8_t j = 0, idx = 0;
            uint8_t c = 0, count = adios_get_stat_set_count(g->vars->type);

            for (c = 0; c < count; c ++)
            {
                while (g->vars->bitmap >> j)
                {
                    if ((g->vars->bitmap >> j) & 1)
                    {
                        if (j == adios_statistic_hist)
                        {
                            struct adios_hist_struct * hist = (struct adios_hist_struct *) (g->vars->stats[c][idx].data);
                            free (hist->breaks);
                            free (hist->frequencies);
                            free (hist);
                        }
                        else
                            free (g->vars->stats[c][idx].data);

                        idx ++;
                    }
                    j ++;
                }
                free (g->vars->stats[c]);
            }
            free (g->vars->stats);		
		}

        if (g->vars->data)
            free (g->vars->data);

        free (g->vars);
        g->vars = vars;
    }

    free (root);

    return 0;
}

void trim_spaces (char * str)
{
    char * t = str, * p = NULL;
    while (*t != '\0')
    {
        if (*t == ' ')
        {
            p = t + 1;
            strcpy (t, p);
        }
        else
            t++;
    }

}

static void tokenize_dimensions (const char * str, char *** tokens, int * count)
{
    if (!str)
    {
        *tokens = 0;
        *count = 0;

        return;
    }

    char * save_str = strdup (str);
    char * t = save_str;
    int i;

    trim_spaces (save_str);

    if (strlen (save_str) > 0)
        *count = 1;
    else
    {
        *tokens = 0;
        *count = 0;
        free (save_str);

        return;
    }

    while (*t)
    {
        if (*t == ',')
            (*count)++;
        t++;
    }

    *tokens = (char **) malloc (sizeof (char **) * *count);
    (*tokens) [0] = strdup (strtok (save_str, ","));
    for (i = 1; i < *count; i++)
    {
        (*tokens) [i] = strdup (strtok (NULL, ","));
    }

    free (save_str);
}

static void cleanup_dimensions (char *** tokens, int * count)
{
    int i;
    for (i = 0; i < *count; i++)
    {
        free ((*tokens) [i]);
    }
    free (*tokens);
    *tokens = 0;
    *count = 0;
}

int adios_common_define_var_characteristics (struct adios_group_struct * g
                                            , const char * var_name
                                            , const char * bin_intervals
                                            , const char * bin_min
                                            , const char * bin_max
                                            , const char * bin_count
                                            )
{
    struct adios_var_struct * root = g->vars, * var;

    var = adios_find_var_by_name (root, var_name, adios_flag_no);

    struct adios_hist_struct * hist;

    if (var->type == adios_complex || var->type == adios_double_complex)
        return 0;

    int i = 0, j = 0;
    while ((var->bitmap >> j) && (j < adios_statistic_hist))
    {
        if ((var->bitmap >> j) & 1)
            i ++;
        j ++;
    }

    hist = var->stats[0][i].data = (struct adios_hist_struct *) malloc (sizeof(struct adios_hist_struct));

    if (!var)
    {
           fprintf (stderr, "config.xml: Didn't find the variable %s for analysis\n", var_name);
        return 0;
    }
    else
    {
        int i;
        if (bin_intervals)
        {
            int count;
            char ** bin_tokens = 0;

            tokenize_dimensions (bin_intervals, &bin_tokens, &count);

            if (!count)
            {
                fprintf (stderr, "config.xml: unable to tokenize break points\n");
                return 0;
            }

            hist->breaks = calloc(count, sizeof(double));

            if(!hist || !hist->breaks)
            {
                fprintf (stderr, "config.xml: unable to allocate memory for histogram break points in "
                        "adios_common_define_var_characteristics\n");
                return 0;
            }

            for(i = 0; i < count; i++)
            {
                hist->breaks[i] = atof(bin_tokens[i]);
                if(i > 0 && (hist->breaks[i] <= hist->breaks[i-1]))
                {
                    fprintf (stderr, "config.xml: break points should be in increasing order "
                        "adios_common_define_var_characteristics\n");
                    return 0;
                }
            }

            hist->num_breaks = count;
            hist->min = hist->breaks[0];

            if(count > 0)
                hist->max = hist->breaks[count - 1];
            else
                hist->max = hist->min;

            var->bitmap = var->bitmap | (1 << adios_statistic_hist);
        }
        else
        {
            if(!bin_max || !bin_min || !bin_count)
            {
                fprintf (stderr, "config.xml: unable to generate break points\n");
                return 0;
            }

            int count = atoi(bin_count);

            if (!count)
            {
                fprintf (stderr, "config.xml: bin count is undefined\n");
                return 0;
            }

            hist->num_breaks = count + 1;
            hist->min = atof(bin_min);
            hist->max = atof(bin_max);
            hist->breaks = calloc(hist->num_breaks, sizeof(double));

            if(!hist || !hist->breaks)
            {
                fprintf (stderr, "config.xml: unable to allocate memory for histogram break points in "
                        "adios_common_define_var_characteristics\n");
                return 0;
            }

            if (hist->min >= hist->max)
            {
                fprintf (stderr, "config.xml: minimum boundary value greater than maximum\n");
                return 0;
            }

            for(i = 0; i < hist->num_breaks; i ++)
                hist->breaks[i] = hist->min + i * (hist->max - hist->min) / count;

            var->bitmap = var->bitmap | (1 << adios_statistic_hist);
        }
    }


    return 1;
}

int adios_common_define_var (int64_t group_id, const char * name
                            ,const char * path, enum ADIOS_DATATYPES type
                            ,const char * dimensions
                            ,const char * global_dimensions
                            ,const char * local_offsets
                            )
{
    struct adios_group_struct * t = (struct adios_group_struct *) group_id;
    struct adios_var_struct * v = (struct adios_var_struct *)
                               malloc (sizeof (struct adios_var_struct));
    char * dim_temp;
    char * g_dim_temp;
    char * lo_dim_temp;
    enum ADIOS_FLAG flag;
    uint8_t i;
    if (dimensions)
        dim_temp = strdup (dimensions);
    else
        dim_temp = 0;
    if (global_dimensions)
        g_dim_temp = strdup (global_dimensions);
    else
        g_dim_temp = 0;
    if (local_offsets)
        lo_dim_temp = strdup (local_offsets);
    else
        lo_dim_temp = 0;

    v->name = strdup (name);
    v->path = strdup (path);
    v->type = type;
    v->dimensions = 0;
    v->is_dim = adios_flag_no;
    v->got_buffer = adios_flag_no;
    v->free_data = adios_flag_no;

    v->data = 0;
    v->write_offset = 0;

    v->data_size = 0;

    v->next = 0;

    // NCSU - Initializing stat related info
    v->stats = 0;
    v->bitmap = 0;

    // Q.L. - Check whether stats are disabled or not
    if (t->stats_on == adios_flag_yes)
    {
        // '1' at the bit location of stat id in adios_bp_v1.h, enables calculation of statistic.
        for (i = 0; i < ADIOS_STAT_LENGTH; i++)
            v->bitmap |= (1 << i);

        // Default values for histogram not yet implemented. Disabling it.
        v->bitmap ^= (1 << adios_statistic_hist);

        // For complex numbers, the set of statistics occur thrice: stat[0] - magnitude, stat[1] - real, stat[2] - imaginary
        if (v->type == adios_complex || v->type == adios_double_complex)
        {
            uint8_t c;
            v->stats = malloc (3 * sizeof(struct adios_stat_struct *));

            for (c = 0; c < 3; c ++)
                v->stats[c] = calloc (ADIOS_STAT_LENGTH, sizeof(struct adios_stat_struct));
        }
        else
        {
            v->stats = malloc (sizeof(struct adios_stat_struct *));
            v->stats[0] = calloc (ADIOS_STAT_LENGTH, sizeof(struct adios_stat_struct));
        }
    }

    // NCSU - End of initializing stat related info

    if (dim_temp && strcmp (dim_temp, ""))
    {
        int dim_count;
        char ** dim_tokens = 0;

        int g_dim_count;
        char ** g_dim_tokens = 0;

        int lo_dim_count;
        char ** lo_dim_tokens = 0;

        int i = 0;

        tokenize_dimensions (dim_temp, &dim_tokens, &dim_count);
        tokenize_dimensions (g_dim_temp, &g_dim_tokens, &g_dim_count);
        tokenize_dimensions (lo_dim_temp, &lo_dim_tokens, &lo_dim_count);

        while (i < dim_count)
        {
            int ret;
            struct adios_dimension_struct * d =
                     (struct adios_dimension_struct *)
                         calloc (1, sizeof (struct adios_dimension_struct));

            if (!d)
            {
                fprintf (stderr, "config.xml: out of memory in "
                                 "adios_common_define_var\n"
                        );

                return 0;
            }
            char * dim = 0;
            char * g_dim = "0";
            char * lo_dim = "0";

            if (i < dim_count)
                dim = dim_tokens [i];
            if (i < g_dim_count)
                g_dim = g_dim_tokens [i];
            if (i < lo_dim_count)
                lo_dim = lo_dim_tokens [i];

            if (!(ret = adios_parse_dimension (dim, g_dim, lo_dim, t, d)))
            {
                free (dim_temp);
                free (g_dim_temp);
                free (lo_dim_temp);
                free (v->name);
                free (v->path);
                free (v);
                cleanup_dimensions (&dim_tokens, &dim_count);
                cleanup_dimensions (&g_dim_tokens, &g_dim_count);
                cleanup_dimensions (&lo_dim_tokens, &lo_dim_count);

                return ret;
            }

            adios_append_dimension (&v->dimensions, d);

            i++;
        }
        cleanup_dimensions (&dim_tokens, &dim_count);
        cleanup_dimensions (&g_dim_tokens, &g_dim_count);
        cleanup_dimensions (&lo_dim_tokens, &lo_dim_count);
    }

    if (dim_temp)
        free (dim_temp);
    if (g_dim_temp)
        free (g_dim_temp);
    if (lo_dim_temp)
        free (lo_dim_temp);

    flag = adios_append_var (&t->vars, v, ++t->member_count);
    if (flag == adios_flag_no)
    {
        t->all_unique_var_names = adios_flag_no;
    }
    t->var_count++;

    return 1;
}

void adios_common_get_group (int64_t * group_id, const char * name)
{
    struct adios_group_list_struct * g = adios_get_groups ();

    *group_id = 0;

    while (g)
    {
        if (!strcasecmp (g->group->name, name))
        {
            *group_id = (int64_t) g->group;

            return;
        }

        g = g->next;
    }

    fprintf (stderr, "adios-group '%s' not found in configuration file\n"
            ,name
            );
}

// *****************************************************************************
static void buffer_write (char ** buffer, uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,const void * data, uint64_t size
                         )
{
    if (*buffer_offset + size > *buffer_size || *buffer == 0)
    {
        char * b = realloc (*buffer, *buffer_offset + size + 1000);
        if (b)
        {
            *buffer = b;
            *buffer_size = (*buffer_offset + size + 1000);
        }
        else
        {
            fprintf (stderr, "Cannot allocate memory in buffer_write.  "
                             "Requested: %llu\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}

static uint16_t adios_calc_var_characteristics_dims_overhead
                                                  (struct adios_var_struct * v)
{
    uint16_t overhead = 0;
    struct adios_dimension_struct * d = v->dimensions;

    overhead += 1; // count
    overhead += 2; // length

    while (d)
    {
        overhead += 8 + 8 + 8; // the dims

        d = d->next;
    }

    return overhead;
}

// NCSU -This function precomputes the amount of overhead consumed by statistics
uint16_t adios_calc_var_characteristics_stat_overhead (struct adios_var_struct * var)
{
    uint16_t i, j, overhead;
    overhead = j = i = 0;

    while (var->bitmap >> j)
    {
        // NCSU - This characteristic is present. It adds to the overhead
        if ((var->bitmap >> j) & 1)
            overhead += adios_get_stat_size(var->stats[0][i ++].data, var->type, j);
        j ++;
    }

    return overhead;
}


static uint16_t adios_calc_var_characteristics_overhead
                                                  (struct adios_var_struct * v)
{
    uint16_t overhead = 0;

    overhead += 1 + 4; // count + length

    switch (v->type)
    {
        case adios_string:   // nothing for strings
            //overhead += 1; // id
            //overhead += 2; // size
            break;

        default:   // the 12 numeric types
            if (v->dimensions)
            {
                overhead += 1; // id for bitmap
                overhead += 4; // value for bitmap

                overhead += 1;  // id for statistics
                // For complex numbers - min, max, avg repeated thrice
                overhead += adios_get_stat_set_count(v->type) * adios_calc_var_characteristics_stat_overhead (v);

                overhead += 1;  // id
                overhead += adios_calc_var_characteristics_dims_overhead (v);
            }
    }

    return overhead;
}

uint16_t adios_calc_var_overhead_v1 (struct adios_var_struct * v)
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

uint32_t adios_calc_attribute_overhead_v1 (struct adios_attribute_struct * a)
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

uint64_t adios_calc_overhead_v1 (struct adios_file_struct * fd)
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
        overhead += adios_calc_var_overhead_v1 (v);

        v = v->next;
    }

    overhead += 2; // attributes count
    overhead += 8; // attributes length

    while (a)
    {
        overhead += adios_calc_attribute_overhead_v1 (a);

        a = a->next;
    }

    return overhead;
}

int adios_write_process_group_header_v1 (struct adios_file_struct * fd
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

static void index_append_process_group_v1 (
                          struct adios_index_process_group_struct_v1 ** root
                         ,struct adios_index_process_group_struct_v1 * item
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

static void index_append_var_v1 (struct adios_index_var_struct_v1 ** root
                                ,struct adios_index_var_struct_v1 * item
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
                        * sizeof (struct adios_index_characteristic_struct_v1)
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
                        * sizeof (struct adios_index_characteristic_struct_v1)
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

static void index_append_attribute_v1
                                (struct adios_index_attribute_struct_v1 ** root
                                ,struct adios_index_attribute_struct_v1 * item
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
                       * sizeof (struct adios_index_characteristic_struct_v1)
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
                        * sizeof (struct adios_index_characteristic_struct_v1)
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
void adios_merge_index_v1 (struct adios_index_process_group_struct_v1 ** p1
                          ,struct adios_index_var_struct_v1 ** v1
                          ,struct adios_index_attribute_struct_v1 ** a1
                          ,struct adios_index_process_group_struct_v1 * p2
                          ,struct adios_index_var_struct_v1 * v2
                          ,struct adios_index_attribute_struct_v1 * a2
                          )
{
    // this will just add it on to the end and all should work fine
    index_append_process_group_v1 (p1, p2);

    // need to do vars attrs one at a time to merge them properly
    struct adios_index_var_struct_v1 * v_temp;
    struct adios_index_attribute_struct_v1 * a_temp;

    while (v2)
    {
        v_temp = v2->next;
        v2->next = 0;
        index_append_var_v1 (v1, v2);
        v2 = v_temp;
    }

    while (a2)
    {
        a_temp = a2->next;
        a2->next = 0;
        index_append_attribute_v1 (a1, a2);
        a2 = a_temp;
    }
}

// sort pg/var indexes by time index
void adios_sort_index_v1 (struct adios_index_process_group_struct_v1 ** p1
                         ,struct adios_index_var_struct_v1 ** v1
                         ,struct adios_index_attribute_struct_v1 ** a1
                         )
{
    struct adios_index_process_group_struct_v1 * p2 = 0, * p1_temp, * p2_temp, * p2_temp_prev;
    struct adios_index_var_struct_v1 * v2 = 0, * v1_temp, * v2_temp, * v2_temp_prev;
    struct adios_index_attribute_struct_v1 * a2 = 0, * a1_temp, * a2_temp, * a2_temp_prev;
    int i, j;

    while (*p1)
    {
        // if new index list is empty
        if (!p2)
        {
            p2 = *p1;
            *p1 = (*p1)->next;
            p2->next = 0;
        }
        else
        {
            p2_temp = p2;
            p2_temp_prev = p2;

            while (p2_temp && (*p1)->time_index >= p2_temp->time_index)
            {
                p2_temp_prev = p2_temp;
                p2_temp = p2_temp->next;
            }

            if (!p2_temp)
            {
                p2_temp_prev->next = *p1;
                *p1 = (*p1)->next;
                p2_temp_prev->next->next = 0;
            }
            else
            {
                p1_temp = (*p1)->next;
                (*p1)->next = p2_temp;
                p2_temp_prev->next = *p1;

                *p1 = p1_temp;
            }

        }
    }

    *p1 = p2;

    v1_temp = *v1;

    while (v1_temp)
    {
        for (i = 0; i < v1_temp->characteristics_count; i++)
        {
            for (j = 0; j < v1_temp->characteristics_count - i - 1; j++)
            {
                if (v1_temp->characteristics[j].time_index > v1_temp->characteristics[j + 1].time_index)
                {
                    uint64_t t_offset;  // beginning of the var or attr entry
                    struct adios_index_characteristic_dims_struct_v1 t_dims;
                    uint16_t t_var_id;
                    void * t_value;
                    uint64_t t_payload_offset;   // beginning of the var or attr payload
                    uint32_t t_file_index;
                    uint32_t t_time_index;

                    uint32_t t_bitmap;

                    struct adios_index_characteristics_stat_struct ** t_stats;

                    t_offset = v1_temp->characteristics[j].offset;
                    t_dims.count = v1_temp->characteristics[j].dims.count;
                    t_dims.dims = v1_temp->characteristics[j].dims.dims;
                    t_var_id = v1_temp->characteristics[j].var_id;
                    t_value = v1_temp->characteristics[j].value;
                    t_payload_offset = v1_temp->characteristics[j].payload_offset;
                    t_file_index = v1_temp->characteristics[j].file_index;
                    t_time_index = v1_temp->characteristics[j].time_index;
                    t_bitmap = v1_temp->characteristics[j].bitmap;
                    t_stats = v1_temp->characteristics[j].stats;
                    
                    v1_temp->characteristics[j].offset = v1_temp->characteristics[j + 1].offset;
                    v1_temp->characteristics[j].dims.count = v1_temp->characteristics[j + 1].dims.count;
                    v1_temp->characteristics[j].dims.dims = v1_temp->characteristics[j + 1].dims.dims;
                    v1_temp->characteristics[j].var_id = v1_temp->characteristics[j + 1].var_id;
                    v1_temp->characteristics[j].value = v1_temp->characteristics[j + 1].value;
                    v1_temp->characteristics[j].payload_offset = v1_temp->characteristics[j + 1].payload_offset;
                    v1_temp->characteristics[j].file_index = v1_temp->characteristics[j + 1].file_index;
                    v1_temp->characteristics[j].time_index = v1_temp->characteristics[j + 1].time_index;
                    v1_temp->characteristics[j].bitmap = v1_temp->characteristics[j + 1].bitmap;
                    v1_temp->characteristics[j].stats = v1_temp->characteristics[j + 1].stats;

                    v1_temp->characteristics[j + 1].offset = t_offset;
                    v1_temp->characteristics[j + 1].dims.count = t_dims.count;
                    v1_temp->characteristics[j + 1].dims.dims = t_dims.dims;
                    v1_temp->characteristics[j + 1].var_id = t_var_id;
                    v1_temp->characteristics[j + 1].value = t_value;
                    v1_temp->characteristics[j + 1].payload_offset = t_payload_offset;
                    v1_temp->characteristics[j + 1].file_index = t_file_index;
                    v1_temp->characteristics[j + 1].time_index = t_time_index;
                    v1_temp->characteristics[j + 1].bitmap = t_bitmap;
                    v1_temp->characteristics[j + 1].stats = t_stats;
                }
            }
        }

        v1_temp = v1_temp->next;
    }

    // no need to sort attributes
}

static void adios_clear_process_groups_index_v1 (
                            struct adios_index_process_group_struct_v1 * root
                           )
{
    while (root)
    {
        struct adios_index_process_group_struct_v1 * temp = root->next;
        if (root->group_name)
            free (root->group_name);
        if (root->time_index_name)
            free (root->time_index_name);
        free (root);
        root = temp;
    }
}

// NCSU - Clears up the statistical data from variable index table
static void adios_clear_vars_index_v1 (struct adios_index_var_struct_v1 * root)
{
    while (root)
    {
        int i;
        struct adios_index_var_struct_v1 * temp = root->next;

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
            if (root->characteristics [i].value)
                free (root->characteristics [i].value);

            // NCSU - Clears up the statistical data, based on bitmap
            if (root->characteristics [i].stats != 0)
            {
                uint8_t j = 0, idx = 0;
                uint8_t c = 0, count = adios_get_stat_set_count(root->type);

                for (c = 0; c < count; c ++)
                {
                    while (root->characteristics [i].bitmap >> j)
                    {
                        if ((root->characteristics [i].bitmap >> j) & 1)
                        {
							if (j == adios_statistic_hist)
							{
								struct adios_index_characteristics_hist_struct * hist = (struct adios_index_characteristics_hist_struct	*) root->characteristics [i].stats[c][idx].data;
								free (hist->breaks);
								free (hist->frequencies);
							}	
							else
                           		free (root->characteristics [i].stats[c][idx].data);
                            idx ++;
                        }
                        j ++;
                    }
                    free (root->characteristics [i].stats [c]);
                }

				free (root->characteristics [i].stats);
            }
        }
        if (root->characteristics)
            free (root->characteristics);

        free (root);
        root = temp;
    }
}

// NCSU - Clears up the statistical data, based on bitmap
static void adios_clear_attributes_index_v1
                                (struct adios_index_attribute_struct_v1 * root)
{
    while (root)
    {
        int i;
        struct adios_index_attribute_struct_v1 * temp = root->next;

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

            // NCSU - Clears up the statistical data, based on bitmap
            if (root->characteristics [i].stats != 0)
            {
                uint8_t j = 0, idx = 0;
                uint8_t c = 0, count = adios_get_stat_set_count(root->type);
                for (c = 0; c < count; c ++)
                {
                    while (root->characteristics [i].bitmap >> j)
                    {
                        if ((root->characteristics [i].bitmap >> j) & 1)
                        {
							if (j == adios_statistic_hist)
	                    	{
								struct adios_index_characteristics_hist_struct * hist = (struct adios_index_characteristics_hist_struct *) root->characteristics [i].stats[c][idx].data;
								free (hist->breaks);
								free (hist->frequencies);
								free (hist);
	                     	}
							else
								free (root->characteristics [i].stats[c][idx].data);

                            idx ++;
                        }
                        j ++;
                    }
                    free (root->characteristics [i].stats [c]);
                }
                free (root->characteristics [i].stats);
            }

            if (root->characteristics [i].value)
                free (root->characteristics [i].value);
        }
    
        if (root->characteristics)
            free (root->characteristics);

        free (root);
        root = temp;
    }
}


void adios_clear_index_v1 (struct adios_index_process_group_struct_v1 * pg_root
                          ,struct adios_index_var_struct_v1 * vars_root
                          ,struct adios_index_attribute_struct_v1 * attrs_root
                          )
{
    adios_clear_process_groups_index_v1 (pg_root);
    adios_clear_vars_index_v1 (vars_root);
    adios_clear_attributes_index_v1 (attrs_root);
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

void adios_copy_var_written (struct adios_var_struct ** root
                            ,struct adios_var_struct * var
                            ,struct adios_file_struct * fd
                            )
{
    struct adios_var_struct * var_new;

    while (root)
    {
        if (!*root)
        {
            var_new = (struct adios_var_struct *) malloc
                                                  (sizeof (struct adios_var_struct));
            //var_new->id = ++fd->group->member_count;
            var_new->id = var->id;
            var_new->parent_id = var->id;
            var_new->name = strdup (var->name);
            var_new->path = strdup (var->path);
            var_new->type = var->type;
            var_new->dimensions = 0;
            var_new->got_buffer = var->got_buffer;
            var_new->is_dim = var->is_dim;
            var_new->write_offset = var->write_offset;
            var_new->stats = 0;
            var_new->free_data = var->free_data;
            var_new->data = 0;
            var_new->data_size = var->data_size;
            var_new->next = 0;

            uint64_t size = adios_get_type_size (var->type, var->data);
            switch (var->type)
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
                    if (var->dimensions)
                    {
                        uint8_t c;
                        uint8_t j;
                        struct adios_dimension_struct * d = var->dimensions;
                        /*
                         *
                         * NOT ALL METHODS TRACK MIN/MAX.  CHECK BEFORE TRYING TO COPY.
                         *
                         */
                       // NCSU Statistics - copy stat to new var struct
                        uint8_t count = adios_get_stat_set_count(var->type);
                        uint8_t idx = 0;
                        uint64_t characteristic_size;

                        var_new->bitmap = var->bitmap;
                        var_new->stats = malloc (count * sizeof(struct adios_stat_struct *));

                        // Set of characteristics will be repeated thrice for complex numbers
                        for (c = 0; c < count; c ++)
                        {
                            var_new->stats[c] = calloc(ADIOS_STAT_LENGTH, sizeof (struct adios_stat_struct));

                            j = idx = 0;
                            while (var->bitmap >> j)
                             {
                                if ((var->bitmap >> j) & 1)
                                {
                                    if (var->stats[c][idx].data != NULL)
                                    {
                                        if (j == adios_statistic_hist)
                                        {
                                            var_new->stats[c][idx].data = (struct adios_hist_struct *) malloc (sizeof(struct adios_hist_struct));

                                            struct adios_hist_struct * var_hist = var->stats[c][idx].data;
                                            struct adios_hist_struct * var_new_hist = var_new->stats[c][idx].data;

                                            var_new_hist->min = var_hist->min;
                                            var_new_hist->max = var_hist->max;
                                            var_new_hist->num_breaks = var_hist->num_breaks;

                                            var_new_hist->frequencies = malloc ((var_hist->num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
                                            memcpy (var_new_hist->frequencies, var_hist->frequencies, (var_hist->num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
                                            var_new_hist->breaks = malloc ((var_hist->num_breaks) * adios_get_type_size(adios_double, ""));
                                            memcpy (var_new_hist->breaks, var_hist->breaks, (var_hist->num_breaks) * adios_get_type_size(adios_double, ""));
                                        }
                                        else
                                        {
                                            characteristic_size = adios_get_stat_size(var->stats[c][idx].data, var->type, j);
                                            var_new->stats[c][idx].data = malloc (characteristic_size);
                                            memcpy (var_new->stats[c][idx].data, var->stats[c][idx].data, characteristic_size);
                                        }

                                        idx ++;
                                    }
                                }
                                j ++;
                            }
                        }

                        c = count_dimensions (var->dimensions);

                        for (j = 0; j < c; j++)
                        {
                            struct adios_dimension_struct * d_new = (struct adios_dimension_struct *)
                                                            malloc (sizeof (struct adios_dimension_struct));
                            // de-reference dimension id
                            d_new->dimension.id = 0;
                            d_new->dimension.rank = get_value_for_dim (fd, &d->dimension);
                            d_new->dimension.time_index = d->dimension.time_index;
                            d_new->global_dimension.id = 0;
                            d_new->global_dimension.rank = get_value_for_dim (fd, &d->global_dimension);
                            d_new->global_dimension.time_index = d->global_dimension.time_index;
                            d_new->local_offset.id = 0;
                            d_new->local_offset.rank = get_value_for_dim (fd, &d->local_offset);
                            d_new->local_offset.time_index = d->local_offset.time_index;
                            d_new->next = 0;

                            adios_append_dimension (&var_new->dimensions, d_new);

                            d = d->next;
                        }
                    }
                    else
                    {
                        var_new->stats = 0;
                        var_new->data = malloc (size);
                        memcpy (var_new->data, var->data, size);
                    }

                    break;

                case adios_string:
                {
                    var_new->data = malloc (size + 1);
                    memcpy (var_new->data, var->data, size);
                    ((char *) (var_new->data)) [size] = 0;

                    break;
                }
            }

            *root = var_new;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}

void adios_build_index_v1 (struct adios_file_struct * fd
                          ,struct adios_index_process_group_struct_v1 ** pg_root
                          ,struct adios_index_var_struct_v1 ** vars_root
                          ,struct adios_index_attribute_struct_v1 ** attrs_root
                          )
{
    struct adios_group_struct * g = fd->group;
    struct adios_var_struct * v = g->vars_written;
    struct adios_attribute_struct * a = g->attributes;
    struct adios_index_process_group_struct_v1 * g_item;

    uint64_t process_group_count = 0;
    uint16_t var_count = 0;

    g_item = (struct adios_index_process_group_struct_v1 *)
                malloc (sizeof (struct adios_index_process_group_struct_v1));
    g_item->group_name = (g->name ? strdup (g->name) : 0L);
    g_item->adios_host_language_fortran = g->adios_host_language_fortran;
    g_item->process_id = g->process_id;
    g_item->time_index_name = (g->time_index_name ? strdup (g->time_index_name) : 0L);
    g_item->time_index = g->time_index;
    g_item->offset_in_file = fd->pg_start_in_file;
    g_item->next = 0;

    // build the groups and vars index
    index_append_process_group_v1 (pg_root, g_item);

    while (v)
    {
        // only add items that were written to the index
        if (v->write_offset != 0)
        {
            struct adios_index_var_struct_v1 * v_index;
            v_index = malloc (sizeof (struct adios_index_var_struct_v1));
            v_index->characteristics = malloc (
                           sizeof (struct adios_index_characteristic_struct_v1)
                          );

            v_index->id = v->id;
            v_index->group_name = (g->name ? strdup (g->name) : 0L);
            v_index->var_name = (v->name ? strdup (v->name) : 0L);
            v_index->var_path = (v->path ? strdup (v->path) : 0L);
            v_index->type = v->type;
            v_index->characteristics_count = 1;
            v_index->characteristics_allocated = 1;
            v_index->characteristics [0].offset = v->write_offset;
            // Find the old var in g->vars.
            // We need this to calculate the correct payload_offset, because that
            // holds the variable references in the dimensions, while v-> contains
            // only numerical values
            struct adios_var_struct * old_var = adios_find_var_by_id (g->vars, v->parent_id);
            v_index->characteristics [0].payload_offset = v->write_offset
                            + adios_calc_var_overhead_v1 (old_var)
                            - strlen (old_var->path)  // take out the length of path defined in XML
                            + strlen (v->path); // add length of the actual, current path of this var
            v_index->characteristics [0].file_index = fd->subfile_index;
            v_index->characteristics [0].time_index = g_item->time_index;

            v_index->characteristics [0].value = 0;
            v_index->characteristics [0].dims.count = 0;

            // NCSU - Initializing stat related info in index
            v_index->characteristics [0].bitmap = 0;
            v_index->characteristics [0].stats = 0;

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

                        // NCSU - Copy statistics from var struct to index
                        uint8_t count = adios_get_stat_set_count(v->type);
                        uint8_t idx = 0;
                        uint64_t characteristic_size;

                        v_index->characteristics [0].bitmap = v->bitmap;
                        v_index->characteristics [0].stats = malloc (count * sizeof(struct adios_index_characteristics_stat_struct *));

                        // Set of characteristics will be repeated thrice for complex numbers
                        for (c = 0; c < count; c ++)
                        {
                            v_index->characteristics [0].stats[c] = calloc(ADIOS_STAT_LENGTH, sizeof (struct adios_index_characteristics_stat_struct));

                            j = idx = 0;
                            while (v_index->characteristics [0].bitmap >> j)
                            {
                                if ((v_index->characteristics [0].bitmap >> j) & 1)
                                {
                                    if (v->stats[c][idx].data != NULL)
                                    {
                                        if (j == adios_statistic_hist)
                                        {
                                            v_index->characteristics [0].stats[c][idx].data = (struct adios_index_characteristics_hist_struct *) malloc (sizeof(struct adios_index_characteristics_hist_struct));

                                            struct adios_hist_struct * v_hist = v->stats[c][idx].data;
                                            struct adios_hist_struct * v_index_hist = v_index->characteristics [0].stats[c][idx].data;

                                            v_index_hist->min = v_hist->min;
                                            v_index_hist->max = v_hist->max;
                                            v_index_hist->num_breaks = v_hist->num_breaks;

                                            v_index_hist->frequencies = malloc ((v_hist->num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
                                            memcpy (v_index_hist->frequencies, v_hist->frequencies, (v_hist->num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
                                            v_index_hist->breaks = malloc ((v_hist->num_breaks) * adios_get_type_size(adios_double, ""));
                                            memcpy (v_index_hist->breaks, v_hist->breaks, (v_hist->num_breaks) * adios_get_type_size(adios_double, ""));
                                        }
                                        else
                                        {
                                            characteristic_size = adios_get_stat_size(v->stats[c][idx].data, v->type, j);
                                            v_index->characteristics [0].stats[c][idx].data = malloc (characteristic_size);
                                            memcpy (v_index->characteristics [0].stats[c][idx].data, v->stats[c][idx].data, characteristic_size);
                                        }

                                        idx ++;
                                    }
                                }
                                j ++;
                            }
                        }
                        // NCSU - End of copy, for statistics

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
                        v_index->characteristics [0].bitmap = 0;
                        v_index->characteristics [0].stats = 0;

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
            index_append_var_v1 (vars_root, v_index);
        }

        v = v->next;
    }

    while (a)
    {
        // only add items that were written to the index
        if (a->write_offset != 0)
        {
            struct adios_index_attribute_struct_v1 * a_index;
            a_index = malloc (sizeof (struct adios_index_attribute_struct_v1));
            a_index->characteristics = malloc (
                           sizeof (struct adios_index_characteristic_struct_v1)
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
            a_index->characteristics [0].payload_offset = a->write_offset + adios_calc_attribute_overhead_v1 (a);
            a_index->characteristics [0].file_index = fd->subfile_index;
            a_index->characteristics [0].time_index = 0;

            // NCSU -,Initializing stat related info in attribute index
            a_index->characteristics [0].bitmap = 0;
            a_index->characteristics [0].stats = 0;

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
            index_append_attribute_v1 (attrs_root, a_index);
        }

        a = a->next;
    }
}

int adios_write_index_v1 (char ** buffer
                         ,uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,uint64_t index_start
                         ,struct adios_index_process_group_struct_v1 * pg_root
                         ,struct adios_index_var_struct_v1 * vars_root
                         ,struct adios_index_attribute_struct_v1 * attrs_root
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

            // add a file index characteristic for all vars
            characteristic_set_count++;
            flag = (uint8_t) adios_characteristic_file_index;
            buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
            index_size += 1;
            var_size += 1;
            characteristic_set_length += 1;

            buffer_write (buffer, buffer_size, buffer_offset
                         ,&vars_root->characteristics [i].file_index, 4
                         );
            index_size += 4;
            var_size += 4;
            characteristic_set_length += 4;

            // add a time index characteristic for all vars
            characteristic_set_count++;
            flag = (uint8_t) adios_characteristic_time_index;
            buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
            index_size += 1;
            var_size += 1;
            characteristic_set_length += 1;

            buffer_write (buffer, buffer_size, buffer_offset
                         ,&vars_root->characteristics [i].time_index, 4
                         );
            index_size += 4;
            var_size += 4;
            characteristic_set_length += 4;

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

                        // NCSU - Adding bitmap
                        characteristic_set_count++;
                        flag = (uint8_t) adios_characteristic_bitmap;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&flag, 1
                                     );
                        index_size += 1;
                        var_size += 1;
                        characteristic_set_length += 1;

                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&vars_root->characteristics [i].bitmap, 4
                                     );
                        index_size += 4;
                        var_size += 4;
                        characteristic_set_length += 4;

                        // NCSU - Adding statistics
                           characteristic_set_count++;
                        flag = (uint8_t) adios_characteristic_stat;
                        buffer_write (buffer, buffer_size, buffer_offset
                                     ,&flag, 1
                                     );
                        index_size += 1;
                        var_size += 1;
                        characteristic_set_length += 1;

                        uint8_t count = adios_get_stat_set_count(vars_root->type);
                        uint8_t idx = 0, c, j;
                        uint64_t characteristic_size;

                        for (c = 0; c < count; c ++)
                        {
                            j = idx = 0;

                            while (vars_root->characteristics [i].bitmap >> j)
                            {
                                if ((vars_root->characteristics [i].bitmap >> j) & 1)
                                {
                                    if (j == adios_statistic_hist)
                                    {
                                        struct adios_hist_struct * hist = vars_root->characteristics [i].stats[c][idx].data;

                                        buffer_write (  buffer, buffer_size, buffer_offset
                                                        , &hist->num_breaks, 4
                                                     );
                                        characteristic_size = 4;

                                        buffer_write (  buffer, buffer_size, buffer_offset
                                                        , &hist->min, 8
                                                     );
                                        characteristic_size += 8;

                                        buffer_write (  buffer, buffer_size, buffer_offset
                                                        , &hist->max, 8
                                                     );
                                        characteristic_size += 8;

                                        buffer_write (  buffer, buffer_size, buffer_offset
                                                        , hist->frequencies, (hist->num_breaks + 1) * 4
                                                     );
                                        characteristic_size += (hist->num_breaks + 1) * 4;

                                        buffer_write (  buffer, buffer_size, buffer_offset
                                                        , hist->breaks, hist->num_breaks * 8
                                                     );
                                        characteristic_size += (hist->num_breaks) * 8;
                                    }
                                    else
                                    {
                                        characteristic_size = adios_get_stat_size(vars_root->characteristics [i].stats[c][idx].data, vars_root->type, j);

                                        buffer_write ( 	buffer, buffer_size, buffer_offset
                                                         ,vars_root->characteristics [i].stats[c][idx].data, characteristic_size
                                                      );

                                    }

                                    index_size += characteristic_size;
                                    var_size += characteristic_size;
                                    characteristic_set_length += characteristic_size;

                                    idx ++;
                                }
                                j ++;
                            }
                        }
                        // NCSU - End of addition statistic to buffer
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

            // add a file index characteristic for all attrs
            characteristic_set_count++;
            flag = (uint8_t) adios_characteristic_file_index;
            buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
            index_size += 1;
            attr_size += 1;
            characteristic_set_length += 1;

            buffer_write (buffer, buffer_size, buffer_offset
                         ,&attrs_root->characteristics [i].file_index, 4
                         );
            index_size += 4;
            attr_size += 4;
            characteristic_set_length += 4;

            // add a time index characteristic for all attrs
            characteristic_set_count++;
            flag = (uint8_t) adios_characteristic_time_index;
            buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
            index_size += 1;
            attr_size += 1;
            characteristic_set_length += 1;

            buffer_write (buffer, buffer_size, buffer_offset
                         ,&attrs_root->characteristics [i].time_index, 4
                         );
            index_size += 4;
            attr_size += 4;
            characteristic_set_length += 4;

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

int adios_write_version_v1 (char ** buffer
                           ,uint64_t * buffer_size
                           ,uint64_t * buffer_offset
                           )
{
    uint32_t test = 1;

    if (!*(char *) &test)
        test = 0x80000000;
    else
        test = 0;

    test += 1;   // first data storage version
    // For the new read API to be able to read back older file format,
    // set this flag
    test |= ADIOS_VERSION_HAVE_TIME_INDEX_CHARACTERISTIC;

    test = htonl (test);

    buffer_write (buffer, buffer_size, buffer_offset, &test, 4);

    return 0;
}

int adios_write_version_flag_v1 (char ** buffer
                                ,uint64_t * buffer_size
                                ,uint64_t * buffer_offset
                                ,uint32_t flag
                                )
{
    uint32_t test = 1;

    if (!*(char *) &test)
        test = 0x80000000;
    else
        test = 0;

    // version number 1 byte, endiness 1 byte,
    // the rest is user-defined options 2 bytes
    test += 1;   // master index file version
    // For the new read API to be able to read back older file format,
    // set this flag
    test |= ADIOS_VERSION_HAVE_TIME_INDEX_CHARACTERISTIC | flag;

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
uint64_t adios_write_dimension_v1 (struct adios_file_struct * fd
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
uint16_t adios_write_dimensions_v1 (struct adios_file_struct * fd
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
        size += adios_write_dimension_v1 (fd, dimensions);

        dimensions = dimensions->next;
    }

    return size;
}

uint16_t adios_write_var_characteristics_dims_v1 (struct adios_file_struct * fd
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

uint16_t adios_write_var_characteristics_v1 (struct adios_file_struct * fd
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

                len = adios_write_var_characteristics_dims_v1 (fd, v);
                index_size += len;
                characteristic_set_length += len;

                // add the bitmap of characteristics
                characteristic_set_count++;
                flag = (uint8_t) adios_characteristic_bitmap;
                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                             ,&flag, 1
                             );
                index_size += 1;
                characteristic_set_length += 1;

                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                             ,&v->bitmap, 4
                             );
                index_size += 4;
                characteristic_set_length += 4;

                // add a stat value characteristic
                characteristic_set_count++;
                flag = (uint8_t) adios_characteristic_stat;
                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                             ,&flag, 1
                             );
                index_size += 1;
                characteristic_set_length += 1;

                uint8_t j, c;
                uint8_t count = adios_get_stat_set_count(v->type);
                uint8_t idx = 0;
                uint64_t characteristic_size;

                for (c = 0; c < count; c ++)
                {
                    j = idx = 0;
                    while (v->bitmap >> j)
                    {
                        if ((v->bitmap >> j) & 1)
                        {
                            if (j == adios_statistic_hist)
                            {
                                struct adios_hist_struct * hist = v->stats[c][idx].data;
                                int32_t num_breaks = hist->num_breaks;
                                // Adding number of bins
                                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                                             ,&hist->num_breaks, 4
                                             );
                                characteristic_size = 4;

                                 // Adding min bin
                                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                                            ,&hist->min, 8
                                             );
                                characteristic_size += 8;

                                 // Adding max bin
                                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                                             ,&hist->max, 8
                                             );
                                characteristic_size += 8;

                                // add a frequencies value characteristic
                                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                                             ,hist->frequencies, 4 * (num_breaks + 1)
                                             );
                                characteristic_size += 4 * (num_breaks + 1);

                                // add the breaks value characteristic
                                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                                             ,hist->breaks, 8 * num_breaks
                                             );
                                characteristic_size += 8 * num_breaks;
                            }
                            else
                            {
                                characteristic_size = adios_get_stat_size(v->stats[c][idx].data, v->type, j);

                                buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset
                                         ,v->stats[c][idx].data, characteristic_size
                                         );
                            }

                            index_size += characteristic_size;
                            characteristic_set_length += characteristic_size;
                            idx ++;
                        }
                        j ++;
                    }
                }
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

int adios_generate_var_characteristics_v1 (struct adios_file_struct * fd
                                    ,struct adios_var_struct * var
                                    )
{
    uint64_t total_size = adios_get_var_size (var, fd->group, var->data);
    uint64_t size = 0;

    if (var->bitmap == 0)
        return ;

    int32_t map[32];
    memset (map, -1, sizeof(map));

#if 1
#define HIST(a) \
{ \
int j = 0, low, high, mid; \
low=0; \
high=hist->num_breaks - 1; \
if (hist->breaks[low] > a) \
    hist->frequencies[0] += 1; \
else if (a >= hist->breaks[high]) \
    hist->frequencies[high + 1] += 1; \
else if (hist->breaks[low] <= a && a < hist->breaks[high]) \
{ \
    while(high-low>=2) \
    { \
        mid=(high+low)/2; \
        if(a >= hist->breaks[mid]) \
        { \
            low=mid; \
        } \
        else \
        { \
            high=mid; \
        } \
    } \
    hist->frequencies[low + 1] += 1; \
} \
}
#endif

#if 1
#define ADIOS_STATISTICS(a,b) \
{\
a * data = (a *) var->data; \
int i, j; \
struct adios_stat_struct * stats = var->stats[0]; \
a * min, * max; \
double * sum, * sum_square; \
uint32_t * cnt; \
struct adios_hist_struct * hist = 0; \
i = j = 0; \
while (var->bitmap >> j) { \
    if ((var->bitmap >> j) & 1)	{\
        map [j] = i; \
        if (j == adios_statistic_hist) \
            ;\
        else \
            stats[i].data = malloc(adios_get_stat_size(NULL, var->type, j)); \
        i ++; \
    } \
    j ++; \
} \
min = (a *) stats[map[adios_statistic_min]].data; \
max = (a *) stats[map[adios_statistic_max]].data; \
sum = (double *) stats[map[adios_statistic_sum]].data; \
sum_square = (double *) stats[map[adios_statistic_sum_square]].data; \
cnt = (uint32_t *) stats[map[adios_statistic_cnt]].data; \
*cnt = 0;\
if (map[adios_statistic_hist] != -1) {\
    hist = (struct adios_hist_struct *) stats[map[adios_statistic_hist]].data; \
    hist->frequencies = calloc ((hist->num_breaks + 1), adios_get_type_size(adios_unsigned_integer, "")); \
} \
int finite = 0; \
size = 0; \
while ((size * b) < total_size) \
{ \
    if (isnan (data [size]) || !isfinite (data [size])) {\
        size ++; \
        continue; \
    }\
    if (!finite) { \
        *min = data [size]; \
        *max = data [size]; \
        *sum = data [size]; \
        *sum_square = (data [size] * data [size]) ; \
        *cnt = *cnt + 1; \
        if (map[adios_statistic_hist] != -1) \
            HIST(data [size]); \
        finite = 1; \
        size ++; \
        continue; \
    } \
    if (data [size] < *min) \
        *min = data [size]; \
    if (data [size] > *max) \
        *max = data [size]; \
    *sum += data [size]; \
    *sum_square += (data [size] * data [size]) ; \
    *cnt = *cnt + 1; \
    if (map[adios_statistic_hist] != -1) \
        HIST(data [size]); \
       size++; \
} \
if (map[adios_statistic_finite] != -1) \
    * ((uint8_t * ) stats[map[adios_statistic_finite]].data) = finite; \
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
            ADIOS_STATISTICS(int8_t,1)

        case adios_unsigned_byte:
            ADIOS_STATISTICS(uint8_t,1)

        case adios_short:
            ADIOS_STATISTICS(int16_t,2)

        case adios_unsigned_short:
            ADIOS_STATISTICS(uint16_t,2)

        case adios_integer:
            ADIOS_STATISTICS(int32_t,4)

        case adios_unsigned_integer:
            ADIOS_STATISTICS(uint32_t,4)

        case adios_long:
            ADIOS_STATISTICS(int64_t,8)

        case adios_unsigned_long:
            ADIOS_STATISTICS(uint64_t,8)

        case adios_real:
            ADIOS_STATISTICS(float,4)

        case adios_double:
            ADIOS_STATISTICS(double,8)

        case adios_long_double:
            ADIOS_STATISTICS(long double,16)

        case adios_complex:
        {
            int i, j, c, count = 3;
            struct adios_stat_struct ** stats = var->stats;
            float * data = var->data;
            i = j = 0;

            while (var->bitmap >> j) {
                if ((var->bitmap >> j) & 1)	{
                    map [j] = i;
                    for (c = 0; c < count; c ++)
                        if (j != adios_statistic_hist)
                            stats[c][i].data = malloc(adios_get_stat_size(NULL, var->type, j));
                    i ++;
                }
                j ++;
            }

            double *min_i, *min_r, *min_m;
            double *max_i, *max_r, *max_m;
            double *sum_i, *sum_r, *sum_m;
            double *sum_square_i, *sum_square_r, *sum_square_m;
            //struct adios_hist_struct *hist, *hist_i, *hist_r, *hist_m;
            uint32_t *cnt_i, *cnt_r, *cnt_m;
            double magnitude;
            uint8_t finite = 0;

            min_m = (double *) stats[0][map[adios_statistic_min]].data;
            min_r = (double *) stats[1][map[adios_statistic_min]].data;
            min_i = (double *) stats[2][map[adios_statistic_min]].data;

            max_m = (double *) stats[0][map[adios_statistic_max]].data;
            max_r = (double *) stats[1][map[adios_statistic_max]].data;
            max_i = (double *) stats[2][map[adios_statistic_max]].data;

            sum_m = (double *) stats[0][map[adios_statistic_sum]].data;
            sum_r = (double *) stats[1][map[adios_statistic_sum]].data;
            sum_i = (double *) stats[2][map[adios_statistic_sum]].data;

            sum_square_m = (double *) stats[0][map[adios_statistic_sum_square]].data;
            sum_square_r = (double *) stats[1][map[adios_statistic_sum_square]].data;
            sum_square_i = (double *) stats[2][map[adios_statistic_sum_square]].data;

            cnt_m = (uint32_t *) stats[0][map[adios_statistic_cnt]].data;
            cnt_r = (uint32_t *) stats[1][map[adios_statistic_cnt]].data;
            cnt_i = (uint32_t *) stats[2][map[adios_statistic_cnt]].data;

            // Histogram is not available for complex numbers, yet.
            /*
            if (map[adios_statistic_hist] != -1) {
                hist_r = (struct adios_hist_struct *) stat[0][map[adios_statistic_hist]].data;
                hist_i = (struct adios_hist_struct *) stat[1][map[adios_statistic_hist]].data;
                hist_m = (struct adios_hist_struct *) stat[2][map[adios_statistic_hist]].data;

                hist_r->frequencies = calloc ((hist->num_breaks + 1), adios_get_type_size(adios_unsigned_integer, ""));
                hist_i->frequencies = calloc ((hist->num_breaks + 1), adios_get_type_size(adios_unsigned_integer, ""));
                hist_m->frequencies = calloc ((hist->num_breaks + 1), adios_get_type_size(adios_unsigned_integer, ""));
            }
            */

            *cnt_r = *cnt_i = *cnt_m = 0;
            *min_r = *min_i = *min_m = INFINITY;
            *max_r = *max_i = *max_m = -INFINITY;
            *sum_r = *sum_i = *sum_m = 0;
            *sum_square_r = *sum_square_i = *sum_square_m = 0;

            while ((size * sizeof(float)) < total_size) {

                magnitude = sqrt((double) data [size] * data [size] + (double) data[size + 1] * data[size + 1]);

                // Both real and imaginary parts have to be finite, else skip calculating the characteristic
                if ( isnan(data [size]) || !isfinite(data [size]) || isnan(data[size + 1]) || !isfinite(data[size + 1]) ) {
                    size += 2;
                    continue;
                }

                finite = 1;

                // Updating the characteristic values
                if (data [size] < *min_r)
                    *min_r = data [size];
                if (data [size + 1] < *min_i)
                    *min_i = data [size + 1];
                if (magnitude < *min_m)
                    *min_m = magnitude;

                if (data [size] > *max_r)
                    *max_r = data [size];
                if (data [size + 1] > *max_i)
                    *max_i = data [size + 1];
                if (magnitude > *max_m)
                    *max_m = magnitude;

                *sum_r += data [size];
                *sum_i += data [size + 1];
                *sum_m += magnitude;

                *sum_square_r += (double) data [size] * data [size];
                *sum_square_i += (double) data [size + 1] * data [size + 1];
                *sum_square_m += magnitude * magnitude;

                *cnt_r = *cnt_r + 1;
                *cnt_i = *cnt_i + 1;
                *cnt_m = *cnt_m + 1;
                // Histogram not available yet
                /*
                if (map[adios_statistic_hist] != -1)
                {
                    hist = hist_r;
                    HIST (data[size]);

                    hist = hist_i;
                    HIST (data[size + 1]);

                    hist = hist_m;
                    HIST (magnitude);
                }
                */

                   size += 2;
            }

            if (map[adios_statistic_finite] != -1)
                for (c = 0; c < count; c ++)
                    * ((uint8_t * ) stats[c][map[adios_statistic_finite]].data) = finite;

            return 0;
        }

        case adios_double_complex:
        {
            int i, j, c, count = 3;
            struct adios_stat_struct ** stats = var->stats;
            double * data = var->data;
            i = j = 0;

            while (var->bitmap >> j) {
                if ((var->bitmap >> j) & 1)	{
                    map [j] = i;
                    for (c = 0; c < count; c ++)
                        if (j != adios_statistic_hist)
                            stats[c][i].data = malloc(adios_get_stat_size(NULL, var->type, j));
                    i ++;
                }
                j ++;
            }

            long double *min_i, *min_r, *min_m;
            long double *max_i, *max_r, *max_m;
            long double *sum_i, *sum_r, *sum_m;
            long double *sum_square_i, *sum_square_r, *sum_square_m;
            //struct adios_hist_struct *hist, *hist_i, *hist_r, *hist_m;
            uint32_t *cnt_i, *cnt_r, *cnt_m;
            long double magnitude;
            uint8_t finite = 0;

            min_m = (long double *) stats[0][map[adios_statistic_min]].data;
            min_r = (long double *) stats[1][map[adios_statistic_min]].data;
            min_i = (long double *) stats[2][map[adios_statistic_min]].data;

            max_m = (long double *) stats[0][map[adios_statistic_max]].data;
            max_r = (long double *) stats[1][map[adios_statistic_max]].data;
            max_i = (long double *) stats[2][map[adios_statistic_max]].data;

            sum_m = (long double *) stats[0][map[adios_statistic_sum]].data;
            sum_r = (long double *) stats[1][map[adios_statistic_sum]].data;
            sum_i = (long double *) stats[2][map[adios_statistic_sum]].data;

            sum_square_m = (long double *) stats[0][map[adios_statistic_sum_square]].data;
            sum_square_r = (long double *) stats[1][map[adios_statistic_sum_square]].data;
            sum_square_i = (long double *) stats[2][map[adios_statistic_sum_square]].data;


            // Histogram not available for complex numbers yet
            /*
            if (map[adios_statistic_hist] != -1) {
                hist_r = (struct adios_hist_struct *) stat[0][map[adios_statistic_hist]].data;
                hist_i = (struct adios_hist_struct *) stat[1][map[adios_statistic_hist]].data;
                hist_m = (struct adios_hist_struct *) stat[2][map[adios_statistic_hist]].data;

                hist_r->frequencies = calloc ((hist->num_breaks + 1), adios_get_type_size(adios_unsigned_integer, ""));
                hist_i->frequencies = calloc ((hist->num_breaks + 1), adios_get_type_size(adios_unsigned_integer, ""));
                hist_m->frequencies = calloc ((hist->num_breaks + 1), adios_get_type_size(adios_unsigned_integer, ""));
            }
            */

            cnt_m = (uint32_t *) stats[0][map[adios_statistic_cnt]].data;
            cnt_r = (uint32_t *) stats[1][map[adios_statistic_cnt]].data;
            cnt_i = (uint32_t *) stats[2][map[adios_statistic_cnt]].data;

            *cnt_r = *cnt_i = *cnt_m = 0;
            *min_r = *min_i = *min_m = INFINITY;
            *max_r = *max_i = *max_m = -INFINITY;
            *sum_r = *sum_i = *sum_m = 0;
            *sum_square_r = *sum_square_i = *sum_square_m = 0;

            while ((size * sizeof(double)) < total_size)
            {
                long double magnitude = sqrt((long double) data [size] * data [size] + (long double) data[size + 1] * data[size + 1]);

                // Both real and imaginary parts have to be finite, else skip calculating the characteristic
                if ( isnan(data [size]) || !isfinite(data [size]) || isnan(data[size + 1]) || !isfinite(data[size + 1]) ) {
                    size += 2;
                    continue;
                }

                finite = 1;

                if (data [size] < *min_r)
                    *min_r = data [size];
                if (data [size + 1] < *min_i)
                    *min_i = data [size + 1];
                if (magnitude < *min_m)
                    *min_m = magnitude;

                if (data [size] > *max_r)
                    *max_r = data [size];
                if (data [size + 1] > *max_i)
                    *max_i = data [size + 1];
                if (magnitude > *max_m)
                    *max_m = magnitude;

                *sum_r += data [size];
                *sum_i += data [size + 1];
                *sum_m += magnitude;

                *sum_square_r += (long double) data [size] * data [size];
                *sum_square_i += (long double) data [size + 1] * data [size + 1];
                *sum_square_m += magnitude * magnitude;

                // Histgram has not available for complex yet
                /*
                if (map[adios_statistic_hist] != -1)
                {
                    hist = hist_r;
                    HIST (data[size]);

                    hist = hist_i;
                    HIST (data[size + 1]);

                    hist = hist_m;
                    HIST (magnitude);
                }
                */

                *cnt_r = *cnt_r + 1;
                *cnt_i = *cnt_i + 1;
                *cnt_m = *cnt_m + 1;

                   size += 2;
            }

            if (map[adios_statistic_finite] != -1)
                for (c = 0; c < count; c ++)
                    * ((uint8_t * ) stats[c][map[adios_statistic_finite]].data) = finite;

            return 0;
        }

        case adios_string:
        {
            var->stats = 0;

            return 0;
        }

        default:
        {
            uint64_t data = adios_unknown;
            var->stats = 0;

            return 0;
        }
    }
}

// data is only there for sizing
uint64_t adios_write_var_header_v1 (struct adios_file_struct * fd
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

    total_size += adios_write_dimensions_v1 (fd, v->dimensions);

    adios_generate_var_characteristics_v1 (fd, v);
    total_size += adios_write_var_characteristics_v1 (fd, v);

    total_size += adios_get_var_size (v, fd->group, v->data); // payload

    buffer_write (&fd->buffer, &fd->buffer_size, &start, &total_size, 8);

    fd->vars_written++;

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return total_size;
}

int adios_write_var_payload_v1 (struct adios_file_struct * fd
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

int adios_write_attribute_v1 (struct adios_file_struct * fd
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

int adios_write_open_vars_v1 (struct adios_file_struct * fd)
{
    fd->vars_written = 0;

    // it is now setup to write the vars and then the attrs on close
    fd->vars_start = fd->offset;

    fd->offset += (2 + 8); // (count + size)

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return 0;
}

int adios_write_close_vars_v1 (struct adios_file_struct * fd)
{
    // close the var area (count and total size) and write the attributes
    uint64_t size = fd->offset - fd->vars_start;
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->vars_start, &fd->vars_written, 2);

    buffer_write (&fd->buffer, &fd->buffer_size, &fd->vars_start, &size, 8);

    return 0;
}

int adios_write_open_attributes_v1 (struct adios_file_struct * fd)
{
    fd->vars_start = fd->offset;   // save the start of attr area for size
    fd->offset += (2 + 8);         // space to write the count and size
    fd->vars_written = 0;

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return 0;
}

int adios_write_close_attributes_v1 (struct adios_file_struct * fd)
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

uint64_t adios_get_var_size (struct adios_var_struct * var
                            ,struct adios_group_struct * group, void * data
                            )
{
    uint64_t size = 0;

    size = adios_get_type_size (var->type, data);

    if (var->dimensions)
    {
        struct adios_dimension_struct * d = var->dimensions;

        while (d)
        {
            // calculate the size for this dimension element
            if (d->dimension.id != 0)
            {
                struct adios_var_struct * dim_var = 0;

                dim_var = adios_find_var_by_id (group->vars, d->dimension.id);

                // first check to make sure all vars are provided
                if (!dim_var)
                {
                    struct adios_attribute_struct * attr = 0;
                    attr = adios_find_attribute_by_id (group->attributes
                                                      ,d->dimension.id
                                                      );
                    if (attr)
                    {
                        if (attr->var)
                        {
                            if (!attr->var->data)
                            {
                                fprintf (stderr, "adios_get_var_size: "
                                                 "sizing of %s failed because "
                                                 "dimension component %s was "
                                                 "not provided\n"
                                        ,var->name, attr->var->name
                                        );

                                return 0;
                            }
                            else
                            {
                                if (!adios_multiply_dimensions (&size, var
                                                               ,attr->var->type
                                                               ,attr->var->data
                                                               )
                                   )
                                {
                                    return 0;
                                }
                            }
                        }
                        else
                        {
                            if (!adios_multiply_dimensions (&size, var
                                                           ,attr->type
                                                           ,attr->value
                                                           )
                               )
                            {
                                return 0;
                            }
                        }
                    }
                    else
                    {
                        fprintf (stderr, "adios_get_var_size: "
                                         "sizing of %s failed because "
                                         "dimension component was not "
                                         "provided\n"
                                ,var->name
                                );

                        return 0;
                    }
                }
                else
                {
                    if (!dim_var->data)
                    {
                        fprintf (stderr, "adios_get_var_size: "
                                         "sizing of %s failed because "
                                         "dimension component %s was not "
                                         "provided\n"
                                ,var->name, dim_var->name
                                );

                        return 0;
                    }
                    else
                    {
                        if (!adios_multiply_dimensions (&size, var
                                                       ,dim_var->type
                                                       ,dim_var->data
                                                       )
                           )
                        {
                            return 0;
                        }
                    }
                }
            }
            else
            {
                if (d->dimension.time_index == adios_flag_no)
                {
                    size *= d->dimension.rank;
                }
                // the time index doesn't take up space...
            }

            d = d->next;
        }
    }

    return size;
}

const char * adios_type_to_string_int (int type)
{
    switch (type)
    {
        case adios_unsigned_byte:    return "unsigned byte";
        case adios_unsigned_short:   return "unsigned short";
        case adios_unsigned_integer: return "unsigned integer";
        case adios_unsigned_long:    return "unsigned long long";

        case adios_byte:             return "byte";
        case adios_short:            return "short";
        case adios_integer:          return "integer";
        case adios_long:             return "long long";

        case adios_real:             return "real";
        case adios_double:           return "double";
        case adios_long_double:      return "long double";

        case adios_string:           return "string";
        case adios_complex:          return "complex";
        case adios_double_complex:   return "double complex";

        default:
        {
            static char buf [50];
            sprintf (buf, "(unknown: %d)", type);
            return buf;
        }
    }
}

const char * adios_file_mode_to_string (int mode)
{
    static char buf [50];

    switch (mode)
    {
        case adios_mode_write:  return "write";
        case adios_mode_read:   return "read";
        case adios_mode_update: return "update";
        case adios_mode_append: return "append";

        default:
            sprintf (buf, "(unknown: %d)", mode);
    }

    return buf;
}

//////////////////////////////////////////////////////////////////////////////
// Queue management code intended for adaptive API use
//////////////////////////////////////////////////////////////////////////////
void list_init (List * list, void (* destroy) (void * data))
{
    list->size = 0;
    list->destroy = destroy;
    list->head = NULL;
    list->tail = NULL;

    return;
}

void list_destroy (List * list)
{
    void * data;

    while (list_size (list) > 0)
    {
        if (list_rem_next (list, NULL, (void **) &data) == 0 && list->destroy != NULL)
        {
            list->destroy (data);
        }
    }

    memset (list, 0, sizeof (List));

    return;
}

int list_ins_next (List * list, ListElmt * element, const void * data)
{
    ListElmt * new_element;

    if ((new_element = (ListElmt *) malloc (sizeof (ListElmt))) == NULL)
        return -1;

    new_element->data = (void *) data;

    if (element == NULL)
    {
        if  (list_size (list) == 0)
            list->tail = new_element;

        new_element->next = list->head;
        list->head = new_element;
    }
    else
    {
        if (element->next == NULL)
            list->tail = new_element;

        new_element->next = element->next;
        element->next = new_element;
    }

    list->size++;

    return 0;
}

int list_rem_next (List * list, ListElmt * element, void ** data)
{
    ListElmt * old_element;

    if (list_size (list) == 0)
        return -1;

    if (element == NULL)
    {
        *data = list->head->data;
        old_element = list->head;
        list->head = list->head->next;

        if (list_size (list) == 1)
            list->tail = NULL;
    }
    else
    {
        if (element->next == NULL)
            return -1;

        *data = element->next->data;
        old_element = element->next;
        element->next = element->next->next;

        if (element->next == NULL)
            list->tail = element;
    }

    free (old_element);

    list->size--;

    return 0;
}

int queue_enqueue (Queue * queue, const void * data)
{
    return list_ins_next (queue, list_tail (queue), data);
}

int queue_dequeue (Queue * queue, void ** data)
{
    return list_rem_next (queue, NULL, data);
}
