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

static struct adios_method_list_struct * adios_methods = 0;
static struct adios_group_list_struct * adios_groups = 0;

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

        if (g->vars->min)
            free (g->vars->min);
        if (g->vars->max)
            free (g->vars->max);
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

static void tokenize_dimensions (char * str, char *** tokens, int * count)
{
    if (!str)
    {
        *tokens = 0;
        *count = 0;

        return;
    }

    trim_spaces (str);

    char * t = str;
    char * save_str = strdup (str);
    int i;

    if (strlen (str) > 0)
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
    v->min = 0;
    v->max = 0;

    v->data_size = 0;

    v->next = 0;

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
                overhead += 1;  // id
                overhead += adios_get_type_size (v->type, ""); // min

                overhead += 1;  // id
                overhead += adios_get_type_size (v->type, ""); // max

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
            var_new->id = ++fd->group->member_count;
            var_new->parent_id = var->id;
            var_new->name = strdup (var->name);
            var_new->path = strdup (var->path);
            var_new->type = var->type;
            var_new->dimensions = 0;
            var_new->got_buffer = var->got_buffer;
            var_new->is_dim = var->is_dim;
            var_new->write_offset = var->write_offset;
            var_new->min = 0;
            var_new->max = 0;
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
                        var_new->min = malloc (size);
                        var_new->max = malloc (size);
                        memcpy (var_new->min, var->min, size);
                        memcpy (var_new->max, var->max, size);
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
                        var_new->min = 0;
                        var_new->max = 0;
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
            // We need this to calculate the correct payload_offset
            struct adios_var_struct * old_var = adios_find_var_by_id (g->vars, v->parent_id);
            v_index->characteristics [0].payload_offset = v->write_offset
                                                    + adios_calc_var_overhead_v1 (old_var);
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

    test = htonl (test);

    buffer_write (buffer, buffer_size, buffer_offset, &test, 4);

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

int adios_generate_var_characteristics_v1 (struct adios_file_struct * fd
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
