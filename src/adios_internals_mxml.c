/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/stat.h> /* struct stat */

// xml parser
#include <mxml.h>

#ifdef _NOMPI
    /* Sequential processes can use the library compiled with -D_NOMPI */
#   include "mpidummy.h"
#else
    /* Parallel applications should use MPI to communicate  */
#   include "mpi.h"
#endif

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "adios_internals_mxml.h"
#include "buffer.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

static enum ADIOS_FLAG adios_host_language_fortran = adios_flag_yes;

struct adios_method_list_struct * adios_methods = 0;
struct adios_group_list_struct * adios_groups = 0;

//extern struct adios_method_list_struct * adios_methods;
//extern struct adios_group_list_struct * adios_groups;

struct adios_transport_struct * adios_transports = 0;
static int adios_transports_initialized = 0;

// this macro makes getting the attributes easier
// fix the bgp bugs
#define GET_ATTR(n,attr,var,en)                              \
if (!strcasecmp (n, attr->name)) {                           \
    if (!var)                                                \
    {                                                        \
        var = attr->value;                                   \
        continue;                                            \
    }                                                        \
    else                                                     \
    {                                                        \
        fprintf (stderr, "xml: duplicate attribute %s on %s (ignored)",n,en); \
        continue;                                            \
    }                                                        \
}

static enum ADIOS_DATATYPES parseType (const char * type, const char * name)
{
    if (   !strcasecmp (type, "byte")
        || !strcasecmp (type, "integer*1")
       )
        return adios_byte;

    if (   !strcasecmp (type, "short")
        || !strcasecmp (type, "integer*2")
       )
        return adios_short;

    if (   !strcasecmp (type, "integer")
        || !strcasecmp (type, "integer*4")
       )
        return adios_integer;

    if (   !strcasecmp (type, "long")
        || !strcasecmp (type, "integer*8")
       )
        return adios_long;

    if (   !strcasecmp (type, "unsigned byte")
        || !strcasecmp (type, "unsigned integer*1")
       )
        return adios_unsigned_byte;

    if (   !strcasecmp (type, "unsigned short")
        || !strcasecmp (type, "unsigned integer*2")
       )
        return adios_unsigned_short;

    if (   !strcasecmp (type, "unsigned integer")
        || !strcasecmp (type, "unsigned integer*4")
       )
        return adios_unsigned_integer;

    if (   !strcasecmp (type, "unsigned long")
        || !strcasecmp (type, "unsigned integer*8")
       )
        return adios_unsigned_long;

    if (   !strcasecmp (type, "real")
        || !strcasecmp (type, "real*4")
        || !strcasecmp (type, "float")
       )
        return adios_real;

    if (   !strcasecmp (type, "real*8")
        || !strcasecmp (type, "double")
        || !strcasecmp (type, "long float")
       )
        return adios_double;

    if (   !strcasecmp (type, "real*16")
        || !strcasecmp (type, "long double")
       )
        return adios_long_double;

    if (!strcasecmp (type, "string"))
        return adios_string;

    if (   !strcasecmp (type, "complex")
        || !strcasecmp (type, "complex*8")
       )
        return adios_complex;

    if (   !strcasecmp (type, "double complex")
        || !strcasecmp (type, "complex*16")
       )
        return adios_double_complex;

    fprintf (stderr, "config.xml: invalid type: %s in var %s\n", type, name);

    return adios_unknown;
}

static enum ADIOS_FLAG parseFlag (const char * attr_name, const char * flag
                                 ,enum ADIOS_FLAG default_value
                                 )
{
    if (!flag)
        return default_value;

    if (!strcasecmp (flag, "yes"))
        return adios_flag_yes;

    if (!strcasecmp (flag, "no"))
        return adios_flag_no;

    fprintf (stderr, "config.xml: %s must have a value of 'yes' or 'no' "
                     "not: %s\n", attr_name, flag
            );

    return adios_flag_unknown;
}


static void adios_append_mesh_item (struct adios_mesh_item_list_struct ** root
                                   ,struct adios_mesh_item_list_struct * item
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

static void adios_append_mesh_var (struct adios_mesh_var_list_struct ** root
                                  ,struct adios_mesh_var_list_struct * var
                                  )
{
    while (root)
    {
        if (!*root)
        {
            *root = var;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}

static void adios_append_mesh_cell_list
                           (struct adios_mesh_cell_list_list_struct ** root
                           ,struct adios_mesh_cell_list_list_struct * cell_list
                           )
{
    while (root)
    {
        if (!*root)
        {
            *root = cell_list;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}


// dimensions is a comma separated list of dimension magnitudes
static int parseMeshUniformDimensions (const char * dimensions
                                      ,struct adios_group_struct * new_group
                                      ,struct adios_mesh_uniform_struct * mesh
                                      )
{
    char * c;  // comma location
    char * d1; // save of strdup
    char * tmp;
    struct adios_mesh_item_list_struct * item = 0;

    if (!dimensions)
    {
        fprintf (stderr, "config.xml: mesh uniform dimensions value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (dimensions);

    c = d1;
    tmp = c;

    while (c && *c)
    {
        if (*c == ',')
        {
            item = (struct adios_mesh_item_list_struct *) malloc
                            (sizeof (struct adios_mesh_item_list_struct));
            item->next = 0;

            if (!item)
            {
                fprintf (stderr, "Out of memory parseMeshUniformDimensions\n");
                free (d1);

                return 0;
            }

            *c = '\0';
            if (adios_int_is_num (tmp))
            {
                item->item.rank = strtod (tmp, 0);
                item->item.var = 0;
            }
            else
            {
                item->item.rank = 0.0;
                item->item.var =
                    adios_find_var_by_name (new_group->vars, tmp
                                           ,new_group->all_unique_var_names
                                           );
                if (!item->item.var)
                {
                    fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                            ,tmp
                            );
                    free (d1);

                    return 0;
                }
            }
            adios_append_mesh_item (&(mesh->dimensions), item);
            tmp = c + 1;
        }
        else
            c++;
    }

    free (d1);

    return 1;
}

static int parseMeshUniformOrigin (const char * origin
                                  ,struct adios_group_struct * new_group
                                  ,struct adios_mesh_uniform_struct * mesh
                                  )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_list_struct * item = 0;

    if (!origin)
    {
        fprintf (stderr, "config.xml: mesh uniform origin value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (origin);

    c = strtok (d1, ",");

    while (c)
    {
        item = (struct adios_mesh_item_list_struct *) malloc
                            (sizeof (struct adios_mesh_item_list_struct));
        item->next = 0;

        if (!item)
        {
            fprintf (stderr, "Out of memory parseMeshUniformOrigin\n");
            free (d1);

            return 0;
        }

        if (adios_int_is_var (c))
        {
            item->item.rank = 0.0;
            item->item.var =
                    adios_find_var_by_name (new_group->vars, c
                                           ,new_group->all_unique_var_names
                                           );
            if (!item->item.var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            item->item.rank = strtod (c, 0);
            item->item.var = 0;
        }

        adios_append_mesh_item (&(mesh->origin), item);

        c = strtok (NULL, ",");
    }

    free (d1);

    return 1;
}

static int parseMeshUniformSpacing (const char * spacing
                                   ,struct adios_group_struct * new_group
                                   ,struct adios_mesh_uniform_struct * mesh
                                   )
{
    char * c;  // comma location
    char * d1; // save of strdup
    char * tmp;
    struct adios_mesh_item_list_struct * item = 0;

    if (!spacing)
    {
        fprintf (stderr, "config.xml: mesh uniform spacing value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (spacing);

    c = d1;
    tmp = c;

    while (c && *c)
    {
        if (*c == ',')
        {
            item = (struct adios_mesh_item_list_struct *) malloc
                            (sizeof (struct adios_mesh_item_list_struct));
            item->next = 0;

            if (!item)
            {
                fprintf (stderr, "Out of memory parseMeshUniformSpacing\n");
                free (d1);

                return 0;
            }

            *c = '\0';
            if (adios_int_is_num (tmp))
            {
                item->item.rank = strtod (tmp, 0);
                item->item.var = 0;
            }
            else
            {
                item->item.rank = 0.0;
                item->item.var =
                    adios_find_var_by_name (new_group->vars, tmp
                                           ,new_group->all_unique_var_names
                                           );
                if (!item->item.var)
                {
                    fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                            ,tmp
                            );
                    free (d1);

                    return 0;
                }
            }
            adios_append_mesh_item (&(mesh->spacing), item);
            tmp = c + 1;
        }
        else
            c++;
    }

    free (d1);

    return 1;
}

static int parseMeshRectilinearDimensions (const char * dimensions
                                          ,struct adios_group_struct * new_group
                                          ,struct adios_mesh_rectilinear_struct * mesh
                                          )
{
    char * c;  // comma location
    char * d1; // save of strdup
    char * tmp;
    struct adios_mesh_item_list_struct * item = 0;

    if (!dimensions)
    {
        fprintf (stderr, "config.xml: mesh rectilinear dimensions value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (dimensions);

    c = d1;
    tmp = c;

    while (c && *c)
    {
        if (*c == ',')
        {
            item = (struct adios_mesh_item_list_struct *) malloc
                            (sizeof (struct adios_mesh_item_list_struct));
            item->next = 0;

            if (!item)
            {
                fprintf (stderr, "Out of memory parseMeshRectilinearDimensions\n");
                free (d1);

                return 0;
            }

            *c = '\0';
            if (adios_int_is_num (tmp))
            {
                item->item.rank = strtod (tmp, 0);
                item->item.var = 0;
            }
            else
            {
                item->item.rank = 0.0;
                item->item.var =
                    adios_find_var_by_name (new_group->vars, tmp
                                           ,new_group->all_unique_var_names
                                           );
                if (!item->item.var)
                {
                    fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                            ,tmp
                            );
                    free (d1);

                    return 0;
                }
            }
            adios_append_mesh_item (&(mesh->dimensions), item);
            tmp = c + 1;
        }
        else
            c++;
    }

    free (d1);

    return 1;
}

static int parseMeshRectilinearCoordinatesMultiVar (const char * coordinates
                                                   ,struct adios_group_struct * new_group
                                                   ,struct adios_mesh_rectilinear_struct * mesh
                                                   )
{
    char * c;  // comma location
    char * d1; // save of strdup
    char * tmp;
    struct adios_mesh_var_list_struct * var = 0;

    if (!coordinates)
    {
        fprintf (stderr, "config.xml: mesh rectilinear coordinates-multi-var "
                         "value required\n"
                );

        return 0;
    }

    d1 = strdup (coordinates);

    c = d1;
    tmp = c;

    while (c && *c)
    {
        if (*c == ',')
        {
            var = (struct adios_mesh_var_list_struct *) malloc
                            (sizeof (struct adios_mesh_var_list_struct));
            var->next = 0;

            if (!var)
            {
                fprintf (stderr, "Out of memory "
                                 "parseMeshRectilinearCoordinatesMultiVar\n"
                        );
                free (d1);

                return 0;
            }

            *c = '\0';
            if (adios_int_is_var (tmp))
            {
                var->var =
                    adios_find_var_by_name (new_group->vars, tmp
                                           ,new_group->all_unique_var_names
                                           );
                if (!var->var)
                {
                    fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                            ,tmp
                            );
                    free (d1);

                    return 0;
                }
            }
            adios_append_mesh_var (&(mesh->coordinates), var);
            tmp = c + 1;
        }
        else
            c++;
    }

    free (d1);

    return 1;
}

static int parseMeshRectilinearCoordinatesSingleVar (const char * coordinates
                                                    ,struct adios_group_struct * new_group
                                                    ,struct adios_mesh_rectilinear_struct * mesh
                                                    )
{
    char * c;  // comma location
    char * d1; // save of strdup
    char * tmp;
    struct adios_mesh_var_list_struct * var = 0;

    if (!coordinates)
    {
        fprintf (stderr, "config.xml: mesh rectilinear coordinates-single-var "
                         "value required\n"
                );

        return 0;
    }

    d1 = strdup (coordinates);

    c = d1;
    tmp = c;

    while (c && *c)
    {
        if (*c == ',')
        {
            var = (struct adios_mesh_var_list_struct *) malloc
                            (sizeof (struct adios_mesh_var_list_struct));
            var->next = 0;

            if (!var)
            {
                fprintf (stderr, "Out of memory "
                                 "parseMeshRectilinearCoordinatesSingleVar\n"
                        );
                free (d1);

                return 0;
            }

            *c = '\0';
            if (adios_int_is_var (tmp))
            {
                var->var =
                    adios_find_var_by_name (new_group->vars, tmp
                                           ,new_group->all_unique_var_names
                                           );
                if (!var->var)
                {
                    fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                            ,tmp
                            );
                    free (d1);

                    return 0;
                }
            }
            adios_append_mesh_var (&(mesh->coordinates), var);
            tmp = c + 1;
        }
        else
            c++;
    }

    free (d1);

    return 1;
}

static int parseMeshStructuredNspace (const char * nspace
                                     ,struct adios_group_struct * new_group
                                     ,struct adios_mesh_structured_struct * mesh
                                     )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_struct * item = 0;

    if (!nspace)
    {
        fprintf (stderr, "config.xml: mesh structured nspace value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (nspace);

    c = strtok (d1, ",");

    while (c)
    {
        item = (struct adios_mesh_item_struct *) malloc
                            (sizeof (struct adios_mesh_item_struct));

        if (!item)
        {
            fprintf (stderr, "Out of memory parseMeshStructuredNspace\n");
            free (d1);

            return 0;
        }

        if (adios_int_is_var (c))
        {
            item->rank = 0.0;
            item->var =
                    adios_find_var_by_name (new_group->vars, c
                                           ,new_group->all_unique_var_names
                                           );
            if (!item->var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            item->rank = strtod (c, 0);
            item->var = 0;
        }

        mesh->nspace = item;

        c = strtok (NULL, ",");
    }

    free (d1);

    return 1;
}

static int parseMeshStructuredDimensions (const char * dimensions
                                         ,struct adios_group_struct * new_group
                                         ,struct adios_mesh_structured_struct * mesh
                                         )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_list_struct * item = 0;

    if (!dimensions)
    {
        fprintf (stderr, "config.xml: mesh structured dimensions value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (dimensions);

    c = strtok (d1, ",");

    while (c)
    {
        item = (struct adios_mesh_item_list_struct *) malloc
                            (sizeof (struct adios_mesh_item_list_struct));
        item->next = 0;

        if (!item)
        {
            fprintf (stderr, "Out of memory parseMeshStructuredDimensions\n");
            free (d1);

            return 0;
        }

        if (adios_int_is_var (c))
        {
            item->item.rank = 0.0;
            item->item.var =
                    adios_find_var_by_name (new_group->vars, c
                                           ,new_group->all_unique_var_names
                                           );
            if (!item->item.var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            item->item.rank = strtod (c, 0);
            item->item.var = 0;
        }

        adios_append_mesh_item (&(mesh->dimensions), item);

        c = strtok (NULL, ",");
    }

    free (d1);

    return 1;
}

static int parseMeshStructuredPointsMultiVar (const char * points
                                             ,struct adios_group_struct * new_group
                                             ,struct adios_mesh_structured_struct * mesh
                                             )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_var_list_struct * var = 0;

    if (!points)
    {
        fprintf (stderr, "config.xml: mesh structured points-multi-var value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (points);

    c = strtok (d1, ",");

    while (c)
    {
        var = (struct adios_mesh_var_list_struct *) malloc
                        (sizeof (struct adios_mesh_var_list_struct));
        var->next = 0;

        if (!var)
        {
            fprintf (stderr, "Out of memory parseMeshStructuredPointsMultiVar\n");
            free (d1);

            return 0;
        }

        if (adios_int_is_var (c))
        {
            var->var =
                    adios_find_var_by_name (new_group->vars, c
                                           ,new_group->all_unique_var_names
                                           );
            if (!var->var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            var->var = 0;
        }

        adios_append_mesh_var (&(mesh->points), var);

        c = strtok (NULL, ",");
    }

    free (d1);

    return 1;
}

static int parseMeshStructuredPointsSingleVar (const char * points
                                              ,struct adios_group_struct * new_group
                                              ,struct adios_mesh_structured_struct * mesh
                                              )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_var_list_struct * var = 0;

    if (!points)
    {
        fprintf (stderr, "config.xml: mesh structured points-single-var value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (points);

    c = strtok (d1, ",");

    while (c)
    {
        var = (struct adios_mesh_var_list_struct *) malloc
                        (sizeof (struct adios_mesh_var_list_struct));
        var->next = 0;

        if (!var)
        {
            fprintf (stderr, "Out of memory parseMeshStructuredPointsSingleVar\n");
            free (d1);

            return 0;
        }

        if (adios_int_is_var (c))
        {
            var->var =
                    adios_find_var_by_name (new_group->vars, c
                                           ,new_group->all_unique_var_names
                                           );
            if (!var->var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            var->var = 0;
        }

        adios_append_mesh_var (&(mesh->points), var);

        c = strtok (NULL, ",");
    }

    free (d1);

    return 1;
}

static int parseMeshUnstructuredPoints (const char * components
                                       ,const char * number_of_points
                                       ,const char * value
                                       ,struct adios_group_struct * new_group
                                       ,struct adios_mesh_unstructured_struct * mesh
                                       )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_var_struct * var = 0;

    if (!components)
    {
        fprintf (stderr, "config.xml: mesh structured points-single-var value "
                         "required\n"
                );

        return 0;
    }

    d1 = strdup (components);

    c = strtok (d1, ",");

    while (c)
    {
        if (adios_int_is_var (c))
        {
            var = adios_find_var_by_name (new_group->vars, c
                                         ,new_group->all_unique_var_names
                                         );
            if (!var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            var = 0;
        }

        mesh->points = var;

        c = strtok (NULL, ",");
    }

    free (d1);

    return 1;
}

static int parseMeshUnstructuredUniformCells (const char * count
                                             ,const char * data
                                             ,const char * type
                                             ,struct adios_group_struct * new_group
                                             ,struct adios_mesh_unstructured_struct * mesh
                                             )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_cell_list_list_struct * cell_list = 0;

    if (!count)
    {
        fprintf (stderr, "config.xml: mesh unstructured uniform-cells "
                         "count value required\n"
                );

        return 0;
    }
    if (!data)
    {
        fprintf (stderr, "config.xml: mesh unstructured uniform-cells "
                         "data value required\n"
                );

        return 0;
    }
    if (!type)
    {
        fprintf (stderr, "config.xml: mesh unstructured uniform-cells "
                         "type value required\n"
                );

        return 0;
    }

    d1 = strdup (count);
    c = strtok (d1, ",");

    cell_list = (struct adios_mesh_cell_list_list_struct *) malloc
                        (sizeof (struct adios_mesh_cell_list_list_struct));
    cell_list->next = 0;

    if (!cell_list)
    {
        fprintf (stderr, "Out of memory parseMeshStructuredPointsSingleVar\n");
        free (d1);

        return 0;
    }

    while (c)
    {
        if (adios_int_is_var (c))
        {
            cell_list->cell_list.count.var =
                    adios_find_var_by_name (new_group->vars, c
                                           ,new_group->all_unique_var_names
                                           );
            cell_list->cell_list.count.rank = 0;
            if (!cell_list->cell_list.count.var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            cell_list->cell_list.count.var = 0;
            cell_list->cell_list.count.rank = strtod (c, 0);
        }

        c = strtok (NULL, ",");
    }
    free (d1);

    d1 = strdup (data);
    c = strtok (d1, ",");
    while (c)
    {
        if (adios_int_is_var (c))
        {
            cell_list->cell_list.data =
                 adios_find_var_by_name (new_group->vars, c
                                        ,new_group->all_unique_var_names
                                        );
            if (!cell_list->cell_list.data)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            cell_list->cell_list.data = 0;
        }

        c = strtok (NULL, ",");
    }
    free (d1);

    d1 = strdup (type);
    c = strtok (d1, ",");
    while (c)
    {
        if (adios_int_is_var (c))
        {
            cell_list->cell_list.type.var =
                     adios_find_var_by_name (new_group->vars, c
                                            ,new_group->all_unique_var_names
                                            );
            cell_list->cell_list.type.rank = 0;
            if (!cell_list->cell_list.type.var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = strtod (c, 0);
        }

        c = strtok (NULL, ",");
    }
    free (d1);

    adios_append_mesh_cell_list (&(mesh->cell_list), cell_list);

    return 1;
}

static int parseMeshUnstructuredMixedCells (const char * count
                                           ,const char * data
                                           ,const char * types
                                           ,struct adios_group_struct * new_group
                                           ,struct adios_mesh_unstructured_struct * mesh
                                           )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_cell_list_list_struct * cell_list = 0;

    if (!count)
    {
        fprintf (stderr, "config.xml: mesh unstructured uniform-cells "
                         "count value required\n"
                );

        return 0;
    }
    if (!data)
    {
        fprintf (stderr, "config.xml: mesh unstructured uniform-cells "
                         "data value required\n"
                );

        return 0;
    }
    if (!types)
    {
        fprintf (stderr, "config.xml: mesh unstructured uniform-cells "
                         "types value required\n"
                );

        return 0;
    }

    d1 = strdup (count);
    c = strtok (d1, ",");

    cell_list = (struct adios_mesh_cell_list_list_struct *) malloc
                        (sizeof (struct adios_mesh_cell_list_list_struct));
    cell_list->next = 0;

    if (!cell_list)
    {
        fprintf (stderr, "Out of memory parseMeshStructuredPointsSingleVar\n");
        free (d1);

        return 0;
    }

    while (c)
    {
        if (adios_int_is_var (c))
        {
            cell_list->cell_list.count.var =
                   adios_find_var_by_name (new_group->vars, c
                                          ,new_group->all_unique_var_names
                                          );
            cell_list->cell_list.count.rank = 0;
            if (!cell_list->cell_list.count.var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            cell_list->cell_list.count.var = 0;
            cell_list->cell_list.count.rank = strtod (c, 0);
        }

        c = strtok (NULL, ",");
    }
    free (d1);

    d1 = strdup (data);
    c = strtok (d1, ",");
    while (c)
    {
        if (adios_int_is_var (c))
        {
            cell_list->cell_list.data =
                   adios_find_var_by_name (new_group->vars, c
                                          ,new_group->all_unique_var_names
                                          );
            if (!cell_list->cell_list.data)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            cell_list->cell_list.data = 0;
        }

        c = strtok (NULL, ",");
    }
    free (d1);

    d1 = strdup (types);
    c = strtok (d1, ",");
    while (c)
    {
        if (adios_int_is_var (c))
        {
            cell_list->cell_list.type.var =
                   adios_find_var_by_name (new_group->vars, c
                                          ,new_group->all_unique_var_names
                                          );
            cell_list->cell_list.type.rank = 0;
            if (!cell_list->cell_list.type.var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,c
                        );
                free (d1);

                return 0;
            }
        }
        else
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = strtod (c, 0);
        }

        c = strtok (NULL, ",");
    }
    free (d1);

    adios_append_mesh_cell_list (&(mesh->cell_list), cell_list);

    return 1;
}

// primary mesh XML parsing
static int parseMeshUniform (mxml_node_t * node
                            ,struct adios_group_struct * new_group
                            ,struct adios_mesh_uniform_struct ** mesh
                            )
{
    mxml_node_t * n;
    int saw_dimensions = 0;
    int saw_origin = 0;
    int saw_spacing = 0;

    for (n = mxmlWalkNext (node, node, MXML_DESCEND)
        ;n
        ;n = mxmlWalkNext (n, node, MXML_DESCEND)
        )
    {
        if (n->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (n->value.element.name, "dimensions"))
        {
            const char * dimensions;

            if (saw_dimensions)
            {
                fprintf (stderr, "config.xml: only one dimensions "
                 "definition allowed per mesh sructured-points\n"
                        );

                return 0;
            }

            saw_dimensions = 1;
            dimensions = mxmlElementGetAttr (n, "value");

            if (!dimensions)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "dimensions required\n"
                        );

                return 0;
            }

            if (!parseMeshUniformDimensions (dimensions, new_group, *mesh))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "origin"))
        {
            const char * value;

            if (saw_origin)
            {
                fprintf (stderr, "config.xml: only one origin "
                                 "definition allowed per mesh uniform\n"
                        );

                return 0;
            }

            saw_origin = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "origin required\n"
                        );

                return 0;
            }

            if (!parseMeshUniformOrigin (value, new_group, *mesh))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "spacing"))
        {
            const char * value;

            if (saw_spacing)
            {
                fprintf (stderr, "config.xml: only one spacing "
                                 "definition allowed per mesh uniform\n"
                        );

                return 0;
            }

            saw_spacing = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "spacing required\n"
                        );

                return 0;
            }

            if (!parseMeshUniformSpacing (value, new_group, *mesh))
                return 0;
        } else
        {
            if (!strncmp (n->value.element.name, "!--", 3)) // a comment
            {
                continue;
            }
        }
    }

    if (!saw_dimensions)
    {
        fprintf (stderr, "config.xml: dimensions required on mesh "
                         "type=uniform\n"
                );

        return 0;
    }

    return 1;
}

static int parseMeshRectilinear (mxml_node_t * node
                               ,struct adios_group_struct * new_group
                               ,struct adios_mesh_rectilinear_struct ** mesh
                               )
{
    mxml_node_t * n;
    int saw_dimensions = 0;
    int saw_coordinates_multi_var = 0;
    int saw_coordinates_single_var = 0;

    for (n = mxmlWalkNext (node, node, MXML_DESCEND)
        ;n
        ;n = mxmlWalkNext (n, node, MXML_DESCEND)
        )
    {
        if (n->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (n->value.element.name, "dimensions"))
        {
            const char * value;

            if (saw_dimensions)
            {
                fprintf (stderr, "config.xml: only one dimensions "
                                 "definition allowed per mesh rectilinear\n"
                        );

                return 0;
            }

            saw_dimensions = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "dimensions required\n"
                        );

                return 0;
            }

            if (!parseMeshRectilinearDimensions (value, new_group, *mesh))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "coordinates-multi-var"))
        {
            const char * value;

            if (saw_coordinates_multi_var || saw_coordinates_single_var)
            {
                fprintf (stderr, "config.xml: only one coordinates "
                                 "definition allowed per mesh rectilinear\n"
                        );

                return 0;
            }

            saw_coordinates_multi_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "coordinates-multi-var required\n"
                        );

                return 0;
            }

            if (!parseMeshRectilinearCoordinatesMultiVar (value, new_group, *mesh))
                return 0;
            (*mesh)->coordinates_single_var = adios_flag_no;
        } else
        if (!strcasecmp (n->value.element.name, "coordinates-single-var"))
        {
            const char * value;

            if (saw_coordinates_single_var || saw_coordinates_multi_var)
            {
                fprintf (stderr, "config.xml: only one coordinates "
                                 "definition allowed per mesh rectilinear\n"
                        );

                return 0;
            }

            saw_coordinates_single_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "coordinates-single-var required\n"
                        );

                return 0;
            }

            if (!parseMeshRectilinearCoordinatesSingleVar (value, new_group, *mesh))
                return 0;
            (*mesh)->coordinates_single_var = adios_flag_yes;
        } else
        {
            if (!strncmp (n->value.element.name, "!--", 3)) // a comment
            {
                continue;
            }
        }
    }

    if (!saw_dimensions)
    {
        fprintf (stderr, "config.xml: dimensions required on mesh "
                         "type=rectilinear\n"
                );

        return 0;
    }
    if (!saw_coordinates_multi_var && !saw_coordinates_single_var)
    {
        fprintf (stderr, "config.xml: coordinates-multi-var or "
                         "coordinates-single-var required on mesh "
                         "type=rectilinear\n"
                );

        return 0;
    }

    return 1;
}

static int parseMeshStructured (mxml_node_t * node
                               ,struct adios_group_struct * new_group
                               ,struct adios_mesh_structured_struct ** mesh
                               )
{
    mxml_node_t * n;
    int saw_nspace = 0;
    int saw_dimensions = 0;
    int saw_points_multi_var = 0;
    int saw_points_single_var = 0;

    for (n = mxmlWalkNext (node, node, MXML_DESCEND)
        ;n
        ;n = mxmlWalkNext (n, node, MXML_DESCEND)
        )
    {
        if (n->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (n->value.element.name, "nspace"))
        {
            const char * value;

            if (saw_nspace)
            {
                fprintf (stderr, "config.xml: only one nspace "
                                 "definition allowed per mesh structured\n"
                        );

                return 0;
            }

            saw_nspace = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "nspace required\n"
                        );

                return 0;
            }

            if (!parseMeshStructuredNspace (value, new_group, *mesh))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "dimensions"))
        {
            const char * value;

            if (saw_dimensions)
            {
                fprintf (stderr, "config.xml: only one dimensions "
                                 "definition allowed per mesh structured\n"
                        );

                return 0;
            }

            saw_dimensions = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "dimensions required\n"
                        );

                return 0;
            }

            if (!parseMeshStructuredDimensions (value, new_group, *mesh))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "points-multi-var"))
        {
            const char * value;

            if (saw_points_multi_var || saw_points_single_var)
            {
                fprintf (stderr, "config.xml: only one points "
                                 "definition allowed per mesh structured\n"
                        );

                return 0;
            }

            saw_points_multi_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "points-multi-var required\n"
                        );

                return 0;
            }

            if (!parseMeshStructuredPointsMultiVar (value, new_group, *mesh))
                return 0;
            (*mesh)->points_single_var = adios_flag_no;
        } else
        if (!strcasecmp (n->value.element.name, "points-single-var"))
        {
            const char * value;

            if (saw_points_multi_var || saw_points_single_var)
            {
                fprintf (stderr, "config.xml: only one points "
                                 "definition allowed per mesh structured\n"
                        );

                return 0;
            }

            saw_points_single_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "points-single-var required\n"
                        );

                return 0;
            }

            if (!parseMeshStructuredPointsSingleVar (value, new_group, *mesh))
                return 0;
            (*mesh)->points_single_var = adios_flag_yes;
        } else
        {
            if (!strncmp (n->value.element.name, "!--", 3)) // a comment
            {
                continue;
            }
        }
    }

    if (!saw_dimensions)
    {
        fprintf (stderr, "config.xml: dimensions required on mesh "
                         "type=structured\n"
                );

        return 0;
    }
    if (!saw_points_multi_var && !saw_points_single_var)
    {
        fprintf (stderr, "config.xml: points-single-var or points-multi-var "
                         "required on mesh type=structured\n"
                );

        return 0;
    }

    return 1;
}

static int parseMeshUnstructured (mxml_node_t * node
                                 ,struct adios_group_struct * new_group
                                 ,struct adios_mesh_unstructured_struct ** mesh
                                 )
{
    mxml_node_t * n;
    int saw_points = 0;
    int saw_cell_set = 0;

    for (n = mxmlWalkNext (node, node, MXML_DESCEND)
        ;n
        ;n = mxmlWalkNext (n, node, MXML_DESCEND)
        )
    {
        if (n->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (n->value.element.name, "points"))
        {
            const char * components;
            const char * number_of_points;
            const char * value;

            if (saw_points)
            {
                fprintf (stderr, "config.xml: only one points "
                                 "definition allowed per mesh unstructured\n"
                        );

                return 0;
            }

            saw_points = 1;
            components = mxmlElementGetAttr (n, "components");
            number_of_points = mxmlElementGetAttr (n, "number-of-points");
            value = mxmlElementGetAttr (n, "value");

            if (!components)
            {
                fprintf (stderr, "config.xml: components attribute on "
                                 "points required\n"
                        );

                return 0;
            }
            if (!number_of_points)
            {
                fprintf (stderr, "config.xml: number-of-points attribute on "
                                 "points required\n"
                        );

                return 0;
            }
            if (!value)
            {
                fprintf (stderr, "config.xml: value attribute on "
                                 "points required\n"
                        );

                return 0;
            }

            if (!parseMeshUnstructuredPoints (components, number_of_points
                                             ,value, new_group, *mesh
                                             )
               )
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "uniform-cells"))
        {
            const char * count;
            const char * data;
            const char * type;

            saw_cell_set = 1;
            count = mxmlElementGetAttr (n, "count");
            data = mxmlElementGetAttr (n, "data");
            type = mxmlElementGetAttr (n, "type");

            if (!count)
            {
                fprintf (stderr, "config.xml: count attribute on "
                                 "uniform-cells required\n"
                        );

                return 0;
            }
            if (!data)
            {
                fprintf (stderr, "config.xml: data attribute on "
                                 "uniform-cells required\n"
                        );

                return 0;
            }
            if (!type)
            {
                fprintf (stderr, "config.xml: type attribute on "
                                 "uniform-cells required\n"
                        );

                return 0;
            }

            if (!parseMeshUnstructuredUniformCells (count, data, type
                                                   ,new_group, *mesh
                                                   )
               )
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "mixed-cells"))
        {
            const char * count;
            const char * data;
            const char * types;

            saw_cell_set = 1;
            count = mxmlElementGetAttr (n, "count");
            data = mxmlElementGetAttr (n, "data");
            types = mxmlElementGetAttr (n, "types");

            if (!count)
            {
                fprintf (stderr, "config.xml: count attribute on "
                                 "mixed-cells required\n"
                        );

                return 0;
            }
            if (!data)
            {
                fprintf (stderr, "config.xml: data attribute on "
                                 "mixed-cells required\n"
                        );

                return 0;
            }
            if (!types)
            {
                fprintf (stderr, "config.xml: types attribute on "
                                 "mixed-cells required\n"
                        );

                return 0;
            }

            if (!parseMeshUnstructuredMixedCells (count, data, types
                                                   ,new_group, *mesh
                                                   )
               )
                return 0;
        } else
        {
            if (!strncmp (n->value.element.name, "!--", 3)) // a comment
            {
                continue;
            }
        }
    }

    if (!saw_points)
    {
        fprintf (stderr, "config.xml: grid-nodes required on mesh "
                         "type=unstructured\n"
                );

        return 0;
    }
    if (!saw_cell_set)
    {
        fprintf (stderr, "config.xml: at least one cell-set required on "
                         "mesh type=unstructured\n"
                );

        return 0;
    }

    return 1;
}

static int validatePath (const struct adios_var_struct * vars
                        ,const char * test_path
                        )
{
    // if it is a default path, it is ok by default
    if (!strcmp (test_path, "/"))
    {
        return 1;
    }

    char * path = strdup (test_path);
    int len = strlen (path);
    char * path_only;
    char * var_only;
    char * last_slash = strrchr (path, '/'); // find the last '/'
    path_only = (char *) malloc (len + 1);
    var_only = (char *) malloc (len + 1);
    memset (path_only, 0, len + 1);
    memset (var_only, 0, len + 1);
    if (last_slash == path + len - 1)  // if it is a trailing '/', remove
    {
        last_slash = '\0';
        last_slash = strrchr (path, '/');
    }
    if (last_slash == 0)
    {
        strcpy (path_only, "/");
        strcpy (var_only, path);
    }
    else
    {
        strncpy (path_only, path, (last_slash - path));
        strncpy (var_only, last_slash + 1, (len - (last_slash - path + 1)));
    }

    while (vars)
    {
        int path_only_len = strlen (path_only);
        int var_path_len = strlen (vars->path);
        int var_name_len = strlen (vars->name);
        char full_path_matches = (!strcasecmp (vars->path, path));
        char path_matches = (!strcasecmp (vars->path, path_only));
        char var_matches = (!strcasecmp (vars->name, var_only));
        char prefix_matches = 0;
        char * path_var;
        path_var = (char *) malloc (var_path_len + var_name_len + 2);
        sprintf (path_var, "%s/%s", vars->path, vars->name);
        char path_var_matches = (!strcasecmp (path_var, path));

        if (var_path_len >= len && path_only_len > 0)
            prefix_matches = (!strncmp (vars->path, path_only, path_only_len));
        if (!prefix_matches)
            prefix_matches = (!strncmp (vars->path, path, len));
        int var_len = strlen (var_only);

        if (   (path_matches && var_matches)
            || (path_matches && var_len == 0)
            || (full_path_matches)
            || (prefix_matches)
            || (path_var_matches)
           )
        {
            free (path);
            free (path_only);
            free (var_only);
            free (path_var);

            return 1;
        }
        vars = vars->next;
        free (path_var);
    }

    // not found
    free (path);
    free (path_only);
    free (var_only);

    return 0;
}

static int parseGroup (mxml_node_t * node)
{
    mxml_node_t * n;

    const char * datagroup_name = 0;
    const char * coordination_comm = 0;
    const char * coordination_var = 0;
    const char * host_language = 0;
    const char * time_index_name = 0;
    const char * stats = 0;

    int64_t      ptr_new_group;
    struct adios_group_struct * new_group;
    enum ADIOS_FLAG host_language_fortran = adios_flag_yes, enable_stats = adios_flag_yes;
    int i;

    for (i = 0; i < node->value.element.num_attrs; i++)
    {
        mxml_attr_t * attr = &node->value.element.attrs [i];

        GET_ATTR("name",attr,datagroup_name,"adios-group")
        // JL: 1-2010
        // Although this is not used, we aer leaving in the retrevial
        // of this to avoid messages from all of the existing XML files.
        // In a few months, once everything has been updated, we can remove
        // this code
        GET_ATTR("coordination-communicator",attr,coordination_comm,"adios-group")
        GET_ATTR("coordination-var",attr,coordination_var,"adios-group")
        GET_ATTR("host-language",attr,host_language,"adios-group")
        GET_ATTR("time-index",attr,time_index_name,"adios-group")
        GET_ATTR("stats",attr,stats,"adios-group")
        fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                         "(ignored)\n"
                ,attr->name
                ,"adios-group"
                );
    }

    if (!datagroup_name)
    {
        fprintf (stderr,
                 "config.xml: name attribute required on adios-group\n");

        return 0;
    }
    if (!host_language)
    {
        host_language_fortran = adios_host_language_fortran;
    }
    else
    {
        if (!strcasecmp (host_language, "Fortran"))
        {
            host_language_fortran = adios_flag_yes;
        }
        else
        {
            if (!strcasecmp (host_language, "C"))
            {
                host_language_fortran = adios_flag_no;
            }
            else
            {
                fprintf (stderr, "config.xml: invalid host-language %s"
                        ,host_language
                        );

                return 0;
            }
        }
    }

    if (!stats)
    {
        enable_stats = adios_flag_yes;
    }
    else
    {
        if (!strcasecmp (stats, "On"))
        {
            enable_stats = adios_flag_yes;
        }
        else if (!strcasecmp (stats, "Off"))
        {
            enable_stats = adios_flag_no;
        }
        else
        {
            fprintf (stderr, "config.xml, invalid stats %s"
                    ,stats
                    );
            return 0;
        }
    }

// fix the bgp bugs 
/*
    adios_common_declare_group ((int64_t *) &new_group, datagroup_name
                               ,host_language_fortran, coordination_comm
                               ,coordination_var, time_index_name
                               );
*/
    adios_common_declare_group (&ptr_new_group, datagroup_name
                               ,host_language_fortran, coordination_comm
                               ,coordination_var, time_index_name
                               ,enable_stats
                               );
     new_group = (struct adios_group_struct *)ptr_new_group;

    for (n = mxmlWalkNext (node, node, MXML_DESCEND)
        ;n
        ;n = mxmlWalkNext (n, node, MXML_NO_DESCEND)
        )
    {
        const char * gb_global_dimensions = "";
        const char * gb_local_offsets = "";

        if (n->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (n->value.element.name, "var"))
        {
            const char * name = 0;
            const char * path = 0;
            const char * type = 0;
            const char * dimensions = 0;
            const char * dimension = 0;
            const char * gread = 0;
            const char * gwrite = 0;
            const char * read_flag = 0;
            enum ADIOS_DATATYPES t1;

            for (i = 0; i < n->value.element.num_attrs; i++)
            {
                mxml_attr_t * attr = &n->value.element.attrs [i];

                GET_ATTR("name",attr,name,"var")
                GET_ATTR("path",attr,path,"var")
                GET_ATTR("type",attr,type,"var")
                GET_ATTR("dimensions",attr,dimensions,"var")
                GET_ATTR("dimension",attr,dimension,"var")
                GET_ATTR("gwrite",attr,gwrite,"var")
                GET_ATTR("gread",attr,gread,"var")
                GET_ATTR("read",attr,read_flag,"var")
                fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"var"
                        );
            }

            if (!name)
                name = "";  // this will catch the error
            if (!path)
                path = "/";
            if (!type)
                type = ""; // this will catch the error
            t1 = parseType (type, name);

            if (!dimensions)
            {
                dimensions = dimension;
                if (!dimensions)
                    dimensions = "";
            }

            if (read_flag)
                parseFlag ("read", read_flag, adios_flag_no);

// fix the bgp bugs
//            if (!adios_common_define_var (*(int64_t *) &new_group, name
            if (!adios_common_define_var (ptr_new_group, name
                                         ,path, t1, dimensions
                                         ,gb_global_dimensions
                                         ,gb_local_offsets
                                         )
               )
            {
                return 0;
            }
        } else
        if (!strcasecmp (n->value.element.name, "global-bounds"))
        {
            mxml_node_t * n1;   // used for global_bounds
            struct adios_global_bounds_struct * new_global_bounds = 0;

            const char * dimensions = 0;
            const char * dimension = 0;
            const char * global_dimensions = 0;
            const char * global_dimension = 0;
            const char * offsets = 0;
            const char * offset = 0;
            const char * local_offsets = 0;
            const char * local_offset = 0;

            for (i = 0; i < n->value.element.num_attrs; i++)
            {
                mxml_attr_t * attr = &n->value.element.attrs [i];

                GET_ATTR("dimensions",attr,dimensions,"var")
                GET_ATTR("dimension",attr,dimension,"var")

                GET_ATTR("global-dimensions",attr,global_dimensions,"var")
                GET_ATTR("global-dimension",attr,global_dimension,"var")

                GET_ATTR("offsets",attr,offsets,"var")
                GET_ATTR("offset",attr,offset,"var")

                GET_ATTR("local-offsets",attr,local_offsets,"var")
                GET_ATTR("local-offset",attr,local_offset,"var")

                fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"global-bounds"
                        );
            }

            if (!dimensions)
            {
                dimensions = (dimension ? dimension : global_dimensions);
                dimensions = (dimensions ? dimensions : global_dimension);
                if (!dimensions)
                {
                    fprintf (stderr, "config.xml: dimensions required on "
                                     "global-bounds\n"
                            );

                    return 0;
                }
            }
            if (!offsets)
            {
                offsets = (offset ? offset : local_offsets);
                offsets = (offsets ? offsets : local_offset);
                fprintf (stderr, "config.xml: offsets required on "
                                 "global-bounds\n"
                        );

                return 0;
            }

            gb_global_dimensions = dimensions;
            gb_local_offsets = offsets;

            for (n1 = mxmlWalkNext (n, n, MXML_DESCEND)
                ;n1
                ;n1 = mxmlWalkNext (n1, n, MXML_DESCEND)
                )
            {
                if (n1->type != MXML_ELEMENT)
                {
                    continue;
                }

                if (!strcasecmp (n1->value.element.name, "var"))
                {
                    const char * name = 0;
                    const char * path = 0;
                    const char * type = 0;
                    const char * dimension = 0;
                    const char * dimensions = 0;
                    const char * gwrite = 0;
                    const char * gread = 0;
                    const char * read_flag = 0;
                    enum ADIOS_DATATYPES t1;

                    for (i = 0; i < n1->value.element.num_attrs; i++)
                    {
                        mxml_attr_t * attr = &n1->value.element.attrs [i];

                        GET_ATTR("name",attr,name,"var")
                        GET_ATTR("path",attr,path,"var")
                        GET_ATTR("type",attr,type,"global-bounds var")
                        GET_ATTR("dimensions",attr,dimensions,"var")
                        GET_ATTR("dimension",attr,dimension,"var")
                        GET_ATTR("gwrite",attr,gwrite,"var")
                        GET_ATTR("gread",attr,gread,"var")
                        GET_ATTR("read",attr,read_flag,"var")
                        fprintf (stderr, "config.xml: unknown attribute '%s' "
                                         "on %s (ignored)\n"
                                ,attr->name
                                ,"var"
                                );
                    }

                    if (!name)
                        name = "";  // this will catch the error
                    if (!path)
                        path = "/";
                    if (!type)
                        type = ""; // this will catch the error
                    t1 = parseType (type, name);
                    if (!dimensions)
                        dimensions = dimension;

                    if (read_flag)
                        parseFlag ("read", read_flag, adios_flag_no);
// fix the bgp bugs
//                    if (!adios_common_define_var (*(int64_t *) &new_group
                    if (!adios_common_define_var (ptr_new_group
                                                 ,name
                                                 ,path, t1, dimensions
                                                 ,gb_global_dimensions
                                                 ,gb_local_offsets
                                                 )
                       )
                    {
                        return 0;
                    }
                } else
                {
                    if (!strncmp (n1->value.element.name, "!--", 3)) // comment
                    {
                        continue;
                    }
                    else
                    {
                        fprintf (stderr, "config.xml: invalid xml element: "
                                         "'%s'\n"
                                ,n1->value.element.name
                                );

                        return 0;
                    }
                }
            }

            gb_global_dimensions = "";
            gb_local_offsets = "";
        } else
        if (!strcasecmp (n->value.element.name, "attribute"))
        {
            const char * name = 0;
            const char * path = 0;
            const char * value = 0;
            const char * type = 0;
            const char * var = 0;
            enum ADIOS_DATATYPES t1;

            for (i = 0; i < n->value.element.num_attrs; i++)
            {
                mxml_attr_t * attr = &n->value.element.attrs [i];

                GET_ATTR("name",attr,name,"var")
                GET_ATTR("path",attr,path,"var")
                GET_ATTR("type",attr,type,"var")
                GET_ATTR("value",attr,value,"var")
                GET_ATTR("var",attr,var,"var")
                fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"attribute"
                        );
            }

            if (!name)
            {
                fprintf (stderr, "config.xml: attribute element requires "
                                 "name\n");

                return 0;
            }
            if (!path)
            {
                fprintf (stderr, "config.xml: attribute element requires "
                                 "path\n");

                return 0;
            }
            if ((!value && !var) || (value && var))
            {
                fprintf (stderr, "config.xml: attriute element '%s' "
                                 "requires either value OR var\n"
                        ,name
                        );

                return 0;
            }
            if (var && type)
            {
                fprintf (stderr, "config.xml: attriute element '%s'. "
                                 "The type of an associated var is part "
                                 "of the associated var element and cannot "
                                 "be provided as part of the attribute "
                                 "element."
                                 "\n"
                        ,name
                        );

                return 0;
            }
            if (!type && value)
            {
                type = "string";
            }
            if (!var)
            {
                t1 = parseType (type, name);
            }
            else
            {
                t1 = adios_unknown;
            }
// fix the bgp bugs
//            if (!adios_common_define_attribute (*(int64_t *) &new_group, name
            if (!adios_common_define_attribute (ptr_new_group, name
                                               ,path, t1, value, var
                                               )
               )
            {
                return 0;
            }
        } else
        if (!strcasecmp (n->value.element.name, "mesh"))
        {
            const char * type;
            const char * time_varying;

            new_group->mesh = (struct adios_mesh_struct *) malloc
                                           (sizeof (struct adios_mesh_struct));
            type = mxmlElementGetAttr (n, "type");
            time_varying = mxmlElementGetAttr (n, "time-varying");

            if (!type)
                type = "";

            if (!strcasecmp (type, "uniform"))
            {
                new_group->mesh->type = ADIOS_MESH_UNIFORM;
                new_group->mesh->uniform =
                    (struct adios_mesh_uniform_struct *)
                         calloc (1, sizeof (struct adios_mesh_uniform_struct));
                parseMeshUniform (n, new_group, &new_group->mesh->uniform);
            } else
            if (!strcasecmp (type, "structured"))
            {
                new_group->mesh->type = ADIOS_MESH_STRUCTURED;
                new_group->mesh->structured =
                    (struct adios_mesh_structured_struct *)
                       calloc (1, sizeof (struct adios_mesh_structured_struct));
                parseMeshStructured (n, new_group
                                    ,&new_group->mesh->structured);
            } else
            if (!strcasecmp (type, "rectilinear"))
            {
                new_group->mesh->type = ADIOS_MESH_RECTILINEAR;
                new_group->mesh->rectilinear =
                    (struct adios_mesh_rectilinear_struct *)
                      calloc (1, sizeof (struct adios_mesh_rectilinear_struct));
                parseMeshRectilinear (n, new_group
                                     ,&new_group->mesh->rectilinear);
            } else
            if (!strcasecmp (type, "unstructured"))
            {
                new_group->mesh->type = ADIOS_MESH_UNSTRUCTURED;
                new_group->mesh->unstructured =
                    (struct adios_mesh_unstructured_struct *)
                     calloc (1, sizeof (struct adios_mesh_unstructured_struct));

                parseMeshUnstructured (n, new_group
                                      ,&new_group->mesh->unstructured);
            } else
            {
                fprintf (stderr, "config.xml: invalid mesh type: '%s'\n"
                        ,type
                        );

                return 0;
            }

            new_group->mesh->time_varying = parseFlag
                               ("time-varying", time_varying, adios_flag_no);

            if (new_group->mesh->time_varying == adios_flag_unknown)
            {
                return 1;
            }
        } else
        if (!strcasecmp (n->value.element.name, "gwrite"))
        {
            const char * src = 0;

            for (i = 0; i < n->value.element.num_attrs; i++)
            {
                mxml_attr_t * attr = &n->value.element.attrs [i];

                GET_ATTR("src",attr,src,"var")
                fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"gwrite"
                        );
            }
            if (!src)
            {
                fprintf (stderr, "config.xml: gwrite element requires "
                                 "src\n");

                return 0;
            }
        } else
        {
            if (!strncmp (n->value.element.name, "!--", 3)) // a comment
            {
                continue;
            }
            else
            {
                fprintf (stderr, "config.xml: invalid xml element: '%s'\n"
                        ,n->value.element.name
                        );

                return 0;
            }
        }
    }

    // now that we have declared the whole group, validate that the
    // paths specified in attributes refer to real things or give
    // a warning
    struct adios_attribute_struct * a = new_group->attributes;
    while (a)
    {
        if (!validatePath (new_group->vars, a->path))
        {
/*
            fprintf (stderr, "config.xml warning: attribute element '%s' "
                             "has path '%s' that does not match "
                             "any var path or name.\n"
                    ,a->name, a->path
                    );
*/
        }

        a = a->next;
    }

    return 1;
}

static int parseAnalysis (mxml_node_t * node)
{
    mxml_node_t * n;

    const char * group = 0;
    const char * var = 0;
    const char * bin_intervals = 0;
    const char * bin_count = 0;
    const char * bin_min = 0;
    const char * bin_max = 0;

    int i;
    int64_t group_id;
    struct adios_group_struct * g;

    for (i = 0; i < node->value.element.num_attrs; i++)
    {
        mxml_attr_t * attr = &node->value.element.attrs [i];

        GET_ATTR("adios-group",attr,group,"analysis")
        GET_ATTR("var",attr,var,"analysis")
        GET_ATTR("break-points",attr,bin_intervals,"analysis")
        GET_ATTR("min",attr,bin_min,"analysis")
        GET_ATTR("max",attr,bin_max,"analysis")
        GET_ATTR("count",attr,bin_count,"analysis")
        fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                         "(ignored)\n"
                ,attr->name
                ,"method"
                );
    }

    if (!var)
    {
        fprintf (stderr, "config.xml: variable name must be given\n");
        return 0;
    }

    if (!group)
    {
        fprintf (stderr, "config.xml: adios-group name must be given\n");
        return 0;
    }

    adios_common_get_group (&group_id, group);
    g = (struct adios_group_struct *) group_id;

    if (!g)
    {
        fprintf (stderr, "config.xml: Didn't find group %s for analysis\n", group);
        return 0;
    }
    if(!adios_common_define_var_characteristics(g, var, bin_intervals, bin_min, bin_max, bin_count))
        return 0;

    return 1;
}

static int parseMethod (mxml_node_t * node)
{
    mxml_node_t * n;

    const char * priority = 0;
    const char * iterations = 0;
    const char * base_path = 0;
    const char * method = 0;
    const char * group = 0;
    const char * parameters = 0;
    int p1;
    int i1;
    int i;

    for (i = 0; i < node->value.element.num_attrs; i++)
    {
        mxml_attr_t * attr = &node->value.element.attrs [i];

        GET_ATTR("priority",attr,priority,"method")
        GET_ATTR("iterations",attr,iterations,"method")
        GET_ATTR("base-path",attr,base_path,"method")
        GET_ATTR("method",attr,method,"method")
        GET_ATTR("group",attr,group,"method")
        fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                         "(ignored)\n"
                ,attr->name
                ,"method"
                );
    }

    // Check for parameters, if they exist
    n = mxmlWalkNext (node, node, MXML_DESCEND);
    if (n != NULL)
    {
        parameters = n->value.text.string;
    }
    else
    {
        parameters = NULL;
    }

    if (!priority)
        p1 = 1;
    else
        p1 = atoi (priority);
    if (!iterations)
        i1 = 1;
    else
        i1 = atoi (iterations);
    if (!parameters)
        parameters = "";
    if (!base_path)
        base_path = "";
    else
    {
        uint16_t len = strlen (base_path);
        if (len > 0 && base_path [len - 1] != '/')
        {
            fprintf (stderr, "config.xml: method %s for group %s base-path "
                             "must end with a '/' character\n"
                    ,method, group
                    );

            return 0;
        }
    }
    if (!group)
        group = "";
    if (!method)
        method = "";

    if (!adios_common_select_method (p1, method, parameters, group
                                    ,base_path, i1
                                    )
       )
    {
        return 0;
    }

    return 1;
}

static int parseBuffer (mxml_node_t * node)
{
    const char * size_MB = 0;
    const char * free_memory_percentage = 0;
    const char * allocate_time = 0;

    int i;

    int size = -1;

    for (i = 0; i < node->value.element.num_attrs; i++)
    {
        mxml_attr_t * attr = &node->value.element.attrs [i];

        GET_ATTR("size-MB",attr,size_MB,"method")
        GET_ATTR("free-memory-percentage",attr,free_memory_percentage,"method")
        GET_ATTR("allocate-time",attr,allocate_time,"method")
        fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                         "(ignored)\n"
                ,attr->name
                ,"buffer"
                );
    }



    if ((!size_MB && !free_memory_percentage) || !allocate_time)
    {
        fprintf (stderr, "config.xml: must define allocate-time and either "
                         "size-MB or free-memory-percentage for "
                         "buffer element\n"
                );

        return 0;
    }
    else
    {
        if (!strcasecmp (allocate_time, "now"))
        {
            adios_buffer_alloc_when_set (ADIOS_BUFFER_ALLOC_NOW);
        }
        else
        {
            if (!strcasecmp (allocate_time, "oncall"))
            {
                adios_buffer_alloc_when_set (ADIOS_BUFFER_ALLOC_LATER);
            }
            else
            {
                fprintf (stderr, "config.xml: buffer allocate-time %s "
                                 "invalid. ('now' or 'oncall')\n"
                        ,allocate_time
                        );

                return 0;
            }
        }

        if (size_MB)
        {
            adios_buffer_alloc_percentage_set (0);
            size = atoi (size_MB);
            if (size_MB == 0)
            {
                fprintf (stderr, "config.xml: buffer size-MB is either 0 or "
                                 "cannot be parsed: %s"
                        ,size_MB
                        );

                return 0;
            }

            if (size < 1)
                size = 1; // we need a minimum 1 MB buffer

            adios_buffer_size_requested_set ((uint64_t) size * 1024 * 1024);
        }
        else
        {
            adios_buffer_alloc_percentage_set (1);
            size = atoi (free_memory_percentage);
            if (size > 0 && size <= 100)
            {
                adios_buffer_size_requested_set ((uint64_t) size);
            }
            else
            {
                fprintf (stderr, "config.xml: buffer free-memory-percentage %s "
                                 "is not an integer between 1 and 100\n"
                        ,free_memory_percentage
                        );

                return 0;
            }
        }

        if (adios_buffer_alloc_when_get() == ADIOS_BUFFER_ALLOC_NOW)
        {
            return adios_set_buffer_size ();
        }
    }

    return 1;
}

int adios_parse_config (const char * config)
{
    FILE * fp = 0;
    mxml_node_t * doc = NULL;
    mxml_node_t * node = NULL;
    mxml_node_t * root = NULL;
    int saw_datagroup = 0;
    int saw_method = 0;
    int saw_buffer = 0;

    if (!adios_transports_initialized)
    {
        adios_transports_initialized = 1;
        adios_init_transports (&adios_transports);
    }

    char * buffer = NULL;
//#if HAVE_MPI
    int buffer_size = 0;
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
//#endif
        fp = fopen (config, "r");
        if (!fp)
        {
            fprintf (stderr, "missing config file %s\n", config);

            return 0;
        }
        struct stat s;
        if (stat (config, &s) == 0)
        {
            buffer = malloc (s.st_size + 1);
            buffer [s.st_size] = 0;
        }
        if (buffer)
        {
            size_t bytes_read = fread (buffer, 1, s.st_size, fp);
            if (bytes_read != s.st_size)
            {
                fprintf (stderr, "error reading config file: %s. Expected %d Got %d\n"
                        ,config, s.st_size, bytes_read );

                return 0;
            }
        }
        else
        {
            fprintf (stderr, "error allocating %d for reading config.\n"
                    ,s.st_size + 1
                    );

            return 0;
        }
        fclose (fp);
//#if HAVE_MPI
        buffer_size = s.st_size;
        MPI_Bcast (&buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast (buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast (&buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        buffer = malloc (buffer_size + 1);
        if (!buffer)
        {
            fprintf (stderr, "cannot allocate %d bytes to receive config file\n"
                    ,buffer_size + 1
                    );

            return 0;
        }
        MPI_Bcast (buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD);
        buffer [buffer_size] = 0;
    }
//#endif

    doc = mxmlLoadString (NULL, buffer, MXML_TEXT_CALLBACK);
    free (buffer);
    buffer = NULL;

    if (!doc)
    {
        fprintf (stderr, "config.xml: unknown error parsing XML "
                         "(probably structural)\n"
                         "Did you remember to start the file with\n"
                         "<?xml version=\"1.0\"?>\n");

        return 0;
    }

    root = doc;

    while (root && root->type != MXML_ELEMENT)
    {
        root = mxmlWalkNext (root, doc, MXML_DESCEND);
    }

    while (!strncmp (root->value.element.name, "!--", 3))
    {
        root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
        root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
    }

    if (strcasecmp (root->value.element.name, "adios-config"))
    {
        if (strncmp (root->value.element.name, "?xml", 4))
        {
            fprintf (stderr, "config.xml: invalid root xml element: %s\n"
                    ,root->value.element.name
                    );

            mxmlRelease (doc);

            return 0;
        }
        else
        {
            while (!strncmp (root->value.element.name, "!--", 3))
            {
                root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
            }

            root = mxmlWalkNext (root, doc, MXML_DESCEND);  // skip ver num
            root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);  // get next
            while (!strncmp (root->value.element.name, "!--", 3))
            {
                root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
                root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
            }
        }
    }
    else
    {
        //printf ("it is adios-config\n");
    }


    if (strcasecmp (root->value.element.name, "adios-config"))
    {
        fprintf (stderr, "config.xml: invalid root xml element: %s\n"
                ,root->value.element.name
                );

        mxmlRelease (doc);

        return 0;
    }
    else
    {
        const char * host_language = 0;
        int i;

        for (i = 0; i < root->value.element.num_attrs; i++)
        {
            mxml_attr_t * attr = &root->value.element.attrs [i];

            GET_ATTR("host-language",attr,host_language,"var")
            fprintf (stderr, "config.xml: unknown attribute '%s' on %s "
                             "(ignored)\n"
                    ,attr->name
                    ,"adios-config"
                    );
        }

        if (!host_language)
        {
            host_language = "Fortran";
        }

        if (!strcasecmp (host_language, "Fortran"))
        {
            adios_host_language_fortran = adios_flag_yes;
        }
        else
        {
            if (!strcasecmp (host_language, "C"))
            {
                adios_host_language_fortran = adios_flag_no;
            }
            else
            {
                fprintf (stderr, "config.xml: invalid host-language %s"
                        ,host_language
                        );

                mxmlRelease (doc);

                return 0;
            }
        }
    }

    for (node = mxmlWalkNext (root, doc, MXML_DESCEND_FIRST)
        ;node
        ;node = mxmlWalkNext (node, root, MXML_NO_DESCEND)
        )
    {
        if (node->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcasecmp (node->value.element.name, "adios-group"))
        {
            if (!parseGroup (node))
                break;
            saw_datagroup = 1;
        }
        else
        {
            if (   !strcasecmp (node->value.element.name, "transport")
                || !strcasecmp (node->value.element.name, "method")
               )
            {
                if (!parseMethod (node))
                    break;
                saw_method = 1;
            }
            else
            {
                if (!strcasecmp (node->value.element.name, "buffer"))
                {
                    if (!parseBuffer (node))
                        break;
                    saw_buffer = 1;
                }
				else
				{
                    if (!strcasecmp (node->value.element.name, "analysis"))
                    {
                        if (!parseAnalysis(node))
                            break;
                    }
                	else
                	{
                	    if (!strncmp (node->value.element.name, "!--", 3))
                	    {
                	        continue;
                	    }
                	    else
                	    {
                	        fprintf (stderr, "config.xml: invalid element: %s\n"
                	                ,node->value.element.name
                	                );

                	        break;
                	    }
					}
                }
            }
        }
    }

    mxmlRelease (doc);

    if (!saw_datagroup)
    {
        fprintf (stderr, "config.xml: must define at least 1 adios-group in "
                         "config.xml\n"
                );

        return 0;
    }
    if (!saw_method)
    {
        fprintf (stderr, "config.xml: must define at least 1 method for "
                         "the adios-group in config.xml\n"
                );

        return 0;
    }
    if (!saw_buffer)
    {
        fprintf (stderr, "config.xml: must define the buffer element in "
                         "config.xml\n"
                );

        return 0;
    }

    return 1;
}

int adios_local_config ()
{
    if (!adios_transports_initialized)
    {
        adios_transports_initialized = 1;
        adios_init_transports (&adios_transports);
    }

    return 1;
}

int adios_common_select_method (int priority, const char * method
                               ,const char * parameters, const char * group
                               ,const char * base_path, int iters
                               )
{
    int64_t group_id;
    struct adios_group_struct * g;
    struct adios_method_struct * new_method;
    int requires_group_comm = 0;

    new_method = (struct adios_method_struct *)
                           malloc (sizeof (struct adios_method_struct));
    
    new_method->m = ADIOS_METHOD_UNKNOWN;
    new_method->base_path = strdup (base_path);
    new_method->method = strdup (method);
    new_method->parameters = strdup (parameters);
    new_method->iterations = iters;
    new_method->priority = priority;
    new_method->method_data = 0;
    new_method->group = 0;

    if (adios_parse_method (method, &new_method->m, &requires_group_comm))
    {
        if (   new_method->m != ADIOS_METHOD_UNKNOWN
            && new_method->m != ADIOS_METHOD_NULL
            && adios_transports [new_method->m].adios_init_fn
           )
        {
            adios_transports [new_method->m].adios_init_fn
                                       (parameters, new_method);
        }
    }
    else
    {
        fprintf (stderr, "config.xml: invalid transport: %s\n", method);

        free (new_method->base_path);
        free (new_method->method);
        free (new_method->parameters);
        free (new_method);

        return 0;
    }

    adios_common_get_group (&group_id, group);
    g = (struct adios_group_struct *) group_id;
    if (!g)
    {
        fprintf (stderr, "config.xml: Didn't find group: %s for transport: %s\n"
                ,group, method
                );

        free (new_method->base_path);
        free (new_method->method);
        free (new_method->parameters);
        free (new_method);

        return 0;
    }
    else
    {
        // JL: 1-2010
        // we no longer require this since we moved the coordiantion
        // communicator to the open call. Leaving code here, just in case.
        // once this has been validated thoroughly, remove this block of code.
        //if (requires_group_comm && !g->group_comm)
        //{
        //    fprintf (stderr, "config.xml: method %s for group %s.  Group does "
        //                     "not have the required coordination-communicator"
        //                     ".\n"
        //            ,method, group
        //            );

        //    free (new_method->base_path);
        //    free (new_method->method);
        //    free (new_method->parameters);
        //    free (new_method);

        //    return 0;
        //}
        adios_add_method_to_group (&g->methods, new_method);
        new_method->group = g;
    }

    adios_append_method (new_method);

    return 1;
}

int adios_common_select_method_by_group_id (int priority, const char * method
                                           ,const char * parameters, int64_t group_id
                                           ,const char * base_path, int iters
                                           )
{
    struct adios_group_struct * g;
    struct adios_method_struct * new_method;
    int requires_group_comm = 0;

    new_method = (struct adios_method_struct *)
                           malloc (sizeof (struct adios_method_struct));

    new_method->m = ADIOS_METHOD_UNKNOWN;
    new_method->base_path = strdup (base_path);
    new_method->method = strdup (method);
    new_method->parameters = strdup (parameters);
    new_method->iterations = iters;
    new_method->priority = priority;
    new_method->method_data = 0;
    new_method->group = 0;

    if (adios_parse_method (method, &new_method->m, &requires_group_comm))
    {
        if (   new_method->m != ADIOS_METHOD_UNKNOWN
            && new_method->m != ADIOS_METHOD_NULL
            && adios_transports [new_method->m].adios_init_fn
           )
        {
            adios_transports [new_method->m].adios_init_fn
                                       (parameters, new_method);
        }
    }
    else
    {
        fprintf (stderr, "config.xml: invalid transport: %s\n", method);

        free (new_method->base_path);
        free (new_method->method);
        free (new_method->parameters);
        free (new_method);

        return 0;
    }

    g = (struct adios_group_struct *) group_id;
    if (!g)
    {
        fprintf (stderr, "config.xml: invalid group id: %llu for transport: %s\n"
                ,group_id, method
                );

        free (new_method->base_path);
        free (new_method->method);
        free (new_method->parameters);
        free (new_method);

        return 0;
    }
    else
    {
        if (requires_group_comm && !g->group_comm)
        {
            fprintf (stderr, "config.xml: method %s for group %s.  Group does "
                             "not have the required coordination-communicator"
                             ".\n"
                    ,method, g->name
                    );

            free (new_method->base_path);
            free (new_method->method);
            free (new_method->parameters);
            free (new_method);

            return 0;
        }
        adios_add_method_to_group (&g->methods, new_method);
        new_method->group = g;
    }

    adios_append_method (new_method);

    return 1;
}

void adios_cleanup ()
{
    adios_transports_initialized = 0;
    if (adios_transports)
        free (adios_transports);
    adios_transports = 0;

    while (adios_methods)
    {
        struct adios_method_list_struct * methods = adios_methods->next;
        if (adios_methods->method->base_path)
            free (adios_methods->method->base_path);
        if (adios_methods->method->method)
            free (adios_methods->method->method);
        if (adios_methods->method->method_data)
            free (adios_methods->method->method_data);
        if (adios_methods->method->parameters)
            free (adios_methods->method->parameters);
        free (adios_methods->method);
        free (adios_methods);
        adios_methods = methods;
    }

    while (adios_groups)
    {
        struct adios_group_list_struct * groups = adios_groups->next;

        if (adios_groups->group->name)
            free (adios_groups->group->name);

        while (adios_groups->group->vars)
        {
            struct adios_var_struct * vars = adios_groups->group->vars->next;

            if (adios_groups->group->vars->name)
                free (adios_groups->group->vars->name);
            if (adios_groups->group->vars->path)
                free (adios_groups->group->vars->path);

            while (adios_groups->group->vars->dimensions)
            {
                struct adios_dimension_struct * dimensions
                                = adios_groups->group->vars->dimensions->next;

                free (adios_groups->group->vars->dimensions);
                adios_groups->group->vars->dimensions = dimensions;
            }

			// NCSU - Clean up stat
            if (adios_groups->group->vars->stats)
			{
				int j, idx;
				int c, count = 1;

				if (adios_groups->group->vars->type == adios_complex || adios_groups->group->vars->type == adios_double_complex)
					count = 3;

				for (c = 0; c < count; c ++)
				{
					j = idx = 0;
					while (adios_groups->group->vars->bitmap >> j)
					{
						if (adios_groups->group->vars->bitmap >> j & 1)
						{
                            if (j == adios_statistic_hist)
                            {
								struct adios_index_characteristics_hist_struct * hist = (struct adios_index_characteristics_hist_struct *) adios_groups->group->vars->stats[c][idx].data;
								free (hist->breaks);
								free (hist->frequencies);
								free (hist);
							}
						    else
								free (adios_groups->group->vars->stats[c][idx].data);
							idx ++;
						}
						j ++;
					}
					free (adios_groups->group->vars->stats[c]);
				}

                free (adios_groups->group->vars->stats);
			}
            if (adios_groups->group->vars->data)
                free (adios_groups->group->vars->data);

            free (adios_groups->group->vars);
            adios_groups->group->vars = vars;
        }

        while (adios_groups->group->attributes)
        {
            struct adios_attribute_struct * attributes
                                        = adios_groups->group->attributes->next;

            if (adios_groups->group->attributes->name)
                free (adios_groups->group->attributes->name);
            if (adios_groups->group->attributes->path)
                free (adios_groups->group->attributes->path);
            if (adios_groups->group->attributes->value)
                free (adios_groups->group->attributes->value);

            free (adios_groups->group->attributes);
            adios_groups->group->attributes = attributes;
        }

        if (adios_groups->group->group_comm)
            free (adios_groups->group->group_comm);
        if (adios_groups->group->group_by)
            free (adios_groups->group->group_by);
        if (adios_groups->group->time_index_name)
            free (adios_groups->group->time_index_name);

        while (adios_groups->group->methods)
        {
            struct adios_method_list_struct * m = adios_groups->group->methods->next;
            free (adios_groups->group->methods);
            adios_groups->group->methods = m;
        }

        if (adios_groups->group->mesh)
        {
            switch (adios_groups->group->mesh->type)
            {
                case ADIOS_MESH_UNIFORM:
                {
                    struct adios_mesh_item_list_struct * i;
                    while (adios_groups->group->mesh->uniform->dimensions)
                    {
                        i = adios_groups->group->mesh->uniform->dimensions->next;
                        free (adios_groups->group->mesh->uniform->dimensions);
                        adios_groups->group->mesh->uniform->dimensions = i;
                    }
                    while (adios_groups->group->mesh->uniform->origin)
                    {
                        i = adios_groups->group->mesh->uniform->origin->next;
                        free (adios_groups->group->mesh->uniform->origin);
                        adios_groups->group->mesh->uniform->origin = i;
                    }
                    while (adios_groups->group->mesh->uniform->spacing)
                    {
                        i = adios_groups->group->mesh->uniform->spacing->next;
                        free (adios_groups->group->mesh->uniform->spacing);
                        adios_groups->group->mesh->uniform->spacing = i;
                    }

                    break;
                }

                case ADIOS_MESH_STRUCTURED:
                {
                    struct adios_mesh_item_list_struct * i;
                    struct adios_mesh_var_list_struct * v;
                    while (adios_groups->group->mesh->structured->dimensions)
                    {
                        i = adios_groups->group->mesh->structured->dimensions->next;
                        free (adios_groups->group->mesh->structured->dimensions);
                        adios_groups->group->mesh->structured->dimensions = i;
                    }
                    while (adios_groups->group->mesh->structured->points)
                    {
                        v = adios_groups->group->mesh->structured->points->next;
                        free (adios_groups->group->mesh->structured->points);
                        adios_groups->group->mesh->structured->points = v;
                    }
                    if (adios_groups->group->mesh->structured->nspace)
                        free (adios_groups->group->mesh->structured->nspace);

                    break;
                }

                case ADIOS_MESH_RECTILINEAR:
                {
                    struct adios_mesh_item_list_struct * i;
                    struct adios_mesh_var_list_struct * v;
                    while (adios_groups->group->mesh->rectilinear->dimensions)
                    {
                        i = adios_groups->group->mesh->rectilinear->dimensions->next;
                        free (adios_groups->group->mesh->rectilinear->dimensions);
                        adios_groups->group->mesh->rectilinear->dimensions = i;
                    }
                    while (adios_groups->group->mesh->rectilinear->coordinates)
                    {
                        v = adios_groups->group->mesh->rectilinear->coordinates->next;
                        free (adios_groups->group->mesh->rectilinear->coordinates);
                        adios_groups->group->mesh->rectilinear->coordinates = v;
                    }

                    break;
                }

                case ADIOS_MESH_UNSTRUCTURED:
                {
                    if (adios_groups->group->mesh->unstructured->components)
                        free (adios_groups->group->mesh->unstructured->components);
                    if (adios_groups->group->mesh->unstructured->points_count)
                        free (adios_groups->group->mesh->unstructured->points_count);
                    while (adios_groups->group->mesh->unstructured->cell_list)
                    {
                        struct adios_mesh_cell_list_list_struct * next
                          = adios_groups->group->mesh->unstructured->cell_list->next;
                        free (adios_groups->group->mesh->unstructured->cell_list);
                        adios_groups->group->mesh->unstructured->cell_list = next;
                    }

                    break;
                }
            }

            free (adios_groups->group->mesh);
        }

        free (adios_groups->group);
        free (adios_groups);
        adios_groups = groups;
    }
}


