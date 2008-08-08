#include <math.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <stdint.h>

// xml parser
#include <mxml.h>

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

enum ADIOS_FLAG adios_host_language_fortran = adios_flag_yes;

// buffer sizing may be problematic.  To get a more accurate picture, check:
// http://chandrashekar.info/vault/linux-system-programs.html
static uint64_t adios_buffer_size_requested = 0;
static uint64_t adios_buffer_size_max = 0;
static uint64_t adios_buffer_size_remaining = 0;
static int adios_buffer_alloc_percentage = 0;  // 1 = yes, 0 = no
static enum ADIOS_BUFFER_ALLOC_WHEN {ADIOS_BUFFER_ALLOC_UNKNOWN
                                    ,ADIOS_BUFFER_ALLOC_NOW
                                    ,ADIOS_BUFFER_ALLOC_LATER
                                    } adios_buffer_alloc_when
                                         = ADIOS_BUFFER_ALLOC_UNKNOWN;

static struct adios_method_list_struct * adios_methods = 0;
static struct adios_group_list_struct * adios_groups = 0;

struct adios_transport_struct * adios_transports = 0;
static int adios_transports_initialized = 0;

// this macro makes getting the attributes easier
#define GET_ATTR(n,attr,var,en)                              \
if (!strcasecmp (n, attr->name))                             \
    if (!var)                                                \
    {                                                        \
        var = attr->value;                                   \
        continue;                                            \
    }                                                        \
    else                                                     \
    {                                                        \
        fprintf (stderr, "config.xml: duplicate attribute "  \
                         n                                   \
                         " on "                              \
                         en                                  \
                         " (ignored)");                      \
        continue;                                            \
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

static int is_var (const char * temp) // 1 == yes, 0 == no
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

static int is_num (char * temp) // 1 == yes, 0 == no
{
    char * extra = 0;

    strtod (temp, &extra);

    if (extra)
        return 0;
    else
        return 1;

    return 0;
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
            if (is_num (tmp))
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

        if (is_var (c))
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
            if (is_num (tmp))
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
            if (is_num (tmp))
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
            if (is_var (tmp))
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
            if (is_var (tmp))
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

        if (is_var (c))
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

        if (is_var (c))
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

        if (is_var (c))
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

        if (is_var (c))
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
        if (is_var (c))
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
        if (is_var (c))
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
        if (is_var (c))
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
        if (is_var (c))
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
        if (is_var (c))
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
        if (is_var (c))
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
        if (is_var (c))
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

        if (var_path_len >= len)
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
    const char * time_index = 0;

    struct adios_group_struct * new_group;
    enum ADIOS_FLAG host_language_fortran = adios_flag_yes;
    int i;

    for (i = 0; i < node->value.element.num_attrs; i++)
    {
        mxml_attr_t * attr = &node->value.element.attrs [i];

        GET_ATTR("name",attr,datagroup_name,"adios-group")
        GET_ATTR("coordination-communicator",attr,coordination_comm,"adios-group")
        GET_ATTR("coordination-var",attr,coordination_var,"adios-group")
        GET_ATTR("host-language",attr,host_language,"adios-group")
        GET_ATTR("time-index",attr,time_index,"adios-group")
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

    adios_common_declare_group ((long long *) &new_group, datagroup_name
                               ,host_language_fortran, coordination_comm
                               ,coordination_var, time_index
                               );

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
            const char * gname = 0;
            enum ADIOS_DATATYPES t1;

            for (i = 0; i < n->value.element.num_attrs; i++)
            {
                mxml_attr_t * attr = &n->value.element.attrs [i];

                GET_ATTR("name",attr,name,"var")
                GET_ATTR("path",attr,path,"var")
                GET_ATTR("type",attr,type,"var")
                GET_ATTR("dimensions",attr,dimensions,"var")
                GET_ATTR("dimension",attr,dimension,"var")
                GET_ATTR("gname",attr,gname,"var")
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

            if (!adios_common_define_var (*(long long *) &new_group, name
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
                    const char * gname = 0;
                    enum ADIOS_DATATYPES t1;

                    for (i = 0; i < n1->value.element.num_attrs; i++)
                    {
                        mxml_attr_t * attr = &n1->value.element.attrs [i];

                        GET_ATTR("name",attr,name,"var")
                        GET_ATTR("path",attr,path,"var")
                        GET_ATTR("type",attr,type,"global-bounds var")
                        GET_ATTR("dimensions",attr,dimensions,"var")
                        GET_ATTR("dimension",attr,dimension,"var")
                        GET_ATTR("gname",attr,gname,"var")
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

                    if (!adios_common_define_var (*(long long *) &new_group
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
            if (!validatePath (new_group->vars, path))
            {
                fprintf (stderr, "config.xml: attribute element '%s' "
                                 "has path '%s' that does not match "
                                 "an existing var path or name.\n"
                        ,name, path
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

            if (!adios_common_define_attribute (*(long long *) &new_group, name
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
            adios_buffer_alloc_when = ADIOS_BUFFER_ALLOC_NOW;
        }
        else
        {
            if (!strcasecmp (allocate_time, "oncall"))
            {
                adios_buffer_alloc_when = ADIOS_BUFFER_ALLOC_LATER;
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
            adios_buffer_alloc_percentage = 0;
            size = atoi (size_MB);
            if (size_MB == 0)
            {
                fprintf (stderr, "config.xml: buffer size-MB is either 0 or "
                                 "cannot be parsed: %s"
                        ,size_MB
                        );

                return 0;
            }

            adios_buffer_size_requested = (  (uint64_t) size
                                           * 1024
                                           * 1024
                                          );
        }
        else
        {
            adios_buffer_alloc_percentage = 1;
            size = atoi (free_memory_percentage);
            if (size > 0 && size <= 100)
            {
                adios_buffer_size_requested = (uint64_t) size;
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

        if (adios_buffer_alloc_when == ADIOS_BUFFER_ALLOC_NOW)
        {
            return adios_set_buffer_size ();
        }
    }

    return 1;
}

int adios_set_buffer_size ()
{
    if (!adios_buffer_size_max) // not called before
    {
        long pagesize;
        long pages;

        pagesize = sysconf (_SC_PAGE_SIZE);
        pages = sysconf (_SC_AVPHYS_PAGES);

        if (adios_buffer_alloc_percentage)
        {
            adios_buffer_size_max =   (pages * pagesize / 100.0)
                                    * adios_buffer_size_requested;
        }
        else
        {
            if (pagesize * pages >= adios_buffer_size_requested)
            {
                // sufficient memory, do nothing
                adios_buffer_size_max = adios_buffer_size_requested;
            }
            else
            {
                fprintf (stderr, "adios_allocate_buffer (): insufficient memory: "
                                 "%llu requested, %llu available.  Using "
                                 "available.\n"
                        ,adios_buffer_size_requested
                        ,(((uint64_t) pagesize) * pages)
                        );
                 adios_buffer_size_max = ((uint64_t) pagesize) * pages;
           }
        }

        adios_buffer_size_remaining = adios_buffer_size_max;

        return 1;
    }
    else
    {
        fprintf (stderr, "adios_allocate_buffer already called. "
                         "No changes made.\n"
                );

        return 0;
    }
}

uint64_t adios_method_buffer_alloc (uint64_t size)
{
    if (adios_buffer_size_remaining >= size)
    {
        adios_buffer_size_remaining -= size;

        return size;
    }
    else
    {
        uint64_t remaining = adios_buffer_size_remaining;

        adios_buffer_size_remaining = 0;

        return remaining;
    }
}

int adios_method_buffer_free (uint64_t size)
{
    if (size + adios_buffer_size_remaining > adios_buffer_size_max)
    {
        fprintf (stderr, "ERROR: attempt to return more bytes to buffer "
                         "pool than were originally available\n"
                );

        adios_buffer_size_remaining = adios_buffer_size_max;

        return 0;
    }
    else
    {
        adios_buffer_size_remaining += size;

        return 1;
    }
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

void adios_parse_dimension (const char * dimension
                           ,const char * global_dimension
                           ,const char * local_offset
                           ,struct adios_group_struct * g
                           ,struct adios_dimension_struct * dim
                           )
{
    if (!dimension)
    {
        fprintf (stderr, "adios_parse_dimension: dimension not provided\n");

        return;
    }

    if (is_var (dimension))
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
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,dimension
                        );

                return;
            }
            else
            {
                if (attr->var)
                {
                    attr->var->is_dim = adios_flag_yes;
                }
                dim->dimension.id = attr->id;
            }
        }
        else
        {
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

        return;
    }

    if (is_var (global_dimension))
    {
        struct adios_var_struct * var = 0;
        dim->global_dimension.rank = 0;
        var = adios_find_var_by_name (g->vars, global_dimension
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
                fprintf (stderr, "config.xml: invalid global-bounds "
                                 "dimension: %s\n"
                        ,global_dimension
                        );

                return;
            }
            else
            {
                dim->global_dimension.id = attr->id;
            }
        }
        else
        {
            var->is_dim = adios_flag_yes;
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

        return;
    }

    if (is_var (local_offset))
    {
        struct adios_var_struct * var = 0;
        dim->local_offset.rank = 0;
        var = adios_find_var_by_name (g->vars, local_offset
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
                fprintf (stderr, "config.xml: invalid var local_offset: %s\n"
                        ,local_offset
                        );

                return;
            }
            else
            {
                dim->local_offset.id = attr->id;
            }
        }
        else
        {
            var->is_dim = adios_flag_yes;
        }
    }
    else
    {
        dim->local_offset.id = 0;
        dim->local_offset.rank = strtol (local_offset, NULL, 10);
    }
}

struct adios_method_list_struct * adios_get_methods ()
{
    return adios_methods;
}

struct adios_group_list_struct * adios_get_groups ()
{
    return adios_groups;
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

    fp = fopen (config, "r");
    if (!fp)
    {
        fprintf (stderr, "missing config file %s\n", config);

        return 0;
    }
    doc = mxmlLoadFile (NULL, fp, MXML_TEXT_CALLBACK);
    fclose (fp);
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
            if (!strcasecmp (node->value.element.name, "method"))
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
                                    ,adios_type_to_string (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = (void *) t;

                            return 1;
                        }
                    case adios_short:
                        if (t < SHRT_MIN || t > SHRT_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = (void *) t;

                            return 1;
                        }
                    case adios_integer:
                        if (t < INT_MIN || t > INT_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = (void *) t;

                            return 1;
                        }
                }
            }
        }
        case adios_long:
        {
            int errno_save = errno;
            long long t = strtoll (value, &end, 10);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "type is %s, value is out of range: '%s'\n"
                        ,adios_type_to_string (type), value
                        );

                return 0;
            }
            else
            {
                *out = (void *) t;

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
                                    ,adios_type_to_string (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = (void *) t;

                            return 1;
                        }
                    case adios_unsigned_short:
                        if (t > USHRT_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = (void *) t;

                            return 1;
                        }
                    case adios_unsigned_integer:
                        if (t > UINT_MAX)
                        {
                            fprintf (stderr, "type is %s, value "
                                             "is out of range: '%s'\n"
                                    ,adios_type_to_string (type), value
                                    );

                            return 0;
                        }
                        else
                        {
                            *out = (void *) t;

                            return 1;
                        }
                }
            }
        }
        case adios_unsigned_long:
        {
            int errno_save = errno;
            unsigned long long t = strtoull (value, &end, 10);
            if (errno != errno_save || (end != 0 && *end != '\0'))
            {
                fprintf (stderr, "type is %s, value is out of range: '%s'\n"
                        ,adios_type_to_string (type), value
                        );

                return 0;
            }
            else
            {
                *out = (void *) t;

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
                        ,adios_type_to_string (type), value
                        );

                return 0;
            }
            else
            {
                // there must be a better way to do this...
                // and get rid of the warning
                *out = (void *) *(uint32_t *) &t;

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
                        ,adios_type_to_string (type), value
                        );

                return 0;
            }
            else
            {
                *out = (void *) *(uint64_t *) &t;

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
                        ,adios_type_to_string (type), value
                        );

                return 0;
            }
            else
            {
                fprintf (stderr, "adios_parse_scalar_string long_double not "
                                 "supported\n"
                        );

                return 1;
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
#if 0
            int errno_save = errno;
            long t = strtol(const char *restrict, char **restrict, int);
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
                            fprintf (stderr, "type is signed byte, value "
                                             "is out of range: '%s'\n
                                    ,value
                                    );

                            return 0;
                        }
                        else
                        {
                            return 1;
                        }
                    case adios_short:
                        if (t < SHRT_MIN || t > SHRT_MAX)
                        {
                            fprintf (stderr, "type is signed short, value "
                                             "is out of range: '%s'\n
                                    ,value
                                    );

                            return 0;
                        }
                        else
                        {
                            return 1;
                        }
                    case adios_int:
                        if (t < INT_MIN || t > INT_MAX)
                        {
                            fprintf (stderr, "type is signed int, value "
                                             "is out of range: '%s'\n
                                    ,value
                                    );

                            return 0;
                        }
                        else
                        {
                            return 1;
                        }
                }
            }
#endif
        }
        case adios_double_complex:
        {
            fprintf (stderr, "adios_double_complex type validation needs to "
                             "be implemented\n");
            return 1;
#if 0
            int errno_save = errno;
            long t = strtol(const char *restrict, char **restrict, int);
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
                            fprintf (stderr, "type is signed byte, value "
                                             "is out of range: '%s'\n
                                    ,value
                                    );

                            return 0;
                        }
                        else
                        {
                            return 1;
                        }
                    case adios_short:
                        if (t < SHRT_MIN || t > SHRT_MAX)
                        {
                            fprintf (stderr, "type is signed short, value "
                                             "is out of range: '%s'\n
                                    ,value
                                    );

                            return 0;
                        }
                        else
                        {
                            return 1;
                        }
                    case adios_int:
                        if (t < INT_MIN || t > INT_MAX)
                        {
                            fprintf (stderr, "type is signed int, value "
                                             "is out of range: '%s'\n
                                    ,value
                                    );

                            return 0;
                        }
                        else
                        {
                            return 1;
                        }
                }
            }
#endif
        }

        case adios_unknown:
        default:
            fprintf (stderr, "unknown type cannot be validated\n");

            return 0;
    }

    return 1;
}

int adios_common_define_attribute (long long group, const char * name
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

            free (attr);

            return 0;
        }
    }

    attr->next = 0;

    adios_append_attribute (&g->attributes, attr, ++g->member_count);

    return 1;
}

void adios_extract_string (char * out, const char * in, int size)
{
    if (in && out)
        strcpy (out, in);
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
}

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
int adios_common_declare_group (long long * id, const char * name
                               ,enum ADIOS_FLAG host_language_fortran
                               ,const char * coordination_comm
                               ,const char * coordination_var
                               ,const char * time_index
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
    g->attributes = 0;
    g->group_by = (coordination_var ? strdup (coordination_var) : 0L);
    g->group_comm = (coordination_comm ? strdup (coordination_comm) : 0L);
    g->time_index = (time_index ? strdup (time_index) : 0L);
    g->process_id = 0;
    g->methods = 0;
    g->mesh = 0;

    *id = (long long) g;

    adios_append_group (g);

    return 1;
}

static void tokenize_dimensions (char * str, char *** tokens, int * count)
{
    char * t = str;
    char * save_str = strdup (str);
    int i;

    if (strlen (str) > 0)
        *count = 1;
    else
    {
        *tokens = 0;
        *count = 0;

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

int adios_common_define_var (long long group_id, const char * name
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

    if (strcmp (dim_temp, ""))
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
            
            adios_parse_dimension (dim, g_dim, lo_dim, t, d);

            adios_append_dimension (&v->dimensions, d);

            i++;
        }
        free (dim_temp);
        free (g_dim_temp);
        free (lo_dim_temp);
    }

    flag = adios_append_var (&t->vars, v, ++t->member_count);
    if (flag == adios_flag_no)
    {
        t->all_unique_var_names = adios_flag_no;
    }
    t->var_count++;

    return 1;
}

int adios_common_select_method (int priority, const char * method
                               ,const char * parameters, const char * group
                               ,const char * base_path, int iters
                               )
{
    long long group_id;
    struct adios_group_struct * g;
    struct adios_method_struct * new_method;

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

    if (adios_parse_method (method, &new_method->m))
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
        fprintf (stderr, "config.xml: invalid method: %s\n", method);

        return 0;
    }

    adios_common_get_group (&group_id, group);
    g = (struct adios_group_struct *) group_id;
    if (!g)
    {
        fprintf (stderr, "config.xml: Didn't find group: %s for method: %s\n"
                ,group, method
                );

        return 0;
    }
    else
    {
        adios_add_method_to_group (&g->methods, new_method);
    }

    adios_append_method (new_method);

    return 1;
}

void adios_common_get_group (long long * group_id, const char * name)
{
    struct adios_group_list_struct * g = adios_get_groups ();

    *group_id = 0;

    while (g)
    {
        if (!strcasecmp (g->group->name, name))
        {
            *group_id = (long long) g->group;

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
                         ,uint64_t * buffer_offset, void * data, uint64_t size
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
        if (d->dimension.id != 0)
            overhead += 2; // member id
        else
            overhead += 8; // value

        overhead += 1; // var flag
        if (d->global_dimension.id != 0)
            overhead += 2; // member id
        else
            overhead += 8; // value

        overhead += 1; // var flag
        if (d->local_offset.id != 0)
            overhead += 2; // member id
        else
            overhead += 8; // value

        d = d->next;
    }
    overhead += adios_get_type_size (v->type, ""); // min
    overhead += adios_get_type_size (v->type, ""); // max

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
    overhead += 2; // coordination comm id
    overhead += 2; // coordination var id
    overhead += 2; // timestep var id

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

    var = adios_find_var_by_name (g->vars, g->group_comm
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

    var = adios_find_var_by_name (g->vars, g->time_index
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

    struct adios_method_list_struct * m = fd->group->methods;
    uint8_t methods_count = 0;
    uint16_t methods_length = 0;
    while (m)
    {
        methods_count++;
        methods_length += 1 + 2 + strlen (m->method->parameters);

        m = m->next;
    }
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &methods_count, 1);
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &methods_length, 2);

    m = fd->group->methods;
    while (m)
    {
        uint16_t len = strlen (m->method->parameters);

        flag = (uint8_t) m->method->m;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &flag, 1);

        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &len, 2);

        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, m->method->parameters, len);

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
                if (  (*root)->entries_count + item->entries_count
                    > (*root)->entries_allocated
                   )
                {
                    int new_items = (item->entries_count == 1)
                                             ? 100 : item->entries_count;
                    (*root)->entries_allocated =   (*root)->entries_count
                                                 + new_items;
                    void * ptr;
                    ptr = realloc ((*root)->entries
                            ,  (*root)->entries_allocated
                             * sizeof (struct adios_index_var_entry_struct_v1)
                            );

                    if (ptr)
                    {
                        (*root)->entries = ptr;
                    }
                    else
                    {
                        fprintf (stderr, "error allocating memory to build "
                                         "var index.  Index aborted\n"
                                );

                        return;
                    }
                }
                memcpy (&(*root)->entries [(*root)->entries_count]
                       ,item->entries
                       ,  item->entries_count
                        * sizeof (struct adios_index_var_entry_struct_v1)
                       );

                (*root)->entries_count += item->entries_count;

                free (item->entries);
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
                          ,struct adios_index_process_group_struct_v1 * p2
                          ,struct adios_index_var_struct_v1 * v2
                          )
{
    // this will just add it on to the end and all should work fine
    index_append_process_group_v1 (p1, p2);

    // need to do vars one at a time to merge them properly
    struct adios_index_var_struct_v1 * v_temp;

    while (v2)
    {
        v_temp = v2->next;
        v2->next = 0;
        index_append_var_v1 (v1, v2);
        v2 = v_temp;
    }
}

static void adios_clear_process_groups_index_v1 (
                            struct adios_index_process_group_struct_v1 * root
                           )
{
    while (root)
    {
        struct adios_index_process_group_struct_v1 * temp = root->next;
        free (root);
        root = temp;
    }
}

static void adios_clear_vars_index_v1 (struct adios_index_var_struct_v1 * root)
{
    while (root)
    {
        struct adios_index_var_struct_v1 * temp = root->next;
        free (root->entries);
        free (root);
        root = temp;
    }
}

void adios_clear_index_v1 (struct adios_index_process_group_struct_v1 * pg_root
                          ,struct adios_index_var_struct_v1 * vars_root
                          )
{
    adios_clear_process_groups_index_v1 (pg_root);
    adios_clear_vars_index_v1 (vars_root);
}

void adios_build_index_v1 (struct adios_file_struct * fd
                       ,struct adios_index_process_group_struct_v1 ** pg_root
                       ,struct adios_index_var_struct_v1 ** vars_root
                       )
{
    struct adios_group_struct * g = fd->group;
    struct adios_var_struct * v = g->vars;
    struct adios_index_process_group_struct_v1 * g_item;

    uint64_t process_group_count = 0;
    uint16_t var_count = 0;

    g_item = (struct adios_index_process_group_struct_v1 *)
                malloc (sizeof (struct adios_index_process_group_struct_v1));
    g_item->group_name = g->name;
    g_item->process_id = g->process_id;
    g_item->timestep = 0;
    g_item->offset_in_file = fd->pg_start_in_file;
    g_item->next = 0;

    // build the groups and vars index
    index_append_process_group_v1 (pg_root, g_item);

    while (v)
    {
        if (v->write_offset != 0)
        {
            struct adios_index_var_struct_v1 * v_index;
            v_index = malloc (sizeof (struct adios_index_var_struct_v1));
            v_index->entries = malloc (
                                 sizeof (struct adios_index_var_entry_struct_v1)
                                 );

            v_index->group_name = g->name;
            v_index->var_name = v->name;
            v_index->var_path = v->path;
            v_index->type = v->type;
            v_index->entries_count = 1;
            v_index->entries_allocated = 1;
            v_index->entries [0].offset = v->write_offset;
            v_index->entries [0].min = v->min;
            v_index->entries [0].max = v->max;

            v_index->next = 0;

            // this fn will either take ownership for free
            index_append_var_v1 (vars_root, v_index);
        }

        v = v->next;
    }
}

uint64_t adios_calc_size_process_group_index_v1 (
                        struct adios_index_process_group_struct_v1 * pg_root
                       )
{
    uint64_t size = 2500;
fprintf (stderr, "How to do adios_calc_size_process_group_index?\n");

    return size;
}

int adios_write_index_v1 (char ** buffer
                         ,uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,uint64_t index_start
                         ,struct adios_index_process_group_struct_v1 * pg_root
                         ,struct adios_index_var_struct_v1 * vars_root
                         )
{
    uint64_t groups_count = 0;
    uint16_t vars_count = 0;

    uint64_t index_size = 0;
    uint64_t pg_index_start = index_start;
    uint64_t vars_index_start = 0;

    // we need to save the offset we will write the count and size
    uint64_t buffer_offset_start = 0; // since we realloc, we can't save a ptr

    // save for the process group index
    buffer_offset_start = *buffer_offset;

    *buffer_offset += (8 + 8); // save space for groups count and index size

    while (pg_root)
    {
        uint16_t len;
        uint16_t group_size = 0;
        uint64_t group_start = *buffer_offset;

        groups_count++;

        *buffer_offset += 2; // save space for the size

        len = strlen (pg_root->group_name);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        group_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset, pg_root->group_name, len);
        index_size += len;
        group_size += len;

        buffer_write (buffer, buffer_size, buffer_offset, &pg_root->process_id, 4);
        index_size += 4;
        group_size += 4;
        buffer_write (buffer, buffer_size, buffer_offset, &pg_root->timestep, 4);
        index_size += 4;
        group_size += 4;
        buffer_write (buffer, buffer_size, buffer_offset, &pg_root->offset_in_file, 8);
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

        vars_count++;

        *buffer_offset += 4; // save space for var length

        len = strlen (vars_root->group_name);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        var_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset, vars_root->group_name, len);
        index_size += len;
        var_size += len;

        len = strlen (vars_root->var_name);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        var_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset, vars_root->var_name, len);
        index_size += len;
        var_size += len;

        len = strlen (vars_root->var_path);
        buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
        index_size += 2;
        var_size += 2;
        buffer_write (buffer, buffer_size, buffer_offset, vars_root->var_path, len);
        index_size += len;
        var_size += len;

        flag = vars_root->type;
        buffer_write (buffer, buffer_size, buffer_offset, &flag, 1);
        index_size += 1;
        var_size += 1;

        buffer_write (buffer, buffer_size, buffer_offset, &vars_root->entries_count, 8);
        index_size += 8;
        var_size += 8;

        for (int i = 0; i < vars_root->entries_count; i++)
        {
            uint64_t size;

            buffer_write (buffer, buffer_size, buffer_offset, &vars_root->entries [i].offset, 8);
            index_size += 8;
            var_size += 8;

            size = adios_get_type_size (vars_root->type
                                       ,&vars_root->entries [i].min
                                       );

            buffer_write (buffer, buffer_size, buffer_offset, &vars_root->entries [i].min, size);
            index_size += size;
            var_size += size;

            buffer_write (buffer, buffer_size, buffer_offset, &vars_root->entries [i].max, size);
            index_size += size;
            var_size += size;
        }

        buffer_write (buffer, buffer_size, &var_start, &var_size, 4);

        vars_root = vars_root->next;
    }

    // vars index count/size prefix
    buffer_write (buffer, buffer_size, &buffer_offset_start, &vars_count, 2);
    buffer_write (buffer, buffer_size, &buffer_offset_start, &index_size, 8);

    // location of the beginning of the indexes (first proc groups then vars)
    buffer_write (buffer, buffer_size, buffer_offset, &pg_index_start, 8);
    buffer_write (buffer, buffer_size, buffer_offset, &vars_index_start, 8);

    return 0;
}

int adios_write_version_v1 (char ** buffer
                           ,uint64_t * buffer_size
                           ,uint64_t * buffer_offset
                           )
{
    uint64_t test = 1;

    if (!*(char *) &test)
        test = 0x80000000;
    else
        test = 0;

    test += 1;   // current version

    test = htonl (test);

    buffer_write (buffer, buffer_size, buffer_offset, &test, 4);

    return 0;
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

static uint16_t calc_dimension_size (struct adios_dimension_struct * dimension)
{
    uint16_t size = 0;

    size += 1; // var (y or n)

    if (dimension->dimension.id == 0)  // it is a number
    {
        size += 8;  // size of value
    }
    else   // it is a var
    {
        size += 2;  // size of var ID
    }

    if (dimension->global_dimension.id == 0)
    {
        size += 8; // default to a rank
    }
    else
    {
        if (dimension->global_dimension.id == 0)  // it is a number
        {
            size += 8;  // size of value
        }
        else   // it is a var
        {
            size += 2;  // size of var ID
        }
    }

    if (dimension->local_offset.id == 0)
    {
        size += 8;  // default to a rank
    }
    else
    {
        if (dimension->local_offset.id == 0)  // it is a number
        {
            size += 8;  // size of value
        }
        else   // it is a var
        {
            size += 2;  // size of var ID
        }
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

    if (dimension->dimension.id == 0)
    {
        var = 'n';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dimension->dimension.rank, 8);
        size += 8;
    }
    else
    {
        var = 'y';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dimension->dimension.id, 2);
        size += 2;
    }

    if (dimension->global_dimension.id == 0)
    {
        var = 'n';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dimension->global_dimension.rank, 8);
        size += 8;
    }
    else
    {
        var = 'y';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dimension->global_dimension.id, 2);
        size += 2;
    }

    if (dimension->local_offset.id == 0)
    {
        var = 'n';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dimension->local_offset.rank, 8);
        size += 8;
    }
    else
    {
        var = 'y';
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var, 1);
        size += 1;
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &dimension->local_offset.id, 2);
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

    total_size += adios_get_type_size (v->type, v->data);     // min
    total_size += adios_get_type_size (v->type, v->data);     // max
    total_size += adios_get_var_size (v, fd->group, v->data); // payload

    buffer_write (&fd->buffer, &fd->buffer_size, &start, &total_size, 8);

    fd->vars_written++;

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return total_size;
}

static void calc_min_max (struct adios_group_struct * group
                         ,struct adios_var_struct * var
                         )
{
    uint64_t total_size = adios_get_var_size (var, group, var->data);
    uint64_t size = 0;

    switch (var->type)
    {
        case adios_byte:
        {
            int8_t * data = (int8_t *) var->data;
            int8_t min = data [0];
            int8_t max = data [0];
            size++;
            while ((size * 1) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;
          
            return;
        }

        case adios_unsigned_byte:
        {
            uint8_t * data = (uint8_t *) var->data;
            uint8_t min = data [0];
            uint8_t max = data [0];
            size++;
            while ((size * 1) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;
          
            return;
        }

        case adios_string:
        {
            if (!var)
            {
                uint8_t data = adios_string;
                var->min = (void *) data;
                var->max = (void *) data;
                return;
            }
            else
            {
                uint8_t data = adios_string;
                var->min = (void *) data;
                var->max = (void *) data;
                return; // strlen ((char *) var);
            }
        }

        case adios_short:
        {
            int16_t * data = (int16_t *) var->data;
            int16_t min = data [0];
            int16_t max = data [0];
            size++;
            while ((size * 2) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;
          
            return;
        }

        case adios_unsigned_short:
        {
            uint16_t * data = (uint16_t *) var->data;
            uint16_t min = data [0];
            uint16_t max = data [0];
            size++;
            while ((size * 2) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;
          
            return;
        }

        case adios_integer:
        {
            int32_t * data = (int32_t *) var->data;
            int32_t min = data [0];
            int32_t max = data [0];
            size++;
            while ((size * 4) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;

            return;
        }

        case adios_unsigned_integer:
        {
            uint32_t * data = (uint32_t *) var->data;
            uint32_t min = data [0];
            uint32_t max = data [0];
            size++;
            while ((size * 4) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;

            return;
        }

        case adios_long:
        {
            int64_t * data = (int64_t *) var->data;
            int64_t min = data [0];
            int64_t max = data [0];
            size++;
            while ((size * 8) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;

            return;
        }

        case adios_unsigned_long:
        {
            uint64_t * data = (uint64_t *) var->data;
            uint64_t min = data [0];
            uint64_t max = data [0];
            size++;
            while ((size * 8) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;

            return;
        }

        case adios_real:
        {
            float * data = (float *) var->data;
            float min = data [0];
            float max = data [0];
            size++;
            while ((size * 4) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) *(uint32_t *) &min;
            var->max = (void *) *(uint32_t *) &max;

            return;
        }

        case adios_double:
        {
            double * data = (double *) var->data;
            double min = data [0];
            double max = data [0];
            size++;
            while ((size * 8) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) *(uint64_t *) &min;
            var->max = (void *) *(uint64_t *) &max;

            return;
        }

        case adios_long_double:
        {
            long double * data = (long double *) var->data;
            long double min = data [0];
            long double max = data [0];
            size++;
            while ((size * 16) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) *(uint64_t *) &min;
            var->max = (void *) *(uint64_t *) &max;

            return;
	}

        case adios_complex:
	{
            uint32_t * data = (uint32_t *) var->data;
            uint32_t min = (void *) data [0];
            uint32_t max = (void *) data [0];
            size++;
            while ((size * 4) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;

            return;
	}

        case adios_double_complex:
	{
            uint32_t * data = (uint32_t *) var->data;
            uint32_t min = (void *) data [0];
            uint32_t max = (void *) data [0];
            size++;
            while ((size * 4) < total_size)
            {
                if (data [size] < min)
                    min = data [size];
                if (data [size] > max)
                    max = data [size];
                size++;
            }
            var->min = (void *) min;
            var->max = (void *) max;

            return;
	}

        default:
	{
            uint64_t data = adios_unknown;
            var->min = (void *) *(uint64_t *) &data;
            var->max = (void *) *(uint64_t *) &data;
            return;
	}
    }
}

int adios_generate_var_characteristics_v1 (struct adios_file_struct * fd
                                          ,struct adios_var_struct * var
                                          )
{
    calc_min_max (fd->group, var);

    return 0;
}

int adios_write_var_characteristics_v1 (struct adios_file_struct * fd
                                       ,struct adios_var_struct * var
                                       )
{
    uint64_t size;

    size = adios_get_type_size (var->type, var->data);

    // write min
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var->min, size);

    // write max
    buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &var->max, size);

    if (fd->bytes_written < fd->offset)
        fd->bytes_written = fd->offset;

    return 0;
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
        buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &a->var->id, 2);
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

        if (a->type == adios_string)
            buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, a->value, t);
        else
            buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &a->value, t);
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

uint64_t adios_get_type_size (enum ADIOS_DATATYPES type, void * var)
{
    switch (type)
    {
        case adios_byte:
        case adios_unsigned_byte:
            return 1;

        case adios_string:
            if (!var)
                return 0;
            else
                return strlen ((char *) var);

        case adios_short:
        case adios_unsigned_short:
            return 2;

        case adios_integer:
        case adios_unsigned_integer:
            return 4;

        case adios_long:
        case adios_unsigned_long:
            return 8;

        case adios_real:
            return 4;

        case adios_double:
            return 8;

        case adios_long_double:
            return 16;

        case adios_complex:
            return 2 * 4;

        case adios_double_complex:
            return 2 * 8;

        default:
            return -1;
    }
}

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
                    ,adios_type_to_string (type)
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
                struct adios_var_struct * var = 0;

                var = adios_find_var_by_id (group->vars, d->dimension.id);

                // first check to make sure all vars are provided
                if (!var)
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
                                                           ,&attr->value
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
                    if (!var->data)
                    {
                        fprintf (stderr, "adios_get_var_size: "
                                         "sizing of %s failed because "
                                         "dimension component %s was not "
                                         "provided\n"
                                ,var->name, var->name
                                );

                        return 0;
                    }
                    else
                    {
                        if (!adios_multiply_dimensions (&size, var
                                                       ,var->type
                                                       ,var->data
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
                size *= d->dimension.rank;
            }

            d = d->next;
        }
    }

    return size;
}

const char * adios_type_to_string (int type)
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
