#include <math.h>
#include <string.h>
#include <unistd.h>

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

// Chen's encoder
#include "bw-utils.h"
#include "br-utils.h"

#include "binpack-general.h"
#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"

enum ADIOS_FLAG adios_host_language_fortran = adios_flag_yes;

MPI_Comm adios_mpi_comm_world;
MPI_Comm adios_mpi_comm_self;
MPI_Info adios_mpi_info;

// buffer sizing may be problematic.  To get a more accurate picture, check:
// http://chandrashekar.info/vault/linux-system-programs.html
static unsigned long long adios_buffer_size_requested = 0;
static unsigned long long adios_buffer_size_max = 0;
static unsigned long long adios_buffer_size_remaining = 0;
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

static int parseType (const char * type, const char * name)
{
    if (!strcmp (type, "byte"))
        return adios_byte;

    if (!strcmp (type, "integer*4") || !strcmp (type, "integer"))
        return adios_integer;

    if (!strcmp (type, "integer*8") || !strcmp (type, "long"))
        return adios_long;

    if (!strcmp (type, "real*4") || !strcmp (type, "real"))
        return adios_real;

    if (!strcmp (type, "complex"))
        return adios_complex;

    if (!strcmp (type, "string"))
        return adios_string;

    if (!strcmp (type, "real*8") || !strcmp (type, "double"))
        return adios_double;

    fprintf (stderr, "config.xml: invalid type: %s in var %s\n", type, name);

    return adios_unknown;
}

static int parseFlag (const char * attr_name, const char * flag
                     ,int default_value
                     )
{
    if (!flag)
        return default_value;

    if (!strcmp (flag, "yes"))
        return adios_flag_yes;

    if (!strcmp (flag, "no"))
        return adios_flag_no;

    fprintf (stderr, "config.xml: %s must have a value of 'yes' or 'no' "
                     "not: %s\n", attr_name, flag
            );

    return adios_flag_unknown;
}

static int is_var (char * temp) // 1 == yes, 0 == no
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
                item->item.var = adios_find_var_by_name (new_group->vars, tmp);
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
            item->item.var = adios_find_var_by_name (new_group->vars, c);
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
                item->item.var = adios_find_var_by_name (new_group->vars, tmp);
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
                item->item.var = adios_find_var_by_name (new_group->vars, tmp);
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
                var->var = adios_find_var_by_name (new_group->vars, tmp);
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
                var->var = adios_find_var_by_name (new_group->vars, tmp);
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
            item->var = adios_find_var_by_name (new_group->vars, c);
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
            item->item.var = adios_find_var_by_name (new_group->vars, c);
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
            var->var = adios_find_var_by_name (new_group->vars, c);
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
            var->var = adios_find_var_by_name (new_group->vars, c);
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
            var = adios_find_var_by_name (new_group->vars, c);
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
            cell_list->cell_list.count.var = adios_find_var_by_name (new_group->vars, c);
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
            cell_list->cell_list.data = adios_find_var_by_name (new_group->vars, c);
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
            cell_list->cell_list.type.var = adios_find_var_by_name (new_group->vars, c);
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
            cell_list->cell_list.count.var = adios_find_var_by_name (new_group->vars, c);
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
            cell_list->cell_list.data = adios_find_var_by_name (new_group->vars, c);
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
            cell_list->cell_list.type.var = adios_find_var_by_name (new_group->vars, c);
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

        if (!strcmp (n->value.element.name, "dimensions"))
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
        if (!strcmp (n->value.element.name, "origin"))
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
        if (!strcmp (n->value.element.name, "spacing"))
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

        if (!strcmp (n->value.element.name, "dimensions"))
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
        if (!strcmp (n->value.element.name, "coordinates-multi-var"))
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
        if (!strcmp (n->value.element.name, "coordinates-single-var"))
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

        if (!strcmp (n->value.element.name, "nspace"))
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
        if (!strcmp (n->value.element.name, "dimensions"))
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
        if (!strcmp (n->value.element.name, "points-multi-var"))
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
        if (!strcmp (n->value.element.name, "points-single-var"))
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

        if (!strcmp (n->value.element.name, "points"))
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
        if (!strcmp (n->value.element.name, "uniform-cells"))
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
        if (!strcmp (n->value.element.name, "mixed-cells"))
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

static int parseGroup (mxml_node_t * node)
{
    mxml_node_t * n;
    const char * datagroup_name;
    const char * coordination_comm;
    const char * coordination_var;
    struct adios_group_struct * new_group;

    datagroup_name = mxmlElementGetAttr (node, "name");
    coordination_comm = mxmlElementGetAttr (node, "coordination-communicator");
    coordination_var = mxmlElementGetAttr (node, "coordination-var");
    if (!datagroup_name)
    {
        fprintf (stderr,
                 "config.xml: name attribute required on adios-group\n");

        return 0;
    }
    adios_common_declare_group ((long long *) &new_group, datagroup_name
                               ,adios_host_language_fortran
                               ,coordination_comm, coordination_var);

    for (n = mxmlWalkNext (node, node, MXML_DESCEND)
        ;n
        ;n = mxmlWalkNext (n, node, MXML_NO_DESCEND)
        )
    {
        if (n->type != MXML_ELEMENT)
        {
            continue;
        }

        if (!strcmp (n->value.element.name, "var"))
        {
            const char * name;
            const char * path;
            const char * type;
            const char * dimensions;
            const char * copy_on_write;
            int t1;
            int c1;

            name = mxmlElementGetAttr (n, "name");
            path = mxmlElementGetAttr (n, "path");
            type = mxmlElementGetAttr (n, "type");
            dimensions = mxmlElementGetAttr (n, "dimensions");
            copy_on_write = mxmlElementGetAttr (n, "copy-on-write");
            if (!name)
                name = "";  // this will catch the error
            if (!path)
                path = "/";
            if (!type)
                type = ""; // this will catch the error
            t1 = parseType (type, name);
            c1 = parseFlag ("copy-on-write", copy_on_write, adios_flag_no);

            if (!dimensions)
            {
                dimensions = mxmlElementGetAttr (n, "dimension");
                if (!dimensions)
                    dimensions = "";
            }

            if (!adios_common_define_var (*(long long *) &new_group, name
                                         ,path, t1, c1, dimensions, 0
                                         )
               )
            {
                return 0;
            }
        } else
        if (!strcmp (n->value.element.name, "global-bounds"))
        {
            mxml_node_t * n1;   // used for global_bounds
            struct adios_global_bounds_struct * new_global_bounds = 0;

            const char * dimensions;
            const char * offsets;

            dimensions = mxmlElementGetAttr (n, "dimensions");
            offsets = mxmlElementGetAttr (n, "offsets");

            if (!dimensions)
            {
                dimensions = mxmlElementGetAttr (n, "global-dimensions");
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
                fprintf (stderr, "config.xml: offsets required on "
                                 "global-bounds\n"
                        );

                return 0;
            }

            if (!adios_common_define_global_bounds (*(long long *) &new_group
                                                   ,dimensions, offsets
                                                   ,&new_global_bounds
                                                   )
               )
            {
                return 0;
            }

            for (n1 = mxmlWalkNext (n, n, MXML_DESCEND)
                ;n1
                ;n1 = mxmlWalkNext (n1, n, MXML_DESCEND)
                )
            {
                if (n1->type != MXML_ELEMENT)
                {
                    continue;
                }

                if (!strcmp (n1->value.element.name, "var"))
                {
                    const char * name;
                    const char * path;
                    const char * type;
                    const char * dimensions;
                    const char * copy_on_write;
                    int t1;
                    int c1;

                    name = mxmlElementGetAttr (n1, "name");
                    path = mxmlElementGetAttr (n1, "path");
                    type = mxmlElementGetAttr (n1, "type");
                    dimensions = mxmlElementGetAttr (n1, "dimensions");
                    copy_on_write = mxmlElementGetAttr (n1, "copy-on-write");
                    if (!name)
                        name = "";  // this will catch the error
                    if (!path)
                        path = "/";
                    if (!type)
                        type = ""; // this will catch the error
                    t1 = parseType (type, name);
                    c1 = parseFlag ("copy-on-write", copy_on_write, adios_flag_no);
                    if (!dimensions)
                        dimensions = mxmlElementGetAttr (n, "dimension");


                    if (!adios_common_define_var (*(long long *) &new_group, name
                                                 ,path, t1, c1, dimensions
                                                 ,new_global_bounds
                                                 )
                       )
                    {
                        return 0;
                    }
                } else
                {
                    if (!strncmp (n1->value.element.name, "!--", 3)) // a comment
                    {
                        continue;
                    }
                    else
                    {
                        fprintf (stderr, "config.xml: invalid xml element: '%s'\n"
                                ,n1->value.element.name
                                );

                        return 0;
                    }
                }
            }
        } else
        if (!strcmp (n->value.element.name, "attribute"))
        {
            const char * name;
            const char * path;
            const char * value;
            const char * type;
            const char * var;

            name = mxmlElementGetAttr (n, "name");
            path = mxmlElementGetAttr (n, "path");
            value = mxmlElementGetAttr (n, "value");
            type = mxmlElementGetAttr (n, "type");
            var = mxmlElementGetAttr (n, "var");

            if (!name)
            {
                fprintf (stderr, "config.xml: attribute element requires name\n");

                return 0;
            }
            if (!path)
            {
                fprintf (stderr, "config.xml: attribute element requires path\n");

                return 0;
            }
            if (value || type || var)
            {
                if (!(   (!value && type && var)
                      || (value && !type && !var)
                     )
                   )
                {
                    fprintf (stderr, "config.xml: attriute element '%s' "
                                     "requires either value OR type and var\n"
                            ,name
                            );

                    return 0;
                }
            }
            else
            {
                fprintf (stderr, "config.xml: attriute element '%s' "
                                 "requires either value OR type and var\n"
                        ,name
                        );

                return 0;
            }
          
#if 0
            if (!value)
            {
                fprintf (stderr, "config.xml: attribute element requires value\n");

                return 0;
            }
#endif
            if (!adios_common_define_attribute (*(long long *) &new_group, name
                                               ,path, value, type, var
                                               )
               )
            {
                return 0;
            }
        } else
        if (!strcmp (n->value.element.name, "mesh"))
        {
            const char * type;
            const char * time_varying;

            new_group->mesh = (struct adios_mesh_struct *) malloc
                                           (sizeof (struct adios_mesh_struct));
            type = mxmlElementGetAttr (n, "type");
            time_varying = mxmlElementGetAttr (n, "time-varying");

            if (!type)
                type = "";

            if (!strcmp (type, "uniform"))
            {
                new_group->mesh->type = ADIOS_MESH_UNIFORM;
                new_group->mesh->uniform =
                    (struct adios_mesh_uniform_struct *)
                         calloc (1, sizeof (struct adios_mesh_uniform_struct));
                parseMeshUniform (n, new_group, &new_group->mesh->uniform);
            } else
            if (!strcmp (type, "structured"))
            {
                new_group->mesh->type = ADIOS_MESH_STRUCTURED;
                new_group->mesh->structured =
                    (struct adios_mesh_structured_struct *)
                         calloc (1, sizeof (struct adios_mesh_structured_struct));
                parseMeshStructured (n, new_group, &new_group->mesh->structured);
            } else
            if (!strcmp (type, "rectilinear"))
            {
                new_group->mesh->type = ADIOS_MESH_RECTILINEAR;
                new_group->mesh->rectilinear =
                    (struct adios_mesh_rectilinear_struct *)
                         calloc (1, sizeof (struct adios_mesh_rectilinear_struct));
                parseMeshRectilinear (n, new_group, &new_group->mesh->rectilinear);
            } else
            if (!strcmp (type, "unstructured"))
            {
                new_group->mesh->type = ADIOS_MESH_UNSTRUCTURED;
                new_group->mesh->unstructured =
                    (struct adios_mesh_unstructured_struct *)
                         calloc (1, sizeof (struct adios_mesh_unstructured_struct));

                parseMeshUnstructured (n, new_group, &new_group->mesh->unstructured);
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

    const char * priority;
    const char * method;
    const char * iterations;
    const char * group;
    const char * parameters;
    const char * base_path;
    int p1;
    int i1;

    priority = mxmlElementGetAttr (node, "priority");
    iterations = mxmlElementGetAttr (node, "iterations");
    base_path = mxmlElementGetAttr (node, "base-path");
    method = mxmlElementGetAttr (node, "method");
    group = mxmlElementGetAttr (node, "group");
    /* Check for parameters, if they exist */
    n = mxmlWalkNext (node, node, MXML_DESCEND);
    if ( n != NULL) 
      {parameters = n->value.text.string;    }
    else
      {parameters = NULL;}
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

    if (!adios_common_select_method (p1, method, parameters, group, base_path, i1))
    {
        return 0;
    }

    return 1;
}

static int parseBuffer (mxml_node_t * node)
{
    const char * size_MB;
    const char * free_memory_percentage;
    const char * allocate_time;

    int size = -1;

    size_MB = mxmlElementGetAttr (node, "size-MB");
    free_memory_percentage = mxmlElementGetAttr (node, "free-memory-percentage");
    allocate_time = mxmlElementGetAttr (node, "allocate-time");

    if ((!size_MB && !free_memory_percentage) || !allocate_time)
    {
        fprintf (stderr, "config.xml: must define allocate-time and either "
                         "size-MB or free-memory-percentage for buffer element\n"
                );

        return 0;
    }
    else
    {
        if (!strcmp (allocate_time, "now"))
        {
            adios_buffer_alloc_when = ADIOS_BUFFER_ALLOC_NOW;
        }
        else
        {
            if (!strcmp (allocate_time, "oncall"))
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

            adios_buffer_size_requested = (  (unsigned long long) size
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
                adios_buffer_size_requested = (unsigned long long) size;
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
                        ,((unsigned long long) (pagesize * pages))
                        );
                 adios_buffer_size_max = (unsigned long long) (pagesize * pages);
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

unsigned long long adios_method_buffer_alloc (unsigned long long size)
{
    if (adios_buffer_size_remaining >= size)
    {
        adios_buffer_size_remaining -= size;

        return size;
    }
    else
    {
        unsigned long long remaining = adios_buffer_size_remaining;

        adios_buffer_size_remaining = 0;

        return remaining;
    }
}

int adios_method_buffer_free (unsigned long long size)
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
                                                 )
{
    int done = 0;

    if (!name)
    {
        done = 1;
        root = 0;
    }

    while (!done && root)
    {
        if (!strcmp (name, root->name))
        {
            done = 1;
        }
        else
        {
            root = root->next;
        }
    }

    return root;
}

struct adios_var_struct * adios_find_attribute_var_by_name (struct adios_attribute_struct * root, const char * name)
{
    int done = 0;
    struct adios_var_struct * v = 0;

    if (!name)
    {
        done = 1;
    }

    while (!done && root)
    {
        if (root->var.name && !strcmp (name, root->var.name))
        {
            done = 1;
            v = &root->var;
        }
        else
        {
            root = root->next;
        }
    }

    return v;
}

#if 0
static int parse_subitem (char * d, struct adios_group_struct * g
                         ,struct adios_dimension_item_struct * bound
                         ,struct adios_dimension_item_struct * use_bound
                         )
{
/*  These have been tried successfully
    char * tests [] =
    {    "5"
        ,"x"
        ,"x(y)"
        ,"1:10"
        ,"-3:5"
        ,"x:10"
        ,"10:x"
        ,"1(5):10"
        ,"1(x):10"
        ,"x(y):10"
        ,"2:10(5)"
        ,"3:10(x)"
        ,"4:x(y)"
        ,"x(y):z(w)"
        ,0
    };
*/
    char * p = d;  // paren location
    int has_paren = 0;

    while (p && *p && !has_paren)
    {
        if (*p == '(')
            has_paren = 1;
        else
            p++;
    }

    if (has_paren)
    {
        char * left = d;
        char * right = p + 1;
        *p = '\0';
        // kill the close paren
        while (*right)
        {
            if (*right == ')')
                *right = '\0';
            else
                right++;
        }
        // put right back at the front of that piece
        right = p + 1;

        if (is_var (left))
        {
            bound->rank = 0;
            bound->var = adios_find_var_by_name (g->vars, left);
            if (!bound->var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,left
                        );
    
                return 0;
            }
        }
        else
        {
            bound->var = 0;
            bound->rank = atoi (left);
        }

        if (is_var (right))
        {
            use_bound->rank = 0;
            use_bound->var = adios_find_var_by_name (g->vars, right);
            if (!use_bound->var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,right
                        );

                return 0;
            }
        }
        else
        {
            use_bound->rank = atoi (right);
            use_bound->var = 0;
        }
    }
    else
    {
        if (is_var (d))
        {
            bound->rank = 0;
            use_bound->rank = 0;
            bound->var = adios_find_var_by_name (g->vars, d);
            if (!bound->var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,d
                        );

                return 0;
            }
            use_bound->var = adios_find_var_by_name (g->vars, d);
            if (!use_bound->var)
            {
                fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                        ,d
                        );

                return 0;
            }
        }
        else
        {
            bound->rank = atoi (d);
            use_bound->rank = bound->rank;
            bound->var = 0;
            use_bound->var = 0;
        }
    }

    return 1;
}
#endif

void adios_parse_dimension (char * dimension, struct adios_group_struct * g
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
        dim->dimension.rank = 0;
        dim->dimension.var = adios_find_var_by_name (g->vars, dimension);
        if (!dim->dimension.var)
        {
            fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                    ,dimension
                    );

            return;
        }
    }
    else
    {
        dim->dimension.var = 0;
        dim->dimension.rank = atoi (dimension);
    }
}

#if 0
// format:
// x(u):y(v):z
// x = allocated lower bound
// y = allocated upper bound
// u = use lower bound
// v = use upper bound
// z = stride
// if stride is to be provided, x is required.
// if stride is to be defaulted to 1, then you should omit the : too.
// u & v are optional and default to x & y respectively.
// x is optional and defaults to 1 for Fortran and 0 for C
//    (exception to optionality: using stride)
void adios_parse_dimension (char * dimension, struct adios_group_struct * g
                           ,struct adios_dimension_struct * dim
                           )
{
    char * d = strdup (dimension);  // since we modify the string
    char * c_temp = d; // use to find the colon locations
    char * c1 = NULL; // first colon location
    char * c2 = NULL; // second colon location
    char ** temp = NULL; // colon location selector

    // find the location of the first and second colon
    while (c_temp && *c_temp && !c2)
    {
        if (!c1)  // not seen first :
        {
            temp = &c1;
        }
        else
        {
            temp = &c2;
        }
 
        if (*c_temp == ':')
            *temp = c_temp;

        c_temp++;
    }

    if (!c1)  // no colons
    {
        dim->lower_bound.rank = 1;
        dim->lower_bound.var = 0;
        dim->use_lower_bound = dim->lower_bound;
        dim->stride.rank = 1;
        dim->stride.var = 0;

        parse_subitem (d, g, &dim->upper_bound, &dim->use_upper_bound);
        if (adios_host_language_fortran == adios_flag_no)
        {
            dim->lower_bound.rank--;
            dim->use_lower_bound.rank--;
            dim->upper_bound.rank--;
            dim->use_upper_bound.rank--;
        }
        //printf ("dimension: %s\n", d);
    }
    else
    {
        char * lower;
        char * upper;
        char * stride;

        if (c2)
        {
            *c2 = '\0';
            c2++;
        }
        *c1 = '\0';
        c1++;

        lower = d;
        upper = c1;
        stride = c2;

        parse_subitem (lower, g, &dim->lower_bound, &dim->use_lower_bound);
        parse_subitem (upper, g, &dim->upper_bound, &dim->use_upper_bound);
        if (stride)
        {
            if (is_var (stride))
            {
                dim->stride.rank = 0;
                dim->stride.var = adios_find_var_by_name (g->vars, stride);
                if (!dim->stride.var)
                {
                    fprintf (stderr, "config.xml: invalid var dimension: %s\n"
                            ,stride
                            );
    
                    return;
                }
            }
            else
            {
                dim->stride.rank = atoi (stride);
                dim->stride.var = 0;
            }
        }
        else
        {
            dim->stride.rank = 1;
            dim->stride.var = 0;
        }
        //printf ("\tleft: %s right: %s stride: %s\n", left, right, stride);
    }

    free (d);
}
#endif

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
        fprintf (stderr, "config.xml: unknown error parsing XML\n"
                         "Did you remember to start the file with\n"
                         "<?xml version=\"1.0\"?>");

        return 0;
    }

    root = mxmlWalkNext (doc, doc, MXML_DESCEND); // get rid of the xml version
    root = mxmlWalkNext (root, doc, MXML_DESCEND); // should be our root tag
    
    while (!strncmp (root->value.element.name, "!--", 3)) // skip comments
    {
        root = mxmlWalkNext (root, doc, MXML_NO_DESCEND);
    }

    if (strcmp (root->value.element.name, "adios-config"))
    {
        fprintf (stderr, "config.xml: invalid root xml element: %s\n"
                ,root->value.element.name
                );

        mxmlRelease (doc);

        return 0;
    }
    else
    {
        const char * host_language = NULL;
        host_language = mxmlElementGetAttr (root, "host-language");
        if (!host_language)
        {
            host_language = "Fortran";
        }

        if (!strcmp (host_language, "Fortran"))
        {
            adios_host_language_fortran = adios_flag_yes;
        }
        else
        {
            if (!strcmp (host_language, "C"))
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

        if (!strcmp (node->value.element.name, "adios-group"))
        {
            if (!parseGroup (node))
                break;
            saw_datagroup = 1;
        }
        else
        {
            if (!strcmp (node->value.element.name, "method"))
            {
                if (!parseMethod (node))
                    break;
                saw_method = 1;
            }
            else
            {
                if (!strcmp (node->value.element.name, "buffer"))
                {
                    if (!parseBuffer (node))
                        break;
                    saw_buffer = 1;
                }
                else
                {
                    if (!strncmp (node->value.element.name, "!--", 3)) // a comment
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

int adios_common_define_attribute (long long group, const char * name
                                  ,const char * path, const char * value
                                  ,const char * type, const char * var
                                  )
{
    struct adios_group_struct * g = (struct adios_group_struct *) group;
    struct adios_attribute_struct * attr = (struct adios_attribute_struct *)
                              malloc (sizeof (struct adios_attribute_struct));

#if 0
    struct adios_var_struct * f;

    if (attr_type == adios_attribute_data)
    {
        char * a = strdup (path);
        char * extracted_name;  // for testing
        char * extracted_path;  // for testing
        int last_slash = -1;
        int len;
        int i;

        // extract name and path from path
        extracted_path = a;
        len = strlen (extracted_path);
        if (extracted_path [len - 1] == '/')
            extracted_path [len - 1] = 0;
        len--;
        for (i = len - 1; i >= 0; i--)
        {
            if (extracted_path [i] == '/')
            {
                last_slash = i;
                break;
            }
        }

        if (last_slash >= 0)
        {
            extracted_name = strdup (extracted_path + last_slash + 1);
            if (last_slash == 0)
                extracted_path [last_slash + 1] = 0;
            else
                extracted_path [last_slash] = 0;
        }
        else
        {
            fprintf (stderr, "config.xml: data item '%s' for attribute '%s' "
                             "not found (no root '/' provided in path: "
                             "last_slash position: %d extracted_path: '%s')\n"
                    ,path, name, last_slash, extracted_path
                    );
            free (a);

            return 0;
        }

        f = adios_find_var_by_name (t->vars, extracted_name);
        if (f)
        {
            if (   strcmp (extracted_path, v->path)
                || strcmp (extracted_name, v->name)
               )
            {
                fprintf (stderr, "config.xml: data item '%s' for attribute '%s' "
                                 "not found (extracted_path: '%s' "
                                 "extracted_name: '%s')\n"
                        ,path, name, extracted_path, extracted_name
                        );
                free (a);
                free (extracted_name);

                return 0;
            }

//            free (a);
//            free (extracted_name);
        }
        else
        {
            fprintf (stderr, "config.xml: data item '%s' for attribute '%s' "
                             "not found\n", path, name
                    );
            free (a);
            free (extracted_name);

            return 0;
        }
    }
    else
    {
        if (attr_type == adios_attribute_group)
        {
            v = t->vars;
            int done = 0;
            int len = strlen (path);
            char * a = strdup (path);
            if (len > 1)
            {
                if (a [len - 1] == '/')
                    a [len - 1] = 0;
                len--;
            }

            while (v && !done)
            {
                // check for trailing '/' on both parts and fix up
                if (!strcmp (v->path, a))
                {
                    done = 1;
                }
                else
                {
                    v = v->next;
                }
            }
            if (!done)
            {
                fprintf (stderr, "config.xml: group '%s' for attribute '%s' "
                                 "not found\n", path, name
                        );

                return 0;
            }
        }
        else
        {
            fprintf (stderr, "invalid attribute type: %d\n", attr_type);

           return 0;
        }
    }
#endif

    attr->name = strdup (name);
    attr->path = strdup (path);

    if (value)
    {
        attr->var.type = adios_string;
        attr->var.data = adios_dupe_data (&attr->var, (void *) value);
    }
    else
    {
        attr->var.name = strdup (var);
        attr->var.data = 0;
        attr->var.type = parseType (type, name);
    }

    attr->next = 0;

    adios_append_attribute (&g->attributes, attr);

    return 1;
}

void * adios_dupe_data (struct adios_var_struct * v, void * data)
{
    struct adios_dimension_struct * d = v->dimensions;
    int element_size = bp_getsize (v->type, data);

    if (v->dimensions)
    {
        int rank = 10;
        int full_write = 1;
        unsigned long long use_count = 1;
        unsigned long long total_count = 1;
        struct adios_bp_dimension_struct * dims =
               (struct adios_bp_dimension_struct *)
                      calloc (rank, sizeof (struct adios_bp_dimension_struct));

        adios_dims_to_bp_dims (v->name, d, v->global_bounds, &rank, dims);
        adios_var_element_count (rank, dims, &use_count, &total_count);

        if (use_count == total_count)
            full_write = 1;
        else
            full_write = 0;

        d = malloc (use_count * element_size); // no strings
        if (!d)
        {
            fprintf (stderr, "cannot allocate %llu byte buffer to duplicate "
                             "array %s\n"
                    ,use_count * element_size
                    ,v->name
                    );
            return 0;
        }

        memcpy (d, data, use_count * element_size);
    }
    else
    {
        d = malloc (bp_getsize (v->type, data));

        if (!d)
        {
            fprintf (stderr, "cannot allocate %d bytes to copy scalar %s\n"
                    ,bp_getsize (v->type, data)
                    ,v->name
                    );

            return 0;
        }

        switch (v->type)
        {
            case adios_byte:
            case adios_integer:
            case adios_long:
            case adios_real:
            case adios_double:
            case adios_string:
                memcpy ((char *) d, (char *) data, element_size);
                break;

            default:
                d = 0;
                break;
        }
    }

    return d;
}

int adios_do_write_var (struct adios_var_struct * v
                       ,void * buf
                       ,unsigned long long buf_size
                       ,unsigned long long buf_start
                       ,unsigned long long * buf_end
                       )
{
    long long size = 0;
    int start = (int) buf_start;
    int end = (int) *buf_end;

    if (v->data)
    {
        if (!v->dimensions)
        {
            size = bcalsize_scalar (v->path, v->name, v->type
                                   ,v->data
                                   );
            if (size + buf_start > buf_size)
            {
                return 1; // overflowed;
            }
            //printf ("Scalar name: %s total size: %d\n", v->name, size);
            bw_scalar (buf, start, &end, v->path, v->name, v->data
                      ,v->type
                      );
        }
        else
        {
            int rank = 10;
            struct adios_bp_dimension_struct * dims = 0;

            dims = (struct adios_bp_dimension_struct *)
                    calloc (rank, sizeof (struct adios_bp_dimension_struct));

            adios_dims_to_bp_dims (v->name, v->dimensions
                                  ,v->global_bounds, &rank, dims
                                  );
            size = bcalsize_dset (v->path, v->name, v->type
                                 ,rank, dims
                                 );

            if (size + buf_start > buf_size)
            {
                return 1; // overflowed
            }
            //printf ("\ttotalsize: %d\n", size);
            bw_dset (buf, start, &end, v->path, v->name, v->data
                    ,v->type, rank, dims
                    );
            free (dims);
        }
        buf_start = start;
        *buf_end = end;
    }
    else
    {
        //printf ("Skipping %s (no data provided)\n"
        //       ,v->name
        //       );
    }

    return 0;
}

int adios_do_write_attribute (struct adios_attribute_struct * a
                             ,void * buf
                             ,unsigned long long buf_size
                             ,unsigned long long buf_start
                             ,unsigned long long * buf_end
                             )
{
    long long size = 0;
    int start = (int) buf_start;
    int end = (int) *buf_end;

    size = bcalsize_attr (a->path, a->name, a->var.data, a->var.type);
    if (size + buf_start > buf_size)
    {
        return 1; // overflowed
    }

    bw_attr (buf, start, &end, a->path, a->name, a->var.data, a->var.type);

    buf_start = start;
    *buf_end = end;

    return 0;
}

void adios_pre_element_fetch (struct adios_bp_element_struct * element
                             ,void ** buffer, long long * buffer_size
                             ,void * private_data
                             )
{
    struct adios_parse_buffer_struct * d =
                  (struct adios_parse_buffer_struct *) private_data;
    struct adios_var_struct * v;

    // if it isn't a data element, skip it
    if (element->tag == SCR_TAG || element->tag == DST_TAG)
    {
        v = adios_find_var_by_name (d->vars, element->name);
        if (!v)
        {
            fprintf (stderr, "Data item %s being read ignored\n", element->name);
        }
    }
    else
    {
        v = 0;
    }

    if (v)
    {
        if (v->data)
        {
            if (v->dimensions)
            {
                int full_read = 1;
                unsigned long long use_count = 1;
                unsigned long long total_count = 1;

                adios_var_element_count (element->ranks, element->dims
                                        ,&use_count, &total_count
                                        );

                if (use_count == total_count)
                    full_read = 1;
                else
                    full_read = 0;

                if (full_read)
                {
                    *buffer = v->data;
                    *buffer_size = element->size;
                }
                else
                {
                    if (d->buffer_len < element->size)
                    {
                        if (d->buffer)
                            free (d->buffer);
                        d->buffer = malloc (element->size);
                        d->buffer_len = element->size;
                    }
                    *buffer = d->buffer;
                    *buffer_size = element->size;
                }
            }
            else
            {
                *buffer = v->data;
                *buffer_size = element->size;
            }
        }
        else
        {
            *buffer = 0;
            *buffer_size = 0;
        }
    }
    else
    {
        *buffer = 0;
        *buffer_size = 0;
    }
}

void adios_post_element_fetch (struct adios_bp_element_struct * element
                              ,void * buffer, long long buffer_size
                              ,void * private_data
                              )
{
    struct adios_parse_buffer_struct * d =
                  (struct adios_parse_buffer_struct *) private_data;
    struct adios_var_struct * v;

    // if it isn't a data element, skip it
    if (element->tag == SCR_TAG || element->tag == DST_TAG)
        v = adios_find_var_by_name (d->vars, element->name);
    else
        v = 0;

    if (v)
    {
        if (v->dimensions)
        {
            int full_read = 1;
            int i = 0;
            unsigned long long use_count = 1;
            unsigned long long total_count = 1;

            adios_var_element_count (element->ranks, element->dims
                                    ,&use_count, &total_count
                                    );

            if (use_count == total_count)
                full_read = 1;
            else
                full_read = 0;

            if (full_read)
            {
                // already wrote into proper buffer
            }
            else // copy in into the proper places
            {
                int element_size = bp_getsize (element->type, v->data);
                int reading = 0;
                char * b = (char *) v->data;
                char * start = 0;
                int position [element->ranks];
                for (i = 0; i < total_count; i++)
                {
                    int use_it = adios_should_use_data (i, element->ranks
                                                       ,element->dims
                                                       ,position
                                                       );

                    if (use_it)
                    {
                        if (reading)
                        {
                            // do nothing
                        }
                        else
                        {
                            reading = 1;
                            start = b;
                        }
                    }
                    else
                    {
                        if (reading)
                        {
                            memcpy (start, buffer, b - start);
                            buffer += b - start;
                            reading = 0;
                            start = 0;
                        }
                        else
                        {
                            // do nothing
                        }
                    }
                    b += element_size;
                }

                if (start) // should be reading something through end
                {
                    memcpy (start, buffer, b - start);
                    buffer += b - start;
                }
            }
        }

        v->data = 0;
    }
    else
    {
        // didn't read anything so nothing to cleanup
    }
}

void adios_parse_buffer (struct adios_file_struct * fd, char * buffer
                        ,long long len
                        )
{
    struct adios_var_struct * v = fd->group->vars;

    unsigned long DATALEN = 300 * 1024;

    long long handle;
    struct adios_bp_element_struct * element;
    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer_len = DATALEN;
    data.buffer = malloc (DATALEN);
    if (!data.buffer)
    {
        fprintf (stderr, "cannot allocate %lu for data buffer\n", DATALEN);

        return;
    }

    handle = br_bopen (buffer, len);
    while (br_get_next_element_specific (handle
                                        ,adios_pre_element_fetch
                                        ,adios_post_element_fetch
                                        ,&data
                                        ,&element
                                        )
          )
    {
        br_free_element (element);
    }

    br_bclose (handle);

    if (data.buffer)
        free (data.buffer);
}

unsigned long long adios_data_size (struct adios_group_struct * g)
{
    unsigned long long size = 0;
    struct adios_var_struct * v;
    struct adios_attribute_struct * a;

    if (!g)
        return -1;

    v = g->vars;
    a = g->attributes;

    while (v)
    {
        if (v->data)
            size += adios_size_of_var (v, v->data);
        v = v->next;
    }

    while (a)
    {
        size += adios_size_of_attribute (a);
        a = a->next;
    }

    return size;
}

unsigned long long adios_size_of_var (struct adios_var_struct * v, void * data)
{
    unsigned long long size = 0;

    if (!v->dimensions)
    {
        size = bcalsize_scalar (v->path, v->name, v->type, data);
    }
    else
    {
        int rank = 10;
        struct adios_bp_dimension_struct * dims =
            (struct adios_bp_dimension_struct *)
               calloc (rank, sizeof (struct adios_bp_dimension_struct));

        adios_dims_to_bp_dims (v->name, v->dimensions, v->global_bounds, &rank, dims);
        size = bcalsize_dset (v->path, v->name, v->type, rank, dims);
        free (dims);
    }
   //printf("name:%s, size:%llu\n",v->name,size);

    return size;
}

unsigned long long adios_size_of_attribute (struct adios_attribute_struct * a)
{
    unsigned long long size = 0;

    size = bcalsize_attr (a->path, a->name, a->var.data, a->var.type);

    return size;
}

void adios_extract_string (char * out, const char * in, int size)
{
    if (in && out)
        strcpy (out, in);
/* for some Fortran implementations, we get a size for a string.
 * for others (like PGI), we don't and it isn't null terminated
 * unless we do it explicitly.  Assume that it is null terminated
 * for now.
 */
/*
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
*/
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

#if 0
void adios_append_global_bounds (struct adios_global_bounds_struct * bounds)
{
    struct adios_global_bounds_struct ** root = &adios_global_bounds;

    while (root)
    {
        if (!*root)
        {
            *root = bounds;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }
}
#endif

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

void adios_append_var (struct adios_var_struct ** root
                      ,struct adios_var_struct * var
                      )
{
    int id = 1;

    while (root)
    {
        if (!*root)
        {
            var->id = id;
            *root = var;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
            id++;
        }
    }
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
                            )
{
    int id = 1;

    while (root)
    {
        if (!*root)
        {
            attribute->var.id = id;
            *root = attribute;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
            id++;
        }
    }
}

void adios_dims_to_bp_dims (char * name
                           ,struct adios_dimension_struct * adios_dims
                           ,struct adios_global_bounds_struct * global
                           ,int * rank
                           ,struct adios_bp_dimension_struct * bp_dims
                           )
{
    int i = *rank;
    struct adios_dimension_struct * d = adios_dims;
    struct adios_dimension_struct * g = 0;
    struct adios_dimension_struct * o = 0;

    // check size of target bp_dimension array
    while (d)
    {
        i--;
        d = d->next;
    }
    if (i < 0)
    {
        fprintf (stderr, "conversion of '%s' failed because too small "
                         "destination array provided.\n"
                ,name
                );
    }
    else
    {
        d = adios_dims;
        *rank = 0;
    }

    if (global)
    {
        g = global->dimensions;
        o = global->offsets;
    }

    // do the conversion
    while (d)
    {
        // first check to make sure all vars are provided
        if (   d->dimension.var
            && !d->dimension.var->data
           )
        {
            fprintf (stderr, "sizing of %s failed because "
                             "dimension component %s was not "
                             "provided\n"
                    ,name, d->dimension.var->name
                    );

            break;
        }

        // calculate the size for this dimension element
        if (d->dimension.var)
        {
            bp_dims [*rank].local_bound =
                              (*(int *) d->dimension.var->data);
        }
        else
        {
            bp_dims [*rank].local_bound = d->dimension.rank;
        }
        if (g)
        {
            if (g->dimension.var)
            {
                bp_dims [*rank].global_bound =
                           (*(int *) g->dimension.var->data);
            }
            else
            {
                bp_dims [*rank].global_bound = g->dimension.rank;
            }
        }
        else
        {
            bp_dims [*rank].global_bound = 0;
        }
        if (o)
        {
            if (o->dimension.var)
            {
                bp_dims [*rank].global_offset =
                           (*(int *) o->dimension.var->data);
            }
            else
            {
                bp_dims [*rank].global_offset = o->dimension.rank;
            }
        }
        else
        {
            bp_dims [*rank].global_offset = 0;
        }

        //printf ("\tdim (%d): %d(%d)[%d]\n", *rank
        //       ,bp_dims [*rank].local_bound
        //       ,bp_dims [*rank].global_bound
        //       ,bp_dims [*rank].global_offset
        //       );
        d = d->next;
        if (g)
            g = g->next;
        if (o)
            o = o->next;
        (*rank)++;
    }
}

///////////////////////////////////////////////////////////////////////////////
// functions to support C & Fortran interface
///////////////////////////////////////////////////////////////////////////////
int adios_common_declare_group (long long * id, const char * name
                               ,enum ADIOS_FLAG host_language_fortran
                               ,const char * coordination_comm
                               ,const char * coordination_var
                               )
{
    struct adios_group_struct * g = (struct adios_group_struct *)
                             malloc (sizeof (struct adios_group_struct));

    g->name = strdup (name);
    g->adios_host_language_fortran = host_language_fortran;
    g->id = 0; // will be set in adios_append_group
    g->var_count = 0;
    g->vars = 0;
    g->attributes = 0;
    g->group_by = (coordination_var ? strdup (coordination_var) : 0L); 
    g->group_comm = (coordination_comm ? strdup (coordination_comm) : 0L);
    g->methods = 0;
    g->mesh = 0;
    
    *id = (long long) g;
    
    adios_append_group (g);
    
    return 1;
}

int adios_common_define_global_bounds (long long group_id
                                      ,const char * dimensions
                                      ,const char * offsets
                                      ,struct adios_global_bounds_struct ** b
                                      )
{
    struct adios_group_struct * t = (struct adios_group_struct *) group_id;

    *b = (struct adios_global_bounds_struct *)
                        malloc (sizeof (struct adios_global_bounds_struct));
    char * dim;
    char * dim_temp;
    char * off;
    char * offsets_temp;

    if (dimensions)
        dim_temp = strdup (dimensions);
    else
        dim_temp = 0;

    if (offsets)
        offsets_temp = strdup (offsets);
    else
        offsets_temp = 0;

    (*b)->dimensions = 0;
    (*b)->offsets = 0;

    if (dim_temp)
    {
        dim = strtok (dim_temp, ",");
        while (dim)
        {
            struct adios_dimension_struct * d =
                     (struct adios_dimension_struct *)
                         calloc (1, sizeof (struct adios_dimension_struct));

            adios_parse_dimension (dim, t, d);

            adios_append_dimension (&(*b)->dimensions, d);

            dim = strtok (NULL, ",");
        }
        free (dim_temp);
    }

    if (offsets_temp)
    {
        off = strtok (offsets_temp, ",");
        while (off)
        {
            struct adios_dimension_struct * d =
                     (struct adios_dimension_struct *)
                         calloc (1, sizeof (struct adios_dimension_struct));

            adios_parse_dimension (off, t, d);

            adios_append_dimension (&(*b)->offsets, d);

            off = strtok (NULL, ",");
        }
        free (offsets_temp);
    }

    return 1;
}

int adios_common_define_var (long long group_id, const char * name
                            ,const char * path, int type
                            ,int copy_on_write
                            ,const char * dimensions
                            ,struct adios_global_bounds_struct * global_bounds
                            )
{
    struct adios_group_struct * t = (struct adios_group_struct *) group_id;
    struct adios_var_struct * v = (struct adios_var_struct *)
                               malloc (sizeof (struct adios_var_struct));
    char * dim;
    char * dim_temp;
    if (dimensions)
        dim_temp = strdup (dimensions);
    else
        dim_temp = 0;

    v->name = strdup (name);
    v->path = strdup (path);
    v->type = type;
    v->dimensions = 0;
    v->global_bounds = global_bounds;
    v->copy_on_write = copy_on_write;
    v->got_buffer = adios_flag_no;
    v->free_data = adios_flag_no;
    v->data = 0;
    v->data_size = 0;
    v->next = 0;

    if (dim_temp)
    {
        dim = strtok (dim_temp, ",");
        while (dim)
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
            adios_parse_dimension (dim, t, d);

            adios_append_dimension (&v->dimensions, d);

            dim = strtok (NULL, ",");
        }
        free (dim_temp);
    }

    adios_append_var (&t->vars, v);
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

    new_method->m = -1;
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
        if (!strcmp (g->group->name, name))
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

#if CHUNK_DATA
static int whole_data_only = 0;
/* copy data elements into the buffer without overflowing it */
static int chunk_data (struct adios_var_struct * f_param, char * buf, int buf_len)
{
                               /* where we are in the current list of vars */
    static struct adios_var_struct * f;
    static int var_offset;   /* what offset in current var to start copy */

    int copied_so_far = 0;
    int var_size;
    struct adios_var_struct * f_start;

    if (f_param)
    {
        f = f_param;
        var_offset = 0;
    }

    f_start = f;

    if (f)
    {
        var_size = adios_size_of_var (f, v->data);
    }
    else
    {
        return -1;
    }

    if (whole_data_only)
    {
        while (f && copied_so_far + var_size <= buf_len)
        {
            memcpy (buf, v->data, var_size);
            copied_so_far += var_size;
            buf += var_size;
            f = v->next;

            if (f)
            {
                var_size = adios_size_of_var (f, v->data);
            }
        }
    }
    else
    {
        while (v && copied_so_far < buf_len)
        {
            int to_copy = 0;
            int buf_left = buf_len - copied_so_far;

            if (buf_left >= (var_size - var_offset))
            {
                to_copy = var_size - var_offset;
            }
            else
            {
                to_copy = buf_left;
            }
            memcpy (buf, ((char *) v->data) + var_offset, to_copy);
            copied_so_far += to_copy;
            buf += to_copy;
            if (var_size == to_copy + var_offset)
            {
                var_offset = 0;
                v = v->next;
                if (v)
                {
                    var_size = adios_size_of_var (v, v->data);
                }
            }
            else
            {
                if (var_offset == 0)
                    var_offset = to_copy;
                else
                    var_offset += to_copy;
                /* copied_so_far should == buf_len */
            }
        }
    }

    if (f == f_start && copied_so_far == 0)
        return -2;
    else
        return copied_so_far;
}
#endif
