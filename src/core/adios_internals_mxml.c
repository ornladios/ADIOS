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

#include "public/adios.h"
#include "public/adios_error.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/adios_internals_mxml.h"
#include "core/buffer.h"
#include "core/adios_logger.h"
#include "core/util.h" // PairStruct*

#ifdef DMALLOC
#include "dmalloc.h"
#endif

static enum ADIOS_FLAG adios_host_language_fortran = adios_flag_yes;

extern struct adios_method_list_struct * adios_methods;
extern struct adios_group_list_struct * adios_groups;

//extern struct adios_method_list_struct * adios_methods;
//extern struct adios_group_list_struct * adios_groups;

struct adios_transport_struct * adios_transports = 0;
static int adios_transports_initialized = 0;


// ADIOS Schema: adding utility functions
void conca_att_nam(char ** returnstr, const char * meshname, char * att_nam);
void conca_numb_att_nam(char ** returnstr, const char * meshname, char * att_nam, char counterstr[5]);

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
        log_warn ("xml: duplicate attribute %s on %s (ignored)",n,en); \
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

    log_error ("config.xml: invalid type: %s in var %s\n", type, name);

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

    log_error ("config.xml: %s must have a value of 'yes' or 'no' "
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

static int parseMeshUniformDimensions (const char * dimensions
                                      ,struct adios_group_struct * new_group
                                      ,struct adios_mesh_uniform_struct * mesh
                                      ,const char * name
                                      )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_list_struct * item = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * dim_att_nam = 0; // dimensions attribute name
    char * getdimsfrom = 0; // dimensions attribute that is a var
    int counter = 0;        // used to create dimX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create dimX attributes

    if (!dimensions)
    {
        log_warn ("config.xml: dimensions value required for"
                         "uniform mesh: %s\n"
                         ,name
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
            log_warn ("Out of memory parseMeshUniformDimensions of mesh: %s\n"
                             ,name
                    );
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
                log_warn ("config.xml: invalid variable %s\n"
                                 "for dimensions of mesh: %s\n"
                                 ,c
                                 ,name
                        );
                free (d1);

                return 0;

            }else{
                // Found variable ==> create a dims attribute for it.
                // char * getdimsfrom = 0;
                getdimsfrom = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                conca_numb_att_nam(&getdimsfrom, name, "dim", counterstr);
                adios_common_define_attribute (p_new_group,getdimsfrom,"/",adios_string,item->item.var->name,"");
                free (getdimsfrom);
                counter++;
            }
        }
        else
        {
            item->item.rank = strtod (c, 0);
            counterstr[0] = '\0';
            snprintf(counterstr, 5, "%d", counter);
            dim_att_nam = 0;
            conca_numb_att_nam(&dim_att_nam, name, "dim", counterstr);
            adios_common_define_attribute (p_new_group,dim_att_nam,"/",adios_double,c,"");
            free (dim_att_nam);
            item->item.var = 0;
            counter++;
        }

        adios_append_mesh_item (&(mesh->dimensions), item);

        c = strtok (NULL, ",");
    }

    char * dims = 0;
    counterstr[0] = '\0';
    snprintf(counterstr, 5, "%d", counter);
    dims = 0;
    conca_att_nam(&dims, name, "ndims");

    adios_common_define_attribute (p_new_group,dims,"/",adios_double,counterstr,"");

    free (dims);

    free (d1);

    return 1;
}

static int parseMeshUniformMaxima (const char * maximum
                                      ,struct adios_group_struct * new_group
                                      ,struct adios_mesh_uniform_struct * mesh
                                      ,const char * name
                                      )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_list_struct * item = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * max_att_nam = 0; // maxima attribute name
    char * getmaxafrom = 0; // maxima attribute name that is a var
    int counter = 0;        // used to create maxX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create maxX attributes

    if (!maximum)
    {
         log_warn ("config.xml: maximum value required"
                         "for uniform mesh: %s\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (maximum);

    c = strtok (d1, ",");

    while (c)
    {
        item = (struct adios_mesh_item_list_struct *) malloc
                            (sizeof (struct adios_mesh_item_list_struct));
        item->next = 0;

        if (!item)
        {
            log_warn ("Out of memory parseMeshUniformMaxima of mesh: %s\n"
                             ,name
                    );
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
                log_warn ("config.xml: invalid variable %s for maximum\n"
                                 "of mesh: %s\n"
                                 ,c
                                 ,name
                        );
                free (d1);

                return 0;
            }else{
                // Found variable ==> create a maxa attribute for it.
                // char * getmaxafrom = 0;
                getmaxafrom = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                conca_numb_att_nam(&getmaxafrom, name, "max", counterstr);
                adios_common_define_attribute (p_new_group,getmaxafrom,"/",adios_string,item->item.var->name,"");
                free (getmaxafrom);
                counter++;
            }
        }
        else
        {
            // Create attributes for each maximum
            item->item.rank = strtod (c, 0);
            counterstr[0] = '\0';
            snprintf(counterstr, 5, "%d", counter);
            max_att_nam = 0;
            conca_numb_att_nam(&max_att_nam, name, "max", counterstr);
            adios_common_define_attribute (p_new_group,max_att_nam,"/",adios_double,c,"");
            free (max_att_nam);
            item->item.var = 0;
            counter++;
        }

        adios_append_mesh_item (&(mesh->maximum), item);

        c = strtok (NULL, ",");
    }

    char * maxa = 0;
    counterstr[0] = '\0';
    snprintf(counterstr, 5, "%d", counter);
    maxa = 0;
    conca_att_nam(&maxa, name, "maxa");
    adios_common_define_attribute (p_new_group,maxa,"/",adios_double,counterstr,"");
    free (maxa);
    free (d1);

    return 1;
}

static int parseMeshUniformOrigin (const char * origin
                                  ,struct adios_group_struct * new_group
                                  ,struct adios_mesh_uniform_struct * mesh
                                  ,const char * name
                                  )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_list_struct * item = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * org_att_nam = 0; // origins attribute name
    char * getorgsfrom = 0; // origins attribute name that is a var
    int counter = 0;        // used to create orgX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create orgX attributes

    if (!origin)
    {
        log_warn ("config.xml: origin value required "
                         "for uniform mesh: %s\n"
                         ,name
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
            log_warn ("Out of memory parseMeshUniformOrigin of mesh: %s\n", name);
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
                log_warn ("config.xml: invalid variable %s for origin\n"
                        "of mesh: %s\n"
                        ,c
                        ,name
                        );
                free (d1);

                return 0;
            }else{
                // Found variable ==> create a orgs attribute for it.
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                getorgsfrom = 0;
                conca_numb_att_nam(&getorgsfrom, name, "org", counterstr);
                adios_common_define_attribute (p_new_group,getorgsfrom,"/",adios_string,item->item.var->name,"");
                free (getorgsfrom);
                counter++;
            }
        }
        else
        {
            // Create attributes for each origin
            item->item.rank = strtod (c, 0);
            item->item.var = 0;
            counterstr[0] = '\0';
            snprintf(counterstr, 5, "%d", counter);
            org_att_nam = 0;
            conca_numb_att_nam(&org_att_nam, name, "org", counterstr);
            adios_common_define_attribute (p_new_group,org_att_nam,"/",adios_double,c,"");
            free (org_att_nam);
            counter++;
        }

        adios_append_mesh_item (&(mesh->origin), item);

        c = strtok (NULL, ",");
    }

    char * orgs = 0;
    counterstr[0] = '\0';
    snprintf(counterstr, 5, "%d", counter);
    orgs = 0;
    conca_att_nam(&orgs, name, "orgs");
    adios_common_define_attribute (p_new_group,orgs,"/",adios_double,counterstr,"");
    free (orgs);

    free (d1);

    return 1;
}

static int parseMeshUniformSpacing (const char * spacing
                                   ,struct adios_group_struct * new_group
                                   ,struct adios_mesh_uniform_struct * mesh
                                   ,const char * name
                                   )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_list_struct * item = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * spa_att_nam = 0; // spacings attribute name
    char * getspasfrom = 0; // spacings attribute name that is a var
    int counter = 0;        // used to create spaX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create spaX attributes if (!spacing)

    if (!spacing)
    {
        log_warn ("config.xml: mesh uniform spacing value "
                "required for mesh: %s\n"
                ,name
                );
        return 0;
    }

    d1 = strdup (spacing);

    c = strtok (d1, ",");

    while (c)
    {
        item = (struct adios_mesh_item_list_struct *) malloc
                            (sizeof (struct adios_mesh_item_list_struct));
        item->next = 0;

        if (!item)
        {
            log_warn ("Out of memory parseMeshUniformSpacing for mesh: %s\n", name);
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
                log_warn ("config.xml: invalid variable: %s for spacing\n"
                                 "of mesh: %s\n"
                                 ,c
                                 ,name
                        );
                free (d1);

                return 0;
            }else{
                // Found variable ==> create a spas attribute for it.
                // char * getspasfrom = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                getspasfrom = 0;
                conca_numb_att_nam(&getspasfrom, name, "spa", counterstr);
                adios_common_define_attribute (p_new_group,getspasfrom,"/",adios_string,item->item.var->name,"");
                free (getspasfrom);
                counter++;
            }
        }
        else
        {
            // Create attributes for each dimension
            item->item.rank = strtod (c, 0);
            counterstr[0] = '\0';
            snprintf(counterstr, 5, "%d", counter);
            spa_att_nam = 0;
            conca_numb_att_nam(&spa_att_nam, name, "spa", counterstr);
            adios_common_define_attribute (p_new_group,spa_att_nam,"/",adios_double,c,"");
            free (spa_att_nam);
            item->item.var = 0;
            counter++;
        }

        adios_append_mesh_item (&(mesh->origin), item);

        c = strtok (NULL, ",");
    }

    char * spas = 0;
    counterstr[0] = '\0';
    snprintf(counterstr, 5, "%d", counter);
    spas = 0;
    conca_att_nam(&spas, name, "spas");
    adios_common_define_attribute (p_new_group,spas,"/",adios_double,counterstr,"");
    free (spas);

    free (d1);

    return 1;
}

static int parseMeshRectilinearDimensions (const char * dimensions
                                          ,struct adios_group_struct * new_group
                                          ,struct adios_mesh_rectilinear_struct * mesh
                                          ,const char * name
                                          )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_list_struct * item = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * dim_att_nam = 0; // dimensions attribute name
    char * getdimsfrom = 0; // dimensions attribute name that is a var
    int counter = 0;        // used to create dimX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create dimX attributes

    if (!dimensions)
    {
        log_warn ("config.xml: dimensions value required"
                         "for rectilinear mesh: %s\n"
                         ,name
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
            log_warn ("Out of memory parseMeshRectilinearDimensions for mesh: %s\n"
                             ,name
                             );
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
                log_warn ("config.xml: invalid variable: %s for dimensions\n"
                                 "of mesh: %s\n"
                                 ,c
                                 ,name
                        );
                free (d1);

                return 0;
            }else
            {
                // Found variable ==> create a dims attribute for it.
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                getdimsfrom = 0;
                conca_numb_att_nam(&getdimsfrom, name, "dim", counterstr);
                adios_common_define_attribute (p_new_group,getdimsfrom,"/",adios_string,item->item.var->name,"");
                free (getdimsfrom);
                counter++;
            }
        }
        else
        {
            // Create attributes for each dimension
            item->item.rank = strtod (c, 0);
            counterstr[0] = '\0';
            snprintf(counterstr, 5, "%d", counter);
            dim_att_nam = 0;
            conca_numb_att_nam(&dim_att_nam, name, "dim", counterstr);
            adios_common_define_attribute (p_new_group,dim_att_nam,"/",adios_double,c,"");
            free (dim_att_nam);
            item->item.var = 0;
            counter++;
        }

        adios_append_mesh_item (&(mesh->dimensions), item);

        c = strtok (NULL, ",");
    }

    char * dims = 0;
    counterstr[0] = '\0';
    snprintf(counterstr, 5, "%d", counter);
    dims = 0;
    conca_att_nam(&dims, name, "dims");
    adios_common_define_attribute (p_new_group,dims,"/",adios_double,counterstr,"");

    free (dims);

    free (d1);

    return 1;
}

static int parseMeshRectilinearCoordinatesMultiVar (const char * coordinates
                                                   ,struct adios_group_struct * new_group
                                                   ,struct adios_mesh_rectilinear_struct * mesh
                                                   ,const char * name
                                                   )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_var_list_struct * var = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * coo_att_nam = 0; // coordinates attribute name
    int counter = 0;        // used to create ptsX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create ptsX attributes

    if (!coordinates)
    {
        log_warn ("config.xml: coordinates-multi-var value required"
                         "for rectilinear mesh: %s\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (coordinates);

    c = strtok (d1, ",");

    while (c)
    {
        var = (struct adios_mesh_var_list_struct *) malloc
                        (sizeof (struct adios_mesh_var_list_struct));
        var->next = 0;

        if (!var)
        {
            log_warn ("Out of memory parseMeshRectilinearCoordinatesMultiVar\n"
                             "for mesh: %s"
                             ,name
                    );
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
                log_warn ("config.xml: invalid variable: %s for coordinates\n"
                                 "of mesh: %s\n"
                                 ,c
                                 ,name
                        );
                free (d1);

                return 0;
            }else
            {
                // Found variable ==> create a coords attribute for it.
                coo_att_nam = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                conca_numb_att_nam(&coo_att_nam, name, "pts", counterstr);
                adios_common_define_attribute (p_new_group,coo_att_nam,"/",adios_string,c,"");
                free (coo_att_nam);
                counter++;
            }
        }
        else
        {
            var->var = 0;
            log_warn ("config.xml: invalid variable: %s for coordinates"
                             "of rectilinear mesh: %s\n"
                             ,c
                             ,name
                    );
            free (d1);
            return 0;
        }

        adios_append_mesh_var (&(mesh->coordinates), var);

        c = strtok (NULL, ",");
    }

    // At this points, coordinates should point to at least 2 variables
    // otherwise let the user know to use the coordinates-single-var tag
    if (counter > 1) {
        char * coords = 0;
        counterstr[0] = '\0';
        snprintf(counterstr, 5, "%d", counter);
        conca_att_nam(&coords, name, "nvars");
        adios_common_define_attribute (p_new_group,coords,"/",adios_double,counterstr,"");
        free (coords);
    } else
    {
        log_warn ("config.xml: coordinates-multi-var expects "
                         "at least 2 variables (%s)\n"
                         ,name
                );
        free (d1);
        return 0;
    }

    free (d1);

    return 1;
}

static int parseMeshRectilinearCoordinatesSingleVar (const char * coordinates
                                                    ,struct adios_group_struct * new_group
                                                    ,struct adios_mesh_rectilinear_struct * mesh
                                                    ,const char * name
                                                    )
{
    char * d1; // save of strdup
    struct adios_mesh_var_list_struct * var = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * coo_att_nam = 0; // coordinates attribute name

    if (!coordinates)
    {
        log_warn ("config.xml: coordinates-single-var value required"
                         "for rectilinear mesh: %s\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (coordinates);

    var = (struct adios_mesh_var_list_struct *) malloc
                    (sizeof (struct adios_mesh_var_list_struct));
    var->next = 0;

    if (!var)
    {
        log_warn ("Out of memory parseMeschRectilinearCoordinatesSingleVar\n");
        free (d1);

        return 0;
     }

     if (adios_int_is_var (d1))
     {
        var->var =
               adios_find_var_by_name (new_group->vars, d1
                                      ,new_group->all_unique_var_names
                                      );
        if (!var->var)
        {
            log_warn ("config.xml: invalid variable: %s for coordinates"
                             "of mesh: %s\n"
                             ,d1
                             ,name
                    );
            free (d1);

            return 0;
        }else
        {
            // Found variable ==> create a nvars attribute for it.
            conca_att_nam(&coo_att_nam, name, "nvars");
            adios_common_define_attribute (p_new_group,coo_att_nam,"/",adios_string,d1,"");
            free (coo_att_nam);
        }
    }
    else
    {
        var->var = 0;
        log_warn ("config.xml: invalid variable: %s for coordinates of mesh: %s\n"
                                 ,d1
                                 ,name
                );
        free (d1);
        return 0;
    }

    adios_append_mesh_var (&(mesh->coordinates), var);

    free (d1);

    return 1;
}

static int parseMeshStructuredNspace (const char * nspace
                                     ,struct adios_group_struct * new_group
                                     ,struct adios_mesh_structured_struct * mesh
                                     ,const char * name
                                     )
{
    char * d1; // save of strdup
    int64_t      p_new_group = (int64_t) new_group;
    struct adios_mesh_item_struct * item = 0;
    char * nsp_att_nam = 0; // nspace attribute name

    if (!nspace)
    {
        log_warn ("config.xml: npsace value required "
                         "for structured mesh: %s\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (nspace);

    item = (struct adios_mesh_item_struct *) malloc
                            (sizeof (struct adios_mesh_item_struct));

    if (!item)
    {
        log_warn ("Out of memory parseMeshStructuredNspace for mesh: %s\n", name);
        free (d1);

        return 0;
    }

    conca_att_nam(&nsp_att_nam, name, "nsp");

    if (adios_int_is_var (nspace))
    {
        item->rank = 0.0;
        item->var =
               adios_find_var_by_name (new_group->vars, nspace
                                       ,new_group->all_unique_var_names
                                       );
        if (!item->var)
        {
           log_warn ("config.xml: invalid variable: %s for nspace of mesh: %s\n"
                            ,nspace
                            ,name
                   );
            free (d1);

            return 0;
        }else
        {
            adios_common_define_attribute (p_new_group,nsp_att_nam,"/",adios_string,item->var->name,"");
            free (nsp_att_nam);
        }
    }
    else
    {
        item->rank = strtod (nspace, 0);
        adios_common_define_attribute (p_new_group,nsp_att_nam,"/",adios_double,nspace,"");
        free (nsp_att_nam);
        item->var = 0;
    }

    mesh->nspace = item;

    free (d1);

    return 1;
}

static int parseMeshStructuredDimensions (const char * dimensions
                                         ,struct adios_group_struct * new_group
                                         ,struct adios_mesh_structured_struct * mesh
                                         ,const char * name
                                         )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_item_list_struct * item = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * dim_att_nam = 0; // dimensions attribute name
    char * getdimsfrom = 0; // dimensions attribute name that is a var
    int counter = 0;        // used to create dimX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create dimX attributes

    if (!dimensions)
    {
        log_warn ("config.xml: dimensions value required"
                         "for structured mesh: %s\n"
                         ,name
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
            log_warn ("Out of memory parseMeshStructuredDimensions of mesh: %s\n"
                             ,name
                    );
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
                log_warn ("config.xml: invalid variable: %s for dimensions\n"
                                 "of mesh: %s\n"
                                 ,c
                                 ,name
                        );

                free (d1);

                return 0;
            } else
            {
                // Found variable ==> create a dims attribute for it.
                // char * getdimsfrom = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                getdimsfrom = 0;
                conca_numb_att_nam(&getdimsfrom, name, "dim", counterstr);
                adios_common_define_attribute (p_new_group,getdimsfrom,"/",adios_string,item->item.var->name,"");
                free (getdimsfrom);
                counter++;
            }
        }
        else
        {
            item->item.rank = strtod (c, 0);
            item->item.var = 0;
            counterstr[0] = '\0';
            snprintf(counterstr, 5, "%d", counter);
            dim_att_nam = 0;
            conca_numb_att_nam(&dim_att_nam, name, "dim", counterstr);
            adios_common_define_attribute (p_new_group,dim_att_nam,"/",adios_double,c,"");
            free (dim_att_nam);
            counter++;

        }

        adios_append_mesh_item (&(mesh->dimensions), item);

        c = strtok (NULL, ",");
    }

    char * dims = 0;
    counterstr[0] = '\0';
    snprintf(counterstr, 5, "%d", counter);
    dims = 0;
    conca_att_nam(&dims, name, "ndims");
    adios_common_define_attribute (p_new_group,dims,"/",adios_double,counterstr,"");

    free (dims);

    free (d1);

    return 1;
}

static int parseMeshStructuredPointsMultiVar (const char * points
                                             ,struct adios_group_struct * new_group
                                             ,struct adios_mesh_structured_struct * mesh
                                             ,const char * name
                                             )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_var_list_struct * var = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * pts_att_nam = 0; // pointss attribute name
    int counter = 0;        // used to create ptsX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create ptsX attributes

    if (!points)
    {
        log_warn ("config.xml: points-multi-var value required"
                         "for structured mesh: %s\n"
                         ,name
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
            log_warn ("Out of memory parseMeshStructuredPointsMultiVar"
                             "of mesh: %s\n"
                             ,name
                    );
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
                log_warn ("config.xml: invalid variable: %s points of mesh: %s\n"
                                 ,c
                                 ,name
                        );
                free (d1);

                return 0;
            }else
            {
                // Found variable ==> create a points attribute for it.
                pts_att_nam = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                conca_numb_att_nam(&pts_att_nam, name, "pts", counterstr);
                adios_common_define_attribute (p_new_group,pts_att_nam,"/",adios_string,c,"");
                free (pts_att_nam);
                counter++;
            }
        }
        else
        {
            var->var = 0;
            log_warn ("config.xml: invalid variable: %s for points of mesh: %s\n"
                             ,c
                             ,name
                    );
            free (d1);

            return 0;
        }

        adios_append_mesh_var (&(mesh->points), var);

        c = strtok (NULL, ",");
    }

    // Define an attribute showing the number of mesh_vars
    // Should be more than one in this multi-var parsing
    if (counter > 1){
        char * pts = 0;
        counterstr[0] = '\0';
        snprintf(counterstr, 5, "%d", counter);
        conca_att_nam(&pts, name, "nvars");
        adios_common_define_attribute (p_new_group,pts,"/",adios_double,counterstr,"");
        free (pts);
    } else
    {
            log_warn ("config.xml: points-multi-var tag for mesh: %s "
                             " expects at least 2 variables\n"
                             ,name
                    );
            free (d1);

            return 0;
    }

    free (d1);

    return 1;
}

static int parseMeshStructuredPointsSingleVar (const char * points
                                              ,struct adios_group_struct * new_group
                                              ,struct adios_mesh_structured_struct * mesh
                                              ,const char * name
                                              )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_var_list_struct * var = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * pts_att_nam = 0; // points attribute name

    if (!points)
    {
        log_warn ("config.xml: points-single-var value required"
                         "for structured mesh: %s\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (points);

    var = (struct adios_mesh_var_list_struct *) malloc
                    (sizeof (struct adios_mesh_var_list_struct));
    var->next = 0;

    if (!var)
    {
        log_warn ("Out of memory parseMeshStructuredPointsSingleVar"
                         "for mesh: %s\n"
                         ,name
                );
        free (d1);

        return 0;
     }

     if (adios_int_is_var (d1))
     {
        var->var =
               adios_find_var_by_name (new_group->vars, d1
                                      ,new_group->all_unique_var_names
                                      );
        if (!var->var)
        {
            log_warn ("config.xml: invalid variable: %s for points of mesh: %s\n"
                             ,d1
                             ,name
                    );
            free (d1);

            return 0;
        }else
        {
            // Found variable ==> create a number of vars attribute for it.
            conca_att_nam(&pts_att_nam, name, "nvars");
            adios_common_define_attribute (p_new_group,pts_att_nam,"/",adios_string,d1,"");
            free (pts_att_nam);
        }
    }
    else
    {
        var->var = 0;
        log_warn ("config.xml: invalid variable: %s for points of mesh: %s\n"
                         ,d1
                         ,name
                );
        free (d1);

        return 0;
    }

    adios_append_mesh_var (&(mesh->points), var);

    free (d1);

    return 1;
}

static int parseMeshUnstructuredNspace (const char * nspace
                                     ,struct adios_group_struct * new_group
                                     ,struct adios_mesh_unstructured_struct * mesh
                                     ,const char * name
                                     )
{
    char * d1; // save of strdup
    int64_t      p_new_group = (int64_t) new_group;
    struct adios_mesh_item_struct * item = 0;
    char * nsp_att_nam = 0; // nspace attribute name

    if (!nspace)
    {
        log_warn ("config.xml: nspace value required"
                         "for unstructured mesh: %s\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (nspace);

    item = (struct adios_mesh_item_struct *) malloc
                            (sizeof (struct adios_mesh_item_struct));

    if (!item)
    {
        log_warn ("Out of memory parseMeshUnstructuredNspace for mesh: %s\n"
                         ,name
                );
        free (d1);

        return 0;
    }

    conca_att_nam(&nsp_att_nam, name, "nsp");

    if (adios_int_is_var (nspace))
    {
        item->rank = 0.0;
        item->var =
               adios_find_var_by_name (new_group->vars, nspace
                                       ,new_group->all_unique_var_names
                                       );
        if (!item->var)
        {
           log_warn ("config.xml: invalid variable: %s for nspace of mesh: %s\n"
                            ,nspace
                            ,name
                   );
            free (d1);

            return 0;
        }else
        {
            adios_common_define_attribute (p_new_group,nsp_att_nam,"/",adios_string,item->var->name,"");
            free (nsp_att_nam);
        }
    }
    else
    {
        item->rank = strtod (nspace, 0);
        adios_common_define_attribute (p_new_group,nsp_att_nam,"/",adios_double,nspace,"");
        free (nsp_att_nam);
        item->var = 0;
    }

    mesh->nspace = item;

    free (d1);

    return 1;
}

static int parseMeshUnstructuredNpoints (const char * npoints
                                     ,struct adios_group_struct * new_group
                                     ,struct adios_mesh_unstructured_struct * mesh
                                     ,const char * name
                                     )
{
    char * d1; // save of strdup
    int64_t      p_new_group = (int64_t) new_group;
    struct adios_mesh_item_struct * item = 0;
    char * npts_att_nam = 0; // npoints attribute name

    if (!npoints)
    {
        log_warn ("config.xml: npoints value required"
                         "for unstructured mesh\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (npoints);

    item = (struct adios_mesh_item_struct *) malloc
                            (sizeof (struct adios_mesh_item_struct));

    if (!item)
    {
        log_warn ("Out of memory parseMeshUnstructuredNpoints for mesh: %s\n"
                         ,name
                );
        free (d1);

        return 0;
    }

    conca_att_nam(&npts_att_nam, name, "npoints");

    if (adios_int_is_var (npoints))
    {
        item->rank = 0.0;
        item->var =
               adios_find_var_by_name (new_group->vars, npoints
                                       ,new_group->all_unique_var_names
                                       );
        if (!item->var)
        {
           log_warn ("config.xml: invalid variable: %s for npoints of mesh: %s\n"
                   ,npoints
                   ,name
                   );
            free (d1);

            return 0;
        }else
        {
            adios_common_define_attribute (p_new_group,npts_att_nam,"/",adios_string,item->var->name,"");
            free (npts_att_nam);
        }
    }
    else
    {
        item->rank = strtod (npoints, 0);
        adios_common_define_attribute (p_new_group,npts_att_nam,"/",adios_double,npoints,"");
        free (npts_att_nam);
        item->var = 0;
    }

    mesh->points_count = item;

    free (d1);

    return 1;
}

static int parseMeshUnstructuredPointsMultiVar (const char * points
                                             ,struct adios_group_struct * new_group
                                             ,struct adios_mesh_unstructured_struct * mesh
                                             ,const char * name
                                             )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_var_list_struct * var = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * pts_att_nam = 0; // pointss attribute name
    int counter = 0;        // used to create ptsX attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create ptsX attributes

    if (!points)
    {
        log_warn ("config.xml: points-multi-var value required"
                         "for unstructured mesh: %s\n"
                         ,name
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
            log_warn ("Out of memory parseMeshUnstructuredPointsMultiVar\n"
                             "for mesh: %s\n"
                             ,name
                    );
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
                log_warn ("config.xml: invalid variable: %s for points of mesh: %s\n"
                        ,c
                        ,name
                        );
                free (d1);

                return 0;
            }else
            {
                // Found variable ==> create a points attribute for it.
                pts_att_nam = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                conca_numb_att_nam(&pts_att_nam, name, "pts", counterstr);
                adios_common_define_attribute (p_new_group,pts_att_nam,"/",adios_string,c,"");
                free (pts_att_nam);
                counter++;
            }
        }
        else
        {
            var->var = 0;
            log_warn ("config.xml: invalid variable: %s for points of mesh: %s\n"
                             ,c
                             ,name
                    );
            free (d1);
            return 0;
        }

        adios_append_mesh_var (&(mesh->points), var);

        c = strtok (NULL, ",");
    }

    // At this point we expect at least 2 "points-multi-var values
    if (counter > 1){
        char * pts = 0;
        counterstr[0] = '\0';
        snprintf(counterstr, 5, "%d", counter);
        conca_att_nam(&pts, name, "nvars");
        adios_common_define_attribute (p_new_group,pts,"/",adios_double,counterstr,"");
        free (pts);
        free (d1);
    } else
    {
        log_warn ("config.xml: points-multi-var tag expects "
                         " at least two variabels. (%s)\n"
                ,name
                );
        free (d1);
        return 0;
    }

    return 1;
}

static int parseMeshUnstructuredPointsSingleVar (const char * points
                                              ,struct adios_group_struct * new_group
                                              ,struct adios_mesh_unstructured_struct * mesh
                                              ,const char * name
                                              )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_var_list_struct * var = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * pts_att_nam = 0; // points attribute name

    if (!points)
    {
        log_warn ("config.xml: points-single-var value required"
                         "for unstructured mesh: %s\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (points);

    var = (struct adios_mesh_var_list_struct *) malloc
                    (sizeof (struct adios_mesh_var_list_struct));
    var->next = 0;

    if (!var)
    {
        log_warn ("Out of memory parseMeshUnstructuredPointsSingleVar"
                         "for mesh: %s\n"
                         ,name
                );
        free (d1);

        return 0;
     }

     if (adios_int_is_var (d1))
     {
        var->var =
               adios_find_var_by_name (new_group->vars, d1
                                      ,new_group->all_unique_var_names
                                      );
        if (!var->var)
        {
            log_warn ("config.xml: invalid variable: %s for points of mesh: %s\n"
                    ,d1
                    ,name
                    );
            free (d1);

            return 0;
        }else
        {
            // Found variable ==> create a number of vars attribute for it.
            conca_att_nam(&pts_att_nam, name, "nvars");
            adios_common_define_attribute (p_new_group,pts_att_nam,"/",adios_string,d1,"");
            free (pts_att_nam);
        }
    }
    else
    {
        var->var = 0;
        log_warn ("config.xml: invalid variable: %s "
                         " for points-multi-var of mesh: %s\n"
                ,d1
                ,name
                );
        free (d1);
        return 0;
    }

    adios_append_mesh_var (&(mesh->points), var);

    free (d1);

    return 1;
}

static int parseMeshUnstructuredUniformCells (const char * count
                                             ,const char * data
                                             ,const char * type
                                             ,struct adios_group_struct * new_group
                                             ,struct adios_mesh_unstructured_struct * mesh
                                             ,const char * name
                                             )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_cell_list_list_struct * cell_list = 0;
    struct adios_mesh_item_struct * item = 0;
    int64_t      p_new_group = (int64_t) new_group;
    char * ncellset_att_nam = 0;  // ncellset attribute
    char * cellcount_att_nam = 0; // single cell count attribute
    char * celldata_att_nam = 0;  // single cell data  attribute
    char * celltype_att_nam = 0;  // single cell type attribute

    item = (struct adios_mesh_item_struct *) malloc
                            (sizeof (struct adios_mesh_item_struct));

    if (!item)
    {
        log_warn ("Out of memory parseMeshUnstructuredUniformCells of mesh: %s\n"
                         ,name
                );
        free (d1);

        return 0;
    }else{
        item->rank = 1;
        conca_att_nam(&ncellset_att_nam,name,"ncsets");
        adios_common_define_attribute (p_new_group,ncellset_att_nam,"/",adios_double,"1","");
        free (ncellset_att_nam);
        item->var = 0;
    }

    // Given that we expect only one variable in this case
    // We should not have any "," in these values... Catch errors
    mesh->cell_set_count = item;

    if (!count)
    {
        log_warn ("config.xml: uniform-cells count value required"
                         "for unstructured mesh: %s\n"
                         ,name
                );

        return 0;
    }
    if (!data)
    {
        log_warn ("config.xml: uniform-cells data value required"
                         "for unstructured mesh: %s\n"
                         ,name
                );

        return 0;
    }
    if (!type)
    {
        log_warn ("config.xml: uniform-cells type value required"
                         "for unstructured mesh: %s\n"
                         ,name
                );

        return 0;
    }

    d1 = strdup (count);

    cell_list = (struct adios_mesh_cell_list_list_struct *) malloc
                        (sizeof (struct adios_mesh_cell_list_list_struct));
    cell_list->next = 0;

    if (!cell_list)
    {
        log_warn ("Out of memory parseMeshUnstructuredUniformCells"
                         "of mesh: %s\n"
                         ,name
                );
        free (d1);

        return 0;
    }

    if (adios_int_is_var (d1))
    {
        cell_list->cell_list.count.var =
                adios_find_var_by_name (new_group->vars, d1
                                       ,new_group->all_unique_var_names
                                       );
        cell_list->cell_list.count.rank = 0;
        if (!cell_list->cell_list.count.var)
        {
            log_warn ("config.xml: invalid variable: %s for 'count' attribute\n"
                             "of mesh: %s\n"
                             ,d1
                             ,name
                    );
            free (d1);

            return 0;
        } else
        {
           conca_att_nam(&cellcount_att_nam, name, "ccount");
           adios_common_define_attribute (p_new_group,cellcount_att_nam,"/",adios_string,cell_list->cell_list.count.var->name,"");
           free (cellcount_att_nam);
        }
    }
    else
    {
        cell_list->cell_list.count.var = 0;
        cell_list->cell_list.count.rank = strtod (d1, 0);
        conca_att_nam(&cellcount_att_nam, name, "ccount");
        adios_common_define_attribute (p_new_group,cellcount_att_nam,"/",adios_double,d1,"");
        free (cellcount_att_nam);
    }

    free (d1);

    d1 = strdup (data);

    if (adios_int_is_var (d1))
    {
        cell_list->cell_list.data =
             adios_find_var_by_name (new_group->vars, d1
                                    ,new_group->all_unique_var_names
                                    );
        if (!cell_list->cell_list.data)
        {
            log_warn ("config.xml: invalid variable: %s for cell data\n"
                             "of mesh: %s\n"
                             ,d1
                             ,name
                    );
            free (d1);

            return 0;
        } else
        {
            conca_att_nam(&celldata_att_nam, name, "cdata");
            adios_common_define_attribute (p_new_group,celldata_att_nam,"/",adios_string,d1,"");
            free (celldata_att_nam);
        }
    }
    else
    {
        cell_list->cell_list.data = 0;
        log_warn ("config.xml: invalid data dimension: %s of mesh: %s\n"
                         "Please provide a valid variable for cell data.\n"
                         ,d1
                         ,name
                );
        free (d1);
        return 0;
    }

    free (d1);

    d1 = strdup (type);

    conca_att_nam(&celltype_att_nam, name, "ctype");

    if (!strcmp(d1,"pt") || !strcmp(d1,"point"))
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = 1; //strtod (c, 0);
        adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
    } else
    if (!strcmp(d1,"line"))
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = 2; //strtod (c, 0);
        adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
    } else
    if (!strcmp(d1,"tri") || !strcmp(d1,"triangle"))
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = 3; //strtod (c, 0);
        adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
    } else
    if (!strcmp(d1,"quad") || !strcmp(d1,"quadrilateral"))
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = 4; //strtod (c, 0);
        adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
    } else
    if (!strcmp(d1,"hex") || !strcmp(d1,"hexahedron"))
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = 5; //strtod (c, 0);
        adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
    } else
    if (!strcmp(d1,"pri") || !strcmp(d1,"prism"))
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = 6; //strtod (c, 0);
        adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
    } else
    if (!strcmp(d1,"tet") || !strcmp(d1,"tetrahedron"))
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = 7; //strtod (c, 0);
        adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
    } else
    if (!strcmp(d1,"pyr") || !strcmp(d1,"pyramid"))
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = 8; //strtod (c, 0);
        adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
    } else
    if (adios_int_is_var (d1))
    {
        cell_list->cell_list.type.var =
            adios_find_var_by_name (new_group->vars, d1
                    ,new_group->all_unique_var_names
                    );
        cell_list->cell_list.type.rank = 0;
        if (!cell_list->cell_list.type.var)
        {
            log_warn ("config.xml: invalid variable: %s for type attribute\n"
                             "of mesh: %s\n"
                             ,d1
                             ,name
                    );
            free (celltype_att_nam);
            free (d1);

            return 0;
        } else
        {
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,d1,"");
        }
    } else
    {
        cell_list->cell_list.type.var = 0;
        cell_list->cell_list.type.rank = strtod (d1, 0);
        if (cell_list->cell_list.type.rank > 0 && cell_list->cell_list.type.rank < 9)
        {
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_double,d1,"");
        } else {
            log_warn ("config.xml: invalid type attribute: %s of mesh: %s.\n"
                             "Please select a number between 1 and 8 for cell types\n"
                             "or: pt, line, tri, quad, hex, pri, tet or pyr.\n"
                             ,d1
                             ,name
                    );
            free (celltype_att_nam);
            free (d1);

            return 0;
        }
    }

    free(celltype_att_nam);

    free (d1);

    adios_append_mesh_cell_list (&(mesh->cell_list), cell_list);

    return 1;
}

static int parseMeshUnstructuredMixedCells (const char * count
                                           ,const char * data
                                           ,const char * types
                                           ,struct adios_group_struct * new_group
                                           ,struct adios_mesh_unstructured_struct * mesh
                                           ,const char * name
                                           )
{
    char * c;  // comma location
    char * d1; // save of strdup
    struct adios_mesh_cell_list_list_struct * cell_list = 0;
    struct adios_mesh_item_struct * item = 0;
    int counter = 0;        // used to create countX, typeX, dataX? attributes
    char counterstr[5] = {0,0,0,0,0}; // used to create countX, typeX, dataX? attributes
    int64_t      p_new_group = (int64_t) new_group;
    char * ncellset_att_nam = 0;  // ncellset attribute
    char * ccounts_att_nam = 0;   // ccountX attributes
    char * cdata_att_nam = 0;     // cdataX attributes
    char * celltype_att_nam = 0;  // ctypeX attributes

    item = (struct adios_mesh_item_struct *) malloc
                            (sizeof (struct adios_mesh_item_struct));

    if (!item)
    {
        log_warn ("Out of memory parseMeshUnstructuredMixedCells of mesh: %s\n"
                         ,name
                );
        free (d1);

        return 0;
    }

    if (!count)
    {
        log_warn ("config.xml: mixed-cells count value required"
                         "for unstructured mesh: %s\n"
                         ,name
                );

        return 0;
    }
    if (!data)
    {
        log_warn ("config.xml: mixed-cells data value required"
                         "for unstructured mesh: %s\n"
                         ,name
                );

        return 0;
    }
    if (!types)
    {
        log_warn ("config.xml: mixed-cellsi type value required "
                         "for unstructured mesh: %s\n"
                         ,name
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
        log_warn ("Out of memory parseMeshUnstructuredMixedCells\n");
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
                log_warn ("config.xml: invalid variable: %s for cell count\n"
                                 "of mesh: %s\n"
                        ,c, name
                        );
                free (d1);

                return 0;
            }else{
                ccounts_att_nam = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                conca_numb_att_nam(&ccounts_att_nam, name, "ccount", counterstr);
                adios_common_define_attribute (p_new_group,ccounts_att_nam,"/",adios_string,c,"");
                free (ccounts_att_nam);
                counter++;
            }
        }
        else
        {
            cell_list->cell_list.count.var = 0;
            cell_list->cell_list.count.rank = strtod (c, 0);
            counterstr[0] = '\0';
            snprintf(counterstr, 5, "%d", counter);
            ccounts_att_nam = 0;
            conca_numb_att_nam(&ccounts_att_nam, name, "ccount", counterstr);
            adios_common_define_attribute (p_new_group,ccounts_att_nam,"/",adios_double,c,"");
            free (ccounts_att_nam);
            counter++;
        }

        c = strtok (NULL, ",");
    }
    free (d1);

    // We should have at least 2 cell sets, otherwise the cells are uniform
    if (counter <= 1){
       log_warn ("config.xml: Please provide at least 2 cell counts of mesh: %s\n"
                        "or use the 'uniform-cells' tag.\n"
                        ,name
               );
        return 0;
    }

    item->rank = (double) counter;
    conca_att_nam(&ncellset_att_nam, name, "ncsets");
    adios_common_define_attribute (p_new_group,ncellset_att_nam,"/",adios_double,counterstr,"");
    free (ncellset_att_nam);
    item->var = 0;

    mesh->cell_set_count = item;

    // From the number of counts expect the same number of data and type items
    int cell_set_count = counter;
    // Reset counter
    counter = 0;

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
                log_warn ("config.xml: invalid variable: %s for data of mesh: %s\n"
                                 ,c
                                 ,name
                        );
                free (d1);

                return 0;
            } else
            {
                cdata_att_nam = 0;
                counterstr[0] = '\0';
                snprintf(counterstr, 5, "%d", counter);
                conca_numb_att_nam(&cdata_att_nam, name, "cdata", counterstr);
                adios_common_define_attribute (p_new_group,cdata_att_nam,"/",adios_string,c,"");
                free (cdata_att_nam);
                counter++;
            }
        }
        else
        {
            cell_list->cell_list.data = 0;
            log_warn ("config.xml: invalid data: %s of mesh: %s\n"
                             "Please provide a valid variable for cell data\n"
                    ,d1
                    ,name                    );
            free (d1);
            return 0;
        }

        c = strtok (NULL, ",");
    }
    free (d1);

    // If the number of data variables does not match the number of counts
    // Generate an error message
    if (counter != cell_set_count){
       log_warn ("config.xml: Please provide at least %d cell data of mesh: %s\n"
                        "or use the 'uniform-cells' tag\n"
                        ,cell_set_count
                        ,name
               );
        return 0;
    }

    // Reset counter
    counter = 0;

    d1 = strdup (types);
    c = strtok (d1, ",");

    while (c)
    {
        celltype_att_nam = 0;
        counterstr[0] = '\0';
        snprintf(counterstr, 5, "%d", counter);
        conca_numb_att_nam(&celltype_att_nam, name, "ctype", counterstr);

        if (!strcmp(c,"pt") || !strcmp(c,"point"))
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = 1; //strtod (c, 0);
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
        } else
        if (!strcmp(c,"line"))
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = 2; //strtod (c, 0);
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
        } else
        if (!strcmp(c,"tri") || !strcmp(c,"triangle"))
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = 3; //strtod (c, 0);
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
        } else
        if (!strcmp(c,"quad") || !strcmp(c,"quadrilateral"))
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = 4; //strtod (c, 0);
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
        } else
        if (!strcmp(c,"hex") || !strcmp(c,"hexahedron"))
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = 5; //strtod (c, 0);
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
        } else
        if (!strcmp(c,"pri") || !strcmp(c,"prism"))
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = 6; //strtod (c, 0);
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
        } else
        if (!strcmp(c,"tet") || !strcmp(c,"tetrahedron"))
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = 7; //strtod (c, 0);
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
        } else
        if (!strcmp(c,"pyr") || !strcmp(c,"pyramid"))
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = 8; //strtod (c, 0);
            adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
        } else
        if (adios_int_is_var (c))
        {
            cell_list->cell_list.type.var =
                   adios_find_var_by_name (new_group->vars, c
                                          ,new_group->all_unique_var_names
                                          );
            cell_list->cell_list.type.rank = 0;
            if (!cell_list->cell_list.type.var)
            {
                log_warn ("config.xml: invalid var for type variable: %s\n"
                                 "of mesh: %s\n"
                                 ,c
                                 ,name
                        );
                free (d1);

                return 0;
            } else
            {
                adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_string,c,"");
            }
        }
        else
        {
            cell_list->cell_list.type.var = 0;
            cell_list->cell_list.type.rank = strtod (c, 0);
            if (cell_list->cell_list.type.rank > 0 && cell_list->cell_list.type.rank < 9)
            {
                adios_common_define_attribute (p_new_group,celltype_att_nam,"/",adios_double,c,"");
            } else {
                log_warn ("config.xml: invalid type attribute: %s of mesh: %s.\n"
                                 "Please select a number between 1 and 8 for cell types\n"
                                 "or: pt, line, tri, quad, hex, pri, tet or pyr.\n"
                                 ,c
                                 ,name
                        );
                free (celltype_att_nam);
                free (d1);

                return 0;
            }
        }

        c = strtok (NULL, ",");
        counter++;
        free (celltype_att_nam);
    }
    free (d1);

    // If the number of data variables does not match the number of counts
    // Generate an error message
    if (counter != cell_set_count){
       log_warn ("config.xml: Please provide at least %d cell types of mesh: %s\n"
                        "or use the 'uniform-cells' tag\n"
                        ,cell_set_count
                        ,name
               );
        return 0;
    }

    adios_append_mesh_cell_list (&(mesh->cell_list), cell_list);

    return 1;
}

// primary mesh XML parsing
static int parseMeshUniform (mxml_node_t * node
                            ,struct adios_group_struct * new_group
                            ,struct adios_mesh_uniform_struct ** mesh
                            ,const char * name
                            )
{
    mxml_node_t * n;
    int saw_dimensions = 0;
    int saw_origin = 0;
    int saw_spacing = 0;
    int saw_maximum = 0;

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
                log_warn ("config.xml: only one dimensions definition "
                                 "allowed per mesh sructured-points (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_dimensions = 1;
            dimensions = mxmlElementGetAttr (n, "value");

            if (!dimensions)
            {
                log_warn ("config.xml: value attribute on "
                                 "dimensions required (%s)\n"
                                 ,name
                        );
                return 0;
            }

            if (!parseMeshUniformDimensions (dimensions, new_group, *mesh, name))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "origin"))
        {
            const char * value;

            if (saw_origin)
            {
                log_warn ("config.xml: only one origin definition "
                                 "allowed per mesh uniform (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_origin = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "origin required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUniformOrigin (value, new_group, *mesh, name))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "spacing"))
        {
            const char * value;

            if (saw_spacing)
            {
                log_warn ("config.xml: only one spacing "
                                 "definition allowed per mesh uniform (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_spacing = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "spacing required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUniformSpacing (value, new_group, *mesh, name))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "maximum"))
        {
            const char * value;

            if (saw_maximum)
            {
                log_warn ("config.xml: only one maximum "
                                 "definition allowed per mesh uniform (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_maximum = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "max required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUniformMaxima (value, new_group, *mesh, name))
                return 0;
        } else
        {
            if (!strncmp (n->value.element.name, "!--", 3)) // a comment
            {
                continue;
            }
        }
    }

    // If nothing is given, simply assume basic uniform plots using
    // the dimensions of the variable, origin=0 and spacing=1
    /*if (!saw_dimensions)
    {
        log_debug ("config.xml: dimensions required on mesh "
                   "type=uniform (%s)\n",name);
        return 0;
    }*/

    return 1;
}

static int parseMeshRectilinear (mxml_node_t * node
                               ,struct adios_group_struct * new_group
                               ,struct adios_mesh_rectilinear_struct ** mesh
                               ,const char * name
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
                log_warn ("config.xml: only one dimensions "
                                 "definition allowed per mesh rectilinear (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_dimensions = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "dimensions required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshRectilinearDimensions (value, new_group, *mesh, name))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "coordinates-multi-var"))
        {
            const char * value;

            if (saw_coordinates_multi_var || saw_coordinates_single_var)
            {
                log_warn ("config.xml: only one coordinates "
                                 "definition allowed per mesh rectilinear (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_coordinates_multi_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "coordinates-multi-var required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshRectilinearCoordinatesMultiVar (value, new_group, *mesh, name))
                return 0;
            (*mesh)->coordinates_single_var = adios_flag_no;
        } else
        if (!strcasecmp (n->value.element.name, "coordinates-single-var"))
        {
            const char * value;

            if (saw_coordinates_single_var || saw_coordinates_multi_var)
            {
                log_warn ("config.xml: only one coordinates "
                                 "definition allowed per mesh rectilinear (%s)\n"
                                 ,name
                        );
                return 0;
            }

            saw_coordinates_single_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "coordinates-single-var required (%s)\n"
                                 ,name
                        );
                return 0;
            }

            if (!parseMeshRectilinearCoordinatesSingleVar (value, new_group, *mesh, name))
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
        log_warn ("config.xml: dimensions required on mesh "
                         "type=rectilinear (%s)\n"
                         ,name
                );
        return 0;
    }
    if (!saw_coordinates_multi_var && !saw_coordinates_single_var)
    {
        log_warn ("config.xml: coordinates-multi-var or "
                         "coordinates-single-var required on mesh "
                         "type=rectilinear (%s)\n"
                         ,name
                );

        return 0;
    }

    return 1;
}

static int parseMeshStructured (mxml_node_t * node
                               ,struct adios_group_struct * new_group
                               ,struct adios_mesh_structured_struct ** mesh
                               ,const char * name
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
                log_warn ("config.xml: only one nspace "
                                 "definition allowed per mesh structured (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_nspace = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "nspace required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshStructuredNspace (value, new_group, *mesh, name))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "dimensions"))
        {
            const char * value;

            if (saw_dimensions)
            {
                log_warn ("config.xml: only one dimensions "
                                 "definition allowed per mesh structured (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_dimensions = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "dimensions required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshStructuredDimensions (value, new_group, *mesh, name))
                return 0;
        } else
        if (!strcasecmp (n->value.element.name, "points-multi-var"))
        {
            const char * value;

            if (saw_points_multi_var || saw_points_single_var)
            {
                log_warn ("config.xml: only one points "
                                 "definition allowed per mesh structured (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_points_multi_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "points-multi-var required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshStructuredPointsMultiVar (value, new_group, *mesh, name))
                return 0;
            (*mesh)->points_single_var = adios_flag_no;
        } else
        if (!strcasecmp (n->value.element.name, "points-single-var"))
        {
            const char * value;

            if (saw_points_multi_var || saw_points_single_var)
            {
                log_warn ("config.xml: only one points "
                                 "definition allowed per mesh structured (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_points_single_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "points-single-var required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshStructuredPointsSingleVar (value, new_group, *mesh, name))
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
        log_warn ("config.xml: dimensions required on mesh "
                         "type=structured (%s)\n"
                         ,name
                );

        return 0;
    }
    if (!saw_points_multi_var && !saw_points_single_var)
    {
        log_warn ("config.xml: points-single-var or points-multi-var "
                         "required on mesh type=structured (%s)\n"
                         ,name
                );

        return 0;
    }
    // nspace should not be required since the single var could be multi-dim
    /*if (saw_points_single_var && !saw_nspace)
    {
        log_warn ("config.xml: points-single-var and nspace "
                         "required on mesh type=structured (%s)\n"
                         ,name
                );

        return 0;
    }*/

    return 1;
}

static int parseMeshUnstructured (mxml_node_t * node
                                 ,struct adios_group_struct * new_group
                                 ,struct adios_mesh_unstructured_struct ** mesh
                                 ,const char * name
                                 )
{
    mxml_node_t * n;
    int saw_points = 0;
    int saw_nspace =0;
    int saw_number_of_points = 0;
    int saw_points_multi_var = 0;
    int saw_points_single_var = 0;
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

        if (!strcasecmp (n->value.element.name, "nspace"))
        {
            const char * value;

            if (saw_nspace)
            {
                log_warn ("config.xml: only one nspace "
                                 "definition allowed per mesh structured (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_nspace = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "nspace required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUnstructuredNspace (value, new_group, *mesh, name))
                return 0;
        }else
        if (!strcasecmp (n->value.element.name, "number-of-points"))
        {
            const char * value;

            if (saw_number_of_points)
            {
                log_warn ("config.xml: only one number-of-points "
                                 "definition allowed per mesh structured (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_number_of_points = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "number-of-points required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUnstructuredNpoints (value, new_group, *mesh, name))
                return 0;
        }else
        if (!strcasecmp (n->value.element.name, "points-multi-var"))
        {
            const char * value;

            if (saw_points_multi_var || saw_points_single_var)
            {
                log_warn ("config.xml: only one points "
                                 "definition allowed per mesh unstructured (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_points_multi_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "points-multi-var required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUnstructuredPointsMultiVar (value, new_group, *mesh, name))
                return 0;
            (*mesh)->points_single_var = adios_flag_no;
        } else
        if (!strcasecmp (n->value.element.name, "points-single-var"))
        {
            const char * value;

            if (saw_points_multi_var || saw_points_single_var)
            {
                log_warn ("config.xml: only one points "
                                 "definition allowed per mesh unstructured (%s)\n"
                                 ,name
                        );

                return 0;
            }

            saw_points_single_var = 1;
            value = mxmlElementGetAttr (n, "value");

            if (!value)
            {
                log_warn ("config.xml: value attribute on "
                                 "points-single-var required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUnstructuredPointsSingleVar (value, new_group, *mesh, name))
                return 0;
            (*mesh)->points_single_var = adios_flag_yes;
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
                log_warn ("config.xml: count attribute on "
                                 "uniform-cells required (%s)\n"
                                 ,name
                        );

                return 0;
            }
            if (!data)
            {
                log_warn ("config.xml: data attribute on "
                                 "uniform-cells required (%s)\n"
                                 ,name
                        );

                return 0;
            }
            if (!type)
            {
                log_warn ("config.xml: type attribute on "
                                 "uniform-cells required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUnstructuredUniformCells (count, data, type
                                                   ,new_group, *mesh
                                                   ,name
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
            types = mxmlElementGetAttr (n, "type");

            if (!count)
            {
                log_warn ("config.xml: count attribute on "
                                 "mixed-cells required (%s)\n"
                                 ,name
                        );

                return 0;
            }
            if (!data)
            {
                log_warn ("config.xml: data attribute on "
                                 "mixed-cells required (%s)\n"
                                 ,name
                        );

                return 0;
            }
            if (!types)
            {
                log_warn ("config.xml: types attribute on "
                                 "mixed-cells required (%s)\n"
                                 ,name
                        );

                return 0;
            }

            if (!parseMeshUnstructuredMixedCells (count, data, types
                                                   ,new_group, *mesh
                                                   ,name
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

    if (!saw_points_multi_var && !saw_points_single_var)
    {
        log_warn ("config.xml: points-single-var or points-multi-var "
                         "required on mesh type=unstructured (%s)\n"
                         ,name
                );

        return 0;
    }
    // nspace should not be required since the single variable could be multi-dim
    /*if (saw_points_single_var && !saw_nspace && !saw_number_of_points)
    {
        log_warn ("config.xml: with points-single-var, nspace or number-of-points "
                         "required on mesh type=unstructured (%s)\n"
                         ,name
                );

        return 0;

    }*/

    if (!saw_cell_set)
    {
        log_warn ("config.xml: at least one cell-set required on "
                         "mesh type=unstructured (%s)\n"
                         ,name
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
        // Although this is not used, we are leaving in the retrevial
        // of this to avoid messages from all of the existing XML files.
        // In a few months, once everything has been updated, we can remove
        // this code
        GET_ATTR("coordination-communicator",attr,coordination_comm,"adios-group")
        GET_ATTR("coordination-var",attr,coordination_var,"adios-group")
        GET_ATTR("host-language",attr,host_language,"adios-group")
        GET_ATTR("time-index",attr,time_index_name,"adios-group")
        GET_ATTR("stats",attr,stats,"adios-group")
        log_warn ("config.xml: unknown attribute '%s' on %s "
                         "(ignored)\n"
                ,attr->name
                ,"adios-group"
                );
    }

    if (!datagroup_name)
    {
        adios_error(err_invalid_attrname,"config.xml: name attribute required on adios-group\n");

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
                adios_error (err_invalid_host_language, "config.xml: invalid host-language %s"
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
            log_error ("config.xml, invalid stats %s"
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
            const char * mesh = 0;
            const char * type = 0;
            const char * dimensions = 0;
            const char * dimension = 0;
            const char * gread = 0;
            const char * gwrite = 0;
            const char * read_flag = 0;
            const char * transform_type = 0; // NCSU ALACRITY-ADIOS
            enum ADIOS_DATATYPES t1;
            char  * mpath = 0;

            for (i = 0; i < n->value.element.num_attrs; i++)
            {
                mxml_attr_t * attr = &n->value.element.attrs [i];

                GET_ATTR("name",attr,name,"var")
                GET_ATTR("mesh",attr,mesh,"var")
                GET_ATTR("path",attr,path,"var")
                GET_ATTR("type",attr,type,"var")
                GET_ATTR("dimensions",attr,dimensions,"var")
                GET_ATTR("dimension",attr,dimension,"var")
                GET_ATTR("gwrite",attr,gwrite,"var")
                GET_ATTR("gread",attr,gread,"var")
                GET_ATTR("read",attr,read_flag,"var")
                GET_ATTR("transform",attr,transform_type,"var") // NCSU ALACRITY-ADIOS
                log_warn ("config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"var"
                        );
            }

            if (!name)
                name = ""; // this will catch the error
            if (!path)
                path = "/";
            if (!type)
                type = ""; // this will catch the error
            if (!mesh)
                mesh = "";

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
                                         ,transform_type // NCSU ALACRITY-ADIOS
                                         )
               )
            {
                return 0;
            }else{
                // Successfully define a variable, so now
                // an attribute for the mesh if it exists.
                if (strcmp(mesh,"")){
                    mpath = malloc(strlen("/adios_schema")+strlen(name)+1);
                    strcpy(mpath,name);
                    strcat(mpath,"/adios_schema");
                    adios_common_define_attribute (ptr_new_group,mpath,path,adios_string,mesh,"");
                }
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

                log_warn ("config.xml: unknown attribute '%s' on %s "
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
                    adios_error (err_global_dim_required, "config.xml: dimensions required on "
                                     "global-bounds\n"
                            );

                    return 0;
                }
            }
            if (!offsets)
            {
                offsets = (offset ? offset : local_offsets);
                offsets = (offsets ? offsets : local_offset);
                adios_error (err_global_offset_required, "config.xml: offsets required on "
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
                    const char * mesh = 0;
                    const char * path = 0;
                    const char * type = 0;
                    const char * dimension = 0;
                    const char * dimensions = 0;
                    const char * gwrite = 0;
                    const char * gread = 0;
                    const char * read_flag = 0;
                    const char * transform_type = 0; // NCSU ALACRITY-ADIOS
                    enum ADIOS_DATATYPES t1;
                    char * mpath = 0;

                    for (i = 0; i < n1->value.element.num_attrs; i++)
                    {
                        mxml_attr_t * attr = &n1->value.element.attrs [i];

                        GET_ATTR("name",attr,name,"var")
                        GET_ATTR("mesh",attr,mesh,"var")
                        GET_ATTR("path",attr,path,"var")
                        GET_ATTR("type",attr,type,"global-bounds var")
                        GET_ATTR("dimensions",attr,dimensions,"var")
                        GET_ATTR("dimension",attr,dimension,"var")
                        GET_ATTR("gwrite",attr,gwrite,"var")
                        GET_ATTR("gread",attr,gread,"var")
                        GET_ATTR("read",attr,read_flag,"var")
                        GET_ATTR("transform",attr,transform_type,"var") // NCSU ALACRITY-ADIOS
                        log_warn ("config.xml: unknown attribute '%s' "
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
                    if (!mesh)
                        mesh = "";

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
                                                 ,transform_type // NCSU ALACRITY-ADIOS
                                                 )
                       )
                    {
                        return 0;
                    }else{
                        // Successfully define a variable, so now
                        // an attribute for the mesh if it exists.
                        if (strcmp(mesh,"")){
                            mpath = malloc(strlen("/adios_schema")+strlen(name)+1);
                            strcpy(mpath,name);
                            strcat(mpath,"/adios_schema");
                            adios_common_define_attribute (ptr_new_group,mpath,path,adios_string,mesh,"");
                         }
                    }
                } else
                {
                    if (!strncmp (n1->value.element.name, "!--", 3)) // comment
                    {
                        continue;
                    }
                    else
                    {
                        log_warn ("config.xml: invalid xml element: "
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
                log_warn ("config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"attribute"
                        );
            }

            if (!name)
            {
                log_warn ("config.xml: attribute element requires "
                                 "name\n");

                return 0;
            }
            if (!path)
            {
                log_warn ("config.xml: attribute element requires "
                                 "path\n");

                return 0;
            }
            if ((!value && !var) || (value && var))
            {
                log_warn ("config.xml: attribute element '%s' "
                                 "requires either value OR var\n"
                        ,name
                        );

                return 0;
            }
            if (var && type)
            {
                log_warn ("config.xml: attribute element '%s'. "
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
            int t_varying;
            const char * name;

            // Get the mesh name
            name = mxmlElementGetAttr (n, "name");
            // Get the mesh type
            type = mxmlElementGetAttr (n, "type");
            // Get the time varying parameter
            time_varying = mxmlElementGetAttr(n, "time-varying");

            if (!type)
                type = "";

            if (!strcmp(time_varying,"yes")){
                t_varying = adios_flag_yes;
            }else if (!strcmp(time_varying,"no")){
                t_varying = adios_flag_no;
            }else{
                t_varying = adios_flag_no;
                // If the user enters anything else than "yes" or "no"
                // Output a warning letting them no that the default ("no"
                // will be use give instead of their value 
                log_warn ("config.xml: the value of the time varying "
                        "attribute can only be 'yes' or 'no'. The "
                        "unrecognize value of '%s' is ignored and "
                        "replaced by 'no'."
                        "\n"
                        ,time_varying
                        );
            }

            char * meshtype = 0;
            char * meshtime = 0;

            conca_att_nam(&meshtype, name, "type");
            conca_att_nam(&meshtime, name, "time");

            if (!strcasecmp (type, "uniform"))
            {
                struct adios_mesh_struct * mes;
                mes = adios_common_define_mesh(ptr_new_group, name,
                                               t_varying, ADIOS_MESH_UNIFORM);

                if (mes) {
                    // Define attribute for the type and time varying characteristics
                    adios_common_define_attribute (ptr_new_group,meshtype,"/",adios_string,type,"");
                    adios_common_define_attribute (ptr_new_group,meshtime,"/",adios_string,time_varying,"");
                    free (meshtype);
                    free (meshtime);

                    // Parse the uniform mesh tags
                    parseMeshUniform (n, new_group, &mes->uniform, name);
                }
            } else
            if (!strcasecmp (type, "structured"))
            {
                struct adios_mesh_struct * mes;
                mes = adios_common_define_mesh(ptr_new_group, name,
                                               t_varying, ADIOS_MESH_STRUCTURED);

                if (mes) {
                    // Define attribute for the type and time varying characteristics
                    adios_common_define_attribute (ptr_new_group,meshtype,"/",adios_string,type,"");
                    adios_common_define_attribute (ptr_new_group,meshtime,"/",adios_string,time_varying,"");
                    free (meshtype);
                    free (meshtime);

                    // Parse the uniform mesh tags
                    parseMeshStructured (n, new_group, &mes->structured, name);
                }
            } else
            if (!strcasecmp (type, "rectilinear"))
            {
                struct adios_mesh_struct * mes;
                mes = adios_common_define_mesh(ptr_new_group, name,
                                               t_varying, ADIOS_MESH_RECTILINEAR);

                if (mes) {
                    // Define attribute for the type and time varying characteristics
                    adios_common_define_attribute (ptr_new_group,meshtype,"/",adios_string,type,"");
                    adios_common_define_attribute (ptr_new_group,meshtime,"/",adios_string,time_varying,"");
                    free (meshtype);
                    free (meshtime);

                    // Parse the uniform mesh tags
                    parseMeshRectilinear (n, new_group, &mes->rectilinear, name);
                }
            } else
            if (!strcasecmp (type, "unstructured"))
            {
                struct adios_mesh_struct * mes;
                mes = adios_common_define_mesh(ptr_new_group, name,
                                               t_varying, ADIOS_MESH_UNSTRUCTURED);

                if (mes) {
                    // Define attribute for the type and time varying characteristics
                    adios_common_define_attribute (ptr_new_group,meshtype,"/",adios_string,type,"");
                    adios_common_define_attribute (ptr_new_group,meshtime,"/",adios_string,time_varying,"");
                    free (meshtype);
                    free (meshtime);

                    // Parse the uniform mesh tags
                    parseMeshUnstructured (n, new_group, &mes->unstructured, name);
                }
            } else
            {
                log_warn ("config.xml: invalid mesh type: '%s'\n"
                        ,type
                        );

                return 0;
            }
        } else
        if (!strcasecmp (n->value.element.name, "gwrite"))
        {
            const char * src = 0;

            for (i = 0; i < n->value.element.num_attrs; i++)
            {
                mxml_attr_t * attr = &n->value.element.attrs [i];

                GET_ATTR("src",attr,src,"var")
                log_warn ("config.xml: unknown attribute '%s' on %s "
                                 "(ignored)\n"
                        ,attr->name
                        ,"gwrite"
                        );
            }
            if (!src)
            {
                log_warn ("config.xml: gwrite element requires "
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
                log_warn ("config.xml: invalid xml element: '%s'\n"
                        ,n->value.element.name
                        );

                return 0;
            }
        }
    }

    // now that we have declared the whole group, validate that the
    // paths specified in attributes refer to real things or give
    // a warning
    /* // This is not required in ADIOS so why do this warnings?
    struct adios_attribute_struct * a = new_group->attributes;
    while (a)
    {
        if (!validatePath (new_group->vars, a->path))
        {
             log_warn ("config.xml warning: attribute element '%s' "
                       "has path '%s' that does not match "
                       "any var path or name.\n",a->name, a->path);
        }

        a = a->next;
    }
    */

    return 1;
}

// concat numbered attribute name strings
void conca_numb_att_nam(char ** returnstr, const char * meshname, char * att_nam, char counterstr[5]) {
    *returnstr = malloc (strlen("adios_schema/") + strlen(meshname) + strlen(att_nam) + strlen(counterstr) + 3);
    strcpy(*returnstr,"adios_schema");
    strcat(*returnstr,"/");
    strcat(*returnstr,meshname);
    strcat(*returnstr,"/");
    strcat(*returnstr,att_nam);
    strcat(*returnstr,counterstr);
}

// concat attribute name strings
void conca_att_nam(char ** returnstr, const char * meshname, char * att_nam) {
    int slength = 0;
    slength = strlen("adios_schema/");
    slength = slength + strlen(meshname);
    slength = slength + 1;
    slength = slength + 1;
    slength = slength + strlen(att_nam);

    *returnstr = malloc (slength);

    strcpy(*returnstr,"adios_schema/");
    strcat(*returnstr,meshname);
    strcat(*returnstr,"/");
    strcat(*returnstr,att_nam);
}

// Append a mesh to a group
enum ADIOS_FLAG adios_append_mesh (struct adios_mesh_struct ** root
                                 ,struct adios_mesh_struct * mesh
                                 ,uint16_t id
                                 )
{
    while (root)
    {
        if (*root && !strcasecmp ((*root)->name, mesh->name))
        {
            return adios_flag_no;
        }
        if (!*root)
        {
            *root = mesh;
            root = 0;
        }
        else
        {
            root = &(*root)->next;
        }
    }

    return adios_flag_yes;
}

// Define a new mesh
struct adios_mesh_struct * adios_common_define_mesh (
        int64_t group_id, const char * name,
        enum ADIOS_FLAG time_varying,
        enum ADIOS_MESH_TYPE type)
{
    struct adios_group_struct * t = (struct adios_group_struct *) group_id;
    struct adios_mesh_struct * m = (struct adios_mesh_struct *)
                               malloc (sizeof (struct adios_mesh_struct));
    enum ADIOS_FLAG flag;

    m->name = strdup (name);
    m->type = type;
    m->time_varying = time_varying;
    m->next = 0;

    if (type == 1){
        m->uniform = (struct adios_mesh_uniform_struct *) calloc (1, sizeof (struct adios_mesh_uniform_struct));
    }else if (type == 3){
        m->rectilinear = (struct adios_mesh_rectilinear_struct *) calloc (1, sizeof (struct adios_mesh_rectilinear_struct));
    }else if (type == 2){
        m->structured = (struct adios_mesh_structured_struct *) calloc (1, sizeof (struct adios_mesh_structured_struct));
    }else if (type == 4){
        m->unstructured = (struct adios_mesh_unstructured_struct *) calloc (1, sizeof (struct adios_mesh_unstructured_struct));
    }else{
        m->uniform = (struct adios_mesh_uniform_struct *) calloc (1, sizeof (struct adios_mesh_uniform_struct));
    }
    flag = adios_append_mesh (&t->meshs, m, t->mesh_count);
    if (flag == adios_flag_no)
    {
        log_warn ("config.xml: unique mesh names required; "
                         "second mesh: %s will be ignored.\n"
                         ,name
                );
        free(m);
        m = 0;
    } else {
        t->mesh_count++;
    }

    return m;
}

// Find a mesh
struct adios_mesh_struct * adios_find_mesh_by_name (struct adios_mesh_struct * root
                                                 ,const char * name
                                                 )
{
    int done = 0;
    struct adios_mesh_struct * mesh = 0;

    if (!name)
        done = 1;

    while (!done && root)
    {
        if (!strcasecmp (name, root->name))
        {
            done = 1;
            mesh = root;
        }
        else
        {
            root = root->next;
        }
    }

    return mesh;
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
        log_warn ("config.xml: unknown attribute '%s' on %s "
                         "(ignored)\n"
                ,attr->name
                ,"method"
                );
    }

    if (!var)
    {
        log_warn ("config.xml: variable name must be given\n");
        return 0;
    }

    if (!group)
    {
        log_warn ("config.xml: adios-group name must be given\n");
        return 0;
    }

    adios_common_get_group (&group_id, group);
    g = (struct adios_group_struct *) group_id;

    if (!g)
    {
        log_warn ("config.xml: Didn't find group %s for analysis\n", group);
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
        log_warn ("config.xml: unknown attribute '%s' on %s "
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
            adios_error (err_invalid_method, "config.xml: method %s for group %s base-path "
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
        log_warn ("config.xml: unknown attribute '%s' on %s "
                         "(ignored)\n"
                ,attr->name
                ,"buffer"
                );
    }



    if ((!size_MB && !free_memory_percentage) || !allocate_time)
    {
        adios_error (err_invalid_buffer_size, "config.xml: must define allocate-time and either "
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
                adios_error (err_invalid_buffer_size, "config.xml: buffer allocate-time %s "
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
                adios_error (err_invalid_buffer_size, "config.xml: buffer size-MB is either 0 or "
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
                adios_error (err_invalid_buffer_size, "config.xml: buffer free-memory-percentage %s "
                                 "is not an integer between 1 and 100\n"
                        ,free_memory_percentage
                        );

                return 0;
            }
        }

        if (adios_buffer_alloc_when_get() == ADIOS_BUFFER_ALLOC_NOW)
        {

// Do not attempt to allocate the buffer when this is being called from adios_lint
#ifndef _INTERNAL
            return adios_set_buffer_size ();
#endif

        }
    }

    return 1;
}

int adios_parse_config (const char * config, MPI_Comm comm)
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
    // NCSU ALACRITY-ADIOS - Initialize transform methods
    adios_transform_init();

    char * buffer = NULL;
//#if HAVE_MPI
    int buffer_size = 0;
    int rank;
    MPI_Comm_rank (comm, &rank);
    if (rank == 0)
    {
//#endif
        fp = fopen (config, "r");
        if (!fp)
        {
            adios_error (err_missing_config_file, "missing config file %s\n", config);

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
                adios_error (err_expected_read_size_mismatch, "error reading config file: %s. Expected %d Got %d\n"
                        ,config, s.st_size, bytes_read );

                return 0;
            }
        }
        else
        {
            adios_error (err_allocating_buffer_size, "error allocating %d for reading config.\n"
                    ,s.st_size + 1
                    );

            return 0;
        }
        fclose (fp);
//#if HAVE_MPI
        buffer_size = s.st_size;
        MPI_Bcast (&buffer_size, 1, MPI_INT, 0, comm);
        MPI_Bcast (buffer, buffer_size, MPI_BYTE, 0, comm);
    }
    else
    {
        MPI_Bcast (&buffer_size, 1, MPI_INT, 0, comm);
        buffer = malloc (buffer_size + 1);
        if (!buffer)
        {
            adios_error (err_allocating_buffer_size, "cannot allocate %d bytes to receive config file\n"
                    ,buffer_size + 1
                    );

            return 0;
        }
        MPI_Bcast (buffer, buffer_size, MPI_BYTE, 0, comm);
        buffer [buffer_size] = 0;
    }
//#endif

    doc = mxmlLoadString (NULL, buffer, MXML_TEXT_CALLBACK);
    free (buffer);
    buffer = NULL;

    if (!doc)
    {
        adios_error (err_invalid_xml_doc, "config.xml: unknown error parsing XML "
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
            adios_error (err_invalid_xml_doc, "config.xml: invalid root xml element: %s\n"
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
        adios_error (err_invalid_xml_doc, "config.xml: invalid root xml element: %s\n"
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
            log_warn ("config.xml: unknown attribute '%s' on %s "
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
                adios_error (err_invalid_host_language, "config.xml: invalid host-language %s"
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
                	        log_warn ("config.xml: invalid element: %s\n"
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
        adios_error (err_no_group_defined, "config.xml: must define at least 1 adios-group in "
                         "config.xml\n"
                );

        return 0;
    }
    if (!saw_method)
    {
        adios_error (err_no_method_defined, "config.xml: must define at least 1 method for "
                         "the adios-group in config.xml\n"
                );

        return 0;
    }
    if (!saw_buffer)
    {
        adios_error (err_no_buffer_defined, "config.xml: must define the buffer element in "
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
    // NCSU ALACRITY-ADIOS - Initialize transform methods
    adios_transform_init();

    return 1;
}

static PairStruct * get_and_preprocess_params (const char * parameters)
{
    PairStruct *params, *p, *prev_p;
    int verbose_level, removeit, save;
    char *end;

    params = text_to_name_value_pairs (parameters);

    /*
    p = params;
    while (p) {
        fprintf(stderr, "-------  Param    name = %s  value = %s\n", p->name, p->value);
        p = p->next;
    }
    */
    prev_p = NULL;
    p = params;
    while (p) {
        /*fprintf(stderr, "Parameter    name = %s  value = %s\n", p->name, p->value);*/
        removeit = 0;
        if (!strcasecmp (p->name, "verbose"))
        {
            if (p->value) {
                errno = 0;
                verbose_level = strtol(p->value, &end, 10);
                if (errno || (end != 0 && *end != '\0')) {
                    log_error ("Invalid 'verbose' parameter passed to read init function: '%s'\n", p->value);
                    verbose_level = 1; // print errors only
                }
            } else {
                verbose_level = 3;  // info level
            }
            adios_verbose_level = verbose_level;
            removeit = 1;
        }
        else if (!strcasecmp (p->name, "quiet"))
        {
            adios_verbose_level = 0; //don't print errors
            removeit = 1;
        }
        else if (!strcasecmp (p->name, "logfile"))
        {
            /*fprintf (stderr,"****************** logfile = %s\n", p->value);*/
            if (p->value) {
                adios_logger_open (p->value, -1);
            }
            removeit = 1;
        }
        else if (!strcasecmp (p->name, "abort_on_error"))
        {
            adios_abort_on_error = 1;
            save = adios_verbose_level;
            adios_verbose_level = 2;
            log_warn ("ADIOS is set to abort on error\n");
            adios_verbose_level = save;
            removeit = 1;
        }
        if (removeit) {
            if (p == params) {
                // remove head
                //fprintf(stderr, "  Remove HEAD  p = %x p->next = %x\n", p, p->next);
                p = p->next;
                params->next = NULL;
                free_name_value_pairs (params);
                params = p;
            } else {
                // remove from middle of the list
                //fprintf(stderr, "  Remove MIDDLE prev = %x p = %x p->next = %x\n", prev_p, p, p->next);
                prev_p->next = p->next;
                p->next = NULL;
                free_name_value_pairs (p);
                p = prev_p->next;
            }
        } else {
            //fprintf(stderr, "  Keep MIDDLE prev = %x p = %x p->next = %x\n", prev_p, p, p->next);
            prev_p = p;
            p = p->next;
        }
    }

    return params;
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
    new_method->parameters = strdup (parameters); // string goes into BP file
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
            PairStruct * params = get_and_preprocess_params (parameters);

            adios_transports [new_method->m].adios_init_fn
                                       (params, new_method);

            free_name_value_pairs (params);
        }
    }
    else
    {
        adios_error (err_invalid_write_method, "config.xml: invalid transport: %s\n", method);

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
        adios_error (err_missing_invalid_group, "config.xml: Didn't find group: %s for transport: %s\n"
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
            PairStruct * params = get_and_preprocess_params (parameters);

            adios_transports [new_method->m].adios_init_fn
                                       (params, new_method);

            free_name_value_pairs (params);
        }
    }
    else
    {
        adios_error (err_invalid_write_method, "config.xml: invalid transport: %s\n", method);

        free (new_method->base_path);
        free (new_method->method);
        free (new_method->parameters);
        free (new_method);

        return 0;
    }

    g = (struct adios_group_struct *) group_id;
    if (!g)
    {
        adios_error (err_missing_invalid_group, "config.xml: invalid group id: %llu for transport: %s\n"
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
            adios_error (err_group_method_mismatch, "config.xml: method %s for group %s.  Group does "
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
            // ADIOS Schema
            // if (adios_groups->group->vars->mesh)
                // free (adios_groups->group->vars->mesh);


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

            // NCSU ALACRITY-ADIOS - Clean transform metadata
            adios_transform_clear_transform_var(adios_groups->group->vars);

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

        // ADIOS Schema
        /*while (adios_groups->group->meshs)
        {
            struct adios_mesh_struct * meshs = adios_groups->group->meshs->next;

            if (adios_groups->group->meshs->name)
                free (adios_groups->group->meshs->name);

            switch (adios_groups->group->meshs->type)
            {
                case ADIOS_MESH_UNIFORM:
                {
                    struct adios_mesh_item_list_struct * i;
                    while (adios_groups->group->meshs->uniform->dimensions)
                    {
                        i = adios_groups->group->meshs->uniform->dimensions->next;
                        free (adios_groups->group->meshs->uniform->dimensions);
                        adios_groups->group->meshs->uniform->dimensions = i;
                    }
                    while (adios_groups->group->meshs->uniform->origin)
                    {
                        i = adios_groups->group->meshs->uniform->origin->next;
                        free (adios_groups->group->meshs->uniform->origin);
                        adios_groups->group->meshs->uniform->origin = i;
                    }
                    while (adios_groups->group->meshs->uniform->spacing)
                    {
                        i = adios_groups->group->meshs->uniform->spacing->next;
                        free (adios_groups->group->meshs->uniform->spacing);
                        adios_groups->group->meshs->uniform->spacing = i;
                    }
                    while (adios_groups->group->meshs->uniform->origin)
                    {
                        i = adios_groups->group->meshs->uniform->origin->next;
                        free (adios_groups->group->meshs->uniform->origin);
                        adios_groups->group->meshs->uniform->origin = i;
                    }
                    break;
                }
                case ADIOS_MESH_STRUCTURED:
                {
                    struct adios_mesh_item_list_struct * i;
                    struct adios_mesh_var_list_struct * v;
                    while (adios_groups->group->meshs->structured->dimensions)
                    {
                        i = adios_groups->group->meshs->structured->dimensions->next;
                        free (adios_groups->group->meshs->structured->dimensions);
                        adios_groups->group->meshs->structured->dimensions = i;
                    }
                    while (adios_groups->group->meshs->structured->points)
                    {
                        v = adios_groups->group->meshs->structured->points->next;
                        free (adios_groups->group->meshs->structured->points);
                        adios_groups->group->meshs->structured->points = v;
                    }
                    if (adios_groups->group->meshs->structured->nspace)
                        free (adios_groups->group->meshs->structured->nspace);

                    break;
                }
                case ADIOS_MESH_RECTILINEAR:
                {
                    struct adios_mesh_item_list_struct * i;
                    struct adios_mesh_var_list_struct * v;
                    while (adios_groups->group->meshs->rectilinear->dimensions)
                    {
                        i = adios_groups->group->meshs->rectilinear->dimensions->next;
                        free (adios_groups->group->meshs->rectilinear->dimensions);
                        adios_groups->group->meshs->rectilinear->dimensions = i;
                    }
                    while (adios_groups->group->meshs->rectilinear->coordinates)
                    {
                        v = adios_groups->group->meshs->rectilinear->coordinates->next;
                        free (adios_groups->group->meshs->rectilinear->coordinates);
                        adios_groups->group->meshs->rectilinear->coordinates = v;
                    }

                    break;
                }
                case ADIOS_MESH_UNSTRUCTURED:
                {
                    struct adios_mesh_var_list_struct * v;
                    if (adios_groups->group->meshs->unstructured->cell_set_count)
                        free (adios_groups->group->meshs->unstructured->cell_set_count);
                    if (adios_groups->group->meshs->unstructured->points_count)
                        free (adios_groups->group->meshs->unstructured->points_count);
                    if (adios_groups->group->meshs->unstructured->nspace)
                        free (adios_groups->group->meshs->unstructured->nspace);
                    while (adios_groups->group->meshs->unstructured->points)
                    {
                        v = adios_groups->group->meshs->unstructured->points->next;
                        free (adios_groups->group->meshs->unstructured->points);
                        adios_groups->group->meshs->unstructured->points = v;
                    }
                    while (adios_groups->group->meshs->unstructured->cell_list)
                    {
                        struct adios_mesh_cell_list_list_struct * next
                          = adios_groups->group->meshs->unstructured->cell_list->next;
                        free (adios_groups->group->meshs->unstructured->cell_list);
                        adios_groups->group->meshs->unstructured->cell_list = next;
                    }

                    break;
                }
            }

            free (adios_groups->group->meshs);
            adios_groups->group->meshs = meshs;
        }*/

        free (adios_groups->group);
        free (adios_groups);
        adios_groups = groups;
    }
}


