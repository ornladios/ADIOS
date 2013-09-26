/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "config.h"
#include "adios.h"
#include "adios_types.h"
#include "adios_bp_v1.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"
#ifdef HAVE_NSSI
#include "nssi_client.h"
#include "adios_nssi_args.h"
#include "adios_nssi_config.h"
#include "nssi_logger.h"
#endif

#include "io_timer.h"


#define NUM_GP 24
void adios_nssi_filter_get_write_buffer(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        uint64_t *size,
        void **buffer,
        struct adios_method_struct *method)
{
}

///////////////////////////
// Datatypes
///////////////////////////
struct adios_nssi_filter_data_struct
{
    char *sm_method;
    char *sm_parameters;

    int      fd;
    MPI_Comm group_comm;
    int      rank;
    int      size;

    struct adios_method_struct *submethod;
};

///* Need a struct to encapsulate var offset info.
// */
//struct var_offset {
//    char      opath[ADIOS_PATH_MAX];
//    char      oname[ADIOS_PATH_MAX];
//    void     *ovalue;
//    uint64_t  osize;
//};
//
///* list of variable offsets */
//List var_offset_list;
//
///* Need a struct to encapsulate var dim info.
// */
//struct var_dim {
//    char      dpath[ADIOS_PATH_MAX];
//    char      dname[ADIOS_PATH_MAX];
//    void     *dvalue;
//    uint64_t  dsize;
//};
//
///* list of variable offsets */
//List var_dim_list;

/* Need a struct to encapsulate anonymous dim info.
 */
struct anonymous_dim {
    char      adpath[ADIOS_PATH_MAX];
    char      adpathname[ADIOS_PATH_MAX];
    void     *advalue;
    uint64_t  adsize;
    struct adios_dimension_item_struct *dim;
};

///* list of anonymous dims */
//List anonymous_dim_list;

/* Need a struct to encapsulate open file info
 */
struct open_file {
    char                           fpath[ADIOS_PATH_MAX];
    char                           fname[ADIOS_PATH_MAX];
    struct adios_nssi_filter_data_struct *md;
    struct adios_file_struct      *f;

    List anonymous_dim_list;
};

/* list of variable offsets */
static List open_file_list;

///////////////////////////
// Global Variables
///////////////////////////
static int adios_nssi_filter_initialized = 0;

static int global_rank=-1;
struct adios_nssi_config nssi_cfg;

//static log_level adios_nssi_filter_debug_level;
static int DEBUG=0;


extern struct adios_transport_struct * adios_transports;
//struct adios_method_struct *submethod=NULL;
struct adios_method_struct *self=NULL;


///////////////////////////
// Function Definitions
///////////////////////////

struct adios_var_struct *vars_deep_copy(struct adios_var_struct *orig)
{
    struct adios_var_struct *new = NULL;
    struct adios_var_struct *current = NULL;

    if (orig) {
        new = current = (struct adios_var_struct *)malloc(sizeof(struct adios_var_struct));
        current->next = NULL;
        while(orig) {
            memcpy(current, orig, sizeof(struct adios_var_struct));

            current->name = strdup (orig->name);
            current->path = strdup (orig->path);
            current->type = orig->type;
            current->got_buffer = adios_flag_no;

            current->write_offset = 0;

            current->stats=0;
            current->bitmap=0;

            if (orig->dimensions) {

                // NCSU Statistics - copy stat to new var struct
                uint8_t count = adios_get_stat_set_count(orig->type);
                uint8_t idx = 0;
                uint64_t characteristic_size;
                uint8_t c;
                uint8_t j;

//                printf("vars_deep_copy(): copying %d stat sets\n", count);

                current->bitmap = orig->bitmap;
                current->stats = malloc (count * sizeof(struct adios_stat_struct *));

                // Set of characteristics will be repeated thrice for complex numbers
                for (c = 0; c < count; c ++)
                {
                    current->stats[c] = calloc(ADIOS_STAT_LENGTH, sizeof (struct adios_stat_struct));

                    j = idx = 0;
                    while (orig->bitmap >> j)
                    {
//                        printf("name(%s) j(%d) ((orig->bitmap >> j) & 1)(%d) orig->stats[%d][%d].data(%p)\n", orig->name, j, ((orig->bitmap >> j) & 1), c, idx, orig->stats[c][idx].data);

                        if ((orig->bitmap >> j) & 1)
                        {
                            if (orig->stats[c][idx].data != NULL) {
                                if (j == adios_statistic_hist)
                                {
                                    current->stats[c][idx].data = (struct adios_hist_struct *) malloc (sizeof(struct adios_hist_struct));

                                    struct adios_hist_struct * orig_hist = orig->stats[c][idx].data;
                                    struct adios_hist_struct * current_hist = current->stats[c][idx].data;

                                    current_hist->min = orig_hist->min;
                                    current_hist->max = orig_hist->max;
                                    current_hist->num_breaks = orig_hist->num_breaks;

                                    current_hist->frequencies = malloc ((orig_hist->num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
                                    memcpy (current_hist->frequencies, orig_hist->frequencies, (orig_hist->num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
                                    current_hist->breaks = malloc ((orig_hist->num_breaks) * adios_get_type_size(adios_double, ""));
                                    memcpy (current_hist->breaks, orig_hist->breaks, (orig_hist->num_breaks) * adios_get_type_size(adios_double, ""));
                                }
                                else
                                {
                                    characteristic_size = adios_get_stat_size(orig->stats[c][idx].data, orig->type, j);
                                    current->stats[c][idx].data = malloc (characteristic_size);
                                    memcpy (current->stats[c][idx].data, orig->stats[c][idx].data, characteristic_size);
                                }

                                idx++;
                            }
                        }
                        j++;
                    }
                }
                // NCSU - End of copy, for statistics
            }

            current->free_data = adios_flag_no;
            current->data = 0;
            current->data_size = 0;

            current->next = NULL;
            orig = orig->next;
            if (orig) {
                current->next = (struct adios_var_struct *)malloc(sizeof(struct adios_var_struct));
                current = current->next;
            }
        }
    }

    return(new);
}

struct adios_attribute_struct *attrs_deep_copy(
        struct adios_attribute_struct *orig,
        struct adios_group_struct *new_group)
{
    struct adios_attribute_struct *new = NULL;
    struct adios_attribute_struct *current = NULL;

    if (orig) {
        new = current = (struct adios_attribute_struct *)malloc(sizeof(struct adios_attribute_struct));
        current->next = NULL;
        while(orig) {
            memcpy(current, orig, sizeof(struct adios_attribute_struct));


            current->name = strdup (orig->name);
            current->path = strdup (orig->path);
            current->type = orig->type;

            if (orig->var == 0) {
                uint64_t size = adios_get_type_size(orig->type, orig->value);
                current->value=malloc(size);
                memcpy(current->value, orig->value, size);
            } else {
                current->value = 0;
                current->type = adios_unknown;
                current->var = adios_find_var_by_id(new_group->vars,
                                                    orig->var->id);
            }

            current->next = NULL;
            orig = orig->next;
            if (orig) {
                current->next = (struct adios_attribute_struct *)malloc(sizeof(struct adios_attribute_struct));
                current = current->next;
            }
        }
    }

    return(new);
}

struct adios_group_struct *group_deep_copy(struct adios_group_struct *orig)
{
    struct adios_group_struct *new = NULL;

    if (orig) {
        new = (struct adios_group_struct *)malloc(sizeof(struct adios_group_struct));
        memcpy(new, orig, sizeof(struct adios_group_struct));

        new->vars = vars_deep_copy(orig->vars);
        new->attributes = attrs_deep_copy(orig->attributes, new);
    }

    return(new);
}



struct adios_method_struct *init_submethod(
        const char * method,
        const char * parameters)
{
    int64_t group_id;
    struct adios_group_struct * g;
    struct adios_method_struct * new_method;
    int requires_group_comm = 0;

    new_method = (struct adios_method_struct *)
                           malloc (sizeof (struct adios_method_struct));

    new_method->m = ADIOS_METHOD_UNKNOWN;
    new_method->base_path = strdup(self->base_path);
    new_method->method = strdup (method);
    new_method->parameters = strdup (parameters);
    new_method->iterations = self->iterations;
    new_method->priority = self->priority;
    new_method->method_data = 0;
    new_method->group = 0;

    if (adios_parse_method (method, &new_method->m, &requires_group_comm))
    {
        if (   new_method->m != ADIOS_METHOD_UNKNOWN
            && new_method->m != ADIOS_METHOD_NULL
            && adios_transports [new_method->m].adios_init_fn
           )
        {
            adios_transports[new_method->m].adios_init_fn
                                       (parameters, new_method);
        }
    }
    else
    {
        fprintf (stderr, "config.xml: invalid transport: %s\n", method);

        free (new_method->method);
        free (new_method->parameters);
        free (new_method);

        return NULL;
    }

    return new_method;
}


static struct anonymous_dim *anonymous_dim_create(const char *path, const char *pathname, void *value, uint64_t size, struct adios_dimension_item_struct *dim)
{
    struct anonymous_dim *ad=calloc(1,sizeof(struct anonymous_dim));

    strcpy(ad->adpath, path);
    strcpy(ad->adpathname, pathname);
    ad->advalue = value;
    ad->adsize  = size;
    ad->dim     = dim;

    return(ad);
}
static void anonymous_dim_free(void *ad)
{
    free(ad);
}
static int anonymous_dim_equal(const struct anonymous_dim *ad1, const struct anonymous_dim *ad2)
{
    if (strcmp(ad1->adpathname, ad2->adpathname) == 0) return TRUE;

    return FALSE;
}
static struct anonymous_dim *anonymous_dim_find(List *ad_list, const char *pathname)
{
    ListElmt *elmt;
    struct anonymous_dim *ad;

    if (DEBUG>3) printf("looking for adpathname(%s)\n", pathname);

    elmt = list_head(ad_list);
    while(elmt) {
        ad = list_data(elmt);
        if (DEBUG>3) printf("comparing to adpathname(%s)\n", ad->adpathname);
        if (strcmp(pathname, ad->adpathname) == 0) {
            if (DEBUG>3) printf("adpathname(%s) matches search\n", ad->adpathname);
            return ad;
        }
        elmt = list_next(elmt);
    }

    return NULL;
}
static void anonymous_dim_printall(List *ad_list)
{
    ListElmt *elmt;
    struct anonymous_dim *ad;

    elmt = list_head(ad_list);
    while(elmt) {
        ad = list_data(elmt);
        if (DEBUG>3) printf("adpathname(%s)\n", ad->adpathname);
        elmt = list_next(elmt);
    }
}

static struct open_file *open_file_create(
        const char *path,
        const char *name,
        struct adios_nssi_filter_data_struct *method_private_data,
        struct adios_file_struct *f)
{
    struct open_file *of=calloc(1,sizeof(struct open_file));

    strcpy(of->fpath, path);
    strcpy(of->fname, name);
    of->md = method_private_data;
    of->f = f;

    list_init(&of->anonymous_dim_list, anonymous_dim_free);


    return(of);
}
static void open_file_free(void *of)
{
    list_destroy(&((struct open_file *)of)->anonymous_dim_list);
    free(of);
}
static int open_file_equal(const struct open_file *of1, const struct open_file *of2)
{
    if ((strcmp(of1->fpath, of2->fpath) == 0) && (strcmp(of1->fname, of2->fname) == 0)) return TRUE;

    return FALSE;
}
static struct open_file *open_file_find(const char *path, const char *name)
{
    ListElmt *elmt;
    struct open_file *of;

    if (DEBUG>3) printf("looking for fpath(%s) fname(%s)\n", path, name);

    elmt = list_head(&open_file_list);
    while(elmt) {
        of = list_data(elmt);
        if (DEBUG>3) printf("comparing to fpath(%s) fname(%s)\n", of->fpath, of->fname);
        if ((strcmp(path, of->fpath) == 0) && (strcmp(name, of->fname) == 0)) {
            if (DEBUG>3) printf("fpath(%s) fname(%s) matches search\n", of->fpath, of->fname);
            return of;
        }
        elmt = list_next(elmt);
    }

    return NULL;
}
static struct open_file *open_file_delete(const char *path, const char *name)
{
    ListElmt *elmt, *prev;
    struct open_file *of;

    if (DEBUG>3) printf("trying to delete fpath(%s) fname(%s)\n", path, name);

    prev = elmt = list_head(&open_file_list);
    while(elmt) {
        of = list_data(elmt);
        if (DEBUG>3) printf("comparing to fpath(%s) fname(%s)\n", of->fpath, of->fname);
        if ((strcmp(path, of->fpath) == 0) && (strcmp(name, of->fname) == 0)) {
            if (DEBUG>3) printf("fpath(%s) fname(%s) matches search\n", of->fpath, of->fname);
            if (list_is_head(&open_file_list, elmt)) {
                list_rem_next(&open_file_list, NULL, &of);
            } else {
                list_rem_next(&open_file_list, prev, &of);
            }
        }
        prev = elmt;
        elmt = list_next(elmt);
    }

    return NULL;
}
static void open_file_printall(void)
{
    ListElmt *elmt;
    struct open_file *of;

    elmt = list_head(&open_file_list);
    while(elmt) {
        of = list_data(elmt);
        if (DEBUG>3) printf("fpath(%s) fname(%s)\n", of->fpath, of->fname);
        elmt = list_next(elmt);
    }
}

//static struct var_offset *var_offset_create(const char *path, const char *name, void *value, uint64_t size)
//{
//    struct var_offset *vo=calloc(1,sizeof(struct var_offset));
//
//    strcpy(vo->opath, path);
//    strcpy(vo->oname, name);
//    vo->ovalue = value;
//    vo->osize  = size;
//
//    return(vo);
//}
//static void var_offset_free(void *vo)
//{
//    free(vo);
//}
//static int var_offset_equal(const struct var_offset *vo1, const struct var_offset *vo2)
//{
//    if ((strcmp(vo1->opath, vo2->opath) == 0) && (strcmp(vo1->oname, vo2->oname) == 0)) return TRUE;
//
//    return FALSE;
//}
//static struct var_offset *var_offset_find(const char *path, const char *name)
//{
//    ListElmt *elmt;
//    struct var_offset *vo;
//
//    if (DEBUG>3) printf("looking for opath(%s) oname(%s)\n", path, name);
//
//    elmt = list_head(&var_offset_list);
//    while(elmt) {
//        vo = list_data(elmt);
//        if (DEBUG>3) printf("comparing to opath(%s) oname(%s)\n", vo->opath, vo->oname);
//        if ((strcmp(path, vo->opath) == 0) && (strcmp(name, vo->oname) == 0)) {
//            if (DEBUG>3) printf("opath(%s) oname(%s) matches search\n", vo->opath, vo->oname);
//            return vo;
//        }
//        elmt = list_next(elmt);
//    }
//
//    return NULL;
//}
//static void var_offset_printall(void)
//{
//    ListElmt *elmt;
//    struct var_offset *vo;
//
//    elmt = list_head(&var_offset_list);
//    while(elmt) {
//        vo = list_data(elmt);
//        if (DEBUG>3) printf("opath(%s) oname(%s)\n", vo->opath, vo->oname);
//        elmt = list_next(elmt);
//    }
//}
//
//static struct var_dim *var_dim_create(const char *path, const char *name, void *value, uint64_t size)
//{
//    struct var_dim *vd=calloc(1,sizeof(struct var_dim));
//
//    strcpy(vd->dpath, path);
//    strcpy(vd->dname, name);
//    vd->dvalue = value;
//    vd->dsize  = size;
//
//    return(vd);
//}
//static void var_dim_free(void *vd)
//{
//    free(vd);
//}
//static int var_dim_equal(const struct var_dim *vd1, const struct var_dim *vd2)
//{
//    if ((strcmp(vd1->dpath, vd2->dpath) == 0) && (strcmp(vd1->dname, vd2->dname) == 0)) return TRUE;
//
//    return FALSE;
//}
//static struct var_dim *var_dim_find(const char *path, const char *name)
//{
//    ListElmt *elmt;
//    struct var_dim *vd;
//
//    if (DEBUG>3) printf("looking for opath(%s) oname(%s)\n", path, name);
//
//    elmt = list_head(&var_dim_list);
//    while(elmt) {
//        vd = list_data(elmt);
//        if (DEBUG>3) printf("comparing to opath(%s) oname(%s)\n", vd->dpath, vd->dname);
//        if ((strcmp(path, vd->dpath) == 0) && (strcmp(name, vd->dname) == 0)) {
//            if (DEBUG>3) printf("opath(%s) oname(%s) matches search\n", vd->dpath, vd->dname);
//            return vd;
//        }
//        elmt = list_next(elmt);
//    }
//
//    return NULL;
//}
//static void var_dim_printall(void)
//{
//    ListElmt *elmt;
//    struct var_dim *vd;
//
//    elmt = list_head(&var_dim_list);
//    while(elmt) {
//        vd = list_data(elmt);
//        if (DEBUG>3) printf("dpath(%s) dname(%s)\n", vd->dpath, vd->dname);
//        elmt = list_next(elmt);
//    }
//}

static void parse_dimension_size(
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_dimension_item_struct *dim,
        size_t *dimsize)
{
    struct adios_var_struct *var_linked = NULL;
    struct adios_attribute_struct *attr_linked;
    if (dim->id) {
        var_linked = adios_find_var_by_id (pvar_root , dim->id);
        if (!var_linked) {
            attr_linked = adios_find_attribute_by_id (patt_root, dim->id);
            if (!attr_linked->var) {
                switch (attr_linked->type) {
                case adios_unsigned_byte:
                    *dimsize = *(uint8_t *)attr_linked->value;
                    break;
                case adios_byte:
                    *dimsize = *(int8_t *)attr_linked->value;
                    break;
                case adios_unsigned_short:
                    *dimsize = *(uint16_t *)attr_linked->value;
                    break;
                case adios_short:
                    *dimsize = *(int16_t *)attr_linked->value;
                    break;
                case adios_unsigned_integer:
                    *dimsize = *(uint32_t *)attr_linked->value;
                    break;
                case adios_integer:
                    *dimsize = *(int32_t *)attr_linked->value;
                    break;
                case adios_unsigned_long:
                    *dimsize = *(uint64_t *)attr_linked->value;
                    break;
                case adios_long:
                    *dimsize = *(int64_t *)attr_linked->value;
                    break;
                default:
                    fprintf (stderr, "Invalid datatype for array dimension on "
                            "var %s: %s\n"
                            ,attr_linked->name
                            ,adios_type_to_string_int (var_linked->type)
                    );
                    break;
                }
            } else {
                var_linked = attr_linked->var;
            }
        }
        if (var_linked && var_linked->data) {
            *dimsize = *(int *)var_linked->data;
        }
    } else {
        if (dim->time_index == adios_flag_yes) {
            *dimsize = 1;
        } else {
            *dimsize = dim->rank;
        }
    }

    return;
}
static void parse_dimension_name(
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_dimension_item_struct *dim,
        char *dimname)
{
    struct adios_var_struct *var_linked = NULL;
    struct adios_attribute_struct *attr_linked;
    if (dim->id) {
        var_linked = adios_find_var_by_id (pvar_root , dim->id);
        if (!var_linked) {
            attr_linked = adios_find_attribute_by_id (patt_root, dim->id);
            if (!attr_linked->var) {
//                strcpy(dimname, attr_linked->name);
                sprintf(dimname, "%s", attr_linked->name);
            } else {
                var_linked = attr_linked->var;
            }
        }
        if (var_linked && var_linked->name) {
//            strcpy(dimname, var_linked->name);
            sprintf(dimname, "%s", var_linked->name);
        }
    } else {
        if (dim->time_index == adios_flag_yes) {
//            strcpy(dimname, group->time_index_name);
            sprintf(dimname, "%s", group->time_index_name);
        } else {
            dimname[0] = '\0';
        }
    }

    return;
}
static int is_anonymous_dimension(
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_dimension_item_struct *dim)
{
    int rc=FALSE;
    struct adios_var_struct       *var_linked = NULL;
    struct adios_attribute_struct *attr_linked = NULL;

    if (dim->id) {
        var_linked = adios_find_var_by_id (pvar_root , dim->id);
        if (!var_linked) {
            attr_linked = adios_find_attribute_by_id (patt_root, dim->id);
            if (!attr_linked->var) {
                rc=FALSE;
            } else {
                var_linked = attr_linked->var;
            }
        }
        if (var_linked && var_linked->name) {
            rc=FALSE;
        } else {
            rc=TRUE;
        }
    } else {
        if (dim->time_index == adios_flag_yes) {
            rc=FALSE;
        } else {
            rc=TRUE;
        }
    }

    return rc;
}


static int gen_anonymous_dim_list(
        struct open_file *of,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root)
{
    struct adios_var_struct *v;
    struct adios_dimension_struct *dims;
    struct anonymous_dim *ad;
    char ad_name[255];
    uint64_t *value;

    v = pvar_root;
    while (v) {
        dims=v->dimensions;
        int ad_idx=0;
        uint64_t vsize = 8;

        while (dims) {
            if (is_anonymous_dimension(pvar_root, patt_root, &dims->dimension) == TRUE) {
                sprintf(ad_name, "%s_%s_dim_%d", /*v->path*/"", v->name, ad_idx);
                if (DEBUG>3) printf("gen: dim_name(%s)\n", ad_name);

                value=calloc(1,sizeof(uint64_t));
                parse_dimension_size(group, pvar_root, patt_root, &dims->dimension, value);

                ad = anonymous_dim_create(v->path, ad_name, value, vsize, &dims->dimension);
                list_ins_next(&of->anonymous_dim_list, list_tail(&of->anonymous_dim_list), ad);
            }

            if (is_anonymous_dimension(pvar_root, patt_root, &dims->global_dimension) == TRUE) {
                sprintf(ad_name, "%s_%s_global_%d", /*v->path*/"", v->name, ad_idx);
                if (DEBUG>3) printf("gen: gdim_name(%s)\n", ad_name);

                value=calloc(1,sizeof(uint64_t));
                parse_dimension_size(group, pvar_root, patt_root, &dims->global_dimension, value);

                ad = anonymous_dim_create(v->path, ad_name, value, vsize, &dims->global_dimension);
                list_ins_next(&of->anonymous_dim_list, list_tail(&of->anonymous_dim_list), ad);
            }

            if (is_anonymous_dimension(pvar_root, patt_root, &dims->local_offset) == TRUE) {
                sprintf(ad_name, "%s_%s_offset_%d", /*v->path*/"", v->name, ad_idx);
                if (DEBUG>3) printf("gen: ad_name(%s)\n", ad_name);

                value=calloc(1,sizeof(uint64_t));
                parse_dimension_size(group, pvar_root, patt_root, &dims->local_offset, value);

                ad = anonymous_dim_create(v->path, ad_name, value, vsize, &dims->local_offset);
                list_ins_next(&of->anonymous_dim_list, list_tail(&of->anonymous_dim_list), ad);
            }

            ad_idx++;
            dims = dims->next;
        }
        v = v->next;
    }
}


void adios_nssi_filter_init(
        const char *parameters,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    int verbose=5;
    char logfile[1024];
    int log_rank;

    struct adios_nssi_filter_data_struct *self_md=NULL;

    if (!adios_nssi_filter_initialized) {
        adios_nssi_filter_initialized = 1;
    }

    MPI_Comm_rank(method->init_comm, &global_rank);

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_init\n", global_rank);

//    sprintf(logfile, "%s.%04d", "adios_nssi_filter_client.log", log_rank);
//    logger_init((log_level)verbose, logfile);

//    logger_init((log_level)verbose, NULL);

    if (self == NULL) {
        self=method;
        self_md=calloc(1, sizeof(struct adios_nssi_filter_data_struct));
        self->method_data=self_md;
    } else {
        self_md=(struct adios_nssi_filter_data_struct *)self->method_data;
    }

    /*
     * initialize sub_method
     */
    char *sm_method=NULL;
    char *sm_parameters=NULL;

    if (DEBUG>3) printf("parameters=%s\n", parameters);

    char *pkey=NULL;
    char *pkey_str="submethod=\"";
    char *pvalue=NULL;
    int   pvalue_len=0;
    pkey=strstr(parameters, pkey_str);
    if (pkey != NULL) {
        pvalue=pkey+strlen(pkey_str);
        char *quote_pos=strchr(pvalue, '"');
        pvalue_len=quote_pos-pvalue+1;
        sm_method=(char *)malloc(pvalue_len+1);
        strncpy(sm_method, pvalue, pvalue_len);
        sm_method[pvalue_len-1]='\0';
    }

    if (DEBUG>3) printf("sm_method=%s\n", sm_method);

    if (sm_method != NULL) {
        pkey=NULL;
        pkey_str="subparameters=\"";
        pvalue=NULL;
        pvalue_len=0;
        pkey=strstr(parameters, pkey_str);
        if (pkey != NULL) {
            pvalue=pkey+strlen(pkey_str);
            char *quote_pos=strchr(pvalue, '"');
            pvalue_len=quote_pos-pvalue+1;
            sm_parameters=(char *)malloc(pvalue_len+1);
            strncpy(sm_parameters, pvalue, pvalue_len);
            sm_parameters[pvalue_len-1]='\0';
        }
    }

    if (DEBUG>3) printf("sm_parameters=%s\n", sm_parameters);

    if (sm_parameters == NULL) {
        sm_parameters="";
    }

    self_md->sm_method = strdup(sm_method);
    self_md->sm_parameters = strdup(sm_parameters);


    list_init(&open_file_list, open_file_free);
//    list_init(&var_offset_list, var_offset_free);
//    list_init(&var_dim_list, var_dim_free);

    return;
}


enum ADIOS_FLAG adios_nssi_filter_should_buffer(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;

    struct open_file *of=NULL;
    struct adios_nssi_filter_data_struct *md=NULL;

    enum ADIOS_FLAG sm_should_buffer=adios_flag_no;
    enum ADIOS_FLAG old_shared_buffer=adios_flag_no;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_should_buffer\n", global_rank);

    if (DEBUG>3) {
        struct adios_var_struct *v=f->group->vars;
        while(v) {
            printf("adios_nssi_filter_should_buffer: fname(%s) vname(%s)\n", f->name, v->name);
            v=v->next;
        }
    }

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "nssi_filter_should_buffer: file(%s, %s) is not open.  FAIL.", method->base_path, f->name);
        return adios_flag_no;
    }
    md=of->md;

    /*
     * call sub_method
     */
    if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
        && md->submethod->m != ADIOS_METHOD_NULL
        && adios_transports[md->submethod->m].adios_should_buffer_fn
       )
    {
        sm_should_buffer = adios_transports[md->submethod->m].adios_should_buffer_fn
                                                                (f, md->submethod);
    }

    if (DEBUG>3) printf("sm_should_buffer==%d\n", sm_should_buffer);
    return sm_should_buffer;
}

int adios_nssi_filter_open(
        struct adios_file_struct *f,
        struct adios_method_struct *method,
        MPI_Comm comm)
{
    int rc=NSSI_OK;

    struct open_file *of=NULL;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("global_rank(%d): enter adios_nssi_filter_open (%s)\n", global_rank, f->name);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        struct adios_nssi_filter_data_struct *self_md=(struct adios_nssi_filter_data_struct *)self->method_data;;

        md             = malloc(sizeof(struct adios_nssi_filter_data_struct));
        md->fd         = -1;
        md->rank       = -1;
        md->size       = 0;
        md->group_comm = comm;

        md->submethod = init_submethod(self_md->sm_method, self_md->sm_parameters);
        md->submethod->group = group_deep_copy(f->group);
        f->group = md->submethod->group;

        of=open_file_create(method->base_path, f->name, md, f);
    } else {
        md=of->md;

        // sanity check
        if (md->fd == -1) {
            if (DEBUG>3) printf("open: %s is open but fd==-1.  sanity check failed.  attempting reopen.\n", f->name);
            open_file_delete(of->fpath, of->fname);
        } else {
            // file already open
            return adios_flag_no;
        }
    }

    if (md->group_comm != MPI_COMM_NULL) {
        if (DEBUG>3) printf("global_rank(%d): adios_nssi_filter_open: get rank and size: group_comm(%p)\n", global_rank, md->group_comm);
        MPI_Comm_rank(md->group_comm, &md->rank);
        MPI_Comm_size(md->group_comm, &md->size);
        if (DEBUG>3) printf("global_rank(%d): adios_nssi_filter_open: size(%d) rank(%d)\n", global_rank, md->size, md->rank);
    } else {
        md->group_comm=MPI_COMM_SELF;
    }
    f->group->process_id = md->rank;


    gen_anonymous_dim_list(
            of,
            f->group,
            f->group->vars,
            f->group->attributes);


//    gen_offset_list(
//            f->group,
//            f->group->vars,
//            f->group->attributes);
//    gen_dim_list(
//            f->group,
//            f->group->vars,
//            f->group->attributes);
//    var_offset_printall();


    /*
     * call sub_method
     */
    if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
        && md->submethod->m != ADIOS_METHOD_NULL
        && adios_transports[md->submethod->m].adios_open_fn
       )
    {
        adios_transports[md->submethod->m].adios_open_fn
                                   (f, md->submethod, md->group_comm);
    }


    list_ins_next(&open_file_list, list_tail(&open_file_list), of);

    return 1;
}

void adios_nssi_filter_start_calculation(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int myrank;

    ListElmt *elmt;
    struct open_file *of;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_start_calc\n", global_rank);

    elmt = list_head(&open_file_list);
    while(elmt) {
        of = list_data(elmt);
        md = of->md;

        /*
         * call sub_method
         */
        if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
            && md->submethod->m != ADIOS_METHOD_NULL
            && adios_transports[md->submethod->m].adios_start_calculation_fn
           )
        {
            adios_transports[md->submethod->m].adios_start_calculation_fn
                                       (md->submethod);
        }

        elmt = list_next(elmt);
    }

    if (DEBUG>3) printf("rank(%d) exit adios_nssi_filter_start_calc\n", global_rank);

    return;
}

void adios_nssi_filter_end_iteration(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int myrank;

    ListElmt *elmt;
    struct open_file *of;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_end_iteration\n", global_rank);

    elmt = list_head(&open_file_list);
    while(elmt) {
        of = list_data(elmt);
        md = of->md;

        /*
         * call sub_method
         */
        if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
            && md->submethod->m != ADIOS_METHOD_NULL
            && adios_transports[md->submethod->m].adios_end_iteration_fn
        )
        {
            adios_transports[md->submethod->m].adios_end_iteration_fn
                                        (md->submethod);
        }

        elmt = list_next(elmt);
    }

    if (DEBUG>3) printf("rank(%d) exit adios_nssi_filter_end_iteration\n", global_rank);

    return;
}

void adios_nssi_filter_stop_calculation(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int remote_rc=NSSI_OK;
    int myrank;

    ListElmt *elmt;
    struct open_file *of;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_stop_calc\n", global_rank);

    elmt = list_head(&open_file_list);
    while(elmt) {
        of = list_data(elmt);
        md = of->md;

        /*
         * call sub_method
         */
        if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
            && md->submethod->m != ADIOS_METHOD_NULL
            && adios_transports[md->submethod->m].adios_stop_calculation_fn
        )
        {
            adios_transports[md->submethod->m].adios_stop_calculation_fn
                                        (md->submethod);
        }

        elmt = list_next(elmt);
    }

    if (DEBUG>3) printf("rank(%d) exit adios_nssi_filter_stop_calc\n", global_rank);

    return;
}

void adios_nssi_filter_write(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *data,
        struct adios_method_struct *method)
{
    static int first_write = 1;

    struct open_file *of=NULL;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_write - var(%s, %s)\n", global_rank, v->path, v->name);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "nssi_filter_write: file(%s, %s) is not open.  FAIL.", method->base_path, f->name);
        return;
    }
    md=of->md;

    /*
     * call sub_method
     */
    if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
        && md->submethod->m != ADIOS_METHOD_NULL
        && adios_transports[md->submethod->m].adios_write_fn
       )
    {
        adios_transports[md->submethod->m].adios_write_fn
                                   (f, v, data, md->submethod);
    }


//    if (f->mode == adios_mode_write || f->mode == adios_mode_append) {
//
//        if (md->rank==0) {
//            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
//            if (DEBUG>3) fprintf(stderr, "write var: %s start!\n", v->name);
//        }
//        uint64_t var_size = adios_get_var_size (v, f->group, data);
//        if (DEBUG>3) printf("vname(%s) vsize(%ld)\n", v->name, var_size);
//        write_var(md->fd,
//                f->group,
//                f->group->vars,
//                f->group->attributes,
//                v,
//                var_size,
//                f->group->adios_host_language_fortran,
//                md->rank,
//                md->size,
//                md->group_comm);
//    } else {
//        if (DEBUG>3) fprintf(stderr, "entering unknown nc4 mode %d!\n", f->mode);
//    }
//    if (md->rank==0) {
//        if (DEBUG>3) fprintf(stderr, "write var: %s end!\n", v->name);
//        if (DEBUG>3) fprintf(stderr, "-------------------------\n");
//    }

    return;
}


void adios_nssi_filter_read(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *buffer,
        uint64_t buffersize,
        struct adios_method_struct *method)
{
    struct open_file *of=NULL;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_read\n", global_rank);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "nssi_filter_read: file(%s, %s) is not open.  FAIL.", method->base_path, f->name);
        return;
    }
    md=of->md;

    /*
     * call sub_method
     */
    if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
        && md->submethod->m != ADIOS_METHOD_NULL
        && adios_transports[md->submethod->m].adios_read_fn
       )
    {
        adios_transports[md->submethod->m].adios_read_fn
                                   (f, v, buffer, buffersize, md->submethod);
    }


//    if(f->mode == adios_mode_read) {
//        v->data = buffer;
//        v->data_size = buffersize;
//
//        if (md->rank==0) {
//            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
//            if (DEBUG>3) fprintf(stderr, "read var: %s! start\n", v->name);
//        }
//        read_var(md->fd,
//                v,
//                md->rank,
//                md->size,
//                md->group_comm);
//        if (md->rank==0) {
//            if (DEBUG>3) fprintf(stderr, "read var: %s! end\n", v->name);
//            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
//        }
//    }

    return;
}

void adios_nssi_filter_close(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    struct adios_attribute_struct * a = f->group->attributes;
    int myrank;

    struct open_file *of=NULL;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("global_rank(%d) enter adios_nssi_filter_close\n", global_rank);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "nssi_filter_close: file(%s, %s) is not open.  FAIL.", method->base_path, f->name);
        return;
    }
    md=of->md;
    myrank=md->rank;

    if (DEBUG>3) printf("myrank(%d) enter adios_nssi_filter_close\n", myrank);

    /*
     * call sub_method
     */
    if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
        && md->submethod->m != ADIOS_METHOD_NULL
        && adios_transports[md->submethod->m].adios_close_fn
       )
    {
        adios_transports[md->submethod->m].adios_close_fn
                                   (f, md->submethod);
    }

    if (f->mode == adios_mode_read) {
        struct adios_var_struct * v = f->group->vars;
        while (v)
        {
            v->data = 0;
            v = v->next;
        }

        if (md->rank==0) {
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
            if (DEBUG>1) fprintf(stderr, "reading done, NSSI_FILTER file is virtually closed;\n");
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
        }
    } else if (f->mode == adios_mode_write || f->mode == adios_mode_append) {
        //fprintf(stderr, "entering nc4 write attribute mode!\n");
        if (md->rank==0) {
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
            if (DEBUG>1) fprintf(stderr, "writing done, NSSI_FILTER file is virtually closed;\n");
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
        }
    }

//    struct adios_group_struct *group_clone=f->group;
//    if (group_clone) {
//        while (group_clone->vars)
//        {
//            struct adios_var_struct *vars = group_clone->vars->next;
//            free (group_clone->vars);
//            group_clone->vars = vars;
//        }
//
//        while (group_clone->attributes)
//        {
//            struct adios_attribute_struct * attributes = group_clone->attributes->next;
//            free (group_clone->attributes);
//            group_clone->attributes = attributes;
//        }
//
//        free(group_clone);
//        f->group=NULL;
//    }

    open_file_delete(method->base_path, f->name);
    free(md);
    of=NULL;
    md=NULL;

    if (DEBUG>3) printf("global_rank(%d) exit adios_nssi_filter_close\n", global_rank);

    return;
}

void adios_nssi_filter_finalize(
        int mype,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    int myrank;

    ListElmt *elmt;
    struct open_file *of;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_finalize\n", global_rank);

    elmt = list_head(&open_file_list);
    while(elmt) {
        of = list_data(elmt);
        md = of->md;

        /*
         * call sub_method
         */
        if (   md->submethod->m != ADIOS_METHOD_UNKNOWN
            && md->submethod->m != ADIOS_METHOD_NULL
            && adios_transports[md->submethod->m].adios_finalize_fn
        )
        {
            adios_transports[md->submethod->m].adios_finalize_fn
                                        (mype, md->submethod);
        }

        elmt = list_next(elmt);
    }

    free_nssi_config(&nssi_cfg);

    if (adios_nssi_filter_initialized) {
        adios_nssi_filter_initialized = 0;
        self=NULL;
    }

    if (DEBUG>3) printf("rank(%d) exit adios_nssi_filter_finalize\n", global_rank);
}


int adios_nssi_filter_is_anon_dim(
        int fd,
        const char *dimname)
{
    struct adios_file_struct *f=(struct adios_file_struct *)fd;
    int myrank;

    struct open_file *of=NULL;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("global_rank(%d) enter BACKDOOR adios_nssi_filter_is_anon_dim\n", global_rank);

    if ((!adios_nssi_filter_initialized) || (self == NULL)) {
        return FALSE;
    }

    of=open_file_find(self->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "nssi_filter_is_anon: file(%s, %s) is not open.  FAIL.", self->base_path, f->name);
        return FALSE;
    }
    md=of->md;
    myrank=md->rank;

    if (anonymous_dim_find(&of->anonymous_dim_list, dimname) != NULL) {
        return TRUE;
    }
    return FALSE;
}
void adios_nssi_filter_set_anon_dim(
        int fd,
        const char *dimname,
        const uint64_t dimvalue)
{
    struct adios_file_struct *f=(struct adios_file_struct *)fd;
    int myrank;

    struct open_file *of=NULL;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("global_rank(%d) enter BACKDOOR adios_nssi_filter_set_anon_dim\n", global_rank);

    if ((!adios_nssi_filter_initialized) || (self == NULL)) {
        return;
    }

    of=open_file_find(self->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "nssi_filter_set_anon: file(%s, %s) is not open.  FAIL.", self->base_path, f->name);
        return;
    }
    md=of->md;
    myrank=md->rank;

    struct anonymous_dim *ad=anonymous_dim_find(&of->anonymous_dim_list, dimname);
    if (ad == NULL) {
        fprintf(stderr, "couldn't find anonymous dimension (%s)\n", dimname);
        return;
    }
    ad->dim->rank=dimvalue;

    return;
}
