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
#include "mpi.h"
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

#ifndef HAVE_NSSI
void adios_nssi_filter_init(
        const char *parameters,
        struct adios_method_struct * method)
{
}
void adios_nssi_filter_finalize(
        int mype,
        struct adios_method_struct * method)
{
}
enum ADIOS_FLAG adios_nssi_filter_should_buffer(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    return adios_flag_unknown;
}
int adios_nssi_filter_open(
        struct adios_file_struct *f,
        struct adios_method_struct *method,
        void * comm)
{
    return -1;
}
void adios_nssi_filter_close(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
}
void adios_nssi_filter_write(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *data,
        struct adios_method_struct *method)
{
}
void adios_nssi_filter_read(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *buffer,
        uint64_t buffersize,
        struct adios_method_struct *method)
{
}
#else

///////////////////////////
// Datatypes
///////////////////////////
struct adios_nssi_filter_data_struct
{
    int      fd;
    MPI_Comm group_comm;
    int      rank;
    int      size;

    void * comm; // temporary until moved from should_buffer to open
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
List open_file_list;

///////////////////////////
// Global Variables
///////////////////////////
static int adios_nssi_filter_initialized = 0;

static int global_rank=-1;
struct adios_nssi_config nssi_cfg;

//static log_level adios_nssi_filter_debug_level;
static int DEBUG=3;


extern struct adios_transport_struct * adios_transports;
struct adios_method_struct *submethod=NULL;
struct adios_method_struct *self=NULL;


///////////////////////////
// Function Declarations
///////////////////////////


///////////////////////////
// Function Definitions
///////////////////////////

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

    new_method->group = self->group;


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
//				strcpy(dimname, attr_linked->name);
                sprintf(dimname, "%s", attr_linked->name);
            } else {
                var_linked = attr_linked->var;
            }
        }
        if (var_linked && var_linked->name) {
//			strcpy(dimname, var_linked->name);
            sprintf(dimname, "%s", var_linked->name);
        }
    } else {
        if (dim->time_index == adios_flag_yes) {
//			strcpy(dimname, group->time_index_name);
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

//static int gen_offset_list(
//        struct adios_group_struct *group,
//        struct adios_var_struct *pvar_root,
//        struct adios_attribute_struct *patt_root)
//{
//    struct adios_var_struct *v;
//    struct adios_dimension_struct *dims;
//    struct var_info *vi;
//    char offset_name[255];
//    uint64_t *value;
//
//    v = pvar_root;
//    while (v) {
//        dims=v->dimensions;
//        int loffs_idx=0;
//        while (dims) {
//            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, offset_name);
//            if (offset_name[0] == '\0') {
//                sprintf(offset_name, "offset_%d", loffs_idx);
//            }
//            if (DEBUG>3) printf("gen: offset_name(%s)\n", offset_name);
//            value=calloc(1,sizeof(uint64_t));
//            parse_dimension_size(group, pvar_root, patt_root, &dims->local_offset, value);
//
//            uint64_t vsize = 4; /* adios_get_var_size(v, group, value); */
//            vi = var_offset_create(v->path, offset_name, value, vsize);
//            list_ins_next(&var_offset_list, list_tail(&var_offset_list), vi);
//
//            loffs_idx++;
//            dims = dims->next;
//        }
//        v = v->next;
//    }
//}
//
//static void create_offset_list_for_var(
//        struct adios_write_args *args,
//        struct adios_var_struct *v,
//        struct adios_group_struct *group,
//        struct adios_var_struct *pvar_root,
//        struct adios_attribute_struct *patt_root)
//{
//    struct adios_dimension_struct *dims;
//    char offset_name[255];
//
//    args->offsets.offsets_len=0;
//    args->offsets.offsets_val=NULL;
//
//    if ((v) && (v->dimensions)) {
//        int local_offset_count=0;
//        dims=v->dimensions;
//        while (dims) {
//            if (dims->dimension.time_index == adios_flag_yes) {
//                dims = dims->next;
//                continue;
//            }
//            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, offset_name);
//            local_offset_count++;
//            dims = dims->next;
//        }
//
//        args->offsets.offsets_len=local_offset_count;
//        args->offsets.offsets_val=calloc(local_offset_count, sizeof(struct adios_var));
//
//        dims=v->dimensions;
//        int loffs_idx=0;
//        while (dims) {
//            if (dims->dimension.time_index == adios_flag_yes) {
//                dims = dims->next;
//                continue;
//            }
//            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, offset_name);
//            if (offset_name[0] == '\0') {
//                sprintf(offset_name, "offset_%d", loffs_idx);
//            }
//            args->offsets.offsets_val[loffs_idx].vpath=strdup(v->path);
//            args->offsets.offsets_val[loffs_idx].vname=strdup(offset_name);
////            struct var_offset *vo=var_offset_find("", offset_name);
////            memcpy(&(args->offsets.offsets_val[loffs_idx].vdata), vo->ovalue, vo->osize);
////            args->offsets.offsets_val[loffs_idx].vdatasize=vo->osize;
////            printf("create: offset_name(%s) offset_value(%lu)\n", offset_name, vo->ovalue);
//            uint64_t value=0;
//            parse_dimension_size(group, pvar_root, patt_root, &dims->local_offset, &value);
//            memcpy(&(args->offsets.offsets_val[loffs_idx].vdata), &value, 4);
//            args->offsets.offsets_val[loffs_idx].vdatasize=4;
//            if (DEBUG>3) printf("create: offset_name(%s) offset_value(%lu)\n", offset_name, value);
//
//            loffs_idx++;
//            dims = dims->next;
//        }
//    }
//}
//
//static int gen_dim_list(
//        struct adios_group_struct *group,
//        struct adios_var_struct *pvar_root,
//        struct adios_attribute_struct *patt_root)
//{
//    struct adios_var_struct *v;
//    struct adios_dimension_struct *dims;
//    struct var_info *vi;
//    char dim_name[255];
//    uint64_t *value;
//
//    v = pvar_root;
//    while (v) {
//        dims=v->dimensions;
//        int dim_idx=0;
//        while (dims) {
//            parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dim_name);
//            if (dim_name[0] == '\0') {
//                sprintf(dim_name, "dim_%d", dim_idx);
//            }
//            if (DEBUG>3) printf("gen: dim_name(%s)\n", dim_name);
//            value=calloc(1,sizeof(uint64_t));
//            parse_dimension_size(group, pvar_root, patt_root, &dims->dimension, value);
//            if (global_rank==0) {
//                if (DEBUG>3) printf(":o(%d)", *value);
//            }
//
//            uint64_t vsize = 4; /* adios_get_var_size(v, group, value); */
//            vi = var_dim_create(v->path, dim_name, value, vsize);
//            list_ins_next(&var_dim_list, list_tail(&var_dim_list), vi);
//
//            dim_idx++;
//            dims = dims->next;
//        }
//        v = v->next;
//    }
//}
//
//static void create_dim_list_for_var(
//        struct adios_write_args *args,
//        struct adios_var_struct *v,
//        struct adios_group_struct *group,
//        struct adios_var_struct *pvar_root,
//        struct adios_attribute_struct *patt_root)
//{
//    struct adios_dimension_struct *dims;
//    char dim_name[255];
//
//    args->dims.dims_len=0;
//    args->dims.dims_val=NULL;
//
//    if ((v) && (v->dimensions)) {
//        int dim_count=0;
//        dims=v->dimensions;
//        while (dims) {
//            if (dims->dimension.time_index == adios_flag_yes) {
//                dims = dims->next;
//                continue;
//            }
//            parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dim_name);
//            dim_count++;
//            dims = dims->next;
//        }
//
//        args->dims.dims_len=dim_count;
//        args->dims.dims_val=calloc(dim_count, sizeof(struct adios_var));
//
//        dims=v->dimensions;
//        int dim_idx=0;
//        while (dims) {
//            if (dims->dimension.time_index == adios_flag_yes) {
//                dims = dims->next;
//                continue;
//            }
//            parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dim_name);
//            if (dim_name[0] == '\0') {
//                sprintf(dim_name, "dim_%d", dim_idx);
//            }
//            args->dims.dims_val[dim_idx].vpath=strdup(v->path);
//            args->dims.dims_val[dim_idx].vname=strdup(dim_name);
////            struct var_dim *vd=var_dim_find("", dim_name);
////            memcpy(&(args->dims.dims_val[dim_idx].vdata), vd->dvalue, vd->dsize);
////            args->dims.dims_val[dim_idx].vdatasize=vd->dsize;
//            uint64_t value=0;
//            parse_dimension_size(group, pvar_root, patt_root, &dims->dimension, &value);
//            memcpy(&(args->dims.dims_val[dim_idx].vdata), &value, 4);
//            args->dims.dims_val[dim_idx].vdatasize=4;
//            if (DEBUG>3) printf("create: dim_name(%s) dvalue(%lu)\n", dim_name, value);
//            if (global_rank==0) {
//                if (DEBUG>3) printf(":o(%d)", value);
//            }
//
//            dim_idx++;
//            dims = dims->next;
//        }
//    }
//}

//static int read_var(
//        int fd,
//        struct adios_var_struct *pvar,
//        int myrank,
//        int nproc,
//        MPI_Comm group_comm)
//{
//    int return_code=0;
//    int i, rc;
//
//    adios_read_args args;
//    adios_read_res  res;
//
//    args.fd       = fd;
//    args.max_read = pvar->data_size;
//    args.vpath = strdup(pvar->path);
//    args.vname = strdup(pvar->name);
//
//    Func_Timer("ADIOS_READ_OP",
//            rc = nssi_call_rpc_sync(&svcs[default_svc],
//            ADIOS_READ_OP,
//            &args,
//            pvar->data,
//            pvar->data_size,
//            &res););
//    if (rc != NSSI_OK) {
//        //log_error(adios_nssi_filter_debug_level, "unable to call remote adios_read");
//        return_code=-2;
//    }
//
//    free(args.vpath);
//    free(args.vname);
//
//    return return_code;
//}
//static int write_var(
//        int fd,
//        struct adios_group_struct *group,
//        struct adios_var_struct *pvar_root,
//        struct adios_attribute_struct *patt_root,
//        struct adios_var_struct *pvar,
//        uint64_t var_size,
//        enum ADIOS_FLAG fortran_flag,
//        int myrank,
//        int nproc,
//        MPI_Comm group_comm)
//{
//    int i;
//    int rc;
//    int return_code=0;
//
//    adios_write_args args;
//    adios_write_res  res;
//
////    var_offset_printall();
//
//    memset(&args, 0, sizeof(adios_write_args));
//    args.fd    = fd;
//    args.vpath = strdup(pvar->path);
//    args.vname = strdup(pvar->name);
//    args.vsize = var_size;
//    args.atype = pvar->type;
//    if (pvar->dimensions) {
//        args.is_scalar = FALSE;
//    } else {
//        args.is_scalar = TRUE;
//    }
//    args.writer_rank=myrank;
//    args.offsets.offsets_len=0;
//    args.offsets.offsets_val=NULL;
//    args.dims.dims_len=0;
//    args.dims.dims_val=NULL;
//    if (pvar->dimensions) {
//        create_offset_list_for_var(
//                 &args,
//                 pvar,
//                 group,
//                 group->vars,
//                 group->attributes);
//        create_dim_list_for_var(
//                 &args,
//                 pvar,
//                 group,
//                 group->vars,
//                 group->attributes);
//     }
//
////    MPI_Barrier(group_comm);
//    Func_Timer("ADIOS_WRITE_OP",
//            rc = nssi_call_rpc_sync(&svcs[default_svc],
//            ADIOS_WRITE_OP,
//            &args,
//            pvar->data,
//            var_size,
//            &res););
//    if (rc != NSSI_OK) {
//        //log_error(adios_nssi_filter_debug_level, "unable to call remote adios_write");
//        return_code=-2;
//    }
////    MPI_Barrier(group_comm);
//
//    free(args.vpath);
//    free(args.vname);
//    free(args.offsets.offsets_val);
//
//    return return_code;
//}


static void adios_var_to_comm_nssi(
        enum ADIOS_FLAG host_language_fortran,
        void *data,
        MPI_Comm *comm)
{
    if (data) {
        int t = *(int *) data;
        if (host_language_fortran == adios_flag_yes) {
            *comm = MPI_Comm_f2c (t);
        } else {
            *comm = *(MPI_Comm *) data;
        }
    } else {
        fprintf (stderr, "coordination-communication not provided. "
                "Using md->group_comm instead\n");
        *comm = MPI_COMM_WORLD;
    }

    return;
}

void adios_nssi_filter_init(
        const char *parameters,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    int verbose=5;
    char logfile[1024];
    int log_rank;

    if (!adios_nssi_filter_initialized) {
        adios_nssi_filter_initialized = 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_init\n", global_rank);

//    MPI_Comm_rank(md->group_comm, &log_rank);
//    sprintf(logfile, "%s.%04d", "adios_nssi_filter_client.log", log_rank);
//    logger_init((log_level)verbose, logfile);

//    logger_init((log_level)verbose, NULL);

    self=method;

    /*
     * initialize sub_method
     */
    char *sm_method=NULL;
    char *sm_parameters=NULL;

    printf("parameters=%s\n", parameters);

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

    printf("sm_method=%s\n", sm_method);

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

    printf("sm_parameters=%s\n", sm_parameters);

    if (sm_parameters == NULL) {
        sm_parameters="";
    }

    submethod=init_submethod(sm_method, sm_parameters);
    if (submethod == NULL) {
        fprintf(stderr, "adios_nssi_filter: invalid submethod.  aborting.\n");
        adios_nssi_filter_initialized = 0;
        self=NULL;

        return;
    }


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

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "file is not open.  FAIL.");
        return adios_flag_no;
    }
    md=of->md;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_should_buffer\n", global_rank);

    /*
     * call sub_method
     */
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_should_buffer_fn
       )
    {
        sm_should_buffer = adios_transports[submethod->m].adios_should_buffer_fn
                                                                (f, submethod);
    }

    return sm_should_buffer;
}

int adios_nssi_filter_open(
        struct adios_file_struct *f,
        struct adios_method_struct *method,
        void *comm)
{
    int rc=NSSI_OK;

    struct open_file *of=NULL;
    struct adios_nssi_filter_data_struct *md=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_open\n", global_rank);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        md             = malloc(sizeof(struct adios_nssi_filter_data_struct));
        md->fd         = -1;
        md->rank       = -1;
        md->size       = 0;
        md->group_comm = MPI_COMM_NULL;

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

    if (DEBUG>3) printf("global_rank(%d): enter adios_nssi_filter_open (%s)\n", global_rank, f->name);

    md->comm = comm;
    if (DEBUG>3) printf("global_rank(%d): adios_nssi_filter_open: setup group_comm\n", global_rank);
    adios_var_to_comm_nssi(f->group->adios_host_language_fortran, md->comm, &md->group_comm);
    if (md->group_comm != MPI_COMM_NULL) {
        if (DEBUG>3) printf("global_rank(%d): adios_nssi_filter_open: get rank and size\n", global_rank);
        MPI_Comm_rank(md->group_comm, &md->rank);
        MPI_Comm_size(md->group_comm, &md->size);
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
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_open_fn
       )
    {
        adios_transports[submethod->m].adios_open_fn
                                   (f, submethod, comm);
    }


    list_ins_next(&open_file_list, list_tail(&open_file_list), of);

    return 1;
}

void adios_nssi_filter_start_calculation(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int myrank;


    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_start_calc\n", global_rank);

    /*
     * call sub_method
     */
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_start_calculation_fn
       )
    {
        adios_transports[submethod->m].adios_start_calculation_fn
                                   (submethod);
    }


    return;
}

void adios_nssi_filter_end_iteration(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int myrank;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_end_iter\n", global_rank);


    /*
     * call sub_method
     */
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_end_iteration_fn
       )
    {
        adios_transports[submethod->m].adios_end_iteration_fn
                                   (submethod);
    }

    return;
}

void adios_nssi_filter_stop_calculation(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int remote_rc=NSSI_OK;
    int myrank;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_stop_calc\n", global_rank);

    /*
     * call sub_method
     */
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_stop_calculation_fn
       )
    {
        adios_transports[submethod->m].adios_stop_calculation_fn
                                   (submethod);
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

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_write\n", global_rank);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "file is not open.  FAIL.");
        return;
    }
    md=of->md;

    /*
     * call sub_method
     */
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_write_fn
       )
    {
        adios_transports[submethod->m].adios_write_fn
                                   (f, v, data, submethod);
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
        fprintf(stderr, "file is not open.  FAIL.");
        return;
    }
    md=of->md;

    /*
     * call sub_method
     */
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_read_fn
       )
    {
        adios_transports[submethod->m].adios_read_fn
                                   (f, v, buffer, buffersize, submethod);
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
        fprintf(stderr, "file is not open.  FAIL.");
        return;
    }
    md=of->md;
    myrank=md->rank;

    if (DEBUG>3) printf("myrank(%d) enter adios_nssi_filter_close\n", myrank);

    /*
     * call sub_method
     */
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_close_fn
       )
    {
        adios_transports[submethod->m].adios_close_fn
                                   (f, submethod);
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

    open_file_delete(method->base_path, f->name);

    md->group_comm = MPI_COMM_NULL;
    md->fd = -1;
    md->rank = -1;
    md->size = 0;

    if (DEBUG>3) printf("global_rank(%d) exit adios_nssi_filter_close\n", global_rank);

    return;
}

void adios_nssi_filter_finalize(
        int mype,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    int myrank;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_filter_finalize\n", global_rank);

    /*
     * call sub_method
     */
    if (   submethod->m != ADIOS_METHOD_UNKNOWN
        && submethod->m != ADIOS_METHOD_NULL
        && adios_transports[submethod->m].adios_finalize_fn
       )
    {
        adios_transports[submethod->m].adios_finalize_fn
                                   (mype, submethod);
    }

    free_nssi_config(&nssi_cfg);

    if (adios_nssi_filter_initialized) {
        adios_nssi_filter_initialized = 0;
        self=NULL;
    }
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
        fprintf(stderr, "file is not open.  FAIL.\n");
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
        fprintf(stderr, "file is not open.  FAIL.\n");
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

#endif
