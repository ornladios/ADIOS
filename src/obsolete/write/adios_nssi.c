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
#include "util.h"
#ifdef HAVE_NSSI
#include "nssi_client.h"
#include "adios_nssi_args.h"
#include "adios_nssi_config.h"
#include "nssi_logger.h"
#endif

#include "io_timer.h"


#define NUM_GP 24
void adios_nssi_get_write_buffer(
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
struct adios_nssi_file_data_struct
{
    int      fd;
    MPI_Comm group_comm;
    int      size;
    int      rank;

    int      svc_index;
    MPI_Comm collective_op_comm;
    int      collective_op_size;
    int      collective_op_rank;

    int8_t use_single_server;
};

struct adios_nssi_method_data_struct
{
    nssi_request start_calc_req;
    int          has_outstanding_req;
    int          default_svc_index;  /* service to use when there is no open file (eg. finalize) */
};

/* Need a struct to encapsulate var offset info.
 */
struct var_offset {
    struct adios_group_struct *ogroup;
    uint16_t  oid;
    char      opath[ADIOS_PATH_MAX];
    char      oname[ADIOS_PATH_MAX];
    void     *ovalue;
    uint64_t  osize;
    uint8_t   is_anonymous;
};

/* list of variable offsets */
static List var_offset_list;

/* Need a struct to encapsulate var dim info.
 */
struct var_dim {
    struct adios_group_struct *dgroup;
    uint16_t  did;
    char      dpath[ADIOS_PATH_MAX];
    char      dname[ADIOS_PATH_MAX];
    void     *dvalue;
    uint64_t  dsize;
    uint8_t   is_anonymous;
};

/* list of variable offsets */
static List var_dim_list;

/* Need a struct to encapsulate open file info
 */
struct open_file {
    char                           fpath[ADIOS_PATH_MAX];
    char                           fname[ADIOS_PATH_MAX];
    struct adios_nssi_file_data_struct *file_data;
    struct adios_file_struct      *f;
    nssi_request                   start_calc_req;
    int                            has_outstanding_req;
//    List                           outstanding_reqs;
};

/* list of variable offsets */
static List open_file_list;

///////////////////////////
// Global Variables
///////////////////////////
static int adios_nssi_initialized = 0;

static char *job_id=NULL;
static int global_rank=-1;
static int global_size=-1;
static nssi_service *svcs;
struct adios_nssi_config nssi_cfg;

//static log_level adios_nssi_debug_level;
static int DEBUG=0;


///////////////////////////
// Function Declarations
///////////////////////////


///////////////////////////
// Function Definitions
///////////////////////////

static struct open_file *open_file_create(const char *path, const char *name, struct adios_nssi_file_data_struct *file_private_data)
{
    struct open_file *of=calloc(1,sizeof(struct open_file));

    strcpy(of->fpath, path);
    strcpy(of->fname, name);
    of->file_data = file_private_data;
    of->has_outstanding_req=FALSE;
//    list_init(&(of->outstanding_reqs), free);

    return(of);
}
static void open_file_free(void *of)
{
//    list_destroy((List *)&(of->outstanding_reqs));
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

static struct var_offset *var_offset_create(const char *path, const char *name, void *value, uint64_t size, uint16_t id, struct adios_group_struct *group)
{
    struct var_offset *vo=calloc(1,sizeof(struct var_offset));

    vo->ogroup = group;
    vo->oid = id;
    strcpy(vo->opath, path);
    strcpy(vo->oname, name);
    vo->ovalue = value;
    vo->osize  = size;

    return(vo);
}
static void var_offset_free(void *vo)
{
    free(vo);
}
static int var_offset_equal(const struct var_offset *vo1, const struct var_offset *vo2)
{
    if ((strcmp(vo1->opath, vo2->opath) == 0) && (strcmp(vo1->oname, vo2->oname) == 0)) return TRUE;

    return FALSE;
}
static struct var_offset *var_offset_find(struct adios_group_struct *group, const char *path, const char *name)
{
    ListElmt *elmt;
    struct var_offset *vo;

    if (DEBUG>2) printf("looking for group(%p) opath(%s) oname(%s)\n", group, path, name);

    elmt = list_head(&var_offset_list);
    while(elmt) {
        vo = list_data(elmt);
        if (DEBUG>2) printf("comparing to group(%p) opath(%s) oname(%s)\n", vo->ogroup, vo->opath, vo->oname);
        if ((group == vo->ogroup) && (strcmp(path, vo->opath) == 0) && (strcmp(name, vo->oname) == 0)) {
            if (DEBUG>2) printf("group(%p) opath(%s) oname(%s) matches search\n", vo->ogroup, vo->opath, vo->oname);
            return vo;
        }
        elmt = list_next(elmt);
    }

    return NULL;
}
static struct var_offset *var_offset_find_by_id(const uint16_t id)
{
    ListElmt *elmt;
    struct var_offset *vo;

    if (DEBUG>2) printf("looking for oid(%d)\n", id);

    elmt = list_head(&var_offset_list);
    while(elmt) {
        vo = list_data(elmt);
        if (DEBUG>2) printf("comparing to oid(%d)\n", vo->oid);
        if (vo->oid == id) {
            if (DEBUG>2) printf("oid(%d) opath(%s) oname(%s) matches search\n", vo->oid, vo->opath, vo->oname);
            return vo;
        }
        elmt = list_next(elmt);
    }

    return NULL;
}
static void var_offset_printall(void)
{
    ListElmt *elmt;
    struct var_offset *vo;

    elmt = list_head(&var_offset_list);
    while(elmt) {
        vo = list_data(elmt);
        if (DEBUG>3) printf("opath(%s) oname(%s)\n", vo->opath, vo->oname);
        elmt = list_next(elmt);
    }
}

static struct var_dim *var_dim_create(const char *path, const char *name, void *value, uint64_t size, uint16_t id, struct adios_group_struct *group)
{
    struct var_dim *vd=calloc(1,sizeof(struct var_dim));

    vd->dgroup = group;
    vd->did = id;
    strcpy(vd->dpath, path);
    strcpy(vd->dname, name);
    vd->dvalue = value;
    vd->dsize  = size;

    return(vd);
}
static void var_dim_free(void *vd)
{
    free(vd);
}
static int var_dim_equal(const struct var_dim *vd1, const struct var_dim *vd2)
{
    if ((strcmp(vd1->dpath, vd2->dpath) == 0) && (strcmp(vd1->dname, vd2->dname) == 0)) return TRUE;

    return FALSE;
}
static struct var_dim *var_dim_find(struct adios_group_struct *group, const char *path, const char *name)
{
    ListElmt *elmt;
    struct var_dim *vd;

    if (DEBUG>2) printf("looking for group(%p) dpath(%s) dname(%s)\n", group, path, name);

    elmt = list_head(&var_dim_list);
    while(elmt) {
        vd = list_data(elmt);
        if (DEBUG>2) printf("comparing to group(%p) dpath(%s) dname(%s)\n", vd->dgroup, vd->dpath, vd->dname);
        if ((group == vd->dgroup) && (strcmp(path, vd->dpath) == 0) && (strcmp(name, vd->dname) == 0)) {
            if (DEBUG>2) printf(" group(%p) dpath(%s) dname(%s) matches search\n", vd->dgroup, vd->dpath, vd->dname);
            return vd;
        }
        elmt = list_next(elmt);
    }

    return NULL;
}
static struct var_dim *var_dim_find_by_id(const uint16_t id)
{
    ListElmt *elmt;
    struct var_dim *vd;

    if (DEBUG>2) printf("looking for did(%d)\n", id);

    elmt = list_head(&var_dim_list);
    while(elmt) {
        vd = list_data(elmt);
        if (DEBUG>2) printf("comparing to did(%d)\n", vd->did);
        if (vd->did == id) {
            if (DEBUG>2) printf("did(%d) dpath(%s) dname(%s) matches search\n", vd->did, vd->dpath, vd->dname);
            return vd;
        }
        elmt = list_next(elmt);
    }

    return NULL;
}
static void var_dim_printall(void)
{
    ListElmt *elmt;
    struct var_dim *vd;

    elmt = list_head(&var_dim_list);
    while(elmt) {
        vd = list_data(elmt);
        if (DEBUG>3) printf("dpath(%s) dname(%s)\n", vd->dpath, vd->dname);
        elmt = list_next(elmt);
    }
}

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


static int gen_offset_list(
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root)
{
    struct adios_var_struct *v;
    struct adios_dimension_struct *dims;
    struct var_info *vi;
    char offset_name[255];
    uint64_t *value;

    v = pvar_root;
    while (v) {
        dims=v->dimensions;
        int loffs_idx=0;
        while (dims) {
            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, offset_name);
            if (offset_name[0] == '\0') {
                sprintf(offset_name, "offset_%d", loffs_idx);
            }
            if (DEBUG>3) printf("gen: offset_name(%s)\n", offset_name);
            value=calloc(1,sizeof(uint64_t));
            parse_dimension_size(group, pvar_root, patt_root, &dims->local_offset, value);

            uint64_t vsize = 4; /* adios_get_var_size(v, group, value); */
            vi = var_offset_create(v->path, offset_name, value, vsize, dims->local_offset.id, group);
            list_ins_next(&var_offset_list, list_tail(&var_offset_list), vi);

            loffs_idx++;
            dims = dims->next;
        }
        v = v->next;
    }
}

static void create_offset_list_for_var(
        struct adios_var_struct *v,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        u_int *offset_count,
        struct adios_var **offset_values)
{
    struct adios_dimension_struct *dims;
    char offset_name[255];

    *offset_count=0;
    *offset_values=NULL;

    if ((v) && (v->dimensions)) {
        int local_offset_count=0;
        dims=v->dimensions;
        while (dims) {
            if (dims->dimension.time_index == adios_flag_yes) {
                dims = dims->next;
                continue;
            }
            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, offset_name);
            local_offset_count++;
            dims = dims->next;
        }

        *offset_count=local_offset_count;
        *offset_values=calloc(local_offset_count, sizeof(struct adios_var));

        dims=v->dimensions;
        int loffs_idx=0;
        while (dims) {
            if (dims->dimension.time_index == adios_flag_yes) {
                dims = dims->next;
                continue;
            }
            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, offset_name);
            if (offset_name[0] == '\0') {
                sprintf(offset_name, "%s_%s_offset_%d", /*v->path*/"", v->name, loffs_idx);
                (*offset_values)[loffs_idx].is_anonymous=TRUE;
            }
            (*offset_values)[loffs_idx].vpath=strdup(v->path);
            (*offset_values)[loffs_idx].vname=strdup(offset_name);
//            struct var_offset *vo=var_offset_find("", offset_name);
//            memcpy(&(args->offsets.offsets_val[loffs_idx].vdata), vo->ovalue, vo->osize);
//            args->offsets.offsets_val[loffs_idx].vdatasize=vo->osize;
//            printf("create: offset_name(%s) offset_value(%lu)\n", offset_name, vo->ovalue);
            uint64_t value=0;
            parse_dimension_size(group, pvar_root, patt_root, &dims->local_offset, &value);
            memcpy(&((*offset_values)[loffs_idx].vdata), &value, 4);
            (*offset_values)[loffs_idx].vdatasize=4;
            if (DEBUG>3) printf("create: offset_name(%s) offset_value(%lu)\n", offset_name, value);

            loffs_idx++;
            dims = dims->next;
        }
    }
}

static int gen_ldim_list(
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root)
{
    struct adios_var_struct *v;
    struct adios_dimension_struct *dims;
    struct var_info *vi;
    char dim_name[255];
    uint64_t *value;

    v = pvar_root;
    while (v) {
        dims=v->dimensions;
        int dim_idx=0;
        while (dims) {
            parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dim_name);
            if (dim_name[0] == '\0') {
                sprintf(dim_name, "dim_%d", dim_idx);
            }
            if (DEBUG>3) printf("gen: dim_name(%s)\n", dim_name);
            value=calloc(1,sizeof(uint64_t));
            parse_dimension_size(group, pvar_root, patt_root, &dims->dimension, value);
            if (global_rank==0) {
                if (DEBUG>3) printf(":o(%d)", *value);
            }

            uint64_t vsize = 4; /* adios_get_var_size(v, group, value); */
            vi = var_dim_create(v->path, dim_name, value, vsize, dims->dimension.id, group);
            list_ins_next(&var_dim_list, list_tail(&var_dim_list), vi);

            dim_idx++;
            dims = dims->next;
        }
        v = v->next;
    }
}

static void create_ldim_list_for_var(
        struct adios_var_struct *v,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        u_int *dim_count,
        struct adios_var **dim_values)
{
    struct adios_dimension_struct *dims;
    char dim_name[255];

    *dim_count=0;
    *dim_values=NULL;

    if ((v) && (v->dimensions)) {
        int local_dim_count=0;
        dims=v->dimensions;
        while (dims) {
            if (dims->dimension.time_index == adios_flag_yes) {
                dims = dims->next;
                continue;
            }
            parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dim_name);
            local_dim_count++;
            dims = dims->next;
        }

        *dim_count=local_dim_count;
        *dim_values=calloc(local_dim_count, sizeof(struct adios_var));

        dims=v->dimensions;
        int dim_idx=0;
        while (dims) {
            if (dims->dimension.time_index == adios_flag_yes) {
                dims = dims->next;
                continue;
            }
            parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dim_name);
            if (dim_name[0] == '\0') {
                sprintf(dim_name, "%s_%s_dim_%d", /*v->path*/"", v->name, dim_idx);
                (*dim_values)[dim_idx].is_anonymous=TRUE;
            }
            (*dim_values)[dim_idx].vpath=strdup(v->path);
            (*dim_values)[dim_idx].vname=strdup(dim_name);
//            struct var_dim *vd=var_dim_find("", dim_name);
//            memcpy(&(args->dims.dims_val[dim_idx].vdata), vd->dvalue, vd->dsize);
//            args->dims.dims_val[dim_idx].vdatasize=vd->dsize;
            uint64_t value=0;
            parse_dimension_size(group, pvar_root, patt_root, &dims->dimension, &value);
            memcpy(&((*dim_values)[dim_idx].vdata), &value, 4);
            (*dim_values)[dim_idx].vdatasize=4;
            if (DEBUG>3) printf("create: dim_name(%s) dvalue(%lu)\n", dim_name, value);
            if (global_rank==0) {
                if (DEBUG>3) printf(":o(%d)", value);
            }

            dim_idx++;
            dims = dims->next;
        }
    }
}

static void create_gdim_list_for_var(
        struct adios_var_struct *v,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        u_int *dim_count,
        struct adios_var **dim_values)
{
    struct adios_dimension_struct *dims;
    char dim_name[255];

    *dim_count=0;
    *dim_values=NULL;

    if ((v) && (v->dimensions)) {
        int global_dim_count=0;
        dims=v->dimensions;
        while (dims) {
            if (dims->dimension.time_index == adios_flag_yes) {
                dims = dims->next;
                continue;
            }
            parse_dimension_name(group, pvar_root, patt_root, &dims->global_dimension, dim_name);
            global_dim_count++;
            dims = dims->next;
        }

        *dim_count=global_dim_count;
        *dim_values=calloc(global_dim_count, sizeof(struct adios_var));

        dims=v->dimensions;
        int dim_idx=0;
        while (dims) {
            if (dims->dimension.time_index == adios_flag_yes) {
                dims = dims->next;
                continue;
            }
            parse_dimension_name(group, pvar_root, patt_root, &dims->global_dimension, dim_name);
            if (dim_name[0] == '\0') {
                sprintf(dim_name, "%s_%s_dim_%d", /*v->path*/"", v->name, dim_idx);
                (*dim_values)[dim_idx].is_anonymous=TRUE;
            }
            (*dim_values)[dim_idx].vpath=strdup(v->path);
            (*dim_values)[dim_idx].vname=strdup(dim_name);
//            struct var_dim *vd=var_dim_find("", dim_name);
//            memcpy(&(args->dims.dims_val[dim_idx].vdata), vd->dvalue, vd->dsize);
//            args->dims.dims_val[dim_idx].vdatasize=vd->dsize;
            uint64_t value=0;
            parse_dimension_size(group, pvar_root, patt_root, &dims->global_dimension, &value);
            memcpy(&((*dim_values)[dim_idx].vdata), &value, 4);
            (*dim_values)[dim_idx].vdatasize=4;
            if (DEBUG>3) printf("create: dim_name(%s) dvalue(%lu)\n", dim_name, value);
            if (global_rank==0) {
                if (DEBUG>3) printf(":o(%d)", value);
            }

            dim_idx++;
            dims = dims->next;
        }
    }
}

static int read_var(
        struct adios_nssi_file_data_struct *file_data,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_var_struct *pvar)
{
    int return_code=0;
    int i, rc;

    adios_read_args args;
    adios_read_res  res;

    memset(&args, 0, sizeof(args));
    memset(&res, 0, sizeof(res));
    args.fd    = file_data->fd;
    args.max_read = pvar->data_size;
    args.vpath = strdup(pvar->path);
    args.vname = strdup(pvar->name);


    if (pvar->dimensions) {
        args.is_scalar = FALSE;
    } else {
        args.is_scalar = TRUE;
    }

    args.offsets.offsets_len=0;
    args.offsets.offsets_val=NULL;
    args.ldims.ldims_len=0;
    args.ldims.ldims_val=NULL;
    args.gdims.gdims_len=0;
    args.gdims.gdims_val=NULL;
    if (pvar->dimensions) {
        create_offset_list_for_var(
                 pvar,
                 group,
                 group->vars,
                 group->attributes,
                 &args.offsets.offsets_len,
                 &args.offsets.offsets_val);
        create_ldim_list_for_var(
                 pvar,
                 group,
                 group->vars,
                 group->attributes,
                 &args.ldims.ldims_len,
                 &args.ldims.ldims_val);
        create_gdim_list_for_var(
                 pvar,
                 group,
                 group->vars,
                 group->attributes,
                 &args.gdims.gdims_len,
                 &args.gdims.gdims_val);
     }

    Func_Timer("ADIOS_READ_OP",
            rc = nssi_call_rpc_sync(&svcs[file_data->svc_index],
            ADIOS_READ_OP,
            &args,
            pvar->data,
            pvar->data_size,
            &res););
    if (rc != NSSI_OK) {
        fprintf(stderr, "NSSI ERROR: ADIOS_READ_OP failed\n");
        return_code=-2;
    }

    free(args.vpath);
    free(args.vname);
    if (pvar->dimensions) {
        for (i=0;i<args.offsets.offsets_len;i++) {
            free(args.offsets.offsets_val[i].vpath);
            free(args.offsets.offsets_val[i].vname);
        }
        free(args.offsets.offsets_val);
        for (i=0;i<args.ldims.ldims_len;i++) {
            free(args.ldims.ldims_val[i].vpath);
            free(args.ldims.ldims_val[i].vname);
        }
        free(args.ldims.ldims_val);
        for (i=0;i<args.gdims.gdims_len;i++) {
            free(args.gdims.gdims_val[i].vpath);
            free(args.gdims.gdims_val[i].vname);
        }
        free(args.gdims.gdims_val);
    }

    return return_code;
}
static int write_var(
        struct adios_nssi_file_data_struct *file_data,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_var_struct *pvar,
        uint64_t var_size,
        enum ADIOS_FLAG fortran_flag)
{
    int i;
    int rc;
    int return_code=0;

    adios_write_args args;
    adios_write_res  res;

//    var_offset_printall();

    memset(&args, 0, sizeof(args));
    memset(&res, 0, sizeof(res));
    args.fd    = file_data->fd;
    args.vpath = strdup(pvar->path);
    args.vname = strdup(pvar->name);
    args.vsize = var_size;
    args.atype = pvar->type;
    if (pvar->dimensions) {
        args.is_scalar = FALSE;
    } else {
        args.is_scalar = TRUE;
    }
    args.writer_rank=file_data->rank;
    args.offsets.offsets_len=0;
    args.offsets.offsets_val=NULL;
    args.ldims.ldims_len=0;
    args.ldims.ldims_val=NULL;
    args.gdims.gdims_len=0;
    args.gdims.gdims_val=NULL;
    if (pvar->dimensions) {
        create_offset_list_for_var(
                 pvar,
                 group,
                 group->vars,
                 group->attributes,
                 &args.offsets.offsets_len,
                 &args.offsets.offsets_val);
        create_ldim_list_for_var(
                 pvar,
                 group,
                 group->vars,
                 group->attributes,
                 &args.ldims.ldims_len,
                 &args.ldims.ldims_val);
        create_gdim_list_for_var(
                 pvar,
                 group,
                 group->vars,
                 group->attributes,
                 &args.gdims.gdims_len,
                 &args.gdims.gdims_val);
     }

    Func_Timer("ADIOS_WRITE_OP",
            rc = nssi_call_rpc_sync(&svcs[file_data->svc_index],
            ADIOS_WRITE_OP,
            &args,
            pvar->data,
            var_size,
            &res););
    if (rc != NSSI_OK) {
        fprintf(stderr, "NSSI ERROR: ADIOS_WRITE_OP failed\n");
        return_code=-2;
    }

    // looking for alternatives to this barrier.
    // if variable writes on the staging server are collective (NC4), then clients must stay in sync.
    // one solution is to cache scalar writes in addition to array writes.
    MPI_Barrier(file_data->group_comm);

    free(args.vpath);
    free(args.vname);
    if (pvar->dimensions) {
        for (i=0;i<args.offsets.offsets_len;i++) {
            free(args.offsets.offsets_val[i].vpath);
            free(args.offsets.offsets_val[i].vname);
        }
        free(args.offsets.offsets_val);
        for (i=0;i<args.ldims.ldims_len;i++) {
            free(args.ldims.ldims_val[i].vpath);
            free(args.ldims.ldims_val[i].vname);
        }
        free(args.ldims.ldims_val);
        for (i=0;i<args.gdims.gdims_len;i++) {
            free(args.gdims.gdims_val[i].vpath);
            free(args.gdims.gdims_val[i].vname);
        }
        free(args.gdims.gdims_val);
    }

    return return_code;
}


void adios_nssi_init(
        const PairStruct *parameters,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    int verbose=5;
    char logfile[1024];
    int log_rank;
    struct adios_nssi_method_data_struct *private;

    if (!adios_nssi_initialized) {
        adios_nssi_initialized = 1;
    }

    MPI_Comm_rank(method->init_comm, &global_rank);
    MPI_Comm_size(method->init_comm, &global_size);

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_init\n", global_rank);

//    MPI_Comm_rank(file_data->group_comm, &log_rank);
//    sprintf(logfile, "%s.%04d", "adios_nssi_client.log", log_rank);
//    logger_init((log_level)verbose, logfile);

//    logger_init((log_level)verbose, NULL);

#ifdef HAVE_PORTALS
    nssi_ptl_init(PTL_IFACE_CLIENT, getpid() + 1000);
    nssi_rpc_init(NSSI_RPC_PTL, NSSI_RPC_XDR);
#endif
#ifdef HAVE_INFINIBAND
    nssi_ib_init(NULL);
    rc = nssi_rpc_init(NSSI_RPC_IB, NSSI_RPC_XDR);
#endif

    /* Register the client operations */
    NSSI_REGISTER_CLIENT_STUB(ADIOS_OPEN_OP, adios_open_args, void, adios_open_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_GROUP_SIZE_OP, adios_group_size_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_READ_OP, adios_read_args, void, adios_read_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_WRITE_OP, adios_write_args, void, adios_write_res);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_END_ITER_OP, adios_end_iter_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_START_CALC_OP, adios_start_calc_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_STOP_CALC_OP, adios_stop_calc_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_CLOSE_OP, adios_close_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(ADIOS_FINALIZE_OP, adios_finalize_args, void, void);

    list_init(&open_file_list, open_file_free);
    list_init(&var_offset_list, var_offset_free);
    list_init(&var_dim_list, var_dim_free);

    private = malloc(sizeof(struct adios_nssi_method_data_struct));
    private->has_outstanding_req=FALSE;
    private->default_svc_index=-1;
    method->method_data = private;

    parse_nssi_config(getenv("ADIOS_NSSI_CONFIG_FILE"), &nssi_cfg);

    return;
}


enum ADIOS_FLAG adios_nssi_should_buffer(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;

    adios_group_size_args args;

    struct open_file *of=NULL;
    struct adios_nssi_file_data_struct *file_data=NULL;

    uint64_t max_data_size=0;

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "file is not open.  FAIL.\n");
        return adios_flag_no;
    }
    file_data=of->file_data;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_should_buffer (write_size_bytes=%lu)\n", global_rank, f->write_size_bytes);

    memset(&args, 0, sizeof(args));
    args.fd = file_data->fd;
    MPI_Reduce(
            &f->write_size_bytes,
            &max_data_size,
            1,
            MPI_UNSIGNED_LONG,
            MPI_SUM,
            0,
            file_data->collective_op_comm);
    if (file_data->collective_op_rank == 0) {
        if (DEBUG > 3) printf("max_data_size==%lu\n", max_data_size);
        args.data_size = max_data_size*file_data->collective_op_size;
        args.data_size = max_data_size;
        if (DEBUG > 3) printf("args.data_size==%lu\n", args.data_size);

        Func_Timer("ADIOS_GROUP_SIZE_OP",
                rc = nssi_call_rpc_sync(&svcs[file_data->svc_index],
                        ADIOS_GROUP_SIZE_OP,
                        &args,
                        NULL,
                        0,
                        NULL););
        if (rc != NSSI_OK) {
            fprintf(stderr, "NSSI ERROR: ADIOS_GROUP_SIZE_OP failed\n");
        }
    }
    // this barrier ensures that none of the clients race ahead and
    // attempt writes before the file is open and prepped on the staging server.
    MPI_Barrier(file_data->group_comm);

    return adios_flag_no;
}

int adios_nssi_open(
        struct adios_file_struct *f,
        struct adios_method_struct *method,
        MPI_Comm comm)
{
    int rc=NSSI_OK;

    adios_open_args args;
    adios_open_res  res;

    struct open_file *of=NULL;
    struct adios_nssi_file_data_struct *file_data=NULL;
    struct adios_nssi_method_data_struct *private=(struct adios_nssi_method_data_struct *)method->method_data;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_open\n", global_rank);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        file_data             = malloc(sizeof(struct adios_nssi_file_data_struct));
        file_data->fd         = -1;
        file_data->group_comm = MPI_COMM_NULL;
        file_data->size       = 0;
        file_data->rank       = 0;

        file_data->svc_index=-1;
        file_data->collective_op_comm=MPI_COMM_NULL;
        file_data->collective_op_size=0;
        file_data->collective_op_rank=-1;

        file_data->use_single_server=FALSE;

        of=open_file_create(method->base_path, f->name, file_data);
    } else {
        // sanity check
        if (file_data->fd == -1) {
            if (DEBUG>3) printf("open: %s is open but fd==-1.  sanity check failed.  attempting reopen.\n", f->name);
            open_file_delete(of->fpath, of->fname);
        } else {
            // file already open
            return adios_flag_no;
        }
    }

    if (DEBUG>3) printf("global_rank(%d): enter adios_nssi_open (%s)\n", global_rank, f->name);

    file_data->group_comm = comm;
    if (DEBUG>3) printf("global_rank(%d): adios_nssi_open: setup group_comm\n", global_rank);
    if (file_data->group_comm != MPI_COMM_NULL) {
        if (DEBUG>3) printf("global_rank(%d): adios_nssi_open: get rank and size\n", global_rank);
        MPI_Comm_rank(file_data->group_comm, &file_data->rank);
        MPI_Comm_size(file_data->group_comm, &file_data->size);
        if (DEBUG>3) printf("global_rank(%d): adios_nssi_open: size(%d) rank(%d)\n", global_rank, file_data->size, file_data->rank);
    } else {
        file_data->group_comm=MPI_COMM_SELF;
    }
    f->group->process_id = file_data->rank;

    if (file_data->size <= nssi_cfg.num_servers) {
        // there are fewer clients than servers.
        // assume file-per-process and use a single server for this file.
        file_data->use_single_server=TRUE;
        if (file_data->size < global_size) {
            // a subset of all clients is writing
            file_data->svc_index = ((global_rank/file_data->size)%nssi_cfg.num_servers);
        } else {
            file_data->svc_index = 0;
        }
    } else {
        file_data->use_single_server=FALSE;
        if ((file_data->size%nssi_cfg.num_servers) > 0) {
            file_data->svc_index = file_data->rank/((file_data->size/nssi_cfg.num_servers)+1);
        } else {
            file_data->svc_index = file_data->rank/(file_data->size/nssi_cfg.num_servers);
        }
    }
    if (private->default_svc_index == -1) {
        private->default_svc_index=file_data->svc_index;
    }

    /* create a new communicator for just those clients, who share a default service. */
    if (DEBUG>3) printf("global_rank(%d): adios_nssi_open: before MPI_Comm_split\n", global_rank);
    MPI_Comm_split(file_data->group_comm, file_data->svc_index, file_data->rank, &file_data->collective_op_comm);
    if (DEBUG>3) printf("global_rank(%d): adios_nssi_open: after MPI_Comm_split\n", global_rank);
    /* find my rank in the new communicator */
    if (DEBUG>3) printf("global_rank(%d): adios_nssi_open: before MPI_Comm_size\n", global_rank);
    MPI_Comm_size(file_data->collective_op_comm, &file_data->collective_op_size);
    if (DEBUG>3) printf("global_rank(%d): adios_nssi_open: before MPI_Comm_rank\n", global_rank);
    MPI_Comm_rank(file_data->collective_op_comm, &file_data->collective_op_rank);

    if (DEBUG>3) printf("global_rank(%d) file_data->rank(%d) file_data->collective_op_rank(%d) default_service(%d)\n", global_rank, file_data->rank, file_data->collective_op_rank, file_data->svc_index);

    svcs=(nssi_service *)calloc(nssi_cfg.num_servers, sizeof(nssi_service));
    /* !global_rank0 has a preferred server for data transfers.  connect to preferred server.
     * connect to other servers on-demand.
     */
    double GetSvcTime=MPI_Wtime();
    if (DEBUG>3) printf("get staging-service: file_data->svc_index(%d) nid(%lld) pid(%llu) hostname(%s) port(%d)\n",
            file_data->svc_index,
            nssi_cfg.nssi_server_ids[file_data->svc_index].nid,
            nssi_cfg.nssi_server_ids[file_data->svc_index].pid,
            nssi_cfg.nssi_server_ids[file_data->svc_index].hostname,
            nssi_cfg.nssi_server_ids[file_data->svc_index].port);
    rc = nssi_get_service(nssi_cfg.nssi_server_ids[file_data->svc_index], -1, &svcs[file_data->svc_index]);
    if (rc != NSSI_OK) {
        fprintf(stderr, "NSSI ERROR: nssi_get_service failed\n");
        return;
    }


    if (job_id==NULL) {
        if (file_data->rank==0) {
            job_id = getenv("PBS_JOBID");
            if (job_id == NULL) {
                fprintf(stderr, "adios_nssi_open: unable to determine job id.  defaulting id to \"UNKNOWN_JOB_ID\".\n");
                job_id = strdup("UNKNOWN_JOB_ID");
            } else {
                int len=strlen(job_id)+36+1;
                job_id=calloc(len,1);

                struct uuid_st;
                extern int uuid_create   (      struct uuid_st **_uuid);
                extern int uuid_destroy  (      struct uuid_st  *_uuid);
                extern int uuid_make     (      struct uuid_st  *_uuid, unsigned int _mode, ...);
                extern int uuid_export   (const struct uuid_st  *_uuid, unsigned int _fmt,       void **_data_ptr, size_t *_data_len);

                struct uuid_st *uuid;
                char *uuid_str=NULL;
                uuid_create(&uuid);
                uuid_make(uuid, 1);
                uuid_export(uuid, 1, &uuid_str, NULL);
                uuid_destroy(uuid);

                sprintf(job_id, "%s.%s", getenv("PBS_JOBID"), uuid_str);

                free(uuid_str);

                MPI_Bcast (job_id, len, MPI_BYTE, 0, file_data->group_comm);
            }
        } else {
            job_id = getenv("PBS_JOBID");
            if (job_id == NULL) {
                fprintf(stderr, "adios_nssi_open: unable to determine job id.  defaulting id to \"UNKNOWN_JOB_ID\".\n");
                job_id = strdup("UNKNOWN_JOB_ID");
            } else {
                int len=strlen(job_id)+36+1;
                job_id=calloc(len,1);

                MPI_Bcast (job_id, len, MPI_BYTE, 0, file_data->group_comm);
            }
        }
    }

    fprintf(stderr, "adios_nssi_open: job_id=%s\n", job_id);




    gen_offset_list(
            f->group,
            f->group->vars,
            f->group->attributes);
    gen_ldim_list(
            f->group,
            f->group->vars,
            f->group->attributes);
    var_offset_printall();

    memset(&args, 0, sizeof(args));
    memset(&res, 0, sizeof(res));
    args.client_id = strdup(job_id);
    args.fname = malloc(sizeof(char) * (strlen(method->base_path) + strlen(f->name) + 1));
    sprintf(args.fname, "%s%s", method->base_path, f->name);
    args.gname = strdup(method->group->name);
    switch(f->mode) {
        case adios_mode_read:
            args.mode = ADIOS_MODE_READ;
            break;
        case adios_mode_write:
            args.mode = ADIOS_MODE_WRITE;
            break;
        case adios_mode_append:
            args.mode = ADIOS_MODE_APPEND;
            break;
        case adios_mode_update:
            args.mode = ADIOS_MODE_UPDATE;
            break;
    }
    args.use_single_server=file_data->use_single_server;

    Func_Timer("ADIOS_OPEN_OP",
            rc = nssi_call_rpc_sync(&svcs[file_data->svc_index],
            ADIOS_OPEN_OP,
            &args,
            NULL,
            0,
            &res););
    if (rc != NSSI_OK) {
        fprintf(stderr, "NSSI ERROR: ADIOS_OPEN_OP failed\n");
    }

    file_data->fd = res.fd;

    list_ins_next(&open_file_list, list_tail(&open_file_list), of);

    free(args.fname);
    free(args.gname);

    return 1;
}

void adios_nssi_start_calculation(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int myrank;

    adios_start_calc_args args;

    ListElmt *of_elmt;
    struct open_file *of=NULL;
    struct adios_nssi_file_data_struct *file_data=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_start_calc\n", global_rank);

    of_elmt = list_head(&open_file_list);
    while(of_elmt) {
        of = list_data(of_elmt);
        if (of == NULL) {
            fprintf(stderr, "file is not open.  FAIL.");
            return;
        }
        file_data=of->file_data;
        myrank=file_data->rank;

        if (of->has_outstanding_req == FALSE) {
            if (file_data->collective_op_rank == 0) {
                memset(&args, 0, sizeof(args));
                args.fd = file_data->fd;
                Func_Timer("ADIOS_START_CALC_OP",
                        rc = nssi_call_rpc(&svcs[file_data->svc_index],
                                ADIOS_START_CALC_OP,
                                &args,
                                NULL,
                                0,
                                NULL,
                                &of->start_calc_req););
                if (rc != NSSI_OK) {
                    fprintf(stderr, "NSSI ERROR: ADIOS_START_CALC_OP failed\n");
                } else {
                    of->has_outstanding_req=TRUE;
                }
            }
        } else {
            if (DEBUG>3) fprintf(stderr, "rank(%d) has an outstanding ADIOS_START_CALC_OP.  skipping.\n", global_rank);
        }

        of_elmt = list_next(of_elmt);
    }

    if (DEBUG>3) printf("rank(%d) exit adios_nssi_start_calc\n", global_rank);

    return;
}

void adios_nssi_end_iteration(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int myrank;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_end_iter\n", global_rank);


//    adios_end_iter_args args;
//
//    struct open_file *of=NULL;
//    struct adios_nssi_file_data_struct *file_data=NULL;
//
//    of=open_file_find(method->base_path, f->name);
//    if (of == NULL) {
//        fprintf(stderr, "file is not open.  FAIL.");
//        return;
//    }
//    file_data=of->file_data;
//
//    myrank=file_data->rank;
//
//    if (file_data->collective_op_rank == 0) {
//        args.fd = file_data->fd;
//        Func_Timer("ADIOS_END_ITER_OP",
//                rc = nssi_call_rpc_sync(&svcs[file_data->svc_index],
//                ADIOS_END_ITER_OP,
//                &args,
//                NULL,
//                0,
//                NULL););
//        if (rc != NSSI_OK) {
//            //log_error(adios_nssi_debug_level, "unable to call remote adios_read");
//        }
//    }

    return;
}

void adios_nssi_stop_calculation(
        struct adios_method_struct * method)
{
    int rc=NSSI_OK;
    int remote_rc=NSSI_OK;
    int myrank;

    adios_stop_calc_args args;

    ListElmt *of_elmt;
    struct open_file *of=NULL;
    struct adios_nssi_file_data_struct *file_data=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_stop_calc\n", global_rank);


    of_elmt = list_head(&open_file_list);
    while(of_elmt) {
        of = list_data(of_elmt);
        if (of == NULL) {
            fprintf(stderr, "file is not open.  FAIL.");
            return;
        }
        file_data=of->file_data;
        myrank=file_data->rank;

        if (of->has_outstanding_req == TRUE) {
            // wait for any async writes to finish
            if (file_data->collective_op_rank == 0) {
                Func_Timer("ADIOS_START_CALC_OP nssi_wait",
                        rc=nssi_wait(&of->start_calc_req, &remote_rc););
                if (rc != NSSI_OK) {
                    fprintf(stderr, "NSSI ERROR: unable to wait for remote ADIOS_START_CALC_OP\n");
                }
                if (remote_rc != NSSI_OK) {
                    fprintf(stderr, "NSSI ERROR: ADIOS_START_CALC_OP failed\n");
                }
                of->has_outstanding_req = FALSE;
            }
        } else {
            if (DEBUG>3) fprintf(stderr, "rank(%d) has no outstanding ADIOS_START_CALC_OP.  skipping.\n", global_rank);
        }

        of_elmt = list_next(of_elmt);
    }


    of_elmt = list_head(&open_file_list);
    while(of_elmt) {
        of = list_data(of_elmt);
        if (of == NULL) {
            fprintf(stderr, "file is not open.  FAIL.");
            return;
        }
        file_data=of->file_data;
        myrank=file_data->rank;

        if (file_data->collective_op_rank == 0) {
            memset(&args, 0, sizeof(args));
            args.fd = file_data->fd;
            Func_Timer("ADIOS_STOP_CALC_OP",
                    rc = nssi_call_rpc_sync(&svcs[file_data->svc_index],
                            ADIOS_STOP_CALC_OP,
                            &args,
                            NULL,
                            0,
                            NULL););
            if (rc != NSSI_OK) {
                fprintf(stderr, "NSSI ERROR: ADIOS_STOP_CALC_OP failed\n");
            }
        }

        of_elmt = list_next(of_elmt);
    }

    if (DEBUG>3) printf("rank(%d) exit adios_nssi_stop_calc\n", global_rank);

    return;
}

void adios_nssi_write(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *data,
        struct adios_method_struct *method)
{
    static int first_write = 1;

    struct open_file *of=NULL;
    struct adios_nssi_file_data_struct *file_data=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_write\n", global_rank);

//    struct var_offset *vo = var_offset_find_by_id(v->id);
    struct var_offset *vo = var_offset_find(f->group, "/", v->name);
    if (vo != NULL) {
        if (DEBUG>2) printf("rank(%d) var(%s, %s) is an offset.  skip write.\n", global_rank, v->path, v->name);
        return;
    }
//    struct var_dim *vd = var_dim_find_by_id(v->id);
    struct var_dim *vd = var_dim_find(f->group, "/", v->name);
    if (vd != NULL) {
        if (DEBUG>2) printf("rank(%d) var(%s, %s) is a dimension.  skip write.\n", global_rank, v->path, v->name);
        return;
    }

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "file is not open.  FAIL.");
        return;
    }
    file_data=of->file_data;

    if (f->mode == adios_mode_write || f->mode == adios_mode_append) {

        if (file_data->rank==0) {
            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
            if (DEBUG>3) fprintf(stderr, "write var: %s start!\n", v->name);
        }
        uint64_t var_size = adios_get_var_size (v, data);
        if (DEBUG>2) printf("rank (%d) adios_nssi_write: vpath(%s) vname(%s) vsize(%ld)\n", global_rank, v->path, v->name, var_size);
        write_var(file_data,
                f->group,
                f->group->vars,
                f->group->attributes,
                v,
                var_size,
                f->group->adios_host_language_fortran);
    } else {
        if (DEBUG>3) fprintf(stderr, "entering unknown nc4 mode %d!\n", f->mode);
    }
    if (file_data->rank==0) {
        if (DEBUG>3) fprintf(stderr, "write var: %s end!\n", v->name);
        if (DEBUG>3) fprintf(stderr, "-------------------------\n");
    }

    return;
}


void adios_nssi_read(
        struct adios_file_struct *f,
        struct adios_var_struct *v,
        void *buffer,
        uint64_t buffersize,
        struct adios_method_struct *method)
{
    struct open_file *of=NULL;
    struct adios_nssi_file_data_struct *file_data=NULL;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_read\n", global_rank);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "file is not open.  FAIL.");
        return;
    }
    file_data=of->file_data;

    if(f->mode == adios_mode_read) {
        v->data = buffer;
        v->data_size = buffersize;

        if (v->is_dim == adios_flag_yes) {
            // this is a dimension variable.  values in the file are unreliable.
            // assume the caller provided a valid value in 'buffer'.
            if (DEBUG>3) fprintf(stderr, "------------------------------\n");
            if (DEBUG>3) fprintf(stderr, "read var: %s! skipping dim var\n", v->name);
            if (DEBUG>3) fprintf(stderr, "------------------------------\n");
            return;
        }

        if (file_data->rank==0) {
            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
            if (DEBUG>3) fprintf(stderr, "read var: %s! start\n", v->name);
        }
        if (DEBUG>3) printf("rank (%d) adios_nssi_read: vname(%s) vsize(%ld)\n", global_rank, v->name, v->data_size);
        read_var(file_data,
                f->group,
                f->group->vars,
                f->group->attributes,
                v);
        if (file_data->rank==0) {
            if (DEBUG>3) fprintf(stderr, "read var: %s! end\n", v->name);
            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
        }
    }

    return;
}

void adios_nssi_close(
        struct adios_file_struct *f,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    struct adios_attribute_struct * a = f->group->attributes;
    int myrank;

    adios_close_args close_args;

    struct open_file *of=NULL;
    struct adios_nssi_file_data_struct *file_data=NULL;

    if (DEBUG>3) printf("global_rank(%d) enter adios_nssi_close\n", global_rank);

    of=open_file_find(method->base_path, f->name);
    if (of == NULL) {
        fprintf(stderr, "file is not open.  FAIL.");
        return;
    }
    file_data=of->file_data;
    myrank=file_data->rank;

    if (DEBUG>3) printf("myrank(%d) enter adios_nssi_close\n", myrank);

    if (f->mode == adios_mode_read) {
        struct adios_var_struct * v = f->group->vars;
        while (v)
        {
            v->data = v->adata = 0;
            v = v->next;
        }

        if (file_data->rank==0) {
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
            if (DEBUG>1) fprintf(stderr, "reading done, NSSI file is virtually closed;\n");
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
        }
    } else if (f->mode == adios_mode_write || f->mode == adios_mode_append) {
        //fprintf(stderr, "entering nc4 write attribute mode!\n");
        if (file_data->rank==0) {
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
            if (DEBUG>1) fprintf(stderr, "writing done, NSSI file is virtually closed;\n");
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
        }
    }

//    // make sure all clients have finsihed I/O before closing
//    MPI_Barrier(file_data->group_comm);
    if (file_data->collective_op_rank == 0) {
        memset(&close_args, 0, sizeof(close_args));
        close_args.fname = malloc(sizeof(char) * (strlen(method->base_path) + strlen(f->name) + 1));
        sprintf(close_args.fname, "%s%s", method->base_path, f->name);
        close_args.fd = file_data->fd;
        if (DEBUG>3) printf("rank(%d) sending ADIOS_CLOSE_OP\n", myrank);
        Func_Timer("ADIOS_CLOSE_OP",
                rc = nssi_call_rpc_sync(&svcs[file_data->svc_index],
                ADIOS_CLOSE_OP,
                &close_args,
                NULL,
                0,
                NULL););
        if (rc != NSSI_OK) {
            fprintf(stderr, "NSSI ERROR: ADIOS_CLOSE_OP failed\n");
        }
        free(close_args.fname);
    }

    open_file_delete(method->base_path, f->name);

    file_data->group_comm = MPI_COMM_NULL;
    file_data->fd = -1;
    file_data->rank = -1;
    file_data->size = 0;

    if (DEBUG>3) printf("global_rank(%d) exit adios_nssi_close\n", global_rank);

    return;
}

void adios_nssi_finalize(
        int mype,
        struct adios_method_struct *method)
{
    int rc=NSSI_OK;
    int myrank;
    adios_finalize_args finalize_args;
    struct adios_nssi_method_data_struct *private=(struct adios_nssi_method_data_struct *)method->method_data;

    if (DEBUG>3) printf("rank(%d) enter adios_nssi_finalize\n", global_rank);

    if (global_rank == 0) {
        if (job_id==NULL) {
            if (DEBUG>2) fprintf(stderr, "adios_nssi_finalize - adios_nssi_open was never called.  skipping finalize.\n");
        } else {
            memset(&finalize_args, 0, sizeof(finalize_args));
            finalize_args.client_id = strdup(job_id);
            if (DEBUG>3) printf("rank(%d) sending ADIOS_FINALIZE_OP\n", myrank);
            Func_Timer("ADIOS_FINALIZE_OP",
                    rc = nssi_call_rpc_sync(&svcs[private->default_svc_index],
                            ADIOS_FINALIZE_OP,
                            &finalize_args,
                            NULL,
                            0,
                            NULL););
            if (rc != NSSI_OK) {
                fprintf(stderr, "NSSI ERROR: ADIOS_FINALIZE_OP failed\n");
            }
            free(finalize_args.client_id);
        }
    }

    free_nssi_config(&nssi_cfg);

    if (adios_nssi_initialized)
        adios_nssi_initialized = 0;
}
