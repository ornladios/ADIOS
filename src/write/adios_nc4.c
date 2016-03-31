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

#include "mpi.h"
#include "public/adios_types.h"
#include "core/adios_bp_v1.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_internals.h"
#include "core/util.h"
#include "netcdf_par.h"
#include "netcdf.h"

#include "nssi/io_timer.h"


#ifndef FALSE
#      define  FALSE   (0)
#endif
#ifndef TRUE
#      define  TRUE    (1)
#endif
#ifndef NULL
#      define  NULL 0
#endif


typedef char nc4_dimname_t[256];

#define NUM_GP 24
void adios_nc4_end_iteration(
        struct adios_method_struct * method)
{
}
void adios_nc4_start_calculation(
        struct adios_method_struct * method)
{
}
void adios_nc4_stop_calculation(
        struct adios_method_struct * method)
{
}
void adios_nc4_get_write_buffer(
        struct adios_file_struct *fd,
        struct adios_var_struct *v,
        uint64_t *size,
        void **buffer,
        struct adios_method_struct *method)
{
}

typedef struct {
    struct adios_dimension_struct *dims;

    nc4_dimname_t  gbdims_dim0_name;
    nc4_dimname_t  gbdims_dim1_name;
    nc4_dimname_t  gbdims_name;
    nc4_dimname_t *nc4_global_dimnames;
    nc4_dimname_t *nc4_local_dimnames;
    nc4_dimname_t *nc4_local_offset_names;

    size_t    *nc4_gbdims;
    size_t    *nc4_globaldims;
    size_t    *nc4_localdims;
    size_t    *nc4_offsets;
    ptrdiff_t *nc4_strides;
    int       *nc4_global_dimids;
    int       *nc4_local_dimids;
    int       *nc4_loffs_dimids;

    ptrdiff_t nc4_gbstrides[2];
    size_t    nc4_gbglobaldims[2];
    size_t    nc4_gblocaldims[2];
    size_t    nc4_gboffsets[2];
    int       nc4_gbglobaldims_dimids[2];
    int       nc4_gbglobaldims_varid;

    enum ADIOS_FLAG has_globaldims;
    enum ADIOS_FLAG has_localdims;
    enum ADIOS_FLAG has_localoffsets;

    enum ADIOS_FLAG has_timedim;
    int             timedim_index;

    int global_dim_count;
    int local_dim_count;
    int local_offset_count;
} deciphered_dims_t;
struct adios_nc4_data_struct
{
    int      fd;
    int      ncid;
    int      root_ncid;
    MPI_Comm group_comm;
    int      rank;
    int      size;
};

#define NC4_PATH_MAX 1024
/* Need a struct to encapsulate open file info
 */
struct open_file {
    char                           fpath[NC4_PATH_MAX];
    char                           fname[NC4_PATH_MAX];
    struct adios_nc4_data_struct *md;
    struct adios_file_struct      *f;
};

/* list of variable offsets */
static List open_file_list;


static int global_rank=-1;
static int DEBUG=0;



///////////////////////////
// Function Declarations
///////////////////////////
int getTypeSize(
        enum ADIOS_DATATYPES type,
        void *val);

// adios_flag determine whether it is dataset or group

int ncd_gen_name(
        char *fullname,
        char *path,
        char *name);
int getNC4TypeId(
        enum ADIOS_DATATYPES type,
        nc_type *nc4_type_id,
        enum ADIOS_FLAG fortran_flag);




static struct open_file *open_file_create(
        const char *path,
        const char *name,
        struct adios_nc4_data_struct *method_private_data,
        struct adios_file_struct *f)
{
    struct open_file *of=calloc(1,sizeof(struct open_file));

    strcpy(of->fpath, path);
    strcpy(of->fname, name);
    of->md = method_private_data;
    of->f = f;


    return(of);
}
static void open_file_free(void *of)
{
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
            if (DEBUG>3) printf("fpath(%s) fname(%s) ncid(%d) matches search\n", of->fpath, of->fname, of->md->ncid);
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
        if (DEBUG>3) printf("fpath(%s) fname(%s) ncid(%d)\n", of->fpath, of->fname, of->md->ncid);
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
    struct adios_attribute_struct *attr_linked;
    if (dim->var) {
        if (dim->var->data) {
            *dimsize = *(int *)dim->var->data;
        }
    }
    else if (dim->attr) 
    {
        attr_linked = dim->attr;
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
                            ,adios_type_to_string_int (attr_linked->type)
                            );
                    break;
            }
        } else {
            if (attr_linked->var->data) {
                *dimsize = *(int *)attr_linked->var->data;
            }
        }
    } else {
        if (dim->is_time_index == adios_flag_yes) {
            *dimsize = NC_UNLIMITED;
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
        nc4_dimname_t dimname)
{
    struct adios_var_struct *var_linked = NULL;
    struct adios_attribute_struct *attr_linked;
    if (dim->var) {
        sprintf(dimname, "%s_dim", dim->var->name);
    } else if (dim->attr) {
        if (!dim->attr->var) {
            sprintf(dimname, "%s_dim", dim->attr->name);
        } else {
            sprintf(dimname, "%s_dim", dim->attr->var->name);
        }
    } else {
        if (dim->is_time_index == adios_flag_yes) {
            sprintf(dimname, "%s_dim", group->time_index_name);
        } else {
            if (dim->rank > 0) {
                sprintf(dimname, "n%lld_dim", dim->rank);
            } else {
                dimname[0] = '\0';
            }
        }
    }

    return;
}

static int write_attribute(
        int ncid,
        int root_group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt,
        enum ADIOS_FLAG fortran_flag,
        int myrank,
        int nproc)
{
    int i, rank = 0;
    int rc;

    nc_type nc4_type_id;
    int varid;
    int attid;

    char attname[255];
    char varname[255];

    struct adios_dimension_struct * dims;
    struct adios_var_struct * var_linked;

    int err_code = 0;


    ncd_gen_name(attname, patt->path, patt->name);
    varid = -1;
    if (strcmp(patt->path,"/")==0) {
        varid = NC_GLOBAL;
        strcpy(attname, patt->name);
    } else {
        ncd_gen_name(varname, patt->path, "");
        rc = nc_inq_varid(ncid, varname, &varid);
        if (rc == NC_NOERR) {
//            if (myrank==0) fprintf (stderr, "NC4 ERROR variable(%s) exists in write_attribute, rc=%d\n", varname, rc);
        } else if (rc == NC_ENOTVAR) {
            if (myrank==0) fprintf (stderr, "NC4 ERROR variable(%s) does not exist in write_attribute, rc=%d\n", varname, rc);
            err_code = -2;
            goto escape;
        } else {
            if (myrank==0) fprintf (stderr, "NC4 ERROR inquiring about variable(%s) in write_attribute, rc=%d\n", varname, rc);
            err_code = -2;
            goto escape;
        }
    }

//    printf("looking for var(%s)\n", last);
//    printf("got varid(%d) for grp_id(%d).variable(%s) in write_attribute\n", varid, ncid, last);
    rc = nc_inq_attid(ncid, varid, attname, &attid);
    if (rc == NC_NOERR) {
        if (myrank==0) fprintf (stderr, "NC4 ERROR attribute(%s) already exists in write_attribute, rc=%d\n", attname, rc);
        err_code = 0;
        goto escape;
    } else if (rc == NC_ENOTATT) {
//        if (myrank==0) fprintf (stderr, "NC4 ERROR attribute(%s) does not exist in write_attribute, rc=%d\n", attname, rc);
    } else {
        if (myrank==0) fprintf (stderr, "NC4 ERROR inquiring about attribute(%s) in write_attribute, rc=%d\n", attname, rc);
        err_code = -2;
        goto escape;
    }

//    printf("patt->type=%d attname : %s\n", patt->type, attname);
    if (patt->type == -1) {
        var_linked = patt->var;
        if (!var_linked || (var_linked && !var_linked->data)) {
            fprintf (stderr, "NC4 ERROR: invalid data in var_linked(%s(%d)) (in attribute write), rc=%d\n"
                    ,var_linked->name, var_linked->id, rc);
            err_code = -2;
            goto escape;
        } else {
            dims = var_linked->dimensions;
        }
        getNC4TypeId(var_linked->type, &nc4_type_id, fortran_flag);
        // Scalar variable as attribute
        if (!dims) {
            rc = nc_put_att(ncid, varid, attname, nc4_type_id, 1, var_linked->data);
            if (rc != NC_NOERR) {
                fprintf (stderr, "NC4 ERROR unable to put attribute(%s) in write_attribute, rc=%d\n", attname, rc);
                err_code = -2;
                goto escape;
            }
        } else {
            fprintf (stderr, "NC4 ERROR multi-dimensional attribute(%s) unsupported for variable(%s) in write_attribute, rc=%d\n", attname, pvar_root->name, rc);
            err_code = -2;
            goto escape;
            //             while (dims) {
            //                ++rank;
            //                dims = dims->next;
            //             }
            //
            //             h5_localdims = (hsize_t *) malloc (rank * sizeof(hsize_t));
            //             dims = var_linked->dimensions;
            //             for ( i = 0; i < rank; i++) {
            //                 if ( dims->dimension.rank == 0 && dims->dimension.id) {
            //                     var_linked = adios_find_var_by_id (pvar_root , dims->dimension.id);
            //                     if ( var_linked) {
            //                         h5_localdims [i] = *(int *)var_linked->data;
            //                     }
            //                 } else {
            //                     h5_localdims [i] = dims->dimension.rank;
            //                 }
            //             }
            //             h5_dataspace_id = H5Screate_simple(rank,h5_localdims, NULL);
            //             h5_attribute_id = H5Aopen_name ( ncid, attname);
            //             if (h5_attribute_id < 0) {
            //                 h5_attribute_id = H5Acreate ( ncid, attname
            //                                          ,h5_type_id,h5_dataspace_id,0);
            //                 if (h5_attribute_id < 0) {
            //                     fprintf (stderr, "NC4 ERROR: getting negative attribute_id "
            //                                      "in write_attribute: %s\n", attname);
            //                     err_code = -2;
            //                  }
            //             }
            //             if (h5_attribute_id > 0) {
            //                 if (myrank == 0 && var_linked->data) {
            //                     H5Awrite ( h5_attribute_id, h5_type_id, var_linked->data);
            //                 }
            //                 H5Aclose ( h5_attribute_id);
            //             }
            //             H5Sclose ( h5_dataspace_id);
            //             free (h5_localdims);
        }
    }
    if (patt->type > 0) {
        getNC4TypeId(patt->type, &nc4_type_id, fortran_flag);
        if (nc4_type_id > 0) {
            if (patt->type == adios_string) {
                rc = nc_put_att_text(ncid, varid, attname, strlen((char *)patt->value), (const char *)patt->value);
                if (rc != NC_NOERR) {
                    fprintf (stderr, "NC4 ERROR unable to put attribute(%s) in write_attribute, rc=%d\n", attname, rc);
                    err_code = -2;
                    goto escape;
                }
            } else {
                rc = nc_put_att(ncid, varid, attname, nc4_type_id, 1, patt->value);
                if (rc != NC_NOERR) {
                    fprintf (stderr, "NC4 ERROR unable to put attribute(%s) in write_attribute, rc=%d\n", attname, rc);
                    err_code = -2;
                    goto escape;
                }
            }
        }
    }
escape:
    return err_code;
}

static int decipher_dims(
        int ncid,
        int root_group,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_var_struct *pvar,
        int myrank,
        int nproc,
        deciphered_dims_t *deciphered_dims)
{
    int i=0;

    struct adios_dimension_struct *dims;
    nc4_dimname_t *nc4_global_dimnames=NULL;
    nc4_dimname_t *nc4_local_dimnames=NULL;
    nc4_dimname_t *nc4_local_offset_names=NULL;
    size_t    *nc4_gbdims;
    size_t    *nc4_globaldims;
    size_t    *nc4_localdims;
    size_t    *nc4_offsets;
    ptrdiff_t *nc4_strides;
    int       *nc4_global_dimids;
    int       *nc4_local_dimids;
    int       *nc4_loffs_dimids;

    int nc4_gbglobaldims_dimids[2];
    int nc4_gbglobaldims_varid;

    enum ADIOS_FLAG has_globaldims=adios_flag_no;
    enum ADIOS_FLAG has_localdims=adios_flag_no;
    enum ADIOS_FLAG has_localoffsets=adios_flag_no;

    enum ADIOS_FLAG has_timedim=adios_flag_no;
    int             timedim_index=-1;

    int global_dim_count=0;
    int local_dim_count=0;
    int local_offset_count=0;

    char dimname[255];
    void * id;

    memset(deciphered_dims, 0, sizeof(deciphered_dims_t));

    dims=pvar->dimensions;
    while (dims) {
        if ((dims->dimension.is_time_index == adios_flag_yes) &&
            (dims->dimension.var == NULL && dims->dimension.attr == NULL)) {
            has_timedim = adios_flag_yes;
            timedim_index = local_dim_count;
            local_dim_count++;
        } else if ((dims->dimension.rank != 0) ||
            (dims->dimension.rank == 0) && 
            (dims->dimension.var != NULL || dims->dimension.attr != NULL)) {
            has_localdims=adios_flag_yes;
            local_dim_count++;
        }
        if ((dims->global_dimension.rank != 0) ||
            (dims->global_dimension.rank == 0) && 
            (dims->global_dimension.var != NULL || dims->global_dimension.attr != NULL)) {
            has_globaldims=adios_flag_yes;
            global_dim_count++;
        }
        if ((dims->local_offset.rank != 0) ||
            (dims->local_offset.rank == 0) && 
            (dims->local_offset.var != NULL || dims->local_offset.attr != NULL)) {
            has_localoffsets=adios_flag_yes;
            local_offset_count++;
        }
        if (DEBUG>3) printf("gdims[%d].rank=%llu; var=%d, attr=%d, is_time_index=%d\n", i, dims->global_dimension.rank, dims->global_dimension.var, dims->global_dimension.attr, dims->global_dimension.is_time_index);
        if (DEBUG>3) printf("ldims[%d].rank=%llu; var=%d, attr=%d, is_time_index=%d\n", i, dims->dimension.rank, dims->dimension.var, dims->dimension.attr, dims->dimension.is_time_index);
        if (DEBUG>3) printf("loffs[%d].rank=%llu; var=%d, attr=%d, is_time_index=%d\n", i, dims->local_offset.rank, dims->local_offset.var, dims->local_offset.attr, dims->local_offset.is_time_index);
        i++;
        dims = dims->next;
    }

    if (DEBUG>3) printf("global_dim_count  ==%d\n", global_dim_count);
    if (DEBUG>3) printf("local_dim_count   ==%d\n", local_dim_count);
    if (DEBUG>3) printf("calculated local_offset_count==%d\n", local_offset_count);
    if ((has_localoffsets == adios_flag_yes) && (local_offset_count < global_dim_count)) {
        if (DEBUG>3) printf("assuming local_offset_count should equal global_dim_count.  FORCING EQUALITY\n");
        local_offset_count = global_dim_count;
    }

    nc4_gbdims     = (size_t *)calloc(local_dim_count * 3, sizeof(size_t));
    nc4_globaldims = nc4_gbdims;
    nc4_localdims  = nc4_gbdims + local_dim_count;
    nc4_offsets    = nc4_gbdims + (2*local_dim_count);
    nc4_strides    = (ptrdiff_t *)calloc(local_dim_count, sizeof(ptrdiff_t));
    nc4_global_dimnames      = (nc4_dimname_t *)calloc(local_dim_count, sizeof(nc4_dimname_t));
    nc4_local_dimnames       = (nc4_dimname_t *)calloc(local_dim_count, sizeof(nc4_dimname_t));
    nc4_local_offset_names   = (nc4_dimname_t *)calloc(local_dim_count, sizeof(nc4_dimname_t));
    nc4_global_dimids     = (int *)calloc(global_dim_count, sizeof(int));
    nc4_local_dimids      = (int *)calloc(local_dim_count, sizeof(int));
    nc4_loffs_dimids      = (int *)calloc(local_offset_count, sizeof(int));

    dims = pvar->dimensions;
    for (i=0;i<global_dim_count;i++) {
        parse_dimension_name(group, pvar_root, patt_root, &dims->global_dimension, dimname);
        //ncd_gen_name(nc4_global_dimnames[i], pvar->path, dimname);
        ncd_gen_name(nc4_global_dimnames[i], "", dimname);
        if (DEBUG>3) printf("global_dimension[%d]->name==%s, ->rank==%llu, ->var==%d, ->attr==%d, is_time_index==%d\n",
                i, nc4_global_dimnames[i], dims->global_dimension.rank, dims->global_dimension.var, dims->global_dimension.attr, dims->global_dimension.is_time_index);
        if (dims) {
            dims = dims -> next;
        }
    }
    dims = pvar->dimensions;
    for (i=0;i<local_dim_count;i++) {
        parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dimname);
        //ncd_gen_name(nc4_local_dimnames[i], pvar->path, dimname);
        ncd_gen_name(nc4_local_dimnames[i], "", dimname);
        if (DEBUG>3) printf("local_dimension[%d]->name ==%s, ->rank==%llu, ->var==%d, ->attr==%d, is_time_index==%d\n",
                i, nc4_local_dimnames[i], dims->dimension.rank, dims->dimension.var, dims->dimension.attr, dims->dimension.is_time_index);
        if (dims) {
            dims = dims -> next;
        }
    }
    dims = pvar->dimensions;
    for (i=0;i<local_offset_count;i++) {
        parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, dimname);
        //ncd_gen_name(nc4_local_offset_names[i], pvar->path, dimname);
        ncd_gen_name(nc4_local_offset_names[i], "", dimname);
        if (DEBUG>3) printf("local_offset[%d]->name    ==%s, ->rank==%llu, ->var==%d, ->attr==%d, is_time_index==%d\n",
                i, nc4_local_offset_names[i], dims->local_offset.rank, dims->local_offset.var, dims->local_offset.attr, dims->local_offset.is_time_index);
        if (dims) {
            dims = dims -> next;
        }
    }
    int global_idx=0;
    int local_idx =0;
    int loffs_idx =0;
    dims = pvar->dimensions;
    if (myrank==0) if (DEBUG>3) printf("timedim_index=%d\n", timedim_index);
    while (dims) {
        /* get the global/local/offset arrays */
        nc4_strides[local_idx] = 1;
        if (timedim_index == local_idx) {
            nc4_globaldims[global_idx] = NC_UNLIMITED;
            nc4_localdims[local_idx] = 1;
            nc4_offsets[loffs_idx] = 0;
            parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dimname);
            //ncd_gen_name(nc4_global_dimnames[global_idx], pvar->path, dimname);
            ncd_gen_name(nc4_global_dimnames[global_idx], "", dimname);
            strcpy(nc4_local_dimnames[local_idx], nc4_global_dimnames[global_idx]);
            strcpy(nc4_local_offset_names[loffs_idx], nc4_global_dimnames[global_idx]);
            if ((global_dim_count < local_dim_count) && (local_idx < local_dim_count)) {
                global_idx++;
                loffs_idx++;
            }
        } else {
            parse_dimension_name(group, pvar_root, patt_root, &dims->dimension, dimname);
            if (dimname[0] == '\0') {
                sprintf(dimname, "local_%d", local_idx);
            }
            //ncd_gen_name(nc4_local_dimnames[local_idx], pvar->path, dimname);
            ncd_gen_name(nc4_local_dimnames[local_idx], "", dimname);
            parse_dimension_size(group, pvar_root, patt_root, &dims->dimension, &nc4_localdims[local_idx]);
        }
        if (myrank==0) {
            if (DEBUG>3) printf("\t%s[%d]: l(%d)", pvar->name, local_idx, nc4_localdims[local_idx]);
        }
        local_idx++;
        if (global_idx < local_dim_count) {
            parse_dimension_name(group, pvar_root, patt_root, &dims->global_dimension, dimname);
            if (dimname[0] == '\0') {
                sprintf(dimname, "global_%d", global_idx);
            }
            //ncd_gen_name(nc4_global_dimnames[global_idx], pvar->path, dimname);
            ncd_gen_name(nc4_global_dimnames[global_idx], "", dimname);
            parse_dimension_size(group, pvar_root, patt_root, &dims->global_dimension, &nc4_globaldims[global_idx]);
            if (myrank==0) {
                if (DEBUG>3) printf(":g(%d)", nc4_globaldims[global_idx]);
            }
            global_idx++;
        }
        if (loffs_idx < local_dim_count) {
            parse_dimension_name(group, pvar_root, patt_root, &dims->local_offset, dimname);
            if (dimname[0] == '\0') {
                sprintf(dimname, "offset_%d", loffs_idx);
            }
            //ncd_gen_name(nc4_local_offset_names[loffs_idx], pvar->path, dimname);
            ncd_gen_name(nc4_local_offset_names[loffs_idx], "", dimname);
            parse_dimension_size(group, pvar_root, patt_root, &dims->local_offset, &nc4_offsets[loffs_idx]);
            if (myrank==0) {
                if (DEBUG>3) printf(":o(%d)", nc4_offsets[loffs_idx]);
            }
            loffs_idx++;
        }
        if (myrank==0) {
            if (DEBUG>3) printf("\n");
        }

        if (dims) {
            dims = dims -> next;
        }
    }
    for (i=0;i<local_dim_count;i++) {
        if (DEBUG>3) printf("nc4_global_dimnames[%d]   ==%s, nc4_globaldims[%d]==%d\n", i, nc4_global_dimnames[i], i, nc4_globaldims[i]);
        if (DEBUG>3) printf("nc4_local_dimnames[%d]    ==%s, nc4_localdims[%d] ==%d\n", i, nc4_local_dimnames[i], i, nc4_localdims[i]);
        if (DEBUG>3) printf("nc4_local_offset_names[%d]==%s, nc4_offsets[%d]   ==%d\n", i, nc4_local_offset_names[i], i, nc4_offsets[i]);
    }

    if ((has_timedim == adios_flag_yes) && (global_dim_count < local_dim_count)) {
        global_dim_count++;
        local_offset_count++;
    }

    deciphered_dims->nc4_gbstrides[0]    = 1;
    deciphered_dims->nc4_gbstrides[1]    = 1;
    deciphered_dims->nc4_gbglobaldims[0] = nproc;
    deciphered_dims->nc4_gbglobaldims[1] = local_dim_count * 3;
    deciphered_dims->nc4_gboffsets[0]    = myrank;
    deciphered_dims->nc4_gboffsets[1]    = 0;
    deciphered_dims->nc4_gblocaldims[0]  = 1;
    deciphered_dims->nc4_gblocaldims[1]  = local_dim_count * 3;

    ncd_gen_name(dimname, pvar->path, pvar->name);
    sprintf(deciphered_dims->gbdims_name, "_%s_gbdims", dimname);
    sprintf(deciphered_dims->gbdims_dim0_name, "_%s_gbdims_dim0", dimname);
    sprintf(deciphered_dims->gbdims_dim1_name, "_%s_gbdims_dim1", dimname);


    /*
     * Copy local scalers and pointers into deciphered_dims
     */
    deciphered_dims->dims=pvar->dimensions;
    deciphered_dims->nc4_global_dimnames=nc4_global_dimnames;
    deciphered_dims->nc4_local_dimnames=nc4_local_dimnames;
    deciphered_dims->nc4_local_offset_names=nc4_local_offset_names;

    deciphered_dims->nc4_gbdims    =nc4_gbdims;
    deciphered_dims->nc4_globaldims=nc4_globaldims;
    deciphered_dims->nc4_localdims =nc4_localdims;
    deciphered_dims->nc4_offsets   =nc4_offsets;
    deciphered_dims->nc4_strides   =nc4_strides;
    deciphered_dims->nc4_global_dimids=nc4_global_dimids;
    deciphered_dims->nc4_local_dimids =nc4_local_dimids;
    deciphered_dims->nc4_loffs_dimids =nc4_loffs_dimids;

    memcpy(deciphered_dims->nc4_gbglobaldims_dimids, nc4_gbglobaldims_dimids, 2*sizeof(int));
    deciphered_dims->nc4_gbglobaldims_varid=nc4_gbglobaldims_varid;

    deciphered_dims->has_globaldims=has_globaldims;
    deciphered_dims->has_localdims=has_localdims;
    deciphered_dims->has_localoffsets=has_localoffsets;

    deciphered_dims->has_timedim   =has_timedim;
    deciphered_dims->timedim_index =timedim_index;

    deciphered_dims->global_dim_count=global_dim_count;
    deciphered_dims->local_dim_count=local_dim_count;
    deciphered_dims->local_offset_count=local_offset_count;
}

static int cleanup_deciphered_dims(
        deciphered_dims_t *deciphered_dims)
{
    if (deciphered_dims->nc4_gbdims             != NULL) free(deciphered_dims->nc4_gbdims);
    if (deciphered_dims->nc4_global_dimnames    != NULL) free(deciphered_dims->nc4_global_dimnames);
    if (deciphered_dims->nc4_local_dimnames     != NULL) free(deciphered_dims->nc4_local_dimnames);
    if (deciphered_dims->nc4_local_offset_names != NULL) free(deciphered_dims->nc4_local_offset_names);
    if (deciphered_dims->nc4_strides            != NULL) free(deciphered_dims->nc4_strides);
    if (deciphered_dims->nc4_global_dimids      != NULL) free(deciphered_dims->nc4_global_dimids);
    if (deciphered_dims->nc4_local_dimids       != NULL) free(deciphered_dims->nc4_local_dimids);
    if (deciphered_dims->nc4_loffs_dimids       != NULL) free(deciphered_dims->nc4_loffs_dimids);
}

static int read_var(
        int ncid,
        int root_group,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_var_struct *pvar,
        enum ADIOS_FLAG fortran_flag,
        int myrank,
        int nproc)
{
    int return_code=0;
    int i, rc;
    struct adios_dimension_struct * dims = pvar->dimensions;
    deciphered_dims_t deciphered_dims;
    char fullname[255];

    nc_type nc4_type_id;
    int nc4_varid;

    memset(&deciphered_dims, 0, sizeof(deciphered_dims_t));

    getNC4TypeId (pvar->type, &nc4_type_id, fortran_flag);
    if (nc4_type_id <=0 )
    {
        fprintf (stderr, "ERROR in getNC4TypeId in read_var!\n");
        return_code=-2;
        goto escape;
    }

    ncd_gen_name(fullname, pvar->path, pvar->name);

    rc = nc_inq_varid(ncid, fullname, &nc4_varid);
    if (rc == NC_ENOTVAR) {
        fprintf(stderr, "NC4 ERROR variable(%s) does not exist in read_var, rc=%d\n", fullname, rc);
        return_code=-2;
        goto escape;
    } else if (rc != NC_NOERR) {
        fprintf(stderr, "NC4 ERROR checking existence of variable(%s) in read_var, rc=%d\n", fullname, rc);
        return_code=-2;
        goto escape;
    }

    if(myrank==0) if (DEBUG>3) printf("\tenter global reading!\n");

    if (!pvar->dimensions) { // begin scalar read
        rc = nc_get_var(ncid, nc4_varid, pvar->data);
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR getting scalar variable(%s) in read_var\n", fullname);
            return_code=-2;
            goto escape;
        }

        return_code=0;
        goto escape;
    } // end scalar write

    if (myrank==0) if (DEBUG>3) printf("read_var deciphering dims\n");
    decipher_dims(ncid,
            root_group,
            group,
            pvar_root,
            patt_root,
            pvar,
            myrank,
            nproc,
            &deciphered_dims);
    dims = pvar->dimensions;


    if (deciphered_dims.has_timedim == adios_flag_no) {

        /* begin reading array with fixed dimensions */

        if (myrank==0) if (DEBUG>3) printf("\treading fixed dimension array var!\n");

        for (i=0;i<deciphered_dims.local_dim_count;i++) {
            if(myrank==0) {
                if (DEBUG>3) printf("\tDIMS var:%s dim[%d]:  %d %d %d\n",fullname
                        ,i, deciphered_dims.nc4_globaldims[i], deciphered_dims.nc4_localdims[i], deciphered_dims.nc4_offsets[i]);
          }
        }

        Func_Timer("nc4_varid par_access", rc = nc_var_par_access(ncid, nc4_varid, NC_COLLECTIVE););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR setting parallel access for scalar variable(%s) in read_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }

//        rc = nc_get_vars(ncid, nc4_varid, deciphered_dims.nc4_offsets, deciphered_dims.nc4_localdims, deciphered_dims.nc4_strides, pvar->data);
        Func_Timer("getvars", rc = nc_get_vara(ncid, nc4_varid, deciphered_dims.nc4_offsets, deciphered_dims.nc4_localdims, pvar->data););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR getting array variable(%s) in read_var\n", fullname);
            return_code=-2;
            goto escape;
        }

        /* end reading array with fixed dimensions */

    } else {

        /* begin reading array with unlimited dimension */

        size_t current_timestep=0;

        for (i=0;i<deciphered_dims.local_dim_count;i++) {
            if(myrank==0) {
                if (DEBUG>3) printf("\tDIMS var:%s dim[%d]:  %d %d %d\n",fullname
                        ,i, deciphered_dims.nc4_globaldims[i], deciphered_dims.nc4_localdims[i], deciphered_dims.nc4_offsets[i]);
          }
        }

        Func_Timer("nc4_varid par_access", rc = nc_var_par_access(ncid, nc4_varid, NC_COLLECTIVE););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR setting parallel access for scalar variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }

        Func_Timer("inqdim", rc = nc_inq_dimid(ncid, deciphered_dims.nc4_local_dimnames[deciphered_dims.timedim_index], &deciphered_dims.nc4_global_dimids[deciphered_dims.timedim_index]););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR inquiring about dimension(%s) for array variable(%s) in write_var, rc=%d\n", deciphered_dims.nc4_local_dimnames[i], fullname, rc);
            return_code=-2;
            goto escape;
        }
        /* get the current timestep */
        Func_Timer("inqdimlen", rc = nc_inq_dimlen(ncid, deciphered_dims.nc4_global_dimids[deciphered_dims.timedim_index], &current_timestep););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR error getting current timestep for array variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }
        if (DEBUG>3) printf("current_timestep==%d\n", current_timestep);

        /* decrement.  dims are 1-based, while offsets are 0-based. */
        deciphered_dims.nc4_offsets[deciphered_dims.timedim_index]=current_timestep-1;
        for (i=0;i<deciphered_dims.local_dim_count;i++) {
            if (DEBUG>3) printf("write_var: deciphered_dims.nc4_offsets[%d]=%lu deciphered_dims.nc4_localdims[%d]=%lu\n",
                    i, deciphered_dims.nc4_offsets[i],
                    i, deciphered_dims.nc4_localdims[i]);
        }
//        rc = nc_get_vars(ncid, nc4_varid, deciphered_dims.nc4_offsets, deciphered_dims.nc4_localdims, deciphered_dims.nc4_strides, pvar->data);
        Func_Timer("getvars", rc = nc_get_vara(ncid, nc4_varid, deciphered_dims.nc4_offsets, deciphered_dims.nc4_localdims, pvar->data););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR getting array variable(%s) in read_var\n", fullname);
            return_code=-2;
            goto escape;
        }

        /* end reading array with unlimited dimension */

    }

escape:
    cleanup_deciphered_dims(&deciphered_dims);

    return return_code;
}

static int write_header(
        int ncid,
        int root_group,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_var_struct *pvar,
        enum ADIOS_FLAG fortran_flag,
        int myrank,
        int nproc)
{
    int i;
    int rc;
    int return_code=0;
    nc_type nc4_type_id;
    int nc4_varid;
    deciphered_dims_t deciphered_dims;
    char fullname[255];

//    int myrank=md->rank;
//    int nproc=md->size;

//    struct adios_var_struct *pvar=fd->group->vars;
//    enum ADIOS_FLAG fortran_flag=fd->group->adios_host_language_fortran;

    memset(&deciphered_dims, 0, sizeof(deciphered_dims_t));

    getNC4TypeId(pvar->type, &nc4_type_id, fortran_flag);
    if(nc4_type_id <= 0) {
        fprintf(stderr, "NC4 ERROR in getH5TypeId in write_header\n");
        return_code=-2;
        goto escape;
    }

    ncd_gen_name(fullname, pvar->path, pvar->name);

    if (!pvar->dimensions) { // begin scalar write
        Func_Timer("inqvar", rc = nc_inq_varid(ncid, fullname, &nc4_varid););
        if (rc == NC_ENOTVAR) {
            if (pvar->type == adios_string) {
                int str_var_dimid=0;
                char str_var_dimname[40];
                sprintf(str_var_dimname, "%s_dim", fullname);
                Func_Timer("defdim", rc = nc_def_dim(ncid, str_var_dimname, strlen((char *)pvar->data)+1, &str_var_dimid););
                if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR defining string variable(%s) dim in write_header, rc=%d\n", fullname, rc);
                    return_code=-2;
                    goto escape;
                }
                Func_Timer("defvar", rc = nc_def_var(ncid, fullname, nc4_type_id, 1, &str_var_dimid, &nc4_varid););
                if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR defining string variable(%s) in write_header, rc=%d\n", fullname, rc);
                    return_code=-2;
                    goto escape;
                }
            } else {
                Func_Timer("defvar", rc = nc_def_var(ncid, fullname, nc4_type_id, 0, NULL, &nc4_varid););
                if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR defining scalar variable(%s) in write_header, rc=%d\n", fullname, rc);
                    return_code=-2;
                    goto escape;
                }
            }
        }

        goto escape;
    } // end scalar write

    if (myrank==0) if (DEBUG>3) printf("write_header deciphering dims\n");
    decipher_dims(ncid,
            root_group,
            group,
            pvar_root,
            patt_root,
            pvar,
            myrank,
            nproc,
            &deciphered_dims);

    enum ADIOS_FLAG var_exists = adios_flag_yes;
    Func_Timer("inqvar", rc = nc_inq_varid(ncid, fullname, &nc4_varid););
    if (rc == NC_ENOTVAR) {
        var_exists=adios_flag_no;
    } else if (rc != NC_NOERR) {
        fprintf(stderr, "NC4 ERROR checking existence of variable(%s) in write_header, rc=%d\n", fullname, rc);
        return_code=-2;
        goto escape;
    }

    if (deciphered_dims.has_timedim == adios_flag_no) {

        /* begin writing array with fixed dimensions */

        if (myrank==0) if (DEBUG>3) printf("\twriting fixed dimension array var!\n");

        if (var_exists == adios_flag_no) {

#ifdef THIS_IS_UNDEFINED /* I do not think we need to define ADIOS local dimension variables in NC4 */
            for (i=0;i<deciphered_dims.local_dim_count;i++) {
                Func_Timer("inqdim", rc = nc_inq_dimid(ncid, deciphered_dims.nc4_local_dimnames[i], &deciphered_dims.nc4_local_dimids[i]););
                if (rc == NC_EBADDIM) {
                    Func_Timer("defdim", rc = nc_def_dim(ncid, deciphered_dims.nc4_local_dimnames[i], deciphered_dims.nc4_localdims[i], &deciphered_dims.nc4_local_dimids[i]););
                    if (rc != NC_NOERR) {
                        fprintf(stderr, "NC4 ERROR defining array dimension(%s) in write_header, rc=%d\n", deciphered_dims.nc4_local_dimnames[i], rc);
                        return_code=-2;
                        goto escape;
                    }
                } else if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR inquiring about dimension(%s) for array variable(%s) in write_header, rc=%d\n", deciphered_dims.nc4_local_dimnames[i], fullname, rc);
                    return_code=-2;
                    goto escape;
                }
            }
#endif

            for (i=0;i<deciphered_dims.global_dim_count;i++) {
                Func_Timer("inqdim", rc = nc_inq_dimid(ncid, deciphered_dims.nc4_global_dimnames[i], &deciphered_dims.nc4_global_dimids[i]););
                if (rc == NC_EBADDIM) {
                    Func_Timer("defdim", rc = nc_def_dim(ncid, deciphered_dims.nc4_global_dimnames[i], deciphered_dims.nc4_globaldims[i], &deciphered_dims.nc4_global_dimids[i]););
                    if (rc != NC_NOERR) {
                        fprintf(stderr, "NC4 ERROR defining array dimension(%s) in write_header, rc=%d\n", deciphered_dims.nc4_global_dimnames[i], rc);
                        return_code=-2;
                        goto escape;
                    }
                } else if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR inquiring about dimension(%s) for array variable(%s) in write_header, rc=%d\n", deciphered_dims.nc4_global_dimnames[i], fullname, rc);
                    return_code=-2;
                    goto escape;
                }
            }

            if (deciphered_dims.has_globaldims == adios_flag_yes) {
                Func_Timer("defvar", rc = nc_def_var(ncid, fullname, nc4_type_id, deciphered_dims.global_dim_count, deciphered_dims.nc4_global_dimids, &nc4_varid););
                if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR defining array variable(%s) with global dims in write_header, rc=%d\n", fullname, rc);
                    return_code=-2;
                    goto escape;
                }
            } else {
                Func_Timer("defvar", rc = nc_def_var(ncid, fullname, nc4_type_id, deciphered_dims.local_dim_count, deciphered_dims.nc4_local_dimids, &nc4_varid););
                if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR defining array variable(%s) with local dims in write_header, rc=%d\n", fullname, rc);
                    return_code=-2;
                    goto escape;
                }
            }
        }

        if (DEBUG>3) printf("got varid(%d) for grp_id(%d).variable(%s) in write_header, rc=%d\n", nc4_varid, ncid, fullname, rc);
        if (DEBUG>3) printf("sizeof(size_t)==%d\n", sizeof(size_t));

        /* end writing array with fixed dimensions */

    } else {

        /* begin writing array with unlimited dimension */

        size_t current_timestep=0;

        if (myrank==0) if (DEBUG>3) printf("\twriting timestep array var!\n");

        if (var_exists == adios_flag_no) {
            /* define the dims and var */
            for (i=0;i<deciphered_dims.global_dim_count;i++) {
                if (DEBUG>3) printf("inq dim name=%s, size=%d\n", deciphered_dims.nc4_global_dimnames[i], deciphered_dims.nc4_globaldims[i]);
                Func_Timer("inqdim", rc = nc_inq_dimid(ncid, deciphered_dims.nc4_global_dimnames[i], &deciphered_dims.nc4_global_dimids[i]););
                if (rc == NC_EBADDIM) {
                    if (DEBUG>3) printf("def dim name=%s, size=%d\n", deciphered_dims.nc4_global_dimnames[i], deciphered_dims.nc4_globaldims[i]);
                    Func_Timer("defdim", rc = nc_def_dim(ncid, deciphered_dims.nc4_global_dimnames[i], deciphered_dims.nc4_globaldims[i], &deciphered_dims.nc4_global_dimids[i]););
                    if (rc != NC_NOERR) {
                        fprintf(stderr, "NC4 ERROR defining array dimension(%s) in write_header, rc=%d\n", deciphered_dims.nc4_global_dimnames[i], rc);
                        return_code=-2;
                        goto escape;
                    }
                } else if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR inquiring about dimension(%s) for array variable(%s) in write_header, rc=%d\n", deciphered_dims.nc4_global_dimnames[i], fullname, rc);
                    return_code=-2;
                    goto escape;
                }
            }
            for (i=0;i<deciphered_dims.local_dim_count;i++) {
                if (DEBUG>3) printf("inq dim name=%s, size=%d\n", deciphered_dims.nc4_local_dimnames[i], deciphered_dims.nc4_localdims[i]);
                Func_Timer("inqdim", rc = nc_inq_dimid(ncid, deciphered_dims.nc4_local_dimnames[i], &deciphered_dims.nc4_local_dimids[i]););
                if (rc == NC_EBADDIM) {
                    if (DEBUG>3) printf("def dim name=%s, size=%d\n", deciphered_dims.nc4_local_dimnames[i], deciphered_dims.nc4_localdims[i]);
                    Func_Timer("defdim", rc = nc_def_dim(ncid, deciphered_dims.nc4_local_dimnames[i], deciphered_dims.nc4_localdims[i], &deciphered_dims.nc4_local_dimids[i]););
                    if (rc != NC_NOERR) {
                        fprintf(stderr, "NC4 ERROR defining array dimension(%s) in write_header, rc=%d\n", deciphered_dims.nc4_global_dimnames[i], rc);
                        return_code=-2;
                        goto escape;
                    }
                } else if (rc != NC_NOERR) {
                    fprintf(stderr, "NC4 ERROR inquiring about dimension(%s) for array variable(%s) in write_header, rc=%d\n", deciphered_dims.nc4_global_dimnames[i], fullname, rc);
                    return_code=-2;
                    goto escape;
                }
            }

            Func_Timer("defvar", rc = nc_def_var(ncid, fullname, nc4_type_id, deciphered_dims.local_dim_count, deciphered_dims.nc4_global_dimids, &nc4_varid););
            if (rc != NC_NOERR) {
                fprintf(stderr, "NC4 ERROR defining array variable(%s) in write_header, rc=%d\n", fullname, rc);
                return_code=-2;
                goto escape;
            }
        }

        /* end writing array with unlimited dimension */

    }

escape:
    cleanup_deciphered_dims(&deciphered_dims);

    return return_code;
}

static int write_var(
        int ncid,
        int root_group,
        struct adios_group_struct *group,
        struct adios_var_struct *pvar_root,
        struct adios_attribute_struct *patt_root,
        struct adios_var_struct *pvar,
        enum ADIOS_FLAG fortran_flag,
        int myrank,
        int nproc)
{
    int i;
    int rc;
    int return_code=0;
    nc_type nc4_type_id;
    int nc4_varid;
    deciphered_dims_t deciphered_dims;
    char fullname[255];

    memset(&deciphered_dims, 0, sizeof(deciphered_dims_t));

    getNC4TypeId(pvar->type, &nc4_type_id, fortran_flag);
    if(nc4_type_id <= 0) {
        fprintf(stderr, "NC4 ERROR in getH5TypeId in write_var\n");
        return_code=-2;
        goto escape;
    }

    ncd_gen_name(fullname, pvar->path, pvar->name);

    Func_Timer("inqvar", rc = nc_inq_varid(ncid, fullname, &nc4_varid););
    if (rc == NC_ENOTVAR) {
        write_header(ncid, root_group, group, pvar_root, patt_root, pvar, fortran_flag, myrank, nproc);
//        return 0;
    }

//    Func_Timer("enddef", rc = nc_enddef(ncid););
//    if (rc != NC_NOERR) {
//        if (myrank==0) fprintf(stderr, "NC4 ERROR ending define mode for scalar variable(%s) in write_var, rc=%d\n", fullname, rc);
////        return_code=-2;
////        goto escape;
//    }

    if (DEBUG>3) printf("rank(%d) write_var: ncid(%lu) varid(%lu) pvar->data=%p\n", global_rank, ncid, nc4_varid, pvar->data);

    if (!pvar->dimensions) { // begin scalar write
        Func_Timer("inqvar", rc = nc_inq_varid(ncid, fullname, &nc4_varid););
        if (rc == NC_ENOTVAR) {
            fprintf(stderr, "NC4 ERROR scalar variable(%s) does not exist in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        } else if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR checking existence of variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }
        Func_Timer("nc4_varid par_access", rc = nc_var_par_access(ncid, nc4_varid, NC_COLLECTIVE););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR setting parallel access for scalar variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }

        Func_Timer("putvar", rc = nc_put_var(ncid, nc4_varid, pvar->data););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR putting scalar variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }
        if (DEBUG>3) printf("groupid=%d\n",ncid);
        if (DEBUG>3) printf("write dataset: name=%s/%s rc=%d myrank=%d\n"
                 , pvar->path,fullname, rc, myrank);

        goto escape;
    } // end scalar write


    if (myrank==0) if (DEBUG>3) printf("write_var deciphering dims\n");
    decipher_dims(ncid,
            root_group,
            group,
            pvar_root,
            patt_root,
            pvar,
            myrank,
            nproc,
            &deciphered_dims);

    enum ADIOS_FLAG var_exists = adios_flag_yes;
    Func_Timer("inqvar", rc = nc_inq_varid(ncid, fullname, &nc4_varid););
    if (rc == NC_ENOTVAR) {
        fprintf(stderr, "NC4 ERROR array variable(%s) does not exist in write_var, rc=%d\n", fullname, rc);
        return_code=-2;
        goto escape;
    } else if (rc != NC_NOERR) {
        fprintf(stderr, "NC4 ERROR checking existence of variable(%s) in write_var, rc=%d\n", fullname, rc);
        return_code=-2;
        goto escape;
    }

    if (deciphered_dims.has_timedim == adios_flag_no) {

        /* begin writing array with fixed dimensions */

        if (myrank==0) if (DEBUG>3) printf("\twriting fixed dimension array var!\n");

        Func_Timer("nc4_varid par_access", rc = nc_var_par_access(ncid, nc4_varid, NC_COLLECTIVE););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR setting parallel access for scalar variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }

        if (DEBUG>3) printf("got varid(%d) for grp_id(%d).variable(%s) in write_var, rc=%d\n", nc4_varid, ncid, fullname, rc);
        if (DEBUG>3) printf("sizeof(size_t)==%d\n", sizeof(size_t));

//        Func_Timer("putvars", rc = nc_put_vars(ncid, nc4_varid, deciphered_dims.nc4_offsets, deciphered_dims.nc4_localdims, deciphered_dims.nc4_strides, pvar->data););
        Func_Timer("putvars", rc = nc_put_vara(ncid, nc4_varid, deciphered_dims.nc4_offsets, deciphered_dims.nc4_localdims, pvar->data););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR putting to array variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }

        /* end writing array with fixed dimensions */

    } else {

        /* begin writing array with unlimited dimension */

        size_t current_timestep=0;

        if (myrank==0) if (DEBUG>3) printf("\twriting timestep array var!\n");

        Func_Timer("nc4_varid par_access", rc = nc_var_par_access(ncid, nc4_varid, NC_COLLECTIVE););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR setting parallel access for scalar variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }

        Func_Timer("inqdim", rc = nc_inq_dimid(ncid, deciphered_dims.nc4_local_dimnames[deciphered_dims.timedim_index], &deciphered_dims.nc4_global_dimids[deciphered_dims.timedim_index]););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR inquiring about dimension(%s) for array variable(%s) in write_var, rc=%d\n", deciphered_dims.nc4_local_dimnames[i], fullname, rc);
            return_code=-2;
            goto escape;
        }
        /* get the current timestep */
        Func_Timer("inqdimlen", rc = nc_inq_dimlen(ncid, deciphered_dims.nc4_global_dimids[deciphered_dims.timedim_index], &current_timestep););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR error getting current timestep for array variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }
        if (DEBUG>3) printf("current_timestep==%d\n", current_timestep);
        /* the next timestep goes after the current.  */
        /* THK: don't increment.  dims are 1-based, while offsets are 0-based. */

        deciphered_dims.nc4_offsets[deciphered_dims.timedim_index]=current_timestep;
        for (i=0;i<deciphered_dims.local_dim_count;i++) {
            if (DEBUG>3) printf("write_var: deciphered_dims.nc4_offsets[%d]=%lu deciphered_dims.nc4_localdims[%d]=%lu\n",
                    i, deciphered_dims.nc4_offsets[i],
                    i, deciphered_dims.nc4_localdims[i]);
        }
//        Func_Timer("putvars", rc = nc_put_vars(ncid, nc4_varid, deciphered_dims.nc4_offsets, deciphered_dims.nc4_localdims, deciphered_dims.nc4_strides, pvar->data););
        Func_Timer("putvars", rc = nc_put_vara(ncid, nc4_varid, deciphered_dims.nc4_offsets, deciphered_dims.nc4_localdims, pvar->data););
        if (rc != NC_NOERR) {
            fprintf(stderr, "NC4 ERROR putting to array variable(%s) in write_var, rc=%d\n", fullname, rc);
            return_code=-2;
            goto escape;
        }

        /* end writing array with unlimited dimension */

    }

escape:
    cleanup_deciphered_dims(&deciphered_dims);

    return return_code;
}



static int adios_nc4_initialized = 0;
void adios_nc4_init(
        const PairStruct *parameters,
        struct adios_method_struct *method)
{
    struct adios_nc4_data_struct *md=NULL;
//    struct adios_nc4_data_struct *md = (struct adios_nc4_data_struct *)method->method_data;

    if (!adios_nc4_initialized) {
        adios_nc4_initialized = 1;

        MPI_Comm_rank(method->init_comm, &global_rank);

        list_init(&open_file_list, open_file_free);
    }


//    method->method_data = malloc(sizeof(struct adios_nc4_data_struct));
//    md = (struct adios_nc4_data_struct *)method->method_data;
//    md->ncid       = -1;
//    md->root_ncid  = -1;
//    md->rank       = -1;
//    md->size       = 0;
}

enum BUFFERING_STRATEGY adios_nc4_should_buffer(
        struct adios_file_struct *fd,
        struct adios_method_struct *method)
{
    int rc=NC_NOERR;

    struct open_file *of=NULL;
    struct adios_nc4_data_struct *md=NULL;
//    struct adios_nc4_data_struct *md = (struct adios_nc4_data_struct *)method->method_data;
    char *name;
    MPI_Info info = MPI_INFO_NULL;

    if (DEBUG>3) printf("enter adios_nc4_should_buffer (%s)\n", fd->name);

    of=open_file_find(method->base_path, fd->name);
    if (of == NULL) {
        fprintf(stderr, "file(%s, %s) is not open.  FAIL.\n", method->base_path, fd->name);
        return adios_flag_no;
    }
    md=of->md;

    if (md->ncid != -1) {
        // file already open
        if (DEBUG>3) printf("adios_nc4_should_buffer: file is already open (fname=%s, ncid=%d)\n", fd->name, md->ncid);
        return adios_flag_no;
    }

    if (md->group_comm != MPI_COMM_NULL) {
        if (DEBUG>3) printf("global_rank(%d): adios_nc4_should_buffer: get rank and size: group_comm(%p)\n", global_rank, md->group_comm);
        MPI_Comm_rank(md->group_comm, &md->rank);
        MPI_Comm_size(md->group_comm, &md->size);
        if (DEBUG>3) printf("global_rank(%d): adios_nc4_should_buffer: size(%d) rank(%d)\n", global_rank, md->size, md->rank);
    } else {
        md->group_comm=MPI_COMM_SELF;
    }
    fd->group->process_id = md->rank;
    name = malloc(strlen(method->base_path) + strlen(fd->name) + 1);
    sprintf(name, "%s%s", method->base_path, fd->name);

    int myrank=md->rank;

    // create a new file. If file exists its contents will be overwritten. //

    MPI_Info_create(&info);
    MPI_Info_set(info,"cb_align","2");
    MPI_Info_set(info,"romio_ds_write","disable");
    MPI_Info_set(info,"romio_cb_write","enable");

    switch (fd->mode) {
        case adios_mode_read:
        {
            Func_Timer("nc_open_par", rc = nc_open_par(name, NC_NOWRITE|NC_MPIIO, md->group_comm, info, &md->ncid););
            if (rc != NC_NOERR) {
                fprintf (stderr, "ADIOS NC4: could not open file(%s) for reading, rc=%d\n", name, rc);
                free (name);
                return adios_flag_no;
            }
            break;
        }
        case adios_mode_write:
        case adios_mode_append:
        {
            Func_Timer("nc_create_par", rc = nc_create_par(name, NC_NOCLOBBER|NC_MPIIO|NC_NETCDF4, md->group_comm, info, &md->ncid););
            if (rc == NC_EEXIST) {
                Func_Timer("nc_open_par", rc = nc_open_par(name, NC_WRITE|NC_MPIIO, md->group_comm, info, &md->ncid););
                if (rc != NC_NOERR) {
                    fprintf (stderr, "ADIOS NC4: could not open file(%s) for writing, rc=%d\n", name, rc);
                    free (name);
                    return adios_flag_no;
                }
            } else if (rc != NC_NOERR) {
                fprintf (stderr, "ADIOS NC4: cannot create file(%s), rc=%d\n", name, rc);
                free (name);
                return adios_flag_no;
            }
            break;
        }
    }

    md->root_ncid = md->ncid;

    free(name);

    return no_buffering;
}

int adios_nc4_open(
        struct adios_file_struct *fd,
        struct adios_method_struct *method,
        MPI_Comm comm)
{
    struct open_file *of=NULL;
    struct adios_nc4_data_struct *md=NULL;
//    struct adios_nc4_data_struct * md = (struct adios_nc4_data_struct *)method->method_data;

    if (DEBUG>3) printf("enter adios_nc4_open (%s)\n", fd->name);

    of=open_file_find(method->base_path, fd->name);
    if (of == NULL) {
        md             = malloc(sizeof(struct adios_nc4_data_struct));
        md->fd         = -1;
        md->ncid       = -1;
        md->root_ncid  = -1;
        md->rank       = -1;
        md->size       = 0;
        md->group_comm = comm;

        of=open_file_create(method->base_path, fd->name, md, fd);
    } else {
        md=of->md;

        // sanity check
        if (md->fd == -1) {
            if (DEBUG>3) printf("open: %s is open but fd==-1.  sanity check failed.  attempting reopen.\n", fd->name);
            open_file_delete(of->fpath, of->fname);
        } else {
            // file already open
            return adios_flag_no;
        }
    }


    if (DEBUG>3) printf("open: fname=%s; fd==%p; ncid=%d\n", fd->name, md->fd, md->ncid);

    list_ins_next(&open_file_list, list_tail(&open_file_list), of);

    open_file_printall();


    return 1;
}

void adios_nc4_write(
        struct adios_file_struct *fd,
        struct adios_var_struct *v,
        void *data,
        struct adios_method_struct *method)
{
    struct open_file *of=NULL;
    struct adios_nc4_data_struct *md=NULL;
//    struct adios_nc4_data_struct * md = (struct adios_nc4_data_struct *)method->method_data;
    static int first_write = 1;

    of=open_file_find(method->base_path, fd->name);
    if (of == NULL) {
        fprintf(stderr, "file(%s, %s) is not open.  FAIL.\n", method->base_path, fd->name);
        return;
    }
    md=of->md;

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append) {
//        if (first_write == 1) {
//            write_header(fd, md);
//            first_write = 0;
//        }

        if (md->rank==0) {
            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
            if (DEBUG>3) fprintf(stderr, "write var: %s start!\n", v->name);
        }
        write_var(md->ncid,
                md->root_ncid,
                fd->group,
                fd->group->vars,
                fd->group->attributes,
                v,
                fd->group->adios_host_language_fortran,
                md->rank,md->size);
    } else {
        if (DEBUG>3) fprintf(stderr, "entering unknown nc4 mode %d!\n", fd->mode);
    }
    if (md->rank==0) {
        if (DEBUG>3) fprintf(stderr, "write var: %s end!\n", v->name);
        if (DEBUG>3) fprintf(stderr, "-------------------------\n");
    }
}


void adios_nc4_read(
        struct adios_file_struct *fd,
        struct adios_var_struct *v,
        void *buffer,
        uint64_t buffersize,
        struct adios_method_struct *method)
{
    struct open_file *of=NULL;
    struct adios_nc4_data_struct *md=NULL;
//    struct adios_nc4_data_struct * md = (struct adios_nc4_data_struct *)method->method_data;

    of=open_file_find(method->base_path, fd->name);
    if (of == NULL) {
        fprintf(stderr, "file(%s, %s) is not open.  FAIL.\n", method->base_path, fd->name);
        return;
    }
    md=of->md;

    if(fd->mode == adios_mode_read) {
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

        if (md->rank==0) {
            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
            if (DEBUG>3) fprintf(stderr, "read var: %s! start\n", v->name);
        }
        read_var(md->ncid,
                md->root_ncid,
                fd->group,
                fd->group->vars,
                fd->group->attributes,
                v,
                fd->group->adios_host_language_fortran,
                md->rank,
                md->size);
        if (md->rank==0) {
            if (DEBUG>3) fprintf(stderr, "read var: %s! end\n", v->name);
            if (DEBUG>3) fprintf(stderr, "-------------------------\n");
        }

    }
}

void adios_nc4_buffer_overflow (struct adios_file_struct * fd, 
                                struct adios_method_struct * method)
{
    // this never happens without shared buffering
}

void adios_nc4_close(
        struct adios_file_struct *fd,
        struct adios_method_struct *method)
{
    struct open_file *of=NULL;
    struct adios_nc4_data_struct *md=NULL;
//    struct adios_nc4_data_struct * md = (struct adios_nc4_data_struct*)method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;
    int myrank;

    of=open_file_find(method->base_path, fd->name);
    if (of == NULL) {
        fprintf(stderr, "file(%s, %s) is not open.  FAIL.\n", method->base_path, fd->name);
        return;
    }
    md=of->md;
    myrank=md->rank;

    if (fd->mode == adios_mode_read) {
        struct adios_var_struct * v = fd->group->vars;
        while (v)
        {
            v->data = v->adata = 0;
            v = v->next;
        }

        if (md->rank==0) {
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
            if (DEBUG>1) fprintf(stderr, "reading done, nc4 file is virtually closed;\n");
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
        }
    } else if (fd->mode == adios_mode_write || fd->mode == adios_mode_append) {
        if (DEBUG>3) fprintf(stderr, "entering nc4 write attribute mode!\n");
        // FIXME: temporarily removed attributes writing and right now,
        // we don't support writing attrs in PHDF5/NC4 methods. 
     
/*
        while(a) {
            if (strcmp(a->path, "/__adios__")) {
                Func_Timer("write_attribute", 
                        write_attribute(md->ncid, md->root_ncid, fd->group->vars, a,
                            fd->group->adios_host_language_fortran,
                            md->rank,
                            md->size););
            }
            a = a->next;
        }
*/
        if (md->rank==0) {
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
            if (DEBUG>1) fprintf(stderr, "writing done, nc4 file is virtually closed;\n");
            if (DEBUG>1) fprintf(stderr, "-------------------------\n");
        }
    }

    Func_Timer("nc_sync", nc_sync(md->ncid););
    Func_Timer("nc_close", nc_close(md->ncid););

    free(of->md);
    open_file_delete(method->base_path, fd->name);
}

void adios_nc4_finalize(
        int mype,
        struct adios_method_struct *method)
{
//    struct adios_nc4_data_struct * md = (struct adios_nc4_data_struct*)method->method_data;
//    int myrank=md->rank;

//    if (md) {
//        Func_Timer("nc_close", nc_close(md->ncid););
//    }
//    md->group_comm = MPI_COMM_NULL;
//    md->ncid = -1;
//    md->root_ncid = -1;
//    md->rank = -1;
//    md->size = 0;

    if (adios_nc4_initialized)
        adios_nc4_initialized = 0;
}

/*
 * Maps bp datatypes to h5 datatypes
 *
 * The Mapping is according to HDF5 Reference Manual
 * (http://hdf.ncsa.uiuc.edu/HDF5/doc1.6/Datatypes.html)
 */
int getNC4TypeId(
        enum ADIOS_DATATYPES type,
        nc_type *nc4_type_id,
        enum ADIOS_FLAG fortran_flag)
{
    int size, status=0;
    switch (type)
    {
    case adios_byte:
        *nc4_type_id = NC_BYTE;
        break;
    case adios_short:
        *nc4_type_id = NC_SHORT;
        break;
    case adios_integer:
        *nc4_type_id = NC_INT;
        break;
    case adios_long:
        *nc4_type_id = NC_INT64;
        break;
    case adios_string:
        *nc4_type_id = NC_CHAR;
        break;
    case adios_real:
        *nc4_type_id = NC_FLOAT;
        break;
    case adios_double:
        *nc4_type_id = NC_DOUBLE;
        break;
    case adios_long_double:
        fprintf(stderr, "ERROR in mapping ADIOS Data Types to NC4: long double not supported yet.\n");
        status = -1;
        break;
    case adios_complex:
    case adios_double_complex:
        fprintf(stderr, "ERROR in mapping ADIOS Data Types to NC4: complex not supported yet.\n");
        status = -1;
        break;
    case adios_unsigned_byte:
        *nc4_type_id = NC_UBYTE;
        break;
    case adios_unsigned_short:
        *nc4_type_id = NC_USHORT;
        break;
    case adios_unsigned_integer:
        *nc4_type_id = NC_UINT;
        break;
    case adios_unsigned_long:
        *nc4_type_id = NC_UINT64;
        break;
    default:
        status = -1;
        break;
    }
    return status;
}

int getTypeSize(
        enum ADIOS_DATATYPES type,
        void *val)
{
    switch (type)
    {
    case adios_byte:
    case adios_unsigned_byte:
        return 1;

    case adios_string:
        return strlen ((char *) val);

    case adios_short:
    case adios_unsigned_short:
        return 2;

    case adios_integer:
    case adios_unsigned_integer:
        return 4;

    case adios_real:
        return 4;

    case adios_long:
    case adios_unsigned_long:
        return 8;

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

int ncd_gen_name(
        char *fullname,
        char *path,
        char *name)
{
    int i;
    char *new_path = strdup(path);

    if (path[0] == '/') {
         new_path=new_path+1;
    }

    for (i = 0; i < strlen(new_path); i++) {
        if ( new_path[i] == '[' || new_path[i] == ']' || new_path[i] == '/' || new_path[i] == '\\')
            new_path[i] = '_';
    }
    if (*new_path != '\0') {
        if (new_path[i-1]!='_') {
            if (strcmp(name,"")) {
                sprintf (fullname, "%s_%s", new_path, name);
            } else {
                strcpy (fullname,new_path);
                fullname [strlen(fullname)] = '\0';
            }
        } else {
            if (strcmp(name,"")) {
                sprintf (fullname, "%s%s", new_path, name);
            } else {
                strcpy (fullname,new_path);
                fullname [strlen(fullname)] = '\0';
            }
        }
    } else {
        strcpy (fullname, name);
    }

    if (DEBUG>3) printf("fullname==%s\n", fullname);
    return 0;
}
