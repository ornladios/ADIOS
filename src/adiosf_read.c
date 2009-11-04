#include "../config.h"
#include <stdlib.h>
#include <string.h>
#include "bp_utils.h"
#include "bp_types.h"
#include "common_read.h"
#include "adios_error.h"
#include "futils.h"
#define BYTE_ALIGN 8

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif

/*********************/
/* FORTRAN INTERFACE */
/*********************/
void FC_FUNC_(adios_errmsg, adios_errmsg) (char *msg, int msg_len)
{
    futils_cstr_to_fstr( adios_get_last_errmsg(), (char *)msg, msg_len);
}

void FC_FUNC_(adios_fopen, ADIOS_FOPEN)
        (int64_t * fp,
         char    * fname,
         void    * fcomm,
         int     * groups_count,
         int     * err,
         int       fname_len)
{
    ADIOS_FILE *afp;
    char *namestr;
    MPI_Comm comm = MPI_Comm_f2c (*((int *) fcomm));
    futils_called_from_fortran_set();

    namestr = futils_fstr_to_cstr(fname, fname_len);
    if (namestr != NULL) {
        afp = common_read_fopen (namestr,comm);
        if (afp != NULL) {
            *groups_count = afp->groups_count;
        } else {
            *groups_count = 0;
        }
        *fp = (int64_t) afp;
        free(namestr);
    } else {
        *fp = (int64_t) NULL;
    }
    *err = -adios_errno;
    if (*err)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

void FC_FUNC_(adios_fclose, ADIOS_FCLOSE) (int64_t * fp, int * err)
{
    ADIOS_FILE *afp = (ADIOS_FILE *) *fp;
    *err = common_read_fclose (afp);
    futils_called_from_fortran_unset();
    if (*err)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

void FC_FUNC_(adios_inq_file, ADIOS_INQ_FILE) 
        (int64_t * fp,
         int     * vars_count,
         int     * attrs_count,
         int     * tstart,
         int     * ntsteps,
         void    * gnamelist,
         int     * err,
         int       gnamelist_len)
{
    ADIOS_FILE *afp = (ADIOS_FILE *) *fp;
    int i;
    if (afp != NULL) {
        *vars_count = afp->vars_count;
        *attrs_count = afp->attrs_count;
        *tstart = afp->tidx_start;
        *ntsteps = afp->ntimesteps;
        *err = 0;
        for (i=0;i<afp->groups_count;i++) {
            futils_cstr_to_fstr( afp->group_namelist[i], (char *)gnamelist+i*gnamelist_len, gnamelist_len);
        }
    } else {
        *vars_count = 0;
        *attrs_count = 0;
        *tstart = 0;
        *ntsteps = 0;
        *err = 1;
    }
}

void FC_FUNC_(adios_gopen, ADIOS_GOPEN) 
        (int64_t * fp,
         int64_t * gp,
         char    * grpname,
         int     * vars_count, 
         int     * attrs_count, 
         int     * err,
         int       grpname_len)
{
    char *namestr;
    ADIOS_GROUP *agp;
    ADIOS_FILE *afp = (ADIOS_FILE *) *fp;

    namestr = futils_fstr_to_cstr(grpname, grpname_len);
    if (namestr != NULL) {
        agp = common_read_gopen (afp, namestr);
        if (agp != NULL) {
            *vars_count = agp->vars_count;
            *attrs_count = agp->attrs_count;
        } else {
            *vars_count = 0;
            *attrs_count = 0;
        }
        *gp = (int64_t)agp;
        free(namestr);
    } else {
        *gp = (int64_t) NULL;
    }
    *err = -adios_errno;
    if (*err)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

void FC_FUNC_(adios_gclose, ADIOS_GCLOSE) (int64_t * gp, int * err)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    *err=common_read_gclose(agp);
    if (*err)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

void FC_FUNC_(adios_inq_group, ADIOS_INQ_GROUP)
        (int64_t * gp, 
         void    * vnamelist, 
         void    * anamelist,
         int     * err, 
         int       vnamelist_len, 
         int       anamelist_len) 
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    int i;
    if (agp != NULL) {
        for (i=0;i<agp->vars_count;i++) {
            futils_cstr_to_fstr( agp->var_namelist[i], (char *)vnamelist+i*vnamelist_len, vnamelist_len);
        } 
        for (i=0;i<agp->attrs_count;i++) {
            futils_cstr_to_fstr( agp->attr_namelist[i], (char *)anamelist+i*anamelist_len, anamelist_len);
        } 
        *err = 0;
    } else {
        *err = 1;
    }
}

void FC_FUNC_(adios_inq_var, ADIOS_INQ_VAR) 
        (int64_t  * gp, 
         char     * varname,
         int      * type,
         int      * ndim,
         uint64_t * dims,
         int      * timedim,
         int      * err,
         int        varname_len)
{
    char *varstr;
    int  i;
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    ADIOS_VARINFO *vi = NULL;

    varstr = futils_fstr_to_cstr(varname, varname_len);
    if (varstr != NULL) {
        vi = common_read_inq_var (agp, varstr);
    }
    if (vi != NULL) {
        *type = vi->type;
        *ndim = vi->ndim;
        *timedim = vi->timedim;
        for (i=0;i<vi->ndim;i++)
            dims[i] = vi->dims[i];
        common_read_free_varinfo(vi);
    } else {
        *type = adios_unknown;
        *ndim = 0;
        *timedim = -1;
    }
    *err = -adios_errno;
    if (*err < 0)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

void FC_FUNC_(adios_read_var, ADIOS_READ_VAR) 
        (int64_t  * gp,
         char     * varname,
         uint64_t * start,
         uint64_t * count,
         void     * data,
         int64_t  * read_bytes,
         int varname_len)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    char *varstr;
    int i;
    varstr = futils_fstr_to_cstr(varname, varname_len);
    if (varstr != NULL) {
        *read_bytes = common_read_read_var (agp, varstr, start, count, data);
        free(varstr);
    } else {
        *read_bytes = -adios_errno;
    }
    if (*read_bytes < 0)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

/* Specific function for each data type */
void FC_FUNC_(adios_read_var_int1, ADIOS_READ_VAR_INT1) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_int2, ADIOS_READ_VAR_INT2) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_int4, ADIOS_READ_VAR_INT4) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_int8, ADIOS_READ_VAR_INT8) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_real4, ADIOS_READ_VAR_REAL4) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_real8, ADIOS_READ_VAR_REAL8) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_char, ADIOS_READ_VAR_CHAR) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_complex8, ADIOS_READ_VAR_COMPLEX8) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_complex16, ADIOS_READ_VAR_COMPLEX16) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_logical1, ADIOS_READ_VAR_LOGICAL1) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_logical2, ADIOS_READ_VAR_LOGICAL2) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_logical4, ADIOS_READ_VAR_LOGICAL4) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }
void FC_FUNC_(adios_read_var_logical8, ADIOS_READ_VAR_LOGICAL8) (int64_t * gp, char * varname, uint64_t * start, uint64_t * count, void * data, int64_t * read_bytes, int varname_len) { FC_FUNC_(adios_read_var, ADIOS_READ_VAR) (gp, varname, start, count, data, read_bytes, varname_len); }

void FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) 
        (int64_t * gp,
         char    * varname,
         void    * value,
         void    * gmin,
         void    * gmax,
         void    * mins,
         void    * maxs,
         int     * err,
         int       varname_len)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    ADIOS_VARINFO *vi = NULL;
    char *varstr;
    int i, size, ntime;

    varstr = futils_fstr_to_cstr(varname, varname_len);
    if (varstr != NULL) {
        vi = common_read_inq_var (agp, varstr);
    }
    if (vi != NULL) {
        size = bp_get_type_size(vi->type, vi->value);
        if (vi->type == adios_string) size++;
        if (vi->timedim > -1)
            ntime = agp->fp->ntimesteps;
        else 
            ntime = 1;
        if (vi->value) memcpy(value, vi->value, size);
        if (vi->gmin) memcpy(gmin, vi->gmin, size);
        if (vi->gmax) memcpy(gmax, vi->gmax, size);
        if (vi->mins) {
            for (i=0; i<ntime; i++)
                memcpy((char *)mins+i*size, (char *)(vi->mins)+i*size, size);
        }
        if (vi->maxs) {
            for (i=0; i<ntime; i++)
                memcpy((char *)maxs+i*size, (char *)(vi->maxs)+i*size, size);
        }
        common_read_free_varinfo(vi);
    }
    *err = -adios_errno;
    if (*err < 0)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

/* Specific function for each data type */
void FC_FUNC_(adios_get_varminmax_int1, ADIOS_GET_VARMINMAX_INT1) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_int2, ADIOS_GET_VARMINMAX_INT2) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_int4, ADIOS_GET_VARMINMAX_INT4) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_int8, ADIOS_GET_VARMINMAX_INT8) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_real4, ADIOS_GET_VARMINMAX_REAL4) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_real8, ADIOS_GET_VARMINMAX_REAL8) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_complex8, ADIOS_GET_VARMINMAX_COMPLEX8) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_complex16, ADIOS_GET_VARMINMAX_COMPLEX16) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_char, ADIOS_GET_VARMINMAX_CHAR) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_logical1, ADIOS_GET_VARMINMAX_LOGICAL1) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_logical2, ADIOS_GET_VARMINMAX_LOGICAL2) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_logical4, ADIOS_GET_VARMINMAX_LOGICAL4) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }
void FC_FUNC_(adios_get_varminmax_logical8, ADIOS_GET_VARMINMAX_LOGICAL8) (int64_t * gp, char * varname, void * value, void * gmin, void * gmax, void * mins, void * maxs, int * err, int varname_len) { FC_FUNC_(adios_get_varminmax, ADIOS_GET_VARMINMAX) (gp, varname, value, gmin, gmax, mins, maxs, err, varname_len); }

void FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) 
        (int64_t * gp,
         char    * attrname,
         void    * attr,
         int     * err,
         int       attrname_len)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    char *attrstr;
    int i;
    void *data;
    int size;
    enum ADIOS_DATATYPES type;
    attrstr = futils_fstr_to_cstr(attrname, attrname_len);
    if (attrstr != NULL) {
        *err = common_read_get_attr (agp, attrstr, &type, &size, &data);
        if (data) {
            memcpy(attr, data, size);
            free(data);
        }
    } else {
        *err = -adios_errno;
    }
    if (*err < 0)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

void FC_FUNC_(adios_inq_attr, ADIOS_INQ_ATTR) 
        (int64_t * gp,
         char    * attrname,
         int     * type,
         int     * size,
         int     * err,
         int       attrname_len)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    char *attrstr;
    int i;
    void *data;
    attrstr = futils_fstr_to_cstr(attrname, attrname_len);
    if (attrstr != NULL) {
        *err = common_read_get_attr (agp, attrstr, (enum ADIOS_DATATYPES *)type, size, &data);
        free(data);
    } else {
        *err = -adios_errno;
    }
    if (*err < 0)
        fprintf(stderr, "Error: %s\n", adios_get_last_errmsg());
}

/* Specific function for each data type */
void FC_FUNC_(adios_get_attr_int1, ADIOS_GET_ATTR_INT1) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_int2, ADIOS_GET_ATTR_INT2) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_int4, ADIOS_GET_ATTR_INT4) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_int8, ADIOS_GET_ATTR_INT8) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_real4, ADIOS_GET_ATTR_REAL4) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_real8, ADIOS_GET_ATTR_REAL8) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_complex8, ADIOS_GET_ATTR_COMPLEX8) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_complex16, ADIOS_GET_ATTR_COMPLEX16) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_char, ADIOS_GET_ATTR_CHAR) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_logical1, ADIOS_GET_ATTR_LOGICAL1) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_logical2, ADIOS_GET_ATTR_LOGICAL2) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_logical4, ADIOS_GET_ATTR_LOGICAL4) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
void FC_FUNC_(adios_get_attr_logical8, ADIOS_GET_ATTR_LOGICAL8) (int64_t * gp, char * attrname, void * attr, int * err, int attrname_len) { FC_FUNC_(adios_get_attr, ADIOS_GET_ATTR) (gp, attrname, attr, err, attrname_len); }
