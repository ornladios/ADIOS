#include "../config.h"
#include <stdlib.h>
#include <string.h>
#include "bp_utils.h"
#include "bp_types.h"
#include "adios_read.h"
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
void FC_FUNC_(adiosf_errmsg, ADIOSF_LASTERRMSG) (char *msg, int msg_len)
{
    futils_cstr_to_fstr( adios_errmsg(), (char *)msg, msg_len);
}

void FC_FUNC_(adiosf_fopen, ADIOSF_FOPEN)
    (int64_t * fp,
                  char * fname,
                  void * fcomm,
                  int * err,
                  int fname_len
                 )
{
    ADIOS_FILE *afp;
    char *namestr;
    MPI_Comm comm = MPI_Comm_f2c (*((int *) fcomm));
    futils_called_from_fortran_set();

    namestr = futils_fstr_to_cstr(fname, fname_len);
    if (namestr != NULL) {
        afp = adios_fopen (namestr,comm);
        *fp = (int64_t) afp;
        free(namestr);
    } else {
        *fp = (int64_t) NULL;
    }
    *err = -adios_errno;
    if (*err)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
}

void FC_FUNC_(adiosf_fclose, ADIOSF_FCLOSE)( int64_t * fp, int * err)
{
    ADIOS_FILE *afp = (ADIOS_FILE *) *fp;
    *err = adios_fclose (afp);
    futils_called_from_fortran_unset();
    if (*err)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
}

void FC_FUNC_(adiosf_inq_file, ADIOSF_INQ_FILE) ( int64_t * fp,
                       int * groups_count,
                       int * vars_count,
                       int * attrs_count,
                       int * tstart,
                       int * ntsteps,
                       void * gnamelist,
                       int * err,
                       int gnamelist_len)
{
    ADIOS_FILE *afp = (ADIOS_FILE *) *fp;
    int i;
    *groups_count = afp->groups_count;
    *vars_count = afp->vars_count;
    *attrs_count = afp->attrs_count;
    *tstart = afp->tidx_start;
    *ntsteps = afp->ntimesteps;
    *err = 0;
    for (i=0;i<*groups_count;i++) {
        futils_cstr_to_fstr( afp->group_namelist[i], (char *)gnamelist+i*gnamelist_len, gnamelist_len);
    }
}

void FC_FUNC_(adiosf_gopen, ADIOSF_GOPEN) 
    ( int64_t * fp,
                    int64_t * gp,
                    char * grpname,
                    int * err,
                    int grpname_len)
{
    char *namestr;
    ADIOS_GROUP *agp;
    ADIOS_FILE *afp = (ADIOS_FILE *) *fp;

    namestr = futils_fstr_to_cstr(grpname, grpname_len);
    if (namestr != NULL) {
        agp = adios_gopen (afp, namestr);
        *gp = (int64_t)agp;
        free(namestr);
    } else {
        *gp = (int64_t) NULL;
    }
    *err = -adios_errno;
    if (*err)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
}

void FC_FUNC_(adiosf_gclose, ADIOSF_GCLOSE)( int64_t * gp, int * err)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    *err=adios_gclose(agp);
    if (*err)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
}

void FC_FUNC_(adiosf_inq_group, ADIOSF_INQ_GROUP)
    (int64_t * gp, int *vcnt, void *vnamelist, int *acnt, void *anamelist,
        int *err, int vnamelist_len, int anamelist_len) 
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    int i;
    *vcnt = agp->vars_count;
    for (i=0;i<*vcnt;i++) {
        futils_cstr_to_fstr( agp->var_namelist[i], (char *)vnamelist+i*vnamelist_len, vnamelist_len);
    } 
    *acnt = agp->attrs_count;
    for (i=0;i<*acnt;i++) {
        futils_cstr_to_fstr( agp->attr_namelist[i], (char *)anamelist+i*anamelist_len, anamelist_len);
    } 
    *err = 0;
}

void FC_FUNC_(adiosf_inq_var, ADIOSF_INQ_VAR) 
    (int64_t  * gp, char * varname,
                     int      * type,
                     int      * ndim,
                     uint64_t * dims,
                     int      * timedim,
                     int      * err,
                     int varname_len)
{
    char *varstr;
    int  i;
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    ADIOS_VARINFO *vi = NULL;

    varstr = futils_fstr_to_cstr(varname, varname_len);
    if (varstr != NULL) {
        vi = adios_inq_var (agp, varstr);
    }
    if (vi != NULL) {
        *type = vi->type;
        *ndim = vi->ndim;
        *timedim = vi->timedim;
        for (i=0;i<vi->ndim;i++)
            dims[i] = vi->dims[i];
        adios_free_varinfo(vi);
    } else {
        *type = adios_unknown;
        *ndim = 0;
        *timedim = -1;
    }
    *err = -adios_errno;
    if (*err < 0)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
}

void FC_FUNC_(adiosf_read_var, ADIOSF_READ_VAR) 
    (int64_t  * gp,
                      char     * varname,
                      uint64_t * start,
                      uint64_t * readsize,
                      void     * data,
                      int64_t  * read_bytes,
                      int varname_len)
{
    /* FIXME: Magically, *gh becomes 0 after the C function call, which causes abort in a next call.
       Temporarily we save its value and reassign it but clearly it must be found out why this is
       happening. */
    int64_t tmp=*gp;
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    char *varstr;
    int i;
    varstr = futils_fstr_to_cstr(varname, varname_len);
    if (varstr != NULL) {
        *read_bytes = adios_read_var (agp, varstr, start, readsize, data);
        free(varstr);
    } else {
        *read_bytes = -adios_errno;
    }
    if (*read_bytes < 0)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
    *gp=tmp;
}

void FC_FUNC_(adiosf_get_varminmax, ADIOSF_GET_VARMINMAX) 
    (int64_t * gp,
                           char    * varname,
                           void    * value,
                           void    * gmin,
                           void    * gmax,
                           void    * mins,
                           void    * maxs,
                           int     * err,
                           int varname_len)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    ADIOS_VARINFO *vi = NULL;
    char *varstr;
    int i, size, ntime;

    varstr = futils_fstr_to_cstr(varname, varname_len);
    if (varstr != NULL) {
        vi = adios_inq_var (agp, varstr);
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
        adios_free_varinfo(vi);
    }
    *err = -adios_errno;
    if (*err < 0)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
}

void FC_FUNC_(adiosf_get_attr, ADIOSF_GET_ATTR) 
    (int64_t * gp
                     ,char * attrname
                     ,void * attr
                     ,int * err
                     ,int attrname_len)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    char *attrstr;
    int i;
    void *data;
    int size;
    enum ADIOS_DATATYPES type;
    attrstr = futils_fstr_to_cstr(attrname, attrname_len);
    if (attrstr != NULL) {
        *err = adios_get_attr (agp, attrstr, &type, &size, &data);
        if (data) {
            memcpy(attr, data, size);
            free(data);
        }
    } else {
        *err = -adios_errno;
    }
    if (*err < 0)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
}

void FC_FUNC_(adiosf_inq_attr, ADIOSF_INQ_ATTR) 
    (int64_t * gp
                     ,char * attrname
                     ,int * type
                     ,int * size
                     ,int * err
                     ,int attrname_len)
{
    ADIOS_GROUP *agp = (ADIOS_GROUP *) *gp;
    char *attrstr;
    int i;
    void *data;
    attrstr = futils_fstr_to_cstr(attrname, attrname_len);
    if (attrstr != NULL) {
        *err = adios_get_attr (agp, attrstr, (enum ADIOS_DATATYPES *)type, size, &data);
        free(data);
    } else {
        *err = -adios_errno;
    }
    if (*err < 0)
        fprintf(stderr, "Error: %s\n", adios_errmsg());
}

