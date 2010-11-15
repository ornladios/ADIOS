/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adios.h"
#include "common_read.h"
#include "adios_read.h"
#include "adios_error.h"
#include "adios_read_hooks.h"
#include "futils.h"
#define BYTE_ALIGN 8

#ifdef DMALLOC
#include "dmalloc.h"
#endif


/* Note: MATLAB reloads the mex64 files each time, so all static variables get the original value.
   Therefore static variables cannot be used to pass info between two Matlab/ADIOS calls */
static struct adios_read_hooks_struct * adios_read_hooks = 0;
static enum ADIOS_READ_METHOD selected_method = ADIOS_READ_METHOD_BP;

struct common_read_internals_struct {
    enum ADIOS_READ_METHOD method;
    struct adios_read_hooks_struct * read_hooks; /* Save adios_read_hooks for each fopen for Matlab */
};


int common_read_set_read_method(enum ADIOS_READ_METHOD method)
{
    int retval = 0;
    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, "Invalid read method (=%d) passed to adios_set_read_method().", (int)method);
        retval = -err_invalid_read_method;
    } else {
        selected_method = method;
    }
    return retval;
}


int common_read_init(MPI_Comm comm)
{
    adios_read_hooks_init (&adios_read_hooks); // init the adios_read_hooks_struct if not yet initialized    
    return adios_read_hooks[selected_method].adios_init_fn (comm);
}

int common_read_finalize()
{
    return adios_read_hooks[selected_method].adios_finalize_fn ();
}

ADIOS_FILE * common_read_fopen (const char * fname, MPI_Comm comm)
{
    ADIOS_FILE * fp;
    struct common_read_internals_struct * internals = 
            (struct common_read_internals_struct *) calloc(1,sizeof(struct common_read_internals_struct));

    adios_errno = 0;
    adios_read_hooks_init (&adios_read_hooks); // init the adios_read_hooks_struct if not yet initialized    

    internals->method = selected_method;
    internals->read_hooks = adios_read_hooks;

    fp = adios_read_hooks[internals->method].adios_fopen_fn (fname, comm);


    // save the method in fp->internal_data
    if (fp)
        fp->internal_data = (void *)internals;
    return fp;
}

int common_read_fclose (ADIOS_FILE *fp) 
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_fclose_fn (fp);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_fclose()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}

void common_read_reset_dimension_order (ADIOS_FILE *fp, int is_fortran)
{
    struct common_read_internals_struct * internals;

    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        internals->read_hooks[internals->method].adios_reset_dimension_order_fn (fp);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_reset_dimension_order()");
    }
}

ADIOS_GROUP * common_read_gopen (ADIOS_FILE *fp, const char * grpname)
{
    struct common_read_internals_struct * internals;
    ADIOS_GROUP * retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_gopen_fn (fp, grpname);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_gopen()");
        retval = NULL;
    }
    return retval;
}

ADIOS_GROUP * common_read_gopen_byid (ADIOS_FILE *fp, int grpid)
{
    struct common_read_internals_struct * internals;
    ADIOS_GROUP * retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_gopen_byid_fn (fp, grpid);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_gopen()");
        retval = NULL;
    }
    //printf("%s: gp=%x, gp->fp=%x, gp->gh=%x\n",__func__, retval, retval->fp, retval->gh);
    return retval;
}
                   
int common_read_gclose (ADIOS_GROUP *gp)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (gp) {
        //printf("%s: gp=%x, gp->fp=%x, gp->gh=%x\n",__func__, gp, gp->fp, gp->gh);
        internals = (struct common_read_internals_struct *) gp->fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_gclose_fn (gp);
    } else {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_gclose()");
        retval = -err_invalid_group_struct;
    }
    return retval;
}

int common_read_get_attr (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type,
                    int * size, void ** data)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (gp) {
        internals = (struct common_read_internals_struct *) gp->fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_get_attr_fn (gp, attrname, type, size, data);
    } else {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_get_attr()");
        retval = -err_invalid_group_struct;
    }
    return retval;
}

int common_read_get_attr_byid (ADIOS_GROUP * gp, int attrid, 
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (gp) {
        internals = (struct common_read_internals_struct *) gp->fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_get_attr_byid_fn (gp, attrid, type, size, data);
    } else {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_get_attr_byid()");
        retval = -err_invalid_group_struct;
    }
    return retval;
}

ADIOS_VARINFO * common_read_inq_var (ADIOS_GROUP *gp, const char * varname) 
{
    struct common_read_internals_struct * internals;
    ADIOS_VARINFO * retval;
    
    adios_errno = 0;
    if (gp) {
        internals = (struct common_read_internals_struct *) gp->fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_inq_var_fn (gp, varname);
    } else {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_inq_var_byid()");
        retval = NULL;
    }
    return retval;
}

ADIOS_VARINFO * common_read_inq_var_byid (ADIOS_GROUP *gp, int varid)
{
    struct common_read_internals_struct * internals;
    ADIOS_VARINFO * retval;
    
    adios_errno = 0;
    if (gp) {
        internals = (struct common_read_internals_struct *) gp->fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_inq_var_byid_fn (gp, varid);
    } else {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_inq_var_byid()");
        retval = NULL;
    }
    return retval;
}

void common_read_free_varinfo (ADIOS_VARINFO *vp)
{
    if (vp) {
        if (vp->dims)   free(vp->dims);
        if (vp->value)  free(vp->value);
        if (vp->gmin && vp->gmin != vp->value)   free(vp->gmin);
        if (vp->gmax && vp->gmax != vp->value)   free(vp->gmax);
        //if (vp->mins)   free(vp->mins);
        //if (vp->maxs)   free(vp->maxs);
        free(vp);
    }
}


int64_t common_read_read_var (ADIOS_GROUP * gp, const char * varname,
                        const uint64_t * start, const uint64_t * count,
                        void * data)
{
    struct common_read_internals_struct * internals;
    int64_t retval;
    
    adios_errno = 0;
    if (gp) {
        internals = (struct common_read_internals_struct *) gp->fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_read_var_fn (gp, varname, start, count, data);
    } else {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        retval = -err_invalid_group_struct;
    }
    return retval;
}

int64_t common_read_read_local_var (ADIOS_GROUP    * gp,
                                    const char     * varname,
                                    int            idx,
                                    const uint64_t * start,
                                    const uint64_t * count,
                                    void           * data)
{
    struct common_read_internals_struct * internals;
    int64_t retval;

    adios_errno = 0;
    if (gp) {
        internals = (struct common_read_internals_struct *) gp->fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_read_local_var_fn (gp, varname, idx, start, count, data);
    } else {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_local_var()");
        retval = -err_invalid_group_struct;
    }
    return retval;
}

int64_t common_read_read_var_byid (ADIOS_GROUP    * gp,
                             int              varid,
                             const uint64_t  * start,
                             const uint64_t  * count,
                             void           * data)
{
    struct common_read_internals_struct * internals;
    int64_t retval;
    
    adios_errno = 0;
    if (gp) {
        internals = (struct common_read_internals_struct *) gp->fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_read_var_byid_fn (gp, varid, start, count, data);
    } else {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var_byid()");
        retval = -err_invalid_group_struct;
    }
    return retval;
}

const char * common_read_type_to_string (enum ADIOS_DATATYPES type)
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

int common_read_type_size(enum ADIOS_DATATYPES type, void *data)
{
    return bp_get_type_size(type, data);
}


void common_read_print_groupinfo (ADIOS_GROUP *gp) 
{
    int i;
    printf ("---------------------------\n");
    printf ("     var information\n");
    printf ("---------------------------\n");
    printf ("    var id\tname\n");
    if (gp->var_namelist) {
        for (i=0; i<gp->vars_count; i++)
            printf("\t%d)\t%s\n", i, gp->var_namelist[i]);
    }
    printf ("---------------------------\n");
    printf ("     attribute information\n");
    printf ("---------------------------\n");
    printf ("    attr id\tname\n");
    if (gp->attr_namelist) {
        for (i=0; i<gp->attrs_count; i++)
            printf("\t%d)\t%s\n", i, gp->attr_namelist[i]);
    }
    return;
}


void common_read_print_fileinfo (ADIOS_FILE *fp) 
{
    int i;
    printf ("---------------------------\n");
    printf ("     group information\n");
    printf ("---------------------------\n");
    printf ("\t# of groups:\t%d\n"
        "\t# of variables:\t%d\n"
        "\t# of attributes:%d\n"
        "\t# of timesteps:\t%d starting from %d\n",
        fp->groups_count,
        fp->vars_count,
        fp->attrs_count,
        fp->ntimesteps,
        fp->tidx_start);
    printf ("\t----------------\n");
    printf ("\tgroup id\tname\n");
    if (fp->group_namelist) {
        for (i=0; i<fp->groups_count; i++)
            printf("\t  %d)\t%s\n", i, fp->group_namelist[i]);
    }
    return;
}


