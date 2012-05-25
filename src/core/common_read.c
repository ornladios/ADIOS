/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "public/adios_error.h"
#include "core/common_read.h"
#include "core/adios_read_hooks.h"
#include "core/futils.h"
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


int common_read_init_method (enum ADIOS_READ_METHOD method,
                             MPI_Comm comm,
                             const char * parameters)
{
    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, 
            "Invalid read method (=%d) passed to adios_read_init_method().", (int)method);
        return -err_invalid_read_method;
    } 
    // init the adios_read_hooks_struct if not yet initialized  
    adios_read_hooks_init (&adios_read_hooks); 

    return adios_read_hooks[method].adios_init_method_fn (comm, parameters);
}


int common_read_finalize_method(enum ADIOS_READ_METHOD method)
{
    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, 
            "Invalid read method (=%d) passed to adios_read_finalize_method().", (int)method);
        return -err_invalid_read_method;
    } 

    return adios_read_hooks[method].adios_finalize_method_fn ();
}


ADIOS_FILE * common_read_open_stream (const char * fname, 
                                      enum ADIOS_READ_METHOD method, 
                                      MPI_Comm comm, 
                                      enum ADIOS_LOCKMODE lock_mode, 
                                      int timeout_msec)
{
    ADIOS_FILE * fp;
    struct common_read_internals_struct * internals; 

    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, 
            "Invalid read method (=%d) passed to adios_read_open_stream().", (int)method);
        return NULL;
    } 

    adios_errno = 0;
    internals = (struct common_read_internals_struct *) 
                    calloc(1,sizeof(struct common_read_internals_struct));
    // init the adios_read_hooks_struct if not yet initialized 
    adios_read_hooks_init (&adios_read_hooks); 

    internals->method = selected_method;
    internals->read_hooks = adios_read_hooks;

    fp = adios_read_hooks[internals->method].adios_open_stream_fn (fname, comm, lock_mode, timeout_msec);

    // save the method in fp->internal_data
    if (fp)
        fp->internal_data = (void *)internals;
    return fp;
}


ADIOS_FILE * common_read_open_file (const char * fname, 
                                    enum ADIOS_READ_METHOD method,
                                    MPI_Comm comm)
{
    ADIOS_FILE * fp;
    struct common_read_internals_struct * internals; 

    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, 
            "Invalid read method (=%d) passed to adios_read_open_file().", (int)method);
        return NULL;
    } 

    adios_errno = 0;
    internals = (struct common_read_internals_struct *) 
                    calloc(1,sizeof(struct common_read_internals_struct));
    // init the adios_read_hooks_struct if not yet initialized 
    adios_read_hooks_init (&adios_read_hooks); 

    internals->method = selected_method;
    internals->read_hooks = adios_read_hooks;

    fp = adios_read_hooks[internals->method].adios_open_file_fn (fname, comm);

    // save the method in fp->internal_data
    if (fp)
        fp->internal_data = (void *)internals;
    return fp;
}

int common_read_close (const ADIOS_FILE *fp) 
{
    struct common_read_internals_struct * internals;
    int retval;

    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_close_fn (fp);
        free (internals);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_read_close()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}

void common_read_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{
    struct common_read_internals_struct * internals;

    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        internals->read_hooks[internals->method].adios_reset_dimension_order_fn (fp, is_fortran);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_reset_dimension_order()");
    }
}


int common_read_advance_step (const ADIOS_FILE *fp, int last, int wait_for_step)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_advance_step_fn (fp, last, wait_for_step);
        free (internals);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_advance_step()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}


void common_read_release_step (const ADIOS_FILE *fp) 
{
    struct common_read_internals_struct * internals;

    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        internals->read_hooks[internals->method].adios_release_step_fn (fp);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_reset_dimension_order()");
    }
}


ADIOS_VARINFO * common_read_inq_var (const ADIOS_FILE *fp, const char * varname) 
{
    struct common_read_internals_struct * internals;
    ADIOS_VARINFO * retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_inq_var_fn (fp, varname);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var()");
        retval = NULL;
    }
    return retval;
}


ADIOS_VARINFO * common_read_inq_var_byid (const ADIOS_FILE *fp, int varid)
{
    struct common_read_internals_struct * internals;
    ADIOS_VARINFO * retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_inq_var_byid_fn (fp, varid);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var_byid()");
        retval = NULL;
    }
    return retval;
}


int common_read_inq_var_stat (const ADIOS_FILE *fp, const ADIOS_VARINFO * varinfo,
                             int per_step_stat, int per_block_stat)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_inq_var_stat_fn (fp, varinfo, per_step_stat, per_block_stat);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var_byid()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}


#define MYFREE(p) {free(p); p=NULL;}
void common_read_free_varinfo (ADIOS_VARINFO *vp)
{
    int i;
    if (vp) {
        if (vp->blockinfo) {
            ADIOS_VARBLOCK *bp = vp->blockinfo;
            for (i=0; i<vp->sum_nblocks; i++) {
                if (bp->start) MYFREE (bp->start);
                if (bp->count) MYFREE (bp->count);
                bp++;
            }
            MYFREE(vp->blockinfo);
        }

        if (vp->statistics) {
            ADIOS_VARSTAT *sp = vp->statistics;
            if (sp->min && sp->min != vp->value)   MYFREE(sp->min);
            if (sp->max && sp->max != vp->value)   MYFREE(sp->max);
            if (sp->avg && sp->avg != vp->value)   MYFREE(sp->avg);
            if (sp->std_dev)                       MYFREE(sp->std_dev);

            if (sp->steps) {
                if (sp->steps->mins)        MYFREE(sp->steps->mins);
                if (sp->steps->maxs)        MYFREE(sp->steps->maxs);
                if (sp->steps->avgs)        MYFREE(sp->steps->avgs);
                if (sp->steps->std_devs)    MYFREE(sp->steps->std_devs);
                MYFREE(sp->steps);
            }

            if (sp->blocks) {
                if (sp->blocks->mins)        MYFREE(sp->blocks->mins);
                if (sp->blocks->maxs)        MYFREE(sp->blocks->maxs);
                if (sp->blocks->avgs)        MYFREE(sp->blocks->avgs);
                if (sp->blocks->std_devs)    MYFREE(sp->blocks->std_devs);
                MYFREE(sp->blocks);
            }

            if (sp->histogram) {
                if (sp->histogram->breaks)        MYFREE(sp->histogram->breaks);
                if (sp->histogram->frequencies)   MYFREE(sp->histogram->frequencies);
                if (sp->histogram->gfrequencies)  MYFREE(sp->histogram->gfrequencies);
                MYFREE(sp->histogram);
            }

            MYFREE(vp->statistics);
        }

        if (vp->dims)    MYFREE(vp->dims);
        if (vp->value)   MYFREE(vp->value);
        if (vp->nblocks) MYFREE(vp->nblocks);

        free(vp);
    }
}


int common_read_schedule_read (const ADIOS_FILE      * fp,
                               const ADIOS_SELECTION * sel,
                               const char            * varname,
                               int                     from_steps,
                               int                     nsteps,
                               void                  * data)

{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_schedule_read_fn (fp, sel, varname, from_steps, nsteps, data);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_schedule_read()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}


int common_read_schedule_read_byid (const ADIOS_FILE      * fp,
        const ADIOS_SELECTION * sel,
        int                     varid,
        int                     from_steps,
        int                     nsteps,
        void                  * data)

{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_schedule_read_byid_fn (fp, sel, varid, from_steps, nsteps, data);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_schedule_read_byid()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}


int common_read_perform_reads (const ADIOS_FILE *fp, int blocking)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_perform_reads_fn (fp, blocking);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_perform_reads()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}


int common_read_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_check_reads_fn (fp, chunk);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_check_reads()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}


void common_read_free_chunk (const ADIOS_VARCHUNK *chunk)
{
    /** Free the memory of a chunk allocated inside adios_check_reads().
     * It only frees the ADIOS_VARCHUNK struct and the ADIOS_SELECTION struct
     * pointed by the chunk. The data pointer should never be freed since
     * that memory belongs to the reading method.
     */
     if (chunk) {
        if (chunk->sel) {
            free(chunk->sel);
            chunk->sel = NULL;
        }
        free(chunk);
     }
}


int common_read_get_attr (const ADIOS_FILE * fp, 
                          const char * attrname, 
                          enum ADIOS_DATATYPES * type,
                          int * size, 
                          void ** data)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_get_attr_fn (fp, attrname, type, size, data);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_read_get_attr()");
        retval = -err_invalid_file_pointer;
    }
    return retval;
}


int common_read_get_attr_byid (const ADIOS_FILE * fp, 
                               int attrid, 
                               enum ADIOS_DATATYPES * type, 
                               int * size, 
                               void ** data)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_get_attr_byid_fn (fp, attrid, type, size, data);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_read_get_attr_byid()");
        retval = -err_invalid_file_pointer;
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


int common_read_get_grouplist (const ADIOS_FILE  *fp, char **group_namelist)
{
}


int common_read_group_view (const ADIOS_FILE  *fp, int groupid)
{
}


void common_read_print_fileinfo (const ADIOS_FILE *fp) 
{
    int i;
    int ngroups;
    char **group_namelist;
    ngroups = common_read_get_grouplist (fp, group_namelist);

    printf ("---------------------------\n");
    printf ("     file information\n");
    printf ("---------------------------\n");
    printf ("  # of groups:     %d\n"
            "  # of variables:  %d\n"
            "  # of attributes: %d\n"
            "  current step:    %d\n"
            "  last step:       %d\n",
            ngroups,
            fp->nvars,
            fp->nattrs,
            fp->current_step,
            fp->last_step);
    printf ("---------------------------\n");
    printf ("     var information\n");
    printf ("---------------------------\n");
    printf ("    var id\tname\n");
    if (fp->var_namelist) {
        for (i=0; i<fp->nvars; i++)
            printf("\t%d)\t%s\n", i, fp->var_namelist[i]);
    }
    printf ("---------------------------\n");
    printf ("     attribute information\n");
    printf ("---------------------------\n");
    printf ("    attr id\tname\n");
    if (fp->attr_namelist) {
        for (i=0; i<fp->nattrs; i++)
            printf("\t%d)\t%s\n", i, fp->attr_namelist[i]);
    }
    printf ("---------------------------\n");
    printf ("     group information\n");
    printf ("---------------------------\n");
    if (group_namelist) {
        for (i=0; i<ngroups; i++)
            printf("\t%d)\t%s\n", i, group_namelist[i]);
    }


    return;
}



