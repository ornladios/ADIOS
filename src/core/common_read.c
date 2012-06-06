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
#include "core/bp_utils.h" // struct namelists_struct
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
    
    /* Group view information *//* Actual method provides the group names */
    int     ngroups;
    char ** group_namelist;
    int   * nvars_per_group;     /* # of variables per each group */
    int   * nattrs_per_group;    /* # of attributes per each group */
    int     group_in_view;       /* 0..ngroups-1: selected group in view,
                                  -1: all groups */
    int     group_varid_offset;  /* offset of var IDs from specific group to full list
                                    if a selected group is in view */
    int     group_attrid_offset;
    int     full_nvars;          /* fp->nvars to save here for a group view */
    char ** full_varnamelist;    /* fp->var_namelist to save here if one group is viewed */
    int     full_nattrs;         /* fp->nvars to save here for a group view */
    char ** full_attrnamelist;   /* fp->attr_namelist to save here if one group is viewed */
};


int common_read_init_method (enum ADIOS_READ_METHOD method,
                             MPI_Comm comm,
                             const char * parameters)
{
    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, 
            "Invalid read method (=%d) passed to adios_read_init_method().", (int)method);
        return err_invalid_read_method;
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
        return err_invalid_read_method;
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

    // save the method and group information in fp->internal_data
    if (fp){
        adios_read_hooks[internals->method].adios_get_groupinfo_fn (fp, &internals->ngroups, 
                &internals->group_namelist, &internals->nvars_per_group, &internals->nattrs_per_group);
        internals->group_in_view = -1;
        internals->group_varid_offset = 0;
        internals->group_attrid_offset = 0;
        fp->internal_data = (void *)internals;
    } else {
        free (internals);
    }
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

    // save the method and group information in fp->internal_data
    if (fp){
        adios_read_hooks[internals->method].adios_get_groupinfo_fn (fp, &internals->ngroups, 
                &internals->group_namelist, &internals->nvars_per_group, &internals->nattrs_per_group);
        internals->group_in_view = -1;
        internals->group_varid_offset = 0;
        internals->group_attrid_offset = 0;
        fp->internal_data = (void *)internals;
    } else {
        free (internals);
    }
    return fp;
}

int common_read_close (ADIOS_FILE *fp) 
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
        retval = err_invalid_file_pointer;
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


int common_read_advance_step (ADIOS_FILE *fp, int last, int wait_for_step)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_advance_step_fn (fp, last, wait_for_step);
        /* Update group information too */
        adios_read_hooks[internals->method].adios_get_groupinfo_fn (fp, &internals->ngroups, 
                &internals->group_namelist, &internals->nvars_per_group, &internals->nattrs_per_group);
        if (internals->group_in_view > -1) {
            /* if we have a group view, we need to update the presented list again */
            /* advance_step updated fp->nvars, nattrs, var_namelist, attr_namelist */
            int groupid = internals->group_in_view;
            internals->group_in_view = -1; // we have the full view at this moment 
            common_read_group_view (fp, groupid);
        }
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_advance_step()");
        retval = err_invalid_file_pointer;
    }
    return retval;
}


void common_read_release_step (ADIOS_FILE *fp) 
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


static int common_read_find_name (int n, char ** namelist, const char *name, int role)
{
    /** Find a string name in a list of names and return the index. 
        Search should work with starting / characters and without.
        Create adios error and return -1 if name is null or
          if name is not found in the list.
        role = 0 for variable search, 1 for attribute search
     */
    int id, nstartpos=0, sstartpos;
    char ** s = namelist;
    char *rolename[2] = { "variable", "attribute" };
    int roleerror[2] = { err_invalid_varname, err_invalid_attrname };

    if (!name) {
        adios_error (roleerror[role!=0], "Null pointer passed as %s name!", rolename[role!=0]);
        return -1;
    }

    // find names with or without beginning /
    if (*name == '/') nstartpos = 1;

    for (id=0; id < n; id++) {
        if (*s[0] == '/') sstartpos = 1;
        else sstartpos = 0;
        //DBG_PRINTF("     check %s, startpos=%d\n", *s, sstartpos);
        if (!strcmp (*s+sstartpos, name+nstartpos))
            break; // found this name
        s++;
    }

    if (id == n) {
        adios_error (roleerror[role!=0], "%s '%s' is not found! One "
                "possible error is to set the view to a specific group and "
                "then try to read a %s of another group. In this case, "
                "reset the group view with adios_group_view(fp,-1).", 
                rolename[role!=0], name, rolename[role!=0]);
        return -1;
    }
    return id;
}


ADIOS_VARINFO * common_read_inq_var (const ADIOS_FILE *fp, const char * varname) 
{
    struct common_read_internals_struct * internals;
    ADIOS_VARINFO * retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        int varid = common_read_find_name (fp->nvars, fp->var_namelist, varname, 0);
        if (varid >= 0) {
            retval = common_read_inq_var_byid (fp, varid);
        } else {
            retval = NULL;
        }
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
        /* Translate varid to varid in global varlist if a selected group is in view */ 
        retval = internals->read_hooks[internals->method].adios_inq_var_byid_fn (fp, varid+internals->group_varid_offset);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var_byid()");
        retval = NULL;
    }
    return retval;
}


int common_read_inq_var_stat (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo,
                             int per_step_stat, int per_block_stat)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_inq_var_stat_fn (fp, varinfo, per_step_stat, per_block_stat);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var_stat()");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

int common_read_inq_var_blockinfo (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_inq_var_blockinfo_fn (fp, varinfo);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var_blockinfo()");
        retval = err_invalid_file_pointer;
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
        int varid = common_read_find_name (fp->nvars, fp->var_namelist, varname, 0);
        if (varid >= 0) {
            retval = common_read_schedule_read_byid (fp, sel, varid, from_steps, nsteps, data);
        } else {
            retval = adios_errno; // adios_errno was set in common_read_find_name
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_schedule_read()");
        retval = err_invalid_file_pointer;
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
        retval = internals->read_hooks[internals->method].adios_schedule_read_byid_fn (fp, sel, varid+internals->group_varid_offset, from_steps, nsteps, data);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_schedule_read_byid()");
        retval = err_invalid_file_pointer;
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
        retval = err_invalid_file_pointer;
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
        retval = err_invalid_file_pointer;
    }
    return retval;
}


void common_read_free_chunk (ADIOS_VARCHUNK *chunk)
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
        int attrid = common_read_find_name (fp->nattrs, fp->attr_namelist, attrname, 1);
        if (attrid > -1) {
            retval = common_read_get_attr_byid (fp, attrid, type, size, data);
        } else {
            retval = adios_errno; // adios_errno was set in common_read_find_name
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_read_get_attr()");
        retval = err_invalid_file_pointer;
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
        retval = internals->read_hooks[internals->method].adios_get_attr_byid_fn (fp, attrid+internals->group_attrid_offset, type, size, data);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_read_get_attr_byid()");
        retval = err_invalid_file_pointer;
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


int common_read_get_grouplist (const ADIOS_FILE  *fp, char ***group_namelist)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->ngroups;
        *group_namelist = internals->group_namelist;
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_get_grouplist()");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

/** Select a subset of variables and attributes to present in ADIOS_FILE struct.
    ADIOS_FILE-> nvars, nattrs, var_namelist, attr_namelist will contain
    only a subset of all variables and attributes.
    internals-> full_* stores the complete lists for reset or change of group
 */
int common_read_group_view (ADIOS_FILE  *fp, int groupid)
{
    struct common_read_internals_struct * internals;
    int retval, i;
    
    adios_errno = 0;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        if (groupid >= 0 && groupid < internals->ngroups) {
            /* 1. save complete list if first done */
            if (internals->group_in_view == -1) {
                internals->full_nvars = fp->nvars;
                internals->full_varnamelist = fp->var_namelist;
                internals->full_nattrs = fp->nattrs;
                internals->full_attrnamelist = fp->attr_namelist;
            }
            /* Set ID offsets for easier indexing of vars/attrs in other functions */
            internals->group_varid_offset = 0;
            internals->group_attrid_offset = 0;
            for (i=0; i<groupid; i++) {
                internals->group_varid_offset += internals->nvars_per_group[i];
                internals->group_attrid_offset += internals->nattrs_per_group[i];
            }
            /* Set view to this group */
            fp->nvars = internals->nvars_per_group[groupid];
            fp->var_namelist = &(internals->full_varnamelist [internals->group_varid_offset]);
            fp->nattrs = internals->nattrs_per_group[groupid];
            fp->attr_namelist = &(internals->full_attrnamelist [internals->group_attrid_offset]);
            internals->group_in_view = groupid;
            retval = 0;

        } else if (groupid == -1) {
            /* Reset to full view */
            fp->nvars  = internals->full_nvars;
            fp->var_namelist  = internals->full_varnamelist;
            fp->nattrs = internals->full_nattrs;
            fp->attr_namelist  = internals->full_attrnamelist;
            internals->group_varid_offset = 0;
            internals->group_attrid_offset = 0;
            internals->group_in_view = -1;
            retval = 0;
        } else {
            adios_error (err_invalid_group, "Invalid group ID in adios_group_view()");
            retval = err_invalid_group;
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_group_view()");
        retval = err_invalid_file_pointer;
    }
    return retval;
}


void common_read_print_fileinfo (const ADIOS_FILE *fp) 
{
    int i;
    int ngroups;
    char **group_namelist;
    ngroups = common_read_get_grouplist (fp, &group_namelist);

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


/**    SELECTIONS   **/ 
ADIOS_SELECTION * common_read_selection_boundingbox (uint64_t ndim, const uint64_t *start, const uint64_t *count)
{   
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_BOUNDINGBOX;
        sel->u.bb.ndim = ndim;
        sel->u.bb.start = (uint64_t *)start;
        sel->u.bb.count = (uint64_t *)count;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for bounding box selection");
    }
    return sel;
}


ADIOS_SELECTION * common_read_selection_points (uint64_t ndim, uint64_t npoints, const uint64_t *points)
{   
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_POINTS;
        sel->u.points.ndim = ndim;
        sel->u.points.npoints = npoints;
        sel->u.points.points = (uint64_t *) points;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for points selection");
    }
    return sel;
}

ADIOS_SELECTION * common_read_selection_writeblock (int index)
{   
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_WRITEBLOCK;
        sel->u.block.index = index;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for writeblock selection");
    }
    return sel;
}

ADIOS_SELECTION * common_read_selection_auto (char *hints)
{   
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_AUTO;
        sel->u.autosel.hints = hints;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for auto selection");
    }
    return sel;
}

void common_read_selection_delete (ADIOS_SELECTION *sel)
{
    free(sel);
}
