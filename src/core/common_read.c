/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <inttypes.h>
#include "public/adios_error.h"
#include "core/adios_logger.h"
#include "core/common_read.h"
#include "core/futils.h"
#include "core/bp_utils.h" // struct namelists_struct
#include "public/adios_schema.h"

// NCSU ALACRITY-ADIOS
#include "adios_read_hooks.h"
#include "transforms/adios_transforms_transinfo.h"
#include "transforms/adios_transforms_hooks_read.h"
#include "transforms/adios_transforms_reqgroup.h"
#include "transforms/adios_transforms_datablock.h"
#define BYTE_ALIGN 8

#ifdef DMALLOC
#include "dmalloc.h"
#endif


/* Note: MATLAB reloads the mex64 files each time, so all static variables get the original value.
   Therefore static variables cannot be used to pass info between two Matlab/ADIOS calls */
static struct adios_read_hooks_struct * adios_read_hooks = 0;

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

    // NCSU ALACRITY-ADIOS - Table of sub-requests issued by transform method
    adios_transform_read_request *transform_reqgroups;
};

// NCSU ALACRITY-ADIOS - Forward declaration/function prototypes
static void common_read_free_blockinfo(ADIOS_VARBLOCK **varblock, int sum_nblocks);



int common_read_init_method (enum ADIOS_READ_METHOD method,
                             MPI_Comm comm,
                             const char * parameters)
{
    PairStruct *params, *p, *prev_p;
    int verbose_level, removeit, save;
    int retval; 
    char *end;

    adios_errno = err_no_error;
    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, 
            "Invalid read method (=%d) passed to adios_read_init_method().\n", (int)method);
        return err_invalid_read_method;
    } 
    // init the adios_read_hooks_struct if not yet initialized  
    adios_read_hooks_init (&adios_read_hooks); 
    // NCSU ALACRITY-ADIOS - Initialize transform methods
    adios_transform_read_init();

    // process common parameters here
    params = text_to_name_value_pairs (parameters);
    p = params;
    prev_p = NULL;
    while (p) {
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
                p = p->next;
                params->next = NULL;
                free_name_value_pairs (params);
                params = p;
            } else {
                // remove from middle of the list
                prev_p->next = p->next;
                p->next = NULL;
                free_name_value_pairs (p);
                p = prev_p->next;
            }
        } else {
            prev_p = p;
            p = p->next;
        }
    }

    // call method specific init 
    retval = adios_read_hooks[method].adios_init_method_fn (comm, params);
    free_name_value_pairs (params);
    return retval;
}


int common_read_finalize_method(enum ADIOS_READ_METHOD method)
{
    adios_errno = err_no_error;
    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, 
            "Invalid read method (=%d) passed to adios_read_finalize_method().\n", (int)method);
        return err_invalid_read_method;
    } 

    return adios_read_hooks[method].adios_finalize_method_fn ();
}


ADIOS_FILE * common_read_open (const char * fname, 
                               enum ADIOS_READ_METHOD method, 
                               MPI_Comm comm, 
                               enum ADIOS_LOCKMODE lock_mode, 
                               float timeout_sec)
{
    ADIOS_FILE * fp;
    struct common_read_internals_struct * internals; 

    if ((int)method < 0 || (int)method >= ADIOS_READ_METHOD_COUNT) {
        adios_error (err_invalid_read_method, 
            "Invalid read method (=%d) passed to adios_read_open().\n", (int)method);
        return NULL;
    } 

    adios_errno = err_no_error;
    internals = (struct common_read_internals_struct *) 
                    calloc(1,sizeof(struct common_read_internals_struct));
    // init the adios_read_hooks_struct if not yet initialized 
    adios_read_hooks_init (&adios_read_hooks); 
    // NCSU ALACRITY-ADIOS - Initialize transform methods
    adios_transform_read_init();

    internals->method = method;
    internals->read_hooks = adios_read_hooks;

    fp = adios_read_hooks[internals->method].adios_open_fn (fname, comm, lock_mode, timeout_sec);

    //read mesh names from attributes for example the var is using a mesh named trimesh, 
    //we have /adios_schema/trimesh/type. We can extract trimesh from the string
    fp->nmeshes = 0;
    fp->mesh_namelist = NULL;
    int i;
    if (fp->attr_namelist)
    {
        char ** tmp = (char **) malloc (sizeof(char*) * fp->nattrs);
        for (i=0; i<fp->nattrs; i++)
        {
            // find "/adios_schema/***/type" attributes for getting the names of meshes
            if (strstr (fp->attr_namelist[i], "/adios_schema/") == fp->attr_namelist[i])   // starts with /adios_schema/
            {
                char *s = fp->attr_namelist[i]+strlen("/adios_schema/");
                char *p = strchr (s, '/');  
                if ( p && 
                     strstr (p, "/type") == p)
                {
                    // retrieve the name of the mesh
                    tmp [ fp->nmeshes ] = (char *) malloc (sizeof(char*) * (size_t)(p-s)+1);
                    memcpy ( tmp[ fp->nmeshes ], s, (size_t)(p-s) ); 
                    tmp[ fp->nmeshes ][(p-s)] = '\0'; 
                    fp->nmeshes++;
                }
            }
        }

        if (fp->nmeshes) 
        {
            fp->mesh_namelist = (char **) realloc (tmp, sizeof (char *) * fp->nmeshes);
            assert (fp->mesh_namelist);
        }

/*
        int c = 0;
        for (i=0; i<fp->nattrs; i++)
        {
            if (strstr (fp->attr_namelist[i], "/adios_schema") && strstr (fp->attr_namelist[i], "/type") &&
                strcspn(fp->attr_namelist[i], "/adios_schema") == 0 && strcspn (fp->attr_namelist[i], "/type") == strlen(fp->attr_namelist[i]-strlen("/type")) )
            {
                fp->mesh_namelist[c] = (char *) malloc (strlen (fp->attr_namelist[i]) + 1);
                assert (fp->mesh_namelist[c]);
                //avoid to copy slash
                strncpy (fp->mesh_namelist[c], fp->attr_namelist[i]+strlen("/adios_schema")+1, strlen(fp->attr_namelist[i])-strlen("/adios_schema")-strlen("/type")-1); 
                c++;
            }
        }
*/

    }

/*        for (i=0; i<fp->nmeshes; i++)
        {
            printf ("mesh: %s\n", fp->mesh_namelist[i]);
        }*/
    
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
            "Invalid read method (=%d) passed to adios_read_open_file().\n", (int)method);
        return NULL;
    } 

    adios_errno = err_no_error;
    internals = (struct common_read_internals_struct *) 
                    calloc(1,sizeof(struct common_read_internals_struct));
    // init the adios_read_hooks_struct if not yet initialized 
    adios_read_hooks_init (&adios_read_hooks); 
    // NCSU ALACRITY-ADIOS - Initialize transform methods
    adios_transform_read_init();

    internals->method = method;
    internals->read_hooks = adios_read_hooks;

    fp = adios_read_hooks[internals->method].adios_open_file_fn (fname, comm);
    
    //read mesh names from attributes for example the var is using a mesh named trimesh, 
    //we have /adios_schema/trimesh/type. We can extract trimesh from the string
    fp->nmeshes = 0;
    fp->mesh_namelist = NULL;
    int i;

    if (fp->attr_namelist)
    {
        char ** tmp = (char **) malloc (sizeof(char*) * fp->nattrs);
        for (i=0; i<fp->nattrs; i++)
        {
            // find "/adios_schema/***/type" attributes for getting the names of meshes
            if (strstr (fp->attr_namelist[i], "/adios_schema/") == fp->attr_namelist[i])   // starts with /adios_schema/
            {
                char *s = fp->attr_namelist[i]+strlen("/adios_schema/");
                char *p = strchr (s, '/');
                if ( p &&
                     strstr (p, "/type") == p)
                {
                    // retrieve the name of the mesh
                    tmp [ fp->nmeshes ] = (char *) malloc (sizeof(char*) * (size_t)(p-s)+1);
                    memcpy ( tmp[ fp->nmeshes ], s, (size_t)(p-s) );
                    tmp[ fp->nmeshes ][(p-s)] = '\0';
                    fp->nmeshes++;
                }
            }
        }

        if (fp->nmeshes)
        {
            fp->mesh_namelist = (char **) realloc (tmp, sizeof (char *) * fp->nmeshes);
            assert (fp->mesh_namelist);
        }

    }

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

// NCSU ALACRITY-ADIOS - Cleanup for read request groups
#define MYFREE(p) {if (p) free((void*)p); (p)=NULL;}
static void clean_up_read_reqgroups(adios_transform_read_request **reqgroups_head) {
    adios_transform_read_request *removed;
    while ((removed = adios_transform_read_request_pop(reqgroups_head)) != NULL) {
        adios_transform_read_request_free(&removed);
    }
}
#undef MYFREE

int common_read_close (ADIOS_FILE *fp) 
{
    struct common_read_internals_struct * internals;
    int retval;

    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        if (internals->group_in_view != -1) {
            // reset from group view before calling the real close
            common_read_group_view (fp, -1);
        }
        int i;
        if (fp->nmeshes) {
            for (i=0; i<fp->nmeshes; i++)
                free(fp->mesh_namelist[i]);
            free(fp->mesh_namelist);
        }
                
        retval = internals->read_hooks[internals->method].adios_close_fn (fp);
        free_namelist (internals->group_namelist, internals->ngroups);
        free (internals->nvars_per_group);
        free (internals->nattrs_per_group);
        // NCSU ALACRITY-ADIOS - Cleanup read request groups
        clean_up_read_reqgroups(&internals->transform_reqgroups);
        free (internals);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_read_close()\n");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

void common_read_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{
    struct common_read_internals_struct * internals;

    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        internals->read_hooks[internals->method].adios_reset_dimension_order_fn (fp, is_fortran);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_reset_dimension_order()\n");
    }
}


int common_read_advance_step (ADIOS_FILE *fp, int last, float timeout_sec)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_advance_step_fn (fp, last, timeout_sec);
        if (!retval) {
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
        }
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_advance_step()\n");
        retval = err_invalid_file_pointer;
    }
    return retval;
}


void common_read_release_step (ADIOS_FILE *fp) 
{
    struct common_read_internals_struct * internals;

    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        internals->read_hooks[internals->method].adios_release_step_fn (fp);
    } else {
        adios_error ( err_invalid_file_pointer, "Invalid file pointer at adios_reset_dimension_order()\n");
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
    enum ADIOS_ERRCODES roleerror[2] = { err_invalid_varname, err_invalid_attrname };

    if (!name) {
        adios_error (roleerror[role!=0], "Null pointer passed as %s name!\n", rolename[role!=0]);
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
                "reset the group view with adios_group_view(fp,-1).\n",
                rolename[role!=0], name, rolename[role!=0]);
        return -1;
    }
    return id;
}


ADIOS_VARINFO * common_read_inq_var (const ADIOS_FILE *fp, const char * varname) 
{
    struct common_read_internals_struct * internals;
    ADIOS_VARINFO * retval;
 
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        int varid = common_read_find_name (fp->nvars, fp->var_namelist, varname, 0);
        if (varid >= 0) {
            retval = common_read_inq_var_byid (fp, varid);
        } else {
            retval = NULL;
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var()\n");
        retval = NULL;
    }
    return retval;
}

// NCSU ALACRITY-ADIOS - For copying original metadata from transform
//   info to inq var info
static void patch_varinfo_with_transform_blockinfo(ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    common_read_free_blockinfo(&vi->blockinfo, vi->sum_nblocks);    // Free blockinfo in varinfo
    vi->blockinfo = ti->orig_blockinfo;                                // Move blockinfo from transinfo to varinfo
    ti->orig_blockinfo = 0;                                            // Delink blockinfo from transinfo
}
static void patch_varinfo_with_transinfo(ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    // First make room for the transform info fields
    free(vi->dims);

    // Now move them
    vi->type = ti->orig_type;
    vi->ndim = ti->orig_ndim;
    vi->global = ti->orig_global;
    vi->dims = ti->orig_dims;

    // Finally, delink them from the transform info so they aren't inadvertently free'd
    ti->orig_dims = 0;

    patch_varinfo_with_transform_blockinfo(vi, ti); // Also move blockinfo if extant
}

// NCSU ALACRITY-ADIOS - Delegate to the 'inq_var_raw_byid' function, then
//   patch the original metadata in from the transform info
ADIOS_VARINFO * common_read_inq_var_byid (const ADIOS_FILE *fp, int varid)
{
    ADIOS_VARINFO *vi;
    ADIOS_TRANSINFO *ti;

    vi = common_read_inq_var_raw_byid(fp, varid);
    if (vi == NULL)
        return NULL;
    
    vi->meshinfo = NULL;

    // NCSU ALACRITY-ADIOS - translate between original and transformed metadata if necessary
    ti = common_read_inq_transinfo(fp, vi); // No orig_blockinfo
    if (ti && ti->transform_type != adios_transform_none) {
        patch_varinfo_with_transinfo(vi, ti);
    }
    common_read_free_transinfo(vi, ti);

    return vi;
}

// NCSU ALACRITY-ADIOS - Renaming of common_read_inq_var_byid, named 'raw'
//   because it is oblivious to the original metadata as stored in TRANSINFO
ADIOS_VARINFO * common_read_inq_var_raw_byid (const ADIOS_FILE *fp, int varid)
{
    struct common_read_internals_struct * internals;
    ADIOS_VARINFO * retval;
    
    adios_errno = err_no_error;
    if (fp) {
        if (varid >= 0 && varid < fp->nvars) {
            internals = (struct common_read_internals_struct *) fp->internal_data;
            /* Translate varid to varid in global varlist if a selected group is in view */ 
            retval = internals->read_hooks[internals->method].adios_inq_var_byid_fn 
                                            (fp, varid+internals->group_varid_offset);
            if (retval) {
                /* Translate real varid to the group varid presented to the user */
                retval->varid = varid;
            }
        } else {
            adios_error (err_invalid_varid, 
                         "Variable ID %d is not valid adios_inq_var_byid(). "
                         "Available 0..%d\n", varid, fp->nvars-1);
            retval = NULL;
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var_byid()\n");
        retval = NULL;
    }
    return retval;
}

// NCSU ALACRITY-ADIOS - common-layer inquiry function to read transform info.
// NOTE: does not follow the normal pattern of adding the info into
//   ADIOS_VARINFO because this information should not be sent to the user;
//   only in rare cases will a user application need this information (like a
//   query engine using an transform-embedded index), in which case that code
//   can dive deeper and access this function. Alternatively, if this use case
//   becomes more common, a simple 'transform raw' API could be added.
ADIOS_TRANSINFO * common_read_inq_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi) {
    if (!fp) {
        adios_error (err_invalid_file_pointer,
                     "Null ADIOS_FILE pointer passed to common_read_inq_transinfo()\n");
        return NULL;
    }
    if (!vi) {
        adios_error (err_invalid_argument,
                     "Null ADIOS_VARINFO pointer passed to common_read_inq_transinfo()\n");
        return NULL;
    }

    struct common_read_internals_struct * internals;
    internals = (struct common_read_internals_struct *) fp->internal_data;

    ADIOS_TRANSINFO *ti = internals->read_hooks[internals->method].adios_inq_var_transinfo_fn(fp, vi);
    return ti;
}

int common_read_inq_trans_blockinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO * ti) {
    if (!fp) {
        adios_error (err_invalid_argument,
                     "Null ADIOS_FILE pointer passed to common_read_inq_trans_blockinfo()\n");
        return 1;
    }
    if (!vi) {
        adios_error (err_invalid_argument,
                     "Null ADIOS_VARINFO pointer passed to common_read_inq_trans_blockinfo()\n");
        return 1;
    }
    if (!ti) {
        adios_error (err_invalid_argument,
                     "Null ADIOS_TRANSINFO pointer passed to common_read_inq_trans_blockinfo()\n");
        return 1;
    }

    struct common_read_internals_struct * internals;
    internals = (struct common_read_internals_struct *) fp->internal_data;
    return internals->read_hooks[internals->method].adios_inq_var_trans_blockinfo_fn(fp, vi, ti);
}



int common_read_inq_var_stat (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo,
                             int per_step_stat, int per_block_stat)
{
    struct common_read_internals_struct * internals;
    int retval;
    int group_varid;
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        if (varinfo) {
            /* Translate group varid presented to the user to the real varid */
            group_varid = varinfo->varid;
            varinfo->varid = varinfo->varid + internals->group_varid_offset;
        }
        retval = internals->read_hooks[internals->method].adios_inq_var_stat_fn (fp, varinfo, per_step_stat, per_block_stat);
        /* Translate back real varid to the group varid presented to the user */
        varinfo->varid = group_varid;
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var_stat()\n");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

// NCSU ALACRITY-ADIOS - Delegate to the 'inq_var_blockinfo_raw' function, then
//   patch the original metadata in from the transform info
int common_read_inq_var_blockinfo (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    ADIOS_TRANSINFO *ti;

    int retval = common_read_inq_var_blockinfo_raw(fp, varinfo);
    if (retval != err_no_error)
        return retval;

    // NCSU ALACRITY-ADIOS - translate between original and transformed metadata if necessary
    ti = common_read_inq_transinfo(fp, varinfo);
    if (ti && ti->transform_type != adios_transform_none) {
        retval = common_read_inq_trans_blockinfo(fp, varinfo, ti);
        if (retval != err_no_error)
            return retval;

        patch_varinfo_with_transform_blockinfo(varinfo, ti);
    }
    common_read_free_transinfo(varinfo, ti);

    return err_no_error;
}

// NCSU ALACRITY-ADIOS - Renaming of common_read_inq_var_blockinfo, named 'raw'
//   because it is oblivious to the original metadata as stored in TRANSINFO
int common_read_inq_var_blockinfo_raw (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    struct common_read_internals_struct * internals;
    int retval;
    int group_varid;
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        if (varinfo) {
            /* Translate group varid presented to the user to the real varid */
            group_varid = varinfo->varid;
            varinfo->varid = varinfo->varid + internals->group_varid_offset;
        }
        retval = internals->read_hooks[internals->method].adios_inq_var_blockinfo_fn (fp, varinfo);
        /* Translate back real varid to the group varid presented to the user */
        varinfo->varid = group_varid;
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_inq_var_blockinfo()\n");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

#define MYFREE(p) {free(p); (p)=NULL;}
// NCSU ALACRITY-ADIOS - Factored this out to use elsewhere
static void common_read_free_blockinfo(ADIOS_VARBLOCK **varblock, int sum_nblocks) {
    if (*varblock) {
    int i;
        ADIOS_VARBLOCK *bp = *varblock;
        for (i = 0; i < sum_nblocks; i++) {
                if (bp->start) MYFREE (bp->start);
                if (bp->count) MYFREE (bp->count);
                bp++;
            }
        MYFREE(*varblock);
        }
}

void common_read_free_varinfo (ADIOS_VARINFO *vp)
{
    int i;
    if (vp) {
        common_read_free_blockinfo(&vp->blockinfo, vp->sum_nblocks);

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
        if (vp->meshinfo) MYFREE(vp->meshinfo);
        free(vp);
    }
}

// NCSU ALACRITY-ADIOS - Free transform info
void common_read_free_transinfo(const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    if (ti) {
        if (ti->orig_dims) MYFREE(ti->orig_dims);
        if (ti->transform_metadata && ti->should_free_transform_metadata)
            MYFREE(ti->transform_metadata);

        common_read_free_blockinfo(&ti->orig_blockinfo, vi->sum_nblocks);

        free(ti);
    }
}
#undef MYFREE

// the function is the same as common_read_find_name
// but no ERROR msg print out
// called by common_read_get_attr_mesh
static int common_read_find_name_mesh (int n, char ** namelist, const char *name, int role)
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
    enum ADIOS_ERRCODES roleerror[2] = { err_invalid_varname, err_invalid_attrname };

    if (!name) {
        adios_error (roleerror[role!=0], "Null pointer passed as %s name!\n", rolename[role!=0]);
        return -1;
    }

    // find names with or without beginning /
    if (*name == '/') nstartpos = 1;

    for (id=0; id < n; id++) {
        if (*s[0] == '/') sstartpos = 1;
        else sstartpos = 0;
        //DBG_PRINTF(" check %s, startpos=%d\n", *s, sstartpos);
        if (!strcmp (*s+sstartpos, name+nstartpos))
            break; // found this name
        s++;
    }
    return id;
}

// the function is the same as common_read_get_attr_byid
// but no ERROR msg print out
// called by common_read_get_attr_mesh
int common_read_get_attr_byid_mesh (const ADIOS_FILE * fp,
                               int attrid,
                               enum ADIOS_DATATYPES * type,
                               int * size,
                               void ** data)
{
    struct common_read_internals_struct * internals;
    int retval;

    adios_errno = err_no_error;
    if (fp) {
        if (attrid >= 0 && attrid < fp->nattrs) {
            internals = (struct common_read_internals_struct *) fp->internal_data;
            retval = internals->read_hooks[internals->method].adios_get_attr_byid_fn (fp, attrid+internals->group_attrid_offset, type, size, data);
        } else {
            retval = err_invalid_attrid;
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_read_get_attr_byid()\n");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

// this function is almost the same as common_read_get_attr
// just to avoid the ERROR msg when some attributes are not found
// for example spacing/maximum are optinal in uniform mesh
int common_read_get_attr_mesh (const ADIOS_FILE * fp,
                            const char * attrname,
                            enum ADIOS_DATATYPES * type,
                            int * size,
                            void ** data)
{
    struct common_read_internals_struct * internals;
    int retval;

    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        int attrid = common_read_find_name_mesh (fp->nattrs, fp->attr_namelist, attrname, 1);
        if (attrid > -1) {
            retval = common_read_get_attr_byid_mesh (fp, attrid, type, size, data);
        } else {
            retval = adios_errno; // adios_errno was set in common_read_find_name
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_read_get_attr()\n");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

int common_read_inq_var_meshinfo (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    enum ADIOS_DATATYPES attr_type;
    int  attr_size;
    int  read_fail = 0;
    void * data = NULL;
    int i;
    int match;

    varinfo->meshinfo = (ADIOS_VARMESH *) malloc (sizeof(ADIOS_VARMESH));
    char * var_name = strdup (fp->var_namelist[varinfo->varid]);
//    printf ("var name is %s\n", var_name);
    char * var_mesh = malloc (strlen(var_name)+strlen("/adios_schema")+1);
    strcpy (var_mesh, var_name);
    strcat (var_mesh, "/adios_schema");
    
    read_fail = common_read_get_attr_mesh (fp, var_mesh, &attr_type, &attr_size, &data); 
    if (read_fail)
    {
        printf ("ERROR: no matching mesh for var %s\n", var_name);
        return 1;
    }
    else
    {
        match = 0;
//        printf ("meshname from attr is %s\n", (char *)data);
        for (i=0; i<fp->nmeshes; i++)
        {
//            printf ("mesh name is %s\n", fp->mesh_namelist[i]);
            if ( !strcmp(fp->mesh_namelist[i], (char *)data))
            {
                match = 1;
                varinfo->meshinfo->meshid = i;
            }
        }
        if (match == 0)
        {
            printf ("ERROR: mesh %s for var %s is not stored in meshlist\n", (char *)data, var_name);
            return 1;
        }
    }

    // point centering or cell centering
    char * data_centering = malloc (strlen(var_mesh)+strlen("/centering")+1);
    strcpy (data_centering, var_mesh);
    strcat (data_centering, "/centering");
    read_fail = common_read_get_attr_mesh (fp, data_centering, &attr_type, &attr_size, &data); 
//    printf ("attr data_centering is %s\n", data_centering);
    free (data_centering);
    free (var_mesh);
    if (read_fail)        // if no attr for centering, check if it is unstrcutured
    {
        varinfo->meshinfo->centering = 0;          // no centering
        char * meshtype = malloc (strlen("/adios_schema/")+strlen(var_name)+strlen("/type")+1);
        strcpy (meshtype, "/adios_schema/");
        strcat (meshtype, fp->mesh_namelist[varinfo->meshinfo->meshid]);      
        strcat (meshtype, "/type");
//        printf ("attr meshtype is %s\n", meshtype);
        data = NULL;
        read_fail = common_read_get_attr_mesh (fp, meshtype, &attr_type, &attr_size, &data);
        if (read_fail)
        {
            printf ("ERROR: cannot get mesh name from attr %s is not available\n", meshtype);
            free (meshtype);
            return 1;
        }
        else
        {
            free (meshtype);
            if (!strcmp((char *)data, "unstructured"))
            {
                printf ("ERROR: centering info of var %s on unstructured mesh %s is required\n", var_name, fp->mesh_namelist[varinfo->meshinfo->meshid]);
                return 1;
            }
        }
    }
    else
    {
        if (!strcmp((char *)data, "point"))
        {
            varinfo->meshinfo->centering = 1;      // point centering
        }
        else if (!strcmp((char *)data, "cell"))
        {
            varinfo->meshinfo->centering = 2;      // cell centering
        }
        else
        {
            printf ("ERROR: centering method of var %s on mesh %s is not supported (point/cell)\n", var_name, fp->mesh_namelist[varinfo->meshinfo->meshid]);
            return 1;
        }
    }

    return 0;
}


int adios_get_uniform_mesh_attr (ADIOS_FILE * fp, ADIOS_MESH *meshinfo, char * attrs)      //attr for origins-num(origins), spacings-num(spacings), maximums-num(maximums)
{
    int i, j;
    enum ADIOS_DATATYPES attr_type;
    int attr_size;
    void * data = NULL;
    int  read_fail = 0;
    bool have_max = 0;
    bool have_spacing = 0;

    char * attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/")+strlen(attrs)+strlen("-num")+1 );
    strcpy (attribute, "/adios_schema/");
    strcat (attribute, meshinfo->name);
    strcat (attribute, "/");
    strcat (attribute, attrs);
    strcat (attribute, "-num");
//    printf("attribute is %s\n", attribute);
    data = NULL;
    read_fail = common_read_get_attr_mesh (fp, attribute, &attr_type, &attr_size, &data);
    free (attribute);
    if (!read_fail)                       //found attributes maximums/spacings/origins
    {
        int num_attr = *(int *)data;
        if (num_attr != meshinfo->uniform->num_dimensions)
        {
            if (!strcmp (attrs,"origins"))
            {
                num_attr = meshinfo->uniform->num_dimensions;
                printf("WARNING: mesh %s number of origins %d does not match number of dimensions %d! Use number of dimensions for origins\n", 
                        meshinfo->name, num_attr, meshinfo->uniform->num_dimensions);
            }
            else if (!strcmp (attrs, "spacings"))
            {
                num_attr = meshinfo->uniform->num_dimensions;
                printf("WARNING: mesh %s number of origins %d does not match number of dimensions %d! Use number of dimensions for spacings\n", 
                        meshinfo->name, num_attr, meshinfo->uniform->num_dimensions);
            }
            else if (!strcmp (attrs, "maximums"))
            {
                if (num_attr < meshinfo->uniform->num_dimensions)
                {
                    printf("ERROR: mesh %s number of maximums %d is less than the number of dimensions %d!\n", 
                            meshinfo->name, num_attr, meshinfo->uniform->num_dimensions);
                    meshinfo->uniform = NULL;
                    return -1; 
                }
                else if (num_attr > meshinfo->uniform->num_dimensions)
                {
                    num_attr = meshinfo->uniform->num_dimensions;
                    printf("WARNING: mesh %s provided number of maximums %d is greater than the number of dimensions %d! Use number of dimentions \
                           for maximums and ignore the rest maximums\n", meshinfo->name, num_attr, meshinfo->uniform->num_dimensions);
                }
            }
        }
        if (!strcmp (attrs,"origins"))
        {
            meshinfo->uniform->origins = (double *) malloc (sizeof(double)*num_attr);
            for (i = 0; i < num_attr; i++ )
                meshinfo->uniform->origins[i] = 0.0;
        }
        else if (!strcmp (attrs, "spacings"))
        {
            have_spacing = 1;
            meshinfo->uniform->spacings = (double *) malloc (sizeof(double)*num_attr);
            for (i = 0; i < num_attr; i++ )
            {
                meshinfo->uniform->spacings[i] = 1.0;
//                printf("meshinfo->uniform->spacings[%d] = %lf\n", i, meshinfo->uniform->spacings[i]);
            }
        }
        else if (!strcmp (attrs, "maximums"))
        {
            have_max = 1;
            meshinfo->uniform->maximums = (double *) malloc (sizeof(double)*num_attr);
            for (i = 0; i < num_attr; i++ )
                meshinfo->uniform->maximums[i] = 0.0;
        }

        for (i = 0; i < num_attr; i++ )
        {
            char * i_buffer;
            int i_digits;
            if (num_attr < 10)
                i_buffer = (char *) malloc (sizeof(char)+1);
            else
            {
                printf("ERROR: mesh %s has more than 10 dimensions!\n", meshinfo->name);
                meshinfo->uniform = NULL;
                return 0; 
            }
            i_digits = sprintf (i_buffer, "%d", i);
            //i_digits to reprent the number in dimensions0/dimensions1/... 
            char * value = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/")+strlen(attrs)-strlen("-num")+i_digits+1 );
            strcpy (value, "/adios_schema/");
            strcat (value, meshinfo->name);
            strcat (value, "/");
            strcat (value, attrs);
            strcat (value, i_buffer);
//            printf("origin_value is %s\n", value);
            free (i_buffer);
            data = NULL;
            read_fail = common_read_get_attr_mesh (fp, value, &attr_type, &attr_size, &data);
            free (value);
            if (read_fail)
            {
                if (!strcmp (attrs, "origins"))
                    printf("WARNING: mesh %s origins[%d] value is not found, use default value 0\n", meshinfo->name, i);   // origins are set to 0 at the beginning of origin processing
                else if (!strcmp (attrs, "maximums"))
                {
                    printf ("ERROR: mesh %s maximum of maximums[%d] is not provided!\n", meshinfo->name, i);
                    meshinfo->uniform = NULL;
                    return 0; 
                }
                else if (!strcmp (attrs, "spacings"))
                    printf("WARNING: mesh %s spacings[%d] value is not found, use default value 1\n", meshinfo->name, i);
            }
            else
            {
                char * pEnd;
                double d1;
                char * tmp_dimensions_value = strdup((char *)data);
                d1 = strtod (tmp_dimensions_value, &pEnd);
                if (d1 != 0.0)
                {
                    if (!strcmp (attrs, "origins"))
                        meshinfo->uniform->origins[i] = d1;                           //
                    else if (!strcmp (attrs, "maximums"))
                        meshinfo->uniform->maximums[i] = d1;
                    else if (!strcmp (attrs, "spacings"))
                        meshinfo->uniform->spacings[i] = d1;
                }
                else
                {
                    int var_march = 0;
                    for (j=0; j<fp->nvars; j++)
                    {
                        char * value_tmp = NULL;
                        value_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                        strcpy (value_tmp, "/");
                        strcat (value_tmp, (char *)data);
                        if (!strcmp (fp->var_namelist[j], value_tmp))
                        {
                            var_march = 1;
                            ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                            if (!strcmp (attrs, "origins"))
                            {
                                if (v->type == adios_real)    //adios_real
                                    meshinfo->uniform->origins[i] = *(float *)v->value;                        //
                                else if (v->type == adios_double)    //adios_double
                                    meshinfo->uniform->origins[i] = *(double *)v->value;
                                else if (v->type == adios_integer || v->type == adios_unsigned_integer)  // adios_integer or adios_unsigned_integer
                                    meshinfo->uniform->origins[i] = *(int *)v->value;
                                else if (v->type == adios_long || v->type == adios_unsigned_long)  // adios_long or adios_unsigned_long
                                    meshinfo->uniform->origins[i] = *(long int *)v->value;
                                else if (v->type == adios_short || v->type == adios_unsigned_short)  // adios_short or adios_unsigned_short
                                    meshinfo->uniform->origins[i] = *(short int *)v->value;
                                else
                                {
                                    meshinfo->uniform = NULL;
                                    printf("ERROR: uniform mesh %s originss support (unsigned) short int, (unsigned) int, (unsigned) long int, float and double\n");
                                    return -1;
                                }
                            }
                            else if (!strcmp (attrs, "maximums"))
                            {
                                if (v->type == adios_real)    //adios_real
                                    meshinfo->uniform->maximums[i] = *(float *)v->value;
                                else if (v->type == adios_double)    //adios_double
                                    meshinfo->uniform->maximums[i] = *(double *)v->value;
                                else if (v->type == adios_integer || v->type == adios_unsigned_integer)  // adios_integer or adios_unsigned_integer
                                    meshinfo->uniform->maximums[i] = *(int *)v->value;
                                else if (v->type == adios_long || v->type == adios_unsigned_long)  // adios_long or adios_unsigned_long
                                    meshinfo->uniform->maximums[i] = *(long int *)v->value;
                                else if (v->type == adios_short || v->type == adios_unsigned_short)  // adios_short or adios_unsigned_short
                                    meshinfo->uniform->maximums[i] = *(short int *)v->value;
                                else
                                {
                                    meshinfo->uniform = NULL;
                                    printf("ERROR: uniform mesh %s maximums support (unsigned) short int, (unsigned) int, (unsigned) long int, float and double\n");
                                    return -1;
                                }
                            }
                            else if (!strcmp (attrs, "spacings"))
                            {
                                if (v->type == adios_real)    //adios_real
                                    meshinfo->uniform->spacings[i] = *(float *)v->value;
                                else if (v->type == adios_double)    //adios_double
                                    meshinfo->uniform->spacings[i] = *(double *)v->value;
                                else if (v->type == adios_integer || v->type == adios_unsigned_integer)  // adios_integer or adios_unsigned_integer
                                    meshinfo->uniform->spacings[i] = *(int *)v->value;
                                else if (v->type == adios_long || v->type == adios_unsigned_long)  // adios_long or adios_unsigned_long
                                    meshinfo->uniform->spacings[i] = *(long int *)v->value;
                                else if (v->type == adios_short || v->type == adios_unsigned_short)  // adios_short or adios_unsigned_short
                                    meshinfo->uniform->spacings[i] = *(short int *)v->value;
                                else 
                                {
                                    meshinfo->uniform = NULL;
                                    printf("ERROR: uniform mesh %s spacings support (unsigned) short int, (unsigned) int, (unsigned) long int, float and double\n");
                                    return -1;
                                }
//                                printf ("meshinfo->uniform->spacings[%d] = %lf\n",i, meshinfo->uniform->spacings[i]);
                            }
                            common_read_free_varinfo (v);
                        }
                        free (value_tmp);
                    }
                    if (!var_march)
                    {
                        if (!strcmp (attrs, "origins"))
                            printf ("WARNING: mesh %s origins%d is set to 0\n", meshinfo->name, i);
                        else if (!strcmp (attrs, "maximums"))
                        {
                            printf ("ERROR: mesh %s maximums%d var %s is not provided!\n", meshinfo->name, i, (char *)data);
                            meshinfo->uniform = NULL;
                            return 0; 
                        }
                        else if (!strcmp (attrs, "spacings"))
                            printf ("WARNING: mesh %s spacings%d var %s is not provided, set this sapcing value to 1\n", meshinfo->name, i, (char *)data);
                    }
                }
            }
        }
    }
    else                //not found attributes maximums/spacings/origins
    {
        if (!strcmp (attrs,"origins"))
        {
            meshinfo->uniform->origins = (double *) malloc (sizeof(double)*meshinfo->uniform->num_dimensions);
            for (i = 0; i < meshinfo->uniform->num_dimensions; i++ )
                meshinfo->uniform->origins[i] = 0.0;
            printf("WARNING: uniform mesh %s does not provide origin info. We use default value 0\n", meshinfo->name);
        }
        else if (!strcmp (attrs, "spacings"))
        {
            have_spacing = 1;
            meshinfo->uniform->spacings = (double *) malloc (sizeof(double)*meshinfo->uniform->num_dimensions);
            for (i = 0; i < meshinfo->uniform->num_dimensions; i++ )
                meshinfo->uniform->spacings[i] = 1.0;
            printf("WARNING: uniform mesh %s does not provide spacing info. We use default value 1\n", meshinfo->name);
        }
        else if (!strcmp (attrs, "maximums"))
        {
//            meshinfo->uniform->maximums = (double *) malloc (sizeof(double)*meshinfo->uniform->num_dimensions);
//            for (i = 0; i < meshinfo->uniform->num_dimensions; i++ )
//                meshinfo->uniform->maximums[i] = 0;
            printf("WARNING: uniform mesh %s does not provide maximum info.\n", meshinfo->name);
        }

    }

    if (have_max || have_spacing)
        return 1;
    else
        return 0;

}


ADIOS_MESH * common_read_inq_mesh_byid (ADIOS_FILE *fp, int meshid)
{
    enum ADIOS_DATATYPES attr_type;
    int attr_size;
    void * data = NULL;
    int  read_fail = 0;
    int i, j;

    ADIOS_MESH * meshinfo = (ADIOS_MESH *) malloc (sizeof(ADIOS_MESH));
    //mesh id
    meshinfo->id = meshid;
    //mesh name
    meshinfo->name = strdup(fp->mesh_namelist[meshinfo->id]);
    //construct /adios_schema/***/time-varying
    //mesh time-varying
    char * time_varying = malloc ( strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/time-varying")+1 );
    strcpy (time_varying, "/adios_schema/");
    strcat (time_varying, meshinfo->name);
    strcat (time_varying, "/time-varying");
    read_fail = common_read_get_attr_mesh (fp, time_varying, &attr_type, &attr_size, &data);
    free (time_varying);
    if (read_fail)
        meshinfo->time_varying = 0;
    else
    {
        if ( !strcmp ((char *)data,"yes") )
            meshinfo->time_varying = 1;  
        else
            meshinfo->time_varying = 0;
    }
    //construct /adios_schema/***/type to get mesh type (uniform, rectilinear, structured, unstructured)
    //mesh type
    char * mesh_attribute = malloc ( strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/type")+1 );
    strcpy (mesh_attribute, "/adios_schema/");
    strcat (mesh_attribute, meshinfo->name);
    strcat (mesh_attribute, "/type");
    common_read_get_attr_mesh (fp, mesh_attribute, &attr_type, &attr_size, &data);
    free (mesh_attribute);
    if ( !strcmp((char *)data, "uniform") )
    {
        bool have_spacing = 0;
        bool have_max = 0;

        meshinfo->type = ADIOS_MESH_UNIFORM;
        meshinfo->uniform = (MESH_UNIFORM * ) malloc (sizeof(MESH_UNIFORM));

        // initialize pointers that might not be set below
        meshinfo->uniform->spacings = NULL;
        meshinfo->uniform->maximums = NULL;
        meshinfo->uniform->origins = NULL;
        
        char * dimension_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/dimensions-num")+1 );
        strcpy (dimension_attribute, "/adios_schema/");
        strcat (dimension_attribute, meshinfo->name);
        strcat (dimension_attribute, "/dimensions-num");
        data = NULL;
        read_fail = common_read_get_attr_mesh (fp, dimension_attribute, &attr_type, &attr_size, &data);
        free (dimension_attribute);
//        printf("dimension_attribute is %s\n", dimension_attribute);
//        printf("dimension is %d\n",*(int *)data);
        if (!read_fail)
        {
            meshinfo->uniform->num_dimensions = *(int *)data;
            meshinfo->uniform->dimensions = (uint64_t *) malloc (sizeof(uint64_t)*meshinfo->uniform->num_dimensions);
            for (i = 0; i < meshinfo->uniform->num_dimensions; i++ )
            {
                meshinfo->uniform->dimensions[i] = 0;
                char * i_buffer;
                int i_digits;
                if (meshinfo->uniform->num_dimensions < 10)
                    i_buffer = (char *) malloc (sizeof(char)+1);
                else
                {
                    printf("ERROR: uniform mesh %s points have more than 10 dims!\n", meshinfo->name);
                    return NULL; 
                }
                i_digits = sprintf (i_buffer, "%d", i);
                //i_digits to reprent the number in dimensions0/dimensions1/... 
                char * dimensions_value = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/dimensions")+i_digits+1 );
                strcpy (dimensions_value, "/adios_schema/");
                strcat (dimensions_value, meshinfo->name);
                strcat (dimensions_value, "/dimensions");
                strcat (dimensions_value, i_buffer);
                free (i_buffer);
                data = NULL;
                read_fail = common_read_get_attr_mesh (fp, dimensions_value, &attr_type, &attr_size, &data);
                free (dimensions_value);
                if (read_fail)
                {
                    printf ("ERROR: mesh %s dimensions[%d] is not provided!\n", meshinfo->name, i);
                    return NULL; 
                }
                else
                {
                    char * pEnd;
                    char * tmp_dimensions_value = strdup((char *)data);
                    uint64_t tmp_value = strtoull (tmp_dimensions_value, &pEnd, 10);
                    if (tmp_value)
                        meshinfo->uniform->dimensions[i] = tmp_value;
                    else
                    {
                        int var_march = 0;
                        for (j=0; j<fp->nvars; j++)
                        {
                            char * dimensions_value_tmp = NULL;
                            dimensions_value_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                            strcpy (dimensions_value_tmp, "/");
                            strcat (dimensions_value_tmp, (char *)data);
//                            printf("dimensions_value_tmp is %s\n", dimensions_value_tmp);
//                            printf("var name is %s\n", fp->var_namelist[j]);
                            if (!strcmp (fp->var_namelist[j], dimensions_value_tmp))
                            {
                                var_march = 1;
                                ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                                meshinfo->uniform->dimensions[i] = *(uint64_t *)v->value;
                                common_read_free_varinfo (v);
                            }
                            free (dimensions_value_tmp);
                        }
                        if (!var_march)
                        {
                            printf("ERROR: mesh %s dimensions%d var %s is not correct!\n", meshinfo->name, i, (char *)data);
                            return NULL; 
                        }
                    }
                }  
            }
        }
        else
        {
                double d1;
            printf ("ERROR: mesh %s dimension is required\n", meshinfo->name);
            return NULL;            
        }

        //start processing origins, origin is optional
        adios_get_uniform_mesh_attr (fp, meshinfo, "origins");
//        for (i = 0; i < meshinfo->uniform->num_dimensions; i++ )
//            printf ("origins[%d] is %lf\n", i, meshinfo->uniform->origins[i]);

        //start processing maximums, maximum is optional 
        int have_maximums;
        have_maximums = adios_get_uniform_mesh_attr (fp, meshinfo, "maximums");
//        for (int i = 0; i < meshinfo->uniform->num_dimensions; i++ )
//            printf ("maximums[%d] is %lf\n", i, meshinfo->uniform->maximums[i]);

        //start processing spacings, spacing is optional 
        int have_spacings;
        have_spacings = adios_get_uniform_mesh_attr (fp, meshinfo, "spacings");
//        for (i = 0; i < meshinfo->uniform->num_dimensions; i++ )
//            printf ("spacings[%d] is %lf\n", i, meshinfo->uniform->spacings[i]);

        //if mesh spacing and maximum are both defined, check if consistant
        if (have_spacings == 1 && have_maximums == 1)
        {
            for (i = 0; i < meshinfo->uniform->num_dimensions; i++)
            {
                if (meshinfo->uniform->spacings[i] != (meshinfo->uniform->maximums[i]-meshinfo->uniform->origins[i])/(double)(meshinfo->uniform->dimensions[i]-1))
                {
//                    printf ( "meshinfo->uniform->origins = %lf,  meshinfo->uniform->dimensions = %lf\n", meshinfo->uniform->origins[i], meshinfo->uniform->dimensions[i]);
                    printf ("ERROR: mesh %s provided uniform mesh spacing and maximum is not consistant provided spacing is %lf, calculated spacing is %lf\n", 
                            meshinfo->name, meshinfo->uniform->spacings[i], (meshinfo->uniform->maximums[i]-meshinfo->uniform->origins[i])/(double)(meshinfo->uniform->dimensions[i]-1));
                    return NULL; 
                }
            }
        }
//end of uniform mesh 
    }
    else if ( !strcmp((char *)data, "rectilinear") )   
    {
        meshinfo->type = ADIOS_MESH_RECTILINEAR;
        meshinfo->rectilinear = (MESH_RECTILINEAR * ) malloc (sizeof(MESH_RECTILINEAR));
        meshinfo->rectilinear->use_single_var = 0;    // default value 0 indicates using multi-var

        // start processing rectilinear mesh dimensions
        char * dimension_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/dimensions-num")+1 );
        strcpy (dimension_attribute, "/adios_schema/");
        strcat (dimension_attribute, meshinfo->name);
        strcat (dimension_attribute, "/dimensions-num");
        data = NULL;
        read_fail = common_read_get_attr_mesh (fp, dimension_attribute, &attr_type, &attr_size, &data);
        free (dimension_attribute);
        if (!read_fail)
        {
            meshinfo->rectilinear->num_dimensions = *(int *)data;
//            printf ("the dimension is %d\n", meshinfo->rectilinear->num_dimensions);
            meshinfo->rectilinear->dimensions = (uint64_t *) malloc (sizeof(uint64_t)*meshinfo->rectilinear->num_dimensions);
            for (i = 0; i < meshinfo->rectilinear->num_dimensions; i++ )
            {
                meshinfo->rectilinear->dimensions[i] = 0;
                char * i_buffer;
                int i_digits;
                if (meshinfo->rectilinear->num_dimensions < 10)
                    i_buffer = (char *) malloc (sizeof(char)+1);
                else
                {
                    printf("ERROR: rectilinear mesh %s points have more than 10 dims!\n", meshinfo->name);
                    return NULL; 
                }
                i_digits = sprintf (i_buffer, "%d", i);
                //i_digits to reprent the number in dimensions0/dimensions1/... 
                char * dimensions_value = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/dimensions")+i_digits+1 );
                strcpy (dimensions_value, "/adios_schema/");
                strcat (dimensions_value, meshinfo->name);
                strcat (dimensions_value, "/dimensions");
                strcat (dimensions_value, i_buffer);
                free (i_buffer);
                data = NULL;
                read_fail = common_read_get_attr_mesh (fp, dimensions_value, &attr_type, &attr_size, &data);
//                printf("dimensions_value is %s\n", dimensions_value);
                free (dimensions_value);
                if (read_fail)
                {
                    printf ("ERROR: mesh %s dimensions[%d] is not provided!\n", meshinfo->name, i);
                    return NULL; 
                }
                else
                {
                    char * pEnd;
                    char * tmp_dimensions_value = strdup((char *)data);
                    uint64_t tmp_value = strtoull (tmp_dimensions_value, &pEnd, 10);
                    if (tmp_value)
                        meshinfo->rectilinear->dimensions[i] = tmp_value;
                    else
                    {
                        int var_march = 0;
                        for (j=0; j<fp->nvars; j++)
                        {
                            char * dimensions_value_tmp = NULL;
                            dimensions_value_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                            strcpy (dimensions_value_tmp, "/");
                            strcat (dimensions_value_tmp, (char *)data);
//                          printf("dimensions_value_tmp is %s\n", dimensions_value_tmp);
//                          printf("var name is %s\n", fp->var_namelist[j]);
                            if (!strcmp (fp->var_namelist[j], dimensions_value_tmp))
                            {
                                var_march = 1;
                                ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                                meshinfo->rectilinear->dimensions[i] = *(uint64_t *)v->value;
                                common_read_free_varinfo (v);
                            }
                            free (dimensions_value_tmp);
                        }
                        if (!var_march)
                        {
                            printf("ERROR: mesh %s dimensions%d var %s is not correct!\n", meshinfo->name, i, (char *)data);
                            return NULL; 
                        }
                    }
                }
            }
        }
        else
        {
            printf ("ERROR: rectilinear mesh %s dimension is required\n", meshinfo->name);
            return NULL; 
        }
//        for (int i = 0; i < meshinfo->rectilinear->num_dimensions; i++ )
//            printf ("dimensions[%d] is %d\n", i, meshinfo->rectilinear->dimensions[i]);

        // start processing rectilinear mesh coordinates
        char * coords_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/coordinates-single-var")+1 );
        strcpy (coords_attribute, "/adios_schema/");
        strcat (coords_attribute, meshinfo->name);
        strcat (coords_attribute, "/coords-single-var");
        read_fail = common_read_get_attr_mesh (fp, coords_attribute, &attr_type, &attr_size, &data);
//        printf("coords_attribute is %s\n", coords_attribute);
        free (coords_attribute);
        if (!read_fail)   //use coordinates-single-var
        {
            meshinfo->rectilinear->coordinates = (char **) malloc (sizeof(char *));
            meshinfo->rectilinear->use_single_var = 1;
            for (i=0; i<fp->nvars; i++)
            {
                char * coords_tmp = NULL;
                int var_march = 0;
                coords_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                strcpy (coords_tmp, "/");
                strcat (coords_tmp, (char *)data);
                if (!strcmp (fp->var_namelist[i], coords_tmp))
                {
                    var_march = 1;
                    ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[i]);
                    // check if mesh dimension matches coordinates dimension
                    if (v->ndim == 0)                      //scalar
                    {    
                        printf ("ERROR: mesh %s coordinates dimension is 0\n", meshinfo->name);
                        return NULL; 
                    }
                    else                                   //vector, more than one dim
                    {
                        if (v->ndim == 1)                  //one dim vector
                        {
                            uint64_t dim_tmp = 1;
                            for (j=0; j<meshinfo->rectilinear->num_dimensions; j++)
                                dim_tmp *= meshinfo->rectilinear->dimensions[j];
                            if (dim_tmp*meshinfo->rectilinear->num_dimensions != v->dims[0])
                            {
                                printf ("ERROR: mesh %s coordinates dimension %"PRIu64" does not match mesh dimension %"PRIu64"\n", 
                                        meshinfo->name, v->dims[0], dim_tmp*meshinfo->rectilinear->num_dimensions);
                                return NULL; 
                            }
                            else
                                meshinfo->rectilinear->coordinates[0] = strdup (fp->var_namelist[i]);
                        }
                        else        //multi dim vector
                        {
                            // check if each mesh dimension matches coordinates dimension 
                            for (j=0; j<v->ndim; j++)
                            {
                                if (meshinfo->rectilinear->dimensions[j] != v->dims[j])
                                {
                                    printf ("ERROR: mesh %s dimension[%d]= %d does not match coordinates dimension[%d]= %d\n",
                                        meshinfo->name, j, meshinfo->rectilinear->dimensions[j], j, v->dims[j] );
                                    return NULL; 
                                }
                            }
                            meshinfo->rectilinear->coordinates[0] = strdup (fp->var_namelist[i]); 
                        }
                    }
                    common_read_free_varinfo (v);
                }
                free (coords_tmp);
            }
//            printf ("coordinates is %s\n", meshinfo->rectilinear->coordinates[0]);
        }
        else                              // use coordinates-multi-var?
        {
            coords_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/coords-multi-var-num")+1 );
            strcpy (coords_attribute, "/adios_schema/");
            strcat (coords_attribute, meshinfo->name);
            strcat (coords_attribute, "/coords-multi-var-num");
            data = NULL;
            read_fail = common_read_get_attr_mesh (fp, coords_attribute, &attr_type, &attr_size, &data);
            free (coords_attribute);
            if (!read_fail)  //found attributes coords-multi-var
            {
                int num_coordinates = *(int *)data;
                if (num_coordinates < meshinfo->rectilinear->num_dimensions)
                {
                    printf("ERROR: mesh %s number of coordinates %d is less than the number of dimensions %d!\n", 
                            meshinfo->name, num_coordinates, meshinfo->rectilinear->num_dimensions);
                    return NULL; 
                }
                if (num_coordinates > meshinfo->rectilinear->num_dimensions)
                {
                    num_coordinates = meshinfo->rectilinear->num_dimensions;
                    printf("WARNING: mesh %s number of coordinates %d is greater than the number of dimensions %d! Use number of dimensions \
                               for maximums and ignore the rest maximums\n", meshinfo->name, num_coordinates, meshinfo->rectilinear->num_dimensions);
                }
                meshinfo->rectilinear->coordinates = (char **) malloc (num_coordinates*sizeof(char *)); 
                for (i=0; i<num_coordinates; i++)
                {
                    char * i_buffer;
                    int i_digits;
                    if (num_coordinates < 10)
                        i_buffer = (char *) malloc (sizeof(char)+1);
                    else
                    {
                        printf("ERROR: rectilinear mesh %s has more than 10 coordinates!\n", meshinfo->name);
                        return NULL; 
                    }
                    i_digits = sprintf (i_buffer, "%d", i);
 //                   printf ("i digit is %d\n", i_digits); 
                    char * coords_var = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/coords-multi-var")+i_digits+1 );
                    strcpy (coords_var, "/adios_schema/");
                    strcat (coords_var, meshinfo->name);
                    strcat (coords_var, "/coords-multi-var");
                    strcat (coords_var, i_buffer);
//                    printf ("coords_var is %s\n", coords_var);
                    free (i_buffer);
                    data = NULL;
                    read_fail = common_read_get_attr_mesh (fp, coords_var, &attr_type, &attr_size, &data);
                    free (coords_var);
                    if (read_fail)
                    {
                        printf ("ERROR: mesh %s coordinate of coordinate[%d] is not provided!\n", meshinfo->name, i);
                        return NULL; 
                    }
                    int var_march = 0;
                    for (j=0; j<fp->nvars; j++)
                    {
                        char * coords_var_tmp;
                        coords_var_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                        strcpy (coords_var_tmp, "/");
                        strcat (coords_var_tmp, (char *)data);
                        if (!strcmp (fp->var_namelist[j], coords_var_tmp))
                        {
                            var_march = 1;
//                            printf ("match found \n");
//                            printf ("var_namelist is %s\n", fp->var_namelist[j]);
//                            printf ("coords_var_tmp is %s\n", coords_var_tmp);
                            ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
//                            printf ("meshinfo->rectilinear->dimensions is %d\n", meshinfo->rectilinear->dimensions[i]);
//                            printf ("dims[0] is %d \n", v->dims[0]);
                            if (meshinfo->rectilinear->dimensions[i] != v->dims[0])
                            {
                                printf ("ERROR: mesh %s dimension[%d] = %d does not match coordinates dimension[%d] = %d\n",
                                    meshinfo->name, i, meshinfo->rectilinear->dimensions[i], i, v->dims[0] );
                                return NULL; 
                            }
                            else
                            {
                                meshinfo->rectilinear->coordinates[i] = strdup (fp->var_namelist[j]);
//                                printf ("coordinates[%d] is %s\n", i, meshinfo->rectilinear->coordinates[i]);
                            }
                            common_read_free_varinfo (v);
                        }
                        free (coords_var_tmp);
                    }
                }
            }
            else                //coords-single-var not found, coords-multi-var not found
            {
                printf ("ERROR: rectilinear mesh %s coordinate is not provided\n", meshinfo->name);
                return NULL; 
            }
//            for (int i = 0; i < meshinfo->rectilinear->num_dimensions; i++ )
//                printf ("coordinates[%d] is %s\n", i, meshinfo->rectilinear->coordinates[i]);
        }
    }
    else if ( !strcmp((char *)data, "structured") )
    {
        meshinfo->type = ADIOS_MESH_STRUCTURED;
        meshinfo->structured = (MESH_STRUCTURED* ) malloc (sizeof(MESH_STRUCTURED));
        meshinfo->structured->use_single_var = 0;        // default value 0 indicates using multi-var   
        meshinfo->structured->nspaces = meshinfo->structured->num_dimensions;   //default spaces = # of dims
 
        // start processing structured mesh dimensions
        char * dimension_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/dimensions-num")+1 );
        strcpy (dimension_attribute, "/adios_schema/");
        strcat (dimension_attribute, meshinfo->name);
        strcat (dimension_attribute, "/dimensions-num");
        data = NULL;
        read_fail = common_read_get_attr_mesh (fp, dimension_attribute, &attr_type, &attr_size, &data);
        free (dimension_attribute);
        if (!read_fail)
        {
            meshinfo->structured->num_dimensions = *(int *)data;
            meshinfo->structured->dimensions = (uint64_t *) malloc (sizeof(uint64_t)*meshinfo->structured->num_dimensions);
            for (i = 0; i < meshinfo->structured->num_dimensions; i++ )
            {
                meshinfo->structured->dimensions[i] = 0;
                char * i_buffer;
                int i_digits;
                if (meshinfo->structured->num_dimensions < 10)
                    i_buffer = (char *) malloc (sizeof(char)+1);
                else
                {
                    printf("ERROR: strctured mesh %s has more than 10 dimensions!\n", meshinfo->name);
                    return NULL; 
                }
                i_digits = sprintf (i_buffer, "%d", i);
                char * dimensions_value = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/dimensions")+i_digits+1 );
                strcpy (dimensions_value, "/adios_schema/");
                strcat (dimensions_value, meshinfo->name);
                strcat (dimensions_value, "/dimensions");
                strcat (dimensions_value, i_buffer);
                free (i_buffer);
                data = NULL;
                read_fail = common_read_get_attr_mesh (fp, dimensions_value, &attr_type, &attr_size, &data);
                if (read_fail)
                {
                    printf ("ERROR: strctured mesh %s dimension of dimensions[%d] is not provided!\n", meshinfo->name, i);
                    return NULL; 
                }
                else
                {
                    char * pEnd;
                    char * tmp_dimensions_value = strdup((char *)data);
                    uint64_t tmp_value = strtoull (tmp_dimensions_value, &pEnd, 10);
                    if (tmp_value)
                        meshinfo->structured->dimensions[i] = tmp_value;
                    else
                    {
                        int var_march = 0;
                        for (j=0; j<fp->nvars; j++)
                        {
                            char * dimensions_value_tmp = NULL;
                            dimensions_value_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                            strcpy (dimensions_value_tmp, "/");
                            strcat (dimensions_value_tmp, (char *)data);
                            if (!strcmp (fp->var_namelist[j], dimensions_value_tmp))
                            {
                                var_march = 1;
                                ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                                meshinfo->structured->dimensions[i] = *(uint64_t *)v->value;
                                common_read_free_varinfo (v);
                            }
                            free (dimensions_value_tmp);
                        }
                        if (!var_march)
                        {
                            printf("ERROR: mesh %s dimensions%d var %s is not correct!\n", meshinfo->name, i, (char *)data);
                            return NULL; 
                        }
                    }
                }
            }
        }
        else
        {
            printf ("ERROR: strctured mesh %s dimension is required\n", meshinfo->name);
            return NULL; 
        }

        // start processing structured mesh points
        char * points_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/coordinates-single-var")+1 );
        strcpy (points_attribute, "/adios_schema/");
        strcat (points_attribute, meshinfo->name);
        strcat (points_attribute, "/points-single-var");
        read_fail = common_read_get_attr_mesh (fp, points_attribute, &attr_type, &attr_size, &data);
        free (points_attribute);
        if (!read_fail)   //use points-single-var
        {
            meshinfo->structured->points = (char **) malloc (sizeof(char *));
            meshinfo->structured->use_single_var = 1;         // modify default value to 1
            for (i=0; i<fp->nvars; i++)
            {
                char * coords_tmp = NULL;
                int var_march = 0;
                coords_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                strcpy (coords_tmp, "/");
                strcat (coords_tmp, (char *)data);
                if (!strcmp (fp->var_namelist[i], coords_tmp))
                {
                    var_march = 1;
                    ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[i]); 
                    if (v->ndim == 0)                      //scalar
                    {
                        printf ("ERROR: mesh %s points dimension is 0\n", meshinfo->name);
                        return NULL; 
                    }
                    else                                   //vector, more than one dim
                    {
                        if (v->ndim == 1)     //if points is 1D array
                        {
                            uint64_t dim_tmp = 1;
                            for (j=0; j<meshinfo->structured->num_dimensions; j++)
                                dim_tmp *= meshinfo->structured->dimensions[j];
                            if (dim_tmp*meshinfo->structured->num_dimensions != v->dims[0])
                            {
                                printf ("ERROR: strctured mesh %s points dimension %"PRIu64" does not match mesh dimension %"PRIu64"\n", 
                                        meshinfo->name, v->dims[0], dim_tmp*meshinfo->structured->num_dimensions);
                                return NULL; 
                            }
                            else
                                meshinfo->structured->points[0] = strdup (fp->var_namelist[i]);
                        }
                        else
                        {
                            // check if each mesh dimension matches points dimension 
                            for (j=0; j<v->ndim; j++)              // if points is multi dim array
                            {
                                if (meshinfo->structured->dimensions[j] != v->dims[j])
                                {
                                    printf ("ERROR: mesh %s dimension[%d]= %"PRIu64" does not match points dimension[%d]= %"PRIu64"\n",
                                        meshinfo->name, j, meshinfo->structured->dimensions[j], j, v->dims[j] );
                                    return NULL; 
                                }
                            }
                            meshinfo->structured->points[0] = strdup (fp->var_namelist[i]);
                        }
                    }
                    common_read_free_varinfo (v);
                }
                free (coords_tmp);
            }
        }
        else                    //use points-multi-var
        {
            points_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/coords-multi-var-num")+1 );
            strcpy (points_attribute, "/adios_schema/");
            strcat (points_attribute, meshinfo->name);
            strcat (points_attribute, "/points-multi-var-num");
            data = NULL;
            read_fail = common_read_get_attr_mesh (fp, points_attribute, &attr_type, &attr_size, &data);
            free (points_attribute);
            if (!read_fail)  //found attributes points-multi-var
            {
                int points_dim = *(int *)data;
                if (points_dim < meshinfo->structured->num_dimensions)
                {
                    printf("ERROR: strctured mesh %s provided points dim %d is less than dims of mesh %d!\n",
                            meshinfo->name, points_dim, meshinfo->structured->num_dimensions);
                    return NULL; 
                }
                if (points_dim > meshinfo->structured->num_dimensions)
                {
                    points_dim = meshinfo->structured->num_dimensions;
                    printf("WARNING: strctured mesh %s points dim %d is greater than mesh dimes %d! Use number of mesh dimensions \
                               for points and ignore the rest points\n", meshinfo->name, points_dim, meshinfo->structured->num_dimensions);
                }
                meshinfo->structured->points = (char **) malloc (points_dim*sizeof(char *)); 
                for (i=0; i<points_dim; i++)
                {
                    char * i_buffer;
                    int i_digits;
                    if (points_dim < 10)
                        i_buffer = (char *) malloc (sizeof(char)+1);
                    else
                    {
                        printf("ERROR: structured mesh %s points have more than 10 dims!\n", meshinfo->name);
                        return NULL; 
                    }
                    i_digits = sprintf (i_buffer, "%d", i);
                    char * points_var = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/points-multi-var")+i_digits+1 );
                    strcpy (points_var, "/adios_schema/");
                    strcat (points_var, meshinfo->name);
                    strcat (points_var, "/points-multi-var");
                    strcat (points_var, i_buffer);
                    free (i_buffer);
                    data = NULL;
                    read_fail = common_read_get_attr_mesh (fp, points_var, &attr_type, &attr_size, &data);
                    free (points_var);
                    if (read_fail)
                    {
                        printf ("ERROR: strctured mesh %s points of dim[%d] is not provided!\n", meshinfo->name, i);
                        return NULL; 
                    }
                    int var_march = 0;
                    for (j=0; j<fp->nvars; j++)
                    {
                        char * points_var_tmp;
                        points_var_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                        strcpy (points_var_tmp, "/");
                        strcat (points_var_tmp, (char *)data);
                        if (!strcmp (fp->var_namelist[j], points_var_tmp))
                        {
                            var_march = 1;
                            ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                            //check if dim of mesh matches point dim
                            if (v->ndim == 1)     //if var is 1D array
                            {
                                uint64_t dim_tmp = 1;
                                int m = 0;
                                for (m=0; m<meshinfo->structured->num_dimensions; m++)
                                    dim_tmp *= meshinfo->structured->dimensions[m];
                                if (dim_tmp != v->dims[0])
                                {
                                    printf ("ERROR: strctured mesh %s points dimension %"PRIu64" does not match mesh dimension %"PRIu64"\n", 
                                            meshinfo->name, v->dims[0], dim_tmp);
                                    return NULL; 
                                }
                                else
                                    meshinfo->structured->points[i] = strdup (fp->var_namelist[j]);
                            }
                            else
                            {
                                // check if each mesh dim matches points dim (var is multi dim array)
                                int m = 0;
                                for (m=0; m<v->ndim; m++)              // if points is multi dim array
                                {
                                    if (meshinfo->structured->dimensions[m] != v->dims[m])
                                    {
                                        printf ("ERROR: mesh %s dimension[%d]= %"PRIu64" does not match points dimension[%d]= %"PRIu64"\n",
                                            meshinfo->name, m, meshinfo->structured->dimensions[m], m, v->dims[m] );
                                        return NULL; 
                                    }
                                }
                                meshinfo->structured->points[i] = strdup (fp->var_namelist[j]);
                            }
                            common_read_free_varinfo (v);
                        }
                        free (points_var_tmp);
                    }  // end of for loop j
                    if (var_march == 0)
                    {
                        printf ("ERROR: structured mesh %s var of points-multi-var%d is not provided!\n", meshinfo->name, i); 
                        return NULL; 
                    }
                }   // end of for loop i
            }  // end of if found attributes points-multi-var
            else
            {
                printf ("ERROR: structured mesh %s point is not provided\n", meshinfo->name);
                return NULL; 
            }

        }  // end of else use points-multi-var

        // start processing structured mesh nspace
        char * mesh_nspace = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/nspace")+1 );
        strcpy (mesh_nspace, "/adios_schema/");
        strcat (mesh_nspace, meshinfo->name);
        strcat (mesh_nspace, "/nspace");
        data = NULL;
        read_fail = common_read_get_attr_mesh (fp, mesh_nspace, &attr_type, &attr_size, &data);
        free (mesh_nspace);
        if (read_fail)
            printf ("WARNING: structured mesh %s nspace is not provided\n", meshinfo->name);
        else
        {
            long int d1;
            char * pEnd;
            d1 = strtol((char *)data, &pEnd, 10);
            if (d1)
                meshinfo->structured->nspaces = d1;
            else
            {
                int var_march = 0;
                char * spaces_var_tmp;
                spaces_var_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                strcpy (spaces_var_tmp, "/");
                strcat (spaces_var_tmp, (char *)data);
//                printf ("spaces_var_tmp is %s\n", spaces_var_tmp);
                for (j=0; j<fp->nvars; j++)
                {
                    if (!strcmp (fp->var_namelist[j], spaces_var_tmp))
                    {
                        var_march = 1;
                        ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                        meshinfo->structured->nspaces = *(int *)v->value;
                        common_read_free_varinfo (v);
                    }
                    free (spaces_var_tmp);
                }
                if (!var_march)
                {
                    printf ("ERROR: structured mesh %s nspace is not correct\n", meshinfo->name);
                    return NULL;
                }
            }
        } // end of structured mesh nspace
    }// end of structured mesh
    else if ( !strcmp((char *)data, "unstructured") )
    {
        meshinfo->type = ADIOS_MESH_UNSTRUCTURED;
        meshinfo->unstructured = (MESH_UNSTRUCTURED* ) malloc (sizeof(MESH_UNSTRUCTURED));
//        meshinfo->unstructured->use_single_var = 0;  // default value 0 indicates using multi-var
        meshinfo->unstructured->nvar_points = 1;
//        meshinfo->unstructured->uniform_cell = 1;   // default value 0 indicates using uniform cell
//        meshinfo->unstrutured->nspaces init
//        meshinfo->unstructured->npoints init

        // start processing unstructured mesh points-single-var/points-multi-var
        char * points_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/coordinates-single-var")+1 );
        strcpy (points_attribute, "/adios_schema/");
        strcat (points_attribute, meshinfo->name);
        strcat (points_attribute, "/points-single-var");
        read_fail = common_read_get_attr_mesh (fp, points_attribute, &attr_type, &attr_size, &data);
        free (points_attribute);
        if (!read_fail)   //use points-single-var
        {
            meshinfo->unstructured->points = (char **) malloc (sizeof(char *));  
//            meshinfo->unstructured->use_single_var = 1;         // modify default value to 1
            int var_march = 0;
            char * coords_tmp = NULL;
            coords_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
            strcpy (coords_tmp, "/");
            strcat (coords_tmp, (char *)data);
//            printf ("coords_tmp is %s\n", coords_tmp);
            for (i=0; i<fp->nvars; i++)
            {
                if (!strcmp (fp->var_namelist[i], coords_tmp))
                {
                    var_march = 1;
                    ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[i]);
                    if (v->ndim == 0)                      //scalar
                    {
                        printf ("ERROR: unstructured mesh %s points dimension is 0\n", meshinfo->name);
                        return NULL; 
                    }
                    else                                   //vector
                    {
                        meshinfo->unstructured->nspaces = v->ndim;                //unstructured mesh nspaces init
                        meshinfo->unstructured->npoints = 1;
                        int j = 0;
                        for (j=0; j<v->ndim; j++)
                            meshinfo->unstructured->npoints *= v->dims[j];         
                        meshinfo->unstructured->npoints /= v->ndim;               //unstructured mesh npoints init
                        meshinfo->unstructured->points[0] = strdup (fp->var_namelist[i]);
                    }
                    common_read_free_varinfo (v);
                }
            }
            if (var_march ==0)
            {
                printf ("unstrctured mesh %s var %s for points-single-var is not correct\n", meshinfo->name, coords_tmp);
            }
            free (coords_tmp);
        }
        else                    //use points-multi-var
        {
            points_attribute = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/coords-multi-var-num")+1 );
            strcpy (points_attribute, "/adios_schema/");
            strcat (points_attribute, meshinfo->name);
            strcat (points_attribute, "/points-multi-var-num");
            data = NULL;
            read_fail = common_read_get_attr_mesh (fp, points_attribute, &attr_type, &attr_size, &data);
            free (points_attribute);
            if (!read_fail)  //found attributes points-multi-var
            {
                meshinfo->unstructured->nvar_points = *(int *)data;
                int points_dim = *(int *)data;
                meshinfo->unstructured->points = (char **) malloc (points_dim*sizeof(char *));
                meshinfo->unstructured->nspaces = points_dim;       //unstructured mesh nspaces init
                int first_dim = 1;
                meshinfo->unstructured->npoints = 1;
                char * first_match_var;
                for (i=0; i<points_dim; i++)
                {
                    char * i_buffer;
                    int i_digits;
                    if (points_dim < 10)
                        i_buffer = (char *) malloc (sizeof(char)+1);
                    else
                    {
                        printf("ERROR: structured mesh %s points have more than 10 dims!\n", meshinfo->name);
                        return NULL; 
                    }
                    i_digits = sprintf (i_buffer, "%d", i);
                    char * points_var = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/points-multi-var")+i_digits+1 );
                    strcpy (points_var, "/adios_schema/");
                    strcat (points_var, meshinfo->name);
                    strcat (points_var, "/points-multi-var");
                    strcat (points_var, i_buffer);
                    free (i_buffer);
                    data = NULL;
                    read_fail = common_read_get_attr_mesh (fp, points_var, &attr_type, &attr_size, &data);
                    free (points_var);
                    if (read_fail)
                    {
                        printf ("ERROR: unstructured mesh %s points of dim[%d] is not provided!\n", meshinfo->name, i);
                        return NULL; 
                    }
                    int var_march = 0;
                    char * points_var_tmp;
                    points_var_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                    strcpy (points_var_tmp, "/");
                    strcat (points_var_tmp, (char *)data);
//                    printf ( "points_var_tmp is %s\n", points_var_tmp);
                    for (j=0; j<fp->nvars; j++)
                    {
//                        printf ("var is %s\n", fp->var_namelist[j]);
                        if (!strcmp (fp->var_namelist[j], points_var_tmp))
                        {
                            var_march = 1;
                            ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                            if (first_dim)
                            {
                                first_match_var = strdup (points_var_tmp);
                                int k = 0;
                                for (k=0; k<v->ndim; k++)
                                    meshinfo->unstructured->npoints *= v->dims[k];
                                first_dim = 0;
                            }
                            else
                            {
                                if (v->ndim == 1)
                                {
                                    if (meshinfo->unstructured->npoints != v->dims[0])
                                    {
                                        printf ("ERROR: unstructured mesh %s points-multi-var%d %"PRIu64" does not match points-multi-var0 %"PRIu64"\n", 
                                            meshinfo->name, i, v->dims[0], meshinfo->unstructured->npoints);
                                        return NULL; 
                                    }
                                }
                                else // v->ndim > 1, check if match for each dim
                                {
                                    ADIOS_VARINFO * v_first = common_read_inq_var(fp, first_match_var);
                                    if (v_first->ndim == 1)
                                    {
                                        int k = 0;
                                        uint64_t var_dim_tmp = 1;
                                        for (k=0; k<v->ndim; k++)
                                            var_dim_tmp *= v->dims[k];
                                        if (var_dim_tmp != meshinfo->unstructured->npoints)
                                        {
                                            printf ("ERROR: unstructured mesh %s points-multi-var%d %"PRIu64" does not match points-multi-var0 %"PRIu64"\n",
                                                meshinfo->name, i, var_dim_tmp, meshinfo->unstructured->npoints);
                                            return NULL; 
                                        }
                                    }
                                    else  //both v_first and v has more than 2 dims
                                    {
                                        if (v_first->ndim != v->ndim)
                                        {
                                            printf ("ERROR: unstructured mesh %s points-multi-var%d dim does not match points-multi-var0 dim\n",
                                                meshinfo->name, i);
                                            return NULL; 
                                        }
                                        else
                                        {
                                            int k = 0; 
                                            for (k=0; k<v->ndim; k++)
                                            {
                                                if (v->dims[k]!=v_first->dims[k])
                                                {
                                                    printf ("ERROR: unstructured mesh %s points-multi-var%d dim%d does not match points-multi-var0 dim%d\n",
                                                        meshinfo->name, i, k, k);
                                                    return NULL; 
                                                }
                                            }
                                        }
                                    }
                                    common_read_free_varinfo (v_first);
                                }
                            }
                            meshinfo->unstructured->points[i] = strdup (fp->var_namelist[j]);
                            common_read_free_varinfo (v);
                        }
                    }  // end of for loop j
                    free (points_var_tmp);
                    if (var_march == 0)
                    {
                        printf ("ERROR: unstructured mesh %s var of points-multi-var%d is not provided!\n", meshinfo->name, i);
                        return NULL; 
                    }
                }   // end of for loop i
            }  // end of if found attributes points-multi-var
            else
            {
                printf ("ERROR: unstructured mesh %s point is not provided\n", meshinfo->name);
                return NULL; 
            }

        }  // end of else use points-multi-var


        // start processing unstructured mesh npoints
        char * num_points = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/npoints")+1 );
        strcpy (num_points, "/adios_schema/");
        strcat (num_points, meshinfo->name);
        strcat (num_points, "/npoints");
        data = NULL;
        read_fail = common_read_get_attr_mesh (fp, num_points, &attr_type, &attr_size, &data);
        free (num_points);
        if (read_fail)
            printf ("WARNING: unstructured mesh %s npoints is not provided, use calculated default npoints %"PRIu64" from points-var\n", 
                    meshinfo->name, meshinfo->unstructured->npoints);
        else
        {
            uint64_t d1;
            char * pEnd;
            d1 = strtoull((char *)data, &pEnd, 10);
            if (d1)
            {
                if (meshinfo->unstructured->npoints != d1)
                {
                    printf ("ERROR: unstructured mesh %s provided npoints %"PRIu64" does not match points calculated from points-var %"PRIu64"\n", 
                            meshinfo->name, d1, meshinfo->unstructured->npoints);
                    return NULL;
                }
            }
            else
            {
                int var_march = 0;
                for (j=0; j<fp->nvars; j++)
                {
                    char * points_var_tmp;
                    points_var_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                    strcpy (points_var_tmp, "/");
                    strcat (points_var_tmp, (char *)data);
                    if (!strcmp (fp->var_namelist[j], points_var_tmp))
                    {
                        var_march = 1;
                        ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                        if (v->type==adios_long || v->type==adios_unsigned_long) //long long or unsigned long long   
                        {
                            if (meshinfo->unstructured->npoints != *(uint64_t *)v->value)
                                printf ("WARNING: provided npoints %"PRIu64" does not match points calculated from points-var %"PRIu64"\n", 
                                    *(uint64_t *)v->value, meshinfo->unstructured->npoints);
                        }
                        else      //int
                        {
                            if (meshinfo->unstructured->npoints != *(int *)v->value)
                                printf ("WARNING: provided npoints %d does not match points calculated from points-var %"PRIu64"\n",
                                    *(int *)v->value, meshinfo->unstructured->npoints);
                        }
                        common_read_free_varinfo (v);
                    }
                    free (points_var_tmp);
                }
                if (!var_march)
                    log_warn ("Unstructured mesh %s var of npoints is not correct. "
                             "We use calculated default npoints %"PRIu64" from points-var\n", 
                            meshinfo->name, meshinfo->unstructured->npoints);
            }
        }

        // start processing unstructured mesh nspaces
        char * mesh_nspace = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/nspace")+1 );
        strcpy (mesh_nspace, "/adios_schema/");
        strcat (mesh_nspace, meshinfo->name);
        strcat (mesh_nspace, "/nspace");
        data = NULL;
        read_fail = common_read_get_attr_mesh (fp, mesh_nspace, &attr_type, &attr_size, &data);
        free (mesh_nspace);
        if (read_fail) {
            log_info ("Unstructured mesh %s nspace is not provided. "
                      "We use points dim %d for nspaces\n", 
                      meshinfo->name,  meshinfo->unstructured->nspaces);
        }
        else
        {
            int d1;
            char * pEnd;
            d1 = strtol((char *)data, &pEnd, 10);
            if (d1)
            {   
                if (meshinfo->unstructured->nspaces > d1) {
                    log_warn ("The provided nspaces %d is less the points dim %d. "
                              "We use points dim %d for nspaces\n",
                             d1, meshinfo->unstructured->nspaces, 
                             meshinfo->unstructured->nspaces);
                } else
                    meshinfo->unstructured->nspaces = d1;
            }
            else
            {
                int var_march = 0;
                for (j=0; j<fp->nvars; j++)
                {
                    char * spaces_var_tmp;
                    spaces_var_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                    strcpy (spaces_var_tmp, "/");
                    strcat (spaces_var_tmp, (char *)data);
                    if (!strcmp (fp->var_namelist[j], spaces_var_tmp))
                    {
                        var_march = 1;
                        ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                        if (meshinfo->unstructured->nspaces > *(int *)v->value) {
                            log_warn ("Unstructured mesh %s: the provided nspaces %d "
                                      "is less than the points dim %d. "
                                      "We use points dim %d for nspaces.\n",
                                      meshinfo->name, *(int *)v->value, 
                                      meshinfo->unstructured->nspaces, 
                                      meshinfo->unstructured->nspaces);
                        } else
                            meshinfo->unstructured->nspaces = *(int *)v->value;
                        common_read_free_varinfo (v);
                    }
                    free (spaces_var_tmp);
                }
                if (!var_march)
                    log_warn ("Unstructured mesh %s var of npoints is not correct, "
                              " use points dim %d for nspaces\n",
                              meshinfo->name, meshinfo->unstructured->nspaces);
            }
        }

        // start processing cells, including cell count, cell data and cell type
        // how many types of cells
        char * num_cell_type = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/ncsets")+1 );
        strcpy (num_cell_type, "/adios_schema/");
        strcat (num_cell_type, meshinfo->name);
        strcat (num_cell_type, "/ncsets");
        data = NULL;
        read_fail = common_read_get_attr_mesh (fp, num_cell_type, &attr_type, &attr_size, &data);
//        printf ("num_cell_type is %s\n", num_cell_type);
        free (num_cell_type);
        if (read_fail)
        {
            adios_error (err_mesh_unstructured_missing_ncsets, 
                        "Unstructured mesh %s ncsets is not provided, user should have "
                        "mixed-cells-count/uniform-cells in xml file\n", meshinfo->name);
            return NULL;
        }
        else
        {
            int d1 = *(int *)data;
            if (d1)
            {
                meshinfo->unstructured->ncsets = d1;
//                if (d1 == 1)
//                    meshinfo->unstructured->uniform_cell = 1;
//                else 
//                    meshinfo->unstructured->uniform_cell = 0;
            }
            else
            {
                adios_error (err_mesh_unstructured_invalid_ncsets, 
                            "Reading unstructured mesh %s ncsets failed\n", 
                            meshinfo->name);
                return NULL;

            }
        }

        // start processing ccount, how many cells for cell type
//        if (meshinfo->unstructured->uniform_cell)        //uniform cells
        if (meshinfo->unstructured->ncsets == 1)
        {
            meshinfo->unstructured->ccounts = (uint64_t *) malloc (sizeof(uint64_t));

            char * num_cells = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/ccount")+1 );
            strcpy (num_cells, "/adios_schema/");
            strcat (num_cells, meshinfo->name);
            strcat (num_cells, "/ccount");
//            printf ("num_cells is %s\n", num_cells);
            data = NULL;
            read_fail = common_read_get_attr_mesh (fp, num_cells, &attr_type, &attr_size, &data);
            free (num_cells);
            if (read_fail)
            {
                adios_error (err_mesh_unstructured_missing_ccount,
                            "unstructured mesh %s number of cells (ccount) is required\n", 
                            meshinfo->name);    
                return NULL;
            }
            else
            {
                uint64_t d1;
                char * pEnd;
                d1 = strtoull((char *)data, &pEnd, 10);  //number of cells
                if (d1)
                    meshinfo->unstructured->ccounts[0] = d1;
                else
                {
                    int var_march = 0;
                    char * ccount_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                    strcpy (ccount_tmp, "/");
                    strcat (ccount_tmp, (char *)data);
                    for (j=0; j<fp->nvars; j++)
                    {
                        if (!strcmp (fp->var_namelist[j], ccount_tmp))
                        {
                            var_march = 1;
                            ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                            if (v->type == adios_unsigned_long || v->type == adios_long)
                                    meshinfo->unstructured->ccounts[0] = *(uint64_t *)v->value;
                                else
                                    meshinfo->unstructured->ccounts[0] = *(int *)v->value;
                            common_read_free_varinfo (v); 
                        }
                    }
                    free (ccount_tmp);
                    if (var_march == 0)
                    {
                        adios_error (err_mesh_unstructured_invalid_ccount,
                                    "unstructured mesh %s var for ccount is invalid\n", 
                                    meshinfo->name);
                        return NULL;
                    }
                }
            }    
        }
        else  //mixed cells
        {
            meshinfo->unstructured->ccounts = (uint64_t *) malloc (sizeof(uint64_t)*meshinfo->unstructured->ncsets);
            int i = 0;
            for (i=0; i<meshinfo->unstructured->ncsets; i++)
            {
                char * i_buffer;
                int i_digits;
                if (meshinfo->unstructured->ncsets < 10)
                    i_buffer = (char *) malloc (sizeof(char)+1);
                else 
                {
                    adios_error (err_mesh_unstructured_invalid_ctypes,
                                "unstructured mesh %s has more than 10 cell types!\n", 
                                meshinfo->name);
                    return NULL;
                }
                i_digits = sprintf (i_buffer, "%d", i);
                char * ccount_var = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/ccount")+i_digits+1 );
                strcpy (ccount_var, "/adios_schema/");
                strcat (ccount_var, meshinfo->name);
                strcat (ccount_var, "/ccount");
                strcat (ccount_var, i_buffer);
                free (i_buffer);
                data = NULL;
                read_fail = common_read_get_attr_mesh (fp, ccount_var, &attr_type, &attr_size, &data);
                free (ccount_var);
                if (read_fail)
                {
                    adios_error (err_mesh_unstructured_missing_ccount,
                                "unstructured mesh %s ccount%d is not provided!\n", 
                                meshinfo->name, i);
                    return NULL; 
                }
                else
                {
                    uint64_t d1;
                    char * pEnd;
                    d1 = strtoull((char *)data, &pEnd, 10);
                    if (d1)
                        meshinfo->unstructured->ccounts[i] = d1;
                    else
                    {
                        int var_march = 0; 
                        char * ccount_mix_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                        strcpy (ccount_mix_tmp, "/");
                        strcat (ccount_mix_tmp, (char *)data);
                        for (j=0; j<fp->nvars; j++)
                        {
                            if (!strcmp (fp->var_namelist[j], ccount_mix_tmp))
                            {
                                var_march = 1;
                                ADIOS_VARINFO * v = common_read_inq_var(fp, fp->var_namelist[j]);
                                if (v->type == adios_long || v->type == adios_unsigned_long)
                                    meshinfo->unstructured->ccounts[i] = *(uint64_t *)v->value;
                                else
                                    meshinfo->unstructured->ccounts[i] = *(int *)v->value;
                                common_read_free_varinfo (v);
                            }
                        }
                        free (ccount_mix_tmp);
                        if (var_march == 0)
                        {
                            adios_error (err_mesh_unstructured_invalid_ccount,
                                        "unstructured mesh %s var for ccount%d is invalid\n", 
                                        meshinfo->name, i);
                            return NULL; 
                        }

                    }
                }
            }

        } // end of ccount
        
        // start processing cdata
//        if (meshinfo->unstructured->uniform_cell)        //uniform cells
        if (meshinfo->unstructured->ncsets == 1)
        {
            meshinfo->unstructured->cdata = (char **) malloc (sizeof(char *));
            char * data_cells = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/cdata")+1 );
            strcpy (data_cells, "/adios_schema/");
            strcat (data_cells, meshinfo->name);
            strcat (data_cells, "/cdata");
            data = NULL;
            read_fail = common_read_get_attr_mesh (fp, data_cells, &attr_type, &attr_size, &data);
            free (data_cells);
            if (read_fail)
            {
                adios_error (err_mesh_unstructured_missing_cdata,
                            "unstructured mesh %s cell data is required\n", 
                            meshinfo->name);
                return NULL; 
            }
            else
            {
                char * cdata_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                strcpy (cdata_tmp, "/");
                strcat (cdata_tmp, (char *)data);
                int var_match = 0;
                for (j=0; j<fp->nvars; j++)
                {
                    if (!strcmp (fp->var_namelist[j], cdata_tmp))
                    {
                        var_match = 1;
                        meshinfo->unstructured->cdata[0] = strdup(fp->var_namelist[j]);
                    }
                }
                free (cdata_tmp);
                if (var_match == 0)
                {
                    adios_error (err_mesh_unstructured_invalid_cdata,
                                "unstructured mesh %s var for cdata is invalid\n", 
                                meshinfo->name);
                    return NULL; 
                }
            }
        }
        else
        {
            meshinfo->unstructured->cdata = (char **) malloc (sizeof(char *)*meshinfo->unstructured->ncsets);
            int i = 0;
            for (i=0; i<meshinfo->unstructured->ncsets; i++)
            {
                char * i_buffer;
                int i_digits;
                if (meshinfo->unstructured->ncsets < 10)
                    i_buffer = (char *) malloc (sizeof(char)+1);
                else
                {
                    adios_error (err_mesh_unstructured_invalid_ctypes,
                                "unstructured mesh %s has more than 10 cell types!\n",
                                meshinfo->name);
                    return NULL; 
                }
                i_digits = sprintf (i_buffer, "%d", i);
                char * cdata_var = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/cdata")+i_digits+1 );
                strcpy (cdata_var, "/adios_schema/");
                strcat (cdata_var, meshinfo->name);
                strcat (cdata_var, "/cdata");
                strcat (cdata_var, i_buffer);
                free (i_buffer);
                data = NULL;
                read_fail = common_read_get_attr_mesh (fp, cdata_var, &attr_type, &attr_size, &data);
                free (cdata_var);
                if (read_fail)
                {
                    adios_error (err_mesh_unstructured_missing_cdata,
                                "unstructured mesh %s cdata%d is not provided!\n", 
                                meshinfo->name, i);
                    return NULL; 
                }
                else
                {
                    int var_match = 0;
                    char * cdata_mix_tmp = (char *) malloc (strlen("/")+strlen((char *)data)+1);
                    strcpy (cdata_mix_tmp, "/");
                    strcat (cdata_mix_tmp, (char *)data);
                    for (j=0; j<fp->nvars; j++)
                    {
                        if (!strcmp (fp->var_namelist[j], cdata_mix_tmp))
                        {
                            var_match = 1;
                            meshinfo->unstructured->cdata[i] = strdup(fp->var_namelist[j]);
                        }
                    }
                    if (var_match == 0)
                    {
                        adios_error (err_mesh_unstructured_invalid_cdata,
                                    "unstructured mesh %s var for cdata%d is not correct\n", 
                                    meshinfo->name, i);
                        return NULL; 
                    }
                }
            }
        }// end of cdata

        // start processing ctypes
//        if (meshinfo->unstructured->uniform_cell)        //uniform cells
        if (meshinfo->unstructured->ncsets == 1)
        {
            meshinfo->unstructured->ctypes = (enum ADIOS_CELL_TYPE *) malloc (sizeof(enum ADIOS_CELL_TYPE));
            char * type_cells = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/ctype")+1 );
            strcpy (type_cells, "/adios_schema/");
            strcat (type_cells, meshinfo->name);
            strcat (type_cells, "/ctype");
            data = NULL;
            read_fail = common_read_get_attr_mesh (fp, type_cells, &attr_type, &attr_size, &data);
            free (type_cells);
            if (read_fail)
            {
                adios_error (err_mesh_unstructured_missing_ctype,
                            "unstructured mesh %s cells type is required\n", 
                            meshinfo->name);
                return NULL; 
            } 
            else
            {
                if (!strcmp((char *)data, "line"))
                    meshinfo->unstructured->ctypes[i] = ADIOS_CELL_LINE;
                else if (!strcmp((char *)data, "triangle"))
                    meshinfo->unstructured->ctypes[i] = ADIOS_CELL_TRI;
                else if (!strcmp((char *)data, "quad"))
                    meshinfo->unstructured->ctypes[i] = ADIOS_CELL_QUAD;
                else if (!strcmp((char *)data, "hex"))
                    meshinfo->unstructured->ctypes[i] = ADIOS_CELL_HEX;
                else if (!strcmp((char *)data, "prism"))
                    meshinfo->unstructured->ctypes[i] = ADIOS_CELL_PRI;
                else if (!strcmp((char *)data, "tet"))
                    meshinfo->unstructured->ctypes[i] = ADIOS_CELL_TET;
                else if (!strcmp((char *)data, "pyr"))
                    meshinfo->unstructured->ctypes[i] = ADIOS_CELL_PYR;
                else
                {
                    adios_error (err_mesh_unstructured_invalid_ctype,
                                "unstructured mesh %s type %s of for ctype%d is invalid. " 
                                "we use line, triangle, quad, hex, prism, tet or tet for cell types. "
                                "please choose to use one of them. ",
                                 meshinfo->name, (char *)data, i);
                    return NULL;
                }
            }
        }
        else
        {
            meshinfo->unstructured->ctypes = (enum ADIOS_CELL_TYPE *) malloc (sizeof(enum ADIOS_CELL_TYPE)*meshinfo->unstructured->ncsets);
            int i = 0;
            for (i=0; i<meshinfo->unstructured->ncsets; i++)
            {
                char * i_buffer;
                int i_digits;
                if (meshinfo->unstructured->ncsets < 10)
                    i_buffer = (char *) malloc (sizeof(char)+1);
                else
                {
                    adios_error (err_mesh_unstructured_invalid_ctypes,
                                "unstructured mesh %s has more than 10 cell types!\n",
                                meshinfo->name);
                    return NULL; 
                }
                i_digits = sprintf (i_buffer, "%d", i);
                char * ctype_mix_var = malloc (strlen("/adios_schema/")+strlen(meshinfo->name)+strlen("/ctype")+i_digits+1 );
                strcpy (ctype_mix_var, "/adios_schema/");
                strcat (ctype_mix_var, meshinfo->name);
                strcat (ctype_mix_var, "/ctype");
                strcat (ctype_mix_var, i_buffer);
                free (i_buffer);
                data = NULL;
                read_fail = common_read_get_attr_mesh (fp, ctype_mix_var, &attr_type, &attr_size, &data);
                free (ctype_mix_var);
                if (read_fail)
                {
                    adios_error (err_mesh_unstructured_missing_ctype,
                                "unstructured mesh %s ctype%d is not provided!\n", 
                                meshinfo->name, i);
                    return NULL; 
                }
                else
                {
                    if (!strcmp((char *)data, "line"))
                        meshinfo->unstructured->ctypes[i] = ADIOS_CELL_LINE;
                    else if (!strcmp((char *)data, "triangle"))
                        meshinfo->unstructured->ctypes[i] = ADIOS_CELL_TRI;
                    else if (!strcmp((char *)data, "quad"))
                        meshinfo->unstructured->ctypes[i] = ADIOS_CELL_QUAD;
                    else if (!strcmp((char *)data, "hex"))
                        meshinfo->unstructured->ctypes[i] = ADIOS_CELL_HEX;
                    else if (!strcmp((char *)data, "prism")) 
                        meshinfo->unstructured->ctypes[i] = ADIOS_CELL_PRI;
                    else if (!strcmp((char *)data, "tet"))
                        meshinfo->unstructured->ctypes[i] = ADIOS_CELL_TET;
                    else if (!strcmp((char *)data, "pyr"))
                        meshinfo->unstructured->ctypes[i] = ADIOS_CELL_PYR;
                    else
                    {
                        adios_error (err_mesh_unstructured_invalid_ctype,
                                    "unstructured mesh %s type %s of for ctype%d is invalid. "
                                    "we use line, triangle, quad, hex, prism, tet or tet for cell types. "
                                    "please choose to use one of them. ", 
                                    meshinfo->name, (char *)data, i);
                        return NULL;
                        
                    }
                }
            }
        }
    }  // end of unstructured mesh
             
//    fp->attr_namelist[i]
//    common_read_get_attr_mesh (f, f->attr_namelist[i], &attr_type, &attr_size, &data);    

    return meshinfo;
}

void common_read_free_meshinfo (ADIOS_MESH * meshinfo)
{
    if(meshinfo)
    {
        int i = 0;
        switch (meshinfo->type) {
            case ADIOS_MESH_UNIFORM:
                {
                    MESH_UNIFORM *bp = meshinfo->uniform;
                    if (bp->dimensions)
                        free (bp->dimensions);
                    if (bp->origins)
                        free (bp->origins);
                    if (bp->spacings) 
                        free (bp->spacings);
                    if (bp->maximums)
                        free (bp->maximums);
                    free (meshinfo->uniform);
                    break;
                }
            case ADIOS_MESH_RECTILINEAR:
                {
                    MESH_RECTILINEAR *bp = meshinfo->rectilinear;
                    if (bp->dimensions)
                        free (bp->dimensions);
                    for (i = 0; i < meshinfo->rectilinear->num_dimensions; i++)
                    {
                        if (bp->coordinates[i])
                            free (bp->coordinates[i]);
                    }
                    free (meshinfo->rectilinear);
                    break;
                }
            case ADIOS_MESH_STRUCTURED:
                {
                    MESH_STRUCTURED *bp = meshinfo->structured;
                    if (bp->dimensions)
                        free (bp->dimensions);
                    for (i = 0; i < meshinfo->structured->num_dimensions; i++)
                    {
                        if (bp->points[i])
                            free (bp->points[i]);
                    }
                    free (meshinfo->unstructured);
                    break;
                }
            case ADIOS_MESH_UNSTRUCTURED:
                {
                    MESH_UNSTRUCTURED *bp = meshinfo->unstructured;
                    if (bp->ccounts)
                        free (bp->ccounts);
                    for (i = 0; i < meshinfo->unstructured->ncsets; i++)
                    {
                        if (bp->cdata[i])
                            free (bp->cdata[i]);
                        if (bp->ctypes[i])
                            free (bp->ctypes[i]);

                    }
                    for (i = 0; i < meshinfo->unstructured->nvar_points; i++)
                    {
                        if (bp->points[i])
                            free (bp->points[i]);
                    }
                    free (meshinfo->unstructured);
                    break;
                }
            default:
                break;
        }
        free (meshinfo);
    }
}

int common_read_schedule_read (const ADIOS_FILE      * fp,
                               const ADIOS_SELECTION * sel,
                               const char            * varname,
                               int                     from_steps,
                               int                     nsteps,
                               const char            * param, // NCSU ALACRITY-ADIOS
                               void                  * data)

{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        int varid = common_read_find_name (fp->nvars, fp->var_namelist, varname, 0);
        if (varid >= 0) {
            retval = common_read_schedule_read_byid (fp, sel, varid, from_steps, nsteps, param /* NCSU ALACRITY-ADIOS */, data);
        } else {
            retval = adios_errno; // adios_errno was set in common_read_find_name
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_schedule_read()\n");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

// NCSU ALACRITY-ADIOS - Modified to delegate to transform method when called
//   on a transformed variable
int common_read_schedule_read_byid (const ADIOS_FILE      * fp,
        const ADIOS_SELECTION * sel,
        int                     varid,
        int                     from_steps,
        int                     nsteps,
        const char            * param, // NCSU ALACRITY-ADIOS
        void                  * data)

{
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_start ("adios_schedule_read");
#endif
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = err_no_error;
    if (fp) {
        if (varid >=0 && varid < fp->nvars) {
            // NCSU ALACRITY-ADIOS - If the variable is transformed, intercept
            //   the read scheduling and schedule our own reads
            ADIOS_VARINFO *raw_varinfo = common_read_inq_var_raw_byid(fp, varid);        // Get the *raw* varinfo
            ADIOS_TRANSINFO *transinfo = common_read_inq_transinfo(fp, raw_varinfo);    // Get the transform info (i.e. original var info)
            //assert(vi);
            //assert(ti);

            // If this variable is transformed, delegate to the transform
            // method to generate subrequests
            // Else, do the normal thing
            if (transinfo && transinfo->transform_type != adios_transform_none) {
                adios_transform_raw_read_request *subreq;
                adios_transform_pg_read_request *pg_reqgroup;
                adios_transform_read_request *new_reqgroup;

                internals = (struct common_read_internals_struct *) fp->internal_data;

#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_start ("adios_transform_generate_read_requests");
#endif
                // Generate the read request group and append it to the list
                new_reqgroup = adios_transform_generate_read_reqgroup(raw_varinfo, transinfo, fp, sel, from_steps, nsteps, param, data);
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_stop ("adios_transform_generate_read_requests");
#endif

                // Proceed to register the read request and schedule all of its grandchild raw
                // read requests ONLY IF a non-NULL reqgroup was returned (i.e., the user's
                // selection intersected at least one PG).
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_start ("adios_transform_submit_read_requests");
#endif
                if (new_reqgroup) {
                    adios_transform_read_request_append(&internals->transform_reqgroups, new_reqgroup);

                    // Now schedule all of the new subrequests
                    retval = 0;
                    for (pg_reqgroup = new_reqgroup->pg_reqgroups; pg_reqgroup; pg_reqgroup = pg_reqgroup->next) {
                        for (subreq = pg_reqgroup->subreqs; subreq; subreq = subreq->next) {
                            retval |= internals->read_hooks[internals->method].adios_schedule_read_byid_fn(
                                            fp, subreq->raw_sel, varid+internals->group_varid_offset, pg_reqgroup->timestep, 1, subreq->data);
                        }
                    }
                }
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_stop ("adios_transform_submit_read_requests");
#endif
            } else {
                // Old functionality
            internals = (struct common_read_internals_struct *) fp->internal_data;
            retval = internals->read_hooks[internals->method].adios_schedule_read_byid_fn (fp, sel, varid+internals->group_varid_offset, from_steps, nsteps, data);
            }
        } else {
            adios_error (err_invalid_varid, 
                         "Variable ID %d is not valid in adios_schedule_read_byid(). "
                         "Available 0..%d\n", varid, fp->nvars-1);
            retval = err_invalid_varid;
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_schedule_read_byid()\n");
        retval = err_invalid_file_pointer;
    }

#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_stop ("adios_schedule_read");
#endif
    return retval;
}

// NCSU ALACRITY-ADIOS - Modified to delegate to transform method to combine
//  read subrequests to answer original requests
int common_read_perform_reads (const ADIOS_FILE *fp, int blocking)
{
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_start ("adios_perform_reads");
#endif
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_perform_reads_fn (fp, blocking);

        // NCSU ALACRITY-ADIOS - If this was a blocking call, consider all read
        //   request groups completed, and reassemble via the transform method.
        //   Otherwise, do nothing.
        if (blocking) {
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_start ("adios_perform_reads_transform");
#endif
            adios_transform_process_all_reads(&internals->transform_reqgroups);
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_stop ("adios_perform_reads_transform");
#endif
        } else {
            // Do nothing; reads will be performed by check_reads
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_perform_reads()\n");
        retval = err_invalid_file_pointer;
    }
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 2)
    timer_stop ("adios_perform_reads");
#endif
    return retval;
}

int common_read_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;

        // NCSU ALACRITY-ADIOS - Handle those VARCHUNKs that correspond to
        //   subrequests; don't return until we get a completed one
        do {
            // Read some chunk
        retval = internals->read_hooks[internals->method].adios_check_reads_fn (fp, chunk);
            if (!*chunk) break; // If no more chunks are available, stop now

            // Process the chunk through a transform method, if necessary
            adios_transform_process_read_chunk(internals->transform_reqgroups, chunk);
        } while (!*chunk); // Keep reading until we have a chunk to return
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_check_reads()\n");
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
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        int attrid = common_read_find_name (fp->nattrs, fp->attr_namelist, attrname, 1);
        if (attrid > -1) {
            retval = common_read_get_attr_byid (fp, attrid, type, size, data);
        } else {
            retval = adios_errno; // adios_errno was set in common_read_find_name
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_read_get_attr()\n");
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
    
    adios_errno = err_no_error;
    if (fp) {
        if (attrid >= 0 && attrid < fp->nattrs) {
            internals = (struct common_read_internals_struct *) fp->internal_data;
            retval = internals->read_hooks[internals->method].adios_get_attr_byid_fn (fp, attrid+internals->group_attrid_offset, type, size, data);
        } else {
            adios_error (err_invalid_attrid, 
                         "Attribute ID %d is not valid in adios_get_attr_byid(). "
                         "Available 0..%d\n", attrid, fp->nattrs-1);
            retval = err_invalid_attrid;
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_read_get_attr_byid()\n");
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
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->ngroups;
        *group_namelist = internals->group_namelist;
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_get_grouplist()\n");
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
    
    adios_errno = err_no_error;
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
            adios_error (err_invalid_group, "Invalid group ID in adios_group_view()\n");
            retval = err_invalid_group;
        }
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to adios_group_view()\n");
        retval = err_invalid_file_pointer;
    }
    return retval;
}

/* internal function to support version 1 time-dimension reads
   called from adios_read_v1.c and adiosf_read_v1.c 
*/
int common_read_is_var_timed (const ADIOS_FILE *fp, int varid)
{
    struct common_read_internals_struct * internals;
    int retval;
    
    adios_errno = err_no_error;
    if (fp) {
        internals = (struct common_read_internals_struct *) fp->internal_data;
        retval = internals->read_hooks[internals->method].adios_is_var_timed_fn (fp, varid+internals->group_varid_offset);
    } else {
        adios_error (err_invalid_file_pointer, "Null pointer passed as file to common_read_is_var_timed()\n");
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
ADIOS_SELECTION * common_read_selection_boundingbox (int ndim, const uint64_t *start, const uint64_t *count)
{   
    adios_errno = err_no_error;
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_BOUNDINGBOX;
        sel->u.bb.ndim = ndim;
        sel->u.bb.start = (uint64_t *)start;
        sel->u.bb.count = (uint64_t *)count;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for bounding box selection\n");
    }
    return sel;
}


ADIOS_SELECTION * common_read_selection_points (int ndim, uint64_t npoints, const uint64_t *points)
{   
    adios_errno = err_no_error;
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_POINTS;
        sel->u.points.ndim = ndim;
        sel->u.points.npoints = npoints;
        sel->u.points.points = (uint64_t *) points;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for points selection\n");
    }
    return sel;
}

ADIOS_SELECTION * common_read_selection_writeblock (int index)
{   
    adios_errno = err_no_error;
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_WRITEBLOCK;
        sel->u.block.index = index;
        // NCSU ALACRITY-ADIOS: Set the writeblock selection to be a full-PG selection by default
        sel->u.block.is_absolute_index = 0;
        sel->u.block.is_sub_pg_selection = 0;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for writeblock selection\n");
    }
    return sel;
}

ADIOS_SELECTION * common_read_selection_auto (char *hints)
{   
    adios_errno = err_no_error;
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_AUTO;
        sel->u.autosel.hints = hints;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for auto selection\n");
    }
    return sel;
}

void common_read_selection_delete (ADIOS_SELECTION *sel)
{
    free(sel);
}
