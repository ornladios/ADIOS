/*=================================================================
 * Copyright 2009 Oak Ridge National Laboratory
 * $Revision: 1.0 $  $Date: 2009/08/05 12:53:41 $
 * Author: Norbert Podhorszki <pnorbert@ornl.gov>
 *=================================================================*/
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "adios_readutil.h"
#include "adios_read.h"
#include "adios_types.h"
#include "mpidummy.h"

static int verbose = 0;
void adios_readutil_setverbosity(int level) { verbose=level; }

int adios_readutil_open( const char *file, int64_t *fh,
                         int *ngroups, adios_readutil_namelist **groupnames,
                         int *tstart, int *tend, char *msg)
{
    int     status, idx;
    int     mpi_comm_dummy;
    BP_FILE_INFO finfo;
    adios_readutil_namelist *nm, *current;
    int64_t gh;

    /* get file handler */
    status = adios_fopen (fh, file, mpi_comm_dummy);
    if (status != 0) {
        snprintf(msg, 256, "Error opening bp file %s\n", file);
        return 1;
    }
    /*
    if (verbose>1) printf("    open: opened file %s, fh=%lld\n",file,*fh);
    status = adios_gopen (*fh, &gh, "pot0");
    if (verbose>1) printf("    open: opened group pot0, gh=%lld\n",file,gh);
    */

    /* get groups info from file */
    adios_init_fileinfo ( &finfo, 1);
    adios_inq_file (*fh, &finfo);
    *ngroups = finfo.groups_count;
    *tstart = finfo.tidx_start;
    *tend = finfo.tidx_stop;
    *groupnames = NULL;
    for (idx=0; idx<finfo.groups_count; idx++) {
        if (verbose>2) printf("    %s\n",finfo.group_namelist[idx]);
        nm = (adios_readutil_namelist *) malloc(sizeof(adios_readutil_namelist));
        nm->name = strdup(finfo.group_namelist[idx]);
        nm->next = NULL;
        if (*groupnames == NULL) { /* insert as head the first name */
            *groupnames = nm;
        } else { /* add name to end of list */ 
            current->next = nm;
        }
        current = nm;
    }


    adios_free_fileinfo(&finfo);
    return 0;
}


int adios_readutil_groupinfo (int64_t fh, const char *groupname, int64_t *gh,
                              int *nvars, adios_readutil_namelist **varnames,
                              int *nattrs, adios_readutil_namelist **attrnames,
                              char *msg)
{
    BP_GROUP_INFO ginfo; 
    int i, status;
    adios_readutil_namelist *nm, *current;

    status = adios_gopen (fh, gh, groupname);
    /* FIXME: adios_gopen does not return status */
    adios_init_groupinfo( &ginfo, 1);
    adios_inq_group (*gh, &ginfo);

    if (verbose>1) printf("    ginfo: nvars=%d, nattrs=%d\n",ginfo.vars_count,ginfo.attrs_count);
    /* add variable names to a linked list */
    *nvars = ginfo.vars_count;
    *varnames = NULL;
    for (i=0; i<ginfo.vars_count; i++) {
        nm = (adios_readutil_namelist *) malloc(sizeof(adios_readutil_namelist));
        nm->name = strdup(ginfo.var_namelist[i]);
        nm->next = NULL;
        if (*varnames == NULL) { /* insert as head the first name */
            *varnames = nm;
        } else { /* add name to end of list */ 
            current->next = nm;
        }
        current = nm;
    }

    /* add attribute names to a linked list */
    *nattrs = ginfo.attrs_count;
    *attrnames = NULL;
    for (i=0; i<ginfo.attrs_count; i++) {
        nm = (adios_readutil_namelist *) malloc(sizeof(adios_readutil_namelist));
        nm->name = strdup(ginfo.attr_namelist[i]);
        nm->next = NULL;
        if (*attrnames == NULL) { /* insert as head the first name */
            *attrnames = nm;
        } else { /* add name to end of list */ 
            current->next = nm;
        }
        current = nm;
    }

    adios_free_groupinfo(&ginfo);
    return 0;
}

int adios_readutil_getvarinfo( int64_t gh, const char *path,
                               int *ndims, int *dims, int *type, int *timed, char *msg)
{
    int retval;
    retval = adios_inq_var (gh, path, type, ndims, timed, dims);
    if (retval == -1) {
        snprintf(msg, 256, "Error: variable %s was not found\n", path);
    } else if (retval < 0) {
        snprintf(msg, 256, "Error within adios_inq_var(). Error code was %d. Check stderr\n", retval);
    }
    return retval;
}

int adios_readutil_getattrinfo( int64_t gh, const char *path,
                                int *type, int *size, char *msg)
{
    int retval, size_in_bytes;
    retval = adios_inq_attr (gh, path, type, &size_in_bytes);
    if (retval == -1) {
        snprintf(msg, 256, "Error: attribute %s was not found\n", path);
    } else if (retval < 0) {
        snprintf(msg, 256, "Error within adios_inq_attr(). Error code was %d. Check stderr\n", retval);
    } else {
        if (*type == adios_string)
            *size = size_in_bytes;
        else
            *size = 1;
    }
    return retval;
}

void adios_readutil_fclose(int64_t fh) { adios_fclose(fh); }

void adios_readutil_gclose(int64_t gh) { adios_gclose(gh); }


int adios_readutil_fopen( const char *file, int64_t *fh, char *msg)
{
    int     status;
    int     mpi_comm_dummy;

    status = adios_fopen (fh, file, mpi_comm_dummy);
    if (status != 0) {
        //*msg = (char *) malloc(256*sizeof(char));
        snprintf(msg, 256, "Error opening bp file %s\n", file);
        return 1;
    }
    return 0;
}


int adios_readutil_gopen_byname( const char *groupname, int64_t fh, int64_t *gh_p, char *msg)
{
    int     status;
    BP_FILE_INFO finfo;
    int     retval, idx;

    /* get groups info from file */
    adios_init_fileinfo ( &finfo, 1);
    adios_inq_file (fh, &finfo);
    
    if ( strlen(groupname)==0 || 
         !strcmp(groupname,"/") ||
         !strcmp(groupname," ")) {

        /* special case: group is empty string or / or ' ': read in first group */
        idx = 0;

    } else {
        /* search for groupname in the available groups */
        for (idx=0; idx<finfo.groups_count; idx++) {
            if ( !strcmp(groupname,finfo.group_namelist[idx]) ) {
                break;
            }
        }
    }
    if (idx >= finfo.groups_count) { 
        //*msg = (char *) malloc(256*sizeof(char));
        snprintf(msg, 256, "Error group \"%s\" not found.\n", groupname);
        retval = 1;
    } else {
        status = adios_gopen (fh, gh_p, finfo.group_namelist[idx]);
        /* FIXME: gopen does not return status */
        retval = 0;
    }
    adios_free_fileinfo(&finfo);
    return retval;
}


int adios_readutil_gopen_byindex( int32_t groupidx, int64_t fh, int64_t *gh_p, char *msg)
{
    int     status;
    int     gidx = (int) groupidx;
    BP_FILE_INFO finfo;
    int     retval;

    /* get groups info from file */
    adios_init_fileinfo ( &finfo, 1);
    adios_inq_file (fh, &finfo);

    /* group index in Matlab starts from 1, in C it starts from 0 */
    if (groupidx > 0 && finfo.groups_count >= groupidx) { /* index is valid */
        status = adios_gopen (fh, gh_p, finfo.group_namelist[groupidx-1]);
        /* FIXME: gopen does not return status */
        retval = 0;
    } else {
        //*msg = (char *) malloc(256*sizeof(char));
        snprintf(msg, 256, "Invalid group index %d. Must be 1..%d\n", groupidx, finfo.groups_count);
        retval = 1;
    }
    adios_free_fileinfo(&finfo);
    return retval;
}


const char *adios_readutil_type_to_string(type) { return bp_type_to_string(type); }


#if 1
int adios_readutil_getdatainfo( int64_t gh, const char *path,
                                int *ndims, int *dims, int *type, char *msg, int *isvar)
{
    BP_GROUP_INFO ginfo; 
    int i, retval, hastimesteps, size_in_bytes;
    int found=0;

    adios_init_groupinfo( &ginfo, 1);
    // get variable info from group
    adios_inq_group (gh, &ginfo);

    if (verbose>1) printf("    ginfo # of vars=%d\n",ginfo.vars_count);
    if (verbose>1) printf("    ginfo # of attrs=%d\n",ginfo.attrs_count);
    /* Search among variables */
    for (i=0; i<ginfo.vars_count; i++) {
        if (!strcmp(path, ginfo.var_namelist[i])) {
            // path matches a variable in the group: get info about it
            *isvar = 1;
            retval = adios_inq_var (gh, path, type, ndims, &hastimesteps, dims);
            if (retval == -1) {
                //*msg = (char *) malloc(256*sizeof(char));
                snprintf(msg, 256, "Error: getting info on variable %s\n", path, retval);
            } else {
                found = 1;
                break;
            }
        }
    }

    if (!found) {
        /* Search among attributes */
        *isvar = 0;
        if (verbose>1) printf("    path %s is not variable, search in attrs\n",path);
        for (i=0; i<ginfo.attrs_count; i++) {
            if (!strcmp(path, ginfo.attr_namelist[i])) {
                // path matches an attribute in the group: get info about it
                retval = adios_inq_attr (gh, path, type, &size_in_bytes);
                if (retval < 0) {
                    //*msg = (char *) malloc(256*sizeof(char));
                    snprintf(msg, 256, "Error: getting info on attribute %s\n", path, retval);
                } else {
                    /* It's an attribute */
                    *ndims = 1;
                    if (*type == adios_string)
                        dims[0] = size_in_bytes;
                    else
                        dims[0] = 1;
                    found=1;
                    break;
                }
            }
        }
    }

    if (!found) {
        if (verbose>1) printf("   path %s is not var not attr\n",path);
        //*msg = (char *) malloc(256*sizeof(char));
        snprintf(msg, 256, "Error: neither variable nor attribute was found with path %s\n", path);
        retval = -1;
        if (verbose>1) printf("   set retval to -1\n");
    }

    if (verbose>1) printf("   free groupinfo...\n");
    adios_free_groupinfo(&ginfo);
    if (verbose>1) printf("   freed\n");
    return retval;
}

#else
int adios_readutil_getdatainfo( int64_t gh, const char *path,
                                int *ndims, int *dims, int *type, char *msg, int *isvar)
{
    int retval, hastimesteps, size_in_bytes;
    retval = adios_inq_var (gh, path, type, ndims, &hastimesteps, dims);
    if (retval == -1) {
        /* variable not found, search for attrs */
        retval = adios_inq_attr (gh, path, type, &size_in_bytes);
        if (retval < 0) {
            //*msg = (char *) malloc(256*sizeof(char));
            snprintf(msg, 256, "Error: neither variable nor attribute was found with path %s\n", path, retval);
        } else {
            /* It's an attribute */
            *ndims = 1;
            if (*type == adios_string)
                dims[0] = size_in_bytes;
            else
                dims[0] = 1;
            *isvar = 0;
        }
    } else if (retval < 0) {
        //*msg = (char *) malloc(256*sizeof(char));
        snprintf(msg, 256, "Error within adios_inq_var(). Error code was %d. Check stderr\n", retval);
    } else {
        /* It's a variable */
        *isvar = 1;
    }
    return retval;
}
#endif

int32_t adios_readutil_calctimestep(int64_t fh, int32_t timestep)
{
    BP_FILE_INFO finfo;
    int32_t retval = timestep;

    /* negative timestep value has to be interpreted */
    if (timestep < 0) {
        /* get groups info from file */
        adios_init_fileinfo ( &finfo, 1);
        adios_inq_file (fh, &finfo);

        /* -1 is the last timestep, -2 is the previous */
        retval = finfo.tidx_stop + timestep + 1;

        adios_free_fileinfo(&finfo);
    }
    
    return retval;
}

int adios_readutil_readdata(int64_t gh, const char *path, int isvar, 
                            const int *offset, const int *count, int timestep,
                            void *data, char *msg)
{
    int64_t bytes_read;
    int retval;


    if (isvar) { /* Read variable */
        bytes_read = adios_get_var (gh, path, data, offset, count, timestep);
    } else {
        bytes_read = adios_get_attr (gh, path, data);
    }
    if (bytes_read < 0) {
        //*msg = (char *) malloc(256*sizeof(char));
        snprintf(msg, 256, "Error reading %s %s\n", (isvar ? "variable" : "attribute"), path);
        retval = (int)bytes_read;
    } else {
        if (verbose>1) printf("    adios api has read in %lld bytes\n",bytes_read);
        retval = 0;
    }
    return retval;
}


int adios_readutil_getdirectattributes(int64_t gh, const char* path, adios_readutil_namelist **attrs)
{
    BP_GROUP_INFO ginfo; /* reused for each group */
    int vr, len, nmatched;
    int *idxs;   /* array of attr indices to be stored at the search phase */
    adios_readutil_namelist *nm, *current;

    adios_init_groupinfo( &ginfo, 1);
    // get variable info from group
    adios_inq_group (gh, &ginfo);
    
    len = strlen(path); 
    nmatched = 0;
    *attrs = NULL;
    if (verbose>1) printf("    Check %d attrs if they match path %s\n",ginfo.attrs_count,path);
    for (vr=0; vr<ginfo.attrs_count; vr++) {
        if (verbose>2) printf("        Check path %s\n",ginfo.attr_namelist[vr]);
        if (!strncmp(path, ginfo.attr_namelist[vr], len)) {
            /* head matches */
            if (ginfo.attr_namelist[vr][len] == '/') {
                /* next char is / */
                if (index(&(ginfo.attr_namelist[vr][len+1]),'/') == NULL ) {
                    /* no more / in path */
                    if (verbose>1) 
                        printf("    Attribute %s is child of path %s\n", 
                                ginfo.attr_namelist[vr], path);
                    /* add to attrlist */
                    nm = (adios_readutil_namelist *) malloc(sizeof(adios_readutil_namelist));
                    nm->name = strdup(ginfo.attr_namelist[vr]);
                    nm->next = NULL;
                    if (*attrs == NULL) { /* insert as head the first name */
                        *attrs = nm;
                    } else { /* add name to end of list */ 
                        current->next = nm;
                    }
                    current = nm;
                    nmatched++;
                }
            }
        }
    }

    adios_free_groupinfo(&ginfo);
    return nmatched;
}

void adios_readutil_freenamelist(adios_readutil_namelist *namelist)
{
    adios_readutil_namelist *nm = namelist;
    adios_readutil_namelist *old;

    while (nm != NULL) {
        if (nm->name != NULL) free(nm->name);
        old = nm;
        nm = nm->next;
        free(old);
    }
}

