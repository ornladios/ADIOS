/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */


/****************************/
/* Read method for BP files */
/****************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "adios.h"
#include "bp_utils.h"
#include "bp_types.h"
#include "adios_types.h"
#include "adios_read.h"
#include "adios_read_hooks.h"
#include "adios_error.h"
#include "futils.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

MPI_File * get_BP_file_handle(struct BP_file_handle * l, uint32_t file_index)
{
    if (!l)
        return 0;

    while (l)
    {
        if (l->file_index == file_index)
            return &l->fh;

        l = l->next;
    }

    return 0;
}

void add_BP_file_handle (struct BP_file_handle ** l, struct BP_file_handle * n)
{
    if (!n)
        return;

    n->next = *l;
    *l = n;
}


void close_all_BP_files (struct BP_file_handle * l)
{
    struct BP_file_handle * n;

    while (l)
    {
        n = l->next;

        MPI_File_close (&l->fh);
        free (l);

        l = n;
    }
}

/* Return 0: if file is little endian, 1 if file is big endian 
 * We know if it is different from the current system, so here
 * we determine the current endianness and report accordingly.
 */
static int adios_read_bp_get_endianness( uint32_t change_endianness )
{
   int LE = 0;
   int BE = !LE;
   int i = 1;
   char *p = (char *) &i;
   int current_endianness;
   if (p[0] == 1) // Lowest address contains the least significant byte
       current_endianness = LE;
   else
       current_endianness = BE;
    if (change_endianness == adios_flag_yes)
        return !current_endianness;
    else
        return current_endianness;
}

int adios_read_bp_init (MPI_Comm comm) { return 0; }
int adios_read_bp_finalize () { return 0; }

ADIOS_FILE * adios_read_bp_fopen (const char * fname, MPI_Comm comm)
{
    int i, rank;    
    struct BP_FILE * fh;
    ADIOS_FILE * fp;
    uint64_t header_size;

    adios_errno = 0;
    fh = (struct BP_FILE *) malloc (sizeof (struct BP_FILE));
    if (!fh) {
        adios_error ( err_no_memory, "Cannot allocate memory for file info.");
        return NULL;
    }

    fh->fname = (fname ? strdup (fname) : 0L);
    fh->sfh = 0;
    fh->comm = comm;
    fh->gvar_h = 0;
    fh->pgs_root = 0;
    fh->vars_root = 0;
    fh->attrs_root = 0;
    fh->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    if (!fh->b) {
        adios_error ( err_no_memory, "Cannot allocate memory for file info.");
        return NULL;
    }
    fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    if (!fp) {
        adios_error ( err_no_memory, "Cannot allocate memory for file info.");
        return NULL;
    }

    adios_buffer_struct_init (fh->b);
    MPI_Comm_rank (comm, &rank);
    if (bp_read_open (fname, comm, fh))
        return NULL;

    /* Only rank=0 process reads the footer and it broadcasts to all other processes */
    if ( rank == 0 ) {
        if (bp_read_minifooter (fh))
            return NULL;
    }
    MPI_Bcast (&fh->mfooter, sizeof(struct bp_minifooter),MPI_BYTE, 0, comm);
    
    header_size = fh->mfooter.file_size-fh->mfooter.pgs_index_offset;

    if ( rank != 0) {
        if (!fh->b->buff) {
            bp_alloc_aligned (fh->b, header_size);
            if (!fh->b->buff)
                return NULL;
            memset (fh->b->buff, 0, header_size);
            fh->b->offset = 0;
        }
    }
    MPI_Barrier(comm);
    MPI_Bcast (fh->b->buff, fh->mfooter.file_size-fh->mfooter.pgs_index_offset, MPI_BYTE, 0 , comm);

    /* Everyone parses the index on its own */
    bp_parse_pgs (fh);
    bp_parse_vars (fh);
    bp_parse_attrs (fh);

    /* fill out ADIOS_FILE struct */
    fp->fh = (uint64_t) fh;
    fp->groups_count = fh->gvar_h->group_count;
    fp->vars_count = fh->mfooter.vars_count;
    fp->attrs_count = fh->mfooter.attrs_count;
    fp->tidx_start = fh->tidx_start;
    fp->ntimesteps = fh->tidx_stop - fh->tidx_start + 1;
    fp->file_size = fh->mfooter.file_size;
    fp->version = fh->mfooter.version;
    fp->endianness = adios_read_bp_get_endianness(fh->mfooter.change_endianness);
    alloc_namelist (&fp->group_namelist,fp->groups_count); 
    for (i=0;i<fp->groups_count;i++) {
        if (!fp->group_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate buffer for %d strings in adios_fopen()", fp->groups_count);
            adios_read_bp_fclose(fp);
            return NULL;
        }
        else  {
            strcpy(fp->group_namelist[i],fh->gvar_h->namelist[i]);
        }
    }
    return fp;
}

/* This function can be called if user places 
   the wrong sequences of dims for a var 
*/   
void adios_read_bp_reset_dimension_order (ADIOS_FILE *fp, int is_fortran)
{
    struct BP_FILE * fh = (struct BP_FILE *)(fp->fh);
    struct bp_index_pg_struct_v1 ** root = &(fh->pgs_root);
    struct bp_minifooter * mh = &(fh->mfooter);
    uint64_t i;

    for (i = 0; i < mh->pgs_count; i++) {
        is_fortran ? ((*root)->adios_host_language_fortran = adios_flag_yes) 
               : ((*root)->adios_host_language_fortran = adios_flag_no);
        root = &(*root)->next;
    }
}

int adios_read_bp_fclose (ADIOS_FILE *fp) 
{
    struct BP_FILE * fh = (struct BP_FILE *) fp->fh;
    struct BP_GROUP_VAR * gh = fh->gvar_h;
    struct BP_GROUP_ATTR * ah = fh->gattr_h;
    struct adios_index_var_struct_v1 * vars_root = fh->vars_root, *vr;
    struct adios_index_attribute_struct_v1 * attrs_root = fh->attrs_root, *ar;
    struct bp_index_pg_struct_v1 * pgs_root = fh->pgs_root, *pr;
    int i,j;
    MPI_File mpi_fh = fh->mpi_fh;

    adios_errno = 0;
    if (fh->mpi_fh) 
        MPI_File_close (&mpi_fh);

    if (fh->sfh)
        close_all_BP_files (fh->sfh);

    if (fh->b) {
        adios_posix_close_internal (fh->b);
        free(fh->b);
    }

    /* Free variable structures */
    /* alloc in bp_utils.c: bp_parse_vars() */
    while (vars_root) {
        vr = vars_root;
        vars_root = vars_root->next;
        for (j = 0; j < vr->characteristics_count; j++) {
            // alloc in bp_utils.c:bp_parse_characteristics() <- bp_get_characteristics_data()
            if (vr->characteristics[j].dims.dims)
                free (vr->characteristics[j].dims.dims);
            if (vr->characteristics[j].value)
                free (vr->characteristics[j].value);
            // NCSU - Clearing up statistics
            if (vr->characteristics[j].stats)
            {
                uint8_t k = 0, idx = 0;
                uint8_t i = 0, count = adios_get_stat_set_count(vr->type);

                while (vr->characteristics[j].bitmap >> k)
                {
                    if ((vr->characteristics[j].bitmap >> k) & 1)
                    {
                        for (i = 0; i < count; i ++)
                        {
                            if (k == adios_statistic_hist)
                            {
                                struct adios_index_characteristics_hist_struct * hist = (struct adios_index_characteristics_hist_struct *) (vr->characteristics [j].stats[i][idx].data);
                                free (hist->breaks);
                                free (hist->frequencies);
                                free (hist);
                            }
                            else
                            free (vr->characteristics[j].stats [i][idx].data);
                        }
                        idx ++;
                    }
                    k ++;
                }

                for (i = 0; i < count; i ++)
                    free (vr->characteristics[j].stats [i]);

                free (vr->characteristics[j].stats);
                vr->characteristics[j].stats = 0;
            }
        }
        if (vr->characteristics) 
            free (vr->characteristics);
        if (vr->group_name) 
            free (vr->group_name);
        if (vr->var_name) 
            free (vr->var_name);
        if (vr->var_path) 
            free (vr->var_path);
        free(vr);
    }

    /* Free attributes structures */
    /* alloc in bp_utils.c bp_parse_attrs() */
    while (attrs_root) {
        ar = attrs_root;
        attrs_root = attrs_root->next;
        for (j = 0; j < ar->characteristics_count; j++) {
            if (ar->characteristics[j].value)
                free (ar->characteristics[j].value);
        }
        if (ar->characteristics) 
            free (ar->characteristics);
        if (ar->group_name) 
            free (ar->group_name);
        if (ar->attr_name) 
            free (ar->attr_name);
        if (ar->attr_path) 
            free (ar->attr_path);
        free(ar);
    }


    /* Free process group structures */
    /* alloc in bp_utils.c bp_parse_pgs() first loop */
    //printf ("pgs: %d\n", fh->mfooter.pgs_count);
    while (pgs_root) {
        pr = pgs_root;
        pgs_root = pgs_root->next;
        //printf("%d\tpg pid=%d addr=%x next=%x\n",i, pr->process_id, pr, pr->next);
        if (pr->group_name)
            free(pr->group_name);
        if (pr->time_index_name)
            free(pr->time_index_name);
        free(pr);
    }

    /* Free variable structures in BP_GROUP_VAR */
    if (gh) {
        for (j=0;j<2;j++) { 
            for (i=0;i<gh->group_count;i++) {
                if (gh->time_index[j][i])
                    free(gh->time_index[j][i]);
            }
            if (gh->time_index[j])
                free(gh->time_index[j]);
        }
        free (gh->time_index);
    
        for (i=0;i<gh->group_count;i++) { 
            if (gh->namelist[i])
                free(gh->namelist[i]);
        }
        if (gh->namelist)
            free (gh->namelist);

        for (i=0;i<fh->mfooter.vars_count;i++) {
            if (gh->var_namelist[i])
                free(gh->var_namelist[i]);
            if (gh->var_offsets[i]) 
                free(gh->var_offsets[i]);
        }
        if (gh->var_namelist)
            free (gh->var_namelist);

        if (gh->var_offsets) 
            free(gh->var_offsets);

        if (gh->var_counts_per_group)
            free(gh->var_counts_per_group);

        if (gh->pg_offsets)
            free (gh->pg_offsets);

        free (gh);
    }

    /* Free attribute structures in BP_GROUP_ATTR */
    if (ah) {
        for (i = 0; i < fh->mfooter.attrs_count; i++) {
            if (ah->attr_offsets[i]) 
                free(ah->attr_offsets[i]);
            if (ah->attr_namelist[i]) 
                free(ah->attr_namelist[i]);
        }
        if (ah->attr_offsets)
            free(ah->attr_offsets);
        if (ah->attr_namelist)
            free(ah->attr_namelist);
        if (ah->attr_counts_per_group) 
            free(ah->attr_counts_per_group);

        free(ah);
    }

    if (fh->fname)
        free (fh->fname);
        
    if (fh)
        free (fh);    

    free_namelist ((fp->group_namelist),fp->groups_count);
    free(fp);
    return 0;
}


ADIOS_GROUP * adios_read_bp_gopen (ADIOS_FILE *fp, const char * grpname)
{
    struct BP_FILE * fh = (struct BP_FILE *) fp->fh;
    int grpid; 

    adios_errno = 0;
    for (grpid=0;grpid<(fh->gvar_h->group_count);grpid++) {
        if (!strcmp(fh->gvar_h->namelist[grpid], grpname))
            break; 
    }
    if (grpid >= fh->gvar_h->group_count) {
        adios_error ( err_invalid_group, "Invalid group name %s", grpname);
        return NULL;
    }
    return adios_read_bp_gopen_byid(fp, grpid);
}

ADIOS_GROUP * adios_read_bp_gopen_byid (ADIOS_FILE *fp, int grpid)
{
    struct BP_FILE * fh = (struct BP_FILE *) fp->fh;
    struct BP_GROUP * gh;
    ADIOS_GROUP * gp;
    int i, offset;

    adios_errno = 0;
    if (grpid < 0 || grpid >= fh->gvar_h->group_count) {
        adios_error ( err_invalid_group, "Invalid group index %d", grpid);
        return NULL;
    }

    gh = (struct BP_GROUP *) malloc(sizeof(struct BP_GROUP));
    if (!gh) {
        adios_error ( err_no_memory, "Could not allocate memory for group info");
        return NULL;
    }

    gp = (ADIOS_GROUP *) malloc(sizeof(ADIOS_GROUP));
    if (!gp) {
        adios_error ( err_no_memory, "Could not allocate memory for group info");
        free(gh);
        return NULL;
    }

    /* set offset index of variables (which is a long list of all vars in all groups) in this group */
    offset = 0;
    for (i=0; i<grpid; i++)
        offset += fh->gvar_h->var_counts_per_group[i];

    /* gh->vars_root will point to the list of vars in this group */
    gh->vars_root = fh->vars_root; 
    for (i=0; i<offset; i++)
        gh->vars_root = gh->vars_root->next;

    gh->group_id = grpid;
    gh->vars_offset = offset;
    gh->vars_count = fh->gvar_h->var_counts_per_group[grpid];

    /* set offset of attributes in this group */
    offset = 0;
    for(i=0;i<grpid;i++)
        offset += fh->gattr_h->attr_counts_per_group[i];

    /* gh->attrs_root will point to the list of vars in this group */
    gh->attrs_root = fh->attrs_root; 
    for (i=0; i<offset; i++)
        gh->attrs_root = gh->attrs_root->next;

    gh->attrs_offset = offset;
    gh->attrs_count = fh->gattr_h->attr_counts_per_group[grpid];

    gh->fh = fh; 

    /* fill out ADIOS_GROUP struct */
    gp->grpid = grpid;
    gp->gh = (uint64_t) gh;
    gp->fp = fp;
    gp->vars_count = gh->vars_count;
    gp->attrs_count = gh->attrs_count;
    
    offset = gh->vars_offset;
    alloc_namelist (&(gp->var_namelist), gp->vars_count);
    for (i=0;i<gp->vars_count;i++) {
        if (!gp->var_namelist[i]) { 
            adios_error (err_no_memory, "Could not allocate buffer for %d strings in adios_gopen()", gp->vars_count);
            adios_read_bp_gclose(gp);
            return NULL;
        }
        else
            strcpy(gp->var_namelist[i], gh->fh->gvar_h->var_namelist[i+offset]);
    }

    offset = gh->attrs_offset;
    alloc_namelist (&(gp->attr_namelist), gp->attrs_count);
    for (i=0;i<gp->attrs_count;i++) {
        if (!gp->attr_namelist[i]) {
            adios_error (err_no_memory, "Could not allocate buffer for %d strings in adios_gopen()", gp->vars_count);
            adios_read_bp_gclose(gp);
            return NULL;
        }
        else {
            strcpy(gp->attr_namelist[i], gh->fh->gattr_h->attr_namelist[i+offset]);
        }
    }

    return gp;
}
                   
int adios_read_bp_gclose (ADIOS_GROUP *gp)
{
    struct BP_GROUP * gh = (struct BP_GROUP *) gp->gh;

    adios_errno = 0;
    if (!gh) {
        adios_error (err_invalid_group_struct, "group handle is NULL!");
        return  err_invalid_group_struct;
    }
    else
        free (gh);

    free_namelist ((gp->var_namelist),gp->vars_count);
    free_namelist ((gp->attr_namelist),gp->attrs_count);
    free(gp);
    return 0;
}



int adios_read_bp_get_attr (ADIOS_GROUP * gp, const char * attrname, enum ADIOS_DATATYPES * type,
                    int * size, void ** data)
{
    // Find the attribute: full path is stored with a starting / 
    // Like in HDF5, we need to match names given with or without the starting /
    // startpos is 0 or 1 to indicate if the argument has starting / or not
    int attrid;
    int vstartpos = 0, fstartpos = 0; 
    struct BP_GROUP * gh = (struct BP_GROUP *)gp->gh;
    int offset;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_get_attr()");
        return adios_errno;
    }
    if (!attrname) {
        adios_error (err_invalid_attrname, "Null pointer passed as attribute name to adios_get_attr()!");
        return adios_errno;
    }

    offset = gh->attrs_offset;

    if (attrname[0] == '/') 
        vstartpos = 1;
    for (attrid=0; attrid<(gp->attrs_count);attrid++) {
        //if (gp->attr_namelist[attrid][0] == '/') 
        if (gh->fh->gattr_h->attr_namelist[attrid+offset][0] == '/') 
            fstartpos = 1;
        //if (!strcmp(gp->attr_namelist[attrid]+fstartpos, attrname+vstartpos))
        if (!strcmp(gh->fh->gattr_h->attr_namelist[attrid+offset]+fstartpos, attrname+vstartpos))
            break; 
    }
    if (attrid >= gp->attrs_count) {
        adios_error ( err_invalid_attrname, "Invalid attribute name %s", attrname);
        return adios_errno;
    }

    return adios_read_bp_get_attr_byid(gp, attrid, type, size, data);
}

int adios_read_bp_get_attr_byid (ADIOS_GROUP * gp, int attrid, 
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    int    i, offset, count;
    struct BP_GROUP * gh;
    struct BP_FILE * fh;
    struct adios_index_attribute_struct_v1 * attr_root;
    struct adios_index_var_struct_v1 * var_root;
    int    file_is_fortran;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_get_attr()");
        return adios_errno;
    }
    gh = (struct BP_GROUP *) gp->gh;
    if (!gh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return adios_errno;
    }
    fh = gh->fh;
    if (!fh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return adios_errno;
    }
    if (attrid < 0 || attrid >= gh->attrs_count) {
        adios_error (err_invalid_attrid, "Invalid attribute id %d (allowed 0..%d)", attrid, gh->attrs_count);
        return adios_errno;
    }

    attr_root = gh->attrs_root; /* need to traverse the attribute list of the group */
    for (i = 0; i < attrid && attr_root; i++)
        attr_root = attr_root->next;
    if (i != attrid) {
        adios_error (err_corrupted_attribute, "Attribute id=%d is valid but was not found in internal data structures!",attrid);
        return adios_errno; 
    }

    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);

    if (attr_root->characteristics[0].value) {
        /* Attribute has its own value */
        *size = bp_get_type_size (attr_root->type, attr_root->characteristics[0].value);
        *type = attr_root->type;
        *data = (void *) malloc (*size);  
        if (*data)
            memcpy(*data, attr_root->characteristics[0].value, *size);
    }
    else if (attr_root->characteristics[0].var_id) {
        /* Attribute is a reference to a variable */
        /* FIXME: var ids are not unique in BP. If a group of variables are written several
           times under different path using adios_set_path(), the id of a variable is always
           the same (should be different). As a temporary fix, we look first for a matching
           id plus path between an attribute and a variable. If not found, then we look for
           a match on the ids only.*/
        var_root = gh->vars_root; 
        while (var_root) {
            if (var_root->id == attr_root->characteristics[0].var_id && 
                !strcmp(var_root->var_path, attr_root->attr_path))
                break;
            var_root = var_root->next;
        }
        if (!var_root) {
            var_root = gh->vars_root; 
            while (var_root) {
                if (var_root->id == attr_root->characteristics[0].var_id)
                    break;
                var_root = var_root->next;
            }
        }

        if (!var_root) {
            adios_error (err_invalid_attribute_reference, 
                   "Attribute %s/%s in group %s is a reference to variable ID %d, which is not found", 
                   attr_root->attr_path, attr_root->attr_name, attr_root->group_name,
                   attr_root->characteristics[0].var_id);
            return adios_errno;
        }

        /* default values in case of error */
        *data = NULL;
        *size = 0;
        *type = attr_root->type;

        /* FIXME: variable and attribute type may not match, then a conversion is needed. */
        /* Cases:
                1. attr has no type, var is byte array     ==> string
                2. attr has no type, var is not byte array ==> var type
                3. attr is string, var is byte array       ==> string
                4. attr type == var type                   ==> var type 
                5. attr type != var type                   ==> attr type and conversion needed 
        */
        /* Error check: attr cannot reference an array in general */
        if (var_root->characteristics[0].dims.count > 0) {
            if ( (var_root->type == adios_byte || var_root->type == adios_unsigned_byte) &&
                 (attr_root->type == adios_unknown || attr_root->type == adios_string) &&
                 (var_root->characteristics[0].dims.count == 1)) {
                 ; // this conversions are allowed
            } else {
                adios_error (err_invalid_attribute_reference, 
                    "Attribute %s/%s in group %s, typeid=%d is a reference to an %d-dimensional array variable "
                    "%s/%s of type %s, which is not supported in ADIOS",
                    attr_root->attr_path, attr_root->attr_name, attr_root->group_name, attr_root->type,
                    var_root->characteristics[0].dims.count,
                    var_root->var_path, var_root->var_name, common_read_type_to_string(var_root->type));
                return adios_errno;
            }
        }

        if ( (attr_root->type == adios_unknown || attr_root->type == adios_string) &&
             (var_root->type == adios_byte || var_root->type == adios_unsigned_byte) &&
             (var_root->characteristics[0].dims.count == 1) ) {
            /* 1D byte arrays are converted to string */
            /* 1. read in variable */
            char varname[512];
            char *tmpdata;
            uint64_t start, count;
            int status;
            start = 0; 
            count = var_root->characteristics[0].dims.dims[0];
            snprintf(varname, 512, "%s/%s", var_root->var_path, var_root->var_name);
            tmpdata = (char *) malloc (count+1);
            if (tmpdata == NULL) {
                adios_error (err_no_memory, 
                      "Cannot allocate memory of %lld bytes for reading in data for attribute %s/%s of group %s.",
                      count, attr_root->attr_path, attr_root->attr_name, attr_root->group_name);
                return adios_errno;
            }

            status = adios_read_bp_read_var (gp, varname, &start, &count, tmpdata);
            
            if (status < 0) {
                char *msg = strdup(adios_get_last_errmsg());
                adios_error ((enum ADIOS_ERRCODES) status, 
                      "Cannot read data of variable %s/%s for attribute %s/%s of group %s: %s",
                      var_root->var_path, var_root->var_name, 
                      attr_root->attr_path, attr_root->attr_name, attr_root->group_name,
                      msg);
                free(tmpdata);
                free(msg);
                return status;
            }

            *type = adios_string;
            if (file_is_fortran) {
                /* Fortran byte array to C string */
                *data = futils_fstr_to_cstr( tmpdata, (int)count); /* FIXME: supports only 2GB strings... */
                *size = strlen( (char *)data );
                free(tmpdata);
            } else {
                /* C array to C string */
                tmpdata[count] = '\0';
                *size = count+1;
                *data = tmpdata;
            }
        } else {
            /* other types are inherited */
            *type = var_root->type;
            *size = bp_get_type_size (var_root->type, var_root->characteristics[0].value);
            *data = (void *) malloc (*size);  
            if (*data)
                memcpy(*data, var_root->characteristics[0].value, *size);
        }
    }

    return 0;
}


/* Reverse the order in an array in place.
   use swapping from Fortran/column-major order to ADIOS-read-api/C/row-major order and back
*/
static void swap_order(int n, uint64_t *array, int *timedim)
{
    int i;
    uint64_t tmp;
    for (i=0; i<n/2; i++) {
        tmp = array[i];
        array[i] = array[n-1-i];
        array[n-1-i] = tmp;
    }
    if (*timedim > -1)
        *timedim = (n-1) - *timedim; // swap the time dimension too
}

/* Look up variable id based on variable name.
   Return index 0..gp->vars_count-1 if found, -1 otherwise
*/
static int adios_read_bp_find_var(ADIOS_GROUP *gp, const char *varname)
{
    // Find the variable: full path is stored with a starting / 
    // Like in HDF5, we need to match names given with or without the starting /
    // startpos is 0 or 1 to indicate if the argument has starting / or not
    int varid;
    int vstartpos = 0, fstartpos; 
    struct BP_GROUP * gh = (struct BP_GROUP *)gp->gh;
    int offset;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group");
        return -1;
    }
    if (!varname) {
        adios_error (err_invalid_varname, "Null pointer passed as variable name!");
        return -1;
    }

    /* Search in gp->fh->gvar_h->var_namelist instead of gp->var_namelist, because that
       one comes back from the user. One idiot (who writes this comment) sorted the
       list in place after gopen and before inq_var.
    */
    offset = gh->vars_offset;

    if (varname[0] == '/') 
        vstartpos = 1;
    for (varid=0; varid<(gp->vars_count);varid++) {
        fstartpos = 0;
        /* if (gp->var_namelist[varid][0] == '/') */
        fstartpos = 0;
        if (gh->fh->gvar_h->var_namelist[varid+offset][0] == '/')
            fstartpos = 1;
        /*if (!strcmp(gp->var_namelist[varid]+fstartpos, varname+vstartpos))*/
        if (!strcmp(gh->fh->gvar_h->var_namelist[varid+offset]+fstartpos, varname+vstartpos))
            break; 
    }
    if (varid >= gp->vars_count) {
        adios_error (err_invalid_varname, "Invalid variable name %s", varname);
        return -1;
    }
    return varid;
}

// NCSU - For custom memory allocation 
#define CALLOC(var, num, sz, comment)\
{\
    var = calloc (num, sz); \
    if (!var)    {\
        adios_error_at_line (err_no_memory, __FILE__, __LINE__, "Could not allocate memory for ", comment, " in common_read_get_characteristics"); \
        return; \
    }\
}

#define MALLOC(var,sz,comment)\
{\
    var = malloc (sz); \
    if (!var)    {\
        adios_error_at_line (err_no_memory, __FILE__, __LINE__, "Could not allocate memory for ", comment, " in common_read_get_characteristics"); \
        return; \
    }\
}\

// NCSU - Reading the statistics
/** Get value and statistics, allocate space for them too */
static void adios_read_bp_get_characteristics (struct adios_index_var_struct_v1 * var_root, ADIOS_VARINFO *vi)
{
    int i, j, c, count = 1;
    int size, sum_size, sum_type;

    vi->value = NULL;

    vi->gmin = vi->gmax = NULL;
    vi->gavg = NULL;
    vi->mins = vi->maxs = NULL;
    vi->avgs = NULL;
    vi->gstd_dev = NULL;
    vi->std_devs = NULL;
    vi->hist = NULL;

    // set value for scalars
    if (var_root->characteristics [0].value) {
        size = bp_get_type_size(var_root->type, var_root->characteristics [0].value);
        vi->value = (void *) malloc (size);
        
        if (vi->value)
           memcpy(vi->value, var_root->characteristics [0].value, size);
        else {
            adios_error_at_line (err_no_memory, __FILE__, __LINE__, "Could not allocate memory for value in common_read_get_characteristics");
            return;
        }
    } else {
        vi->value = NULL;
    }

    int npgs = var_root->characteristics_count, timestep, ntimes = -1;
    uint64_t gcnt = 0, * cnts;

    double *gsum = NULL, *gsum_square = NULL;
    double **sums = NULL, **sum_squares = NULL;

    int16_t map[32];
    memset (map, -1, sizeof(map));

    // Bitmap shows which statistical information has been calculated
    i = j = 0;
    while (var_root->characteristics[0].bitmap >> j)
    {
        if ((var_root->characteristics[0].bitmap >> j) & 1)
            map [j] = i ++;
        j ++;
     }

    if(vi->timedim >= 0)
    {    
        ntimes = vi->dims[0];

        if (map[adios_statistic_min] != -1)
        {
            MALLOC(vi->mins, ntimes * sizeof(void *), "minimum per timestep");
            for (i = 0; i < ntimes; i++)
                vi->mins[i] = 0;
        }

        if (map[adios_statistic_max] != -1)
        {
            MALLOC(vi->maxs, ntimes * sizeof(void *), "maximum per timestep");
            for (i = 0; i < ntimes; i++)
                vi->maxs[i] = 0;
        }

        if (map[adios_statistic_sum] != -1)
        {
            MALLOC(sums, ntimes * sizeof(double *), "summation per timestep");
            MALLOC(vi->avgs, ntimes * sizeof(double *), "average per timestep");

            for (i = 0; i < ntimes; i++)
                sums[i] = vi->avgs[i] = 0;

            CALLOC(cnts, ntimes, sizeof(uint64_t), "count of elements per timestep");
        }

        if (map[adios_statistic_sum_square] != -1)
        {
            MALLOC(sum_squares, ntimes * sizeof(double *), "summation per timestep");
            MALLOC(vi->std_devs, ntimes * sizeof(double *), "standard deviation per timestep");

            for (i = 0; i < ntimes; i++)
                vi->std_devs[i] = sum_squares[i] = 0;
        }
    }

    if (map[adios_statistic_hist] != -1 && (var_root->characteristics[0].stats[0][map[adios_statistic_hist]].data))
    {
        struct adios_index_characteristics_stat_struct * stats = var_root->characteristics[0].stats[0];
        struct adios_index_characteristics_hist_struct * hist = stats[map[adios_statistic_hist]].data;
        int num_breaks = hist->num_breaks;

        MALLOC(vi->hist, sizeof(struct ADIOS_HIST), "histogram");
        MALLOC(vi->hist->breaks, num_breaks * sizeof(double), "break points of histogram");    
        MALLOC(vi->hist->gfrequencies, (num_breaks + 1) * sizeof(uint32_t), "global frequencies of histogram");

        vi->hist->num_breaks = hist->num_breaks;
        vi->hist->min = hist->min;
        vi->hist->max = hist->max;

        memcpy(vi->hist->breaks, hist->breaks, num_breaks * sizeof(double));
        CALLOC(vi->hist->gfrequencies, (num_breaks + 1), bp_get_type_size(adios_unsigned_integer, ""), "global frequency");

        if (ntimes > 0)
        {
            MALLOC(vi->hist->frequenciess, (ntimes * sizeof(int32_t *)), "frequencies for timesteps");
            for(i = 0; i < ntimes; i++)
                CALLOC(vi->hist->frequenciess[i], (num_breaks + 1), bp_get_type_size(adios_unsigned_integer, ""), "frequency at timestep");
        }
    }

    size = bp_get_type_size (var_root->type, "");    
    sum_size = bp_get_type_size (adios_double, "");

    if (var_root->type == adios_complex || var_root->type == adios_double_complex)
    {
        int type;
        count = 3;
        timestep = 0;

        if (var_root->type == adios_complex)
            type = adios_double;
        else
            type = adios_long_double;

        // Only a double precision returned for all complex values
        size = bp_get_type_size (adios_double, "");    

           for (i=0; i<var_root->characteristics_count; i++)
        {
            if (ntimes > 0)
                timestep = var_root->characteristics[i].time_index - 1;

            if (!var_root->characteristics[i].stats)
                continue;

            struct adios_index_characteristics_stat_struct ** stats = var_root->characteristics[i].stats;

            if ((map[adios_statistic_finite] != -1) && (* ((uint8_t *) stats[0][map[adios_statistic_finite]].data) == 0))
                continue;

            if (map[adios_statistic_min] != -1 && stats[0][map[adios_statistic_min]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double(type, stats[c][map[adios_statistic_min]].data);

                if(!vi->gmin) {
                    MALLOC (vi->gmin, count * size, "global minimum")
                    for (c = 0; c < count; c ++)
                           ((double * ) vi->gmin)[c] = data[c]; 

                } else {
                    for (c = 0; c < count; c ++)
                        if (data[c] < ((double *) vi->gmin)[c])
                               ((double * ) vi->gmin)[c] = data[c]; 
                }

                if (ntimes > 0) {
                    if(!vi->mins[timestep]) {
                        MALLOC (vi->mins[timestep], count * size, "minimum per timestep")
                        for (c = 0; c < count; c ++)
                               ((double **) vi->mins)[timestep][c] = data[c]; 

                    } else {
                        for (c = 0; c < count; c ++)
                            if (data[c] < ((double **) vi->mins)[timestep][c])
                                   ((double **) vi->mins)[timestep][c] = data[c]; 
                    }
                }
            }

            if (map[adios_statistic_max] != -1 && stats[0][map[adios_statistic_max]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double(type, stats[c][map[adios_statistic_max]].data);

                if(!vi->gmax) {
                    MALLOC (vi->gmax, count * size, "global minimum")
                    for (c = 0; c < count; c ++)
                        ((double * ) vi->gmax)[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        if (data[c] > ((double *) vi->gmax)[c])
                            ((double * ) vi->gmax)[c] = data[c];
                }

                if (ntimes > 0) {
                    if(!vi->maxs[timestep]) {
                        MALLOC (vi->maxs[timestep], count * size, "minimum per timestep")
                        for (c = 0; c < count; c ++)
                            ((double **) vi->maxs)[timestep][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            if (data[c] > ((double **) vi->maxs)[timestep][c])
                                ((double **) vi->maxs)[timestep][c] = data[c];
                    }
                }
            }

            if (map[adios_statistic_sum] != -1 && stats[0][map[adios_statistic_sum]].data)
            {    
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double(type, stats[c][map[adios_statistic_sum]].data);

                if(!gsum) {
                    MALLOC(gsum, count * sum_size, "global summation")
                    for (c = 0; c < count; c ++)
                           gsum[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        gsum[c] = gsum[c] + data[c];
                }

                if (ntimes > 0) {
                    if(!sums[timestep]) {
                        MALLOC(sums[timestep], count * sum_size, "summation per timestep")
                        for (c = 0; c < count; c ++)
                            sums[timestep][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            sums[timestep][c] = sums[timestep][c] + data[c];
                    }
                }
            }

            if (map[adios_statistic_sum_square] != -1 && stats[0][map[adios_statistic_sum_square]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double(type, stats[c][map[adios_statistic_sum_square]].data);

                if(!gsum_square) {
                    MALLOC(gsum_square, count * sum_size, "global summation of squares")
                    for (c = 0; c < count; c ++)
                        gsum_square[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                           gsum_square[c] = gsum_square[c] + data[c]; 
                }

                if (ntimes > 0) {
                    if(!sum_squares[timestep]) {
                        MALLOC(sum_squares[timestep], count * sum_size, "summation of square per timestep")
                        for (c = 0; c < count; c ++)
                            sum_squares[timestep][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            sum_squares[timestep][c] = sum_squares[timestep][c] + data[c]; 
                    }
                }
            }

            if (map[adios_statistic_cnt] != -1 && stats[0][map[adios_statistic_cnt]].data)
            {
                if (ntimes > 0)
                    cnts[timestep] += * ((uint32_t *) stats[0][map[adios_statistic_cnt]].data);
                gcnt += * (uint32_t *) stats[0][map[adios_statistic_cnt]].data;
            }
        }

        if(ntimes > 0 && vi->gmin && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            // min, max, summation exists only for arrays
            // Calculate average / timestep

            for(timestep = 0; timestep < ntimes; timestep ++) {
                MALLOC(vi->avgs[timestep], count * sum_size, "average per timestep")
                for (c = 0; c < count; c ++)
                    vi->avgs[timestep][c] = sums[timestep][c] / cnts[timestep];

                MALLOC(vi->std_devs[timestep], count * sum_size, "standard deviation per timestep")
                for (c = 0; c < count; c ++)
                    vi->std_devs[timestep][c] = sqrt((sum_squares[timestep][c] / cnts[timestep]) - (vi->avgs[timestep][c] * vi->avgs[timestep][c]));

                free (sums[timestep]);
                free (sum_squares[timestep]);
            }
        }

        // Calculate global average
        if(vi->gmin && gsum && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            MALLOC(vi->gavg, count * sum_size, "global average")

            if(gcnt > 0)
                for (c = 0; c < count; c ++)
                    vi->gavg[c] = gsum[c] / gcnt;
            else
                for (c = 0; c < count; c ++)
                    vi->gavg[c] = gsum[c];

            MALLOC(vi->gstd_dev, count * sum_size, "global average")
            if(vi->gavg && gcnt > 0)
                for (c = 0; c < count; c ++)
                    vi->gstd_dev[c] = sqrt(gsum_square[c] / gcnt - (vi->gavg[c] * vi->gavg[c]));
            else
                for (c = 0; c < count; c ++)
                    vi->gstd_dev[c] = 0;
        }
    }
    else
    {
        timestep = 0;
           for (i=0; i<var_root->characteristics_count; i++)
        {
            if (ntimes > 0)
                timestep = var_root->characteristics[i].time_index - 1;
                //timestep = i / (npgs / ntimes);

            if (!var_root->characteristics[i].stats)
                continue;

            struct adios_index_characteristics_stat_struct * stats = var_root->characteristics[i].stats[0];
            struct adios_index_characteristics_hist_struct * hist = stats[map[adios_statistic_hist]].data;

            if (map[adios_statistic_finite] != -1 && (* ((uint8_t *) stats[map[adios_statistic_finite]].data) == 0))
                continue;

            if (map[adios_statistic_min] != -1 && stats[map[adios_statistic_min]].data)
            {
                if(!vi->gmin) {
                    MALLOC (vi->gmin, size, "global minimum")
                       memcpy(vi->gmin, stats[map[adios_statistic_min]].data, size);

                } else if (adios_lt(var_root->type, stats[map[adios_statistic_min]].data, vi->gmin)){
                       memcpy(vi->gmin, stats[map[adios_statistic_min]].data, size);
                }

                if (ntimes > 0) {
                    if(!vi->mins[timestep]) {
                        MALLOC (vi->mins[timestep], size, "minimum per timestep")
                        memcpy(vi->mins[timestep], stats[map[adios_statistic_min]].data, size);

                    } else if (adios_lt(var_root->type, stats[map[adios_statistic_min]].data, vi->mins[timestep])) {
                        memcpy(vi->mins[timestep], stats[map[adios_statistic_min]].data, size);
                    }
                }
            }

            if (map[adios_statistic_max] != -1 && stats[map[adios_statistic_max]].data)
            {
                if(!vi->gmax) {
                    MALLOC (vi->gmax, size, "global maximum")
                    memcpy(vi->gmax, stats[map[adios_statistic_max]].data, size);

                } else if (adios_lt(var_root->type, vi->gmax, stats[map[adios_statistic_max]].data))
                       memcpy(vi->gmax, stats[map[adios_statistic_max]].data, size);
                
                if (ntimes > 0) {
                    if(!vi->maxs[timestep]) {
                        MALLOC (vi->maxs[timestep], size, "maximum per timestep")
                        memcpy(vi->maxs[timestep], stats[map[adios_statistic_max]].data, size);

                    } else if (adios_lt(var_root->type, vi->maxs[timestep], stats[map[adios_statistic_max]].data)) {
                        memcpy(vi->maxs[timestep], stats[map[adios_statistic_max]].data, size);    
                    }
                }
            }
    
            if (map[adios_statistic_sum] != -1 && stats[map[adios_statistic_sum]].data)
            {    
                if(!gsum) {
                    MALLOC(gsum, sum_size, "global summation")
                       memcpy(gsum, stats[map[adios_statistic_sum]].data, sum_size);

                } else {
                    *gsum = *gsum + * ((double *) stats[map[adios_statistic_sum]].data);
                }

                if (ntimes > 0) {
                    if(!sums[timestep]) {
                        MALLOC(sums[timestep], sum_size, "summation per timestep")
                        memcpy(sums[timestep], stats[map[adios_statistic_sum]].data, sum_size);

                    } else {
                        *sums[timestep] = *sums[timestep] + * ((double *) stats[map[adios_statistic_sum]].data);
                    }
                }
            }

            if (map[adios_statistic_sum_square] != -1 && stats[map[adios_statistic_sum_square]].data)
            {
                if(!gsum_square) {
                    MALLOC(gsum_square, sum_size, "global summation of squares")
                    memcpy(gsum_square, stats[map[adios_statistic_sum_square]].data, sum_size);

                } else {
                       *gsum_square = *gsum_square + * ((double *) stats[map[adios_statistic_sum_square]].data);
                }

                if (ntimes > 0) {
                    if(!sum_squares[timestep]) {
                        MALLOC(sum_squares[timestep], sum_size, "summation of square per timestep")
                        memcpy(sum_squares[timestep], stats[map[adios_statistic_sum_square]].data, sum_size);

                    } else {
                        *sum_squares[timestep] = *sum_squares[timestep] + * ((double *) stats[map[adios_statistic_sum_square]].data);
                    }
                }
            }

            if(map[adios_statistic_hist] != -1 && stats[map[adios_statistic_hist]].data)
            {
                for(j = 0; j <= vi->hist->num_breaks; j++)
                {    
                    uint32_t freq = hist->frequencies[j];
                    vi->hist->gfrequencies[j] += freq;
                    if (ntimes > 0)
                        vi->hist->frequenciess[timestep][j] += freq;
                }
            }

            if (map[adios_statistic_cnt] != -1 && stats[map[adios_statistic_cnt]].data)
            {
                if (ntimes > 0)
                    cnts[timestep] += * (uint32_t *) stats[map[adios_statistic_cnt]].data;
                gcnt += * (uint32_t *) stats[map[adios_statistic_cnt]].data;
            }
        }

        if(ntimes > 0 && vi->gmin && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            // min, max, summation exists only for arrays
            // Calculate average / timestep

            for(timestep = 0; timestep < ntimes; timestep ++) {
                MALLOC(vi->avgs[timestep], sum_size, "average per timestep")
                *(vi->avgs[timestep]) = *(sums[timestep]) / cnts[timestep];

                MALLOC(vi->std_devs[timestep], sum_size, "standard deviation per timestep")
                *(vi->std_devs[timestep]) = sqrt(*(sum_squares[timestep]) / cnts[timestep] - ((*(vi->avgs[timestep]) * (*(vi->avgs[timestep])))));

                free (sums[timestep]);
                free (sum_squares[timestep]);
            }
        }

        // Calculate global average
        if(vi->gmin && gsum && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            MALLOC(vi->gavg, sum_size, "global average")
            if(gcnt > 0)
                *vi->gavg = *gsum / gcnt;
            else
                vi->gavg = gsum;

            MALLOC(vi->gstd_dev, sum_size, "global average")
            if(vi->gavg && gcnt > 0)
                *vi->gstd_dev = sqrt(*gsum_square / gcnt - ((*(vi->gavg)) * (*(vi->gavg))));
            else
                *vi->gstd_dev = 0;
        }
    }

    if (!vi->value && vi->gmin) {
        vi->value = vi->gmin; // arrays have no value but we assign here the minimum
    }

    if(!vi->gmin) {
        vi->gmin = vi->value; // scalars have value but not min
    }
    if(!vi->gmax) {
        vi->gmax = vi->value; // scalars have value but not max
    }

    if (sums && gsum) {
        free (sums);
        free (gsum);
    }

    if (sum_squares && gsum_square) {
        free (sum_squares);
        free (gsum_square);
    }    
}

/* get local and global dimensions and offsets from a variable characteristics 
   return: 1 = it is a global array, 0 = local array
*/
static int adios_read_bp_get_dimensioncharacteristics(struct adios_index_characteristic_struct_v1 *ch, 
                                               uint64_t *ldims, uint64_t *gdims, uint64_t *offsets)
{
    int is_global=0; // global array or just an array written by one process?
    int ndim = ch->dims.count;
    int k;

    for (k=0; k < ndim; k++) {
        ldims[k]   = ch->dims.dims[k*3];
        gdims[k]   = ch->dims.dims[k*3+1];
        offsets[k] = ch->dims.dims[k*3+2];
        is_global = is_global || gdims[k];
    }
    return is_global;
}

static void adios_read_bp_get_dimensions (struct adios_index_var_struct_v1 *var_root, int ntsteps, int file_is_fortran,
                                  int *ndim, uint64_t **dims, int *timedim)
{
    int i, k;
    int is_global; // global array or just an array written by one process?
    uint64_t ldims[32];
    uint64_t gdims[32];
    uint64_t offsets[32];

    /* Get dimension information */    
    *ndim = var_root->characteristics [0].dims.count; 
    *dims = NULL;
    *timedim = -1;
    if (*ndim == 0) {
        /* done with this scalar variable */
        return ;
    }

    *dims = (uint64_t *) malloc (sizeof(uint64_t) * (*ndim));
    memset(*dims,0,sizeof(uint64_t)*(*ndim));

    is_global = adios_read_bp_get_dimensioncharacteristics( &(var_root->characteristics[0]),
                                                    ldims, gdims, offsets);

    if (!is_global) {
        /* local array */
        for (i=0; i < *ndim; i++) {
            /* size of time dimension is always registered as 1 for an array */
            if (ldims[i] == 1 && var_root->characteristics_count > 1) {
                *timedim = i;
                (*dims)[i] = ntsteps;
            } else {
                (*dims)[i] = ldims[i];
            }
            gdims[i] = ldims[i];
        }         
    }         
    else {
        /* global array:
           time dimension: ldims=1, gdims=0
           in C array, it can only be the first dim
           in Fortran array, it can only be the last dim
           (always the slowest changing dim)
        */
        if (gdims[*ndim-1] == 0)
        {
            if (!file_is_fortran) {
                /* first dimension is the time (C array)
                 * ldims[0] = 1 but gdims does not contain time info and 
                 * gdims[0] is 1st data dimension and 
                 * gdims is shorter by one value than ldims in case of C.
                 * Therefore, gdims[*ndim-1] = 0 if there is a time dimension. 
                 */
                *timedim = 0;
                // error check
                if (*ndim > 1 && ldims[0] != 1) {
                    fprintf(stderr,"ADIOS Error: this is a BP file with C ordering but we didn't find"
                            "an array to have time dimension in the first dimension. l:g:o = (");
                    for (i=0; i < *ndim; i++) {
                        fprintf(stderr,"%llu:%llu:%llu%s", ldims[i], gdims[i], offsets[i], (i<*ndim-1 ? ", " : "") );
                    }
                    fprintf(stderr, ")\n");
                }
                (*dims)[0] = ntsteps;
                for (i=1; i < *ndim; i++) 
                    (*dims)[i]=gdims[i-1];
            } else {
                // last dimension is the time (Fortran array)
                *timedim = *ndim - 1;

                if (*ndim > 1 && ldims[*ndim-1] != 1) {
                    fprintf(stderr,"ADIOS Error: this is a BP file with Fortran array ordering but we didn't find"
                            "an array to have time dimension in the last dimension. l:g:o = (");
                    for (i=0; i < *ndim; i++) {
                        fprintf(stderr,"%llu:%llu:%llu%s", ldims[i], gdims[i], offsets[i], (i<*ndim-1 ? ", " : "") );
                    }
                    fprintf(stderr, ")\n");
                }
                for (i=0; i < *ndim-1; i++) 
                    (*dims)[i]=gdims[i];
                (*dims)[*timedim] = ntsteps;
            }
        } else {
            // no time dimenstion
            for (i=0; i < *ndim; i++) 
                 (*dims)[i]=gdims[i];
        }
    }
}

ADIOS_VARINFO * adios_read_bp_inq_var (ADIOS_GROUP *gp, const char * varname) 
{
    int varid = adios_read_bp_find_var(gp, varname);
    if (varid < 0)
        return NULL;
    return adios_read_bp_inq_var_byid(gp, varid);
}

ADIOS_VARINFO * adios_read_bp_inq_var_byid (ADIOS_GROUP *gp, int varid)
{
    struct BP_GROUP      * gh;
    struct BP_FILE       * fh;
    ADIOS_VARINFO * vi;
    int file_is_fortran;
    struct adios_index_var_struct_v1 * var_root;
    int i,k;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_inq_var()");
        return NULL;
    }
    gh = (struct BP_GROUP *) gp->gh;
    if (!gh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return NULL;
    }
    fh = gh->fh;
    if (!fh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return NULL;
    }
    if (varid < 0 || varid >= gh->vars_count) {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return NULL;
    }
    vi = (ADIOS_VARINFO *) malloc(sizeof(ADIOS_VARINFO));
    if (!vi) {
        adios_error ( err_no_memory, "Could not allocate memory for variable info");
        return NULL;
    }


    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
    
    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i=0; i<varid && var_root; i++) {
        var_root = var_root->next;
    }

    if (i!=varid) {
        adios_error (err_corrupted_variable, "Variable id=%d is valid but was not found in internal data structures!",varid);
        return NULL; 
    }


    vi->varid = varid;

    vi->type = var_root->type;
    if (!var_root->characteristics_count) {
        adios_error (err_corrupted_variable, "Variable %s does not have information on dimensions", 
              gp->var_namelist[varid]);
        free(vi);
        return NULL;
    }
    vi->characteristics_count = var_root->characteristics_count;

    /* Get value or min/max */

    adios_read_bp_get_dimensions (var_root, fh->tidx_stop - fh->tidx_start + 1, file_is_fortran, 
                          &(vi->ndim), &(vi->dims), &(vi->timedim));
    
    if (file_is_fortran != futils_is_called_from_fortran()) {
        /* If this is a Fortran written file and this is called from C code,  
           or this is a C written file and this is called from Fortran code ==>
           We need to reverse the order of the dimensions */
        swap_order(vi->ndim, vi->dims, &(vi->timedim));
        /*printf("File was written from %s and read now from %s, so we swap order of dimensions\n",
                (file_is_fortran ? "Fortran" : "C"), (futils_is_called_from_fortran() ? "Fortran" : "C"));*/
    }
    
    adios_read_bp_get_characteristics (var_root, vi);

    return vi;
}

void adios_read_bp_free_varinfo (ADIOS_VARINFO *vp)
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

#define MPI_FILE_READ_OPS                           \
        bp_realloc_aligned(fh->b, slice_size);      \
        fh->b->offset = 0;                          \
                                                    \
        MPI_File_seek (fh->mpi_fh                   \
                      ,(MPI_Offset)slice_offset     \
                      ,MPI_SEEK_SET                 \
                      );                            \
                                                    \
        MPI_File_read (fh->mpi_fh                   \
                      ,fh->b->buff                  \
                      ,slice_size                   \
                      ,MPI_BYTE                     \
                      ,&status                      \
                      );                            \
        fh->b->offset = 0;                          \

//We also need to be able to read old .bp which doesn't have 'payload_offset'
#define MPI_FILE_READ_OPS1                                                                  \
        MPI_File_seek (fh->mpi_fh                                                           \
                      ,(MPI_Offset) var_root->characteristics[start_idx + idx].offset       \
                      ,MPI_SEEK_SET);                                                       \
        MPI_File_read (fh->mpi_fh, fh->b->buff, 8, MPI_BYTE, &status);                      \
        tmpcount= *((uint64_t*)fh->b->buff);                                                \
                                                                                            \
        bp_realloc_aligned(fh->b, tmpcount + 8);                                            \
        fh->b->offset = 0;                                                                  \
                                                                                            \
        MPI_File_seek (fh->mpi_fh                                                           \
                      ,(MPI_Offset) (var_root->characteristics[start_idx + idx].offset)     \
                      ,MPI_SEEK_SET);                                                       \
        MPI_File_read (fh->mpi_fh, fh->b->buff, tmpcount + 8, MPI_BYTE, &status);           \
        fh->b->offset = 0;                                                                  \
        adios_parse_var_data_header_v1 (fh->b, &var_header);                                \


// To read subfiles
#define MPI_FILE_READ_OPS2                                                                  \
        bp_realloc_aligned(fh->b, slice_size);                                              \
        fh->b->offset = 0;                                                                  \
                                                                                            \
        MPI_File * sfh;                                                                     \
        sfh = get_BP_file_handle (fh->sfh                                                   \
                                 ,var_root->characteristics[start_idx + idx].file_index     \
                                 );                                                         \
        if (!sfh)                                                                           \
        {                                                                                   \
            int err;                                                                        \
            char * ch, * name_no_path, * name;                                              \
            struct BP_file_handle * new_h =                                                 \
                  (struct BP_file_handle *) malloc (sizeof (struct BP_file_handle));        \
            new_h->file_index = var_root->characteristics[start_idx + idx].file_index;      \
            new_h->next = 0;                                                                \
            if (ch = strrchr (fh->fname, '/'))                                              \
            {                                                                               \
                name_no_path = malloc (strlen (ch + 1) + 1);                                \
                strcpy (name_no_path, ch + 1);                                              \
            }                                                                               \
            else                                                                            \
            {                                                                               \
                name_no_path = malloc (strlen (fh->fname) + 1);                             \
                strcpy (name_no_path, fh->fname);                                           \
            }                                                                               \
                                                                                            \
            name = malloc (strlen (fh->fname) + 5 + strlen (name_no_path) + 1 + 10 + 1);    \
            sprintf (name, "%s.dir/%s.%d", fh->fname, name_no_path, new_h->file_index);     \
                                                                                            \
            err = MPI_File_open (fh->comm                                                   \
                                ,name                                                       \
                                ,MPI_MODE_RDONLY                                            \
                                ,(MPI_Info)MPI_INFO_NULL                                    \
                                ,&new_h->fh                                                 \
                                );                                                          \
           if (err != MPI_SUCCESS)                                                          \
           {                                                                                \
               fprintf (stderr, "can not open file %S\n", name);                            \
               return -1;                                                                   \
           }                                                                                \
                                                                                            \
           add_BP_file_handle (&fh->sfh                                                     \
                              ,new_h                                                        \
                              );                                                            \
           sfh = &new_h->fh;                                                                \
                                                                                            \
           free (name_no_path);                                                             \
           free (name);                                                                     \
        }                                                                                   \
                                                                                            \
        MPI_File_seek (*sfh                                                                 \
                      ,(MPI_Offset)slice_offset                                             \
                      ,MPI_SEEK_SET                                                         \
                      );                                                                    \
        MPI_File_read (*sfh                                                                 \
                      ,fh->b->buff                                                          \
                      ,slice_size                                                           \
                      ,MPI_BYTE                                                             \
                      ,&status                                                              \
                      );                                                                    \
        fh->b->offset = 0;                                                                  \


int64_t adios_read_bp_read_var (ADIOS_GROUP * gp, const char * varname,
                        const uint64_t * start, const uint64_t * count,
                        void * data)
{
    struct BP_GROUP      * gh;
    struct BP_FILE       * fh;
    int varid, has_subfile;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        return -adios_errno;
    }

    gh = (struct BP_GROUP *) gp->gh;
    if (!gh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return -adios_errno;
    }

    fh = gh->fh;
    if (!fh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return -adios_errno;
    }

    varid = adios_read_bp_find_var(gp, varname);
    if (varid < 0 || varid >= gh->vars_count) {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return -adios_errno;
    }

    return adios_read_bp_read_var_byid(gp, varid, start, count, data);
}

/***********************************************
 * This routine is to read in data in a 'local *
 * array fashion (as opposed to global array)  *
 *     Q. Liu, 11/2010                         *
 ***********************************************/
int64_t adios_read_bp_read_local_var (ADIOS_GROUP * gp, const char * varname,
                                      int vidx, const uint64_t * start,
                                      const uint64_t * count, void * data)
{
    struct BP_GROUP      * gh;
    struct BP_FILE       * fh;
    struct adios_index_var_struct_v1 * var_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    int    i,j,k, t, varid, start_idx, idx;
    int    ndim, ndim_notime, has_subfile, file_is_fortran;
    uint64_t size, * dims;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t datasize, nloop, dset_stride,var_stride, total_size=0, items_read;
    uint64_t count_notime[32], start_notime[32];
    int timedim = -1, temp_timedim, is_global = 0, size_of_type;
    uint64_t slice_offset, slice_size, tmpcount = 0;
    uint64_t datatimeoffset = 0; // offset in data to write a given timestep
    MPI_Status status;

    adios_errno = 0;
    if (!gp)
    {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        return -adios_errno;
    }

    gh = (struct BP_GROUP *) gp->gh;
    if (!gh)
    {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return -adios_errno;
    }

    fh = gh->fh;
    if (!fh)
    {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return -adios_errno;
    }

    varid = adios_read_bp_find_var(gp, varname);
    if (varid < 0 || varid >= gh->vars_count)
    {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return -adios_errno;
    }

    /* Check if file is written out by Fortran or C */
    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);

    /* Check whether we need to handle subfiles */
    has_subfile = fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE;

    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i = 0; i< varid && var_root; i++)
    {
        var_root = var_root->next;
    }

    if (i != varid)
    {
        adios_error (err_corrupted_variable, 
                     "Variable id=%d is valid but was not found in internal data structures!",
                     varid);
        return -adios_errno; 
    }

    if (vidx < 0 || vidx >= var_root->characteristics_count)
    {
        adios_error (err_out_of_bound, "idx=%d is out of bound", vidx);
    }

    ndim = var_root->characteristics [vidx].dims.count;

    /* count_notime/start_notime are working copies of count/start */
    for (i = 0; i < ndim; i++)
    {
        count_notime[i] = count[i];
        start_notime[i] = start[i];
    }

    ndim_notime = ndim;

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if (futils_is_called_from_fortran())
    {
        timedim = -1;
        swap_order(ndim_notime, count_notime, &timedim);
        swap_order(ndim_notime, start_notime, &timedim);
    }
    
    /* items_read = how many data elements are we going to read in total */
    items_read = 1;
    for (i = 0; i < ndim_notime; i++)
        items_read *= count_notime[i];

    size_of_type = bp_get_type_size (var_root->type, var_root->characteristics [vidx].value);

    /* READ A SCALAR VARIABLE */
    if (ndim_notime == 0)
    {
        slice_size = size_of_type;
        start_idx = 0; // OPS macros below need it
        idx = vidx; // OPS macros below need it

        if (var_root->type == adios_string)
        {
            // strings are stored without \0 in file
            // size_of_type here includes \0 so decrease by one
            size_of_type--;
        }

        /* Old BP files don't have payload_offset characteristic */
        if (var_root->characteristics[vidx].payload_offset > 0)
        {
            slice_offset = var_root->characteristics[vidx].payload_offset;

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }
        }
        else
        {
            slice_offset = 0;
            MPI_FILE_READ_OPS1
        }

        memcpy((char *)data + total_size, fh->b->buff + fh->b->offset, size_of_type);

        if (fh->mfooter.change_endianness == adios_flag_yes)
            change_endianness((char *)data + total_size
                             ,size_of_type
                             ,var_root->type
                             );

        if (var_root->type == adios_string)
        {
            // add \0 to the end of string
            // size_of_type here is the length of string
            // FIXME: how would this work for strings written over time?
            ((char*)data + total_size)[size_of_type] = '\0';
        }

        total_size += size_of_type;

        return total_size;
    } /* READ A SCALAR VARIABLE END HERE */

    /* READ AN ARRAY VARIABLE */
    uint64_t write_offset = 0;
    int npg = 0;
    tmpcount = 0;
    int flag;
    datasize = 1;
    nloop = 1;
    var_stride = 1;
    dset_stride = 1;
    uint64_t payload_size = size_of_type;

    /* To get ldims for the index vidx */
    adios_read_bp_get_dimensioncharacteristics( &(var_root->characteristics[vidx]),
                                                ldims, gdims, offsets);

    /* Again, a Fortran written file has the dimensions in Fortran order we need to swap here */
    /* Only local dims are needed for reading local vars */ 
    if (file_is_fortran)
    {
        i=-1;
        swap_order(ndim, ldims, &(i));
    }

    /*
    printf("ldims   = "); for (j = 0; j < ndim; j++) printf("%d ",ldims[j]); printf("\n");
    printf("count_notime   = "); for (j = 0; j < ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
    printf("start_notime   = "); for (j = 0; j < ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
    */        

    for (j = 0; j < ndim_notime; j++)
    {
        payload_size *= ldims [j];
    
        if ( (start_notime[j] > ldims[j]) 
            || (start_notime[j] + count_notime[j] > ldims[j]))
        {
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j+1, count_notime[j], start_notime[j], ldims[j] - 1);
                    return -adios_errno;
        }
    }

    /* determined how many (fastest changing) dimensions can we read in in one read */
    int break_dim =  ndim_notime - 1;
    while (break_dim > -1)
    {
        if (start_notime[break_dim] == 0 && ldims[break_dim] == count_notime[break_dim])
        {
            datasize *= ldims[break_dim];
        }
        else
            break;
        
        break_dim--;
    }
    
    slice_offset = 0;
    slice_size = 0;
    /* Note: MPI_FILE_READ_OPS  - for reading single BP file.
     *       MPI_FILE_READ_OPS2 - for reading those with subfiles.
     *       MPI_FILE_READ_OPS1 - for reading old version of BP files
     *                            which don't contain "payload_offset"
     * Whenever to use OPS macro, start_idx and idx variable needs to be
     * properly set.
     */
    
    start_idx = 0;
    idx = vidx;

    if (break_dim <= 0) 
    {
        /* The slowest changing dimensions should not be read completely but
           we still need to read only one block */
   
        uint64_t size_in_dset = count_notime[0];
        uint64_t offset_in_dset = start_notime[0];

        slice_size = (break_dim == -1 ? datasize * size_of_type : size_in_dset * datasize * size_of_type);
    
        if (var_root->characteristics[start_idx + idx].payload_offset > 0)
        {
            slice_offset = var_root->characteristics[start_idx + idx].payload_offset 
                         + offset_in_dset * datasize * size_of_type;

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }
        }
        else
        {
            slice_offset = 0;
            MPI_FILE_READ_OPS1
        }

        memcpy ((char *)data, fh->b->buff + fh->b->offset, slice_size);
        if (fh->mfooter.change_endianness == adios_flag_yes)
        {
            change_endianness((char *)data + write_offset, slice_size, var_root->type);
        }
    }
    else 
    {
        uint64_t stride_offset = 0;
        uint64_t * size_in_dset, * offset_in_dset, * offset_in_var;
        uint64_t start_in_payload, end_in_payload, s;
        uint64_t var_offset;
        uint64_t dset_offset;

        size_in_dset = (uint64_t *) malloc (8 * ndim_notime);
        offset_in_dset = (uint64_t *) malloc (8 * ndim_notime);
        offset_in_var = (uint64_t *) malloc (8 * ndim_notime);
 
        if (size_in_dset == 0 || offset_in_dset == 0 || offset_in_var == 0)
        {
             adios_error (err_no_memory, "Malloc failed in %s at %d\n"
                         , __FILE__, __LINE__
                         );
             return -adios_errno;
        }

        for (i = 0; i < ndim_notime ; i++)
        {
            size_in_dset[i] = count_notime[i];
            offset_in_dset[i] = start_notime[i];
            offset_in_var[i] = 0;
        }
 
        datasize = 1;
        var_stride = 1;
        for (i = ndim_notime - 1; i >= break_dim; i--)
        {
            datasize *= size_in_dset[i];
            dset_stride *= ldims[i];
            var_stride *= count_notime[i];
        }

        /* Calculate the size of the chunk we are trying to read in */
        start_in_payload = 0;
        end_in_payload = 0;
        s = 1;
        for (i = ndim_notime - 1; i >= 0; i--)
        {
            start_in_payload += s * offset_in_dset[i] * size_of_type;
            end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
            s *= ldims[i];
        }
        slice_size = end_in_payload - start_in_payload + 1 * size_of_type;
 
        if (var_root->characteristics[start_idx + idx].payload_offset > 0)
        {
            slice_offset =  var_root->characteristics[start_idx + idx].payload_offset
                          + start_in_payload;
            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }
 
            for ( i = 0; i < ndim_notime ; i++)
            {
                offset_in_dset[i] = 0;
            }
        }
        else
        {
            slice_offset =  start_in_payload;
            MPI_FILE_READ_OPS1
        }

        var_offset = 0;
        dset_offset = 0;
        for (i = 0; i < ndim_notime ; i++)
        {
            var_offset = offset_in_var[i] + var_offset * count_notime[i];
            dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
        }

        copy_data (data
                  ,fh->b->buff + fh->b->offset
                  ,0
                  ,break_dim
                  ,size_in_dset
                  ,ldims
                  ,count_notime
                  ,var_stride
                  ,dset_stride
                  ,var_offset
                  ,dset_offset
                  ,datasize
                  ,size_of_type 
                  );

        free (size_in_dset);
        free (offset_in_dset);
        free (offset_in_var);
    }
    
    total_size += items_read * size_of_type;

    return total_size;
}

// The purpose of keeping this function is to be able
// to read in old BP files. Can be deleted later on.
int64_t adios_read_bp_read_var_byid1 (ADIOS_GROUP    * gp,
                             int              varid,
                             const uint64_t  * start,
                             const uint64_t  * count,
                             void           * data)
{
    struct BP_GROUP      * gh;
    struct BP_FILE       * fh;
    int file_is_fortran;
    struct adios_index_var_struct_v1 * var_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    int    i,j,k, idx, timestep;
    int    start_time, stop_time;
    int    pgoffset, pgcount, next_pgoffset,start_idx, stop_idx;
    int    ndim, ndim_notime;  
    uint64_t size;
    uint64_t *dims;
    uint64_t ldims[32];
    uint64_t gdims[32];
    uint64_t offsets[32];
    uint64_t datasize, nloop, dset_stride,var_stride, total_size=0, items_read;
    uint64_t count_notime[32], start_notime[32];
    MPI_Status status;
    int timedim = -1, temp_timedim, timedim_c;
    int rank;
    int is_global = 0;
    int size_of_type;
    uint64_t slice_offset;
    uint64_t slice_size;
    uint64_t tmpcount = 0;
    uint64_t datatimeoffset = 0; // offset in data to write a given timestep

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        return -adios_errno;
    }
    gh = (struct BP_GROUP *) gp->gh;
    if (!gh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return -adios_errno;
    }
    fh = gh->fh;
    if (!fh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return -adios_errno;
    }
    if (varid < 0 || varid >= gh->vars_count) {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return -adios_errno;
    }
    
    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
    
    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i=0; i<varid && var_root; i++) {
        var_root = var_root->next;
    }

    if (i!=varid) {
        adios_error (err_corrupted_variable, "Variable id=%d is valid but was not found in internal data structures!",varid);
        return -adios_errno; 
    }

    /* Get dimensions and flip if caller != writer language */
    adios_read_bp_get_dimensions (var_root, fh->tidx_stop - fh->tidx_start + 1, file_is_fortran, 
                          &ndim, &dims, &timedim);

    /* Here the cases in which .bp written from Fortran and C are considered separately.
       1) bp written from Fortran */
    if (file_is_fortran)
    {
        /* Get the timesteps we need to read */
        if (timedim > -1) 
        {
            if (timedim != ndim - 1)
            {
                adios_error (err_no_data_at_timestep,"Variable (id=%d) has wrong time dimension index",
                      varid);
                return -adios_errno;
            }
            if (futils_is_called_from_fortran())
            {
                start_time = start[timedim] + fh->tidx_start;
                stop_time = start_time + count[timedim] - 1;
            }
            else
            {
                start_time = start[0] + fh->tidx_start;
                stop_time = start_time + count[0] - 1;
            }
        }
        else 
        {
            /* timeless variable, but we still need to handle the case that
               var is not written in the first few timesteps. 
               This happens in Pixie3D.  */
            for (i = 0; i < fh->mfooter.time_steps; i++)
            {
                pgoffset = fh->gvar_h->time_index[0][gh->group_id][i];
                if (i < fh->mfooter.time_steps - 1)
                    next_pgoffset = fh->gvar_h->time_index[0][gh->group_id][i + 1];
                else
                    next_pgoffset = -1;

                if (fh->gvar_h->pg_offsets[pgoffset] < var_root->characteristics[0].offset
                && (i == fh->mfooter.time_steps - 1 
                   ||fh->gvar_h->pg_offsets[next_pgoffset] > var_root->characteristics[0].offset)
                )
                {
                    start_time = fh->tidx_start + i;
                    stop_time = start_time;
                    break;
                }
            }
        }

        /* flip dims and timedim to C order */
        swap_order(ndim, dims, &timedim);

        /* Take out the time dimension from start[] and count[] */
        /* if we have time dimension */
        if (timedim > -1)
        {
            j = 0;
            if (futils_is_called_from_fortran())
                temp_timedim = ndim - 1;
            else
                temp_timedim = 0;

            for (i = 0; i < temp_timedim; i++)
            {
                count_notime[j] = count[i];
                start_notime[j] = start[i];
                j++;
            }
            i++; // skip timedim
            for (; i < ndim; i++)
            {
                count_notime[j] = count[i];
                start_notime[j] = start[i];
                j++;
            }
            ndim_notime = ndim-1;
        }
        else
        /* if we don't have time dimension */
        {
            for (i = 0; i < ndim; i++)
            {
                count_notime[i] = count[i];
                start_notime[i] = start[i];
            }
            ndim_notime = ndim;
        }
    }
    /* 2) .bp written by C */
    else
    {
        /* Get the timesteps we need to read */
        if (timedim > -1) 
        {
            /* timedim has to be the 1st dimension. To be extended to handle 
               the cases timedim at any dimension */
            if (timedim != 0)
            {
                adios_error (err_no_data_at_timestep,"Variable (id=%d) has wrong time dimension",
                      varid);
                return -adios_errno;
            }

            if (futils_is_called_from_fortran())
            {
                start_time = start[ndim - 1] + fh->tidx_start;
                stop_time = start_time + count[ndim -1] - 1;
            }
            else
            {
                start_time = start[0] + fh->tidx_start;
                stop_time = start_time + count[0] - 1;
            }

            start_time = start[timedim] + fh->tidx_start;
            stop_time = start_time + count[timedim] - 1;
        }
        else 
        {
            /* timeless variable */
            start_time = fh->tidx_start;
            stop_time = fh->tidx_start;
        }

        /* No need to flip dims, timedim as they are already in C order. */
        //swap_order(ndim, dims, &timedim);

        /* Take out the time dimension from start[] and count[] */
        if (timedim == -1) /* timeless variable */ 
        {
            for (i = 0; i < ndim; i++) 
            {
                count_notime[i] = count[i];
                start_notime[i] = start[i];
            }
            ndim_notime = ndim;
        }
        /* if we have time dimension */
        else
        {
            j = 0;
            if (futils_is_called_from_fortran())
                temp_timedim = ndim - 1;
            else
                temp_timedim = 0;

            for (i = 0; i < temp_timedim; i++)
            {
                count_notime[j] = count[i];
                start_notime[j] = start[i];
                j++;
            }
            i++; // skip timedim
            for (; i < ndim; i++)
            {
                count_notime[j] = count[i];
                start_notime[j] = start[i];
                j++;
            }
            ndim_notime = ndim - 1;
        }
    }

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if ( futils_is_called_from_fortran()) {
        swap_order(ndim_notime, count_notime, &timedim);
        swap_order(ndim_notime, start_notime, &timedim);
    }
    
    /* items_read = how many data elements are we going to read in total (per timestep) */
    items_read = 1;
    for (i = 0; i < ndim_notime; i++)
        items_read *= count_notime[i];
    
    MPI_Comm_rank(gh->fh->comm, &rank);

    size_of_type = bp_get_type_size (var_root->type, var_root->characteristics [0].value);

    /* For each timestep, do reading separately (they are stored in different sets of process groups */
    for (timestep = start_time; timestep <= stop_time; timestep++) {

        // pgoffset = the starting offset for the given time step
        // pgcount  = number of process groups of that time step
        pgoffset = fh->gvar_h->time_index[0][gh->group_id][timestep - fh->tidx_start];
        pgcount = fh->gvar_h->time_index[1][gh->group_id][timestep - fh->tidx_start];

        start_idx = -1;
        for (i=0;i<var_root->characteristics_count;i++) {
            if (   (  var_root->characteristics[i].offset > fh->gvar_h->pg_offsets[pgoffset])
                && (  (pgoffset + pgcount == fh->mfooter.pgs_count) 
                    ||(  var_root->characteristics[i].offset < fh->gvar_h->pg_offsets[pgoffset + 1]))
               ) 
            {
                start_idx = i;
                break;
            }
        }
/*
printf ("var_root->characteristics_count = %d\n", var_root->characteristics_count);
printf ("pg offset 0 = %lld\n", fh->gvar_h->pg_offsets[pgoffset]);
printf ("pg offset 1 = %lld\n", fh->gvar_h->pg_offsets[pgoffset + 1]);
printf ("var offset 3 = %lld\n", var_root->characteristics[3].offset);
printf ("var offset 2 = %lld\n", var_root->characteristics[2].offset);
printf ("var offset 1 = %lld\n", var_root->characteristics[1].offset);
printf ("pgcount = %lld\n", pgcount);
*/
        for (i=var_root->characteristics_count-1;i>-1;i--) {
            if (   (  var_root->characteristics[i].offset > fh->gvar_h->pg_offsets[pgoffset])
                && (  (pgoffset + pgcount == fh->mfooter.pgs_count)
                    ||(  var_root->characteristics[i].offset < fh->gvar_h->pg_offsets[pgoffset + pgcount]))
               )
            {
                stop_idx = i;
                break;
            }
        }

        if (start_idx<0) {
            adios_error (err_no_data_at_timestep,"Variable (id=%d) has no data at %d time step",
                varid, timestep);
            return -adios_errno;
        }

        if (ndim_notime == 0) {
            /* READ A SCALAR VARIABLE */

            slice_size = size_of_type;
            idx = 0; // macros below need it

            if (var_root->type == adios_string) {
                // strings are stored without \0 in file
                // size_of_type here includes \0 so decrease by one
                size_of_type--;
            }

            if (var_root->characteristics[start_idx+idx].payload_offset > 0) {
                slice_offset = var_root->characteristics[start_idx+idx].payload_offset;
                MPI_FILE_READ_OPS
            } else {
                slice_offset = 0;
                MPI_FILE_READ_OPS1
            }

            memcpy((char *)data+total_size, fh->b->buff + fh->b->offset, size_of_type);
            if (fh->mfooter.change_endianness == adios_flag_yes) {
                change_endianness((char *)data+total_size, size_of_type, var_root->type);
            }

            if (var_root->type == adios_string) {
                // add \0 to the end of string
                // size_of_type here is the length of string
                // FIXME: how would this work for strings written over time?
                ((char*)data+total_size)[size_of_type] = '\0';
            }

            total_size += size_of_type;
            continue;
        }

        /* READ AN ARRAY VARIABLE */
        //int * idx_table = (int *) malloc (sizeof(int) * pgcount);
        //int * idx_table = (int *) malloc (sizeof(int) * (var_root->characteristics_count - start_idx));
        int * idx_table = (int *) malloc (sizeof(int) * (stop_idx - start_idx + 1));

        uint64_t write_offset = 0;
        int npg = 0;
        tmpcount = 0;
        if (pgcount > var_root->characteristics_count)
            pgcount = var_root->characteristics_count;

        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++) {
            int flag;
            datasize = 1;
            nloop = 1;
            var_stride = 1;
            dset_stride = 1;
            idx_table[idx] = 1;
            uint64_t payload_size = size_of_type;
    
            /* Each pg can have a different sized array, so we need the actual dimensions from it */
            is_global = adios_read_bp_get_dimensioncharacteristics( &(var_root->characteristics[start_idx + idx]),
                                                            ldims, gdims, offsets);
            if (!is_global) {
                // we use gdims below, which is 0 for a local array; set to ldims here
                for (j = 0; j< ndim; j++) {
                    gdims[j]=ldims[j];
                }
            }

            /* Again, a Fortran written file has the dimensions in Fortran order we need to swap here */
            //if (file_is_fortran != futils_is_called_from_fortran()) {
            if (file_is_fortran ) {
                i=-1;
                swap_order(ndim, gdims,   &(i)); // i is dummy 
                swap_order(ndim, ldims,   &(i));
                swap_order(ndim, offsets, &(i));
            }
            
            /* take out the time dimension */
            /* For C, gdims and offset are one size shorter because the timedim part is missing,
               so we take it out only for fortran files
            */
            if (timedim > -1) {
                for (i = timedim; i < ndim-1; i++) {
                    ldims[i] = ldims[i+1];
                    if (file_is_fortran) {
                        gdims[i] = gdims[i+1];
                        offsets[i] = offsets[i+1];
                    }
                }
            }
            /*
            printf("ldims   = "); for (j = 0; j<ndim; j++) printf("%d ",ldims[j]); printf("\n");
            printf("gdims   = "); for (j = 0; j<ndim; j++) printf("%d ",gdims[j]); printf("\n");
            printf("offsets = "); for (j = 0; j<ndim; j++) printf("%d ",offsets[j]); printf("\n");
            printf("count_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
            printf("start_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
            */
            for (j = 0; j < ndim_notime; j++) {
    
                payload_size *= ldims [j];
    
                if ( (count_notime[j] > gdims[j]) 
                  || (start_notime[j] > gdims[j]) 
                  || (start_notime[j] + count_notime[j] > gdims[j])){
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j+1, count_notime[j], start_notime[j], gdims[j] - 1);
                    return -adios_errno;
                }
    
                /* check if there is any data in this pg and this dimension to read in */
                flag = (offsets[j] >= start_notime[j] 
                        && offsets[j] < start_notime[j] + count_notime[j])
                    || (offsets[j] < start_notime[j]
                        && offsets[j] + ldims[j] > start_notime[j] + count_notime[j]) 
                    || (offsets[j] + ldims[j] > start_notime[j] 
                        && offsets[j] + ldims[j] <= start_notime[j] + count_notime[j]);
                idx_table [idx] = idx_table[idx] && flag;
            }
            
            if ( !idx_table[idx] ) {
                continue;
            }
            ++npg;

            /* determined how many (fastest changing) dimensions can we read in in one read */
            int hole_break; 
            for (i = ndim_notime - 1; i > -1; i--) {
                if (offsets[i] == start_notime[i] && ldims[i] == count_notime[i]) {
                    datasize *= ldims[i];
                }
                else
                    break;
            }
    
            hole_break = i;
            slice_offset = 0;
            slice_size = 0;

            if (hole_break == -1) {
                /* The complete read happens to be exactly one pg, and the entire pg */
                /* This means we enter this only once, and npg=1 at the end */
                /* This is a rare case. FIXME: cannot eliminate this? */
                slice_size = payload_size;

                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset = var_root->characteristics[start_idx + idx].payload_offset;
                    MPI_FILE_READ_OPS
                } else {
                    slice_offset = 0;
                    MPI_FILE_READ_OPS1
                }
    
                memcpy( (char *)data, fh->b->buff + fh->b->offset, slice_size);
                if (fh->mfooter.change_endianness == adios_flag_yes) {
                    change_endianness(data, slice_size, var_root->type);
                }
            }
            else if (hole_break == 0) 
            {
                /* The slowest changing dimensions should not be read completely but
                   we still need to read only one block */
                int isize;
                uint64_t size_in_dset = 0;
                uint64_t offset_in_dset = 0;
                uint64_t offset_in_var = 0;
    
                isize = offsets[0] + ldims[0];
                if (start_notime[0] >= offsets[0]) {
                    // head is in
                    if (start_notime[0]<isize) {
                        if (start_notime[0] + count_notime[0] > isize)
                            size_in_dset = isize - start_notime[0];
                        else
                            size_in_dset = count_notime[0];
                        offset_in_dset = start_notime[0] - offsets[0];
                        offset_in_var = 0;
                    }
                }
                else {
                    // middle is in
                    if (isize < start_notime[0] + count_notime[0])
                        size_in_dset = ldims[0];
                    else
                    // tail is in
                        size_in_dset = count_notime[0] + start_notime[0] - offsets[0];
                    offset_in_dset = 0;
                    offset_in_var = offsets[0] - start_notime[0];
                }
    
                slice_size = size_in_dset * datasize * size_of_type;
                write_offset = offset_in_var * datasize * size_of_type;

                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset = var_root->characteristics[start_idx + idx].payload_offset 
                                 + offset_in_dset * datasize * size_of_type;
                    MPI_FILE_READ_OPS
                } else {
                    slice_offset = 0;
                    MPI_FILE_READ_OPS1
                }
    
                memcpy ((char *)data + write_offset, fh->b->buff + fh->b->offset, slice_size);
                if (fh->mfooter.change_endianness == adios_flag_yes) {
                    change_endianness((char *)data + write_offset, slice_size, var_root->type);
                }
    
                //write_offset +=  slice_size;
            }
            else 
            {

                uint64_t stride_offset = 0;
                int isize;
                uint64_t size_in_dset[10];
                uint64_t offset_in_dset[10];
                uint64_t offset_in_var[10];
                memset(size_in_dset, 0 , 10 * 8);
                memset(offset_in_dset, 0 , 10 * 8);
                memset(offset_in_var, 0 , 10 * 8);
                int hit = 0;
                for ( i = 0; i < ndim_notime ; i++) {
                    isize = offsets[i] + ldims[i];
                    if (start_notime[i] >= offsets[i]) {
                        // head is in
                        if (start_notime[i]<isize) {
                            if (start_notime[i] + count_notime[i] > isize)
                                size_in_dset[i] = isize - start_notime[i];
                            else
                                size_in_dset[i] = count_notime[i];
                            offset_in_dset[i] = start_notime[i] - offsets[i];
                            offset_in_var[i] = 0;
                            hit = 1 + hit * 10;
                        }
                        else
                            hit = -1;
                    }
                    else {
                        // middle is in
                        if (isize < start_notime[i] + count_notime[i]) {
                            size_in_dset[i] = ldims[i];
                            hit = 2 + hit * 10;
                        }
                        else {
                            // tail is in
                            size_in_dset[i] = count_notime[i] + start_notime[i] - offsets[i];
                            hit = 3 + hit * 10;
                        }
                        offset_in_dset[i] = 0;
                        offset_in_var[i] = offsets[i] - start_notime[i];
                    }
                }
    
                datasize = 1;
                var_stride = 1;
    
                for ( i = ndim_notime-1; i >= hole_break; i--) {
                    datasize *= size_in_dset[i];
                    dset_stride *= ldims[i];
                    var_stride *= count_notime[i];
                }
    
                uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
                for (i = ndim_notime - 1; i > -1; i--) {
                    start_in_payload += s * offset_in_dset[i] * size_of_type;
                    end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
                    s *= ldims[i];
                }
    
                slice_size = end_in_payload - start_in_payload + 1 * size_of_type;
    
                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset =  var_root->characteristics[start_idx + idx].payload_offset
                                  + start_in_payload;
                    MPI_FILE_READ_OPS
    
                    for ( i = 0; i < ndim_notime ; i++) {
                        offset_in_dset[i] = 0;
                    }
                } else {
                    slice_offset =  start_in_payload;
                    MPI_FILE_READ_OPS1
                }

                uint64_t var_offset = 0;
                uint64_t dset_offset = 0;
                for ( i = 0; i < hole_break; i++) {
                    nloop *= size_in_dset[i];
                }
    
                for ( i = 0; i < ndim_notime ; i++) {
                    var_offset = offset_in_var[i] + var_offset * count_notime[i];
                    dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
                }

                copy_data (data
                          ,fh->b->buff + fh->b->offset
                          ,0
                          ,hole_break
                          ,size_in_dset
                          ,ldims
                          ,count_notime
                          ,var_stride
                          ,dset_stride
                          ,var_offset
                          ,dset_offset
                          ,datasize
                          ,size_of_type 
                          );

            }
        }  // end for (idx ... loop over pgs
    
        free (idx_table);

        total_size += items_read * size_of_type;
        // shift target pointer for next read in
        data = (char *)data + (items_read * size_of_type);

    } // end for (timestep ... loop over timesteps

    free (dims);

    return total_size;
}

// Search for the start var index.
static int get_var_start_index (struct adios_index_var_struct_v1 * v, int t)
{
    int i = 0;

    while (i < v->characteristics_count) {
        if (v->characteristics[i].time_index == t) {
            return i;
        }

        i++;
    }

    return -1;
}

// Search for the stop var index
static int get_var_stop_index (struct adios_index_var_struct_v1 * v, int t)
{
    int i = v->characteristics_count - 1;

    while (i > -1) {
        if (v->characteristics[i].time_index == t) {
            return i;
        }

        i--;      
    }

    return -1;
}

int64_t adios_read_bp_read_var_byid2 (ADIOS_GROUP    * gp,
                                      int            varid,
                                      const uint64_t * start,
                                      const uint64_t * count,
                                      void           * data)
{
    struct BP_GROUP      * gh;
    struct BP_FILE       * fh;
    struct adios_index_var_struct_v1 * var_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    int    i,j,k, idx, t;
    int    start_time, stop_time;
    int    start_idx, stop_idx;
    int    ndim, ndim_notime, has_subfile, file_is_fortran;
    uint64_t size, * dims;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t datasize, nloop, dset_stride,var_stride, total_size=0, items_read;
    uint64_t count_notime[32], start_notime[32];
    int timedim = -1, temp_timedim, is_global = 0, size_of_type;
    uint64_t slice_offset, slice_size, tmpcount = 0;
    uint64_t datatimeoffset = 0; // offset in data to write a given timestep
    MPI_Status status;

    gh = (struct BP_GROUP *) gp->gh;
    fh = gh->fh;

    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);
    has_subfile = fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE;

    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i = 0; i< varid && var_root; i++) {
        var_root = var_root->next;
    }

    if (i!=varid) {
        adios_error (err_corrupted_variable, 
               "Variable id=%d is valid but was not found in internal data structures!",
               varid);
        return -adios_errno; 
    }

    /* Get dimensions and flip if caller != writer language */
    adios_read_bp_get_dimensions (var_root, fh->tidx_stop - fh->tidx_start + 1, file_is_fortran, 
                          &ndim, &dims, &timedim);

    /* In a Fortran written files, dimensions are in reversed order for C */
    //if ( file_is_fortran != futils_is_called_from_fortran() ) 
    if (file_is_fortran) 
        swap_order(ndim, dims, &timedim);

    /* Take out the time dimension from start[] and count[] */
    if (timedim == -1) {
        /* For timeless var, we still search from fh->tidx_start to fh->tidx_stop
           to handle the situation that some variables are dumped out in selected timesteps
        */
        start_time = fh->tidx_start;
        stop_time = fh->tidx_stop;

        for (i = 0; i < ndim; i++) {
             count_notime[i] = count[i];
             start_notime[i] = start[i];
        }
        ndim_notime = ndim;
    } else {
        j = 0;
        if (futils_is_called_from_fortran())
            temp_timedim = ndim - 1;
        else
            temp_timedim = 0;

        start_time = start[temp_timedim] + fh->tidx_start;
        stop_time = start_time + count[temp_timedim] - 1;

        for (i = 0; i < temp_timedim; i++) {
             count_notime[j] = count[i];
             start_notime[j] = start[i];
             j++;
        }
        i++; // skip timedim
        for (; i < ndim; i++) {
             count_notime[j] = count[i];
             start_notime[j] = start[i];
             j++;
        }
        ndim_notime = ndim-1;
    }

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if ( futils_is_called_from_fortran()) {
        swap_order(ndim_notime, count_notime, &timedim);
        swap_order(ndim_notime, start_notime, &timedim);
    }
    
    /* items_read = how many data elements are we going to read in total (per timestep) */
    items_read = 1;
    for (i = 0; i < ndim_notime; i++)
        items_read *= count_notime[i];
    
    size_of_type = bp_get_type_size (var_root->type, var_root->characteristics [0].value);

    /* For each timestep, do reading separately (they are stored in different sets of process groups */
    for (t = start_time; t <= stop_time; t++) {
        start_idx = get_var_start_index(var_root, t);
        stop_idx = get_var_stop_index(var_root, t);

        if (start_idx < 0 || stop_idx < 0) {
            adios_error (err_no_data_at_timestep,"Variable (id=%d) has no data at %d time step",
                varid, t);
//            return -adios_errno;
            continue;
        }

        if (ndim_notime == 0) {
            /* READ A SCALAR VARIABLE */
            slice_size = size_of_type;
            idx = 0; // macros below need it

            if (var_root->type == adios_string) {
                // strings are stored without \0 in file
                // size_of_type here includes \0 so decrease by one
                size_of_type--;
            }

            if (var_root->characteristics[start_idx+idx].payload_offset > 0) {
                slice_offset = var_root->characteristics[start_idx+idx].payload_offset;
                if (!has_subfile) {
                    MPI_FILE_READ_OPS
                } else {
                    MPI_FILE_READ_OPS2
                }
            } else {
                slice_offset = 0;
                MPI_FILE_READ_OPS1
            }

            memcpy((char *)data+total_size, fh->b->buff + fh->b->offset, size_of_type);
            if (fh->mfooter.change_endianness == adios_flag_yes) {
                change_endianness((char *)data+total_size, size_of_type, var_root->type);
            }

            if (var_root->type == adios_string) {
                // add \0 to the end of string
                // size_of_type here is the length of string
                // FIXME: how would this work for strings written over time?
                ((char*)data+total_size)[size_of_type] = '\0';
            }

            total_size += size_of_type;
            
            if (timedim == -1)
                break;
            else
                continue;
        }

         /* READ AN ARRAY VARIABLE */
        int * idx_table = (int *) malloc (sizeof(int) * (stop_idx - start_idx + 1));

        uint64_t write_offset = 0;
        int npg = 0;
        tmpcount = 0;
        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++) {
            int flag;
            datasize = 1;
            nloop = 1;
            var_stride = 1;
            dset_stride = 1;
            idx_table[idx] = 1;
            uint64_t payload_size = size_of_type;
    
            /* Each pg can have a different sized array, so we need the actual dimensions from it */
            is_global = adios_read_bp_get_dimensioncharacteristics( &(var_root->characteristics[start_idx + idx]),
                                                            ldims, gdims, offsets);
            if (!is_global) {
                // we use gdims below, which is 0 for a local array; set to ldims here
                for (j = 0; j< ndim; j++) {
                    gdims[j]=ldims[j];
                }
                // we need to read only the first PG, not all, so let's prevent a second loop
                stop_idx = start_idx;
            }

            /* Again, a Fortran written file has the dimensions in Fortran order we need to swap here */
            //if (file_is_fortran != futils_is_called_from_fortran()) {
            if (file_is_fortran) {
                i=-1;
                swap_order(ndim, gdims,   &(i)); // i is dummy 
                swap_order(ndim, ldims,   &(i));
                swap_order(ndim, offsets, &(i));
            }
            
            /* take out the time dimension */
            /* For C, gdims and offset are one size shorter because the timedim part is missing,
               so we take it out only for fortran files
            */
            if (timedim > -1) {
                for (i = timedim; i < ndim-1; i++) {
                    ldims[i] = ldims[i+1];
                    if (file_is_fortran) {
                        gdims[i] = gdims[i+1];
                        offsets[i] = offsets[i+1];
                    }
                }
            }

            
            printf("ldims   = "); for (j = 0; j<ndim; j++) printf("%d ",ldims[j]); printf("\n");
            printf("gdims   = "); for (j = 0; j<ndim; j++) printf("%d ",gdims[j]); printf("\n");
            printf("offsets = "); for (j = 0; j<ndim; j++) printf("%d ",offsets[j]); printf("\n");
            printf("count_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
            printf("start_notime   = "); for (j = 0; j<ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
            
                
            for (j = 0; j < ndim_notime; j++) {
    
                payload_size *= ldims [j];
    
                if ( (count_notime[j] > gdims[j]) 
                  || (start_notime[j] > gdims[j]) 
                  || (start_notime[j] + count_notime[j] > gdims[j])){
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j+1, count_notime[j], start_notime[j], gdims[j] - 1);
                    return -adios_errno;
                }
    
                /* check if there is any data in this pg and this dimension to read in */
                flag = (offsets[j] >= start_notime[j] 
                        && offsets[j] < start_notime[j] + count_notime[j])
                    || (offsets[j] < start_notime[j]
                        && offsets[j] + ldims[j] > start_notime[j] + count_notime[j]) 
                    || (offsets[j] + ldims[j] > start_notime[j] 
                        && offsets[j] + ldims[j] <= start_notime[j] + count_notime[j]);
                idx_table [idx] = idx_table[idx] && flag;
            }
            
            if ( !idx_table[idx] ) {
                continue;
            }
            ++npg;
    
            /* determined how many (fastest changing) dimensions can we read in in one read */
            int hole_break; 
            for (i = ndim_notime - 1; i > -1; i--) {
                if (offsets[i] == start_notime[i] && ldims[i] == count_notime[i]) {
                    datasize *= ldims[i];
                }
                else
                    break;
            }
    
            hole_break = i;
            slice_offset = 0;
            slice_size = 0;

            if (hole_break == -1) {
                /* The complete read happens to be exactly one pg, and the entire pg */
                /* This means we enter this only once, and npg=1 at the end */
                /* This is a rare case. FIXME: cannot eliminate this? */
                slice_size = payload_size;
    
                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset = var_root->characteristics[start_idx + idx].payload_offset;
                    if (!has_subfile) {
                        MPI_FILE_READ_OPS
                    } else {
                        MPI_FILE_READ_OPS2
                    }
                } else {
                    slice_offset = 0;
                    MPI_FILE_READ_OPS1
                }
    
                memcpy( (char *)data, fh->b->buff + fh->b->offset, slice_size);
                if (fh->mfooter.change_endianness == adios_flag_yes) {
                    change_endianness(data, slice_size, var_root->type);
                }
            }
            else if (hole_break == 0) 
            {
                /* The slowest changing dimensions should not be read completely but
                   we still need to read only one block */
                int isize;
                uint64_t size_in_dset = 0;
                uint64_t offset_in_dset = 0;
    
                isize = offsets[0] + ldims[0];
                if (start_notime[0] >= offsets[0]) {
                    // head is in
                    if (start_notime[0]<isize) {
                        if (start_notime[0] + count_notime[0] > isize)
                            size_in_dset = isize - start_notime[0];
                        else
                            size_in_dset = count_notime[0];
                        offset_in_dset = start_notime[0] - offsets[0];
                    }
                }
                else {
                    // middle is in
                    if (isize < start_notime[0] + count_notime[0])
                        size_in_dset = ldims[0];
                    else
                    // tail is in
                        size_in_dset = count_notime[0] + start_notime[0] - offsets[0];
                    offset_in_dset = 0;
                }
    
                slice_size = size_in_dset * datasize * size_of_type;
    
                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset = var_root->characteristics[start_idx + idx].payload_offset 
                                 + offset_in_dset * datasize * size_of_type;
                    if (!has_subfile) {
                        MPI_FILE_READ_OPS
                    } else {
                        MPI_FILE_READ_OPS2
                    }

                } else {
                    slice_offset = 0;
                    MPI_FILE_READ_OPS1
                }
    
                memcpy ((char *)data + write_offset, fh->b->buff + fh->b->offset, slice_size);
                if (fh->mfooter.change_endianness == adios_flag_yes) {
                    change_endianness((char *)data + write_offset, slice_size, var_root->type);
                }
    
                write_offset +=  slice_size;
            }
            else 
            {

                uint64_t stride_offset = 0;
                int isize;
                uint64_t size_in_dset[10];
                uint64_t offset_in_dset[10];
                uint64_t offset_in_var[10];
                memset(size_in_dset, 0 , 10 * 8);
                memset(offset_in_dset, 0 , 10 * 8);
                memset(offset_in_var, 0 , 10 * 8);
                int hit = 0;
                for ( i = 0; i < ndim_notime ; i++) {
                    isize = offsets[i] + ldims[i];
                    if (start_notime[i] >= offsets[i]) {
                        // head is in
                        if (start_notime[i]<isize) {
                            if (start_notime[i] + count_notime[i] > isize)
                                size_in_dset[i] = isize - start_notime[i];
                            else
                                size_in_dset[i] = count_notime[i];
                            offset_in_dset[i] = start_notime[i] - offsets[i];
                            offset_in_var[i] = 0;
                            hit = 1 + hit * 10;
                        }
                        else
                            hit = -1;
                    }
                    else {
                        // middle is in
                        if (isize < start_notime[i] + count_notime[i]) {
                            size_in_dset[i] = ldims[i];
                            hit = 2 + hit * 10;
                        }
                        else {
                            // tail is in
                            size_in_dset[i] = count_notime[i] + start_notime[i] - offsets[i];
                            hit = 3 + hit * 10;
                        }
                        offset_in_dset[i] = 0;
                        offset_in_var[i] = offsets[i] - start_notime[i];
                    }
                }
    
                datasize = 1;
                var_stride = 1;
    
                for ( i = ndim_notime-1; i >= hole_break; i--) {
                    datasize *= size_in_dset[i];
                    dset_stride *= ldims[i];
                    var_stride *= count_notime[i];
                }
    
                uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
                for (i = ndim_notime - 1; i > -1; i--) {
                    start_in_payload += s * offset_in_dset[i] * size_of_type;
                    end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
                    s *= ldims[i];
                }
    
                slice_size = end_in_payload - start_in_payload + 1 * size_of_type;
    
                if (var_root->characteristics[start_idx + idx].payload_offset > 0) {
                    slice_offset =  var_root->characteristics[start_idx + idx].payload_offset
                                  + start_in_payload;
                    if (!has_subfile) {
                        MPI_FILE_READ_OPS
                    } else {
                        MPI_FILE_READ_OPS2
                    }
 
                    for ( i = 0; i < ndim_notime ; i++) {
                        offset_in_dset[i] = 0;
                    }
                } else {
                    slice_offset =  start_in_payload;
                    MPI_FILE_READ_OPS1
                }
    
                uint64_t var_offset = 0;
                uint64_t dset_offset = 0;
                for ( i = 0; i < hole_break; i++) {
                    nloop *= size_in_dset[i];
                }
    
                for ( i = 0; i < ndim_notime ; i++) {
                    var_offset = offset_in_var[i] + var_offset * count_notime[i];
                    dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
                }
    
                copy_data (data
                          ,fh->b->buff + fh->b->offset
                          ,0
                          ,hole_break
                          ,size_in_dset
                          ,ldims
                          ,count_notime
                          ,var_stride
                          ,dset_stride
                          ,var_offset
                          ,dset_offset
                          ,datasize
                          ,size_of_type 
                          );
            }
        }  // end for (idx ... loop over pgs
    
        free (idx_table);
    
        total_size += items_read * size_of_type;
        // shift target pointer for next read in
        data = (char *)data + (items_read * size_of_type);

        if (timedim == -1)
            break;
    } // end for (timestep ... loop over timesteps

    free (dims);

    return total_size;
}

int64_t adios_read_bp_read_var_byid (ADIOS_GROUP    * gp,
                                     int            varid,
                                     const uint64_t  * start,
                                     const uint64_t  * count,
                                     void            * data)
{
    struct BP_GROUP      * gh;
    struct BP_FILE       * fh;
    int has_time_index_characteristic;

    adios_errno = 0;
    if (!gp) {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        return -adios_errno;
    }

    gh = (struct BP_GROUP *) gp->gh;
    if (!gh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return -adios_errno;
    }

    fh = gh->fh;
    if (!fh) {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return -adios_errno;
    }

    has_time_index_characteristic = fh->mfooter.version & ADIOS_VERSION_HAVE_TIME_INDEX_CHARACTERISTIC;
    if (!has_time_index_characteristic) {
        // read older file format. Can be deleted later on.
        return adios_read_bp_read_var_byid1(gp, varid, start, count, data);
    } else {
        return adios_read_bp_read_var_byid2(gp, varid, start, count, data);
    }
}

// NCSU - Timer series analysis, correlation
double adios_stat_cor (ADIOS_VARINFO * vix, ADIOS_VARINFO * viy, char * characteristic, uint32_t time_start, uint32_t time_end, uint32_t lag)
{
    int i,j;

    double avg_x = 0.0, avg_y = 0.0, avg_lag = 0.0;
    double var_x = 0.0, var_y = 0.0, var_lag = 0.0;
    double cov = 0;

    if (vix == NULL)
    {
        fprintf(stderr, "Variable not defined\n");
        return 0;
    }

    // If the vix and viy are not time series objects, return.
    if ((vix->timedim < 0) && (viy->timedim < 0))
    {             
        fprintf(stderr, "Covariance must involve timeseries data\n");
        return 0;
    }                                                                    

    uint32_t min = vix->dims[0] - 1;
    if (viy && (min > viy->dims[0] - 1))
        min = viy->dims[0] - 1;         
    
    if(time_start == 0 && time_end == 0) 
    { //global covariance
        if(viy == NULL) {
            fprintf(stderr, "Must have two variables for global covariance\n");
            return 0;
        }                                                                          

        // Assign vix to viy, and calculate covariance
        viy = vix;
        time_start = 0;
        time_end = min;
    }
    // Check the bounds of time
    if (    (time_start >= 0) && (time_start <= min)
            &&      (time_end >= 0)   && (time_end <= min)
            &&  (time_start <= time_end))
    {
        if(viy == NULL) //user must want to run covariance against itself
        {
            if(! (time_end+lag) > min)
            {                                                                        
                fprintf(stderr, "Must leave enough timesteps for lag\n");
                return 0;
            }

            if (strcmp(characteristic, "average") == 0 || strcmp(characteristic, "avg") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (adios_double, vix->avgs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (adios_double, vix->avgs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (adios_double, vix->avgs[i]); 
                    double val_lag = bp_value_to_double (adios_double, vix->avgs[i + lag]); 
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1); 
                    var_lag += (val_lag - avg_lag) * (val_lag - avg_lag) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_lag - avg_lag) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "standard deviation") == 0 || strcmp(characteristic, "std_dev") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (adios_double, vix->std_devs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (adios_double, vix->std_devs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (adios_double, vix->std_devs[i]);
                    double val_lag = bp_value_to_double (adios_double, vix->std_devs[i + lag]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_lag += (val_lag - avg_lag) * (val_lag - avg_lag) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_lag - avg_lag) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "minimum") == 0 || strcmp(characteristic, "min") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (vix->type, vix->mins[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (vix->type, vix->mins[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (vix->type, vix->mins[i]); 
                    double val_lag = bp_value_to_double (vix->type, vix->mins[i + lag]); 
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1); 
                    var_lag += (val_lag - avg_lag) * (val_lag - avg_lag) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_lag - avg_lag) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "maximum") == 0 || strcmp(characteristic, "max") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (vix->type, vix->maxs[i]); 
                    double val_lag = bp_value_to_double (vix->type, vix->maxs[i + lag]); 
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1); 
                    var_lag += (val_lag - avg_lag) * (val_lag - avg_lag) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_lag - avg_lag) / (time_end - time_start + 1);
                }
            }
            else
            {
                fprintf (stderr, "Unknown characteristic\n");
                return 0;
            }
            return cov / (sqrt (var_x) * sqrt (var_lag));
        }
        else
        {
            if (strcmp(characteristic, "average") == 0 || strcmp(characteristic, "avg") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(adios_double, vix->avgs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(adios_double, viy->avgs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (adios_double, vix->avgs[i]); 
                    double val_y = bp_value_to_double (adios_double, viy->avgs[i]); 
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1); 
                    var_y += (val_y - avg_y) * (val_y - avg_y) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_y - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "standard deviation") == 0 || strcmp(characteristic, "std_dev") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(adios_double, vix->std_devs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(adios_double, viy->std_devs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (adios_double, vix->std_devs[i]);
                    double val_y = bp_value_to_double (adios_double, viy->std_devs[i]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_y += (val_y - avg_y) * (val_y - avg_y) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_y - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "minimum") == 0 || strcmp(characteristic, "min") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(vix->type, vix->mins[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(viy->type, viy->mins[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (vix->type, vix->mins[i]); 
                    double val_y = bp_value_to_double (viy->type, viy->mins[i]); 
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1); 
                    var_y += (val_y - avg_y) * (val_y - avg_y) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_y - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "maximum") == 0 || strcmp(characteristic, "max") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(vix->type, viy->maxs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (vix->type, vix->maxs[i]); 
                    double val_y = bp_value_to_double (viy->type, viy->maxs[i]); 
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1); 
                    var_y += (val_y - avg_y) * (val_y - avg_y) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_y - avg_y) / (time_end - time_start + 1);
                }
            }
            else
            {
                fprintf (stderr, "Unknown characteristic\n");
                return 0;
            }
            return cov / (sqrt (var_x) * sqrt (var_y));
        }
    }
    else
    {
        fprintf (stderr, "Time values out of bounds\n");
        return 0;
    }
}

// NCSU - Time series analysis, covariance
//covariance(x,y) = sum(i=1,..N) [(x_1 - x_mean)(y_i - y_mean)]/N
double adios_stat_cov (ADIOS_VARINFO * vix, ADIOS_VARINFO * viy, char * characteristic, uint32_t time_start, uint32_t time_end, uint32_t lag)
{
    int i,j;

    double avg_x = 0.0, avg_y = 0.0, avg_lag = 0.0;
    double cov = 0;

    if (vix == NULL)
    {
        fprintf(stderr, "Variable not defined\n");
        return 0;
    }

    // If the vix and viy are not time series objects, return.
    if ((vix->timedim < 0) && (viy->timedim < 0))
    {             
        fprintf(stderr, "Covariance must involve timeseries data\n");
        return 0;
    }                                                                    

    uint32_t min = vix->dims[0] - 1;
    if (viy && (min > viy->dims[0] - 1))
        min = viy->dims[0] - 1;         
    
    if(time_start == 0 && time_end == 0) 
    { //global covariance
        if(viy == NULL) {
            fprintf(stderr, "Must have two variables for global covariance\n");
            return 0;
        }                                                                          

        // Assign vix to viy, and calculate covariance
        viy = vix;
        time_start = 0;
        time_end = min;
    }
    // Check the bounds of time
    if (    (time_start >= 0) && (time_start <= min)
            &&      (time_end >= 0)   && (time_end <= min)
            &&  (time_start <= time_end))
    {
        if(viy == NULL) //user must want to run covariance against itself
        {
            if(! (time_end+lag) > min)
            {                                                                        
                fprintf(stderr, "Must leave enough timesteps for lag\n");
                return 0;
            }

            if (strcmp(characteristic, "average") == 0 || strcmp(characteristic, "avg") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (adios_double, vix->avgs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (adios_double, vix->avgs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                    cov += (bp_value_to_double (adios_double, vix->avgs[i]) - avg_x) * (bp_value_to_double (adios_double, vix->avgs[i+lag]) - avg_lag) / (time_end - time_start + 1);
            }
            else if (strcmp(characteristic, "standard deviation") == 0 || strcmp(characteristic, "std_dev") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (adios_double, vix->std_devs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (adios_double, vix->std_devs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                    cov += (bp_value_to_double (adios_double, vix->std_devs[i]) - avg_x) * (bp_value_to_double (adios_double, vix->std_devs[i+lag]) - avg_lag) / (time_end - time_start + 1);
            }
            else if (strcmp(characteristic, "minimum") == 0 || strcmp(characteristic, "min") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (vix->type, vix->mins[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (vix->type, vix->mins[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                    cov += (bp_value_to_double (vix->type, vix->mins[i]) - avg_x) * (bp_value_to_double (vix->type, vix->mins[i+lag]) - avg_lag) / (time_end - time_start + 1);
            }
            else if (strcmp(characteristic, "maximum") == 0 || strcmp(characteristic, "max") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (vix->type, vix->maxs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                    cov += (bp_value_to_double (vix->type, vix->maxs[i]) - avg_x) * (bp_value_to_double (vix->type, vix->maxs[i+lag]) - avg_lag) / (time_end - time_start + 1);
            }
            else
            {
                fprintf (stderr, "Unknown characteristic\n");
                return 0;
            }
        }
        else
        {
            if (strcmp(characteristic, "average") == 0 || strcmp(characteristic, "avg") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(adios_double, vix->avgs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(adios_double, viy->avgs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    cov += (bp_value_to_double(adios_double, vix->avgs[i]) - avg_x) * (bp_value_to_double(adios_double, viy->avgs[i]) - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "standard deviation") == 0 || strcmp(characteristic, "std_dev") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(adios_double, vix->std_devs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(adios_double, viy->std_devs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    cov += (bp_value_to_double(adios_double, vix->std_devs[i]) - avg_x) * (bp_value_to_double(adios_double, viy->std_devs[i]) - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "minimum") == 0 || strcmp(characteristic, "min") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(vix->type, vix->mins[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(viy->type, viy->mins[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    cov += (bp_value_to_double(vix->type, vix->mins[i]) - avg_x) * (bp_value_to_double(viy->type, viy->mins[i]) - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "maximum") == 0 || strcmp(characteristic, "max") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(vix->type, viy->maxs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    cov += (bp_value_to_double(vix->type, vix->maxs[i]) - avg_x) * (bp_value_to_double(viy->type, viy->maxs[i]) - avg_y) / (time_end - time_start + 1);
                }
            }
            else
            {
                fprintf (stderr, "Unknown characteristic\n");
                return 0;
            }
        }
    }
    else
    {
        fprintf (stderr, "Time values out of bounds\n");
    }
    return cov;
}
