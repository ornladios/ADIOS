/*
    read_xpmem.c       
    Goal: to create evpath io connection layer in conjunction with 
    write/adios_xpmem.g

*/
// system libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/queue.h>
#include <sys/socket.h>
#include <sys/times.h>
#include <netinet/in.h>
#include <sys/time.h>
#include <sys/uio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>
#include <assert.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

// local libraries
#include "config.h"
#include "public/adios.h"
#include "public/adios_types.h"
#include "public/adios_read_v2.h"
#include "core/adios_read_hooks.h"
#include "core/adios_logger.h"
#include "core/common_read.h"
#include "public/adios_error.h"
#include "core/bp_types.h"
#include "core/adios_endianness.h"

//include the bp utils header to get access to
//all the usable functions

#include "core/bp_types.h"
#include "core/bp_utils.h"

#include <xpmem.h>

#include "public/adios_xpmem.h"
#include "read/read_xpmem.h"

inline BP_PROC * GET_BP_PROC (const ADIOS_FILE * fp)
{
	return (BP_PROC*)((xpmem_read_file*)fp->fh)->bp;
}

inline BP_FILE * GET_BP_FILE (const ADIOS_FILE * fp)
{
    return (BP_FILE *) ((xpmem_read_file *) fp->fh)->fh;
}


static int map_req_varid (const ADIOS_FILE * fp, int varid)
{
    BP_PROC * p = GET_BP_PROC (fp);

    return p->varid_mapping[varid];
}

/********** Core ADIOS Read functions. **********/

/*
 * Gathers basic MPI information; sets up reader CM.
 */
int
adios_read_xpmem_init_method (MPI_Comm comm, PairStruct* params)
{
	fp = (xpmem_read_data*)malloc(sizeof(xpmem_read_data));
	
	char *buffer;
	char *index;
	
	//read the segid
	memset(fp, 0, sizeof(xpmem_read_data));

	read_segid(&fp->data_segid, "xpmem.data");
	buffer = attach_segid(fp->data_segid, share_size, &fp->data_apid);

	while(buffer == NULL)
	{
		sleep(1);
		read_segid(&fp->data_segid, "xpmem.data");
		buffer = attach_segid(fp->data_segid, share_size, &fp->data_apid);
	}
		
	fprintf(stderr, "opened the segments %p \n", buffer);

	fp->pg = (shared_data*)buffer;

	log_info("read init completed\n");
	
	return 0;
}

ADIOS_FILE*
adios_read_xpmem_open_file(const char * fname, MPI_Comm comm)
{
	log_info("read open file\n");
	
	ADIOS_FILE *af = (ADIOS_FILE*)malloc(sizeof(ADIOS_FILE));
	if(!af){
		adios_error (err_no_memory, 
		             "Cannot allocate memory for file info.\n");
		return NULL;
	}

	xpmem_read_file *f = (xpmem_read_file*)malloc(sizeof(xpmem_read_file));

	//allocate the fake bp file
	f->fh = (BP_FILE*) malloc(sizeof(BP_FILE));
	f->bp = (BP_PROC*) malloc(sizeof(BP_PROC));

	//zero out the structures
	memset(f->bp, 0, sizeof(BP_PROC));
	memset(f->fh, 0, sizeof(BP_FILE));

	//set the read data to the file

	f->bp->fh = f->fh;
	f->fp = fp;

	//initialize some of the stuff in the BP_FILE
	
    af->fh = (uint64_t)f; 
	af->path = strdup("xpmem.data");
	
	//at this point af->fh is the xpmem_read_file
	//af->fh->fh is the BP_FILE
	//af->fp is the xpmem_read_data
	


	log_info("spinning on version\n");

	//we will spin until the file is version is updated
	//incidating that we have some data here
	
	while(f->fp->pg->version == 0)
	    adios_nanosleep(0, 100000000);

	//now the buffer has some data
	f->fp->data = (char*)malloc(f->fp->pg->size);
	f->fp->dsize = f->fp->pg->size;
	

	//copy the data into a non-shared buffer of the right size
	//this is equivalent to the pg
	memcpy(f->fp->data, f->fp->pg->buffer, f->fp->dsize);
	memcpy(&f->fp->debug, f->fp->pg, sizeof(shared_data));

	//read the data
	xp_read_open(f->fh, f->fp);

	xp_seek_to_step(af, -1, show_hidden_attrs);

    af->endianness =  bp_get_endianness (f->fh->mfooter.change_endianness);
	af->version =  f->fh->mfooter.version & ADIOS_VERSION_NUM_MASK;
	af->current_step = 0;
	af->last_step = 0;

	//just check to make sure that readcount is 0
	log_debug("xpmem readcount = %d\n",
	          f->fp->pg->readcount);

	//now fp->index contains the index
	//index size is fp-index->size
	//copy the buffer out
	
	return af;  
}

/*
 * Still have work to do here.  
 * Change it so that we can support the timeouts and lock_modes.
 */
/*
 * Sets up local data structure for series of reads on an adios file
 * - create evpath graph and structures
 * -- create evpath control stone (outgoing)
 * -- create evpath data stone (incoming)
 * -- rank 0 dumps contact info to file
 * -- create connections using contact info from file
 */
ADIOS_FILE*
adios_read_xpmem_open(const char * fname,
			 MPI_Comm comm,
			 enum ADIOS_LOCKMODE lock_mode,
			 float timeout_sec)
{

	log_info("in read_open. xpmem treats everything like a file right now");
	return adios_read_xpmem_open_file(fname, comm);
}

int adios_read_xpmem_finalize_method ()
{
	
    return 0;
}

void adios_read_xpmem_release_step(ADIOS_FILE *fp)
{

	xpmem_read_file *xf = (xpmem_read_file*)fp->fh;
	xpmem_read_data *xd = xf->fp;

	// if(fp->pg->version != 0)
	// 	fp->pg->version = 0;
	// if(fp->index->version != 0)
	// 	fp->index->version = 0;

	if(xd->pg->readcount < 1)
		xd->pg->readcount = 1;

}

int 
adios_read_xpmem_advance_step(ADIOS_FILE *fp, int last, float timeout_sec) 
{
	xpmem_read_file *xf = (xpmem_read_file*)fp->fh;
	xpmem_read_data *xd = xf->fp;

	// if(fp->pg->version != 0)
	// 	fp->pg->version = 0;
	// if(fp->index->version != 0)
	// 	fp->index->version = 0;

	if(xd->pg->readcount < 1)
		xd->pg->readcount = 1;
    return 0;
}

int adios_read_xpmem_close(ADIOS_FILE * fp)
{
	BP_PROC * p = GET_BP_PROC (fp);
	BP_FILE * fh = GET_BP_FILE (fp);
	xpmem_read_file *xf = (xpmem_read_file*)fp->fh;
	xpmem_read_data *xd = xf->fp;

	if (p->fh)
	{
		bp_close (fh);
		p->fh = 0;
	}

	if (p->varid_mapping)
	{
		free (p->varid_mapping);
		p->varid_mapping = 0;
	}

	if (p->local_read_request_list)
	{
		list_free_read_request (p->local_read_request_list);
		p->local_read_request_list = 0;
	}

	free (p);

	if (fp->var_namelist)
	{
		free_namelist (fp->var_namelist, fp->nvars);
		fp->var_namelist = 0;
	}

	if (fp->attr_namelist)
	{
		free_namelist (fp->attr_namelist, fp->nattrs);
		fp->attr_namelist = 0;
	}

	if (fp->path)
	{
		free (fp->path);
		fp->path = 0;
	}
	// internal_data field is taken care of by common reader layer
	free (fp);

	

	return 0;
}

ADIOS_FILE *adios_read_xpmem_fopen(const char *fname, MPI_Comm comm) {
   return 0;
}

int adios_read_xpmem_is_var_timed(const ADIOS_FILE* fp, int varid) { return 0; }

void adios_read_xpmem_get_groupinfo(
    const ADIOS_FILE *fp, 
    int *ngroups, 
    char ***group_namelist, 
    uint32_t **nvars_per_group, 
    uint32_t **nattrs_per_group) 
{
	xpmem_read_file *xf = (xpmem_read_file*)fp->fh;
	BP_FILE *fh = xf->fh;
	BP_PROC *p = xf->bp;
	
    int i, j, offset;

    * ngroups = fh->gvar_h->group_count;

    *group_namelist = (char **) malloc (sizeof (char *) * fh->gvar_h->group_count);
    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        (*group_namelist)[i] = malloc (strlen (fh->gvar_h->namelist[i]) + 1);
        assert ((*group_namelist)[i]);

        memcpy ((*group_namelist)[i], fh->gvar_h->namelist[i], strlen (fh->gvar_h->namelist[i]) + 1);
    }

    * nvars_per_group = (uint32_t *) malloc (fh->gvar_h->group_count * sizeof (uint32_t));
    assert (* nvars_per_group);

    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        (* nvars_per_group)[i] = fh->gvar_h->var_counts_per_group[i];
    }

    * nattrs_per_group = (uint32_t *) malloc (fh->gattr_h->group_count * sizeof (uint32_t));
    assert (* nattrs_per_group);

    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        offset = 0;
        for (j = 0; j < i; j++)
        {
            offset += fh->gattr_h->attr_counts_per_group[j];
        }

        (* nattrs_per_group)[i] = 0;
        for (j = 0; j < fh->gattr_h->attr_counts_per_group[i]; j++)
        {
            if (!show_hidden_attrs && strstr (fh->gattr_h->attr_namelist[offset + j], "__adios__"))
            {
            }
            else
            {
                (* nattrs_per_group)[i] ++;
            }
        }
    }

    return;
}


int adios_read_xpmem_check_reads(const ADIOS_FILE* fp, ADIOS_VARCHUNK** chunk) { log_debug( "xpmem:adios function check reads\n"); return 0; }

int adios_read_xpmem_perform_reads(const ADIOS_FILE *fp, int blocking)
{
	BP_PROC * p = GET_BP_PROC (fp);
	read_request * r;
	ADIOS_VARCHUNK * chunk;

	/* 1. prepare all reads */
	// check if all user memory is provided for blocking read
	if (blocking)
	{
		r = p->local_read_request_list;
		while (r)
		{
			if (!r->data)
			{
				adios_error (err_operation_not_supported,
				             "Blocking mode at adios_perform_reads() requires that user "
				             "provides the memory for each read request. Request for "
				             "variable %d was scheduled without user-allocated memory\n",
				             r->varid);
				return err_operation_not_supported;
			}

			r = r->next;
		}
	}
	else
	{
		adios_error(err_operation_not_supported, "xpmem does not support non-blocking reads\n");
		return adios_errno;
	}

	//read the variable
	while (p->local_read_request_list)
    {
        chunk = read_var (fp, p->local_read_request_list);

        // remove head from list
        r = p->local_read_request_list;
        p->local_read_request_list = p->local_read_request_list->next;
        common_read_selection_delete (r->sel);
        r->sel = NULL;
        free(r);

        common_read_free_chunk (chunk);
    }

	return 0;
}

int
adios_read_xpmem_inq_var_blockinfo(const ADIOS_FILE* fp,
				      ADIOS_VARINFO* varinfo)
{ log_debug( "xpmem:adios function inq var block info\n"); return 0; }

int
adios_read_xpmem_inq_var_stat(const ADIOS_FILE* fp,
				 ADIOS_VARINFO* varinfo,
				 int per_step_stat,
				 int per_block_stat)
{ log_debug( "xpmem:adios function inq var stat\n"); return 0; }


int 
adios_read_xpmem_schedule_read_byid(const ADIOS_FILE *fp,
				       const ADIOS_SELECTION *sel,
				       int varid,
				       int from_steps,
				       int nsteps,
				       void *data)
{
	log_debug("xpmem:adios function schedule read byid\n");
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);

    read_request * r;
    ADIOS_SELECTION * nullsel = 0;
    struct adios_index_var_struct_v1 * v;
    int i, ndim, ns, file_is_fortran, mapped_varid;
    uint64_t * dims = 0;

    mapped_varid = p->varid_mapping[varid];
    v = bp_find_var_byid (fh, mapped_varid);
    file_is_fortran = is_fortran_file (fh);

    r = (read_request *) malloc (sizeof (read_request));
    assert (r);

    if (!sel)
    {
        bp_get_and_swap_dimensions (fp, v, file_is_fortran,
                                    &ndim, &dims,
                                    &ns,
                                    file_is_fortran != futils_is_called_from_fortran()
                                    );

        nullsel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
        assert (nullsel);

        nullsel->type = ADIOS_SELECTION_BOUNDINGBOX;
        nullsel->u.bb.ndim = ndim;
        nullsel->u.bb.start = (uint64_t *) malloc (nullsel->u.bb.ndim * 8);
        assert (nullsel->u.bb.start);
        nullsel->u.bb.count = (uint64_t *) malloc (nullsel->u.bb.ndim * 8);
        assert (nullsel->u.bb.count);

        for (i = 0; i < nullsel->u.bb.ndim; i++)
        {
            nullsel->u.bb.start[i] = 0;
            nullsel->u.bb.count[i] = dims[i];
        }

        free (dims);
    }

    /* copy selection since we don't want to operate on user memory.
     */
    r->sel = (!nullsel ? copy_selection (sel) : nullsel);
    r->varid = mapped_varid;
    if (!p->streaming)
    {
        r->from_steps = from_steps;
        r->nsteps = nsteps;
    }
    else
    {
        r->from_steps = 0;
        r->nsteps = 1;
    }

    r->data = data;
    r->datasize = get_req_datasize (fp, r, v);
    r->priv = 0;
    r->next = 0;

    list_insert_read_request_next (&p->local_read_request_list, r);

    return 0;
}


int 
adios_read_xpmem_get_attr_byid (const ADIOS_FILE *fp, int attrid,
				   enum ADIOS_DATATYPES *type,
				   int *size, void **data)
{
    int i;
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);
    struct adios_index_attribute_struct_v1 * attr_root;
    struct adios_index_var_struct_v1 * var_root, * v1;
    int file_is_fortran, last_step = fp->last_step, show_hidden_attrs;
    uint64_t k, attr_c_index, var_c_index;

    adios_errno = 0;

    show_hidden_attrs = 0;
    for (i = 0; i < fp->nattrs; i++)
    {
        if (strstr (fp->attr_namelist[i], "__adios__"))
        {
            show_hidden_attrs = 1;
            break;
        }
    }

    attr_root = fh->attrs_root; /* need to traverse the attribute list of the group */
    i = 0;

    if (show_hidden_attrs)
    {
        while (i < attrid && attr_root)
        {
            i++;
            attr_root = attr_root->next;
        }
    }
    else
    {
        while (i < attrid && attr_root)
        {
            if (strstr (attr_root->attr_path, "__adios__"))
            {
            }
            else
            {
                i++;
            }

            attr_root = attr_root->next;
        }

        while (attr_root && strstr (attr_root->attr_path, "__adios__"))
        {
            attr_root = attr_root->next;
        }
    }

    assert (attr_root);

    if (i != attrid)
    {
        adios_error (err_corrupted_attribute, "Attribute id=%d is valid but was not found in internal data structures!\n",attrid);
        return adios_errno;
    }

    /* Look for the last step because some of the hidden attributes, such as last update time,
     * make sense for the most recent value. 07/2011 - Q.Liu
     */

    attr_c_index = -1;
    for (k = 0; k < attr_root->characteristics_count; k++)
    {
        if (attr_root->characteristics[k].time_index - 1 == last_step)
        {
            attr_c_index = k;
            break;
        }
    }

    if (attr_c_index == -1)
    {
        log_debug ("adios_read_xpmemget_attr_byid: cannot find step : %d\n", last_step);
        attr_c_index = 0;
    }

    file_is_fortran = is_fortran_file (fh);

    // check the last version
    if (attr_root->characteristics[attr_c_index].value)
    {
        /* Attribute has its own value */
        *size = bp_get_type_size (attr_root->type, attr_root->characteristics[attr_c_index].value);
        *type = attr_root->type;
        *data = (void *) malloc (*size);
        assert (*data);

        memcpy(*data, attr_root->characteristics[attr_c_index].value, *size);
    }
    else if (attr_root->characteristics[attr_c_index].var_id)
    {
        /* Attribute is a reference to a variable */
        /* FIXME: var ids are not unique in BP. If a group of variables are written several
           times under different path using adios_set_path(), the id of a variable is always
           the same (should be different). As a temporary fix, we look first for a matching
           id plus path between an attribute and a variable. If not found, then we look for
           a match on the ids only.*/
        var_root = fh->vars_root;
        while (var_root)
        {
            if (var_root->id == attr_root->characteristics[attr_c_index].var_id
               && !strcmp(var_root->var_path, attr_root->attr_path)
               && !strcmp(var_root->group_name, attr_root->group_name)
               )
                break;
            var_root = var_root->next;
        }

        if (!var_root)
        {
            var_root = fh->vars_root;
            while (var_root)
            {
                if (var_root->id == attr_root->characteristics[attr_c_index].var_id
                   && !strcmp(var_root->group_name, attr_root->group_name))
                    break;
                var_root = var_root->next;
            }
        }

        if (!var_root)
        {
            var_root = fh->vars_root;
            while (var_root)
            {
                if (var_root->id == attr_root->characteristics[attr_c_index].var_id)
                    break;
                var_root = var_root->next;
            }
        }

        if (!var_root)
        {
            adios_error (err_invalid_attribute_reference,
                   "Attribute %s/%s in group %s is a reference to variable ID %d, which is not found\n",
                   attr_root->attr_path, attr_root->attr_name, attr_root->group_name,
                   attr_root->characteristics[attr_c_index].var_id);
            return adios_errno;
        }

        /* default values in case of error */
        *data = NULL;
        *size = 0;
        *type = attr_root->type;

        var_c_index = -1;
        for (k = 0; k < var_root->characteristics_count; k++)
        {
            if (var_root->characteristics[k].time_index - 1 == last_step)
            {
                var_c_index = k;
                break;
            }
        }

        if (var_c_index == -1)
        {
            var_c_index = 0;
            log_debug ("adios_read_xpmem_get_attr_byid: cannot find step : %d\n", last_step);
        }
        /* FIXME: variable and attribute type may not match, then a conversion is needed. */
        /* Cases:
                1. attr has no type, var is byte array     ==> string
                2. attr has no type, var is not byte array ==> var type
                3. attr is string, var is byte array       ==> string
                4. attr type == var type                   ==> var type
                5. attr type != var type                   ==> attr type and conversion needed
        */
        /* Error check: attr cannot reference an array in general */
        if (var_root->characteristics[var_c_index].dims.count > 0)
        {
            if ( (var_root->type == adios_byte || var_root->type == adios_unsigned_byte) &&
                 (attr_root->type == adios_unknown || attr_root->type == adios_string) &&
                 (var_root->characteristics[var_c_index].dims.count == 1))
            {
                 ; // this conversions are allowed
            }
            else
            {
                adios_error (err_invalid_attribute_reference,
                    "Attribute %s/%s in group %s, typeid=%d is a reference to an %d-dimensional array variable "
                    "%s/%s of type %s, which is not supported in ADIOS\n",
                    attr_root->attr_path, attr_root->attr_name, attr_root->group_name, attr_root->type,
                    var_root->characteristics[var_c_index].dims.count,
                    var_root->var_path, var_root->var_name, common_read_type_to_string(var_root->type));
                return adios_errno;
            }
        }

        if ( (attr_root->type == adios_unknown || attr_root->type == adios_string) &&
             (var_root->type == adios_byte || var_root->type == adios_unsigned_byte) &&
             (var_root->characteristics[var_c_index].dims.count == 1) )
        {
            /* 1D byte arrays are converted to string */
            /* 1. read in variable */
            char varname[512];
            char *tmpdata;
            ADIOS_VARCHUNK *vc;
            read_request * r;
            uint64_t start, count;
            int varid = 0;
            v1 = fh->vars_root;
            while (v1 && v1 != var_root)
            {
                v1 = v1->next;
                varid++;
            }

            start = 0;
            count = var_root->characteristics[var_c_index].dims.dims[0];
            snprintf(varname, 512, "%s/%s", var_root->var_path, var_root->var_name);
            tmpdata = (char *) malloc (count+1);
            assert (tmpdata);

            r = (read_request *) malloc (sizeof (read_request));
            assert (r);

            r->sel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
            r->sel->type = ADIOS_SELECTION_BOUNDINGBOX;
            r->sel->u.bb.ndim = 1;
            r->sel->u.bb.start = &start;
            r->sel->u.bb.count = &count;
            r->varid = varid;
            r->from_steps = fp->last_step;
            r->nsteps = 1;
            r->data = tmpdata;
            r->datasize = count;
            r->priv = 0;
            r->next = 0;

            vc = read_var_bb (fp, r);

            free (r->sel);
            free (r);

            if (vc == 0)
            {
                char *msg = strdup(adios_get_last_errmsg());
                adios_error ((enum ADIOS_ERRCODES) adios_errno,
                      "Cannot read data of variable %s/%s for attribute %s/%s of group %s: %s\n",
                      var_root->var_path, var_root->var_name,
                      attr_root->attr_path, attr_root->attr_name, attr_root->group_name,
                      msg);
                free(tmpdata);
                free(msg);
                return adios_errno;
            }

            *type = adios_string;
            if (file_is_fortran)
            {
                /* Fortran byte array to C string */
	            *data = (void*)futils_fstr_to_cstr( tmpdata, (int)count); /* FIXME: supports only 2GB strings... */
                *size = strlen( (char *)data );
                free(tmpdata);
            }
            else
            {
                /* C array to C string */
                tmpdata[count] = '\0';
                *size = count+1;
                *data = tmpdata;
            }

            free (vc->sel);
            free (vc);
        }
        else
        {
            /* other types are inherited */
            *type = var_root->type;
            *size = bp_get_type_size (var_root->type, var_root->characteristics[var_c_index].value);
            *data = (void *) malloc (*size);
            assert (*data);
            memcpy(*data, var_root->characteristics[var_c_index].value, *size);
        }
    }

    return 0;
}

ADIOS_VARINFO* 
adios_read_xpmem_inq_var_byid (const ADIOS_FILE * fp, int varid)
{
	ADIOS_VARINFO *varinfo;
	int mapped_id = map_req_varid(fp, varid);
	adios_errno = 0;

	varinfo = bp_inq_var_byid(fp, mapped_id);
	varinfo->varid = varid;
	
    return varinfo;
}

void 
adios_read_xpmem_free_varinfo (ADIOS_VARINFO *adiosvar)
{
    //log_debug( "debug: adios_read_xpmem_free_varinfo\n");
    fprintf(stderr, "adios_read_xpmem_free_varinfo called\n");
    return;
}


ADIOS_TRANSINFO* 
adios_read_xpmem_inq_var_transinfo(const ADIOS_FILE *gp, const ADIOS_VARINFO *vi)
{    
    return NULL;
}


int 
adios_read_xpmem_inq_var_trans_blockinfo(const ADIOS_FILE *gp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti)
{
    adios_error(err_operation_not_supported, "Xpmem does not yet support transforms: trans_blockinfo.\n");
    return (int64_t)0;
}

void 
adios_read_xpmem_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{
    //log_debug( "debug: adios_read_xpmem_reset_dimension_order\n");
    adios_error(err_invalid_read_method, "adios_read_xpmem_reset_dimension_order is not implemented.");
}

