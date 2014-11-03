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
	af->current_step = 0;
	af->last_step = 0;
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

	log_info("in read_open");
	return adios_read_xpmem_open_file(fname, comm);
}

int adios_read_xpmem_finalize_method ()
{
    return 0;
}

void adios_read_xpmem_release_step(ADIOS_FILE *adiosfile)
{

	xpmem_read_file *f = (xpmem_read_file*)adiosfile->fh;
	xpmem_read_data *fp = f->fp;

	// if(fp->pg->version != 0)
	// 	fp->pg->version = 0;
	// if(fp->index->version != 0)
	// 	fp->index->version = 0;

	if(fp->pg->readcount < 1)
		fp->pg->readcount = 1;

}

int 
adios_read_xpmem_advance_step(ADIOS_FILE *adiosfile, int last, float timeout_sec) 
{
	xpmem_read_file *f = (xpmem_read_file*)adiosfile->fh;
	xpmem_read_data *fp = f->fp;

	// if(fp->pg->version != 0)
	// 	fp->pg->version = 0;
	// if(fp->index->version != 0)
	// 	fp->index->version = 0;

	if(fp->pg->readcount < 1)
		fp->pg->readcount = 1;
    return 0;
}

int adios_read_xpmem_close(ADIOS_FILE * fp)
{
    return 0;
}

ADIOS_FILE *adios_read_xpmem_fopen(const char *fname, MPI_Comm comm) {
   return 0;
}

int adios_read_xpmem_is_var_timed(const ADIOS_FILE* fp, int varid) { return 0; }

void adios_read_xpmem_get_groupinfo(
    const ADIOS_FILE *adiosfile, 
    int *ngroups, 
    char ***group_namelist, 
    uint32_t **nvars_per_group, 
    uint32_t **nattrs_per_group) 
{
	xpmem_read_file *xf = (xpmem_read_file*)adiosfile->fh;
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

int adios_read_xpmem_perform_reads(const ADIOS_FILE *adiosfile, int blocking)
{
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
adios_read_xpmem_schedule_read(const ADIOS_FILE *adiosfile,
			const ADIOS_SELECTION * sel,
			const char * varname,
			int from_steps,
			int nsteps,
			void * data)
{
	log_debug("xpmem:adios function schedule read\n");
    return 0;
}

int 
adios_read_xpmem_get_attr (int *gp, const char *attrname,
                                 enum ADIOS_DATATYPES *type,
                                 int *size, void **data)
{
	log_debug("xpmem:adios function get attr\n");
    return adios_errno;
}

int 
adios_read_xpmem_get_attr_byid (const ADIOS_FILE *adiosfile, int attrid,
				   enum ADIOS_DATATYPES *type,
				   int *size, void **data)
{
	log_debug("xpmem:adios function get attr by id\n");
    return adios_errno;
}

ADIOS_VARINFO* 
adios_read_xpmem_inq_var(const ADIOS_FILE * adiosfile, const char* varname)
{
	log_debug("xpmem:adios function inq var\n");
    return NULL;
}

ADIOS_VARINFO* 
adios_read_xpmem_inq_var_byid (const ADIOS_FILE * fp, int varid)
{
	ADIOS_VARINFO *varinfo;

	adios_errno = 0;

	varinfo = bp_inq_var_byid(fp, varid);

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
adios_read_xpmem_reset_dimension_order (const ADIOS_FILE *adiosfile, int is_fortran)
{
    //log_debug( "debug: adios_read_xpmem_reset_dimension_order\n");
    adios_error(err_invalid_read_method, "adios_read_xpmem_reset_dimension_order is not implemented.");
}

