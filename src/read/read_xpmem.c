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

// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define BYTE_ALIGN 8
#define MINIFOOTER_SIZE 28


#define BUFREAD8(b,var)  do{var = (uint8_t) *(b->buff + b->offset);     \
		b->offset += 1;}while(0);

#define BUFREAD16(b,var) do{var = *(uint16_t *) (b->buff + b->offset);  \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_16(var); \
                         b->offset += 2;}while(0);

#define BUFREAD32(b,var) do{var = *(uint32_t *) (b->buff + b->offset);  \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_32(var); \
                         b->offset += 4;}while(0);

#define BUFREAD64(b,var) do{var = *(uint64_t *) (b->buff + b->offset);  \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_64(var); \
                         b->offset += 8;}while(0);


typedef struct _xpmem_read_data
{
	xpmem_segid_t data_segid;
	xpmem_apid_t data_apid;
	shared_data *pg;
	char *data;
	uint64_t dsize;
	shared_data debug;
}xpmem_read_data;

typedef struct _xpmem_read_file
{
	uint32_t timestep;
	struct adios_bp_buffer_struct_v1 *b;
	struct bp_minifooter mfooter; 

	BP_FILE *fh;
	xpmem_read_data *fp;

}xpmem_read_file;

xpmem_read_data* fp = NULL;
static int show_hidden_attrs = 0; // don't show hidden attr by default

/* variety of functions taken from bp_utils to make them work with
 * shared memory style method
 */

static int print_shared(shared_data *pg)
{
		fprintf(stderr, "version %u size = %ll readcount = %u offset = %u reader = %u buffer=%p\n", pg->version,
	        pg->size,
	        pg->readcount,
	        pg->offset,
	        pg->reader,
	        pg->buffer);
}

static int xp_seek_to_step (ADIOS_FILE * fp, int tostep, int show_hidden_attrs)
{
	xpmem_read_file *xf = (xpmem_read_file*)fp->fh;
	BP_FILE *fh = xf->fh;

	int j, k, t, allstep;
    struct adios_index_var_struct_v1 * var_root = fh->vars_root;
    struct adios_index_attribute_struct_v1 * attr_root;
    uint64_t i;
    int lenpath, lenname;

    /* Streaming starts with step 0. However, time index in BP file
     * starts with 1. If 'tostep' is -1, that means we want to get all steps.
     * If not, we seek to the specified step.
     */

    //TIME FIX
    // if (tostep == -1)
    // {
    //     allstep = 1;
    // }
    // else
    // {
    //     allstep = 0;
    //     t = get_time (var_root, tostep);
    // }

    allstep = 1;
    
    /* Prepare vars */
    fp->nvars = 0;
//    var_root = fh->vars_root;

    while (var_root)
    {
        for (i = 0; i < var_root->characteristics_count; i++)
        {
            if (allstep || (!allstep && var_root->characteristics[i].time_index == t))
            {
                fp->nvars++;
                break;
            }
        }

        var_root = var_root->next;
    }

    fp->var_namelist = (char **) malloc (sizeof (char *) * fp->nvars);

    var_root = fh->vars_root;
    j = 0;
    k = 0;
    while (var_root)
    {
        for (i = 0; i < var_root->characteristics_count; i++)
        {
            if (allstep || (!allstep && var_root->characteristics[i].time_index == t))
            {
                /* From 1.6, relative and full path (starts with /) are handled separately in search */
                // Full name of variable: concatenate var_path and var_name
                lenpath = strlen(var_root->var_path);
                lenname = strlen(var_root->var_name);
                if (lenpath > 0) {
                    fp->var_namelist [j] = (char *) malloc (lenname + lenpath + 1 + 1);
                                                                    // extra / and ending \0
                    strcpy(fp->var_namelist[j], var_root->var_path);
                    if (var_root->var_path[lenpath-1] != '/') {
                        fp->var_namelist[j][lenpath] = '/';
                        lenpath++;
                    }
                    strcpy(&(fp->var_namelist[j][lenpath]), var_root->var_name);
                }
                else {
                    fp->var_namelist[j] = (char *) malloc (lenname+1); 
                    strcpy(fp->var_namelist[j], var_root->var_name);
                }
                //printf ("Seek to step: Variable %d full path is [%s]\n", j, fp->var_namelist[j]);

                j++;

                break;
            }
        }

        k++;
        var_root = var_root->next;
    }

    /* Prepare attrs */
    fp->nattrs = 0;
    attr_root = fh->attrs_root;

    while (attr_root)
    {
        if (!show_hidden_attrs && strstr (attr_root->attr_path, "__adios__"))
        {
        }
        else
        {
            for (i = 0; i < attr_root->characteristics_count; i++)
            {
                if (allstep || (!allstep && attr_root->characteristics[i].time_index == t))
                {
                    fp->nattrs++;
                    break;
                }
            }
        }

        attr_root = attr_root->next;
    }

    fp->attr_namelist = (char **) malloc (sizeof (char *) * fp->nattrs);

    attr_root = fh->attrs_root;
    j = 0;
    while (attr_root)
    {
        if (!show_hidden_attrs && strstr (attr_root->attr_path, "__adios__"))
        {
        }
        else
        {
            for (i = 0; i < attr_root->characteristics_count; i++)
            {
                if (allstep || (!allstep && attr_root->characteristics[i].time_index == t))
                {
                    // Full name of attribute: concatenate attr_path and attr_name
                    lenpath = strlen(attr_root->attr_path);
                    lenname = strlen(attr_root->attr_name);
                    if (lenpath > 0) {
                        fp->attr_namelist [j] = (char *) malloc (lenname + lenpath + 1 + 1);
                                                                    // extra / and ending \0
                        strcpy(fp->attr_namelist[j], attr_root->attr_path);
                        if (attr_root->attr_path[lenpath-1] != '/') {
                            fp->attr_namelist[j][lenpath] = '/';
                            lenpath++;
                        }
                        strcpy(&(fp->attr_namelist[j][lenpath]), attr_root->attr_name);
                    }
                    else {
                        fp->attr_namelist[j] = (char *) malloc (lenname+1); 
                        strcpy(fp->attr_namelist[j], attr_root->attr_name);
                    }
                    //printf ("Seek to step: Attribute %d full path is [%s], path=[%s], name=[%s]\n", 
                    //        j, fp->attr_namelist[j], attr_root->attr_path, attr_root->attr_name);
                    j++;

                    break;
                }
            }
        }

        attr_root = attr_root->next;
    }

    fp->current_step = tostep;

    return 0;
}


static int xp_read_open(BP_FILE *fh, xpmem_read_data *f)
{
	//we don't need MPI_Comm or file name here
	struct adios_bp_buffer_struct_v1 *b;
	struct bp_minifooter * mh = &fh->mfooter;
	uint64_t attrs_end;
	uint64_t footer_size;	

	fh->sfh = NULL;
	fh->b = malloc(sizeof(struct adios_bp_buffer_struct_v1));
	fh->mfooter.file_size = f->dsize;
	fh->mpi_fh = 0;
	fh->fname = strdup("xpmem");
	fh->comm = 0;
	b = fh->b;
	
	log_info("file size = %llu\n", fh->b->file_size);

	adios_buffer_struct_init(fh->b);
	fh->b->file_size = f->dsize;

	//now we read the minifooter

	//allocate memory if we don't have it already
	if(!fh->b->buff)
	{
		bp_alloc_aligned(fh->b, MINIFOOTER_SIZE);
		memset(fh->b->buff, 0, MINIFOOTER_SIZE);
		fh->b->offset = 0;
	}

	memcpy(b->buff, &f->data[f->dsize - MINIFOOTER_SIZE],
	       MINIFOOTER_SIZE);

	b->offset = MINIFOOTER_SIZE - 4;

	adios_parse_version(b, &mh->version);
	mh->change_endianness = b->change_endianness;

	log_info("mh->version = %d\n", mh->version);

	//reset the offset for b
	b->offset = 0;

	//now we identify where each of the offsets are

	//pg index offset
	BUFREAD64(b, b->pg_index_offset);
	mh->pgs_index_offset = b->pg_index_offset;
    if (b->pg_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. PG index offset (%lld) > file size (%lld)\n",
                b->pg_index_offset, b->file_size);
        return -1;
    }

    //vars_index_offset
    BUFREAD64(b, b->vars_index_offset)
    mh->vars_index_offset = b->vars_index_offset;
    // validity check  
    if (b->vars_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Variable index offset (%lld) > file size (%lld)\n",
                b->vars_index_offset, b->file_size);
        return -1;
    }
    if (b->vars_index_offset < b->pg_index_offset) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Variable index offset (%lld) < PG index offset (%lld)\n",
                b->vars_index_offset, b->pg_index_offset);
        return -1;
    }


    //attrs_index_offset
    BUFREAD64(b, b->attrs_index_offset)
    mh->attrs_index_offset = b->attrs_index_offset;
    // validity check  
    if (b->attrs_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Attribute index offset (%lld) > file size (%lld)\n",
                b->attrs_index_offset, b->file_size);
        return -1;
    }
    
    if (b->attrs_index_offset < b->vars_index_offset) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Attribute index offset (%lld) < Variable index offset (%lld)\n",
                b->attrs_index_offset, b->vars_index_offset);
        return -1;
    }

    //set the values in the buffer
    b->end_of_pgs = b->pg_index_offset;
    b->pg_size = b->vars_index_offset - b->pg_index_offset;
    b->vars_size = b->attrs_index_offset - b->vars_index_offset;
    
    attrs_end = b->file_size - MINIFOOTER_SIZE;
    b->attrs_size = attrs_end - b->attrs_index_offset;

    log_debug("offsets are %llu, %llu %llu\n",
              b->pg_index_offset, b->vars_index_offset,
              b->attrs_index_offset);

    //now figure out the rest of the footer

    //footer starts at pgs_index_offset and goes to end of file
    footer_size =  mh->file_size - mh->pgs_index_offset;
    bp_realloc_aligned (b, footer_size);

    //copy the footer into b starting from pgs_index_offset
    //to end of buffer
    
    memcpy(b->buff, &f->data[mh->pgs_index_offset], footer_size);

    //reset offset so we can read from buffer
    b->offset = 0;	

    bp_parse_pgs(fh);
    bp_parse_vars(fh);
    bp_parse_attrs(fh);
    
	
	return 0;
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
    f->fh->gvar_h = 0;
    f->fh->pgs_root = 0;
    f->fh->vars_root = 0;
    f->fh->attrs_root = 0;
    f->fh->vars_table = 0;
    
	//set the read data to the file
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
adios_read_xpmem_schedule_read_byid(const ADIOS_FILE *adiosfile,
				       const ADIOS_SELECTION *sel,
				       int varid,
				       int from_steps,
				       int nsteps,
				       void *data)
{
	log_debug("xpmem:adios function schedule read byid\n");
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
adios_read_xpmem_inq_var_byid (const ADIOS_FILE * adiosfile, int varid)
{
	log_debug("xpmem:adios function inq var byid\n");
	return NULL;
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
    adios_error(err_operation_not_supported, "Xpmem does not yet support transforms: var_transinfo.\n");
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

