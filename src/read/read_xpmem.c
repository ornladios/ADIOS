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

#include <xpmem.h>

#include "public/adios_xpmem.h"

// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define BYTE_ALIGN 8
#define MINIFOOTER_SIZE 28


#define BUFREAD8(b,var)  var = (uint8_t) *(b->buff + b->offset); \
                         b->offset += 1;

#define BUFREAD16(b,var) var = *(uint16_t *) (b->buff + b->offset); \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_16(var); \
                         b->offset += 2;

#define BUFREAD32(b,var) var = *(uint32_t *) (b->buff + b->offset); \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_32(var); \
                         b->offset += 4;

#define BUFREAD64(b,var) var = *(uint64_t *) (b->buff + b->offset); \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_64(var); \
                         b->offset += 8;


typedef struct _xpmem_read_data
{
	xpmem_segid_t data_segid;
	xpmem_segid_t index_segid;
	xpmem_apid_t data_apid;
	xpmem_apid_t index_apid;

	shared_data *pg;
	shared_data *index;
	
}xpmem_read_data;

typedef struct _xpmem_read_file
{
	uint32_t timestep;
	struct adios_bp_buffer_struct_v1 *b;
	struct bp_minifooter mfooter; 

	char *index;
	uint64_t isize;
	char *data;
	uint64_t dsize;
	
	xpmem_read_data *fp;	
}xpmem_read_file;

xpmem_read_data* fp = NULL;


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
	do
	{
		read_segid(&fp->data_segid, "xpmem.data");
		buffer = attach_segid(fp->data_segid, share_size, &fp->data_apid);
	}while(buffer == NULL);
	
	do
	{
		read_segid(&fp->index_segid, "xpmem.index");
		index = attach_segid(fp->index_segid, index_share_size, &fp->index_apid);
	}while(index == NULL);

	
	fprintf(stderr, "opened the segments %p %p \n", buffer, index);
	fp->pg = (shared_data*)buffer;
	fp->index = (shared_data*)index;

	fprintf(stderr, "version %u size = %ll readcount = %u offset = %u reader = %u buffer=%p\n", fp->pg->version,
	        fp->pg->size,
	        fp->pg->readcount,
	        fp->pg->offset,
	        fp->pg->reader,
	        fp->pg->buffer);

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
	f->fp = fp;
	f->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));

	struct bp_minifooter * mh = &f->mfooter;
    struct adios_bp_buffer_struct_v1 *b = f->b;
    uint64_t attrs_end;


	af->fh = (uint64_t)fp;
	af->current_step = 0;

	adios_buffer_struct_init(f->b);

	log_info("spinning on version\n");
	
	while(f->fp->pg->version == 0)
	    adios_nanosleep(0, 100000000);

	while(f->fp->index->version == 0)
	    adios_nanosleep(0, 100000000);

	//now the buffer has some data


	
	//just check to make sure that readcount is 0
	log_debug("xpmem readcount = %d, %d\n",
	          f->fp->pg->readcount, f->fp->index->readcount);

	//now fp->index contains the index
	//index size is fp-index->size
	//copy the buffer out
	
	f->index = (char*)malloc(f->fp->index->size);
	f->data = (char*)malloc(f->fp->pg->size);

	f->isize = f->fp->index->size;
	f->dsize = f->fp->pg->size;
	
	memcpy(f->index, f->fp->index->buffer, f->isize);
	memcpy(f->data, f->fp->pg->buffer, f->dsize);
	
	f->b->file_size = f->dsize;
	f->mfooter.file_size = f->dsize;

	if(!f->b->buff)
	{
		bp_alloc_aligned(f->b, MINIFOOTER_SIZE);
		memset(f->b->buff, 0, MINIFOOTER_SIZE);
		f->b->offset = 0;
	}

	memcpy(f->b->buff, &f->data[f->dsize - MINIFOOTER_SIZE],
	       MINIFOOTER_SIZE);
	b->offset = MINIFOOTER_SIZE - 4;
	adios_parse_version(b, &mh->version);
	mh->change_endianness = b->change_endianness;

	log_info("mh->version = %d\n", mh->version);

	b->offset = 0;
	
	BUFREAD64(b, b->pg_index_offset);
	mh->pgs_index_offset = b->pg_index_offset;
    if (b->pg_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. PG index offset (%lld) > file size (%lld)\n",
                b->pg_index_offset, b->file_size);
        return NULL;
    }
    BUFREAD64(b, b->vars_index_offset)
    mh->vars_index_offset = b->vars_index_offset;
    // validity check  
    if (b->vars_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Variable index offset (%lld) > file size (%lld)\n",
                b->vars_index_offset, b->file_size);
        return NULL;
    }
    if (b->vars_index_offset < b->pg_index_offset) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Variable index offset (%lld) < PG index offset (%lld)\n",
                b->vars_index_offset, b->pg_index_offset);
        return NULL;
    }


    BUFREAD64(b, b->attrs_index_offset)
    mh->attrs_index_offset = b->attrs_index_offset;
    // validity check  
    if (b->attrs_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Attribute index offset (%lld) > file size (%lld)\n",
                b->attrs_index_offset, b->file_size);
        return NULL;
    }
    
    if (b->attrs_index_offset < b->vars_index_offset) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Attribute index offset (%lld) < Variable index offset (%lld)\n",
                b->attrs_index_offset, b->vars_index_offset);
        return NULL;
    }

    b->end_of_pgs = b->pg_index_offset;
    b->pg_size = b->vars_index_offset - b->pg_index_offset;
    b->vars_size = b->attrs_index_offset - b->vars_index_offset;
    attrs_end = b->file_size - MINIFOOTER_SIZE;
    b->attrs_size = attrs_end - b->attrs_index_offset;

    log_debug("offsets are %llu, %llu %llu\n",
              b->pg_index_offset, b->vars_index_offset,
              b->attrs_index_offset);
    
    uint64_t footer_size = mh->file_size - mh->pgs_index_offset;
    bp_realloc_aligned (b, footer_size);

    memcpy(b->buff, &f->data[mh->pgs_index_offset], footer_size);
    b->offset = 0;
    
	   
	
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
	ADIOS_FILE *af = (ADIOS_FILE *)malloc(sizeof(ADIOS_FILE));
	if(af == NULL)
	{
		adios_error(err_no_memory, "adios_read_xpmem_open");
		return NULL;
	}

	xpmem_read_file *f = (xpmem_read_file*)malloc(sizeof(xpmem_read_file));
	f->fp = fp;

	af->fh = (uint64_t)fp;
	af->current_step = 0;
	af->last_step = 0;

	while(fp->pg->version != 0)
	    adios_nanosleep(0, 100000000);

	while(fp->index->version != 0)
	    adios_nanosleep(0, 100000000);

	//now the buffer has some data

	//just check to make sure that readcount is 0
	log_debug("xpmem readcount = %d, %d\n",
	          fp->pg->readcount, fp->index->readcount);
	

    return af;
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
	if(fp->index->readcount < 1)
		fp->index->readcount = 1;

	
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
	if(fp->index->readcount < 1)
		fp->index->readcount = 1;

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
	log_warn("not supported");
}

int adios_read_xpmem_check_reads(const ADIOS_FILE* fp, ADIOS_VARCHUNK** chunk) { log_debug( "xpmem:adios function check reads\n"); return 0; }

int adios_read_xpmem_perform_reads(const ADIOS_FILE *adiosfile, int blocking)
{
    return 0;
}

int
adios_read_xpmem_inq_var_blockinfo(const ADIOS_FILE* fp,
				      ADIOS_VARINFO* varinfo)
{ /*log_debug( "xpmem:adios function inq var block info\n");*/ return 0; }

int
adios_read_xpmem_inq_var_stat(const ADIOS_FILE* fp,
				 ADIOS_VARINFO* varinfo,
				 int per_step_stat,
				 int per_block_stat)
{ /*log_debug( "xpmem:adios function inq var stat\n");*/ return 0; }


int 
adios_read_xpmem_schedule_read_byid(const ADIOS_FILE *adiosfile,
				       const ADIOS_SELECTION *sel,
				       int varid,
				       int from_steps,
				       int nsteps,
				       void *data)
{   
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
    return 0;
}

int 
adios_read_xpmem_get_attr (int *gp, const char *attrname,
                                 enum ADIOS_DATATYPES *type,
                                 int *size, void **data)
{
    return adios_errno;
}

int 
adios_read_xpmem_get_attr_byid (const ADIOS_FILE *adiosfile, int attrid,
				   enum ADIOS_DATATYPES *type,
				   int *size, void **data)
{
    return adios_errno;
}

ADIOS_VARINFO* 
adios_read_xpmem_inq_var(const ADIOS_FILE * adiosfile, const char* varname)
{
    return NULL;
}

ADIOS_VARINFO* 
adios_read_xpmem_inq_var_byid (const ADIOS_FILE * adiosfile, int varid)
{
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
