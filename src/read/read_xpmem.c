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

#include <xpmem.h>

#include "public/adios_xpmem.h"

// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif


typedef struct _xpmem_read_file
{
	xpmem_segid_t data_segid;
	xpmem_segid_t index_segid;
	xpmem_apid_t data_apid;
	xpmem_apid_t index_apid;

	shared_data *pg;
	shared_data *index;
	
}xpmem_read_file, xpmem_read_data;

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
	read_segid(&fp->data_segid, "xpmem.data");
	read_segid(&fp->index_segid, "xpmem.index");

	buffer = attach_segid(fp->data_segid, share_size, &fp->data_apid);
	index = attach_segid(fp->index_apid, index_share_size, &fp->index_apid);

	fp->pg = (shared_data*)buffer;
	fp->index = (shared_data*)index;
	
	

    return 0;
}

ADIOS_FILE*
adios_read_xpmem_open_file(const char * fname, MPI_Comm comm)
{
    adios_error (err_operation_not_supported,
                 "Xpmem staging method does not support file mode for reading. "
                 "Use adios_read_open() to open a staged dataset.\n");
    return NULL;
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
	ADIOS_FILE *fp = (ADIOS_FILE *)malloc(sizeof(ADIOS_FILE));
	if(fp == NULL)
	{
		adios_error(err_no_memory, "adios_read_xpmem_open");
		return NULL;
	}
	
    return fp;
}

int adios_read_xpmem_finalize_method ()
{
    return 0;
}

void adios_read_xpmem_release_step(ADIOS_FILE *adiosfile) {
}

int 
adios_read_xpmem_advance_step(ADIOS_FILE *adiosfile, int last, float timeout_sec) 
{
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
