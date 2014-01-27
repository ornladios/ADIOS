/*
  read_flexpath.c       
  Goal: to create evpath io connection layer in conjunction with 
  write/adios_flexpath.c
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

// evpath libraries
#include <ffs.h>
#include <atl.h>
//#include <gen_thread.h>
#include <evpath.h>

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
#include "public/adios_error.h"

// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif

/********** Core ADIOS Read functions. **********/

/*
 * Gathers basic MPI information; sets up reader CM.
 */
int
adios_read_icee_init_method (MPI_Comm comm, PairStruct* params)
{     
    adios_error (err_operation_not_supported, "No support yet\n");
    return 0;
}

ADIOS_FILE*
adios_read_icee_open_file(const char * fname, MPI_Comm comm)
{
    adios_error (err_operation_not_supported, "No support yet\n");
    return NULL;
}

ADIOS_FILE*
adios_read_icee_open(const char * fname,
			 MPI_Comm comm,
			 enum ADIOS_LOCKMODE lock_mode,
			 float timeout_sec)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return NULL;
}

int 
adios_read_icee_finalize_method ()
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

void 
adios_read_icee_release_step(ADIOS_FILE *adiosfile) 
{
    adios_error (err_operation_not_supported, "Not impleted yet");
}

int 
adios_read_icee_advance_step(ADIOS_FILE *adiosfile, int last, float timeout_sec) 
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

int 
adios_read_icee_close(ADIOS_FILE * fp)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

ADIOS_FILE *
adios_read_icee_fopen(const char *fname, MPI_Comm comm) 
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return NULL;
}

int 
adios_read_icee_is_var_timed(const ADIOS_FILE* fp, int varid) 
{  
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0; 
}

void 
adios_read_icee_get_groupinfo(const ADIOS_FILE *fp, int *ngroups, char ***group_namelist, uint32_t **nvars_per_group, uint32_t **nattrs_per_group) 
{
    adios_error (err_operation_not_supported, "Not impleted yet");
}

int 
adios_read_icee_check_reads(const ADIOS_FILE* fp, ADIOS_VARCHUNK** chunk) 
{ 
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0; 
}

int adios_read_icee_perform_reads(const ADIOS_FILE *adiosfile, int blocking)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

int
adios_read_icee_inq_var_blockinfo(const ADIOS_FILE* fp,
				      ADIOS_VARINFO* varinfo)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0; 
}

int
adios_read_icee_inq_var_stat(const ADIOS_FILE* fp,
				 ADIOS_VARINFO* varinfo,
				 int per_step_stat,
				 int per_block_stat)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0; 
}


int 
adios_read_icee_schedule_read_byid(const ADIOS_FILE *adiosfile,
				       const ADIOS_SELECTION *sel,
				       int varid,
				       int from_steps,
				       int nsteps,
				       void *data)
{   
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

int 
adios_read_icee_schedule_read(const ADIOS_FILE *adiosfile,
			const ADIOS_SELECTION * sel,
			const char * varname,
			int from_steps,
			int nsteps,
			void * data)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

int 
adios_read_icee_get_attr (int *gp, const char *attrname,
                                 enum ADIOS_DATATYPES *type,
                                 int *size, void **data)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

int 
adios_read_icee_get_attr_byid (const ADIOS_FILE *adiosfile, int attrid,
				   enum ADIOS_DATATYPES *type,
				   int *size, void **data)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

ADIOS_VARINFO* 
adios_read_icee_inq_var(const ADIOS_FILE * adiosfile, const char* varname)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return NULL;
}

ADIOS_VARINFO* 
adios_read_icee_inq_var_byid (const ADIOS_FILE * adiosfile, int varid)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return NULL;
}

void 
adios_read_icee_free_varinfo (ADIOS_VARINFO *adiosvar)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return;
}


ADIOS_TRANSINFO* 
adios_read_icee_inq_var_transinfo(const ADIOS_FILE *gp, 
                                      const ADIOS_VARINFO *vi)
{    
    adios_error (err_operation_not_supported, "Not impleted yet");
    return NULL;
}


int 
adios_read_icee_inq_var_trans_blockinfo(const ADIOS_FILE *gp, 
                                            const ADIOS_VARINFO *vi, 
                                            ADIOS_TRANSINFO *ti)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return 0;
}

void 
adios_read_icee_reset_dimension_order (const ADIOS_FILE *adiosfile, 
                                           int is_fortran)
{
    adios_error (err_operation_not_supported, "Not impleted yet");
    return;
}
