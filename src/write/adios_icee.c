
/*
    adios_flexpath.c
    uses evpath for io in conjunction with read/read_flexpath.c
*/


#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <pthread.h>

// xml parser
#include <mxml.h>

// add by Kimmy 10/15/2012
#include <sys/types.h>
#include <sys/stat.h>
// end of change

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include "public/adios_mpi.h"
#include "public/adios_error.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"
#include "core/adios_logger.h"

#if HAVE_FLEXPATH==1

// // evpath libraries
#include <evpath.h>
#include <cod.h>
#include <sys/queue.h>

// Initializes flexpath write local data structures
extern void 
adios_icee_init(const PairStruct *params, struct adios_method_struct *method) 
{
    adios_error (err_operation_not_supported, "No support yet\n");
    return;
}

extern int 
adios_icee_open(struct adios_file_struct *fd, 
		    struct adios_method_struct *method, 
		    MPI_Comm comm) 
{    
    adios_error (err_operation_not_supported, "No support yet\n");
    return 0;	
}

//  writes data to multiqueue
extern void
adios_icee_write(
    struct adios_file_struct *fd, 
    struct adios_var_struct *f, 
    void *data, 
    struct adios_method_struct *method) 
{
}

extern void 
adios_icee_close(struct adios_file_struct *fd, struct adios_method_struct *method) 
{
}

// wait until all open files have finished sending data to shutdown
extern void 
adios_icee_finalize(int mype, struct adios_method_struct *method) 
{
}

// provides unknown functionality
extern enum ADIOS_FLAG 
adios_icee_should_buffer (struct adios_file_struct * fd,struct adios_method_struct * method) 
{
    return adios_flag_no;
}

// provides unknown functionality
extern void 
adios_icee_end_iteration(struct adios_method_struct *method) 
{
}

// provides unknown functionality
extern void 
adios_icee_start_calculation(struct adios_method_struct *method) 
{
}

// provides unknown functionality
extern void 
adios_icee_stop_calculation(struct adios_method_struct *method) 
{
}

// provides unknown functionality
extern void 
adios_icee_get_write_buffer(struct adios_file_struct *fd, 
				struct adios_var_struct *v, 
				uint64_t *size, 
				void **buffer, 
				struct adios_method_struct *method) 
{
}

// should not be called from write, reason for inclusion here unknown
void 
adios_icee_read(struct adios_file_struct *fd, 
		    struct adios_var_struct *f, 
		    void *buffer, 
		    uint64_t buffer_size, 
		    struct adios_method_struct *method) 
{
}

#else // print empty version of all functions (if HAVE_FLEXPATH == 0)

void 
adios_icee_read(struct adios_file_struct *fd, 
		    struct adios_var_struct *f, 
		    void *buffer, 
		    struct adios_method_struct *method) 
{
}

extern void 
adios_icee_get_write_buffer(struct adios_file_struct *fd, 
				struct adios_var_struct *f, 
				unsigned long long *size, 
				void **buffer, 
				struct adios_method_struct *method) 
{
}

extern void 
adios_icee_stop_calculation(struct adios_method_struct *method) 
{
}

extern void 
adios_icee_start_calculation(struct adios_method_struct *method) 
{
}

extern void 
adios_icee_end_iteration(struct adios_method_struct *method) 
{
}

#endif
