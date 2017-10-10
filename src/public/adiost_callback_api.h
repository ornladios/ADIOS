/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef __ADIOST_CALLBACK_API_H__
#define __ADIOST_CALLBACK_API_H__

/****************************************************************************
 * System include files
 ****************************************************************************/

#include <stdint.h>
#ifndef ADIOST_EXPORT 
#define ADIOST_EXPORT 
#endif
#ifndef ADIOST_WEAK
#define ADIOST_WEAK
#endif

#include "adios_types.h"
#include "adios_read.h"

/****************************************************************************
 * iteration macros - to be expanded by other macros in multiple places
 ****************************************************************************/

#define FOREACH_ADIOST_INQUIRY_FN(macro)  \
    macro (adiost_set_callback)           \
    macro (adiost_get_callback)

/* For each event, specify the callback type and the enumeration value. */
#define FOREACH_ADIOST_EVENT(macro) \
macro(adiost_event_thread,                  adiost_thread_callback_t,      1) \
macro(adiost_event_open,                    adiost_open_callback_t,        3) \
macro(adiost_event_close,                   adiost_file_callback_t,        5) \
macro(adiost_event_write,                   adiost_write_callback_t,       10) \
macro(adiost_event_read,                    adiost_file_callback_t,        12) \
macro(adiost_event_advance_step,            adiost_advance_step_callback_t,        14) \
macro(adiost_event_group_size,              adiost_group_size_callback_t,  51) \
macro(adiost_event_transform,               adiost_file_callback_t,        52) \
macro(adiost_event_fp_send_finalize_msg,    adiost_file_callback_t,       102) \
macro(adiost_event_fp_send_read_msg,        adiost_file_callback_t,       105) \
macro(adiost_event_fp_add_var_to_read_msg,  adiost_file_callback_t,       106) \
macro(adiost_event_fp_copy_buffer,          adiost_file_callback_t,       205) \
macro(adiost_event_read_init_method,        adiost_read_init_method_callback_t, 300) \
macro(adiost_event_read_finalize_method,    adiost_read_finalize_method_callback_t, 301) \
macro(adiost_event_read_open,               adiost_read_open_callback_t, 302) \
macro(adiost_event_read_open_file,          adiost_read_open_file_callback_t, 303) \
macro(adiost_event_release_step,            adiost_file_callback_t, 304) \
macro(adiost_event_inq_var,                 adiost_inq_var_callback_t, 305) \
macro(adiost_event_inq_var_byid,            adiost_inq_var_byid_callback_t, 306) \
macro(adiost_event_free_varinfo,            adiost_free_varinfo_callback_t, 307) \
macro(adiost_event_inq_var_stat,            adiost_inq_var_stat_callback_t, 308) \
macro(adiost_event_inq_var_blockinfo,       adiost_inq_var_blockinfo_callback_t, 309) \
macro(adiost_event_selection_boundingbox,   adiost_selection_boundingbox_callback_t, 310) \
macro(adiost_event_selection_points,        adiost_selection_points_callback_t, 311) \
macro(adiost_event_selection_writeblock,    adiost_selection_writeblock_callback_t, 312) \
macro(adiost_event_selection_auto,          adiost_selection_auto_callback_t, 313) \
macro(adiost_event_selection_delete,        adiost_selection_delete_callback_t, 314) \
macro(adiost_event_schedule_read,           adiost_schedule_read_callback_t, 315) \
macro(adiost_event_schedule_read_byid,      adiost_schedule_read_byid_callback_t, 316) \
macro(adiost_event_perform_reads,           adiost_perform_reads_callback_t, 317) \
macro(adiost_event_check_reads,             adiost_check_reads_callback_t, 318) \
macro(adiost_event_free_chunk,              adiost_free_chunk_callback_t, 319) \
macro(adiost_event_get_attr,                adiost_get_attr_callback_t, 320) \
macro(adiost_event_get_attr_byid,           adiost_get_attr_byid_callback_t, 321) \
macro(adiost_event_type_to_string,          adiost_type_to_string_callback_t, 322) \
macro(adiost_event_type_size,               adiost_type_size_callback_t, 323) \
macro(adiost_event_get_grouplist,           adiost_get_grouplist_callback_t, 324) \
macro(adiost_event_group_view,              adiost_group_view_callback_t, 325) \
macro(adiost_event_stat_cov,                adiost_stat_cov_callback_t, 326) \
macro(adiost_event_inq_mesh_byid,           adiost_inq_mesh_byid_callback_t, 328) \
macro(adiost_event_free_meshinfo,           adiost_free_meshinfo_callback_t, 329) \
macro(adiost_event_inq_var_meshinfo,        adiost_inq_var_meshinfo_callback_t, 330) \
macro(adiost_event_library_shutdown,        adiost_callback_t,            999) \

#endif // #ifdef __ADIOST_CALLBACK_API_H__

typedef enum {
#define adiost_event_macro(event, callback, eventid) event = eventid,
    FOREACH_ADIOST_EVENT(adiost_event_macro)
    #undef adiost_event_macro
} adiost_event_t;

/*---------------------
 * set callback results
 *---------------------*/
typedef enum {
    adiost_set_result_registration_success            = 0,
    adiost_set_result_registration_error              = 1
} adiost_set_result_t;

typedef enum {
    adiost_event_enter,
    adiost_event_exit,
    adiost_event
} adiost_event_type_t;

/****************************************************************************
 * Callback signature types
 ****************************************************************************/

/* initialization */
typedef void (*adiost_interface_fn_t)(void);

typedef adiost_interface_fn_t (*adiost_function_lookup_t)(
    const char *                      /* entry point to look up       */
);

/* Events: adios_thread */
typedef void (*adiost_thread_callback_t)(
    adiost_event_type_t type,
    ADIOS_FILE * file_descriptor,
    const char * thread_name
);

/* Events: adios_open */
typedef void (*adiost_open_callback_t)(
    adiost_event_type_t type,
    int64_t file_descriptor,
    const char * group_name,
    const char * file_name,
    const char * mode
);

/* Events: adios_close, adios_read_begin...adios_write_end */
typedef void (*adiost_file_callback_t)(
    adiost_event_type_t type,
    int64_t file_descriptor);

/* Events: adios_write */
typedef void (*adiost_write_callback_t)(
    adiost_event_type_t type,
    int64_t file_descriptor, 
    const char * name, 
    enum ADIOS_DATATYPES data_type, 
    const int ndims, 
    const char * dims, 
    const void * value);

/* Events: adios_group_size */
typedef void (*adiost_group_size_callback_t)(
    adiost_event_type_t type,
    int64_t file_descriptor,
    uint64_t data_size,
    uint64_t total_size
);

/* Events: adiost_event_library_shutdown */
typedef void (*adiost_callback_t)(void);

/* ----------- Events for Pooky/Pookie/Pukie! ------------- */

/* Events: adios_read_init_method */
typedef void (*adiost_read_init_method_callback_t)(
    adiost_event_type_t type,
    enum ADIOS_READ_METHOD method, 
    MPI_Comm comm, 
    const char * parameters
);

/* Events: adios_read_finalize_method */
typedef void (*adiost_read_finalize_method_callback_t)(
    adiost_event_type_t type,
    enum ADIOS_READ_METHOD method 
);

/* Events: adios_read_open */
typedef void (*adiost_read_open_callback_t)(
    adiost_event_type_t type,
    enum ADIOS_READ_METHOD method, 
    MPI_Comm comm, 
	enum ADIOS_LOCKMODE lock_mode,
    float timeout_sec,
	ADIOS_FILE * file_descriptor
);

/* Events: adios_read_open_file */
typedef void (*adiost_read_open_file_callback_t)(
    adiost_event_type_t type,
    const char * fname,
    enum ADIOS_READ_METHOD method, 
    MPI_Comm comm,
	ADIOS_FILE * file_descriptor
);

/* Events: adios_advance_step */
typedef void (*adiost_advance_step_callback_t)(
    adiost_event_type_t type,
    ADIOS_FILE *fp,
    int last,
    float timeout_sec
);

/* Events: adios_read_inq_var */
typedef void (*adiost_inq_var_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
    const char * varname,
	ADIOS_VARINFO * varinfo
);

/* Events: adios_read_inq_var_byid */
typedef void (*adiost_inq_var_byid_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
    int varid,
	ADIOS_VARINFO * varinfo
);

/* Events: adios_read_free_varinfo */
typedef void (*adiost_free_varinfo_callback_t)(
    adiost_event_type_t type,
	ADIOS_VARINFO * varinfo
);

/* Events: adios_read_inq_var_stat */
typedef void (*adiost_inq_var_stat_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	ADIOS_VARINFO * varinfo,
	int per_prep_stat,
	int per_block_stat
);

/* Events: adios_read_inq_var_blockinfo */
typedef void (*adiost_inq_var_blockinfo_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	ADIOS_VARINFO * varinfo
);

/* Events: adios_read_selection_boundingbox */
typedef void (*adiost_selection_boundingbox_callback_t)(
    adiost_event_type_t type,
    uint64_t ndim,
    const uint64_t *start,
    const uint64_t *count,
	ADIOS_SELECTION * selection
);

/* Events: */
typedef void (*adiost_selection_points_callback_t)(
    adiost_event_type_t type,
    uint64_t ndim,
    uint64_t npoints,
    const uint64_t *points,
	ADIOS_SELECTION * container,
	int free_points_on_delete,
	ADIOS_SELECTION * selection
);

/* Events: */
typedef void (*adiost_selection_writeblock_callback_t)(
    adiost_event_type_t type,
    int index,
	ADIOS_SELECTION * selection
);

/* Events: */
typedef void (*adiost_selection_auto_callback_t)(
    adiost_event_type_t type,
    char * hints,
	ADIOS_SELECTION * selection
);

/* Events: */
typedef void (*adiost_selection_delete_callback_t)(
    adiost_event_type_t type,
	ADIOS_SELECTION * selection
);

/* Events: */
typedef void (*adiost_schedule_read_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	const ADIOS_SELECTION * selection,
	const char * varname,
	int from_steps,
	int nsteps,
	const char * param,
	void * data
);

/* Events: */
typedef void (*adiost_schedule_read_byid_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	const ADIOS_SELECTION * selection,
	int varid,
	int from_steps,
	int nsteps,
	const char * param,
	void * data
);

/* Events: */
typedef void (*adiost_perform_reads_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	int blocking
);

/* Events: */
typedef void (*adiost_check_reads_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	ADIOS_VARCHUNK **chunk
);

/* Events: */
typedef void (*adiost_free_chunk_callback_t)(
    adiost_event_type_t type,
	ADIOS_VARCHUNK *chunk
);

/* Events: */
typedef void (*adiost_get_attr_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	const char * attrname,
	enum ADIOS_DATATYPES * datatypes,
    int * size,
	void **data
);

/* Events: */
typedef void (*adiost_get_attr_byid_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	int attrid,
	enum ADIOS_DATATYPES * datatypes,
    int * size,
	void **data
);

/* Events: */
typedef void (*adiost_type_to_string_callback_t)(
    adiost_event_type_t type,
    const char * name
);

/* Events: */
typedef void (*adiost_type_size_callback_t)(
    adiost_event_type_t type,
    void * dadta,
	int size
);

/* Events: */
typedef void (*adiost_get_grouplist_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	char ***group_namelist
);

/* Events: */
typedef void (*adiost_group_view_callback_t)(
    adiost_event_type_t type,
    ADIOS_FILE *fp,
	int groupid
);

/* Events: */
typedef void (*adiost_stat_cov_callback_t)(
    adiost_event_type_t type,
    ADIOS_VARINFO * vix,
	ADIOS_VARINFO * viy,
	char * characteristic,
	uint32_t time_start,
	uint32_t time_end,
	uint32_t lag,
	double correlation
);

/* Events: */
typedef void (*adiost_inq_mesh_byid_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	int meshid,
	ADIOS_MESH * mesh
);

/* Events: */
typedef void (*adiost_free_meshinfo_callback_t)(
    adiost_event_type_t type,
	ADIOS_MESH * mesh
);

/* Events: */
typedef void (*adiost_inq_var_meshinfo_callback_t)(
    adiost_event_type_t type,
    const ADIOS_FILE *fp,
	ADIOS_VARINFO * varinfo
);

/****************************************************************************
 * ADIOST API
 ****************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif

#define ADIOST_API_FNTYPE(fn) fn##_t

#define ADIOST_API_FUNCTION(return_type, fn, args)  \
    typedef return_type (*ADIOST_API_FNTYPE(fn)) args

/****************************************************************************
 * Initialization functions
 ****************************************************************************/

ADIOST_API_FUNCTION(void, adiost_initialize, (
    adiost_function_lookup_t adiost_fn_lookup,
    const char *runtime_version,
    unsigned int adiost_version
));

/* initialization interface - to be defined by tool */
ADIOST_EXPORT adiost_initialize_t ADIOST_WEAK adiost_tool(void);

/* Registering a callback function */
ADIOST_API_FUNCTION(int, adiost_set_callback, (
    adiost_event_t event,
    adiost_callback_t callback
));

/* Getting a callback function */
ADIOST_API_FUNCTION(int, adiost_get_callback, (
    adiost_event_t event,
    adiost_callback_t *callback
));

#ifdef  __cplusplus
}
#endif

