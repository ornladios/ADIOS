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
macro(adiost_event_advance_step,            adiost_file_callback_t,        14) \
macro(adiost_event_group_size,              adiost_group_size_callback_t,  51) \
macro(adiost_event_transform,               adiost_file_callback_t,        52) \
macro(adiost_event_fp_send_finalize_msg,    adiost_file_callback_t,       102) \
macro(adiost_event_fp_send_read_msg,        adiost_file_callback_t,       105) \
macro(adiost_event_fp_add_var_to_read_msg,  adiost_file_callback_t,       106) \
macro(adiost_event_fp_copy_buffer,          adiost_file_callback_t,       205) \
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
    int64_t file_descriptor,
    adiost_event_type_t type,
    const char * thread_name
);

/* Events: adios_open */
typedef void (*adiost_open_callback_t)(
    int64_t file_descriptor,
    adiost_event_type_t type,
    const char * group_name,
    const char * file_name,
    const char * mode
);

/* Events: adios_close, adios_read_begin...adios_write_end */
typedef void (*adiost_file_callback_t)(
    int64_t file_descriptor, 
    adiost_event_type_t type);

/* Events: adios_write */
typedef void (*adiost_write_callback_t)(
    int64_t file_descriptor, 
    adiost_event_type_t type,
    const char * name, 
    enum ADIOS_DATATYPES data_type, 
    const int ndims, 
    const char * dims, 
    const void * value);

/* Events: adios_group_size */
typedef void (*adiost_group_size_callback_t)(
    int64_t file_descriptor,
    adiost_event_type_t type,
    uint64_t data_size,
    uint64_t total_size
);

/* Events: adiost_event_library_shutdown */
typedef void (*adiost_callback_t)(void);

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

