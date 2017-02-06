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

/****************************************************************************
 * iteration macros - to be expanded by other macros in multiple places
 ****************************************************************************/

#define FOREACH_ADIOST_INQUIRY_FN(macro)  \
    macro (adiost_enumerate_state)        \
    macro (adiost_set_callback)           \
    macro (adiost_get_callback)           \
    macro (adiost_get_state)

/* For each state, specify the state type and the enumeration value. */
#define FOREACH_ADIOST_STATE(macro)                                                               \
    macro (adiost_state_first, 0x71)          /* initial enumeration state */                     \
    macro (adiost_state_idle, 0x10)            /* waiting for work */                             \
    macro (adiost_state_overhead, 0x20)        /* overhead excluding wait states */               \
    macro (adiost_state_undefined, 0x70)       /* undefined thread state */

/* For each event, specify the callback type and the enumeration value. */
#define FOREACH_ADIOST_EVENT(macro) \
macro(adiost_event_open_begin,         adiost_open_callback_t,        1) \
macro(adiost_event_open_end,           adiost_file_callback_t,        2) \
macro(adiost_event_close_begin,        adiost_file_callback_t,        3) \
macro(adiost_event_close_end,          adiost_file_callback_t,        4) \
macro(adiost_event_write_begin,        adiost_file_callback_t,        10) \
macro(adiost_event_write_end,          adiost_file_callback_t,        11) \
macro(adiost_event_read_begin,         adiost_file_callback_t,        12) \
macro(adiost_event_read_end,           adiost_file_callback_t,        13) \
macro(adiost_event_advance_step_begin, adiost_file_callback_t,        14) \
macro(adiost_event_advance_step_end,   adiost_file_callback_t,        15) \
macro(adiost_event_group_size,         adiost_group_size_callback_t,  50) \
macro(adiost_event_library_shutdown,   adiost_callback_t,             99) \

#endif // #ifdef __ADIOST_CALLBACK_API_H__

typedef enum {
#define adiost_state_macro(state, code) state = code,
    FOREACH_ADIOST_STATE(adiost_state_macro)
    #undef adiost_state_macro
} adiost_state_t;

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

/****************************************************************************
 * Callback signature types
 ****************************************************************************/

/* initialization */
typedef void (*adiost_interface_fn_t)(void);

typedef adiost_interface_fn_t (*adiost_function_lookup_t)(
    const char *                      /* entry point to look up       */
);

/* Events: adios_open */
typedef void (*adiost_open_callback_t)(
    int64_t file_descriptor,
    const char * group_name,
    const char * file_name,
    const char * mode
);

/* Events: adios_close, adios_read_begin...adios_write_end */
typedef void (*adiost_file_callback_t)(int64_t file_descriptor);

/* Events: adios_group_size */
typedef void (*adiost_group_size_callback_t)(
    int64_t file_descriptor,
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
adiost_initialize_t adiost_tool(void);

/* Initialization modes - required? */
typedef enum opt_init_mode_e {
    adiost_init_mode_never  = 0,
    adiost_init_mode_false  = 1,
    adiost_init_mode_true   = 2,
    adiost_init_mode_always = 3
} adiost_init_mode_t;

/* Error codes for when registering callbacks. */
typedef enum adiost_set_callback_rc_e {  /* non-standard */
    adiost_set_callback_error      = 0,
    adiost_has_event_no_callback   = 1,
    adiost_no_event_no_callback    = 2,
    adiost_has_event_may_callback  = 3,
    adiost_has_event_must_callback = 4,
} adiost_set_callback_rc_t;

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

