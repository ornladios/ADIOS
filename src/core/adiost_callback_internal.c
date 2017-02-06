/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include "adiost_callback_internal.h"
#include "adios_version.h"
#include <stdio.h>

#define adiost_get_callback_success 1
#define adiost_get_callback_failure 0

#define no_tool_present 0

#define ADIOST_API_ROUTINE static

#ifndef ADIOST_STR_MATCH
#define ADIOST_STR_MATCH(haystack, needle) (!strcasecmp(haystack, needle))
#endif

/* static/global variables for this file */
adiost_callbacks_t adiost_callbacks;
static adiost_initialize_t adiost_initialize_fn = NULL;
static adiost_interface_fn_t adiost_fn_lookup(const char *s);
int adios_tool_enabled = 0;
const char adiost_enabled_env_var[] = {"ADIOS_TOOL_ENABLED"};

adiost_state_info_t adiost_state_info[] = {
#define adiost_state_macro(state, code) { # state, state },
    FOREACH_ADIOST_STATE(adiost_state_macro)
    #undef adiost_state_macro
};

/* Weak function definitions */

__attribute__ (( weak )) adiost_initialize_t adiost_initialize() {
    return NULL;
}

/* Pre-initialization. */

void adiost_pre_init(void) {
    // prevent calling this function more than once
    static int adiost_pre_initialized = 0;
    if (adiost_pre_initialized) return;
    adiost_pre_initialized = 1;

    // check whether we are enabled or disabled
    const char *adiost_env_var = getenv(adiost_enabled_env_var);
    enum adios_tool_setting_e tool_setting = adiost_error;

    // convert the input to our internal enumeration
    if (adiost_env_var == NULL || strcmp(adiost_env_var, "")) {
        tool_setting = adiost_unset;
    } else if (strcmp(adiost_env_var, "disabled") == 0) {
        tool_setting = adiost_disabled;
    } else if (strcmp(adiost_env_var, "enabled") == 0) {
        tool_setting = adiost_enabled;
    }

    // validate the input
    printf("%s: %s = %d\n", __func__, adiost_enabled_env_var, tool_setting);
    switch(tool_setting) {
        case adiost_disabled:
            break;
        case adiost_unset:
        case adiost_enabled:
            adiost_initialize_fn = adiost_tool();
            // if initialization is successful, we are enabled
            if (adiost_initialize_fn) {
                adios_tool_enabled = 1;
            }
            break;
        case adiost_error:
            fprintf(stderr, "Warning: %s has invalid value '%s'.\n", 
                adiost_enabled_env_var, adiost_env_var);
            fprintf(stderr, "Legal values are NULL, 'enabled', 'disabled'.\n");
            break;
    }
    printf("%s: adiost_enabled = %d\n", __func__, adios_tool_enabled);
}

void adiost_post_init(void) {
    // prevent calling this function more than once
    static int adiost_post_initialized = 0;
    if (adiost_post_initialized) return;
    adiost_post_initialized = 1;

    // initialize the tool if specified
    if (adios_tool_enabled) {
        adiost_initialize_fn(adiost_fn_lookup, ADIOS_VERSION, ADIOST_VERSION);
    }
}

void adiost_finalize(void) {
    if (adios_tool_enabled) {
        // call the registered shutdown callback function
        if (adiost_callbacks.adiost_callback(adiost_event_library_shutdown)) {
            adiost_callbacks.adiost_callback(adiost_event_library_shutdown)();
        }
    }
    // prevent any further callback API attempts
    adios_tool_enabled = 0;
}

/* Define all callbacks using macro expansion */

/*****************************************************************************
 * interface operations
 ****************************************************************************/

/*****************************************************************************
 * state
 ****************************************************************************/

ADIOST_API_ROUTINE int adiost_enumerate_state(int current_state, int *next_state,
                                          const char **next_state_name)
{
    const static int len = sizeof(adiost_state_info) / sizeof(adiost_state_info_t);
    int i = 0;

    for (i = 0; i < len - 1; i++) {
        if (adiost_state_info[i].state_id == current_state) {
            *next_state = adiost_state_info[i+1].state_id;
            *next_state_name = adiost_state_info[i+1].state_name;
            return 1;
        }
    }

    return 0;
}


/*****************************************************************************
 * callbacks
 ****************************************************************************/

ADIOST_API_ROUTINE int adiost_set_callback(adiost_event_t evid, adiost_callback_t cb)
{
    switch (evid) {

#define adiost_event_macro(event_name, callback_type, event_id)                  \
    case event_name:                                                           \
        adiost_callbacks.adiost_callback(event_name) = (callback_type) cb;     \
        return adiost_set_result_registration_success;

    FOREACH_ADIOST_EVENT(adiost_event_macro)

#undef adiost_event_macro

    default: return adiost_set_result_registration_error;
    }
}

ADIOST_API_ROUTINE int adiost_get_callback(adiost_event_t evid, adiost_callback_t *cb)
{
    switch (evid) {

#define adiost_event_macro(event_name, callback_type, event_id)                  \
    case event_name:                                                           \
        if (1) {                    \
            adiost_callback_t mycb =                                             \
                (adiost_callback_t) adiost_callbacks.adiost_callback(event_name);    \
            if (mycb) {                                                        \
                *cb = mycb;                                                    \
                return adiost_get_callback_success;                              \
            }                                                                  \
        }                                                                      \
        return adiost_get_callback_failure;

    FOREACH_ADIOST_EVENT(adiost_event_macro)

#undef adiost_event_macro

    default: return adiost_get_callback_failure;
    }
}

ADIOST_API_ROUTINE adiost_state_t adiost_get_state(void)
{
    adiost_state_t thread_state = adiost_state_undefined;
    return thread_state;
}

/*****************************************************************************
 * API inquiry for tool
 ****************************************************************************/

static adiost_interface_fn_t adiost_fn_lookup(const char *s)
{

#define adiost_interface_fn(fn) \
    if (strcmp(s, #fn) == 0) return (adiost_interface_fn_t) fn;

    FOREACH_ADIOST_INQUIRY_FN(adiost_interface_fn)

    //FOREACH_ADIOST_PLACEHOLDER_FN(adiost_interface_fn)

    return (adiost_interface_fn_t) 0;
}

