/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

// Tell the header file that this we want to declare the adiost_tool()
// definition as "weak".
#define ADIOST_INTERNAL
#include "adiost_callback_internal.h"

#include "adios_version.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#define adiost_get_callback_success 1
#define adiost_get_callback_failure 0
#define no_tool_present 0
#define ADIOST_API_ROUTINE static

/* The do { ... } while (0) idiom ensures that the code acts like a statement
 * (function call). The unconditional use of the code ensures that the compiler
 * always checks that your debug code is valid — but the optimizer will remove
 * the code when DEBUG is 0. */
#define DEBUG 0
#define debug_print(...) do { \
    if (DEBUG) { \
        fprintf(stderr, __VA_ARGS__); \
        fflush(stderr);\
    } \
} while (0)

/* static/global variables for this file */
adiost_callbacks_t adiost_callbacks;
static adiost_initialize_t adiost_initialize_fn = NULL;
static adiost_interface_fn_t adiost_fn_lookup(const char *s);
int adios_tool_enabled = 0;
const char adiost_enabled_env_var[] = {"ADIOS_TOOL"};

/* forward declaration of the weak (default) tool */

extern adiost_initialize_t default_adiost_tool(void);

/* function pointer to hold the tool function. */
adiost_initialize_t (*my_adiost_tool)(void) = NULL;

/* Pre-initialization. */

void adiost_pre_init(void) {
    // prevent calling this function more than once
    static int adiost_pre_initialized = 0;
    if (adiost_pre_initialized) return;
    adiost_pre_initialized = 1;

    // check whether we are enabled or disabled
    char *adiost_env_var = (char *)getenv(adiost_enabled_env_var);
    adios_tool_setting_t tool_setting = adiost_error;

    // convert the input to our internal enumeration
    if (adiost_env_var == NULL || strlen(adiost_env_var) == 0) {
        tool_setting = adiost_unset;
    } else if (strcmp(adiost_env_var, "disabled") == 0) {
        tool_setting = adiost_disabled;
    } else if (strcmp(adiost_env_var, "enabled") == 0) {
        tool_setting = adiost_enabled;
    }

	// if a tool function is defined, assign our internal pointer to it.
    // for clang, we always have a weak definition, so prevent compiler warning.
#if defined(__clang__)
	if (adiost_tool() != NULL) {
#else
	if ((adiost_tool != NULL) && (adiost_tool() != NULL)){
#endif
	    my_adiost_tool = &adiost_tool;
	} else {
	    my_adiost_tool = &default_adiost_tool;
	}

    // validate the input
    debug_print("%s: %s = %d\n", __func__, adiost_enabled_env_var, tool_setting);
    switch(tool_setting) {
        case adiost_disabled:
            break;
        case adiost_unset:
        case adiost_enabled:
            if (my_adiost_tool) {
            	adiost_initialize_fn = my_adiost_tool();
            	// if initialization is successful, we are enabled
            	if (adiost_initialize_fn) {
                	adios_tool_enabled = 1;
            	}
            }
            break;
        case adiost_error:
            fprintf(stderr, "Warning: %s has invalid value '%s'.\n", 
                adiost_enabled_env_var, adiost_env_var);
            fprintf(stderr, "Legal values are NULL, 'enabled', 'disabled'.\n");
            break;
    }
    debug_print("%s: adiost_enabled = %d\n", __func__, adios_tool_enabled);
}

/* Post-initialization */

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

/* For shutting down the tool support */

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

/*****************************************************************************
 * interface operations
 ****************************************************************************/

/*****************************************************************************
 * callbacks
 ****************************************************************************/

ADIOST_API_ROUTINE int adiost_set_callback(adiost_event_t evid, adiost_callback_t cb)
{
    switch (evid) {

/* Define all callbacks using macro expansion */

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

/* Define all callbacks using macro expansion */

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

char * adiost_build_dimension_string(struct adios_var_struct *v, int * ndims) { 
    // I hate to use a fixed length string, but...
    char tmpstr[1024] = {0};
    *ndims = 0;
    if (v->dimensions == NULL) {
        return strdup("");
    }
    char dims[256] = {0};
    char global_dims[256] = {0};
    char local_offsets[256] = {0};
    struct adios_dimension_struct * tmp = v->dimensions;
    char delimiter = '[';
    while (tmp != NULL) {
        *ndims = *ndims + 1;
        // just a regular old number? Get its value.
        if (tmp->dimension.rank > 0) {
#if defined(__clang__)
            sprintf(dims, "%s%c%llu", dims, delimiter, tmp->dimension.rank);
#else
            sprintf(dims, "%s%c%lu", dims, delimiter, tmp->dimension.rank);
#endif
        // another ADIOS variable? Get its name.
        } else if (tmp->dimension.var != NULL) {
            sprintf(dims, "%s%c%s", dims, delimiter, tmp->dimension.var->name);
        // another ADIOS attribute? Get its name.
        } else if (tmp->dimension.attr != NULL) {
            sprintf(dims, "%s%c%s", dims, delimiter, tmp->dimension.attr->name);
        }
        // just a regular old number? Get its value.
        if (tmp->global_dimension.rank > 0) {
#if defined(__clang__)
            sprintf(global_dims, "%s%c%llu", global_dims, delimiter, tmp->global_dimension.rank);
#else
            sprintf(global_dims, "%s%c%lu", global_dims, delimiter, tmp->global_dimension.rank);
#endif
        // another ADIOS variable? Get its name.
        } else if (tmp->global_dimension.var != NULL) {
            sprintf(global_dims, "%s%c%s", global_dims, delimiter, tmp->global_dimension.var->name);
        // another ADIOS attribute? Get its name.
        } else if (tmp->global_dimension.attr != NULL) {
            sprintf(global_dims, "%s%c%s", global_dims, delimiter, tmp->global_dimension.attr->name);
        }
        // just a regular old number? Get its value.
        if (tmp->local_offset.rank > 0) {
#if defined(__clang__)
            sprintf(local_offsets, "%s%c%llu", local_offsets, delimiter, tmp->local_offset.rank);
#else
            sprintf(local_offsets, "%s%c%lu", local_offsets, delimiter, tmp->local_offset.rank);
#endif
        // another ADIOS variable? Get its name.
        } else if (tmp->local_offset.var != NULL) {
            sprintf(local_offsets, "%s%c%s", local_offsets, delimiter, tmp->local_offset.var->name);
        // another ADIOS attribute? Get its name.
        } else if (tmp->local_offset.attr != NULL) {
            sprintf(local_offsets, "%s%c%s", local_offsets, delimiter, tmp->local_offset.attr->name);
        }
        // move on to the next dimension?
        delimiter = ',';
        tmp = tmp->next;
    }
    delimiter = ']';
	if (strlen(dims) > 0) {
        sprintf(dims, "%s%c", dims, delimiter);
	} else {
        sprintf(dims, "[]");
	}
	if (strlen(global_dims) > 0) {
        sprintf(global_dims, "%s%c", global_dims, delimiter);
	} else {
        sprintf(global_dims, "[]");
	}
	if (strlen(local_offsets) > 0) {
        sprintf(local_offsets, "%s%c", local_offsets, delimiter);
	} else {
        sprintf(local_offsets, "[]");
	}
    // build the whole thing
    sprintf(tmpstr, "%s;%s;%s", dims, global_dims, local_offsets);
    return strdup(tmpstr);
}

// Create a weak definition for the tool, some compiler/systems require it.
// For example, Clang on OSX requires this symbol to be defined.
// To ensure that the tool can instantiate its own definition, only create a
// weak definition when required, such as with clang.
#if defined(__clang__)
ADIOST_EXPORT ADIOST_WEAK_PRE adiost_initialize_t adiost_tool(void) ADIOST_WEAK_POST
{
    return NULL;
}
#endif
