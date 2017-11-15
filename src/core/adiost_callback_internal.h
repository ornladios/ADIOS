/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef __ADIOST_CALLBACK_INTERNAL_H__
#define __ADIOST_CALLBACK_INTERNAL_H__

/* ADIOS Event Callback API
 *
 * This pair of source files (callback_api.h/c) is used to define a callback
 * API for any tool that wants to be notified when ADIOS events happen. For
 * example, a performance tool may want to know when staged writes happen
 * and be notified how many bytes are transferred, in order to estimate 
 * effective bandwidth.
 */

#include <stddef.h>
// we want our implementation of the tool to be weak, and exported.
#define ADIOST_EXPORT __attribute__((visibility("default")))
#if defined(__clang__)
#define ADIOST_WEAK_PRE 
#define ADIOST_WEAK_POST __attribute__((weak_import))
#else
#define ADIOST_WEAK_PRE __attribute__ (( weak )) 
#define ADIOST_WEAK_POST 
#endif
#include "public/adiost_callback_api.h"
#include "adios_internals.h"

/* Definitions */

#define ADIOST_VERSION 20170202 // 2017, February 02
#define ADIOST_EXTERN extern "C"
extern int adios_tool_enabled;

typedef enum adios_tool_setting_e {
    adiost_error = 0,
    adiost_unset,
    adiost_disabled,
    adiost_enabled
} adios_tool_setting_t;

/* management of registered function callbacks */
#define adiost_callback(e) e ## _callback

typedef struct adiost_callbacks_s {
#define adiost_event_macro(event, callback, eventid) callback adiost_callback(event);

    FOREACH_ADIOST_EVENT(adiost_event_macro)

    #undef adiost_event_macro
} adiost_callbacks_t;


/* Global variables */

extern adiost_callbacks_t adiost_callbacks;

/* internal callback API functions */

void adiost_pre_init(void);
void adiost_post_init(void);
void adiost_finalize(void);
char * adiost_build_dimension_string(struct adios_var_struct *v, int * ndims);

/* These variadic macros save us from having to add 4+ lines of source code
 * each time that we want to make an event callback. */

#define ADIOST_CALLBACK(__event, ...) \
    if (adios_tool_enabled && \
        adiost_callbacks.adiost_callback(__event)) { \
          adiost_callbacks.adiost_callback(__event)( \
            adiost_event, ##__VA_ARGS__); \
    } \

#define ADIOST_CALLBACK_ENTER(__event, ...) \
    if (adios_tool_enabled && \
        adiost_callbacks.adiost_callback(__event)) { \
          adiost_callbacks.adiost_callback(__event)( \
            adiost_event_enter, ##__VA_ARGS__); \
    } \

#define ADIOST_CALLBACK_EXIT(__event, ...) \
    if (adios_tool_enabled && \
        adiost_callbacks.adiost_callback(__event)) { \
          adiost_callbacks.adiost_callback(__event)( \
            adiost_event_exit, ##__VA_ARGS__); \
    } \

#define ADIOST_CALLBACK_WRITE_ENTER(__event, __fd, __v) \
    if (adios_tool_enabled && \
        adiost_callbacks.adiost_callback(__event)) { \
          int __ndims = 0; \
          char * __tmp = adiost_build_dimension_string(__v, &__ndims); \
          adiost_callbacks.adiost_callback(__event)(adiost_event_enter, (int64_t)__fd, \
            __v->name, __v->type, __ndims, __tmp, __v->data); \
    } \

#define ADIOST_CALLBACK_WRITE_EXIT(__event, __fd, __v) \
    if (adios_tool_enabled && \
        adiost_callbacks.adiost_callback(__event)) { \
          int __ndims = 0; \
          char * __tmp = adiost_build_dimension_string(__v, &__ndims); \
          adiost_callbacks.adiost_callback(__event)(adiost_event_exit, (int64_t)__fd, \
            __v->name, __v->type, __ndims, __tmp, __v->data); \
    } \

#endif // #ifndef __ADIOST_CALLBACK_INTERNAL_H__
