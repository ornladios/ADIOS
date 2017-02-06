/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS Event Callback API - Default tool implementation
 *
 * This source file is a default implementation of the ADIOS callback API
 * for tools, or ADIOST. This is a simple tool that will add performance
 * data to the ADIOS metadata.
 */

#include "adiost_callback_api.h"
#include <stdio.h>
#define ADIOST_EXTERN 
#define DEBUG_PRINT printf("In %s!\n", __func__); fflush(stdout);
#define DEBUG_PRINT_FD printf("file_descriptor: %d!\n", file_descriptor); fflush(stdout);

ADIOST_EXTERN void my_adios_open_begin ( int64_t file_descriptor, 
    const char * group_name, const char * file_name, const char * mode) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    printf("group_name: %s!\n", group_name); fflush(stdout);
    printf("file_name: %s!\n", file_name); fflush(stdout);
    printf("mode: %s!\n", mode); fflush(stdout);
}

ADIOST_EXTERN void my_adios_open_end(int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
}

ADIOST_EXTERN void my_adios_close_begin(int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
}

ADIOST_EXTERN void my_adios_close_end(int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
}

ADIOST_EXTERN void my_adios_write_begin( int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
}

ADIOST_EXTERN void my_adios_write_end(int64_t file_descriptor) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
} 

ADIOST_EXTERN void my_adios_read_begin( int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
}

ADIOST_EXTERN void my_adios_read_end(int64_t file_descriptor) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
} 

ADIOST_EXTERN void my_adios_group_size(int64_t file_descriptor, 
    uint64_t data_size, uint64_t total_size) {
    DEBUG_PRINT
    printf("data size: %d!\n", data_size); fflush(stdout);
    printf("total size: %d!\n", total_size); fflush(stdout);
}

ADIOST_EXTERN void my_adios_finalize(void) {
    DEBUG_PRINT
}

// This macro is for checking that the function registration worked.
#define CHECK(EVENT,FUNCTION,NAME) \
    printf("Registering ADIOST callback %s...",NAME); \
    fflush(stderr); \
    if (adiost_fn_set_callback(EVENT, (adiost_callback_t)(FUNCTION)) != \
                    adiost_set_result_registration_success) { \
        printf("\n\tFailed to register ADIOST callback %s!\n",NAME); \
        fflush(stderr); \
    } else { \
        printf("success.\n"); \
    } \

ADIOST_EXTERN void __default_adiost_initialize (
    adiost_function_lookup_t adiost_fn_lookup,
    const char *runtime_version,
    unsigned int adiost_version) {
    adiost_set_callback_t adiost_fn_set_callback = adiost_fn_lookup("adiost_set_callback");
    fprintf(stderr,"Registering ADIOS tool events..."); fflush(stderr);
    CHECK(adiost_event_open_begin,       my_adios_open_begin,  "adios_open_begin");
    CHECK(adiost_event_open_end,         my_adios_open_end,    "adios_open_end");
    CHECK(adiost_event_close_begin,      my_adios_close_begin, "adios_close_begin");
    CHECK(adiost_event_close_end,        my_adios_close_end,   "adios_close_end");
    CHECK(adiost_event_write_begin,      my_adios_write_begin, "adios_write_begin");
    CHECK(adiost_event_write_end,        my_adios_write_end,   "adios_write_end");
    CHECK(adiost_event_read_begin,       my_adios_read_begin,  "adios_read_begin");
    CHECK(adiost_event_read_end,         my_adios_read_end,    "adios_read_end");
    CHECK(adiost_event_group_size,       my_adios_group_size,  "adios_group_size");
    CHECK(adiost_event_library_shutdown, my_adios_finalize,    "adios_finalize");
}

adiost_initialize_t adiost_tool() { return __default_adiost_initialize; }

