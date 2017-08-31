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

#include "config.h"
// we want our implementation of the tool to be weak, and exported.
#define ADIOST_EXPORT __attribute__((visibility("default")))
#define ADIOST_WEAK __attribute__ (( weak )) 

#include "public/adiost_callback_api.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include "adios_clock.h"
#define ADIOST_EXTERN 
#define DEBUG_PRINT printf("In %s!\n", __func__); fflush(stdout);
#define DEBUG_PRINT_FD printf("file_descriptor: %" PRId64 "!\n", file_descriptor); fflush(stdout);
#define ONE_BILLION 1000000000
#define ONE_BILLIONF 1000000000.0

#ifndef _GNU_SOURCE 
#define _GNU_SOURCE 
#endif // _GNU_SOURCE 
#include <errno.h> /* to get program name */

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

/* to get the program name from glibc */
extern char *program_invocation_name;
extern char *program_invocation_short_name;

/* Enumeration of timer indices. */
typedef enum adiost_timer_index {
    adiost_thread_timer = 0,
    adiost_open_timer,
    adiost_close_timer,
    adiost_open_to_close_timer,
    adiost_read_timer,
    adiost_write_timer,
    adiost_advance_step_timer,
    adiost_group_size_timer,
    adiost_transform_timer,
    adiost_fp_send_read_msg_timer,
    adiost_fp_send_finalize_msg_timer,
    adiost_fp_add_var_to_read_msg_timer,
    adiost_fp_copy_buffer_timer,
    adiost_last_timer_unused
} adiost_timer_index_t;

/* Enumeration of counter indices */
typedef enum adiost_counter_index {
    adiost_data_bytes_counter = 0,
    adiost_total_bytes_counter,
    adiost_last_counter_unused
} adiost_counter_index_t;

/* Array of timers for all timed events. This is a static
 * array, limited to the number of events. */
static uint64_t adiost_timers_accumulated[adiost_last_timer_unused] = {0ULL};
static uint64_t adiost_timers_count[adiost_last_timer_unused] = {0ULL};
static struct timespec adiost_timers_start_time[adiost_last_timer_unused];
static uint64_t adiost_counters_count[adiost_last_counter_unused] = {0ULL};
static uint64_t adiost_counters_accumulated[adiost_last_counter_unused] = {0ULL};

void __timer_start(adiost_timer_index_t index) {
    adios_clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(adiost_timers_start_time[index]));
}

/* Subtract the ‘struct timeval’ values X and Y,
   storing the result in RESULT.
   Return the total difference in nanoseconds. */

uint64_t timespec_subtract (struct timespec *end, struct timespec *start)
{
    struct timespec result;
    /* Perform the carry for the later subtraction by updating start. */
    if (end->tv_nsec < start->tv_nsec) {
        uint64_t nsec = (start->tv_nsec - end->tv_nsec) / ONE_BILLION + 1;
        start->tv_nsec -= ONE_BILLION * nsec;
        start->tv_sec += nsec;
    }
    if (end->tv_nsec - start->tv_nsec > ONE_BILLION) {
        uint64_t nsec = (end->tv_nsec - start->tv_nsec) / ONE_BILLION;
        start->tv_nsec += ONE_BILLION * nsec;
        start->tv_sec -= nsec;
    }

    /* Compute the time remaining to wait.
        tv_nsec is certainly positive. */
    result.tv_sec = end->tv_sec - start->tv_sec;
    result.tv_nsec = end->tv_nsec - start->tv_nsec;
    return (result.tv_sec * ONE_BILLION) + result.tv_nsec;
}

void __timer_stop(adiost_timer_index_t index) {
    struct timespec end_time;
    adios_clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    uint64_t diff = timespec_subtract(&end_time, &(adiost_timers_start_time[index]));
    adiost_timers_accumulated[index] = adiost_timers_accumulated[index] + diff;
    adiost_timers_count[index] = adiost_timers_count[index] + 1;
}

ADIOST_EXTERN void my_thread(int64_t file_descriptor, char * name, adiost_event_type_t type) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_thread_timer);
    } else {
        __timer_stop(adiost_thread_timer);
    }
}

ADIOST_EXTERN void my_open ( int64_t file_descriptor, adiost_event_type_t type,
    const char * group_name, const char * file_name, const char * mode) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        debug_print("group_name: %s!\n", group_name);
        debug_print("file_name: %s!\n", file_name);
        debug_print("mode: %s!\n", mode);
        __timer_start(adiost_open_to_close_timer);
        __timer_start(adiost_open_timer);
    } else {
        __timer_stop(adiost_open_timer);
    }
}

ADIOST_EXTERN void my_close(int64_t file_descriptor, adiost_event_type_t type) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_close_timer);
    } else {
        __timer_stop(adiost_close_timer);
        __timer_stop(adiost_open_to_close_timer);
    }
}

ADIOST_EXTERN void my_write( int64_t file_descriptor, adiost_event_type_t type, const char * name, enum ADIOS_DATATYPES data_type, const int ndims, const char * dims, const void * value) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_write_timer);
    } else {
        __timer_stop(adiost_write_timer);
    }
} 

ADIOST_EXTERN void my_read( int64_t file_descriptor, adiost_event_type_t type) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_read_timer);
    } else {
        __timer_stop(adiost_read_timer);
    }
} 

ADIOST_EXTERN void my_advance_step( int64_t file_descriptor,
    adiost_event_type_t type) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_advance_step_timer);
    } else {
        __timer_stop(adiost_advance_step_timer);
    }
} 

ADIOST_EXTERN void my_group_size(int64_t file_descriptor, 
    adiost_event_type_t type, uint64_t data_size, uint64_t total_size) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_group_size_timer);
    } else {
        debug_print("data size: %" PRIu64 "!\n", data_size); fflush(stdout);
        adiost_counters_accumulated[adiost_data_bytes_counter] = 
            adiost_counters_accumulated[adiost_data_bytes_counter] + data_size;
        adiost_counters_count[adiost_data_bytes_counter] = 
            adiost_counters_count[adiost_data_bytes_counter] + 1;
        debug_print("total size: %" PRIu64 "!\n", total_size); fflush(stdout);
        adiost_counters_accumulated[adiost_total_bytes_counter] = 
            adiost_counters_accumulated[adiost_total_bytes_counter] + total_size;
        adiost_counters_count[adiost_total_bytes_counter] = 
            adiost_counters_count[adiost_total_bytes_counter] + 1;
        __timer_stop(adiost_group_size_timer);
    }
} 

ADIOST_EXTERN void my_transform( int64_t file_descriptor,
        adiost_event_type_t type) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_transform_timer);
    } else {
        __timer_stop(adiost_transform_timer);
    }
} 

ADIOST_EXTERN void my_fp_send_read_msg(int64_t file_descriptor,
        adiost_event_type_t type) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_fp_send_read_msg_timer);
    } else {
        __timer_stop(adiost_fp_send_read_msg_timer);
    }
} 

ADIOST_EXTERN void my_fp_send_finalize_msg(int64_t file_descriptor,
        adiost_event_type_t type) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_fp_send_finalize_msg_timer);
    } else {
        __timer_stop(adiost_fp_send_finalize_msg_timer);
    }
} 

ADIOST_EXTERN void my_fp_add_var_to_read_msg(int64_t file_descriptor,
        adiost_event_type_t type) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_fp_add_var_to_read_msg_timer);
    } else {
        __timer_stop(adiost_fp_add_var_to_read_msg_timer);
    }
} 

ADIOST_EXTERN void my_fp_copy_buffer(int64_t file_descriptor,
        adiost_event_type_t type) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    if (type == adiost_event_enter) {
        __timer_start(adiost_fp_copy_buffer_timer);
    } else {
        __timer_stop(adiost_fp_copy_buffer_timer);
    }
} 

/* This function is for printing a timer */
void print_timer(const char *_name, adiost_timer_index_t index) {
    if (adiost_timers_count[index] > 0ULL) {
        debug_print("%s: %s, %" PRIu64 " calls, %3.9f seconds\n",
            program_invocation_short_name, _name,
            adiost_timers_count[index],
            ((double)adiost_timers_accumulated[index])/ONE_BILLIONF);
    }
}

/* This function is for printing a counter */
void print_counter(const char *_name, adiost_counter_index_t index) {
    if (adiost_counters_count[index] > 0ULL) {
        debug_print("%s: %s, %" PRIu64 " calls, %" PRIu64 " bytes\n",
            program_invocation_short_name, _name,
            adiost_counters_count[index],
            adiost_counters_accumulated[index]);
    }
}

ADIOST_EXTERN void my_finalize(void) {
    DEBUG_PRINT
    print_timer("adios_thread", adiost_thread_timer);
    print_timer("adios_open", adiost_open_timer);
    print_timer("adios_close", adiost_close_timer);
    print_timer("adios_open_to_close", adiost_open_to_close_timer);
    print_timer("adios_group_size", adiost_group_size_timer);
    print_timer("adios_advance_step", adiost_advance_step_timer);
    print_timer("adios_transform", adiost_transform_timer);
    print_timer("adios_read", adiost_read_timer);
    print_timer("adios_write", adiost_write_timer);
    print_timer("flexpath send open msg", adiost_fp_send_read_msg_timer);
    print_timer("flexpath send finalize msg", adiost_fp_send_finalize_msg_timer);
    print_timer("flexpath send var msg", adiost_fp_add_var_to_read_msg_timer);
    print_timer("flexpath copy buffer", adiost_fp_copy_buffer_timer);
    print_counter("adios data written", adiost_data_bytes_counter);
    print_counter("adios total written", adiost_total_bytes_counter);
}

// This function is for checking that the function registration worked.
#define CHECK(EVENT,FUNCTION,NAME) \
    debug_print("Registering ADIOST callback %s...",NAME); \
    if (adiost_fn_set_callback(EVENT, (adiost_callback_t)(FUNCTION)) != \
                    adiost_set_result_registration_success) { \
        debug_print("\n\tFailed to register ADIOST callback %s!\n",NAME); \
    } else { \
        debug_print("success.\n"); \
    } \

ADIOST_EXTERN void default_adiost_initialize (adiost_function_lookup_t adiost_fn_lookup,
    const char *runtime_version, unsigned int adiost_version) {

    adiost_set_callback_t adiost_fn_set_callback = 
        (adiost_set_callback_t)adiost_fn_lookup("adiost_set_callback");

    if (!getenv("ADIOST")) return;
    debug_print("Registering ADIOS tool events...\n");
    CHECK(adiost_event_thread,       my_thread,        "adios_thread");
    CHECK(adiost_event_open,         my_open,          "adios_open");
    CHECK(adiost_event_close,        my_close,         "adios_close");
    CHECK(adiost_event_write,        my_write,         "adios_write");
    CHECK(adiost_event_read,         my_read,          "adios_read");
    CHECK(adiost_event_advance_step, my_advance_step,  "adios_advance_step");
    CHECK(adiost_event_group_size,   my_group_size,    "adios_group_size");
    CHECK(adiost_event_transform,    my_transform,     "adios_transform");
    CHECK(adiost_event_fp_send_read_msg, 
        my_fp_send_read_msg, "adios_fp_send_read_msg");
    CHECK(adiost_event_fp_send_finalize_msg, 
        my_fp_send_finalize_msg, "adios_fp_send_finalize_msg");
    CHECK(adiost_event_fp_add_var_to_read_msg, 
        my_fp_add_var_to_read_msg, "adiost_fp_add_var_to_read_msg");
    CHECK(adiost_event_fp_copy_buffer, 
        my_fp_copy_buffer, "adios_fp_copy_buffer");
    CHECK(adiost_event_library_shutdown, my_finalize, "adios_finalize");
}

/* Weak function definitions, declarations */

extern adiost_initialize_t adiost_tool(void) { return default_adiost_initialize; }


