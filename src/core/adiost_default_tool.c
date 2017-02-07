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
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#define ADIOST_EXTERN 
#define DEBUG_PRINT //printf("In %s!\n", __func__); fflush(stdout);
#define DEBUG_PRINT_FD //printf("file_descriptor: %d!\n", file_descriptor); fflush(stdout);
#define ONE_BILLION 1000000000
#define ONE_BILLIONF 1000000000.0

#ifndef _GNU_SOURCE 
#define _GNU_SOURCE 
#endif // _GNU_SOURCE 
#include <errno.h> /* to get program name */

/* to get the program name from glibc */
extern char *program_invocation_name;
extern char *program_invocation_short_name;


/* Enumeration of timer indices. */
enum adiost_timer_index {
    adiost_open_timer = 0,
    adiost_close_timer,
    adiost_open_to_close_timer,
    adiost_read_timer,
    adiost_write_timer,
    adiost_advance_step_timer,
    adiost_group_size_timer,
    adiost_transform_timer,
    adiost_last_timer_unused
};

/* Enumeration of counter indices */
enum adiost_counter_index {
    adiost_data_bytes = 0,
    adiost_total_bytes,
    adiost_last_counter_unused
};

/* Array of timers for all timed events. This is a static
 * array, limited to the number of events. */
static uint64_t adiost_timers_accumulated[adiost_last_timer_unused] = {0ULL};
static uint64_t adiost_timers_count[adiost_last_timer_unused] = {0ULL};
static struct timespec adiost_timers_start_time[adiost_last_timer_unused];
static uint64_t adiost_counters_count[adiost_last_counter_unused] = {0ULL};
static uint64_t adiost_counters_accumulated[adiost_last_counter_unused] = {0ULL};

void __timer_start(enum adiost_timer_index index) {
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(adiost_timers_start_time[index]));
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

void __timer_stop(enum adiost_timer_index index) {
    struct timespec end_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    uint64_t diff = timespec_subtract(&end_time, &(adiost_timers_start_time[index]));
    adiost_timers_accumulated[index] = adiost_timers_accumulated[index] + diff;
    adiost_timers_count[index] = adiost_timers_count[index] + 1;
}

ADIOST_EXTERN void my_adios_open_begin ( int64_t file_descriptor, 
    const char * group_name, const char * file_name, const char * mode) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    printf("group_name: %s!\n", group_name); fflush(stdout);
    printf("file_name: %s!\n", file_name); fflush(stdout);
    printf("mode: %s!\n", mode); fflush(stdout);
    __timer_start(adiost_open_to_close_timer);
    __timer_start(adiost_open_timer);
}

ADIOST_EXTERN void my_adios_open_end(int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_stop(adiost_open_timer);
}

ADIOST_EXTERN void my_adios_close_begin(int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_start(adiost_close_timer);
}

ADIOST_EXTERN void my_adios_close_end(int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_stop(adiost_close_timer);
    __timer_stop(adiost_open_to_close_timer);
}

ADIOST_EXTERN void my_adios_write_begin( int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_start(adiost_write_timer);
}

ADIOST_EXTERN void my_adios_write_end(int64_t file_descriptor) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_stop(adiost_write_timer);
} 

ADIOST_EXTERN void my_adios_read_begin( int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_start(adiost_read_timer);
}

ADIOST_EXTERN void my_adios_read_end(int64_t file_descriptor) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_stop(adiost_read_timer);
} 

ADIOST_EXTERN void my_adios_advance_step_begin( int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_start(adiost_advance_step_timer);
}

ADIOST_EXTERN void my_adios_advance_step_end(int64_t file_descriptor) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_stop(adiost_advance_step_timer);
} 

ADIOST_EXTERN void my_adios_group_size_begin(int64_t file_descriptor) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_start(adiost_group_size_timer);
} 

ADIOST_EXTERN void my_adios_group_size_end(int64_t file_descriptor, 
    uint64_t data_size, uint64_t total_size) {
    DEBUG_PRINT
    printf("data size: %d!\n", data_size); fflush(stdout);
    adiost_counters_accumulated[adiost_data_bytes] = adiost_counters_accumulated[adiost_data_bytes] + data_size;
    adiost_counters_count[adiost_data_bytes] = adiost_counters_count[adiost_data_bytes] + 1;
    printf("total size: %d!\n", total_size); fflush(stdout);
    adiost_counters_accumulated[adiost_total_bytes] = adiost_counters_accumulated[adiost_total_bytes] + total_size;
    adiost_counters_count[adiost_total_bytes] = adiost_counters_count[adiost_total_bytes] + 1;
    __timer_stop(adiost_group_size_timer);
}

ADIOST_EXTERN void my_adios_transform_begin( int64_t file_descriptor) {
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_start(adiost_transform_timer);
}

ADIOST_EXTERN void my_adios_transform_end(int64_t file_descriptor) { 
    DEBUG_PRINT
    DEBUG_PRINT_FD
    __timer_stop(adiost_transform_timer);
} 

ADIOST_EXTERN void my_adios_finalize(void) {
    DEBUG_PRINT
    if (adiost_timers_count[adiost_open_timer] > 0ULL) {
        printf("%s: adios_open, %u calls, %3.9f seconds\n", 
            program_invocation_short_name,
            adiost_timers_count[adiost_open_timer],
            ((double)adiost_timers_accumulated[adiost_open_timer])/ONE_BILLIONF);
    }
    if (adiost_timers_count[adiost_close_timer] > 0ULL) {
        printf("%s: adios_close, %u calls, %3.9f seconds\n", 
            program_invocation_short_name,
            adiost_timers_count[adiost_close_timer],
            ((double)adiost_timers_accumulated[adiost_close_timer])/ONE_BILLIONF);
    }
    if (adiost_timers_count[adiost_open_to_close_timer] > 0ULL) {
        printf("%s: adios_open_to_close, %u calls, %3.9f seconds\n", 
            program_invocation_short_name,
            adiost_timers_count[adiost_open_to_close_timer],
            ((double)adiost_timers_accumulated[adiost_open_to_close_timer])/ONE_BILLIONF);
    }
    if (adiost_timers_count[adiost_group_size_timer] > 0ULL) {
        printf("%s: adios_group_size, %u calls, %3.9f seconds\n", 
            program_invocation_short_name,
            adiost_timers_count[adiost_group_size_timer],
            ((double)adiost_timers_accumulated[adiost_group_size_timer])/ONE_BILLIONF);
    }
    if (adiost_timers_count[adiost_advance_step_timer] > 0ULL) {
        printf("%s: adios_advance_step, %u calls, %3.9f seconds\n", 
            program_invocation_short_name,
            adiost_timers_count[adiost_advance_step_timer],
            ((double)adiost_timers_accumulated[adiost_advance_step_timer])/ONE_BILLIONF);
    }
    if (adiost_timers_count[adiost_transform_timer] > 0ULL) {
        printf("%s: adios_transform, %u calls, %3.9f seconds\n", 
            program_invocation_short_name,
            adiost_timers_count[adiost_transform_timer],
            ((double)adiost_timers_accumulated[adiost_transform_timer])/ONE_BILLIONF);
    }
    if (adiost_timers_count[adiost_read_timer] > 0ULL) {
        printf("%s: adios_read, %u calls, %3.9f seconds\n", 
            program_invocation_short_name,
            adiost_timers_count[adiost_read_timer],
            ((double)adiost_timers_accumulated[adiost_read_timer])/ONE_BILLIONF);
    }
    if (adiost_timers_count[adiost_write_timer] > 0ULL) {
        printf("%s: adios_write, %u calls, %3.9f seconds\n", 
            program_invocation_short_name,
            adiost_timers_count[adiost_write_timer],
            ((double)adiost_timers_accumulated[adiost_write_timer])/ONE_BILLIONF);
    }
    if (adiost_counters_count[adiost_data_bytes] > 0ULL) {
        printf("%s: adios data written, %u calls, %u bytes\n", 
            program_invocation_short_name,
            adiost_counters_count[adiost_data_bytes],
            adiost_counters_accumulated[adiost_data_bytes]);
    }
    if (adiost_counters_count[adiost_total_bytes] > 0ULL) {
        printf("%s: adios total written, %u calls, %u bytes\n", 
            program_invocation_short_name,
            adiost_counters_count[adiost_total_bytes],
            adiost_counters_accumulated[adiost_total_bytes]);
    }
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

ADIOST_EXTERN void __default_adiost_initialize (adiost_function_lookup_t adiost_fn_lookup,
    const char *runtime_version, unsigned int adiost_version) {

    adiost_set_callback_t adiost_fn_set_callback = 
        (adiost_set_callback_t)adiost_fn_lookup("adiost_set_callback");

    fprintf(stderr,"Registering ADIOS tool events..."); fflush(stderr);
    CHECK(adiost_event_open_begin,         my_adios_open_begin,          "adios_open_begin");
    CHECK(adiost_event_open_end,           my_adios_open_end,            "adios_open_end");
    CHECK(adiost_event_close_begin,        my_adios_close_begin,         "adios_close_begin");
    CHECK(adiost_event_close_end,          my_adios_close_end,           "adios_close_end");
    CHECK(adiost_event_write_begin,        my_adios_write_begin,         "adios_write_begin");
    CHECK(adiost_event_write_end,          my_adios_write_end,           "adios_write_end");
    CHECK(adiost_event_read_begin,         my_adios_read_begin,          "adios_read_begin");
    CHECK(adiost_event_read_end,           my_adios_read_end,            "adios_read_end");
    CHECK(adiost_event_advance_step_begin, my_adios_advance_step_begin,  "adios_advance_step_begin");
    CHECK(adiost_event_advance_step_end,   my_adios_advance_step_end,    "adios_advance_step_end");
    CHECK(adiost_event_group_size_begin,   my_adios_group_size_begin,    "adios_group_size_begin");
    CHECK(adiost_event_group_size_end,     my_adios_group_size_end,      "adios_group_size_end");
    CHECK(adiost_event_transform_begin,    my_adios_transform_begin,     "adios_transform_begin");
    CHECK(adiost_event_transform_end,      my_adios_transform_end,       "adios_transform_end");
    CHECK(adiost_event_library_shutdown,   my_adios_finalize,            "adios_finalize");
}

adiost_initialize_t adiost_tool() { return __default_adiost_initialize; }

