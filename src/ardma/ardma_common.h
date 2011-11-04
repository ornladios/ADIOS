/*
   Put here all definitions that all implementations of 
   ardma_server.h and ardma_client.h can use. 
*/

#ifndef ARDMA_COMMON_H
#define ARDMA_COMMON_H

#include <stdint.h> /* uint64_t */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define ARDMA_CONNECTION_FILENAME "rdma_connections"

/*
extern int ardma_verbose_level; // must be defined in main file, and must be set 
#define ardma_logger(verbose_level, ...) if (ardma_verbose_level >= verbose_level) fprintf (stderr, __VA_ARGS__); 
#define log_error(...) ardma_logger(0, __VA_ARGS__)
#define log_warn(...) ardma_logger(1, __VA_ARGS__)
#define log_info(...) ardma_logger(2, __VA_ARGS__)
#define log_debug(...) ardma_logger(3, __VA_ARGS__)
*/

#define  TIMEOUT_IN_MS 1000    /* ms */

static inline void * alloc_memory_aligned(uint64_t size)
{
    void * p;
    int fail_alloc = posix_memalign (&p, sysconf (_SC_PAGE_SIZE), size);
    if (fail_alloc)
        return 0;
    return p;
}

#endif
