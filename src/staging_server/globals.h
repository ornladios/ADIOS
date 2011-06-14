#ifndef __STAGING_GLOBALS_H__
#define __STAGING_GLOBALS_H__

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#if HAVE_PTHREAD
#   include <pthread.h>
#endif
#include <mpi.h>

#define MAXPATHLEN 255

#include "mymutex.h"   // pthreads mutex/cond.var or empty mutex/cond.var with pthreads syntax
#include "queue.h"
//#include "events.h"
#include "ardma_logger.h" // logging functions only from libardma. 

typedef int bool;
#define false 0
#define true  1

#define MIN(a,b) (a <= b ? a : b)
#define MAX(a,b) (a >= b ? a : b)

#ifdef __MAIN__
#define EXT
#else
#define EXT extern
#endif

#define BYTES_512MB 536870912
#define BYTES_1GB   1073741824
#define BYTES_1GB   1073741824
#define BYTES_2GB   2147483648
static const uint64_t MAX_BLOCK_SIZE = BYTES_512MB; // do not write more than that at once
                                                    // unless user specifies it with --maxblock

/* User arguments */
EXT uint64_t user_max_memory_allowed;  // max buffer size to pull client data
EXT uint64_t user_max_block_size;      // pull/write this size in total (unless one client data is bigger)
EXT int verbose;                       // verbosity level for logging
EXT bool logfile_separate_ranks;       // create separate log file per rank
EXT char * logpath;                    // log file, default=stderr

// for each client 
struct globals_client_data {
    // set by transport, used by worker
    uint32_t rank;          // client rank
    uint32_t nid;           // node id of client process 
    uint64_t pg_size;       // size of data block (process group buffer)
    uint64_t idx_size;      // size of data block (index buffer)
    uint64_t idx_offset;    // offset in local rdma buffer for the pulled index 

    // set and used by worker only
    uint64_t pg_offset;     // offset in local rdma buffer for the data
    uint64_t file_offset;   // offset in output file
    int      order_idx;     // 0..gd.nc-1, the order of pulling of this client's data 

    enum { 
        PULLSTAT_INIT, 
        PULLSTAT_REQUESTED, 
        PULLSTAT_DONE 
    }       status;        // PG pulling status
};

struct global_data_struct {   
#if HAVE_PTHREAD
    pthread_t worker_thread;
#endif
    int  nc;                   // number of clients (expected to connect to this server process)
    int  nc_current;           // current number of clients connected to this server process
    int  nc_total;             // total number of clients to all server processes
    int  lrank;                // the lowest rank of clients connecting to this server process
    bool application_finished; // true: application exited so all threads 
                               // should exit after finishing remaining work

    queue_t *wtq;              // Worker to Transport event Queue
    queue_t *twq;              // Transport to Worker event Queue

    MPI_Comm mpi_comm;
    int mpi_rank;
    int mpi_size;
    struct precedence_struct * prec_writer_indexing; // ordering of writer and indexing thread for MPI calls

    uint64_t requested_rdma_buffer_size; // how much mem to use for rdma pull 
                                         // [set in worker_global_init()]
    char *rdma_buffer;         // buffer allocated in ARDMA, size told by worker_global_init
    uint64_t rdma_bufsize;     // actual size of allocated and registered memory
    uint64_t index_size;       // rdma_buffer[0..index_size-1] occupied by indices pulled from clients

    bool need_index; // should pull index buffers [set in worker_global_init()]
    bool terminate;  // true: all threads should terminate asap

    struct globals_client_data * clientdata; // array of nc elements
    int *order;               // ordering of clients for pulling order, created by worker

};

EXT struct global_data_struct gd;


#undef EXT
#endif
