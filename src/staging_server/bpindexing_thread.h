#ifndef __BPINDEXING_THREAD_H__
#define __BPINDEXING_THREAD_H__

#include <stdint.h> // uint64_t
#include "globals.h"

/* bpindexing thread is started from bpworker thread and
   it blocks on a condition variable.

   Requests: do indexing, no parameter is given (use global variables)
             exit 
   Response is the pointer/size for local index plus
      pointer/size for global index on rank 0
*/

enum INDEXING_REQUEST_TYPE {
    INDEXING_REQUEST_INDEX,     // do indexing (use gd.rdma_buffer, gd.index_size as args)
    INDEXING_REQUEST_EXIT      // just exit
};

struct indexing_data_struct {
    bool indexing_completed; // index thread sets to true when indexing is done
    enum INDEXING_REQUEST_TYPE indexing_request;

    uint64_t laddr;  /* memory address of local index */
    uint64_t lsize;  /* size of local index */
    // rank 0 only:
    uint64_t gaddr;  /* memory address of global index*/
    uint64_t gsize;  /* size of local index */
};

// mutexed variables
extern struct indexing_data_struct indexing_data;


/* indexing thread */
#ifdef HAVE_PTHREAD
    pthread_t       indexing_thread;
    pthread_mutex_t indexing_mutex;
    pthread_cond_t  indexing_cv;
    void * bpindexing_thread_main (void *arg);
    // tell indexing thread to start working
#   define bpindexing_request_doindex() { \
        pthread_mutex_lock(&indexing_mutex); \
        indexing_data.indexing_request = INDEXING_REQUEST_INDEX; \
        indexing_data.indexing_completed = false; \
        pthread_cond_signal(&indexing_cv); \
        pthread_mutex_unlock(&indexing_mutex); \
    }
    // tell indexing thread to exit
#   define bpindexing_request_doexit() { \
        pthread_mutex_lock(&indexing_mutex); \
        indexing_data.indexing_request = INDEXING_REQUEST_EXIT; \
        pthread_cond_signal(&indexing_cv); \
        pthread_mutex_unlock(&indexing_mutex); \
    }

#else
    void bpindexing_doindex (void); 
#endif


#endif
