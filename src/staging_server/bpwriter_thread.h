#ifndef __BPWRITER_THREAD_H__
#define __BPWRITER_THREAD_H__

#include <stdint.h> // uint64_t
#include "globals.h"
#include "mymutex.h"
#include "queue.h"

/* bpwriter thread is started from bpworker thread */

enum WRITER_REQUEST_TYPE {
    /* Worker to Writer thread queue elements (in wtw queue) */
    WRITER_REQUEST_OPEN,        // open a file
    WRITER_REQUEST_WRITE,       // write a block of data 
    WRITER_REQUEST_WRITEIDX,    // write metadata block(s) (and free block(s))
    WRITER_REQUEST_CLOSE,       // close the file
    WRITER_REQUEST_EXIT,        // exit thread
};

struct writer_request_open {
    char * path;            // path from simulation client
     // Note: path is the original path from the client app
     //       this function will create path.dir/path.<mpi_id>
};

struct writer_request_write {
    uint64_t addr;  /* memory address of data */
    uint64_t size;
    uint64_t offset; /* offset in output file */
    // after completing the request, set *pso=weo if pso!=NULL
    uint64_t *pso;
    uint64_t weo;
};

struct writer_request_writeidx {
    uint64_t laddr;  /* memory address of data */
    uint64_t lsize;
    uint64_t loffset; /* offset in output file */
    // rank 0 writes the global index too
    uint64_t gaddr;  /* memory address of data */
    uint64_t gsize;
};

struct writer_request {
    enum WRITER_REQUEST_TYPE type;
    union {
        struct writer_request_open     open;
        struct writer_request_write    write;
        struct writer_request_writeidx writeidx;
        // close and exit has no data struct
    } request;
};

extern queue_t *woq; // Worker to Writer queue


/*
enum WRITER_RESPONSE_TYPE {
    WRITER_RESPONSE_OPEN,
    WRITER_RESPONSE_WRITE,
    WRITER_RESPONSE_CLOSE
};

extern queue_t *wrq; // Writer to Worker queue

struct writer_info_struct {
    pthread_mutex_t mutex;
    int             n_outstanding_requests; // = 0: writer is done with all work
};
extern struct writer_info_struct writer_info;
*/

/*
static inline void WRITE_REQUESTS_INCREASE() 
{
    pthread_mutex_lock(&writer_info->mutex); 
    writer_info.n_outstanding_requests++; 
    pthread_mutex_unlock(&writer_info->mutex);
}

static inline void WRITE_REQUESTS_DECREASE() 
{
    pthread_mutex_lock(&writer_info->mutex); 
    writer_info.n_outstanding_requests--; 
    pthread_mutex_unlock(&writer_info->mutex);
}

static inline bool WRITE_REQUESTS_ISEMPTY() 
{
    bool b;
    pthread_mutex_lock(&writer_info->mutex); 
    b = writer_info.n_outstanding_requests--; 
    pthread_mutex_unlock(&writer_info->mutex);
    return b;
}
*/


/* Writer thread */
#ifdef HAVE_PTHREAD
pthread_t writer_thread;
void * bpwriter_thread_main (void *arg);
#else
void bpwriter_init (void);
int  bpwriter_main (void); // 0: no events, 1: one event processed
void bpwriter_finalize (void);
#endif


#endif
