/* Worker thread for BP writing staging server.
   It pulls ADIOS PGs from client and writes into
   files as is, effectively like the POSIX method does
*/
#include <stdlib.h> // qsort
#include <string.h> // memset
#include <time.h>   // nanosleep


#include "worker_thread.h"
#include "bpwriter_thread.h"
#include "bpindexing_thread.h"
#include "globals.h"
#include "mem.h"
#include "events.h"
#include "mymutex.h"
#include "precedence.h"   // pthread cond.var to have strict order of access to MPI Collectives in
                          // the two sub-threads (writer and indexing)

/* Worker should tell
   1. size of rdma buffer used for client data
   2. if index data will be used (will be pulled by transport thread automatically)
   This is called in the initial main thread, 
   NOT from the worker thread!
 */
void worker_global_init (void)
{
    uint64_t available_memory = mem_get_available (user_max_memory_allowed);
        log_info("rank %d: available memory = %lld bytes\n", gd.mpi_rank, available_memory);

    gd.requested_rdma_buffer_size = available_memory*0.9; // Use 90% to pull data
    gd.need_index = 1;
}

static enum {
    WORKER_READY,    // ready to process a new request (may write index to previous file)
    WORKER_PULLING,  // pulling data (and also working on it)
    WORKER_WORKING   // still writing data (but not pulling more & requested all writes)
} worker_status;

#if HAVE_PTHREAD
int worker_main (void *arg);
void worker_init (void);
void * worker_thread_main (void *arg)
{
    int num_event;
    struct timespec treq = {.tv_sec=0, .tv_nsec=2000000}; // 2ms sleeptime
    struct timespec trem;
    worker_init();
    while (!gd.application_finished || worker_status != WORKER_READY) {
        num_event = worker_main (arg);
        if (num_event == 0)
            nanosleep(&treq, &trem);
        if (gd.terminate)
            break;
    }
    if (gd.terminate) {
        log_debug("rank %d: worker thread: exit due to failure.\n", gd.mpi_rank);
    } else {
        log_debug("rank %d: worker thread: application has exited and no more work to do.\n", gd.mpi_rank);
    }

    worker_finalize();
    pthread_exit(NULL);
    return NULL; // just to avoid compiler warning
}
#endif

/* Local variables */
static int np; // number of PGs to pull at once (and write out)
//static int *order; // ordering of clients for pulling order
static struct cyclic_buffer {
    char *buf;           // = gd.rdma_buffer + gd.idx_size
    uint64_t size;       // size of buffer, gd.rdma_bufsize - gd.idx_size
    uint64_t totalsize;  // total size of currently requested PGs (becomes file offset of next PG)
    uint64_t co;         // current offset to put next block to it 
    uint64_t eob;        // end of used buffer (a block may not fit into the end)
    uint64_t pso, peo;   // protected zone (containing all blocks under processing)
                         // pso will be modified (moved ahead) by another thread!
                         // We use no mutex as we only read it in this thread and 
                         // it is moved only forward in the cyclic buffer by the other thread.

    // counters/pointers for the next batch of blocks to be worked on
    uint64_t wso, weo;   // Worker's start and end offset (to be or being written; untouchable region)
    uint64_t woff;       // Current output file offset (Note: wso goes around in cyclic buffer)
    int      wn;         // n'th in 'order' is the first client which is to be written/processed
                         //   (we have sent pull requests for clients n..n+np-1 and will
                         //   start processing when all of them completed)
    int      wnopr;      // number of outstanding (unfinished) pull requests for 
                         //   clients n...n+np-1 (=0: start writing/processing)

    int      nextpull;   // 0..gd.nc-1, the next client in 'order' to pull data from.
                         //   It goes ahead beyond n+np-1, as one new pull is issued
                         //   at each completed pull
    int      nopr;       // number of outstanding (unfinished) pull requests for all clients
} cb; // to handle the rdma buffer in cycle

static bool worker_inited = false;
void worker_init (void)
{
    if (worker_inited)
        return;
    woq = queue_init(QUEUE_UNBOUNDED);
    if (!woq) {
        log_error("rank %d: ERROR: could not create worker-writer queue\n", gd.mpi_rank);
        MPI_Abort(gd.mpi_comm, 1);
    }

    gd.prec_writer_indexing = precedence_create(1); // Writer = 0 and Indexing = 1 in precedence order; 
    if (!gd.prec_writer_indexing) {
        log_error("ERROR: rank %d: Cannot create precedence structure\n", gd.mpi_rank);
        MPI_Abort(gd.mpi_comm, 1);
    }

#if HAVE_PTHREAD
    int rc;
    /* Start bpwriter thread. Will join in finalize */
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    log_debug("rank %d: create writer thread\n", gd.mpi_rank);
    rc = pthread_create(&writer_thread, &attr, bpwriter_thread_main, NULL);
    if (rc) {
        log_error("ERROR: rank %d: Cannot create writer thread, err code = %d\n", gd.mpi_rank, rc);
        MPI_Abort(gd.mpi_comm, 1);
    }
    /* Start bpindexing thread. Will join in finalize */
    
    log_debug("rank %d: create indexing thread\n", gd.mpi_rank);
    rc = pthread_create(&indexing_thread, &attr, bpindexing_thread_main, NULL);
    if (rc) {
        log_error("ERROR: rank %d: Cannot create indexing thread, err code = %d\n", gd.mpi_rank, rc);
        MPI_Abort(gd.mpi_comm, 1);
    }
    
    pthread_attr_destroy(&attr);
#else
    bpwriter_init();
#endif
    gd.order = NULL;
    worker_inited = true;
    worker_status = WORKER_READY;
}

void worker_finalize (void)
{
    if (gd.order)
        free(gd.order);
#if HAVE_PTHREAD
    void *status;
    int rc;

    /* Request writer thread to exit */
    struct writer_request * wr =
        (struct writer_request*) malloc (sizeof(struct writer_request));
    if (wr) {
        wr->type = WRITER_REQUEST_EXIT;
        enqueue (woq, wr);
    }

    // wake up indexing thread and tell to exit
    bpindexing_request_doexit();
    
    // join writer thread
    rc = pthread_join(writer_thread, &status);
    if (rc) {
        log_error("ERROR: rank %d: Cannot join writer thread, err code = %d\n", 
                gd.mpi_rank, rc);
    } else {
        log_debug("rank %d: Joined writer thread.\n", gd.mpi_rank);
    }

    // join indexing thread
    rc = pthread_join(indexing_thread, &status);
    if (rc) {
        log_error("ERROR: rank %d: Cannot join indexing thread, err code = %d\n", 
                gd.mpi_rank, rc);
    } else {
        log_debug("rank %d: Joined indexing thread.\n", gd.mpi_rank);
    }
   
#else  
    // loop until there are requests to be processed in writer
    while (!bpwriter_main());
    bpwriter_finalize();
    queue_destroy (woq);
#endif
}


/* local functions */
static int calc_pull_number(uint64_t maxpgsize);
static int * calc_pull_order(void);

int worker_main (void *arg)
{
    if (!worker_inited) // non-threaded env: transport calls this function only
        worker_init(); 

    int num_event = 0, i;
    bool call_writer_main = true; // call writer after event handling (non-threaded env)
    struct event *e = dequeue(gd.twq);
    if (e != NULL) {

        switch (e->type) {

        case EVENT_CLIENT_REQUEST:

            /* FIXME: this should not be handled until worker_status != WORKER_READY
                it should be postponed and actually the auto-index-pulling by 
                the transport thread should be postponed too!!! */
    
            /* One client request = all client processes, connected to this server,
                sent a request (same file/timestep), and the transport thread has
                pulled indices from them if gd.need_index==1 */
            log_info("rank %d: request from clients nc=%d: timestep=%d index-size=%lld path=%s\n",
                    gd.mpi_rank, 
                    e->event.client_request.number_of_clients,
                    e->event.client_request.timestep,
                    gd.index_size, 
                    e->event.client_request.path);
            uint64_t maxpgsize=0;
            /*Note: gd.nc = e->event.client_request.number_of_clients */
            for (i=0; i<e->event.client_request.number_of_clients; i++) {
                log_debug("rank %d: info client %d: rank=%d, nid=%d, idx_size=%lld, idx_offset=%lld"
                        "  pg_size=%lld\n",
                        gd.mpi_rank, 
                        i,
                        gd.clientdata[i].rank,
                        gd.clientdata[i].nid,
                        gd.clientdata[i].idx_size,
                        gd.clientdata[i].idx_offset,
                        gd.clientdata[i].pg_size
                );
                if (maxpgsize < gd.clientdata[i].pg_size)
                    maxpgsize = gd.clientdata[i].pg_size;
            }

            /*Sanity check: gd.nc = e->event.client_request.number_of_clients */
            if (gd.nc != e->event.client_request.number_of_clients) {
                log_error("rank %d: RDMA ERROR: number of connected clients != "
                          "number of clients that has sent a request (%d and %d, resp.).\n"
                          "Behavior of staging code is undefined from here\n",
                          gd.mpi_rank, gd.nc, e->event.client_request.number_of_clients);
            }

            np = calc_pull_number (maxpgsize);
            if (maxpgsize > gd.rdma_bufsize - gd.index_size) {
                log_error("rank %d: RDMA Error: the largest block from one process "
                          "(%lld bytes) does not fit into available RDMA buffer (%lld bytes)\n",
                          gd.mpi_rank, maxpgsize, gd.rdma_bufsize - gd.index_size);
                np = -1;
            }
            if (np == 1 && maxpgsize > (gd.rdma_bufsize - gd.index_size)/2) {
                log_warn("rank %d: RDMA Warning: the largest block from one process "
                          "(%lld bytes) occupies more than half of available RDMA buffer "
                          "(%lld bytes). Data pull and file writing will not overlap.\n",
                          gd.mpi_rank, maxpgsize, (gd.rdma_bufsize - gd.index_size)/2);
            }

            log_info("rank %d: pulling strategy: %d pulls, at once %d, largest block %lld bytes\n", 
                     gd.mpi_rank, gd.nc, np, maxpgsize);

            gd.order = calc_pull_order();
            for (i=0; i<gd.nc; i++) {
                gd.clientdata[ gd.order[i] ].order_idx = i;
                gd.clientdata[i].status = PULLSTAT_INIT;
            }

            /* Init cyclic buffer */
            cb.buf  = gd.rdma_buffer  + gd.index_size;
            cb.size = gd.rdma_bufsize - gd.index_size;
            cb.totalsize = 0; 
            cb.co = cb.wso = cb.weo = cb.woff = 0; 
            cb.pso = cb.peo = 0;
            cb.eob = cb.size;
            worker_status = WORKER_PULLING;

            /* Enqueue first np pull requests*/
            for (i=0; i<np; i++) {
                struct event *ep = (struct event *) malloc (sizeof(struct event));
                struct event_pull_request *r = &ep->event.pull_request;
                struct globals_client_data *cd = &gd.clientdata[gd.order[i]];
                ep->type  = EVENT_PULL_REQUEST;
                r->rank   = cd->rank;
                r->size   = cd->pg_size;
                r->offset = gd.index_size+cb.weo; // offset in gd.rdma_buffer

                log_debug("rank %d: enqueue pull request, client=%d, size=%d, offset=%lld\n",
                          gd.mpi_rank, r->rank, r->size, r->offset);
                enqueue (gd.wtq, ep);

                cd->status = PULLSTAT_REQUESTED;
                cd->pg_offset = cb.weo;
                cd->file_offset = cb.totalsize;
                cb.weo += cd->pg_size;
                cb.totalsize += cd->pg_size;
                cb.co = cb.weo;
            }

            cb.wn = 0;   /* first pull is order[0] until we get the first np pulls done */
            cb.wnopr = np; /* number of outstanding pull request for 0..np-1 clients */
            cb.nextpull = np; /* next client is order[np] */
            cb.nopr = np; /* number of outstanding pull request for all clients */

            /* Start index processing thread here. */
            bpindexing_request_doindex(); // => indexing_data.indexing_completed=false

            /* Request opening the output file */
            struct writer_request * wr = 
                   (struct writer_request*) malloc (sizeof(struct writer_request));
            if (wr) {
                wr->type = WRITER_REQUEST_OPEN;
                wr->request.open.path = strdup(e->event.client_request.path);
                enqueue (woq, wr);
            }

            /* postpone calling writer for next worker cycle.
               First, pull events should be processed in transport cycle */
            call_writer_main = false; 
            break;

        case EVENT_PULL_COMPLETED:
            log_debug("rank %d: pull completed for client %d. Order=%d\n",
                    gd.mpi_rank, e->event.pull_completed.rank, 
                    gd.clientdata[e->event.pull_completed.rank-gd.lrank].order_idx
                    );

            int cid = e->event.pull_completed.rank - gd.lrank;
            int idx = gd.clientdata[cid].order_idx;
            gd.clientdata[cid].status = PULLSTAT_DONE;
            cb.nopr--;
            if (cid < cb.wn+np) {
                // one pull done for that np requests that we want to process next
                cb.wnopr--;
            }

            while (cb.wnopr == 0) {
                /* np contiguous pull request done: request writing */
                struct writer_request * wr = 
                    (struct writer_request*) malloc (sizeof(struct writer_request));
                struct writer_request * wr2 = NULL;
                if (wr) {
                    cb.peo = cb.weo; // add this batch to the protected zone

                    if (np == 1) {
                        /* special case, one block is pulled ->
                           pso = peo = weo would designate non-protected (empty) and
                           all-protected (filled with one block being written) 
                           => set pso to 0 instead of weo
                        */
                        if (cb.pso == cb.peo && cb.peo == cb.weo) {
                            /*
                            if (cb.pso != 0 && cb.pso != cb.weo) {
                                // sanity check
                                log_error("rank %d: Code error for special case for np=1 handling of pso. "
                                        "Assumption did not hold: 0 != cb.pso=%lld "
                                        "!= cb.weo=%lld, cb.peo=%lld, nc=%d, nextpull=%d\n",
                                        gd.mpi_rank, cb.pso, cb.weo, cb.peo, gd.nc, cb.nextpull);
                            }
                            */
                            log_debug("rank %d: special case for np=1 pso=peo=weo=%lld -> set pso to 0\n",
                                       gd.mpi_rank, cb.pso);
                            cb.pso = 0;
                        }
                    } 

                    wr->type = WRITER_REQUEST_WRITE;
                    wr->request.write.addr = (uint64_t)cb.buf + cb.wso;
                    wr->request.write.offset = cb.woff;
                    wr->request.write.pso = &cb.pso;
                    if (cb.weo < cb.wso) {
                        wr->request.write.size = cb.eob - cb.wso;
                        wr->request.write.weo = cb.eob;
                    } else {
                        wr->request.write.size = cb.weo - cb.wso;
                        wr->request.write.weo = cb.weo;
                    }
                    log_debug("rank %d: write request: addr=%lld size=%lld offset=%lld "
                            "wso=%lld weo=%lld eob=%lld\n",
                            gd.mpi_rank, wr->request.write.addr, wr->request.write.size, 
                            wr->request.write.offset, cb.wso, cb.weo, cb.eob);
                    enqueue (woq, wr);
                    cb.woff += wr->request.write.size;

                    if (cb.weo < cb.wso) {
                        /* Add second half of split buffer as additional request */
                        wr2 = (struct writer_request*) malloc (sizeof(struct writer_request));
                        if (wr2) {
                            wr2->type = WRITER_REQUEST_WRITE;
                            wr2->request.write.addr = (uint64_t)cb.buf;
                            wr2->request.write.offset = cb.woff;
                            wr2->request.write.size = cb.weo;
                            wr2->request.write.pso = &cb.pso;
                            wr2->request.write.weo = cb.weo;
                            log_debug("rank %d: write request: addr=%lld size=%lld offset=%lld "
                                    "wso=%lld weo=%lld eob=%lld\n",
                                    gd.mpi_rank, wr2->request.write.addr, wr2->request.write.size, 
                                    wr2->request.write.offset, cb.wso, cb.weo, cb.eob);
                            enqueue (woq, wr2);
                            cb.woff += wr2->request.write.size;
                        } else {
                            log_error("rank %d: ERROR: Cannot allocate memory for a message of %d bytes. "
                                    "Undefined behavior from now...\n",
                                    gd.mpi_rank, sizeof(struct writer_request));
                        }
                    }
                    // cb.pso will be moved forward by writer thread when finished writing
                } else {
                    log_error("rank %d: ERROR: Cannot allocate memory for a message of %d bytes. "
                              "Undefined behavior from now...\n",
                              gd.mpi_rank, sizeof(struct writer_request));
                }

                /* Update cb.w* counters for next batch */
                if (cb.wn+np >= gd.nc) {
                    // this was the last batch 
                    cb.wso = cb.weo = -1;
                    cb.wn = gd.nc;
                    worker_status = WORKER_WORKING; // no more pulling, all writes requested
                    log_debug("rank %d: worker status = Working (pulling is done)\n", gd.mpi_rank);
                    break; // break out the while loop
                } else {
                    cb.wn = cb.wn+np;
                    if (cb.weo + gd.clientdata[gd.order[cb.wn]].pg_size < cb.size) {
                        cb.wso = cb.weo; // next pg fits after cb.weo
                    } else {
                        cb.wso = 0;  // jump to beginning of buffer
                    }
                    cb.weo = cb.wso; // init weo, increment in loop below

                    /* Calculate weo and wnopr for the next batch */
                    int nextwn = MIN (cb.wn+np, gd.nc);
                    for (i = cb.wn; i < nextwn; i++) { 
                        struct globals_client_data *cd = &gd.clientdata[gd.order[i]];
                        if (cd->status != PULLSTAT_DONE) {
                            // this client is not yet pulled
                            cb.wnopr++;
                        }

                        if (cb.weo + cd->pg_size < cb.size) {
                            cb.weo += cd->pg_size;
                        } else {
                            // split batch in cyclic buffer
                            cb.eob = cb.weo;      
                            cb.weo = cd->pg_size; 
                            // this batch will make two write request
                            // [wso..eob] and [0..weo]
                        }
                    }
                    //FIXME: Why was this here????
                    // if (!cb.wnopr && !cb.nopr)
                    //    worker_status = WORKER_WORKING; // no more pulling
                    
                }

            } // while (cb.wnopr == 0)

            break;

        default:
            break;
        }
        free(e);
        num_event++;
    }

    if (worker_status == WORKER_WORKING && indexing_data.indexing_completed) {
        // pulling is done, all write requests are sent and indexing is done
        // --> send index write request, then file close request
        // and we are ready (although writer thread is still working...)

        /* Note: we don't use mutex on indexing_data, since the indexing_completed
           flag ensures that the other thread is not doing anything with it */

        struct writer_request * wr; 

        // send index write request 
        wr = (struct writer_request*) malloc (sizeof(struct writer_request));
        if (wr) {
            wr->type = WRITER_REQUEST_WRITEIDX;
            wr->request.writeidx.laddr = indexing_data.laddr;
            wr->request.writeidx.lsize = indexing_data.lsize;
            wr->request.writeidx.loffset = cb.woff;
            log_debug("rank %d: local-index write request: addr=%lld size=%lld offset=%lld\n",
                    gd.mpi_rank, wr->request.writeidx.laddr, wr->request.writeidx.lsize, 
                    wr->request.writeidx.loffset);
            if (gd.mpi_rank == 0) {
                // add global index write request on rank 0
                wr->request.writeidx.gaddr = indexing_data.gaddr;
                wr->request.writeidx.gsize = indexing_data.gsize;
                log_debug("rank %d: global-index write request: addr=%lld size=%lld\n",
                          gd.mpi_rank, wr->request.writeidx.gaddr, wr->request.writeidx.gsize);
            }
            enqueue (woq, wr);
        }

        // send close request
        wr = (struct writer_request*) malloc (sizeof(struct writer_request));
        if (wr) {
            wr->type = WRITER_REQUEST_CLOSE;
            log_debug("rank %d: write request: close file\n", gd.mpi_rank);
            enqueue (woq, wr);
        }

        worker_status = WORKER_READY;
        log_debug("rank %d: worker status = Ready (writing may still be going on)\n", gd.mpi_rank);

        /* send acknowledgement to all clients */
        struct event *ep = (struct event *) malloc (sizeof(struct event));
        struct event_ack_request *r = &ep->event.ack_request;
        ep->type  = EVENT_ACK_REQUEST;
        r->status = 0;

        log_debug("rank %d: enqueue ack request\n", gd.mpi_rank);
        enqueue (gd.wtq, ep);
    }

    /* Fill up pull requests up to np clients (gd.order[cb.nextpull]) 
       if there is more space in the buffer and something to pull */
    if (worker_status != WORKER_READY) {
        while (cb.nopr < np && cb.nextpull < gd.nc) {
            /* NOTE: cb.pso can change any time (by writer thread)
                     but it cannot bother this code
            */
            struct globals_client_data *cd = &gd.clientdata[gd.order[cb.nextpull]];
            log_debug("rank %d: filling loop next=%d co=%lld pso=%lld nextsize=%lld\n",
                    gd.mpi_rank, cb.nextpull, cb.co, cb.pso, cd->pg_size);

            uint64_t nextco = cb.co + cd->pg_size;
            /* break the loop if we have no more space to pull in */
            // 1. cyclic buffer, co <= pso and block does not fit before the protected zone, 
            //    co <= pso < co + next block's size < buffer size
            if (cb.co <= cb.pso && cb.pso < nextco && nextco < cb.size) 
                break;
            // 2. cyclic buffer, peo < co and block does not fit either at the end 
            //    of buffer or from the beginning but before the protected zone
            if (nextco > cb.size && cb.pso < cd->pg_size)
                break;

            // Note: wso < co < weo  ||   weo < wso < co  || co < weo < wso
            // are okay, this means we still need to request pull for the batch 
            // considered for the next write request
            struct event *ep = (struct event *) malloc (sizeof(struct event));
            struct event_pull_request *r = &ep->event.pull_request;

            // do we need to go around in cyclic buffer?
            if (cb.co + cd->pg_size > cb.size) {
                //cb.eob = cb.co; // set correctly before use above in loop for cb.w* update
                cb.co = 0;
            }

            ep->type  = EVENT_PULL_REQUEST;
            r->rank   = cd->rank;
            r->size   = cd->pg_size;
            r->offset = gd.index_size+cb.co; // offset in gd.rdma_buffer

            log_debug("rank %d: enqueue pull request, client=%d, size=%d, offset=%lld "
                      "co=%lld\n",
                    gd.mpi_rank, r->rank, r->size, r->offset, cb.co);
            enqueue (gd.wtq, ep);

            cd->status = PULLSTAT_REQUESTED;
            cd->pg_offset = cb.co;
            cd->file_offset = cb.totalsize;
            cb.co += cd->pg_size;
            cb.totalsize += cd->pg_size;
            cb.nextpull++;
            cb.nopr++;
        }
    }

#if ! HAVE_PTHREAD
    if (call_writer_main) { 
        num_event += bpwriter_main();
    }
#endif

    return num_event;
}


/* Calculate how many PGs to pull at once */
static int calc_pull_number(uint64_t maxpgsize)
{
    /* h  = half of available rdma buffersize 
       b  = min (h, max block size we intend to use)
       s  = max (b, max PG size we have to pull in one)
    */
    uint64_t h = (gd.rdma_bufsize - gd.index_size) / 2;
    uint64_t b = MIN (user_max_block_size, h);
    uint64_t s = MAX (maxpgsize, b);
    int n = s / maxpgsize;
    if (n > gd.nc)
        n = gd.nc;
    if (n > user_max_client_pulls)
        n = user_max_client_pulls;
    log_debug("rank %d: calc_pull_number h=%lld, b=%lld, s=%lld, n=%d\n", gd.mpi_rank, h,b,s,n)
    return n;
}

struct nodestr {
    int idx; // idx of gd.clientdata
    int nid; // node id of gd.clientdata[idx]
};

/* compare function for nodestr (compares nid) */
static int nodecmp (const void *p1, const void *p2)
{
    return ( ((struct nodestr *)p1)->nid - ((struct nodestr *)p2)->nid );
}

#define PRINT_ARRAY(a, n, name) log_debug("rank %d: %s=[", gd.mpi_rank, name); \
        for (i=0; i<n; i++) { \
            log_debug("%d ", a[i]); \
        }\
        log_debug("]\n");

/* Determine strategy for ordering of pulls.
   Give an ordering of gd.clientdata[] to maximize the number of different nodes 
   in any interval in gd.clientdata[].
*/
static int * calc_pull_order(void)
{
    struct nodestr *n; // node ids sorted incrementally plus index back to gd.clientdata
    int * p; // moving pointers to separate nodes in n
    int * p1; // pointers to separate nodes in n (fixed pointers to compare p to)
    int * o; // the order for pulling strategy 
    int i, k;
    int prev;  
    int nodes; // number of different nodes

    n = (struct nodestr *) malloc (gd.nc*sizeof(struct nodestr));
    o = (int*) malloc (gd.nc*sizeof(int));

    /* init arrays 
       E.g: 
           n.idx=[0 1 2 3 4 5 6 ]
           n.nid=[1406 1407 1395 1402 1406 1407 1395 ]
    */
    for (i=0; i<gd.nc; i++) {
        n[i].idx = i;
        n[i].nid = gd.clientdata[i].nid;
    }
    memset (o, -1, gd.nc*sizeof(int));

    /* sort n.nid in ascending order of n.nid 
       E.g: 
           n.idx=[2 6 3 0 4 1 5 ]
           n.nid=[1395 1395 1402 1406 1406 1407 1407 ]
    */
    qsort (n, gd.nc, sizeof(struct nodestr), nodecmp);

    
    log_debug("rank %d: n.idx=[", gd.mpi_rank);
    for (i=0; i<gd.nc; i++) {
        log_debug("%d ", n[i].idx);
    }
    log_debug("]\n");

    log_debug("rank %d: n.nid=[", gd.mpi_rank);
    for (i=0; i<gd.nc; i++) {
        log_debug("%d ", n[i].nid);
    }
    log_debug("]\n");
    

    /* first round on p and o: get first list of separate nodes 
       E.g: 
            nodes=4
            p=[0 2 3 5 7 7 7 ]     7=gd.nc=invalid, 4 different nodes on idx 0,2,3,5
            o=[0 2 3 5 -1 -1 -1 ]
    */
    prev=-1;
    nodes=0; // number of different nodes
    for (i=0; i<gd.nc; i++) {
        if (n[i].nid != prev) {
            o[nodes] = i; // = p[nodes]
            prev = n[i].nid;
            nodes++;
        }
    }

    p = (int*) malloc ((nodes+1)*sizeof(int));
    p1 = (int*) malloc ((nodes+1)*sizeof(int));
    for (i=0; i<nodes; i++) {
        p[i] = o[i];
        p1[i] = o[i];
    }
    p[nodes] = p1[nodes] = gd.nc; // bigger than any id; for comparisons
    /* 
    PRINT_ARRAY(p1,nodes,"p1");
    PRINT_ARRAY(p,nodes,"p ");
    PRINT_ARRAY(o,nodes,"o ");
    */
    /* determine next rounds of nodes until all clients are processed */
    k = nodes;
    while (k < gd.nc) {
        
        /* increment p */
        for (i=0; i<nodes; i++) {
            if (p[i]+1 < p1[i+1]) { // note: p is large enough to test p[i+1]
                // we have more of this nid
                p[i]++;
                o[k] = p[i];
                k++;
            //} else {
            //    p[i] = gd.nc; // no more of this nid
            }
        }
        /*
        PRINT_ARRAY(p1,nodes,"p1");
        PRINT_ARRAY(p,nodes,"p ");
        PRINT_ARRAY(o,k,"o ");
        */
    }

    /* result is o ordering of n.idx to maximize number of nodes per pull requests
       E.g.
            p=[1 7 4 6 7 7 7 ]
            o=[0 2 3 5 1 4 6 ]
    */
    

    /* Translate ordering back to gd.clientdata indices */
    for (i=0; i<gd.nc; i++) 
        o[i] = n[o[i]].idx;

    log_debug("rank %d: calc_pull_order nodes=%d o=[", gd.mpi_rank, nodes);
    for (i=0; i<gd.nc; i++) {
        log_debug("%d ", o[i]);
    }
    log_debug("]\n");

    free(p);
    free(p1);
    free(n);
    return o;
}



