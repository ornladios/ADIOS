#include "globals.h"
#include <unistd.h>
#ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 199309 
#endif
#include <time.h>   // nanosleep
#include <string.h>   // memset

#include "ardma_server.h"

#if ! HAVE_PTHREAD
#include "worker_thread.h"
#endif
#include "mem.h"
#include "events.h"


static struct ardma_server_connection *asc;
static struct ardma_memory *mem;
//static int terminate; // = 1: should exit due failure

int check_worker_events(void);

/* This is the main thread, so do not call pthread_exit!!! */
void transport_thread_main (void *arg)
{
    int num_event;
    struct timespec treq = {.tv_sec=0, .tv_nsec=2000000}; // 2ms sleeptime
    struct timespec trem;
    log_info("rank %d: initialize rdma server\n",gd.mpi_rank);
    asc = ardma_server_init(gd.mpi_comm, verbose); 
    if (asc) {

        log_debug("rank %d: asc=%x, number of clients currently=%d, expected=%d, lowest rank=%d\n",
                gd.mpi_rank,asc,asc->nc_current,asc->nc,asc->lrank);
        // When number of clients is 0, we finish work (it is -1 in the beginning!).
        // Should return checking events as fast as possible
        while (!gd.terminate && gd.nc_current) {
            num_event = ardma_server_check_events(asc);
            //log_debug("rank %d: nc current = %d\n", gd.mpi_rank, asc->nc_current);
#if ! HAVE_PTHREAD
            // call worker main if not running as a separate thread
            num_event += worker_main (NULL);
#endif
            num_event += check_worker_events();
            if (num_event == 0)
                nanosleep(&treq, &trem); 
        }
    }

    gd.application_finished = true;
    if (gd.terminate) {
        log_info("rank %d: transport thread: exit due to failure.\n", gd.mpi_rank);
    } else {
        log_info("rank %d: transport thread: application has exited and no more work to do.\n", gd.mpi_rank);
    }
}


/* Called from main thread after all threads joined */
void transport_finalize(void) {
    ardma_server_deregister_memory(mem);
    ardma_server_finalize(&asc);
    log_debug("rank %d: transport thread finalized\n", gd.mpi_rank);
#if ! HAVE_PTHREAD
    // call worker finalize if not running as a separate thread
    worker_finalize();
#endif
}

/*
    Callbacks from ardma layer when checking for events 
*/

void ardma_server_cb_connect (struct ardma_server_connection * asc, 
            int nc, int nc_total, int crank)
{
    /* At first connection, allocate the buffer.
       RDMA layer is not fully initalized until the first connection, 
       so we register memory after that */
    static int first_connect = 1;
    log_debug("rank %d: callback: a client connected\n", gd.mpi_rank);
    if (first_connect) {
        mem = ardma_server_register_memory (gd.requested_rdma_buffer_size); 
        if (!mem) {
            gd.terminate = true;
            log_error("rank %d: Failed to register memory, start termination\n", gd.mpi_rank);
        } else {
            gd.rdma_buffer = mem->buf;
            gd.rdma_bufsize = mem->size;
        }
        gd.nc_current = 0;
        gd.nc = nc;
        gd.nc_total = nc_total;
        gd.lrank = asc->lrank;
        gd.clientdata = (struct globals_client_data *) malloc (nc*sizeof(struct globals_client_data));
        if (!gd.clientdata) {
            gd.terminate = true;
            log_error("rank %d: Failed to allocate %d bytes for client data. start termination\n", 
                        gd.mpi_rank, nc*sizeof(struct globals_client_data));
        }
        memset(gd.clientdata, 0, nc*sizeof(struct globals_client_data));
    }
    first_connect = 0;
    gd.nc_current++;
}

void ardma_server_cb_disconnect (struct ardma_server_connection * asc, int crank)
{
    gd.nc_current--;
}


// variables for pulling nc indices from clients for 'nc' client requests and then just
// send 1 request to worker.
static int nrequests = 0;
static int nrequests_pulledidx = 0;
static char tr_path[MAXPATHLEN];
static int tr_timestep;


void ardma_server_cb_request (struct ardma_server_connection * asc, int crank, int nodeid,
        uint64_t pg_size, uint64_t idx_size, 
        int timestep, char * path)
{
    /* 
    log_debug("rank %d: client rank=%lld request:\n"
            "  PG  (size=%lld, addr=%lld, key=%u)\n"
            "  IDX (size=%lld, addr=%lld, key=%u)\n"
            "  timestep = %d\n"
            "  path     = %s\n",
            gd.mpi_rank, crank,
            pg_size, pg_addr, *(uint32_t*)pg_rkey,
            idx_size, idx_addr, *(uint32_t*)idx_rkey,
            timestep, path);
    */
     
    struct globals_client_data d = {
        .rank = crank, 
        .pg_size = pg_size,
        .idx_size = idx_size,
        .idx_offset = gd.index_size,
        .nid = nodeid
    };
    memcpy(&gd.clientdata[crank-asc->lrank], &d, sizeof(struct globals_client_data));

    strncpy(tr_path, path, MAXPATHLEN);
    tr_timestep = timestep;
    nrequests++;
    if (gd.need_index) {
        ardma_server_pull_index (asc, crank, mem, (uint64_t)gd.index_size);

        gd.index_size += idx_size;
    } else {
        ardma_server_cb_pulled_index (asc, crank);
    }
}

void ardma_server_cb_pulled_index (struct ardma_server_connection * asc, int crank)
{
    nrequests_pulledidx++;

    log_debug("rank %d: cb_pulled_index, nc=%d, nrequests=%d, pulled=%d, path=%s\n", 
                gd.mpi_rank, gd.nc, nrequests, nrequests_pulledidx, tr_path);

    if (nrequests == nrequests_pulledidx && nrequests == gd.nc) {
        // all requests and if needed all indicies are here
        // send on request to worker now
        struct event *e = (struct event *) malloc (sizeof(struct event));
        struct event_client_request *r = &e->event.client_request;

        e->type = EVENT_CLIENT_REQUEST;
        r->number_of_clients = gd.nc;
        strncpy(r->path, tr_path, MAXPATHLEN);
        r->timestep = tr_timestep;

        /*
        log_debug("rank %d: transport sends request to worker, nc=%d:\n"
                "  timestep = %d\n"
                "  path     = %s\n"
                "  index size = %lld\n",
                gd.mpi_rank,
                e->event.client_request.number_of_clients,
                e->event.client_request.timestep,
                e->event.client_request.path,
                gd.index_size);
        */

        enqueue (gd.twq, e);

        nrequests = 0;
        nrequests_pulledidx = 0;
    }
}

void ardma_server_cb_pulled_data (struct ardma_server_connection * asc, int crank)
{
    // send on completion to worker now
    struct event *e = (struct event *) malloc (sizeof(struct event));
    struct event_pull_completed *r = &e->event.pull_completed;

    e->type = EVENT_PULL_COMPLETED;
    r->rank = (uint32_t)crank;
    log_debug("rank %d: transport sends pull completion to worker, crank=%d\n",
              gd.mpi_rank, e->event.pull_completed.rank);
    enqueue (gd.twq, e);
}

void ardma_server_cb_failure (struct ardma_server_connection * asc)
{
    log_error("rank %d: Failure in ARDMA layer. Start terminating\n", gd.mpi_rank);
    gd.terminate = true;
}


/**************************/
/* Handle Worker requests */
/**************************/

int check_worker_events(void)
{
    int nevents=0;
    struct event *e = dequeue(gd.wtq); 
    while (e != NULL) { 
        switch (e->type) { 
        case EVENT_PULL_REQUEST:
            log_debug("rank %d: Worker requested pull: client=%d, size=%d, offset=%lld\n",
                    gd.mpi_rank, e->event.pull_request.rank, e->event.pull_request.size,
                    e->event.pull_request.offset);

            ardma_server_pull_data (asc, e->event.pull_request.rank, mem, 
                                    e->event.pull_request.offset, e->event.pull_request.size);
            break;

        case EVENT_ACK_REQUEST:
            log_debug("rank %d: Worker requested sending acknowledgment, status=%d\n",
                    gd.mpi_rank, e->event.ack_request.status);

            ardma_server_send_acknowledgement (asc, e->event.ack_request.status);
            break;

        default:
            log_error("rank %d: RDMA ERROR: Undefined event from Worker thread to Transport thread (%d).\n", 
                      gd.mpi_rank, (int)e->type );
            break;
        }
        free(e);
        e = dequeue(gd.wtq);
        nevents++;
    }
    return nevents;
}



