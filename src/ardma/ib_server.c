/* 
   ADIOS RDMA Server library for Infiniband

   Server applications should use these functions only to 
   get data from staging methods. 

*/

#include <arpa/inet.h>
#include <sys/socket.h>
#include <netdb.h>
#include <ifaddrs.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <poll.h>
#include <rdma/rdma_cma.h>

#include "ardma_server.h"
#include "ardma_common.h"
#include "ib_common.h"

enum COMPLETION_WAIT {
    COMPL_READY,            // not waiting for any send completion events
    COMPL_PULL_INDEX,       // waiting for index pull completion
    COMPL_PULL_DATA,        // waiting for data pull completion
    COMPL_ACK               // waiting for completion of acknowledgement message 
};

const int DEFAULT_BLOCKSIZE = 1024*1024*2; // chop up large buffer pulls into separate blocks
static int blocksize; // a variable for block size, but currently not user-modifiable

// variables, one per connected client
// we use this struct in work requests and get as id->context in WC completions
struct connection {
    struct rdma_cm_id *id;  // this is id->context->id, equals to top id
    int    connected;       // = 1 if connection is established
    int    crank;           // client's rank
    uint64_t pg_addr;       // remote address of buffer PG
    uint64_t pg_size;       // size of buffer PG
    uint32_t pg_rkey;       // RDMA rkey for pulling buffer PG
    uint64_t pg_cksum;      // hash of PG buffer (sent by client)
    uint64_t idx_addr;      // remote address of buffer IDX
    uint64_t idx_size;      // size of buffer IDX 
    uint32_t idx_rkey;      // RDMA rkey for pulling buffer IDX
    enum COMPLETION_WAIT waitingfor; // to identify send completion events
    uint64_t pg_local_addr; // address of PG in local memory after pull started
    int     needs_ack;      // 1: acknowledgement should be sent, 0: ack was already sent
    // pulling PG in blocks:
    int      nblocks;       // number of blocks to pull (is incremented in ardma_server_pull_block())
    uint64_t pg_offset;     // offset of next block to be pulled
    uint64_t lkey;          // local target's RDMA key
};

// server connection variables (one instance per server)
static struct   sockaddr_in          addr;
static struct   rdma_cm_event      * event = NULL;
static struct   rdma_cm_id         * listener = NULL;
static struct   rdma_event_channel * ec = NULL;
static uint16_t port;
static struct   ibv_pd             * pd = NULL;           // protection domain 
static struct   ibv_cq             * cq = NULL;           // completion queue
static struct   ibv_comp_channel   * comp_channel = NULL;
static struct   ardma_message      * recvbuf = NULL;      // buffer for max 'asc->nc' client messages
static struct   ibv_mr             * recvbuf_mr = NULL;   // MR for recv message buffer
static struct   ardma_message      * sendbuf = NULL;      // buffer for max 'asc->nc' server messages to clients
static struct   ibv_mr             * sendbuf_mr = NULL;   // MR for send message buffer
static struct   connection         * conns = NULL;        // connection array
            
// MPI variables of the RDMA server processes 
static MPI_Comm scomm; 
static int rank;
static int size;

FILE *ardma_logf;
int ardma_verbose_level;

// internal functions
static int post_receive (struct ibv_qp *qp, uint64_t wrid, uint64_t addr, 
                         uint32_t size, uint32_t lkey, int rank, int crank);
static int on_connect_request(struct rdma_cm_id *id, 
                              const void * private_data, uint8_t private_data_len,
                              struct ardma_server_connection *asc); 
static int on_connection(struct rdma_cm_id *id, struct ardma_server_connection *asc);
static int on_disconnect(struct rdma_cm_id *id, struct ardma_server_connection *asc);
static int on_cm_event(struct rdma_cm_event *event, struct ardma_server_connection *asc);
static int ardma_server_send_ack1 (struct ardma_server_connection *asc, 
                                   int idx, // client idx, 0..asc->nc-1 
                                   int status);
static int ardma_server_pull_block (struct ardma_server_connection *asc,
                                    struct connection * conn);


/* Initialize ardma 
   Create listening ports and write a file containing the connection info clients can use
*/
struct ardma_server_connection * ardma_server_init (MPI_Comm comm, int verbose_level)
{
    struct ardma_server_connection * asc = 
            (struct ardma_server_connection *) malloc (sizeof(struct ardma_server_connection *));
    memset (asc, 0, sizeof(struct ardma_server_connection));
    asc->nc_current = -1;

    ardma_verbose_level = verbose_level;
    blocksize = DEFAULT_BLOCKSIZE;

    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;

    /* 1. get my Infiniband IP address */
    struct ifaddrs *ifaddr, *ifa;
    int s;
    char host[NI_MAXHOST], *ibip = NULL;

    if (getifaddrs(&ifaddr) == -1) {
        log_error("%s Error: getifaddrs() failed: %s\n", __func__, strerror(errno));
        free(asc);
        return NULL;
    }

    for (ifa = ifaddr; ifa != NULL; ifa = ifa->ifa_next) {
        if (ifa->ifa_addr->sa_family == AF_INET) {
            s = getnameinfo(ifa->ifa_addr, sizeof(struct sockaddr_in),
                    host, NI_MAXHOST, NULL, 0, NI_NUMERICHOST);
            if (s != 0) {
                log_warn("getnameinfo() failed: %s\n", gai_strerror(s));
            }
            log_debug("<Interface>: %s \t <Address> %s\n", ifa->ifa_name, host);
            if (!strncmp(ifa->ifa_name,"ib",2)) {
                ibip = host;
                break;
            }
        }
    }

    if (!ibip) {
        log_error("%s Error: cannot determine IP address of Infiniband device\n", __func__);
        free(asc);
        return NULL;
    }

    /* 2. create RDMA event channel and listener on an arbitrary port */
#define RETURN_NOW(x) log_error("%s, rank %d Error: " #x " failed\n", __func__, rank); \
        free(asc); return NULL;

#define TEST_NZ(x) if ( (x)) {RETURN_NOW(x)};
#define TEST_Z(x)  if (!(x)) {RETURN_NOW(x)};

    TEST_Z(ec = rdma_create_event_channel());

    /* change to non-blocking mode for the completion channel */
    int flags = fcntl(ec->fd, F_GETFL);
    int rc = fcntl(ec->fd, F_SETFL, flags | O_NONBLOCK);
    if (rc < 0) {
        log_error ("rank %d: RDMA Error: Failed to change file descriptor of "
                "event channel to non-blocking mode\n", rank);
        free(asc);
        return NULL;
    }

    TEST_NZ(rdma_create_id(ec, &listener, NULL, RDMA_PS_TCP));
    TEST_NZ(rdma_bind_addr(listener, (struct sockaddr *)&addr));
    TEST_NZ(rdma_listen(listener, 10)); /* ???is backlog not truncated to SOMAXCONN???? */
    port = ntohs(rdma_get_src_port(listener));

    log_debug("listening on Infiniband IP %s port %d.\n", ibip, port);

    /* 3. Collect IP:port from all server process and write into connection file */ 
    MPI_Comm_dup(comm, &scomm);
    MPI_Comm_rank (scomm, &rank);
    MPI_Comm_size (scomm, &size);

    char * ips;
    int  * ports;
    int  port_as_int = (int)port;
    if (rank == 0) {
        ips = (char *) malloc (NI_MAXHOST * size * sizeof(char));
        ports = (int *) malloc (size * sizeof(int));
    }

    MPI_Gather( ibip, NI_MAXHOST, MPI_CHAR, ips, NI_MAXHOST, MPI_CHAR, 0, scomm); 
    MPI_Gather( &port_as_int, 1, MPI_INT, ports, 1, MPI_INT, 0, scomm); 

    if (rank == 0) {
        int i;
        FILE * f;
        f = fopen(ARDMA_CONNECTION_FILENAME,"w");
        if (f) {
            fprintf(f, "%d\n",size);
            for (i=0; i<size; i++) {
                log_info("rank %d: listens on Infiniband IP %s port %d.\n", i, ips+i*NI_MAXHOST, ports[i]);
                fprintf(f, "%s:%d\n", ips+i*NI_MAXHOST, ports[i]);
            }
            fclose(f);
        } else {
            log_error("%s Error: cannot open file %s. %s\n",
                      __func__, ARDMA_CONNECTION_FILENAME, strerror(errno));
            free(asc);
            asc = NULL;
        }

        free(ips);
        free(ports);
    }
    return asc;
}

/* Allocate and register memory to rdma layer. 
Return: ardma_memory struct pointer to the allocated buffer, NULL on error
 */
struct ardma_memory * ardma_server_register_memory (uint64_t size) 
{
    char *buf = alloc_memory_aligned(size);
    if (!buf) {
        log_error("rank %d: RDMA Error: Failed to allocate %lld bytes of memory\n", rank, size);
        return NULL;
    }

    struct ardma_memory * am = (struct ardma_memory *) malloc (sizeof(struct ardma_memory));
    if (!am) {
        log_error("rank %d: RDMA Error: Failed to allocate %lld bytes of memory for ardma_memory struct\n", rank, sizeof(struct ardma_memory));
        return NULL;
    }

    struct ibv_mr * mr = ibv_reg_mr(pd, buf, size, IBV_ACCESS_LOCAL_WRITE);

    if (!mr) {
        log_error("rank %d: RDMA Error: Failed to register %lld bytes of memory to RDMA\n", rank, size);
        free (buf);
        free (am);
        return NULL;
    }
    log_debug("rank %d: registered memory buf=%x, size=%lld mr=%x lkey=%ld\n", 
              rank, buf, size, mr, mr->lkey);

    am->buf = buf;
    am->size = size;
    am->rdma_data = (void *)mr;
    return am;


}

/* Deregister memory. Return 0 on success, !=0 on error */
int ardma_server_deregister_memory (struct ardma_memory *mem) 
{
    if (!mem) return 0;
    if (mem->buf)
        free(mem->buf);
    if (mem->rdma_data)
        ibv_dereg_mr ((struct ibv_mr *)mem->rdma_data);
    free(mem);
    return 0;
}



/* Poll for events.
   First, poll on completion channel for CQ completions and process one.
   If nothing found, poll on RDMA_CM event channel for connection events and process one.
*/
int ardma_server_check_events (struct ardma_server_connection * asc)
{
    struct pollfd p = {
        .fd      = 0,
        .events  = POLLIN,
        .revents = 0 
    };  
    int ms_timeout = 0;
    int rc = 0;
    int num_events = 0;

    /* 
        poll async for a completion event
    */
    if (comp_channel) {
        p.fd = comp_channel->fd;
        rc = poll(&p, 1, ms_timeout);
    } else
        rc = 0;

    if (rc > 0) {
        // get event(s) out from the completion queue
        void * ev_ctx;
        struct ibv_cq *ev_cq;

        log_debug("rank %d: poll on comp channel = %d\n", rank, rc);

        do { // one step loop so that we can jump out

            if (ibv_get_cq_event(comp_channel, &ev_cq, &ev_ctx)) {
                log_error("rank %d: RDMA Error: Failed to get cq_event\n", rank);
                ardma_server_cb_failure(asc);
                continue;  
            }   
            //log_debug("rank %d: Done with ibv_get_cq_event\n", rank);

            ibv_ack_cq_events(ev_cq, 1);
            //log_debug("rank %d: Done with ibv_ack_cq_events\n", rank);

            if (ibv_req_notify_cq(cq, 0)) {
                log_error("rank %d: RDMA Error: Failed to request notification on CQ\n", rank);
                ardma_server_cb_failure(asc);
                continue;  
            }   
            //log_debug("rank %d: Done with ibv_req_notify_cq\n", rank);

            struct ibv_wc wc;
            int ne;
            do {
                ne = ibv_poll_cq(cq, 1, &wc);
                log_debug("rank %d: Done with ibv_poll_cq, ne=%d\n", rank, ne);

                /* there may have been an extra event on comp_channel with no completion in the CQ */
                if (ne == 0) continue;

                /********************/
                /* RECV COMPLETIONS */
                /********************/
                if (wc.opcode & IBV_WC_RECV) {
                    // a client has sent request, index = wc.wr_id
                    log_debug("rank %d: Message from client %d received, wr_id=%lld, slid=%d, ne=%d\n", rank,asc->lrank+(int)wc.wr_id, wc.wr_id, wc.slid, ne);
                    if (wc.status != IBV_WC_SUCCESS) {
                        log_error("rank %d: RDMA Error: Receive failed, status = %d: %s\n",
                                rank, wc.status, ibv_wc_status_str(wc.status));
                        ardma_server_cb_failure(asc);
                        continue;
                    }

                    struct ardma_message * m = &recvbuf[wc.wr_id];

                    if (m->type != ARDMA_MSG_REQUEST) {
                        log_error("rank %d: RDMA Error: Received message is not a request message. type = %d\n",
                                rank, (int)m->type);
                        ardma_server_cb_failure(asc);
                        continue;
                    }
                    
                    struct connection * conn = &conns[wc.wr_id];
                    conn->pg_addr  = (uint64_t) ntohll (m->msg.request.pg_addr);
                    conn->pg_size  = (uint64_t) ntohll (m->msg.request.pg_size);
                    conn->pg_rkey  = (uint32_t) ntohl  (m->msg.request.pg_rkey);
                    conn->pg_cksum = (uint64_t) ntohll (m->msg.request.pg_cksum);
                    conn->idx_addr = (uint64_t) ntohll (m->msg.request.idx_addr);
                    conn->idx_size = (uint64_t) ntohll (m->msg.request.idx_size);
                    conn->idx_rkey = (uint32_t) ntohl  (m->msg.request.idx_rkey);

                    /* 
                    log_debug("rank %d: client %lld request:\n"
                            "  rank = %d, rank calculated = %d\n"
                            "  PG  (size=%lld, addr=%lld, key=%u)\n"
                            "  IDX (size=%lld, addr=%lld, key=%u)\n"
                            "  timestep = %d\n"
                            "  path     = %s\n"
                            "  checksum = %llu\n",
                            rank, wc.wr_id, 
                            ntohl  (m->msg.request.rank), 
                            asc->lrank+(int)wc.wr_id,
                            ntohll (m->msg.request.pg_size), 
                            ntohll (m->msg.request.pg_addr), 
                            ntohl  (m->msg.request.pg_rkey),
                            ntohll (m->msg.request.idx_size), 
                            ntohll (m->msg.request.idx_addr), 
                            ntohl  (m->msg.request.idx_rkey),
                            ntohl  (m->msg.request.timestep), 
                            m->msg.request.path,
                            ntohll (m->msg.request.pg_cksum)
                            );
                    */

                    ardma_server_cb_request (asc, 
                            (int)    ntohl  (m->msg.request.rank),
                            (int)    wc.slid,
                                     conn->pg_size, 
                                     conn->idx_size, 
                            (int)    ntohl  (m->msg.request.timestep),
                            (char *) m->msg.request.path
                            );

                    /* pre-post the next receive for this client's next request */
                    int poststat = post_receive (conn->id->qp, wc.wr_id, 
                            (uint64_t) &recvbuf[wc.wr_id], 
                            sizeof(struct ardma_message), 
                            recvbuf_mr->lkey, rank, conn->crank);

                    if (poststat) return poststat;


                /*************************/
                /* RDMA_READ COMPLETIONS */
                /*************************/
                } else if (wc.opcode == IBV_WC_RDMA_READ) {

                    // a pull completed, wc.wr_id points to connection struct
                    struct connection * conn = (struct connection *)wc.wr_id;
                    uint64_t cksum = 1;
                    log_debug("rank %d: RDMA READ completion about client %d, waitfor=%d, wr_id=%lld, ne=%d\n", 
                              rank, conn->crank, conn->waitingfor, wc.wr_id, ne);
                    if (wc.status != IBV_WC_SUCCESS) {
                        log_error("rank %d: RDMA Error: Completion failed, waitfor=%d, status = %d: %s\n",
                                rank, conn->waitingfor, wc.status, ibv_wc_status_str(wc.status));
                        ardma_server_cb_failure(asc);
                        /*
                        if (conn->waitingfor == COMPL_PULL_DATA) {
                            // send failure ack to client now
                            ardma_server_send_ack1(asc, conn->crank - asc->lrank, 1);
                            conn->needs_ack = 0;
                            // Note: this send would fail anyway if the read already failed
                        }
                        */
                        continue;
                    }

                    switch (conn->waitingfor) {
                    case COMPL_PULL_INDEX:
                        /* Index pull completed */
                        ardma_server_cb_pulled_index (asc, conn->crank);
                        conn->waitingfor = COMPL_READY;
                        break;
                    case COMPL_PULL_DATA:
                        /* data block pull completed */
                        if (conn->pg_offset < conn->pg_size) {
                            // pull next block
                            ardma_server_pull_block (asc, conn);
                        } else {
                            cksum = ardma_calc_checksum ((char *)conn->pg_local_addr, conn->pg_size);
                            if (cksum != conn->pg_cksum) {
                                log_error("rank %d: RDMA Error: Checksum error of pulled PG of client %d: "
                                        "client sent = %llu, calculated locally = %llu\n",
                                        rank, conn->crank, conn->pg_cksum, cksum);
                                /*
                                   int i;
                                   log_debug("rank %d: pulled PG = [", rank);
                                   for (i=0; i<8; i++) 
                                   log_debug("%hhu ",(char) *((char*)(conn->pg_local_addr+i)));
                                   log_debug("... ");
                                   for (i=8; i>0; i--) 
                                   log_debug("%hhu ",(char) *((char*)conn->pg_local_addr+conn->pg_size-i));
                                   log_debug("]\n");
                                 */
                            }

                            /* report back */
                            ardma_server_cb_pulled_data (asc, conn->crank);
                            conn->waitingfor = COMPL_READY;
                            conn->needs_ack = 1;
                        }
                        break;
                    default:
                        log_error("rank %d: RDMA Error: Unexpected send completion event: "
                                  "client %d, waitfor=%d\n",
                                  rank, conn->crank, conn->waitingfor);
                        ardma_server_cb_failure(asc);
                    }

                /*************************/
                /* SEND COMPLETIONS */
                /*************************/
                } else if (wc.opcode == IBV_WC_SEND) {

                    // a send (acknowledgement) completed, wc.wr_id points to connection struct
                    struct connection * conn = (struct connection *)wc.wr_id;
                    log_debug("rank %d: Completion of Send to client %d, waitfor=%d, wr_id=%lld, ne=%d\n", 
                              rank, conn->crank, conn->waitingfor, wc.wr_id, ne);
                    if (wc.status != IBV_WC_SUCCESS) {
                        log_error("rank %d: RDMA Error: Completion of sending a message to client %d "
                                  "failed, waitfor=%d, status = %d: %s\n",
                                rank, conn->crank, conn->waitingfor, wc.status, ibv_wc_status_str(wc.status));
                        ardma_server_cb_failure(asc);
                        continue;
                    }

                    /* we send only acknowledgements, so some sanity check here */
                    if (conn->waitingfor != COMPL_ACK) {
                        log_error("rank %d: RDMA Error: Unexpected Send completion event for client %d. "
                                  "Server is not waiting for completion of acknowledgement message to "
                                  "this client. Waitingfor=%d\n",
                                  rank, conn->crank, conn->waitingfor);
                    }

                    conn->waitingfor = COMPL_READY;

                } else {
                    log_debug("rank %d: Unexpected completion event wr_id=%lld opcode=%d status=%d\n", 
                            rank, wc.wr_id, wc.opcode, wc.status);
                }


            } while(ne); 
        } while(0);


    } else if (rc < 0) {
        log_error ("rank %d: RDMA Error: polling completion channel failed. errno = %d : %s\n",
                rank, errno, strerror(errno));
        ardma_server_cb_failure(asc);
    }


    /*
        poll async for a connection event 
    */
    int rc2;
    p.fd = ec->fd;
    rc2 = poll(&p, 1, ms_timeout);
    if (rc2 > 0) {
        int r = rdma_get_cm_event(ec, &event);
        if (!r) {
            struct rdma_cm_event event_copy;
            struct ardma_connect_data cd;
            memcpy (&event_copy, event, sizeof(*event));

            // for connection requests we have to copy the private data too
            if (event->event == RDMA_CM_EVENT_CONNECT_REQUEST) {
                memcpy(&cd, event->param.conn.private_data, sizeof(struct ardma_connect_data) );
                event_copy.param.conn.private_data = &cd; 
            }

            // mystery: if acknowledge it after on_cm_event, server will hang 
            rdma_ack_cm_event (event);

            on_cm_event(&event_copy, asc);

            if (event_copy.event == RDMA_CM_EVENT_ESTABLISHED) {
                ardma_server_cb_connect(asc, asc->nc, asc->appsize, 
                                        ((struct connection *)event_copy.id->context)->crank);
            } else if (event_copy.event == RDMA_CM_EVENT_DISCONNECTED) {
                ardma_server_cb_disconnect(asc, ((struct connection *)event_copy.id->context)->crank);
            }
        }

    } else if (rc2 < 0) {
        log_error ("rank %d: RDMA Error: polling cm event channel failed. errno = %d : %s\n",
                rank, errno, strerror(errno));
        ardma_server_cb_failure(asc);
    }

    return num_events;
}

static int ardma_server_send_ack1 (struct ardma_server_connection *asc, 
                                   int idx,
                                   int status)
{
    /* prepare message for send */
    struct connection * conn = &conns[idx];

    struct ardma_ack_message ack = {
        .status = status
    };
    sendbuf[idx].type = ARDMA_MSG_REQUEST;
    memcpy (&sendbuf[idx].msg.ack, &ack, sizeof(struct ardma_ack_message));

    struct ibv_sge list = {
        .addr   = (uintptr_t) &sendbuf[idx], // nth slot will contain ack msg for this client
        .length = sizeof(struct ardma_message),
        .lkey   = sendbuf_mr->lkey
    };

    struct ibv_send_wr wr = {
        .wr_id      = (uint64_t) conn,
        .next       = NULL,
        .sg_list    = &list,
        .num_sge    = 1,
        .opcode     = IBV_WR_SEND,
        .send_flags = IBV_SEND_SIGNALED,
    };
    struct ibv_send_wr *bad_wr = NULL;

    //log_debug("wr=%x, bad_wr=%x sge_addr=%x, sge_len=%d sge_key=%d\n",
    //            &wr,bad_wr, wr.sg_list->addr, wr.sg_list->length, wr.sg_list->lkey);

    errno = 0;
    if (ibv_post_send(conn->id->qp, &wr, &bad_wr)) {
        log_error("rank %d: RDMA Error: post_send to send acknowledgement to client %d failed: (errno=%d) %s\n",
                rank, conn->crank, errno, strerror(errno));
        log_debug("wr=%x, bad_wr=%x\n",&wr,bad_wr);
        return 1;
    }

    conn->waitingfor = COMPL_ACK;
    return 0;
}

/* Send acknowledgement message to all clients that the staging completed, 
   i.e. a new request can be sent.
 */
int ardma_server_send_acknowledgement (struct ardma_server_connection *asc, int status)
{
    int i;
    struct connection * conn;
    for (i=0; i<asc->nc; i++) {
        conn = &conns[i];
        if (conn->needs_ack) {
            // send message to those clients only, who received no failure 
            // acknowledgement before
            ardma_server_send_ack1(asc, i, status);
        }
    }
}


static int ardma_server_post_send (uint64_t laddr, uint64_t size, uint32_t lkey, 
        uint64_t raddr, uint32_t rkey, 
        struct connection * conn )
{
    /* post a send for pulling some data from this client */
    struct ibv_sge list = {
        .addr   = (uintptr_t) (laddr), // local address
        .length = size, 
        .lkey   = lkey
    };
    struct ibv_send_wr wr = {
        .wr_id               = (uint64_t) conn,
        .opcode              = IBV_WR_RDMA_READ,
        .send_flags          = IBV_SEND_SIGNALED,
        .wr.rdma.rkey        = rkey,
        .wr.rdma.remote_addr = raddr,
        .sg_list             = &list,
        .num_sge             = 1,
        .next                = NULL
    };
    struct ibv_send_wr *bad_wr;

    if (ibv_post_send(conn->id->qp, &wr, &bad_wr)) {
            log_error("rank %d: RDMA Error: Couldn't post send for RDMA pull from memory of client %d: %s\n",
                      rank, conn->crank, strerror(errno));
            return 1;
    }

}

/* Pull the index from client */
int ardma_server_pull_index (struct ardma_server_connection *asc, int crank, 
                             struct ardma_memory *mem, uint64_t offset)
{
    int idx = crank - asc->lrank;  // 0..asc->nc-1
    struct connection * conn = &conns[idx];
    struct ibv_mr * mr = (struct ibv_mr *)mem->rdma_data;
    int retval;
 
    log_debug("rank %d: pull index from client %d local rank=%d idx size=%lld "
              "raddr=%lld laddr=%lld\n", 
              rank, crank, idx, conn->idx_size, conn->idx_addr, (uint64_t)mem->buf+offset );

    if (offset+conn->idx_size > mem->size) {
        log_error("rank %d: RDMA Error: index does not fit into server's memory, "
                  "client %d idx size=%lld\n mem size=%lld offset=%lld\n", 
              rank, crank, conn->idx_size, (uint64_t)mem->buf, offset );
        /* send acknowledgement of failure to client */
        ardma_server_send_ack1(asc, idx, 1);
        conn->needs_ack = 0;
        return 1;
    }

    retval = ardma_server_post_send ((uint64_t)mem->buf+offset, conn->idx_size, mr->lkey,
                                     conn->idx_addr, conn->idx_rkey, conn);
    if (!retval) {
        conn->waitingfor = COMPL_PULL_INDEX;
    } else { /* send acknowledgement of failure to client */
        ardma_server_send_ack1(asc, idx, 1);
        conn->needs_ack = 0;
    }

    return retval;
}

/* Pull the next block of the PG data block of a client. 
   Offsets are saved in conn struct itself and are incremented
   in this function.
*/
static int ardma_server_pull_block (struct ardma_server_connection *asc,
                                    struct connection * conn)
{
    int idx = conn->crank - asc->lrank;  // 0..asc->nc-1
    int retval;
    uint64_t pullsize = conn->pg_size - conn->pg_offset;
    if (pullsize > blocksize) 
        pullsize = blocksize;

    log_debug("rank %d: pull block: client=%d size=%lld raddr=%lld laddr=%lld\n", 
              rank, idx, pullsize, conn->pg_addr+conn->pg_offset, 
              conn->pg_local_addr+conn->pg_offset);

    retval = ardma_server_post_send (conn->pg_local_addr+conn->pg_offset, pullsize, conn->lkey,
                                     conn->pg_addr+conn->pg_offset, conn->pg_rkey, conn);
    if (!retval) {
        conn->waitingfor = COMPL_PULL_DATA;
        conn->pg_offset += pullsize;
        conn->nblocks++;
    } else { /* send acknowledgement of failure to client */
        ardma_server_send_ack1(asc, idx, 1);
        conn->needs_ack = 0;
    }
    return retval;
}

/* Pull the data block from client (unscheduled)*/
int ardma_server_pull_data (struct ardma_server_connection *asc, int crank, 
                            struct ardma_memory *mem, uint64_t offset, uint64_t size)
{
    int idx = crank - asc->lrank;  // 0..asc->nc-1
    struct connection * conn = &conns[idx];
    int retval;
 
    log_debug("rank %d: pull data from client %d local rank=%d pgsize=%lld "
              "raddr=%lld laddr=%lld\n", 
              rank, crank, idx, conn->pg_size, conn->pg_addr, (uint64_t)mem->buf+offset );

    /* Sanity check */
    if (size != conn->pg_size) {
        log_error("rank %d: RDMA Error: Requested size for pulling data != the known size "
                  "of the data block in client (rank %d) (%lld != %lld).\n"
                  "Behavior is undefined from this point\n",
                  rank, crank, size, conn->pg_size);
    }

    if (offset+conn->pg_size > mem->size) {
        log_error("rank %d: RDMA Error: data does not fit into server's memory, "
                  "client %d pg size=%lld\n mem size=%lld offset=%lld\n", 
              rank, crank, conn->pg_size, (uint64_t)mem->buf, offset );
        /* send acknowledgement of failure to client */
        ardma_server_send_ack1(asc, idx, 1);
        conn->needs_ack = 0;
        return 1;
    }

    /* Pull only first block if buffer is big */
    conn->pg_offset = 0;
    conn->nblocks = 0;
    conn->lkey  = ((struct ibv_mr *)mem->rdma_data)->lkey;
    conn->pg_local_addr = (uint64_t)mem->buf+offset;  

    retval = ardma_server_pull_block (asc, conn);

    return retval;
}

/* Finalize staging server */
int ardma_server_finalize(struct ardma_server_connection ** asc_p)
{
    struct ardma_server_connection * asc = *asc_p;
    int i;

    if (conns) {
        // destroy each QP per client if they had not disconnected (should have)
        for (i==0; i<asc->nc; i++) {
            if (conns[i].id != 0) {
                log_debug("rank %d: transport: destroy rdma qp & id client=%d\n",
                          rank, conns[i].crank);

                rdma_destroy_qp (conns[i].id);
                rdma_destroy_id (conns[i].id);
            }
        }
        free(conns);
        conns = NULL;
    }

    rdma_destroy_id(listener);
    rdma_destroy_event_channel(ec);

    if (asc) free(asc);
    asc = NULL;

    if (sendbuf) free (sendbuf);
    sendbuf = NULL;

    if (sendbuf_mr) ibv_dereg_mr (sendbuf_mr);
    sendbuf_mr = NULL;
    
    if (recvbuf) free (recvbuf);
    recvbuf = NULL;

    if (recvbuf_mr) ibv_dereg_mr (recvbuf_mr);
    recvbuf_mr = NULL;

    return 0;
}

/***********************************************************************/
/*                          INTERNAL FUNCTIONS                         */
/***********************************************************************/

int post_receive (struct ibv_qp *qp, uint64_t wrid, uint64_t addr, 
                  uint32_t size, uint32_t lkey, int rank, int crank)
{
    /* Pre-post a receive for this client */
    struct ibv_sge list = {
        .addr   = addr, 
        .length = size,
        .lkey   = lkey
    };
    struct ibv_recv_wr wr = {
        .wr_id      = wrid,
        .sg_list    = &list,
        .num_sge    = 1,
    };
    struct ibv_recv_wr *bad_wr;

    if (ibv_post_recv(qp, &wr, &bad_wr)) {
            log_error("rank %d: RDMA Error: Couldn't post receive on the message "
                      "buffer for client %d: %s\n",
                      rank, crank, strerror(errno));
            return 1;
    }
    return 0;
}

int on_connect_request(struct rdma_cm_id *id, 
                       const void * private_data, uint8_t private_data_len,
                       struct ardma_server_connection *asc)
{
    struct rdma_conn_param cm_params;
    struct ardma_connect_data data;
    int crank;

    log_debug("rank %d: received connection request, id=%x\n", rank, id);

    // get data that client has sent in connection request
    memcpy(&data, private_data, sizeof(struct ardma_connect_data) );
    crank = ntohl(data.rank);
    asc->nc = ntohl(data.nc);
    asc->appsize  = ntohl(data.size);
    asc->lrank  = ntohl(data.lrank);
    log_debug("rank %d: client %d sent nclients=%d appsize=%d lrank=%d\n", 
              rank, crank, asc->nc, asc->appsize, asc->lrank);

    /* Create verbs objects now that we know which device to use */
    if (!pd) {
        errno = 0;
        // at very first connection only...
        pd = ibv_alloc_pd(id->verbs);
        if (!pd) {
            log_error("rank %d: RDMA Error: Couldn't create PD (protection domain): %s\n",
                    rank, strerror(errno));
            return 1;
        }

        comp_channel = ibv_create_comp_channel(id->verbs);
        if (!comp_channel) {
            log_error("rank %d: RDMA Error: Couldn't create completion channel: %s\n",
                    rank, strerror(errno));
            return 1;
        }

        /* change to non-blocking mode for the completion channel */
        int flags = fcntl(comp_channel->fd, F_GETFL);
        int rc = fcntl(comp_channel->fd, F_SETFL, flags | O_NONBLOCK);
        if (rc < 0) {
            log_error ("rank %d: RDMA Error: Failed to change file descriptor of "
                    "completion event channel to non-blocking mode\n", rank);
            return 1;
        }

        cq = ibv_create_cq(id->verbs, asc->nc+1, NULL, comp_channel, 0);
        if (!cq) {
            log_error("rank %d: RDMA Error: Couldn't create CQ (completion queue): %s\n",
                    rank, strerror(errno));
            return 1;
        }

        if (ibv_req_notify_cq(cq, 0))
            return 1;

        // create connections array
        log_debug("rank %d: Allocate %d bytes for connection variables for %d clients\n", 
                   rank, asc->nc * sizeof(struct connection), asc->nc);
        conns  = (struct connection *)malloc (asc->nc*sizeof(struct connection));
        if (!conns) {
            log_error("rank %d: RDMA Error: Couldn't allocate connections array with size %d\n",
                      rank, asc->nc * sizeof(struct connection));
            return 1;
        }

        // allocate the receive/send message buffers for 'asc->nc' messages now
        log_debug("rank %d: Allocate %d bytes for messages from %d clients\n", 
                   rank, 2*asc->nc * sizeof(struct ardma_message), asc->nc);
        recvbuf = (struct ardma_message *) malloc (asc->nc * sizeof(struct ardma_message));
        if (!recvbuf) {
            log_error("rank %d: RDMA Error: Couldn't allocate recv message buffer with size %d\n",
                      rank, asc->nc * sizeof(struct ardma_message));
            return 1;
        }

        sendbuf = (struct ardma_message *) malloc (asc->nc * sizeof(struct ardma_message));
        if (!sendbuf) {
            log_error("rank %d: RDMA Error: Couldn't allocate send message buffer with size %d\n",
                      rank, asc->nc * sizeof(struct ardma_message));
            return 1;
        }

        // register memory
        recvbuf_mr = ibv_reg_mr(pd, recvbuf, asc->nc * sizeof(struct ardma_message), IBV_ACCESS_LOCAL_WRITE);
        if (!recvbuf_mr) {
            log_error("rank %d: RDMA Error: Couldn't register recv message buffer with size %d: %s\n",
                      rank, asc->nc * sizeof(struct ardma_message), strerror(errno));
            return 1;
        }

        sendbuf_mr = ibv_reg_mr(pd, sendbuf, asc->nc * sizeof(struct ardma_message), IBV_ACCESS_LOCAL_WRITE);
        if (!sendbuf_mr) {
            log_error("rank %d: RDMA Error: Couldn't register send message buffer with size %d: %s\n",
                      rank, asc->nc * sizeof(struct ardma_message), strerror(errno));
            return 1;
        }

    }

    /* Create queue pair */
    struct ibv_qp_init_attr qp_attr;
    memset(&qp_attr, 0, sizeof(qp_attr));

    qp_attr.send_cq = cq;
    qp_attr.recv_cq = cq;
    qp_attr.qp_type = IBV_QPT_RC;

    qp_attr.cap.max_send_wr = 10;
    qp_attr.cap.max_recv_wr = 10;
    qp_attr.cap.max_send_sge = 1;
    qp_attr.cap.max_recv_sge = 1;

    if (rdma_create_qp(id, pd, &qp_attr))
        return 1;

    /* Modify queue pair */
    struct ibv_qp_attr qpa;
    qpa.retry_cnt = 7; 
    qpa.timeout = 28;
    ibv_modify_qp (id->qp, &qpa, IBV_QP_RETRY_CNT | IBV_QP_TIMEOUT);    


    /* Create connection struct for this client */
    int idx = crank - asc->lrank;  // 0..asc->nc-1
    struct connection * conn = &conns[idx];
    memset (conn, 0, sizeof(struct connection));
    id->context = conn;
    conn->id = id;
    conn->crank = crank;
    conn->connected = 0; // will be connected in on_connection()
    conn->waitingfor = COMPL_READY;
 
    log_debug("rank %d: client %d idx=%d conn=%x id=%x\n", 
              rank, crank, idx, conn, conn->id );

    /* Pre-post a receive for this client */
#if 1 
    int poststat = post_receive (id->qp, idx, 
                        (uint64_t) &recvbuf[idx], // nth slot will contain this clients message
                        sizeof(struct ardma_message), 
                        recvbuf_mr->lkey, rank, crank);

    if (poststat) return poststat;
#else
    struct ibv_sge list = {
        .addr   = (uintptr_t) &recvbuf[idx], // nth slot will contain this clients message
        .length = sizeof(struct ardma_message),
        .lkey   = recvbuf_mr->lkey
    };
    struct ibv_recv_wr wr = {
        .wr_id      = idx,
        .sg_list    = &list,
        .num_sge    = 1,
    };
    struct ibv_recv_wr *bad_wr;

    if (ibv_post_recv(id->qp, &wr, &bad_wr)) {
            log_error("rank %d: RDMA Error: Couldn't post receive on the message buffer for client %d: %s\n",
                      rank, crank, strerror(errno));
            return 1;
    }
#endif

    memset (&cm_params, 0, sizeof(struct rdma_conn_param));
    cm_params.responder_resources = 1;
    cm_params.initiator_depth     = 1;
    cm_params.rnr_retry_count     = 7; 
    cm_params.private_data        = NULL;
    cm_params.private_data_len    = 0;

    if (rdma_accept(id, &cm_params)) {
            log_error("rank %d: RDMA Error: Accepting a connection failed: %s\n",
                    rank, strerror(errno));
        return 1;
    }

    return 0;
}

int on_connection(struct rdma_cm_id *id, struct ardma_server_connection *asc)
{
    log_debug("rank %d: a connection has been established.\n", rank);
    if (asc->nc_current == -1)
        asc->nc_current = 1; // very first client
    else
        asc->nc_current++;

    ((struct connection *)id->context)->connected = 1;

    return 0;
}

int on_disconnect(struct rdma_cm_id *id, struct ardma_server_connection *asc)
{
    log_debug("rank %d: client disconnected from server.\n", rank);

    struct connection *conn = (struct connection *)id->context;
    rdma_destroy_qp(conn->id);
    rdma_destroy_id(conn->id);
    conn->id = NULL;

    asc->nc_current--;
    return 0;
}

int on_cm_event(struct rdma_cm_event *event, struct ardma_server_connection *asc)
{
    int r = 1; // error by default

    switch (event->event) {
        case RDMA_CM_EVENT_CONNECT_REQUEST:
            r = on_connect_request (event->id, event->param.conn.private_data, 
                                    event->param.conn.private_data_len, asc);
            break;
        case RDMA_CM_EVENT_ESTABLISHED:
            r = on_connection (event->id, asc);
            break;
        case RDMA_CM_EVENT_DISCONNECTED:
            r = on_disconnect (event->id, asc);
            break;
        default: 
            log_error("rank %d: unexpected RDMA event received: %s.\n",
                      rank, rdma_event_str(event->event));
            break;
    }

    return r;
}

