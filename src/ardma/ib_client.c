/* 
   RDMA Client library for Infiniband
   Implements the functions in ardma_client.h 
*/

#include "ardma_client.h" 
#include "ardma_common.h"
#include "ib_common.h"

#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <arpa/inet.h>
#include <poll.h>

//#include <infiniband/arch.h>  // htonll included by ib_common.h
#include <rdma/rdma_cma.h>


FILE *ardma_logf;
int ardma_verbose_level;
static int connection_established = 0; // 1: when connection to server is made

enum REQUEST_STATE {
    STATE_INIT,              // no outstanding request (e.g. server sent completion message)
    STATE_REQUEST_SENDING,   // client has posted send request for pulling its data
    STATE_REQUEST_SENT,      // ibv_post_send's completion event has been processed
                             //   it does NOT mean the data has been pulled
    STATE_FAILED             // previous request for pull failed (the pull failed)
};

/* Handler for a connection between a client and a server */
struct ardma_client_connection {
    MPI_Comm comm;
    int rank;
    int size;

    // RDMA Connection Manager variables
    char *server_addr;
    char *server_port;
    struct rdma_cm_id *cmid;  // cmid->context will point to the container struct
    struct rdma_event_channel *ec;

    // for cascaded connection to the server
    int lrank; // the lowest rank which connects to the same server process
    int hrank; // the highest rank which connects to the same server process

    // for sending messages
    struct ibv_pd *pd;            // protection domain
    struct ibv_cq *cq;            // completion queue
    struct ibv_comp_channel *comp_channel;
    //struct ibv_qp *qp;            // queue pair to server
    struct ibv_mr *send_mr;        // message buffer region to send message to server
    struct ardma_message *sendbuf; // memory buffer for the send_mr
    struct ibv_mr *recv_mr;        // message buffer region to get message from server
    struct ardma_message *recvbuf; // memory buffer for the recv_mr
    struct ibv_mr *pg_mr;         // memory region that server can pull 
    struct ibv_mr *idx_mr;        // buffer that server can pull 
    void * pg_buf;                // memory buffer for pg_mr
    void * idx_buf;               // memory buffer for idx_mr

    enum REQUEST_STATE state;     // state of request sent/completed cycle
};

static int on_addr_resolved(struct rdma_cm_id *id, struct ardma_client_connection *acc);
static int on_connection(struct rdma_cm_id *id, struct ardma_client_connection *acc);
static int on_disconnect(struct rdma_cm_id *id, struct ardma_client_connection *acc);
static int on_event(struct rdma_cm_event *event, struct ardma_client_connection *acc);
static int on_route_resolved(struct rdma_cm_id *id, struct ardma_client_connection *acc);
static void calc_index_boundary (int M, int N, int R, int * sip, int * qp, int * sp);


#include <ifaddrs.h>
static void print_myhost(int rank) {
    /* 1. get my Infiniband IP address */
    struct ifaddrs *ifaddr, *ifa;
    int s;
    char host[NI_MAXHOST], *ibip = NULL;

    if (getifaddrs(&ifaddr) == -1) { 
        log_error("%s Error: getifaddrs() failed: %s\n", __func__, strerror(errno));
        return;
    }

    for (ifa = ifaddr; ifa != NULL; ifa = ifa->ifa_next) {
        if (ifa->ifa_addr->sa_family == AF_INET) {
            s = getnameinfo(ifa->ifa_addr, sizeof(struct sockaddr_in),
                    host, NI_MAXHOST, NULL, 0, NI_NUMERICHOST);
            if (s != 0) {
                log_warn("getnameinfo() failed: %s\n", gai_strerror(s));
            }
            //log_debug("<Interface>: %s \t <Address> %s\n", ifa->ifa_name, host);
            if (!strncmp(ifa->ifa_name,"ib",2)) {
                ibip = host;
                break;
            }
        }
    }
    log_debug("rank %d: host = %s\n", rank, ibip);
}

/* Initialize RDMA */
int ardma_client_init (int verbose_level)
{
    ardma_verbose_level = verbose_level;
    ardma_logger_init(NULL, verbose_level, -1);
    return 0;
}

/* Make a connection to a server process 
   The mpi rank and size is needed to determine to which 
   server to connect to. 
   Strategy: if there are M server processes and N app processes
   then rank R will connect to server with rank (int) ((float)M/N) * R
   Connection info is in a file written by the server
   Format: 
      int number_of_servers\n
      string IP:port\n
      string IP:port\n
      ... repeated number_of_servers times
*/
struct ardma_client_connection * ardma_client_connect (MPI_Comm comm)
{
    /* connection info is in a file created by the server. */
    struct ardma_client_connection * acc = 
          (struct ardma_client_connection *) malloc (sizeof(struct ardma_client_connection));
    memset (acc, 0, sizeof(struct ardma_client_connection));
#define CLEN 128
    char *conns = 0;
    int server_size = 0;

    connection_established = 0;

    MPI_Comm_dup(comm, &acc->comm);
    MPI_Comm_rank (acc->comm, &acc->rank);
    MPI_Comm_size (acc->comm, &acc->size);

    print_myhost(acc->rank);
    if (!acc->rank) {
        // root process reads in the connection file and then sends 
        // everyone its connection info
        int i, n;
        FILE * f;
        char line[256];
        f = fopen(ARDMA_CONNECTION_FILENAME,"r");
        if (f) {
            if (fgets (line, sizeof(line), f)) {
                server_size = atoi (line);
                log_debug("number of servers based on first line in file = %d.\n", server_size);         
                if (server_size > 0 && server_size < 10000) { // make it sane

                    conns = (char *) malloc (server_size * CLEN * sizeof(char));
                    memset (conns, 0, server_size * CLEN * sizeof(char));
                
                    n = 0;
                    while ( n < server_size && fgets(line, sizeof(line), f)) {
                        strncpy (conns+n*CLEN, line, CLEN);
                        log_debug("server %d listens on %s", n, conns+n*CLEN);         
                        n++;
                    }
                    if (n != server_size) {
                        log_error("%s Error: first line in file %s says there are %d servers "
                                "but we found %d connection lines. Will use the latter.\n",
                                __func__, ARDMA_CONNECTION_FILENAME, server_size, n);
                        server_size = n;
                    }
                } else {
                    log_error("%s Error: first line in file %s says there are %d servers.\n",
                      __func__, ARDMA_CONNECTION_FILENAME, server_size);
                      server_size = 0;
                }
            } else {
                log_error("%s Error: cannot read first line from file %s.\n",
                   __func__, ARDMA_CONNECTION_FILENAME);
            }
            fclose(f);
        } else {
            log_error("%s Error: cannot open file %s. %s\n",
                      __func__, ARDMA_CONNECTION_FILENAME, strerror(errno));
        }

        // send info to everyone
        MPI_Bcast (&server_size, 1, MPI_INT, 0, acc->comm); 
        if (server_size)
            MPI_Bcast (conns, server_size*CLEN, MPI_CHAR, 0, acc->comm); 

    } else { 
        // other processes get from root how many servers we have
        MPI_Bcast (&server_size, 1, MPI_INT, 0, acc->comm); 
        if (server_size) {
            conns = (char *) malloc (server_size * CLEN * sizeof(char));
            MPI_Bcast (conns, server_size*CLEN, MPI_CHAR, 0, acc->comm); 
        }
    }

    if (server_size) {
        // Calculate the server index to whom this process connects plus
        //   the lowest and highest rank which connects to the same server.
        // Rank R will connect to server with rank (int) ((float)M/N) * R
        int idx; // = ((float)server_size/acc->size) * acc->rank;
        calc_index_boundary (server_size, acc->size, acc->rank, &idx, &acc->lrank, &acc->hrank);
        char s_addr[CLEN];
        char s_port[CLEN];
        sscanf (conns+idx*CLEN, "%[^:]:%s", s_addr, s_port);
        log_info("rank %d will connect to server %d at %s:%s\n", acc->rank, idx, s_addr, s_port);
        acc->server_addr = strdup(s_addr);
        acc->server_port = strdup(s_port);
    } else {
        free (acc);
        acc = NULL;
    }

    if (conns) free(conns);

    // return back now if failed above
    if (!acc) return NULL;

    // Now try to connect to RDMA server
    struct addrinfo *addr;
    struct rdma_cm_event *event = NULL;

#define RETURN_NOW(x) log_error("%s, rank %d Error: " #x " failed\n", __func__, acc->rank); \
                      free(acc); \
                      return NULL;

#define TEST_NZ(x) if ( (x)) {RETURN_NOW(x)};
#define TEST_Z(x)  if (!(x)) {RETURN_NOW(x)};

    TEST_NZ(getaddrinfo(acc->server_addr, acc->server_port, NULL, &addr));

    TEST_Z(acc->ec = rdma_create_event_channel());
    TEST_NZ(rdma_create_id(acc->ec, &(acc->cmid), (void *)acc, RDMA_PS_TCP));
    TEST_NZ(rdma_resolve_addr(acc->cmid, NULL, addr->ai_addr, TIMEOUT_IN_MS));

    freeaddrinfo(addr);

    while ( !connection_established && 
            rdma_get_cm_event(acc->ec, &event) == 0) {

        struct rdma_cm_event event_copy;

        memcpy (&event_copy, event, sizeof(*event));
        rdma_ack_cm_event (event);

        if (on_event(&event_copy, acc))
            break; // exit on non-zero return value
    }

    if (connection_established) {
        acc->state = STATE_INIT; //ready to send request to server
    } else {
        // failed to get connected
        if (acc->cmid->qp) {
            rdma_destroy_qp (acc->cmid);
            acc->cmid->qp = NULL;
        }
        rdma_destroy_id (acc->cmid);
        rdma_destroy_event_channel (acc->ec);
        free(acc);
        acc = NULL;
    } 

    return acc;
}


/* Disconnect from staging server */
int ardma_client_disconnect (struct ardma_client_connection **acc_p)
{
    struct ardma_client_connection * acc = *acc_p;
    /* FIXME: send disconnect message to server ? or will the destroy create an event there?*/
    log_debug("finalize acc %x\n", acc);
    if (acc) {
        log_info("rank %d will disconnect from server %s:%s\n", 
                 acc->rank, acc->server_addr, acc->server_port);
        if (acc->server_addr)
            free(acc->server_addr);
        if (acc->server_port)
            free(acc->server_port);
        if (connection_established) {
            rdma_destroy_qp (acc->cmid);
            rdma_destroy_id (acc->cmid);
            rdma_destroy_event_channel (acc->ec);
            connection_established = 0;
        }
        free(acc);
        *acc_p = NULL;
    }
    return 0;
}


void * ardma_client_register_pg_memory (struct ardma_client_connection *acc,
                                        uint64_t size, void * addr)
{
    if (!acc)
        return NULL;

    if (addr) {
        if (addr != acc->pg_buf) {
            log_error("rank %d: Code logic error: called %s with a non-NULL pointer, "
                      "which is not the PG buffer.\n", acc->rank, __func__);
        }
        ibv_dereg_mr (acc->pg_mr);
        free (addr);
    }

    acc->pg_buf = alloc_memory_aligned(size); 
    if (!acc->pg_buf) {
        log_error("rank %d: RDMA Error: Failed to allocate %lld bytes of memory\n", acc->rank, size);
        return NULL;
    }

    acc->pg_mr = ibv_reg_mr(acc->pd, acc->pg_buf, size, 
                                IBV_ACCESS_LOCAL_WRITE |
                                IBV_ACCESS_REMOTE_READ);

    if (!acc->pg_mr) {
        log_error("rank %d: RDMA Error: Failed to register %lld bytes of memory to RDMA\n", acc->rank, size);
        free (acc->pg_buf);
        return NULL;
    }
    log_debug("rank %d: --- pg_buf=%x, pg_mr=%x lkey=%ld\n", acc->rank, acc->pg_buf, acc->pg_mr, acc->pg_mr->lkey);
    return acc->pg_buf;
}

void * ardma_client_register_index_memory (struct ardma_client_connection *acc,
                                         uint64_t size, void * addr)
{
    if (!acc)
        return NULL;

    if (addr) {
        if (addr != acc->idx_buf) {
            log_error("rank %d: Code logic error: called %s with a non-NULL pointer, "
                      "which is not the index buffer.\n", acc->rank, __func__);
        }
        ibv_dereg_mr (acc->idx_mr);
        free (addr);
    }

    acc->idx_buf = alloc_memory_aligned(size); 
    if (!acc->idx_buf) {
        log_error("rank %d: RDMA Error: Failed to allocate %lld bytes of memory\n", acc->rank, size);
        return NULL;
    }

    acc->idx_mr = ibv_reg_mr(acc->pd, acc->idx_buf, size, 
                                IBV_ACCESS_LOCAL_WRITE |
                                IBV_ACCESS_REMOTE_READ);

    if (!acc->idx_mr) {
        log_error("rank %d: RDMA Error: Failed to register %lld bytes of memory to RDMA\n", acc->rank, size);
        free (acc->idx_buf);
        return NULL;
    }
    log_debug("rank %d: --- idx_buf=%x, idx_mr=%x, lkey=%ld\n", acc->rank, acc->idx_buf, acc->idx_mr, acc->idx_mr->lkey);
    return acc->idx_buf;
}

int ardma_client_deregister_pg_memory(struct ardma_client_connection *acc, void * addr)
{
    if (!acc)
        return -1;
    if (addr != acc->pg_buf) {
        log_error("rank %d: Code logic error: called %s with a pointer, "
                "which is not the PG buffer.\n", acc->rank, __func__);
    }
    ibv_dereg_mr (acc->pg_mr);
    free (acc->pg_buf);
    return 0;
}

int ardma_client_deregister_index_memory(struct ardma_client_connection *acc, void * addr)
{
    if (!acc)
        return -1;
    if (addr != acc->idx_buf) {
        log_error("rank %d: Code logic error: called %s with a pointer, "
                "which is not the index buffer.\n", acc->rank, __func__);
    }
    ibv_dereg_mr (acc->idx_mr);
    free (acc->idx_buf);
    return 0;
}


//Send the request to staging server to pull the completed buffers
int ardma_client_send_request(struct ardma_client_connection *acc, 
                              uint64_t pg_size, 
                              uint64_t idx_size,
                              char *path,
                              int timestep)
{
    /* First, post a receive for the server's completion message,
       which is expected much later but make sure here we are ready to receive */
    struct ibv_sge list = {
        .addr   = (uintptr_t) acc->recvbuf,
        .length = sizeof(struct ardma_message),
        .lkey   = acc->recv_mr->lkey
    };

    struct ibv_recv_wr rwr = {
        .wr_id      = 1, 
        .sg_list    = &list,
        .num_sge    = 1,
    };
    struct ibv_recv_wr *bad_rwr;

    if (ibv_post_recv(acc->cmid->qp, &rwr, &bad_rwr)) {
        log_error("rank %d: RDMA Error: Couldn't post receive on the message buffer: %s\n",
                acc->rank, strerror(errno));
        return 1;
    }

    uint64_t cksum = ardma_calc_checksum (acc->pg_buf, pg_size);
    //uint64_t cksum2 = ardma_calc_checksum (acc->pg_buf, pg_size);

    log_debug("rank %d: Send request:\n"
            "  PG  (size=%lld, addr=%lld, key=%u)\n"
            "  IDX (size=%lld, addr=%lld, key=%u)\n"
            "  timestep = %u\n"
            "  path     = %s\n"
            "  checksum = %llu\n",
            acc->rank, 
            pg_size, acc->pg_buf, acc->pg_mr->lkey, 
            idx_size, acc->idx_buf, acc->idx_mr->lkey, 
            timestep, path, cksum);

    /*
    int i;
    log_debug("rank %d: PG = [", acc->rank);
    for (i=0; i<8; i++) 
        log_debug("%hhu ",(char) *((char*)(acc->pg_buf+i)));
    log_debug("... ");
    for (i=8; i>0; i--)
        log_debug("%hhu ",(char) *((char*)acc->pg_buf+pg_size-i));
    log_debug("]\n");
    */

    /* prepare message for send */
    struct ardma_request_message req = {
        .rank     = htonl(acc->rank),
        .pg_size  = htonll(pg_size),
        .pg_addr  = htonll((uintptr_t) acc->pg_buf),
        .pg_rkey  = htonl(acc->pg_mr->lkey),
        .pg_cksum = htonll(cksum),
        .idx_size = htonll(idx_size), 
        .idx_addr = htonll((uintptr_t) acc->idx_buf),
        .idx_rkey = htonl(acc->idx_mr->lkey),
        .timestep = htonl(timestep),
    };
    strncpy (req.path, path, MAXPATHLEN);

    acc->sendbuf->type = ARDMA_MSG_REQUEST;
    memcpy (&acc->sendbuf->msg.request, &req, sizeof(struct ardma_request_message));

    list.addr   = (uintptr_t) acc->sendbuf;
    list.length = sizeof(struct ardma_message);
    list.lkey   = acc->send_mr->lkey;

    struct ibv_send_wr wr = {
        .wr_id      = 0,
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
    if (ibv_post_send(acc->cmid->qp, &wr, &bad_wr)) {
        log_error("rank %d: RDMA Error: post_send failed: (errno=%d) %s\n", 
            acc->rank, errno, strerror(errno));
        log_debug("wr=%x, bad_wr=%x\n",&wr,bad_wr);
        return 1;
    }

    acc->state = STATE_REQUEST_SENDING; //ready to send request to server

    return 0;
}


/* check if previous staging has completed */
enum ARDMA_STAGING_STATUS ardma_client_staging_status (struct ardma_client_connection *acc)
{
    if (!acc)
        return ARDMA_STAGING_FAILED;

    /* Check if completion event for the last send exists and
       if server's acknowledgement message has arrived */

    /* poll the channel for an event without blocking */
    struct pollfd p = {
        .fd      = acc->comp_channel->fd,
        .events  = POLLIN,
        .revents = 0 
    };

    int ms_timeout = 0; // do not wait at all, since event should be there already

    int rc;
    do {
        rc = poll(&p, 1, ms_timeout);
        if (rc < 0) {
            log_error ("rank %d: RDMA Error: polling completion channel failed. errno = %d : %s\n",
                    acc->rank, errno, strerror(errno));
            acc->state = STATE_FAILED;
            continue; // exit lop
        }
        else if (rc == 0) continue; // no event, exit loop

        // get the event(s) out from the completion queue
        void * ev_ctx;
        struct ibv_cq *ev_cq;
        if (ibv_get_cq_event(acc->comp_channel, &ev_cq, &ev_ctx)) {
            log_error("rank %d: RDMA Error: Failed to get cq_event\n");
            acc->state = STATE_FAILED;
            continue; 
        }

        if (ibv_req_notify_cq(acc->cq, 0)) {
            acc->state = STATE_FAILED;
            continue; 
        }

        ibv_ack_cq_events(ev_cq, 1);

        /* Empty the CQ: poll all of the completions from the CQ (if any exist) */
        int ne;
        do {
            struct ibv_wc wc;
            ne = ibv_poll_cq(acc->cq, 1, &wc);
            if (ne < 0) {
                log_error("rank %d: RDMA Error: Failed to poll completions from the CQ\n", acc->rank);
                acc->state = STATE_FAILED;
                ne = rc = 0; // exit from both loops
                continue;
            }

            /* there may be an extra event with no completion in the CQ */
            if (ne == 0)
                continue;

            if (wc.opcode == IBV_WC_SEND) {
                // client's send completed
                acc->state = STATE_REQUEST_SENT;
                log_debug("rank %d: Request send completed\n", acc->rank);
            } else if (wc.opcode == IBV_WC_RECV) {
                // server sent message about completed pull
                acc->state = STATE_INIT;
                log_debug("rank %d: Server acknowledged pulling data\n", acc->rank);
            } else {
                log_error("rank %d: RDMA Error: Unexpected completion event opcode=%d\n", 
                            acc->rank, wc.opcode);
                acc->state = STATE_FAILED;
            }
                
            if (wc.status != IBV_WC_SUCCESS) {
                log_error("rank %d: RDMA Error: Completion event failed: status = %d: %s\n", 
                            acc->rank, wc.status, ibv_wc_status_str(wc.status));
                acc->state = STATE_FAILED;
            }

        } while (ne);

    } while (rc > 0);


    enum ARDMA_STAGING_STATUS s = ARDMA_STAGING_READY;
    switch (acc->state) {
        case STATE_INIT:
            s = ARDMA_STAGING_READY;
            log_debug ("rank %d: Previous send has completed.\n", acc->rank);
            break;
        case STATE_REQUEST_SENDING:
        case STATE_REQUEST_SENT:
            s = ARDMA_STAGING_INPROGRESS;
            log_debug ("rank %d: Previous send has not completed yet.\n", acc->rank);
            break;
        case STATE_FAILED:
            s = ARDMA_STAGING_FAILED;
            log_debug ("rank %d: Previous staging failed.\n", acc->rank);
            break;
    }

    return s; 
}



/***********************************************************************/
/*                          INTERNAL FUNCTIONS                         */
/***********************************************************************/
int on_addr_resolved(struct rdma_cm_id *id, struct ardma_client_connection *acc)
{
    log_debug("rank %d: address resolved. resolve route...\n", acc->rank);

    /* Create verbs objects now that we know which device to use */
    acc->pd = ibv_alloc_pd(id->verbs);
    if (!acc->pd)
        return 1;

    acc->comp_channel = ibv_create_comp_channel(id->verbs);
    if (!acc->comp_channel)
        return 1;

    /* change to non-blocking mode for the completion channel */
    int flags = fcntl(acc->comp_channel->fd, F_GETFL);
    int rc = fcntl(acc->comp_channel->fd, F_SETFL, flags | O_NONBLOCK);
    if (rc < 0) {
       log_error ("rank %d: RDMA Error: Failed to change file descriptor of "
                  "completion event channel to non-blocking mode\n", acc->rank);
        return 1;
    }


    acc->cq = ibv_create_cq(id->verbs, 2, NULL, acc->comp_channel, 0);
    if (!acc->cq)
        return 1;

    if (ibv_req_notify_cq(acc->cq, 0))
        return 1;

    /* allocate small message memory to send requests and register for rdma */
    acc->sendbuf = (struct ardma_message *) malloc (sizeof(struct ardma_message));
    if (!acc->sendbuf) {
        log_error("rank %d: RDMA Error: Couldn't allocate receive message buffer with size %d\n",
                acc->rank, sizeof(struct ardma_message));
        return 1;
    }

    log_debug("rank %d: --- sendbuf=%x, size=%d\n", acc->rank, acc->sendbuf, sizeof(struct ardma_message));
    acc->send_mr = ibv_reg_mr(acc->pd, acc->sendbuf, 
                                sizeof(struct ardma_message), 
                                IBV_ACCESS_LOCAL_WRITE);
    if (!acc->send_mr) {
        log_error("rank %d: RDMA Error: Couldn't register send message buffer with size %d: %s\n",
                acc->rank, sizeof(struct ardma_message), strerror(errno));
        return 1;
    }

    /* allocate small message memory to receive acknowledgement and register for rdma */
    acc->recvbuf = (struct ardma_message *) malloc (sizeof(struct ardma_message));
    if (!acc->recvbuf) {
        log_error("rank %d: RDMA Error: Couldn't allocate receive message buffer with size %d\n",
                acc->rank, sizeof(struct ardma_message));
        return 1;
    }

    log_debug("rank %d: --- recvbuf=%x, size=%d\n", acc->rank, acc->recvbuf, sizeof(struct ardma_message));
    acc->recv_mr = ibv_reg_mr(acc->pd, acc->recvbuf, 
                                sizeof(struct ardma_message), 
                                IBV_ACCESS_LOCAL_WRITE);
    if (!acc->recv_mr) {
        log_error("rank %d: RDMA Error: Couldn't register receive message buffer with size %d: %s\n",
                acc->rank, sizeof(struct ardma_message), strerror(errno));
        return 1;
    }

    /* Create queue pair */
    struct ibv_qp_init_attr qp_attr;
    memset(&qp_attr, 0, sizeof(qp_attr));

    qp_attr.send_cq = acc->cq;
    qp_attr.recv_cq = acc->cq;
    qp_attr.qp_type = IBV_QPT_RC;

    qp_attr.cap.max_send_wr = 10;
    qp_attr.cap.max_recv_wr = 10;
    qp_attr.cap.max_send_sge = 1;
    qp_attr.cap.max_recv_sge = 1;

    if (rdma_create_qp(id, acc->pd, &qp_attr))
        return 1;

    //sprintf(get_local_message_region(id->context), 
    //        "msg from client with rank %d", acc->rank);

    if (rdma_resolve_route(id, TIMEOUT_IN_MS)) {
        log_error("%s, rank %d Error: Route resolution failed\n", __func__, acc->rank);  
        return 1;
    };

    return 0;
}

int on_route_resolved(struct rdma_cm_id *id, struct ardma_client_connection *acc)
{
    struct rdma_conn_param cm_params;
    struct ardma_connect_data connect_data; // sent with rdma_connect()
    int r=0;
    MPI_Status stat;
    MPI_Request req;

    log_debug("rank %d: route resolved. connect to server...\n", acc->rank);

    memset(&cm_params, 0, sizeof(cm_params));
    cm_params.initiator_depth = cm_params.responder_resources = 1;
    cm_params.rnr_retry_count = 7; /* infinite retry */
    connect_data.rank = htonl(acc->rank);
    connect_data.size = htonl(acc->size);
    connect_data.nc   = htonl(acc->hrank - acc->lrank + 1); 
    connect_data.lrank = htonl(acc->lrank);
    cm_params.private_data = &connect_data;
    cm_params.private_data_len = sizeof(struct ardma_connect_data);
    //log_debug("rank %d: packed private data rank=%d size=%d nc=%d\n", acc->rank,
    //    connect_data.rank, connect_data.size, connect_data.nc);
    log_debug("rank %d: packed private data rank=%d size=%d nc=%d\n", acc->rank,
        ntohl(*(int*)cm_params.private_data), 
        ntohl(*(((int*)cm_params.private_data)+1)), 
        ntohl(*(((int*)cm_params.private_data)+2))
        );

    /* cascaded connection requests to a server */
    if (acc->rank > acc->lrank) { 
        // wait for token from rank-1
        MPI_Recv (&r, 1, MPI_INT, acc->rank-1, acc->rank-1/*TAG*/, acc->comm, &stat);
        /* corresponding send is done in on_connection() or 
           below in case of failed rdma_connect */
    }

    // try to connect only if rank-1 succeeded
    if (!r)
        r = rdma_connect(id, &cm_params); 

    if (r) {
        if (acc->rank < acc->hrank) { 
            // send token (r) to next rank now
            MPI_Isend (&r, 1, MPI_INT, acc->rank+1, acc->rank/*TAG*/, acc->comm, &req);
        }
        log_error("%s, rank %d Error: RDMA connect failed\n", __func__, acc->rank);  
        return 1;
    };

    return 0;
}


int on_connection(struct rdma_cm_id *id, struct ardma_client_connection *acc)
{
    int r = 0;
    MPI_Request req;

    if (acc->rank < acc->hrank) {
        // send token (r) to next rank now
        MPI_Isend (&r, 1, MPI_INT, acc->rank+1, acc->rank/*TAG*/, acc->comm, &req);
    }

    connection_established = 1;
    log_debug("rank %d: connection established.\n", acc->rank);

    return 0;
}

int on_disconnect(struct rdma_cm_id *id, struct ardma_client_connection *acc)
{
    log_debug("rank %d: server disconnected from me.\n", acc->rank);

    rdma_destroy_qp (id);
    rdma_destroy_id (id);
    rdma_destroy_event_channel (acc->ec);
    return 1; /* exit event loop */
}

int on_event(struct rdma_cm_event *event, struct ardma_client_connection *acc)
{
    int r = 1; // 1: bad, 0: okay

    if (event->event == RDMA_CM_EVENT_ADDR_RESOLVED)
        r = on_addr_resolved(event->id, acc);
    else if (event->event == RDMA_CM_EVENT_ROUTE_RESOLVED)
        r = on_route_resolved(event->id, acc);
    else if (event->event == RDMA_CM_EVENT_ESTABLISHED)
        r = on_connection(event->id, acc);
    else if (event->event == RDMA_CM_EVENT_DISCONNECTED)
        r = on_disconnect(event->id, acc);
    else
        log_error("rank %d: unexpected RDMA event received: %s.\n", 
                  acc->rank, rdma_event_str(event->event));

    return r;
}

/* Calculate which lowest (Q) and highest rank (S) has the same server index (SI) to
   connect to than this rank (R).
    Input:
        M = server size
        N = this app's size
        R = rank of this process
    Output:
        sip = server index this process has to connect
        lp  = lowest rank with same server index
        hp  = highest rank with same server index
 */
void calc_index_boundary (int M, int N, int R, int * sip, int * lp, int * hp)
{
    // 1. rank R will connect to server with rank (int) ((float)M/N) * R
    float f = (float)M/N;
    int si = f * R;
    *sip = si;

    // 2. search for lowest rank with same si
    int l = R - N/M - 1; // below this they connect to other server processes
    if (si == 0) {
        l = 0; // for server index 0, rank = 0 is the lowest    
    } else {
        for (; l<R; l++) {
            if ((int)(f*l) == si) break;
        }
    }
    *lp = l;

    // 3. search for highest rank with same si
    int h = R + N/M + 1; // below this they connect to other server processes
    if (si == M-1) {
        h = N-1; // for server index M-1, rank = N-1 is the highest    
    } else {
        for (; h>R; h--) {
            if ((int)(f*h) == si) break;
        }
    }
    *hp = h;
}

