#ifndef __EVENTS_H__
#define __EVENTS_H__

/* Define events passed between Worker and Transport threads */
#include "globals.h" // MAXPATHLEN

enum EVENT_TYPE {
    /* Worker to Transport events (in wtq queue) */
    EVENT_PULL_REQUEST,        // pull data from a client
    EVENT_ACK_REQUEST,         // send acknowledgement to all clients

    /* Transport to Worker events (in twq queue) */
    EVENT_CLIENT_REQUEST,        // client sends request to server to pull data
    EVENT_PULL_COMPLETED,        // client sends request to server to pull data
};


/* I. Worker to Transport events (in wtq queue) */
struct event_pull_request { 
    uint32_t rank;          // client rank
    uint64_t size;          // size of data block 
    uint64_t offset;        // local memory address to fill (offset in registered memory)
                            //   i.e. offset from gd.rdma_buffer which is = mem.buf
};

struct event_ack_request { 
    int status;
};

/* II. Transport to Worker events (in twq queue) */
// Client sent request (1 request to Worker when requests from all clients are here)
struct event_client_request { 
    int      number_of_clients; // number of clients sent the request
    char     path[MAXPATHLEN]; // name of output file
    uint32_t timestep;         // timestep of output
};

struct event_pull_completed { 
    uint32_t rank;          // client rank
};

/* III. Union type of all events */
struct event {
    enum EVENT_TYPE type;
    union {
        struct event_client_request client_request;
        struct event_pull_request   pull_request;
        struct event_pull_completed pull_completed;
        struct event_ack_request    ack_request;
    } event;
};



#endif
