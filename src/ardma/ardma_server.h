/* 
   RDMA Server library

   Server applications should use this header and these functions only to 
   get data from staging methods. 

*/

#ifndef __ADIOS_RDMA_SERVER_H__
#define __ADIOS_RDMA_SERVER_H__

#include <stdint.h> /* uint64_t */
#include <mpi.h>
#include "ardma_logger.h"

/*
// return values for ardma_server_check_events()
enum ARDMA_EVENT { 
    ARDMA_NO_EVENT,              // nothing happened
    ARDMA_CLIENT_CONNECTED,      // a client connected
    ARDMA_CLIENT_DISCONNECTED,   // a client disconnected
    ARDMA_CLIENT_REQUEST,        // a client requested its data to be pulled
    ARDMA_OTHER_EVENT,           // event that you do not need to worry about
    ARDMA_EVENT_HANDLER_FAILURE, // returned if rdma layer below aborted (reason to terminate)
};
*/

struct ardma_server_connection {
    int nc_current;  // current # of clients, -1 if never had any, 0 if there were but now all gone

    /* (each client will send the following values at connect) */
    int nc;          // number of expected clients 
    int appsize;     // total number of clients (to all server processes)
    int lrank;       // lowest rank connecting to this server process 
                     // i.e. rank-lrank will go from 0..nc-1

};

struct ardma_memory {
    char * buf;
    uint64_t size;
    void * rdma_data; // opaque pointer for rdma layer
};


/* Initialize ardma 
   Create listening ports and write a file containing the connection info clients can use
*/
struct ardma_server_connection * ardma_server_init (MPI_Comm comm, int verbose_level);

/* Finalize staging server */
int ardma_server_finalize(struct ardma_server_connection ** asc_p);

/* Allocate and register memory to rdma layer. 
   Return: ardma_memory struct pointer to the allocated buffer, NULL on error
*/
struct ardma_memory * ardma_server_register_memory (uint64_t size);

/* Deregister memory. Return 0 on success, !=0 on error */
int ardma_server_deregister_memory (struct ardma_memory *mem);


/* Check for events (non-blocking). For different events, different 
   callback functions will be called. 
   Return value: number of events
*/
int ardma_server_check_events (struct ardma_server_connection * asc);

/* Pull the index from client into memory pointed by mem->buf at offset 'offset' 
   Return: 0 on success of initiating the pull request, !=0 on error
   See ardma_server_cb_pulled_index() for notification of completion.
*/
int ardma_server_pull_index (struct ardma_server_connection *asc, int crank, 
                             struct ardma_memory *mem, uint64_t offset);

/* Pull the data from client into memory pointed by mem->buf at offset 'offset' 
   Return: 0 on success of initiating the pull request, !=0 on error
   See ardma_server_cb_pulled_data() for notification of completion.
   Note: size is known in the library too, so it is used only for sanity checks.
*/
int ardma_server_pull_data (struct ardma_server_connection *asc, int crank, 
                            struct ardma_memory *mem, uint64_t offset, uint64_t size);

/* Send acknowledgement message to all clients that the staging completed.
   i.e. a new request can be sent.
   status=0 means success, =1 means staging failed.
*/
int ardma_server_send_acknowledgement (struct ardma_server_connection *asc, int status);

/* Callback functions for ardma_server_check_events.
   Your application must define this functions 
*/

// Callback for client connections 
// Very first is interesting for server, when the ARDMA layer is 
// fully initialized so it can register memory at this time
void ardma_server_cb_connect (struct ardma_server_connection * asc,
        int nc, int nc_total, int crank);

void ardma_server_cb_disconnect (struct ardma_server_connection * asc, int crank);

// Callback for client requests (to pull its data)
void ardma_server_cb_request (struct ardma_server_connection * asc, 
                              int crank, int nodeid, 
                              uint64_t pg_size, uint64_t idx_size,
                              int timestep, char * path);

// Callback if a failure happened in ARDMA layer (reason to terminate)
void ardma_server_cb_failure (struct ardma_server_connection * asc);

// Callback when ardma_server_pull_index has accomplished the pull
void ardma_server_cb_pulled_index (struct ardma_server_connection * asc, int crank);

// Callback when ardma_server_pull_data has accomplished the pull
void ardma_server_cb_pulled_data (struct ardma_server_connection * asc, int crank);

/* For testing purposes */
uint64_t ardma_server_checksum (char * buf, uint64_t size);

#endif
