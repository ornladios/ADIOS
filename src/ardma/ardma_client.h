/* 
   RDMA Client library

   Applications should use this header and these functions only to 
   stage their output to a staging server

*/

#ifndef __ADIOS_RDMA_CLIENT_H__
#define __ADIOS_RDMA_CLIENT_H__

#include <stdint.h> /* uint64_t */
#include <mpi.h>
#include "ardma_logger.h"

/* Opaque handler for a connection between a client and a server */
struct ardma_client_connection; 


/* Initialize ardma */
int ardma_client_init (int verbose_level);

/* Make a connection to a server process 
*/
struct ardma_client_connection * ardma_client_connect (MPI_Comm comm);

/* Disconnect from staging server.
   It frees the connection structure, so it must not be used afterwards. */
int ardma_client_disconnect (struct ardma_client_connection ** acc_p);

/* Allocate and register memory to rdma layer. 
   One allocates the large Process Group buffer, 
   the other the small metadata index buffer.
   If the address passed is non-NULL, it is deregistered first and freed 
     before allocating and registering a new buffer. 
   Return: the pointer to the allocated buffer, NULL on error
*/
void * ardma_client_register_pg_memory (struct ardma_client_connection *acc,
                                        uint64_t size, void * addr);
void * ardma_client_register_index_memory (struct ardma_client_connection *acc,
                                         uint64_t size, void * addr);

/* deregister from DMA and free buffer */
int ardma_client_deregister_pg_memory(struct ardma_client_connection *acc, void * addr);
int ardma_client_deregister_index_memory(struct ardma_client_connection *acc, void * addr);

// Send the request to staging server to pull the completed buffers
int ardma_client_send_request (struct ardma_client_connection *acc,
                               uint64_t pg_size,
                               uint64_t idx_size,
                               char *path,
                               int timestep);

enum ARDMA_STAGING_STATUS {
    ARDMA_STAGING_READY,       // init state or previous staging has been completed
    ARDMA_STAGING_INPROGRESS,  // still in progress...
    ARDMA_STAGING_FAILED,      // error occured in previous staging
};

/* Check if previous staging has completed */
enum ARDMA_STAGING_STATUS ardma_client_staging_status (struct ardma_client_connection *acc);

/* For testing purposes */
uint64_t ardma_client_checksum (char * buf, uint64_t size);

#endif
