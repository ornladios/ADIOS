/* Portals server for ARDMA/ADIOS
 */

#include "ardma_server.h"
#include "ardma_common.h"
#include "portals_common.h"

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

enum CONN_STATE 
{
  UNINITIALIZED,
  CONNECTING,
  CONNECTED,
  DISCONNECTED,
  REQUEST_PENDING,
  REQUEST_SERVICING,
  FAILURE
};
	

typedef struct connection_
{
	
  einfo client;
  enum CONN_STATE state;
	
  //anything else we need to know?	
}conn;


typedef struct requests_
{
  req_conn *client;  //link back to the client
  req_data *r; //the actual transfer request from the client
  msg *pg;
  msg *idx;
}requests;



typedef struct sinfo_
{
  unsigned long long maxmem;
  unsigned long long mem;

  einfo self;
  List *conn_list; //list of all the connections
  List *requests; //list of all the requests we have recieved
  List *pending;  //list of all the transfers we have issued
  List *complete; //list of all the completed requests we have to deal with	
}sinfo;

	
ninfo local_info;
MPI_Comm scomm;
sinfo *server;


/* Initialize ardma 
   create connection file
*/

struct ardma_server_connection * ardma_server_init (MPI_Comm comm, int verbose_level)
{
    struct ardma_server_connection * asc = 
	(struct ardma_server_connection*)malloc(sizeof(struct ardma_server_connection));
    int rank = 0, size = 0;
    int retval = 0;
    int *pids, *nids, *pts, *matchs;
	
    ardma_verbose_level = verbose_level;



    portal_init_common(&local_info);
	
    //create the array for the connection information. 
    //this is what we will check in check_event for connections.
    //each connection thats created will get use the 
    //data request event queue and put transfer requests in there
    //disconnect requests also to the connection request queue
    req_conn *conn_list = (req_conn*) malloc(sizeof(req_conn)*LISTSIZE);
    ptl_md_t connmd;
    ptl_handle_md_t connmd_h;
    connmd.start = conn_list;
    connmd.length = sizeof(req_conn)*LISTSIZE;
    connmd.threshold = PTL_MD_THRESH_INF;
    connmd.max_size = sizeof(req_conn)*LISTSIZE;
    connmd.options = PTL_MD_OP_PUT | PTL_MD_TRUNCATE | PTL_MD_MAX_SIZE;
    connmd.eq_handle = local_info.eqh;
    connmd.user_ptr = (void*)asc;
	
    retval = PtlMDAttach(local_info.meh, connmd, 
			 PTL_RETAIN, &connmd_h);
    if(retval != PTL_OK)
    {
	log_error("Portal Error %d %s\n", retval, PtlErrorStr(retval));
	free(asc);
	return NULL;
    }
	
	
	

    /* 3. Collect rdma info from all server process and write into connection file */ 
    scomm = comm;
    MPI_Comm_rank (scomm, &rank);
    MPI_Comm_size (scomm, &size);

    if(rank == 0)
    {
	//create the buffers to hold the connection information
	pids = (int*)malloc(sizeof(int) * size);
	nids = (int*)malloc(sizeof(int) * size);
	pts = (int*) malloc(sizeof(int) * size);
	matchs = (int*) malloc(sizeof(int) * size);
		
    }
	
    MPI_Gather( &local_info.pid.pid, 1, MPI_INT, pids, 1, MPI_INT, 0, scomm);
    MPI_Gather( &local_info.pid.nid, 1, MPI_INT, nids, 1, MPI_INT, 0, scomm);
    MPI_Gather( &local_info.index, 1, MPI_INT, pts, 1, MPI_INT, 0, scomm);
    MPI_Gather( &local_info.match, 1, MPI_INT, matchs, 1, MPI_INT, 0, scomm);

    if(rank == 0)
    {
	int i;
	FILE *f;
	f = fopen(ARDMA_CONNECTION_FILENAME, "w+");
	if(f)
	{
	    for(i = 0; i < size; i++)
	    {
		fprintf(f, "%d %d %d %d\n", 
			pids[i], nids[i], 
			pts[i], matchs[i]);
		log_debug("%d\t %d %d %d %d\n",
			  i, pids[i], nids[i], 
			  pts[i], matchs[i]);
 
	    }
			
	    fflush(f);
	    fclose(f);
			
	    free(pids);
	    free(nids);
	    free(pts);
	    free(matchs);			
	}
		
    }
	

    asc->nc_current = -1;

    server = (sinfo*) malloc(sizeof(sinfo));
	
    //allocate and initialize the lists
    server->conn_list = (List*)malloc(sizeof(List));	
    server->requests = (List*)malloc(sizeof(List));
    server->pending = (List*)malloc(sizeof(List));
    server->complete = (List*)malloc(sizeof(List));

    list_init(server->requests, free);
    list_init(server->conn_list, free);
    list_init(server->pending, free);
    list_init(server->complete, free);
	
	
    return asc;
}

/* reset the initialization
 */

int ardma_server_finalize(struct ardma_server_connection **asc_p)
{
}


/* look at the network side and get the events on the network
 * look at the worker side and get the worker events
 */

int ardma_server_check_events( struct ardma_server_connection *asc)
{
    int retval = 0;
    ptl_event_t event;
    int which = 0;
    
    //check the event queue for remote events
    do
    {
	//
	retval = PtlEQPoll(&local_info.eqh, 1, 100,
			   &event, &which);
	if(retval == PTL_OK)
	{
	    log_info("successfully pulled event\n\t[%d, (%d,%d), %d]", 
		     event.type, event.initiator.pid, 
		     event.initiator.nid, event.pt_index);
	    
	    //we recieved event - now check the match bits
	    //we can have connection event MATCH_CONN
	    //or data transfer request MATCH_DATA
	    
	    if(event.match_bits == CONN_MATCH)
	    {
		log_info("recieved connection request\n");
		
	    }
	    else if(event.match_bits == DATA_MATCH)
	    {
		log_info("recieved data request\n");

	    }
	    else
	    {
		log_error("Unrecognized recognized\n");
		continue;
	    }
	    
	}
	else if(retval == PTL_EQ_EMPTY)
	{
	    log_info("event queue empty\n");
	    break;
	}
	else
	{
	    //actual error
	    log_error("Error in PtlEQPoll: %s\n", PtlErrorStr(retval));
	    break;
	}
    }while(retval != PTL_EQ_EMPTY);

}

void ardma_server_cb_request (struct ardma_server_connection * asc, 
                              int crank, int nodeid, 
                              uint64_t pg_size, uint64_t idx_size,
                              int timestep, char * path)

{

}

							
void ardma_server_cb_failure( struct ardma_server_connection *asc)
{


}
