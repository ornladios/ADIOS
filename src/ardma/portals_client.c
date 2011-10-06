/* 
   RDMA Client library for Portals
   Implements the functions in ardma_client.h 
*/

#include "ardma_client.h" 
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





/*
 * External functions
 */

int ardma_verbose_level =  2;
static int connection_established = 0; // 1: when connection to server is made
static int is_initialized = 0;

enum REQUEST_STATE {
  STATE_INIT,              // no outstanding request (e.g. server sent completion message)
  STATE_REQUEST_SENDING,   // client has posted send request for pulling its data
  STATE_REQUEST_SENT,      // ibv_post_send's completion event has been processed
  //   it does NOT mean the data has been pulled
  STATE_FAILED             // previous request for pull failed (the pull failed)
};


ninfo local_info;


typedef struct cinfo_
{
  //maximum and in use memory for this client
  unsigned long long maxmem;
  unsigned long long mem;

  //one client can have potentially multiple connections
  //one connection is by nature associated with only one client

  //while one client can talk to many servers in the first iteration 
  //we won't worry about that


  //the msgs that the client can send
  msg request;
  msg pg;
  msg idx;
  msg conn;
	

  //addressing the client
  einfo self;

  //MPI communicator and info 
  MPI_Comm comm;
  int rank;
  int size;

}cinfo, clientinfo;


typedef struct sinfo_
{
  einfo server;
  int rank;
	
}sinfo, serverinfo;


List *server_list;

	

typedef struct ardma_client_connection {

  //the end points
  sinfo *server;
  cinfo *client;
	
  //for cascaded connection to the server
  int lrank;
  int hrank;

  //connection memory
  unsigned long long mem;

} client_connection;


int ardma_client_init (int verbose_level)
{

    ardma_verbose_level = verbose_level;
    if(is_initialized)
	return 0;
	
	
	
    is_initialized = 1;

    portal_init_common(&local_info);
	

    //at this point local_info contains 
    //the portal information specific to this process
    //we will now set a conditional variable so we can't have this function duplicated

    server_list = (List*) malloc(sizeof(List));
	
    list_init(server_list, free);
	

    return 0;
	
}

static sinfo *readConnectionFile(char *filename)
{
    char line[1024];
    char *buffer;
	
    int pid;
    int nid;
    int pt;
    int match;
    int count = 0;
	

    FILE *f;
    f = fopen(filename, "r");
    if(f)
    {
	while((buffer = fgets(line, 1024, f))!= NULL)
	{
	    log_debug("line%d\t:%s\n", count, line);
	    sscanf(line, "%d %d %d %d",
		   &pid, &nid, &pt, &match);
	    sinfo *server  = (sinfo*)malloc(sizeof(sinfo));
	    server->rank = count++;
	    server->server.id.pid = pid;
	    server->server.id.nid = nid;
	    server->server.pt = pt;
	    server->server.match = match;
				
	    list_ins(server_list, (void*)server);
				
	}			
	sinfo *headserver = (sinfo*)list_head(server_list);
	return headserver;
			
    }
    else
    {
	log_error("%s Error cannot open file", strerror(errno));
	return NULL;
    }
				
	
}


static sinfo *getConnectionInfo(MPI_Comm comm)
{
    int rank;
    int size;
	
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    char line[1024];
    char *buffer;
	
    int pid;
    int nid;
    int pt;
    int match;
    int count = 0;
	

    if(!rank)
    {
	FILE *f;
	f = fopen(ARDMA_CONNECTION_FILENAME, "r");
	if(f)
	{
	    while((buffer = fgets(line, 1024, f))!= NULL)
	    {
		log_debug("line%d\t:%s\n", count, line);
		sscanf(line, "%d %d %d %d",
		       &pid, &nid, &pt, &match);
		sinfo *server  = (sinfo*)malloc(sizeof(sinfo));
		server->rank = count++;
		server->server.id.pid = pid;
		server->server.id.nid = nid;
		server->server.pt = pt;
		server->server.match = match;
				
		list_ins(server_list, (void*)server);
				
	    }
			
	    //we can broadcast the connection data to all the clients now
	    MPI_Bcast(&count, 1, MPI_INT, 0, comm);
	    if(count)
	    {
		//TODO		
		
	    }
	    //lets use the first server
			
	    sinfo *headserver = (sinfo*)list_head(server_list);
	    return headserver;
	}
	else
	{
	    log_error("%s Error cannot open file", strerror(errno));
	    return NULL;
	}
		
		
    }
	
}


	
struct ardma_client_connection * ardma_client_connect (MPI_Comm comm)
{
    sinfo *server = NULL;
    ptl_event_t event;
    req_conn conn_msg;
    cinfo *client = (cinfo*)malloc(sizeof(cinfo));
    int retval = 0;
	
    if(comm)
    {
	//not null
	server = getConnectionInfo(comm);		
	if(server == NULL)
	{
	    log_error("Unable to get connection info on comm\n");
	    return NULL;			
	}
		
    }
    else
    {
	//just read from param file
	server = readConnectionFile("param");
		
    }
	
    if(server == NULL)
    {
	log_error("bad connection information - remote is NULL\n");
	return NULL;
		
    }
	
	
    //init the connection structure 
	
	
    client->comm = comm; //MPI_dup comm maybe?
    MPI_rank(comm, &client->rank);
    MPI_size(comm, &client->size);
	
//    memcpy(&server->server, remote, sizeof(ptlinfo));
	
    client->self.id.pid = local_info.pid.pid;
    client->self.id.nid = local_info.pid.nid;
    client->self.ac = local_info.ac;	
    client->self.match = server->server.match;
    client->self.pt = local_info.index;
	
    //connection protocol
    //client sends connection message to server
    //client polls the queue and waits for reply from server
    //reply contains the match entry for the request buffer
    //on the server side

    client->conn.md.start = &conn_msg;
    client->conn.md.length = sizeof(conn_msg);
    client->conn.md.threshold = PTL_MD_THRESH_INF;
    client->conn.md.max_size = sizeof(conn_msg);
    client->conn.md.options = PTL_MD_OP_PUT | PTL_MD_OP_GET | 
	PTL_MD_MANAGE_REMOTE | PTL_MD_TRUNCATE | 
	PTL_MD_MAX_SIZE;
    client->conn.md.user_ptr = (void*)client;
    client->conn.md.eq_handle = local_info.eqh;

    memcpy(&conn_msg.info, &client->self, sizeof(einfo));
	
	
    retval = PtlMDAttach(local_info.meh,client->conn.md,
			 PTL_RETAIN, &client->conn.h);
    if(retval != PTL_OK)
    {
	log_error("PtlMDAttach: %s\n", PtlErrorStr(retval));
	return NULL;
    }

    log_info("sending connection request to server\n");
	
    retval = PtlPut(client->conn.h,
		    PTL_NO_ACK_REQ,
		    server->server.id, //remote id
		    server->server.pt, //remote index
		    server->server.ac, //remote ac but not used
		    server->server.match, //remote match id that we read from file
		    0, 0); //don't need offset or header dat afor now
    if(retval != PTL_OK)
    {
	log_error("PtlPut failed %s\n", PtlErrorStr(retval));
	return NULL;
		
    }
	
	
    log_info("connection request sent successfully\nwaiting for response\n"); 
	
    do
    {
	
	retval = PtlEQWait(client->conn.md.eq_handle,
			   &event);
	if(retval != PTL_OK)
	{
	    log_warn("PtlEQGet didn't return correctly %s\n", PtlErrorStr(retval));
	    //don't know for sure what to do here - for safety reasons
	    // I think we should return null for now
	    return NULL;
		
	}

	log_debug("event: type=%d\ninitiator:(%d,%d)\n\
               uid = %d\t jid = %d\n			  \
               pt_index = %d\t match = %d\n		  \
               rlength = %d\t mlength = %d\n	  \
               offset = %d\t md_handle = %d\n	  \
               md = %d\t hdr_data = %ull\n		  \
               seq = %d\t sequence = %d\n", 
		  event.type, event.initiator.nid, event.initiator.pid,
		  event.uid, event.jid, event.pt_index, event.match_bits,
		  event.rlength, event.mlength, event.offset, event.md_handle,
		  event.md, event.hdr_data, event.link, event.sequence);

	if(is_end_event(event))
	{
	    log_debug("event is end of sequence of events\n");
	    if(event.type = PTL_EVENT_PUT_END)
	    {
		log_debug("server has sent a response");
		if(event.md.start != client->conn.md.start)
		{
		    log_error("the MD isn't right %p instead of %p\n",
			      event.md.start, client->conn.md.start);
					
					
		}
		break;
	    }
			
	}
		
	
    }while(1);


    //we got the return event from server
    //at this point we can set up the request MD 
    //that we'll use to send all subsequent requests to 
    //the server

    log_info("size = %d, rank = %d\n", conn_msg.size, conn_msg.rank);
    log_info("id = %d, match = %d\n", conn_msg.id, conn_msg.info.match);
	
	
    req_data *r = (req_data*)malloc(sizeof(req_data));
	
    memset(r, 0, sizeof(req_data));
	
    client->request.md.start = (void*)r;
    client->request.md.length = sizeof(req_data);
    client->request.md.threshold = PTL_MD_THRESH_INF;
    client->request.md.max_size = sizeof(req_data);
    client->request.md.options = PTL_MD_TRUNCATE;
    client->request.md.user_ptr = (void*)client;
    client->request.buffer = (void*)r;
    client->request.size = sizeof(req_data);
    client->request.match = conn_msg.info.match;
	

    //now we can easily send requests to the server
	
    //construct the connection structure and return
	
    struct ardma_client_connection * acc = 
	(struct ardma_client_connection*) malloc(sizeof(struct ardma_client_connection));
	
    memset(acc, 0, sizeof(*acc));
	
    acc->client = client;
    acc->server = server;

				
    return acc;
}

	
/* Disconnect from staging server */
int ardma_client_disconnect (struct ardma_client_connection **acc_p)
{
    return 0;
}

void * ardma_client_register_pg_memory (struct ardma_client_connection *acc,
                                        uint64_t size, void * addr)
{
    return 0;
}

void * ardma_client_register_index_memory (struct ardma_client_connection *acc,
                                         uint64_t size, void * addr)
{
    return 0;
}

int ardma_client_deregister_pg_memory(struct ardma_client_connection *acc, void * addr)
{
    return 0;
}

int ardma_client_deregister_index_memory(struct ardma_client_connection *acc, void * addr)
{
    return 0;
}

int ardma_client_send_request(struct ardma_client_connection *acc, 
                              uint64_t pg_size, 
                              uint64_t idx_size,
                              char *path,
                              int timestep)
{
    return 0;
}

/* check if previous staging has completed */
int ardma_client_operation_completed(struct ardma_client_connection *acc)
{

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


/*
 * local static functions
 */

static ptl_md_t initmd()
{
    ptl_md_t md;
    return md;

}

static int destroymd()
{
    return -13;
}


