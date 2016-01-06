#ifndef _FLEXPATH_H
#define _FLEXPATH_H

#include "core/adios_logger.h"
#include <container_globals.h>
#include <errno.h>


#ifdef _NOMPI
#include <mpirelay_client.h>
#endif


#define CONTACT_LENGTH 1024

#define READER_CONTACT_FILE "reader_info.txt"
#define WRITER_CONTACT_FILE "writer_info.txt"
#define READER_READY_FILE "reader_ready.txt"
#define WRITER_READY_FILE "writer_ready.txt"
#define FP_SENDER_RANK "fp_rank_num"
#define FP_DST_REPLICA_ID "replica_id"
#define FP_DST_RANK "fp_dst_rank"
#define FP_DIM_ATTR_NAME "fp_dim"
#define FP_NDIMS_ATTR_NAME "fp_ndims"

#define CLOSE_MSG 0
#define OPEN_MSG 1
#define ACK_MSG 2
#define INIT_MSG 3
#define EOS_MSG 4            

#define FP_FORTRAN_MODE 1
#define FP_C_MODE 0

//adios_logger(4,1, __VA_ARGS__);

typedef enum {FORMAT, DATA, EVGROUP, STEP } Flush_type;

typedef struct _update_step_msg
{
    int replica_id; // from which replica the msg came.
    int process_id; // from which process_id the message came
    int step;       // the last step the writer has written to this replica.
    int step_count;
    int finalized;
    int condition;
}update_step_msg;

/*
 * Contains the offset information for a variable for all writers.
 * offsets_per_rank is == ndims.
 */
typedef struct _offset_struct
{
    int replica_id;
    int offsets_per_rank;
    int total_offsets;
    uint64_t *local_dimensions;
    uint64_t *local_offsets;
    uint64_t *global_dimensions;
} offset_struct;

typedef struct _var 
{
    char *name;
    int noffset_structs;
    offset_struct * offsets;    
} global_var, *global_var_ptr;

typedef struct _evgroup 
{        
    int replica_id;
    int process_id;
    int condition;
    int num_vars;
    int step;
    global_var *vars;
} evgroup, *evgroup_ptr;

typedef struct _op_msg
{
    int replica_id;
    int process_id;
    char *file_name;
    int type; //4 = end_of_stream, 3 = init, 2 = ack, 1 = open, 0 = close,
    int step;
    int condition;
} op_msg, *op_msg_ptr;
 
typedef struct flush_msg_ 
{
    int replica_id;
    int process_id;
    int step;
    Flush_type type;
    int condition;
    int id;
} Flush_msg, *Flush_msg_ptr;

typedef struct var_msg_ 
{
    int replica_id;
    int process_id;
    int step;
    char* var_name;
    int condition;
} Var_msg, *Var_msg_ptr;

static FMField update_step_msg_field_list[]=
{
    {"replica_id", "integer", sizeof(int), FMOffset(update_step_msg*, replica_id)},
    {"process_id", "integer", sizeof(int), FMOffset(update_step_msg*, process_id)},
    {"step", "integer", sizeof(int), FMOffset(update_step_msg*, step)},
    {"step_count", "integer", sizeof(int), FMOffset(update_step_msg*, step_count)},
    {"finalized", "integer", sizeof(int), FMOffset(update_step_msg*, finalized)},
    {"condition", "integer", sizeof(int), FMOffset(update_step_msg*, condition)},
    {NULL, NULL, 0, 0}
};

static FMField offset_struct_field_list[]=
{
    {"replica_id", "integer", sizeof(int), FMOffset(offset_struct*, replica_id)},
    {"offsets_per_rank", "integer", sizeof(int), FMOffset(offset_struct*, offsets_per_rank)},
    {"total_offsets", "integer", sizeof(int), FMOffset(offset_struct*, total_offsets)},
    {"local_dimensions", "integer[total_offsets]", sizeof(uint64_t), FMOffset(offset_struct*, local_dimensions)},
    {"local_offsets", "integer[total_offsets]", sizeof(uint64_t), FMOffset(offset_struct*, local_offsets)},
    {"global_dimensions", "integer[offsets_per_rank]", sizeof(uint64_t), FMOffset(offset_struct*, global_dimensions)},
    {NULL, NULL, 0, 0}
};

static FMField global_var_field_list[]=
{
    {"name", "string", sizeof(char*), FMOffset(global_var_ptr, name)},
    {"noffset_structs", "integer", sizeof(int), FMOffset(global_var_ptr, noffset_structs)},
    {"offsets", "offset_struct[noffset_structs]", sizeof(offset_struct), FMOffset(global_var_ptr, offsets)},
    {NULL, NULL, 0, 0}
};

static FMField evgroup_field_list[]=
{
    {"relica_id", "integer", sizeof(int), FMOffset(evgroup_ptr, replica_id)},
    {"process_id", "integer", sizeof(int), FMOffset(evgroup_ptr, process_id)},
    {"condition", "integer", sizeof(int), FMOffset(evgroup_ptr, condition)},
    {"num_vars", "integer", sizeof(int), FMOffset(evgroup_ptr, num_vars)},
    {"step", "integer", sizeof(int), FMOffset(evgroup_ptr, step)},
    {"vars", "global_var[num_vars]", sizeof(global_var), FMOffset(evgroup_ptr, vars)},
    {NULL, NULL, 0, 0}
};

static FMField flush_field_list[] =
{   
    {"replica_id", "integer", sizeof(int), FMOffset(Flush_msg_ptr, replica_id)},
    {"process_id", "integer", sizeof(int), FMOffset(Flush_msg_ptr, process_id)},
    {"step", "integer", sizeof(int), FMOffset(Flush_msg_ptr, step)},
    {"type", "integer", sizeof(Flush_type), FMOffset(Flush_msg_ptr, type)},
    {"condition", "integer", sizeof(int), FMOffset(Flush_msg_ptr, condition)},
    {"id", "integer", sizeof(int), FMOffset(Flush_msg_ptr, id)},
    {NULL, NULL, 0, 0}
};

static FMField var_field_list[] =
{
    {"replica_id", "integer", sizeof(int), FMOffset(Var_msg_ptr, replica_id)},
    {"process_id", "integer", sizeof(int), FMOffset(Var_msg_ptr, process_id)},
    {"step", "integer", sizeof(int), FMOffset(Var_msg_ptr, step)},
    {"var_name", "string", sizeof(char*), FMOffset(Var_msg_ptr, var_name)},
    {NULL, NULL, 0, 0}
};


static FMField op_file_field_list[] =
{
    {"replica_id", "integer", sizeof(int), FMOffset(op_msg_ptr, replica_id)},
    {"process_id", "integer", sizeof(int), FMOffset(op_msg_ptr, process_id)},
    {"file_name", "string", sizeof(char*), FMOffset(op_msg_ptr, file_name)},
    {"type", "integer", sizeof(int), FMOffset(op_msg_ptr, type)},
    {"step", "integer", sizeof(int), FMOffset(op_msg_ptr, step)},
    {"condition", "integer", sizeof(int), FMOffset(op_msg_ptr, condition)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec update_step_msg_format_list[]=
{
    {"update_step_msg", update_step_msg_field_list, sizeof(update_step_msg), NULL},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec offset_struct_format_list[] =
{
    {"offset_struct", offset_struct_field_list, sizeof(offset_struct), NULL},
    {NULL, NULL, 0, 0}
};


static FMStructDescRec evgroup_format_list[] =
{   
    {"evgroup", evgroup_field_list, sizeof(evgroup), NULL},
    {"offset_struct", offset_struct_field_list, sizeof(offset_struct), NULL},
    {"global_var", global_var_field_list, sizeof(global_var), NULL},
    {NULL,NULL,0,NULL}
};

static FMStructDescRec flush_format_list[] =
{   
    {"flush", flush_field_list, sizeof(Flush_msg), NULL},
    {NULL,NULL,0,NULL}
};
 
static FMStructDescRec var_format_list[] =
{
    {"varMsg", var_field_list, sizeof(Var_msg), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescRec data_format_list[] =
{
    {"anonymous", NULL, 0, NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescRec op_format_list[] =
{
    {"op_msg", op_file_field_list, sizeof(op_msg), NULL},
    {NULL, NULL, 0, NULL}
};

static char *getFixedName(char *name);


// mode 0 = control contact, mode 1 = data contact
#ifdef _NOMPI
/* char* */
/* gather_contacts2(MPIRelay_client *client, char *endpoint, int root, int myrank); */

/* char* */
/* gather_hostnames2(MPIRelay_client *client, int root, int myrank); */


static char*
gather_contacts2(MPIRelay_client *client, char *endpoint, int root, int myrank)
{
    int worldsize = MPIRelay_client_size(client);
    char *recvbuf = NULL;
    int tmpsize = sizeof(char) * CONTACT_LENGTH;
    if (myrank == root) {
	fprintf(stderr, "here\n");
        recvbuf = (char*)calloc(tmpsize, worldsize);
    }
   
    //gather data endpoints and control endpoints from the other ranks.
    MPIRelay_gather(client, endpoint, tmpsize, MPIChar, recvbuf, tmpsize, MPIChar, 0);
    /* MPI_Gather(endpoint, tmpsize, MPI_CHAR, recvbuf, */
    /* 	       tmpsize, MPI_CHAR, 0, comm); */
    return recvbuf;
}

static char*
gather_hostnames2(MPIRelay_client *client, int root, int myrank)
{
    int worldsize = MPIRelay_client_size(client);
    char hostname[HOSTNAME_SIZE];
    memset(&hostname, '\0', HOSTNAME_SIZE);
    int err = gethostname(hostname, HOSTNAME_SIZE);
    fprintf(stderr, "replica hostname:%s\n", hostname);
    
    if (err == -1) {
	fprintf(stderr, "Error %d getting hostname for rank:%d\n",
		errno, myrank);
	return NULL;
    }
    char *recvbuf = NULL;
    int tmpsize = HOSTNAME_SIZE * sizeof(char);
    if (myrank == root) {
	recvbuf = (char*)calloc(tmpsize, worldsize);	
    }

    MPIRelay_gather(client, hostname, tmpsize, MPIChar, recvbuf, tmpsize, MPIChar, root);
    /* MPI_Gather(hostname, tmpsize, MPI_CHAR, */
    /* 	       recvbuf, tmpsize, MPI_CHAR, */
    /* 	       root, comm); */

    return recvbuf;
}

#else

static inline char*
gather_contacts(MPI_Comm comm, char *endpoint, int root, int myrank)
{
    int worldsize;
    MPI_Comm_size(comm, &worldsize);
    fprintf(stderr, "rank %d sending %s\n", myrank, endpoint);
    char *recvbuf = NULL;
    int tmpsize = sizeof(char) * CONTACT_LENGTH;
    if (myrank == root) {
	fprintf(stderr, "here\n");
        recvbuf = (char*)calloc(tmpsize, worldsize);
    }
   
    //gather data endpoints and control endpoints from the other ranks.
    MPI_Gather(endpoint, tmpsize, MPI_CHAR, recvbuf,
	       tmpsize, MPI_CHAR, 0, comm);
    return recvbuf;
}

static inline char*
gather_hostnames(MPI_Comm comm, int root, int myrank)
{
    int worldsize;
    MPI_Comm_size(comm, &worldsize);
    char hostname[HOSTNAME_SIZE];
    memset(&hostname, '\0', HOSTNAME_SIZE);
    int err = gethostname(hostname, HOSTNAME_SIZE);
    fprintf(stderr, "replica hostname:%s\n", hostname);
    
    if (err == -1) {
	fprintf(stderr, "Error %d getting hostname for rank:%d\n",
		errno, myrank);
	return NULL;
    }
    char *recvbuf = NULL;
    int tmpsize = HOSTNAME_SIZE * sizeof(char);
    if (myrank == root) {
	recvbuf = (char*)calloc(tmpsize, worldsize);	
    }
    MPI_Gather(hostname, tmpsize, MPI_CHAR,
	       recvbuf, tmpsize, MPI_CHAR,
	       root, comm);

    return recvbuf;
}
#endif

#endif
