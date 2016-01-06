/*
    adios_flexpath.c
    uses evpath for io in conjunction with read/read_flexpath.c
*/


#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
//#include "adios_mpi.h"
//#include <mpi.h>
#include <pthread.h>
#include <stdarg.h>
// xml parser
#include <mxml.h>

// add by Kimmy 10/15/2012
#include <sys/types.h>
#include <sys/stat.h>
// end of change

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include "public/adios_mpi.h"
#include "public/adios_error.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"
#include "core/adios_logger.h"
#include "core/globals.h"

#if HAVE_FLEXPATH==1

// // evpath libraries
#include <data_store.h>
#include <data_store_globals.h>

#include <evpath.h>
#include <cod.h>
#include "core/flexpath.h"
#include <sys/queue.h>
//#ifndef _NOMPI
#include <container_replica.h>
#include <containers.h>
//#endif
/************************* Structure and Type Definitions ***********************/
// used for messages in the control queue
typedef enum {
    VAR=0,
    STEP_MSG,
    DATA_FLUSH,
    OPEN,
    CLOSE,
    INIT,
    EVGROUP_BUFFER,
    DATA_BUFFER
} FlexpathMessageType;

// TODO: check to see if replica is still "active".
char *router_func = "{\n\
        int port = -1;\n\
        port = attr_ivalue(event_attrs, \"fp_dst_rank\");\n\
        return port;\n\
     }";

// creates multiqueue function to handle ctrl messages for given bridge stones
// port 0 is always going to terminal stone that has registered handlers.
// CNT: Have to move evgroup flush stuff out of multi_queue stone. Multiqueue
// stone now might actually not serve much of a purpose. Perhaps just a standard
// router stone. Need to move flush step messages out of here too.
char *multiqueue_action = "{\n\
    int found = 0;\n\
    int my_rank = -1;\n\
    attr_list mine;\n\
    if (EVcount_varMsg()>0) {\n\
        EVdiscard_and_submit_varMsg(0, 0);\n\
    }\n\
    if (EVcount_update_step_msg() > 0) {\n\
        update_step_msg *msg = EVdata_update_step_msg(0);\n\
        mine = EVget_attrs_update_step_msg(0);\n\
        found = attr_ivalue(mine, \"replica_id\");\n\
        EVsubmit_attr(found, msg, mine);\n\
        EVdiscard_update_step_msg(0);\n\
    }\n\
    if (EVcount_op_msg()>0) {\n\
        op_msg *msg = EVdata_op_msg(0);\n\
        if (msg->type == 2) {\n\
            mine = EVget_attrs_op_msg(0);\n\
            found = attr_ivalue(mine, \"replica_id\");\n\
            EVsubmit_attr(found, msg, mine);\n\
            EVdiscard_op_msg(0);\n\
        } else {\n\
            EVdiscard_and_submit_op_msg(0,0);\n\
        }\n\
    }\n\
    if (EVcount_evgroup()>0) {\n\
        evgroup *msg = EVdata_evgroup(0);\n\
        mine = EVget_attrs_evgroup(0);\n\
        found = attr_ivalue(mine, \"replica_id\");\n\
        EVsubmit_attr(found, msg, mine);\n\
        EVdiscard_evgroup(0);\n\
    }\n\
    if (EVcount_flush()>0) {\n\
            EVdiscard_and_submit_flush(0,0);\n\
    }\n\
    if (EVcount_anonymous()>0) {\n\
        mine = EVget_attrs_anonymous(0);\n\
        found = attr_ivalue(mine, \"replica_id\");\n\
        EVdiscard_and_submit_anonymous(found,0);\n\
    }\n\
 }";

// maintains variable dimension information
typedef struct _flexpath_var_list
{
    char* var_name;
    struct _flexpath_var_list* dimensions;
    struct _flexpath_var_list* next;
    int rank;
} Flexpath_varlist;

// used to sanitize names
typedef struct _flexpath_name_table
{
    char *originalName;
    char *mangledName;
    LIST_ENTRY(_flexpath_name_table) entries;
} FlexpathNameTable;

// used to sanitize names
typedef struct _flexpath_alt_name
{
    char *name;
    FMField *field;
    LIST_ENTRY(_flexpath_alt_name) entries;
} Flexpath_alt_name;

// used to sanitize names
typedef struct _flexpath_dim_names
{
    char *name;
    LIST_HEAD(alts, _flexpath_alt_name) altList;
    LIST_ENTRY(_flexpath_dim_names) entries;
} Flexpath_dim_names;

// structure for file data (metadata and buffer)
typedef struct _flexpath_fm_structure
{
    FMStructDescRec *format;
    int size;
    unsigned char *buffer;
    FMFormat ioFormat;
    attr_list attrList;
    LIST_HEAD(tableHead, _flexpath_name_table) nameList;
    LIST_HEAD(dims, _flexpath_dim_names) dimList;
} Flexpath_fm;

// maintains connection information
typedef struct _flexpath_stone
{
    int myNum;
    int theirNum;
    int step;
    int opened;
    int created;
    int condition;
    char *contact;
} Flexpath_stone;

typedef struct _Flexpath_writer_replica
{
    char *container_name;
    int replica_id;

    attr_list attrs;
    int num_bridges;    
    Flexpath_stone *bridges;
    EVstone router_stone;
    EVaction faction;
    char *filter; // pass in NULL when creating router action spec.

    int open_count;
    int writer_step;
    int reader_step;
    int active;
    Flexpath_varlist *asked_vars;

    DataStore *datads;
    DataStore *evgroupds;
    DataStore *stepds;

    pthread_mutex_t open_mutex;
    pthread_mutex_t data_mutex;
    pthread_mutex_t evgroup_mutex;
    pthread_mutex_t step_mutex;
    pthread_cond_t data_condition; //fill
    pthread_cond_t evgroup_condition;
    pthread_cond_t step_condition;

    evgroup *gp;
    int opened; // HACK. FIX LATER.
    struct _Flexpath_writer_replica *next;
}Flexpath_writer_replica;

// information used per each flexpath file
typedef struct _flexpath_write_file_data
{
    char *name;
    int host_language;
    // MPI stuff
    MPI_Comm mpiComm;
    int rank;
    int size;

    int replica_id;
    // EVPath stuff
    char *data_endpoint;
    EVstone multiStone;
    EVstone sinkStone;
    EVaction multi_action;

    EVsource data_source;
    EVsource evgroup_source;
    EVsource op_source;
    EVsource step_source;

    attr_list attrs;

    pthread_mutex_t repmutex;
    Flexpath_writer_replica *replicas;
    Flexpath_writer_replica *current;
    //Flexpath_stone *bridges;
    //int num_bridges;

    // server state
    int maxQueueSize;
    int finalized; // have we finalized?
    int use_ctrl_thread;

    DataStore *ctrlds;

    Flexpath_fm* fm;
    Flexpath_varlist *format_vars;
    pthread_t ctrl_thr_id;

    int globalCount;
    // for maintaining open file list
    struct _flexpath_write_file_data* next;
} Flexpath_writer_file;

typedef struct _flexpath_write_data
{
    int tunnel;
    int port;
    int rank;
    Flexpath_writer_file* openFiles;
//#ifndef _NOMPI
    Replica_context *ctx;
//#endif
    CManager cm;
} Flexpath_writer_local;

/************************* Global Variable Declarations *************************/

static Flexpath_writer_local local;
int global = 0;
/**************************** Function Definitions *********************************/

inline void
fp_writer_log(const char *log, const char *format, ...)
{
    //if (local.rank != 1990) return;

    char *env, *tmp = NULL;
    tmp = getenv("FP_DEBUG");

    if (tmp) {
        env = strdup(tmp);
        va_list arguments;
        va_start(arguments, format);
        char header[128];
        char buffer[2056], msg[2056];
        memset(header, 0, 128);
        memset(buffer, 0, 2056);
        memset(msg, 0, 2056);

        sprintf(header, "WRITER:%s", log);
        vsprintf(buffer, format, arguments);
        sprintf(msg, "%s:%s", header, buffer);

        if (strcmp("ALL", env) == 0) {
            fprintf(stderr, "%s\n", msg);
        } else {
            char* env_tok;
            env_tok = strtok(env, ",");
            while (env_tok) {
                if (strcmp(env_tok, log)==0) {
                    fprintf(stderr, "%s\n", msg);
                }
                env_tok = strtok(NULL, ",");
            }

        }
        free(env);
    }
}

static void
reverse_dims(uint64_t *dims, int len)
{
  int i;
  for (i = 0; i<(len/2); i++) {
    uint64_t tmp = dims[i];
    int end = len-1-i;
    //printf("%d %d\n", dims[i], dims[end]);
    dims[i] = dims[end];
    dims[end] = tmp;
  }
}


static char*
resolve_path_name(char *path, char *name)
{
    char *fullname = NULL;
    if (name) {
        if (path) {
            if (strcmp(path, "")) {
                fullname = malloc(strlen(path) + strlen(name) + 2);
                strcpy(fullname, path);
                strcat(fullname, "/");
                strcat(fullname, name);
                return fullname;
            }
        }
        fullname = malloc(strlen(name)+1);
        strcpy(fullname, name);
        return fullname;
    }
    return NULL;
}

static
double dgettimeofday( void )
{
#ifdef HAVE_GETTIMEOFDAY
    double timestamp;
    struct timeval now;
    gettimeofday(&now, NULL);
    timestamp = now.tv_sec + now.tv_usec* 1.0e-6 ;
    return timestamp;
#else
    return -1;
#endif
}

/* char extern_string[] = "double dgettimeofday(); \n"; */
/* cod_extern_entry externs[] = { */
/*     {"dgettimeofday", (void *) dgettimeofday},	        // 0 */
/*     {(void *) 0, (void *) 0} */
/* }; */


// add an attr for each dimension to an attr_list
static void
set_attr_dimensions(char *var_name, char *altName, int numDims, attr_list attrs)
{
    char atomName[200] = "";
    char dimNum[10];
    strcat(atomName, FP_DIM_ATTR_NAME);
    strcat(atomName, "_");
    strcat(atomName, var_name);
    strcat(atomName, "_");
    sprintf(dimNum, "%d", numDims);
    strcat(atomName, dimNum);
    atom_t dimAtom = attr_atom_from_string(atomName);
    add_string_attr(attrs, dimAtom, altName);
    atomName[0] = '\0';
    strcat(atomName, FP_NDIMS_ATTR_NAME);
    strcat(atomName, "_");
    strcat(atomName, altName);

    atom_t ndimsAtom = attr_atom_from_string(atomName);
    add_int_attr(attrs, ndimsAtom, 0);
}

static attr_list
set_dst_replica_atom(attr_list attrs, int value)
{
    atom_t repid_atom = attr_atom_from_string("replica_id");
    int repid;
    if (!get_int_attr(attrs, repid_atom, &repid)) {
	add_int_attr(attrs, repid_atom, value);
    }
    set_int_attr(attrs, repid_atom, value);
    return attrs;
}

static attr_list
set_flush_id_atom(attr_list attrs, int value)
{
    atom_t dst_atom = attr_atom_from_string("fp_flush_id");
    int dst;
    if (!get_int_attr(attrs, dst_atom, &dst)) {
        add_int_attr(attrs, dst_atom, value);
    }
    set_int_attr(attrs, dst_atom, value);
    return attrs;
}

// sets a size atom
attr_list
set_size_atom(attr_list attrs, int value)
{
    atom_t dst_atom = attr_atom_from_string("fp_size");
    int size;
    if (!get_int_attr(attrs, dst_atom, &size)) {
        add_int_attr(attrs, dst_atom, value);
    }
    set_int_attr(attrs, dst_atom, value);
    return attrs;
}

// sets a dst rank atom
static attr_list
set_dst_rank_atom(attr_list attrs, int value)
{
    atom_t dst_atom = attr_atom_from_string("fp_dst_rank");
    int dst;
    if (!get_int_attr(attrs, dst_atom, &dst)) {
        add_int_attr(attrs, dst_atom, value);
    }
    set_int_attr(attrs, dst_atom, value);
    return attrs;
}

// sets a dst condition atom
static attr_list
set_dst_condition_atom(attr_list attrs, int condition)
{
    atom_t dst_atom = attr_atom_from_string("fp_dst_condition");
    int dst;
    if (!get_int_attr(attrs, dst_atom, &dst)) {
	add_int_attr(attrs, dst_atom, condition);
    }
    set_int_attr(attrs, dst_atom, condition);
    return attrs;
}

// free format packets once EVPath is finished with them

static void
free_evgroup(evgroup *msg)
{
    int num_vars = msg->num_vars;
    int i;
    for (i=0; i<num_vars; i++) {
    	free(msg->vars[i].offsets);
    }
    free(msg);
}

static void
update_step_msg_free(void *eventData, void *clientData)
{
        update_step_msg *msg = (update_step_msg*)eventData;
        free(msg);
}

// free data packets once EVPath is finished with them
static void
data_free(void* eventData, void* clientData)
{
    Flexpath_writer_file* file_data = (Flexpath_writer_file*)clientData;
    FMfree_var_rec_elements(file_data->fm->ioFormat, eventData);
    free(eventData);
}

// free op packets once EVPath is finished with them
static void
op_free(void* eventData, void* clientData)
{
    fp_writer_log("OP", "freeing an op message\n");
    op_msg *op = (op_msg*) eventData;
    if (op->file_name) {
        free(op->file_name);
    }
    free(op);
}


// add new var to a var list
static Flexpath_varlist*
add_var(Flexpath_varlist *queue, char *var_name, Flexpath_varlist *dims, int rank)
{
    if (queue) {
        queue->next=add_var(queue->next, var_name, dims, rank);
        return queue;
    } else {
        queue = malloc(sizeof(Flexpath_varlist));
        queue->var_name = strdup(var_name);
        queue->dimensions = dims;
        queue->next = NULL;
        queue->rank = rank;
        return queue;
    }
}

// free a var list
static void
free_vars(Flexpath_varlist* queue)
{
    if (queue) {
        free_vars(queue->next);
        free(queue->var_name);
        free(queue);
    }
}

// search a var list
static Flexpath_varlist*
queue_contains(Flexpath_varlist *queue, const char* name, int rank)
{
    int compare_rank = 0;
    if (rank >= 0 ) {
        compare_rank = 1;
    }
    Flexpath_varlist* tmp = queue;
    while (tmp) {
        if (strcmp(tmp->var_name, name)==0) {
            if (compare_rank) {
                if (tmp->rank == rank) {
                    return tmp;
                }
            } else {
                return tmp;
            }
        }
        tmp = tmp->next;
    }
    return NULL;
}

// returns a name with the dimension prepended
static char*
get_alt_name(const char *name, char *dimName)
{
    int len = strlen(name) + strlen(dimName) + 2;
    char *newName = malloc(sizeof(char) * len);
    strcpy(newName, dimName);
    strcat(newName, "_");
    strcat(newName, name);
    return newName;
}

// lookup a dimensions real name
static Flexpath_alt_name*
find_alt_name(Flexpath_fm *currentFm, char *dimName, const char *var_name)
{
    char *altName = get_alt_name(var_name, dimName);
    Flexpath_dim_names *d = NULL;

    // move to dim name in fm dim name list
    for (d = currentFm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
        if (!strcmp(d->name, dimName)) {
	    break;
	}
    }

    // if reached end of list - create list with current dim name at head
    if (d == NULL) {
        d = malloc(sizeof(Flexpath_dim_names));
        d->name = dimName;
        LIST_INIT(&d->altList);
        LIST_INSERT_HEAD(&currentFm->dimList, d, entries);
    }

    // create Flexpath_alt_name structure and field with alternative name in it
    Flexpath_alt_name *a = malloc(sizeof(Flexpath_alt_name));
    a->name = altName;
    FMField *field = malloc(sizeof(FMField));
    a->field = field;
    field->field_name = strdup(altName);
    // TO FIX: Should really check datatype (another paramater?)
    field->field_type = strdup("integer");
    field->field_size = sizeof(int);
    field->field_offset = -1;
    LIST_INSERT_HEAD(&d->altList, a, entries);
    return a;
}

// populates offsets array
static int
get_var_offsets(struct adios_var_struct *v,
		struct adios_group_struct *g,
		uint64_t **offsets,
		uint64_t **local_dimensions,
		uint64_t **global_dimensions)
{
    struct adios_dimension_struct * dim_list = v->dimensions;

    int ndims = 0;
    while (dim_list) {
        ndims++;
        dim_list = dim_list->next;
    }
    dim_list = v->dimensions;

    if (ndims) {
        uint64_t *local_offsets = malloc(sizeof(uint64_t) * ndims);
        uint64_t *local_sizes = malloc(sizeof(uint64_t) * ndims);
	uint64_t *global_sizes = malloc(sizeof(uint64_t) * ndims);
        int n = 0;
        while (dim_list) {
            local_sizes[n] = (uint64_t)adios_get_dim_value(&dim_list->dimension);
            local_offsets[n] = (uint64_t)adios_get_dim_value(&dim_list->local_offset);
	    global_sizes[n] = (uint64_t)adios_get_dim_value(&dim_list->global_dimension);
            dim_list=dim_list->next;
            n++;
        }
        *offsets = local_offsets;
        *local_dimensions = local_sizes;
	*global_dimensions = global_sizes;
    } else {
        *offsets = NULL;
        *local_dimensions = NULL;
	*global_dimensions = NULL;
    }
    return ndims;
}


// sets a field based on data type
static void
set_field(int type, FMFieldList* field_list_ptr, int fieldNo, int* size)
{
    FMFieldList field_list = *field_list_ptr;
    switch (type) {
    case adios_unknown:
	fprintf(stderr, "set_field: Bad Type Error\n");
	break;

    case adios_unsigned_integer:
	field_list[fieldNo].field_type = strdup("unsigned integer");
	field_list[fieldNo].field_size = sizeof(unsigned int);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(unsigned int);
	break;

    case adios_unsigned_long:
	field_list[fieldNo].field_type = strdup("unsigned long");
	field_list[fieldNo].field_size = sizeof(unsigned long);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(unsigned long);

    case adios_integer:
	field_list[fieldNo].field_type = strdup("integer");
	field_list[fieldNo].field_size = sizeof(int);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(int);
	break;

    case adios_real:
	field_list[fieldNo].field_type = strdup("float");
	field_list[fieldNo].field_size = sizeof(float);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(float);
	break;

    case adios_string:
	field_list[fieldNo].field_type = strdup("string");
	field_list[fieldNo].field_size = sizeof(char *);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(unsigned char *);
	break;

    case adios_double:
	field_list[fieldNo].field_type = strdup("float");
	field_list[fieldNo].field_size = sizeof(double);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(double);
	break;

    case adios_byte:
	field_list[fieldNo].field_type = strdup("char");
	field_list[fieldNo].field_size = sizeof(char);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(char);
	break;

    default:
	fprintf(stderr, "set_field: Unknown Type Error\n");
	break;
    }
    *field_list_ptr = field_list;
}

// find a field in a given field list
static FMField*
internal_find_field(char *name, FMFieldList flist)
{
    FMField *f = flist;
    while (f->field_name != NULL && strcmp(f->field_name, name)) {
	f++;
    }
    return f;
}

static char*
get_dim_name (struct adios_dimension_item_struct *d)
{
    char *vname = NULL;
    if (d->var) {
        vname = resolve_path_name(d->var->path, d->var->name);
    } else if (d->attr) {
        if (d->attr->var)
            vname = d->attr->var->name;
        else
            vname = d->attr->name;
    }
    // else it's a number value, so there is no name
    return vname;
}

static Flexpath_writer_replica*
find_replica_byid(Flexpath_writer_replica *list, int id)
{
    while (list) {
	if (list->replica_id == id)
	    return list;
	list = list->next;
    }
    return NULL;
}

static void
build_bridge(Flexpath_writer_replica *rep, int rank)
{
    Flexpath_stone *reader = &rep->bridges[rank];
    /* printf("\n\nbuilding bridge for writer_id: %d rank: %d contact: %d:%s\n\n", */
    /*        rep->replica_id, rank, reader->theirNum, reader->contact); */
    reader->created = 1;
    reader->myNum = EVcreate_bridge_action(local.cm,
					   attr_list_from_string(reader->contact),
					   reader->theirNum);
    EVaction_set_output(local.cm,
			rep->router_stone,
			rep->faction,
			rank,
			reader->myNum);
}

//#ifndef _NOMPI
static Flexpath_writer_replica*
setup_replica(Flexpath_writer_file *file_data, Replica_info_msg *msg)
{
    Flexpath_writer_replica *rep = malloc(sizeof(Flexpath_writer_replica));
    memset(rep, 0, sizeof(Flexpath_writer_replica));
    rep->writer_step = -1;
    if(msg->state == OFFLINE) {
        printf("\t\tREPLICA OFFLINE!\n");
        rep->active = 0;
    }        
    else
        rep->active = 1;
    rep->next = NULL;
    rep->datads = DSinit(file_data->maxQueueSize);
    
    Flexpath_stone *bridges = malloc(sizeof(Flexpath_stone) * msg->num_readers);
    char temp[CONTACT_LENGTH] = "";
    int their_stone;

    int i;
    for (i = 0; i<msg->num_readers; i++) {
	sscanf(msg->readers[i].endpoint, "%d:%s", &their_stone, temp);
	//printf("rank: %d reader[%d].endpoint: %s\n", file_data->rank, i, msg->readers[i].endpoint);
	bridges[i].created = 0;
	bridges[i].opened = 0;
	bridges[i].step = 0;
	bridges[i].theirNum = their_stone;
	bridges[i].contact = strdup(temp);
    }
    rep->replica_id = msg->replica_id;
    rep->bridges = bridges;
    rep->router_stone = EValloc_stone(local.cm);
    rep->filter = create_router_action_spec(NULL, router_func);
    rep->faction = EVassoc_immediate_action(local.cm, rep->router_stone, rep->filter, NULL);
    rep->num_bridges = msg->num_readers;

    EVaction_set_output(local.cm,
			file_data->multiStone,
			file_data->multi_action,
			rep->replica_id,
			rep->router_stone);


    pthread_mutex_init(&rep->open_mutex, NULL);
    pthread_mutex_init(&rep->data_mutex, NULL);
    pthread_mutex_init(&rep->evgroup_mutex, NULL);
    pthread_mutex_init(&rep->step_mutex, NULL);

    pthread_cond_init(&rep->data_condition, NULL);
    pthread_cond_init(&rep->evgroup_condition, NULL);
    pthread_cond_init(&rep->step_condition, NULL);
    return rep;
}
//#endif

static void
next_replica_rr(Flexpath_writer_file *file_data)
{
    pthread_mutex_lock(&file_data->repmutex);
    if (!file_data->current->next) {
        file_data->current = file_data->replicas;
    }
    else {
	file_data->current = file_data->current->next;
    }
    while (file_data->current->active == 0) {
        file_data->current = file_data->current->next;
        if (!file_data->current) {
            fprintf(stderr,
                    "ERROR: no replicas active! writer_replica: %d, writer_rank: %d\n",
                    file_data->replica_id, file_data->rank);
               
        }
    }
    pthread_mutex_unlock(&file_data->repmutex);
}

static void
add_replica_tolist(Flexpath_writer_replica **list, Flexpath_writer_replica *rep)
{
    if (!(*list)) {
	*list = rep;
    } else {
	Flexpath_writer_replica *tmp = *list;
	while (tmp && tmp->next) {
	    tmp = tmp->next;
	}
	tmp->next = rep;
    }
}

//really a stupid way of balancing queues. need something much better
static void
balance_queues(Flexpath_writer_replica *reps, Flexpath_writer_replica *target, int max)
{
    if (max < 6)
        return;

    int current = 0;
    Flexpath_writer_replica *longest = NULL;
    while (reps) {
        int length = DSget_count(reps->datads);
        if (length > current) {
            current = length;
            longest = reps;
        }
        reps = reps->next;
    }
    if (6 > current)
        return;

    int i;
    for (i = 0; i<3; i++) {
        DataStore_obj *obj = DSget_tail(longest->datads);
        target->writer_step++;
        longest->writer_step--;
        /* printf("queue_balancing. Taking obj %d from replica %d and giving to %d with obj id: %d\n", */
        /*        obj->id, longest->replica_id, target->replica_id, target->writer_step); */
        obj->id = target->writer_step;
        DSadd_obj(target->datads, obj, target->writer_step, obj->tag);
    }    
}

static void
process_decrease_msg(Flexpath_writer_file *file_data, Container_decrease_msg *msg)
{
    int i;
    fprintf(stderr, "process decrease function called.\n");
    pthread_mutex_lock(&file_data->repmutex);
    printf("decrease msg count: %d\n", msg->count);
    for (i = 0; i<msg->count; i++) {
        Flexpath_writer_replica *r = find_replica_byid(file_data->replicas, msg->replica_ids[i]);
        if (!r) {
            fprintf(stderr,
                    "ERROR: Replica: %d not found for decrease by writer_replica: %d\n",
                    i, file_data->replica_id);
        }
        else {
            fprintf(stderr, "setting replica: %d as inactive.\n", r->replica_id);
            r->active = 0;
        }
    }
    pthread_mutex_unlock(&file_data->repmutex);
}

static void
process_update_contact(Flexpath_writer_file *file_data, Update_contact_msg *msg)
{
    int i;
    // CNT: Make a function to convert Replica_info_msg to Flexpath_writer_replica
    for (i = 0; i<msg->num_replicas; i++) {
	//printf("\t\trank: %d processing reader_id: %d\n", file_data->rank, i);
	Replica_info_msg *rep = &msg->replicas[i];
	Flexpath_writer_replica *node = setup_replica(file_data, rep);
	node->attrs = attr_copy_list(file_data->attrs);

	// CNT: setup router stone here. don't add port for all nodes, though.
	pthread_mutex_lock(&file_data->repmutex);
	add_replica_tolist(&file_data->replicas, node);
        balance_queues(file_data->replicas, node, file_data->maxQueueSize);
	pthread_mutex_unlock(&file_data->repmutex);
    }
}


// construct an fm structure based off the group xml file
static Flexpath_fm*
set_format(struct adios_group_struct *t,
	   struct adios_var_struct *fields,
	   Flexpath_writer_file *file_data)
{
    FMStructDescRec *format = malloc(sizeof(FMStructDescRec)*2);
    memset(format, 0, sizeof(FMStructDescRec)*2);

    Flexpath_fm *currentFm = malloc(sizeof(Flexpath_fm));
    memset(currentFm, 0, sizeof(Flexpath_fm));

    LIST_INIT(&currentFm->nameList);
    LIST_INIT(&currentFm->dimList);
    currentFm->format = format;
    format->format_name = strdup(t->name);

    if (t->hashtbl_vars->size(t->hashtbl_vars) == 0) {
	adios_error(err_invalid_group, "set_format: No Variables In Group\n");
	return NULL;
    }

    FMFieldList field_list = malloc(sizeof(FMField) * \
				    ((int)t->hashtbl_vars->size(t->hashtbl_vars) + 1));
    if (field_list == NULL) {
	adios_error(err_invalid_group,
		    "set_format: Field List Memory Allocation Failed.\n");
	return NULL;
    }

    int fieldNo = 0;
    int altvarcount = 0;

    // for each type look through all the fields
    struct adios_var_struct *f;
    for (f = t->vars; f != NULL; f = f->next, fieldNo++) {
	char *fullname = resolve_path_name(f->path, f->name);

	// use the mangled name for the field.
	field_list[fieldNo].field_name = fullname;
        if (fullname!=NULL) {
            int num_dims = 0;
            char atom_name[200] = "";
            Flexpath_varlist *dims=NULL;
            if (f->dimensions) {
                struct adios_dimension_struct *adim = f->dimensions;

                // attach appropriate attrs for dimensions
                for (; adim != NULL; adim = adim->next) {
                    num_dims++;

                    char *vname = get_dim_name(&adim->dimension);
                    if (vname) {
			//char *name = find_fixed_name(currentFm, vname);
			char *aname = get_alt_name(fullname, vname);
			dims=add_var(dims, strdup(aname), NULL, 0);
			set_attr_dimensions(fullname, aname, num_dims, file_data->attrs);
		    }
                    char *gname = get_dim_name(&adim->global_dimension);
		    if (gname) {
			file_data->globalCount++;
			//char *name = find_fixed_name(currentFm, gname);
			char *aname = get_alt_name(fullname, gname);
			dims=add_var(dims, strdup(aname), NULL, 0);
			set_attr_dimensions(fullname, aname, num_dims, file_data->attrs);
		    }
                }
            }
            // attach ndims attr
            strcat(atom_name, FP_NDIMS_ATTR_NAME);
            strcat(atom_name, "_");
            strcat(atom_name, fullname);
            atom_t ndims_atom = attr_atom_from_string(strdup(atom_name));
            add_int_attr(file_data->attrs, ndims_atom, num_dims);
            file_data->format_vars = add_var(file_data->format_vars, fullname, dims, 0);
        }
	// if its a single field
	if (!f->dimensions) {
	    // set the field type size and offset approrpriately
	    set_field(f->type, &field_list, fieldNo, &currentFm->size);
	} else {
	    //it's a vector!
	    struct adios_dimension_struct *d = f->dimensions;
            #define DIMSIZE 10240
	    #define ELSIZE 256
            char dims[DIMSIZE] = "";
	    char el[ELSIZE] = "";
	    int v_offset=-1;

	    //create the textual representation of the dimensions
	    for (; d != NULL; d = d->next) {
                char *vname = get_dim_name(&d->dimension);
                if (vname) {

		    //char *name = find_fixed_name(currentFm, vname);
		    Flexpath_alt_name *a = find_alt_name(currentFm,
						       vname,
						       field_list[fieldNo].field_name);
		    altvarcount++;
		    snprintf(el, ELSIZE, "[%s]", a->name);
		    v_offset = 0;
		} else {
		    snprintf(el, ELSIZE, "[%llu]", d->dimension.rank);
		    v_offset *= d->dimension.rank;
		}
		strncat(dims, el, DIMSIZE);
	    }
	    v_offset *= -1;

	    while (currentFm->size % 8 != 0) {
		currentFm->size ++;
	    }

	    switch (f->type) {
	    case adios_unknown:
		fprintf(stderr, "set_format: Bad Type Error\n");
		fieldNo--;
		break;

	    case adios_integer:
		field_list[fieldNo].field_type =
		    malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "integer%s", dims);
		field_list[fieldNo].field_size = sizeof(int);

		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  }
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(int));  }
		break;

	    case adios_unsigned_integer:
		field_list[fieldNo].field_type =
		    malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "unsigned integer%s", dims);
		field_list[fieldNo].field_size = sizeof(unsigned int);

		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  }
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(unsigned int));  }
		break;

	    case adios_unsigned_long:
		field_list[fieldNo].field_type =
		    malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "unsigned long%s", dims);
		field_list[fieldNo].field_size = sizeof(unsigned long);

		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  }
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(unsigned long));  }
		break;

	    case adios_real:
		field_list[fieldNo].field_type =
		    malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "float%s", dims);
		field_list[fieldNo].field_size = sizeof(float);
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  }
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(float));  }
		break;

	    case adios_string:
		field_list[fieldNo].field_type = strdup("string");
		field_list[fieldNo].field_size = sizeof(char);
		field_list[fieldNo].field_offset = currentFm->size;
		currentFm->size += sizeof(void *);
		break;

	    case adios_double:
		field_list[fieldNo].field_type =
		    malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "float%s", dims);
		field_list[fieldNo].field_size = sizeof(double);
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  }
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(double));  }
		break;

	    case adios_byte:
		field_list[fieldNo].field_type =
		    malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255, "char%s",
			 dims);
		field_list[fieldNo].field_size = sizeof(char);
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  }
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(char));  }
		break;

	    default:
		adios_error(err_invalid_group,
			    "set_format: Unknown Type Error %d: name: %s\n",
			    f->type, field_list[fieldNo].field_name);
		fieldNo--;
		return NULL;
		//break;
	    }
	}

	fp_writer_log("FORMAT","field: %s, %s, %d, %d\n",
		     field_list[fieldNo].field_name,
		     field_list[fieldNo].field_type,
		     field_list[fieldNo].field_size,
		     field_list[fieldNo].field_offset);
    }

    Flexpath_dim_names *d = NULL;
    field_list = realloc(field_list,
			 sizeof(FMField) * (altvarcount +		\
					    (int)t->hashtbl_vars->size(t->hashtbl_vars) \
					    + 1));

    for (d = currentFm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
	Flexpath_alt_name *a = NULL;
	for (a = d->altList.lh_first; a != NULL; a = a->entries.le_next) {
	    a->field->field_offset = currentFm->size;
	    currentFm->size += sizeof(int);
	    memcpy(&field_list[fieldNo], a->field, sizeof(FMField));
	    fieldNo++;
	}
    }

    for (; fieldNo<(t->hashtbl_vars->size(t->hashtbl_vars) + 1+altvarcount); fieldNo++) {
	field_list[fieldNo].field_type = NULL;
	field_list[fieldNo].field_name = NULL;
	field_list[fieldNo].field_offset = 0;
	field_list[fieldNo].field_size = 0;
    }

    format->field_list = field_list;
    currentFm->format->struct_size = currentFm->size;

    currentFm->buffer = malloc(currentFm->size);
    memset(currentFm->buffer, 0, currentFm->size);

    return currentFm;
}

// copies buffer zeroing out arrays that havent been asked for
static void*
copy_buffer(void *buffer,
	    int rank,
	    Flexpath_writer_file *file_data,
	    Flexpath_writer_replica *rep)
{
    char *temp = malloc(file_data->fm->size);
    memcpy(temp, buffer, file_data->fm->size);
    FMField *f = file_data->fm->format->field_list;
    while (f->field_name != NULL) {
        Flexpath_varlist *tmp = queue_contains(rep->asked_vars, f->field_name, rank);
        if (!tmp) {
            Flexpath_varlist *a = queue_contains(file_data->format_vars, f->field_name, -1);
            if ((a)  && (a->dimensions != NULL)) {
                Flexpath_varlist *dim = a->dimensions;
                while (dim) {
                    FMField *f2 = file_data->fm->format->field_list;
                    while (f2->field_name != NULL) {
                        if (strcmp(f2->field_name, dim->var_name)==0) {
                            break;
                        }
                        f2++;
                    }
                    if (f2->field_name != NULL) {
                        memset(&temp[f2->field_offset], 0, f2->field_size);
                    }
                    dim = dim->next;
                }
                memset(&temp[f->field_offset], 0, f->field_size);
            }
        }
        f++;
    }
    return temp;
}

static void
process_step_flush(Flexpath_writer_file *file_data,
		   Flexpath_writer_replica *rep,
		   Flush_msg *flushmsg)
{
    fp_writer_log("DETAILS", 
		  "w_rank: %d, entering process_step_flush step: %d reader_id: %d reader_rank: %d\n",
		  file_data->rank,
		  file_data->current->writer_step,
		  flushmsg->replica_id, flushmsg->process_id);
    
    file_data->attrs = set_dst_replica_atom(file_data->attrs, flushmsg->replica_id);
    file_data->attrs = set_dst_rank_atom(file_data->attrs, flushmsg->process_id);
    file_data->attrs = set_dst_condition_atom(file_data->attrs, flushmsg->condition);

    if (!rep->bridges[flushmsg->process_id].created) {
	build_bridge(rep, flushmsg->process_id);
    }

    update_step_msg *msg = malloc(sizeof(update_step_msg));
    msg->step = rep->writer_step;
    msg->process_id = file_data->rank;
    msg->replica_id = file_data->replica_id;
    msg->step_count = DSget_count(rep->datads);
    int finalized;
    if(rep->active)
        finalized = file_data->finalized;
    else
        finalized = 1; // tricking reader to think we're finalized so it ends up exiting.
    msg->finalized = finalized;

    fp_writer_log("DETAILS",
                  "w_rank: %d sending step: %d to reader_id: %d with queue_count: %d\n",
                  file_data->rank,
                  rep->writer_step,
                  rep->replica_id,
                  msg->step_count);

    EVsubmit_general(file_data->step_source, msg, update_step_msg_free, file_data->attrs);
}

static void
process_evgroup_flush(Flexpath_writer_file *file_data,
		      Flexpath_writer_replica *rep,
		      Flush_msg *flushmsg)
{
    fp_writer_log("DETAILS",
                  "w_rank: %d entering process_evgroup_flush step: %d reader_id: %d reader_rank: %d\n",
		  file_data->rank,
		  file_data->current->writer_step,
		  flushmsg->replica_id,
		  flushmsg->process_id);
    file_data->attrs = set_dst_replica_atom(file_data->attrs, flushmsg->replica_id);
    file_data->attrs = set_dst_rank_atom(file_data->attrs, flushmsg->process_id);
    file_data->attrs = set_dst_condition_atom(file_data->attrs, flushmsg->condition);

    fp_writer_log("DETAILS",
                  "w_rank: %d look evgroup for reader_id: %d reader_rank: %d step: %d\n", 
                  file_data->rank,
                  rep->replica_id,
                  flushmsg->process_id,
                  flushmsg->step);   
    DataStore_obj *node = DSpeek_obj_blocking(rep->datads, flushmsg->step);
    fp_writer_log("DETAILS",
                  "w_rank: %d found evgroup for reader_id: %d reader_rank: %d step: %d\n", 
                  file_data->rank,
                  rep->replica_id,
                  flushmsg->process_id,
                  flushmsg->step);

    evgroup *msg = node->metadata;
    msg->process_id = file_data->rank;
    msg->replica_id = file_data->replica_id;
    EVsubmit_general(file_data->evgroup_source, msg, NULL, file_data->attrs);
    //free_evgroup(msg);
    //free(node);
}

static void
process_data_flush(Flexpath_writer_file *file_data,
		   Flexpath_writer_replica *rep,
		   Flush_msg *flushmsg,
		   void *data)
{
    fp_writer_log("DETAILS",
                  "w_rank: %d process data_flush for reader_id: %d reader_rank: %d step: %d\n",
		   file_data->rank,
                  flushmsg->replica_id, flushmsg->process_id, flushmsg->step);
    void *temp = copy_buffer(data,
                             flushmsg->process_id,
			     file_data,
			     rep);

    file_data->attrs = set_dst_replica_atom(file_data->attrs, flushmsg->replica_id);
    file_data->attrs = set_dst_rank_atom(file_data->attrs, flushmsg->process_id);
    file_data->attrs = set_dst_condition_atom(file_data->attrs, flushmsg->condition);
    file_data->attrs = set_flush_id_atom(file_data->attrs, flushmsg->id);

    EVsubmit_general(file_data->data_source, temp, NULL, file_data->attrs);
    fp_writer_log("DETAILS",
                  "w_rank: %d leaving process data_flush for reader_id: %d reader_rank: %d step: %d\n",
                  file_data->rank,
                  flushmsg->replica_id,
                  flushmsg->process_id,
                  flushmsg->step);
}

static void
process_var_msg(Flexpath_writer_file *file_data, Var_msg *varmsg)
{
    fp_writer_log("DETAILS", 
		  "\t\tw_rank: %d entering process_var_msg step: %d reader_id: %d reader_rank: %d\n",
		  file_data->rank,
                  file_data->current->writer_step,
                  varmsg->replica_id,
                  varmsg->process_id);
    
    pthread_mutex_lock(&file_data->repmutex);
    Flexpath_writer_replica *rep = find_replica_byid(file_data->replicas, varmsg->replica_id);
    pthread_mutex_unlock(&file_data->repmutex);
    rep->asked_vars = add_var(rep->asked_vars,
			      strdup(varmsg->var_name),
			      NULL,
			      varmsg->process_id);
    fp_writer_log("DETAILS", 
		  "\t\tw_rank: %d leaving process_var_msg step: %d reader_id: %d reader_rank: %d\n",
		  file_data->rank,
                  file_data->current->writer_step,
                  varmsg->replica_id, 
		  varmsg->process_id);
}

static void
process_open_msg(Flexpath_writer_file *file_data, Flexpath_writer_replica *rep, op_msg *open)
{

    fp_writer_log("DETAILS",
                  "w_rank: %d process open msg for reader_id: %d reader_rank: %d step: %d\n",
                  file_data->rank,
                  open->replica_id,
                  open->process_id,
                  open->step);
    rep->bridges[open->process_id].step = open->step;
    rep->bridges[open->process_id].condition = open->condition;

    if (!rep->bridges[open->process_id].created) {
	build_bridge(rep, open->process_id);
    }

    if (open->step == rep->reader_step) {
	pthread_mutex_lock(&rep->open_mutex);
	rep->open_count++;

	rep->bridges[open->process_id].opened = 1;
	pthread_mutex_unlock(&rep->open_mutex);

	op_msg *ack = malloc(sizeof(op_msg));
	ack->file_name = strdup(file_data->name);
	ack->replica_id = file_data->replica_id;
	ack->process_id = file_data->rank;
	ack->step = rep->reader_step;
	ack->type = 2;

	file_data->attrs = set_dst_replica_atom(file_data->attrs, open->replica_id);
	file_data->attrs = set_dst_rank_atom(file_data->attrs, open->process_id);
        file_data->attrs = set_dst_condition_atom(file_data->attrs, open->condition);

	EVsubmit_general(file_data->op_source, ack, op_free, file_data->attrs);
    }
    else if (open->step < rep->reader_step) {
        log_error(
            "\t\t\tw_rank: %d Past Open reader_id: %d reader_rank: %d reader_step: %d open_msg step: %d\n",
            file_data->rank,
            open->replica_id,
            open->process_id,
            rep->reader_step,
            open->step);
    }
    else {
        log_error(
            "\t\t\t w_rank: %d op future from reader_id: %d reader_rank: %d reader_step: %d open_msg step: %d\n",
	    file_data->rank,
            open->replica_id,
            open->process_id,
            rep->reader_step,
            open->step);
    }
    fp_writer_log("DETAILS",
                  "w_rank: %d leave process_open for reader_id: %d reader_rank: %d step: %d\n",
                  file_data->rank,
                  open->replica_id,
                  open->process_id,
                  open->step);
}

static void
process_close_msg(Flexpath_writer_file *file_data, op_msg *close)
{
    fp_writer_log("DETAILS",
                  "w_rank: %d processing close msg for reader_id: %d reader_rank: %d reader_step: %d\n",
		  file_data->rank, close->replica_id, close->process_id, close->step);
    pthread_mutex_lock(&file_data->repmutex);
    Flexpath_writer_replica *rep = find_replica_byid(file_data->replicas, close->replica_id);

    pthread_mutex_unlock(&file_data->repmutex);

    pthread_mutex_lock(&rep->open_mutex);
    rep->open_count--;
    rep->bridges[close->process_id].opened=0;
    rep->bridges[close->process_id].condition = close->condition;
    pthread_mutex_unlock(&rep->open_mutex);

    if (rep->open_count==0) {        
        DataStore_obj *node = DSget_obj(rep->datads, close->step);
	FMfree_var_rec_elements(file_data->fm->ioFormat, node->data);

	free_evgroup(node->metadata);
	free(node);
	rep->reader_step++;
	fp_writer_log("DETAILS", 
		      "w_rank: %d dropping data for reader_id: %d reader_rank: %d reader_step: %d\n",
		      file_data->rank, close->replica_id, close->process_id, close->step);
    }

    op_msg *ack = malloc(sizeof(op_msg));
    ack->file_name = strdup(file_data->name);
    ack->process_id = file_data->rank;
    ack->replica_id = file_data->replica_id;
    ack->step = rep->reader_step;
    ack->type = 2;

    file_data->attrs = set_dst_replica_atom(file_data->attrs, close->replica_id);
    file_data->attrs = set_dst_rank_atom(file_data->attrs, close->process_id);
    file_data->attrs = set_dst_condition_atom(file_data->attrs, close->condition);

    EVsubmit_general(file_data->op_source,
		     ack,
		     op_free,
		     file_data->attrs);

    fp_writer_log("DETAILS",
                  "w_rank: %d leaving process close msg for reader_id: %d reader_rank: %d reader_step: %d\n",
		  file_data->rank, close->replica_id, close->process_id, close->step);
}

static int
var_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    Flexpath_writer_file* file_data = (Flexpath_writer_file*) client_data;
    Var_msg* msg = (Var_msg*) vevent;
    EVtake_event_buffer(cm, vevent);    
    DSadd_raw(file_data->ctrlds, msg, NULL, msg->step, VAR);
    return 0;
}

static int
flush_handler(CManager cm, void* vevent, void* client_data, attr_list attrs)
{
    Flexpath_writer_file* file_data = (Flexpath_writer_file*) client_data;
    Flush_msg* msg = (Flush_msg*) vevent;
    int err = EVtake_event_buffer(cm, vevent);
    int type;
    char *logmsg ="";
    if (msg->type == DATA) {
	type = DATA_FLUSH;
	logmsg = "DATA_FLUSH";
    }
    else if (msg->type == EVGROUP) {
	logmsg = "EVGROUP_FLUSH";
	type = EVGROUP_BUFFER;
    }
    else if (msg->type == STEP) {
	logmsg = "STEP_FLUSH";
	type = STEP_MSG;
    }

    DSadd_raw(file_data->ctrlds, msg, NULL, msg->step, type);    
    return 0;
}

static int
op_handler(CManager cm, void* vevent, void* client_data, attr_list attrs)
{
    Flexpath_writer_file *file_data = client_data;
    op_msg *msg = vevent;
    fp_writer_log("HANDLER", 
                  "w_rank: %d entering op_handler from reader_id: %d, reader_rank: %d.\n",
		  file_data->replica_id, file_data->rank, msg->replica_id, msg->process_id);
    EVtake_event_buffer(cm, vevent);
    int type;
    char *logmsg = "";
    if (msg->type == OPEN_MSG) {
	type = OPEN;
	logmsg = "OPEN";
    }
    else if (msg->type == CLOSE_MSG) {
	logmsg = "CLOSE";
	type = CLOSE;
    }
    else
	return 0;

    DSadd_raw(file_data->ctrlds, msg, NULL, msg->step, type);
    
    fp_writer_log(
	"HANDLER", "w_rank: %d leaving op_handler.\n",
        file_data->rank);
    return 0;
}

char*
tag_to_string(FlexpathMessageType tag)
{
    char ret[256];
    switch (tag) {
    case VAR:
        return "var";
    case STEP_MSG:
        return "step";
    case DATA_FLUSH:
        return "data_flush";
    case OPEN:
        return "open";
    case CLOSE:
        return "close";
    case INIT:
        return "init";
    case EVGROUP_BUFFER:
        return "evgroup_buffer";
    case DATA_BUFFER:
        return "DATA_BUFFER";
    default:
        memset(ret, '\0', 256);
        sprintf(ret, "tag %d not found", (int)tag);
        return strdup(ret);
    }    
}

// processes messages from control queue
static void
control_thread(void *arg)
{
    Flexpath_writer_file *file_data = (Flexpath_writer_file*)arg;
    int rank = file_data->rank;
    DataStore_obj *ctrlmsg;
    DataStore_obj *data_node;
  
    while (1) {	      
        ctrlmsg = DSget_head_blocking(file_data->ctrlds);
	if (ctrlmsg) {
	    char *tagstr = tag_to_string((FlexpathMessageType)ctrlmsg->tag);
            fp_writer_log("DETAILS",
                          "w_rank: %d ctrlthread: ctrlmsg->tag: %s\n",
                          rank, tagstr);
	    if (ctrlmsg->tag == VAR) {
		Var_msg *varmsg = (Var_msg*) ctrlmsg->data;
		process_var_msg(file_data, varmsg);
		EVreturn_event_buffer(local.cm, ctrlmsg->data);
	    }
	    else if (ctrlmsg->tag == STEP_MSG) {
		Flush_msg *flushmsg = (Flush_msg*)ctrlmsg->data;
		pthread_mutex_lock(&file_data->repmutex);
		Flexpath_writer_replica *rep = find_replica_byid(file_data->replicas,
								 flushmsg->replica_id);

		pthread_mutex_unlock(&file_data->repmutex);
	
		if (!rep) {
                    DSadd_obj(file_data->ctrlds, ctrlmsg, ctrlmsg->id, ctrlmsg->tag);
		    CMsleep(local.cm, 1);

		}
		else {
		    process_step_flush(file_data, rep, flushmsg);
		    EVreturn_event_buffer(local.cm, flushmsg);
                    free(ctrlmsg);
		}
	    }
	    else if (ctrlmsg->tag == EVGROUP_BUFFER) {
		Flush_msg *flushmsg = (Flush_msg*)ctrlmsg->data;
		pthread_mutex_lock(&file_data->repmutex);
		Flexpath_writer_replica *rep = find_replica_byid(file_data->replicas,
								 flushmsg->replica_id);
		pthread_mutex_unlock(&file_data->repmutex);
		process_evgroup_flush(file_data, rep, flushmsg);
		EVreturn_event_buffer(local.cm, flushmsg);
	    }
	    else if (ctrlmsg->tag == DATA_FLUSH) {
		Flush_msg *flushmsg = (Flush_msg*)ctrlmsg->data;
		pthread_mutex_lock(&file_data->repmutex);

		Flexpath_writer_replica *rep = find_replica_byid(file_data->replicas,
								 flushmsg->replica_id);
		pthread_mutex_unlock(&file_data->repmutex);               
                data_node = DSpeek_head_blocking(rep->datads);
		process_data_flush(file_data, rep, flushmsg, data_node->data);
		EVreturn_event_buffer(local.cm, flushmsg);
                free(ctrlmsg);

	    }
	    else if (ctrlmsg->tag == OPEN) {
		op_msg *open = (op_msg*)ctrlmsg->data;
		pthread_mutex_lock(&file_data->repmutex);

		Flexpath_writer_replica *rep = find_replica_byid(file_data->replicas,
								 open->replica_id);
		pthread_mutex_unlock(&file_data->repmutex);
		if (!rep) {
                    DSadd_obj(file_data->ctrlds, ctrlmsg, ctrlmsg->id, ctrlmsg->tag);
		    CMsleep(local.cm, 1);

		}
		else {
		    process_open_msg(file_data, rep, open);
		    EVreturn_event_buffer(local.cm, open);
		    free(ctrlmsg);
		}
	    }
	    else if (ctrlmsg->tag == CLOSE) {
		op_msg* close = (op_msg*) ctrlmsg->data;
		process_close_msg(file_data, close);
		EVreturn_event_buffer(local.cm, close);
                free(ctrlmsg);
	    }
	    else{
		log_error("control_thread: Unrecognized Control Message\n");
	    }
            fp_writer_log("DETAILS",
                          "w_rank: %d ctrlthread processed ctrlmsg type %s\n",
                          rank, tagstr);
	}
    }
    return;
}

// adds an open file handle to global open file list
static void
add_open_file(Flexpath_writer_file *newFile)
{
    Flexpath_writer_file *last = local.openFiles;
    while (last && last->next) {
        last = last->next;
    }
    if (last) {
        last->next = newFile;
    } else {
        local.openFiles = newFile;
    }
}

// searches for an open file handle
Flexpath_writer_file*
find_open_file(char* name)
{
    Flexpath_writer_file* file = local.openFiles;
    while (file && strcmp(file->name, name)) {
        file = file->next;
    }
    return file;
}

/* void */
/* writer_contact_callback(Replica_info_msg *msg, int num_replicas) */
/* { */
/*     // pull reader contact info. */
/*     // lock list */
/*     // create bridges */
/*     // add them to list */
/*     // unlock list. */
/* } */

static void
setup_stones(Flexpath_writer_file *file_data)
{
    file_data->use_ctrl_thread = 1;
    //file_data->control_queue = NULL;
    file_data->ctrlds = DSinit(DS_NO_QLIMIT);

    //pthread_mutex_init(&file_data->control_mutex, NULL);
    //pthread_cond_init(&file_data->control_condition, NULL);

    int i=0;
    local.rank = file_data->rank;
    //file_data->globalCount = 0;

    file_data->multiStone = EValloc_stone(local.cm);
    file_data->sinkStone = EValloc_stone(local.cm);

    char *string_list;
    file_data->data_endpoint = malloc(sizeof(char)*CONTACT_LENGTH);
    string_list = attr_list_to_string(CMget_contact_list(local.cm));
    sprintf(file_data->data_endpoint, "%d:%s", file_data->multiStone, string_list);
    free(string_list);

    FMStructDescList queue_list[] = {flush_format_list,
				     var_format_list,
				     op_format_list,
				     evgroup_format_list,
				     data_format_list,
				     update_step_msg_format_list,
				     NULL};
    char *q_action_spec = create_multityped_action_spec(queue_list,
							multiqueue_action);
    file_data->multi_action = EVassoc_multi_action(local.cm,
						   file_data->multiStone,
						   q_action_spec,
						   NULL);
    file_data->data_source = EVcreate_submit_handle(local.cm,
						    file_data->multiStone,
						    file_data->fm->format);

    file_data->op_source = EVcreate_submit_handle(local.cm,
						  file_data->multiStone,
						  op_format_list);

    file_data->evgroup_source = EVcreate_submit_handle(local.cm,
						       file_data->multiStone,
						       evgroup_format_list);

    file_data->step_source = EVcreate_submit_handle(local.cm,
						    file_data->multiStone,
						    update_step_msg_format_list);

    EVassoc_terminal_action(local.cm, file_data->sinkStone,
			    var_format_list, var_handler, file_data);
    EVassoc_terminal_action(local.cm, file_data->sinkStone,
			    op_format_list, op_handler, file_data);
    EVassoc_terminal_action(local.cm, file_data->sinkStone,
			    flush_format_list, flush_handler, file_data);

    //link multiqueue to sink
    EVaction_set_output(local.cm,
			file_data->multiStone,
			file_data->multi_action,
			0,
 			file_data->sinkStone);

    FMContext my_context = create_local_FMcontext();
    file_data->fm->ioFormat = register_data_format(my_context, file_data->fm->format);

    pthread_create(&file_data->ctrl_thr_id, NULL, (void*)&control_thread, file_data);
    //EVadd_standard_routines(local.cm, extern_string, externs);
}

static Container_decrease_msg*
exchange_decrease_msg(Flexpath_writer_file *file_data, Container_decrease_msg *op)
{
    if (file_data->rank != 0) {
        op = calloc(1, sizeof(Container_decrease_msg));
    }
    MPI_Bcast(&op->count, 1, MPI_INT, 0, file_data->mpiComm);
    if (file_data->rank != 0) {
        op->replica_ids = calloc(op->count, sizeof(int));
    }
    MPI_Bcast(op->replica_ids, op->count, MPI_INT, 0, file_data->mpiComm);
    return op;
}

//#ifndef _NOMPI
static Update_contact_msg*
exchange_update_msg(Flexpath_writer_file *file_data, Update_contact_msg *op)
{
    if (file_data->rank != 0) {
	op = malloc(sizeof(Update_contact_msg));
    }

    MPI_Bcast(&op->num_replicas, 1, MPI_INT, 0, file_data->mpiComm);
    if (file_data->rank != 0) {
	op->replicas = calloc(op->num_replicas, sizeof(Replica_info_msg));
    }
    
    int i;
    for (i = 0; i < op->num_replicas; i++) {
	MPI_Bcast(&op->replicas[i].replica_id, 1, MPI_INT, 0, file_data->mpiComm);
	MPI_Bcast(&op->replicas[i].num_readers, 1, MPI_INT, 0, file_data->mpiComm);
	int num_readers = op->replicas[i].num_readers;
	int bufsize = sizeof(char) * num_readers * CONTACT_LENGTH;
	char *dstream_readers = malloc(bufsize);
	memset(dstream_readers, 0, bufsize);

	if(file_data->rank == 0) {
	    int pos = 0;
	    int j;
	    for (j = 0; j<num_readers; j++) {
		//printf("writer has: %s\n", op->replicas[i].readers[j].endpoint);
		strncpy(&dstream_readers[pos], op->replicas[i].readers[j].endpoint, CONTACT_LENGTH);
		pos += CONTACT_LENGTH;
	    }
	}
	MPI_Bcast(dstream_readers, bufsize, MPI_BYTE, 0, file_data->mpiComm);

	if (file_data->rank != 0) {
	    op->replicas[i].readers = malloc(sizeof(Contact_info) * num_readers);
	    int j;
	    int pos = 0;
	    for (j = 0; j<num_readers; j++) {
		op->replicas[i].readers[j].endpoint = strndup(&dstream_readers[pos], CONTACT_LENGTH);
		pos += CONTACT_LENGTH;
	    }
	}
    }
    return op;
}

// Initializes flexpath write local data structures
extern void
adios_flexpath_init(const PairStruct *params, struct adios_method_struct *method)
{
    setenv("CMSelfFormats", "1", 1);
    // global data structure creation
    local.rank = -1;
    local.openFiles = NULL;
//#ifndef _NOMPI
    local.ctx = NULL;
//#endif
    // setup CM
    local.cm = CManager_create();
    atom_t CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    char * transport = getenv("CMTransport");

    if (transport == NULL) {
	while (CMlisten(local.cm) == 0) {
	    fprintf(stderr,
		    "writer error: unable to initialize connection manager.\n");
	    //exit(1);
	}
    
    } else {
	attr_list listen_list = create_attr_list();
	add_attr(listen_list, CM_TRANSPORT, Attr_String, (attr_value)strdup(transport));
	CMlisten_specific(local.cm, listen_list);
    }

    // configuration setup
    //gen_pthread_init();
    setenv("CMSelfFormats", "1", 1);

    // fork communications thread
    int forked = CMfork_comm_thread(local.cm);
    if (!forked) {
	fprintf(stderr, "Wrtier error forking comm thread\n");
    }
}

//#endif

int
get_and_process_operation(Flexpath_writer_file *file_data)
{
    int optype = REPLICA_NO_OP;
    void *op = NULL;
    if (file_data->rank == 0) {
        op = replica_get_operation(local.ctx, &optype);
    }

    /* The purpose of this block below is so that we don't let
     * writer proceed with different lists of downstream readers.
     * That would be very bad.	 */
    MPI_Bcast(&optype, 1, MPI_INT, 0, file_data->mpiComm);

    fp_writer_log("BIGBIG","before\n");
    if (optype == REPLICA_ADD_CONTACTS) {
        //printf("should not\n");
        Update_contact_msg *downstream = exchange_update_msg(file_data, (Update_contact_msg*)op);
        //Update_contact_msg *downstream = op;
        process_update_contact(file_data, downstream);
        MPI_Barrier(file_data->mpiComm);
    }
    else if (optype == REPLICA_REMOVE_CONTACTS) {
        Container_decrease_msg *downstream = exchange_decrease_msg(file_data,
                                                                   (Container_decrease_msg*)op);
        fprintf(stderr, "processing decrease message.\n");
        process_decrease_msg(file_data, downstream);
    }
    else if (optype == REPLICA_PAUSE) {

    }
    else {
        // not valid control message at this point. what to do?
    }
    return optype;
}

extern int
adios_flexpath_open(struct adios_file_struct *fd,
		    struct adios_method_struct *method,
		    MPI_Comm comm)
{
    if ( fd == NULL || method == NULL) {
        fprintf(stderr, "open: Bad input parameters\n");
        return -1;
    }
    // file creation
    Flexpath_writer_file *file_data = find_open_file(method->group->name);
    global++;
    if (file_data) {
	fp_writer_log("DETAILS", 
		      "adios_flexpath_open w_rank: %d for reader: %d step: %d reader_step: %d\n",		       
		      file_data->rank, 
		      file_data->current->replica_id, 
		      file_data->current->writer_step+1, 
		      file_data->current->reader_step);
        // CNT
        // stream already open
        // just check for new replicas that were added via the callback.
        // put in MPI_AllReduce call here to make sure everyone has
        // gotten the updated replica.
//#ifndef _NOMPI

        get_and_process_operation(file_data);

	//printf( "after\n");
        file_data->current->writer_step++;
	fp_writer_log("DETAILS", 
		      "flexpath_writer_open leave w_rank: %d for reader: %d step: %d\n",
		      file_data->rank, 
		      file_data->current->replica_id,
                      file_data->current->writer_step);
	fp_writer_log("FUNC", "leaving writer_open when file already exists.\n");

        return 0;
//#endif
    }

    fp_writer_log("BIGBUG", "opening with new file.\n");

    file_data = malloc(sizeof(Flexpath_writer_file));
    memset(file_data, 0, sizeof(Flexpath_writer_file));
    file_data->maxQueueSize = 1;
    file_data->attrs = create_attr_list();
    file_data->current = NULL;

    if (method->parameters) {
        sscanf(method->parameters,"QUEUE_SIZE=%d;", &file_data->maxQueueSize);
    }
    fp_writer_log("BIGBUG", "welp in writer open.\n");
    MPI_Comm_dup(comm, &file_data->mpiComm);
    MPI_Comm_rank((file_data->mpiComm), &file_data->rank);
    MPI_Comm_size((file_data->mpiComm), &file_data->size);

    struct adios_group_struct *t = method->group;
    if (t == NULL) {
	adios_error(err_invalid_group, "Invalid group.\n");
	return err_invalid_group;
    }

    file_data->host_language = t->adios_host_language_fortran;
    
    struct adios_var_struct *fields = t->vars;
    if (fields == NULL) {
	adios_error(err_invalid_group, "Group has no variables.\n");
	return err_invalid_group;
    }

    file_data->fm = set_format(t, fields, file_data);
    file_data->name = strdup(method->group->name);
    setup_stones(file_data);
    atom_t rank_atom = attr_atom_from_string(FP_SENDER_RANK);
    add_int_attr(file_data->attrs, rank_atom, file_data->rank);

    int name_set;
    int id_set;
    char *cname = (char*)globals_adios_get_container_name(&name_set);
    int replica_id = globals_adios_get_application_id(&id_set);

    if (!id_set) {
	replica_id = 1;
    }
    if (!name_set) {
	cname = strdup("lammps");
    }
    fp_writer_log("BIGBUG", "in writer open %s %d\n", cname, replica_id);
    file_data->replica_id = replica_id;

    char data_contact_info[CONTACT_LENGTH];
    memset(&data_contact_info[0], 0, CONTACT_LENGTH);
    // if not using tunneling, do normal shit
    if (!local.tunnel) {
	char *string_list;
	string_list = attr_list_to_string(CMget_contact_list(local.cm));
	sprintf(&data_contact_info[0], "%d:%s", file_data->multiStone, string_list);
	free(string_list);
    }
    else { // else trick the downstream reader
	attr_list string_list = create_attr_list();
	add_int_attr(string_list, attr_atom_from_string("IP_PORT"), local.port-1000);
	add_string_attr(string_list, attr_atom_from_string("IP_HOST"), "localhost");
	add_int_attr(string_list, attr_atom_from_string("IP_ADDR"), 0x7f000001);
	char *tunnel_string = attr_list_to_string(string_list);
	sprintf(&data_contact_info[0], "%d:%s", file_data->multiStone, tunnel_string);

    }
    MPI_Barrier(file_data->mpiComm);
    //printf("cntrank: %d: %s\n", file_data->rank, &data_contact_info[0]);

    fp_writer_log("BIGBUG", "writer init replica.\n");
    fprintf(stderr, "replica_size: %d proc_id: %d\n", file_data->size, file_data->rank);
    local.ctx = replica_init(local.cm,
			     file_data->size,
			     file_data->rank,
			     replica_id,
			     &data_contact_info[0],
			     cname,
			     WRITER_REPLICA);
#ifdef _NOMPI
    char *allcontacts; // fill these out
    char *allhosts;
#endif
    
#ifndef _NOMPI
    char *allcontacts =
	gather_contacts(file_data->mpiComm, &data_contact_info[0], 0, file_data->rank);
    char *allhosts =
	gather_hostnames(file_data->mpiComm, 0, file_data->rank);
#endif
    if (file_data->rank == 0) {
	report_new_replica(local.ctx, allcontacts, allhosts);
    }
    
    int opcode = REPLICA_NO_OP;
    while (opcode == REPLICA_NO_OP) {
        opcode = get_and_process_operation(file_data);
    }
    
//#endif
    file_data->current = file_data->replicas;
    file_data->current->writer_step++;
    add_open_file(file_data);

    //MPI_Barrier((file_data->mpiComm));

    //generate multiqueue function that sends formats or all data based on flush msg
    fp_writer_log("FUNC", "leaving writer_open for new file..\n");
    return 0;
}




//  writes data to multiqueue
extern void
adios_flexpath_write(
    struct adios_file_struct *fd,
    struct adios_var_struct *f,
    const void *data,
    struct adios_method_struct *method)
{
    Flexpath_writer_file* file_data = find_open_file(method->group->name);
    Flexpath_fm* fm = file_data->fm;

    if (fm == NULL)
    {
	log_error("adios_flexpath_write: something has gone wrong with format registration: %s\n",
		  f->name);
	return;
    }

    FMFieldList flist = fm->format->field_list;
    FMField *field = NULL;
    char *fullname = resolve_path_name(f->path, f->name);
    field = internal_find_field(fullname, flist);

    if (field != NULL) {
	//scalar quantity
	if (!f->dimensions) {
	    if (data) {
		//why wouldn't it have data?
		memcpy(&fm->buffer[field->field_offset], data, field->field_size);

		//scalar quantities can have Flexpath_alt_names also so assign those
		if (field->field_name != NULL) {

		    Flexpath_dim_names *d = NULL;
		    for (d = fm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
			if (!strcmp(d->name, field->field_name)) {
			    //matches
			    //check if there are Flexpath_alt_names
			    Flexpath_alt_name *a = NULL;
			    for (a = d->altList.lh_first; a != NULL; a = a->entries.le_next) {
				/* fprintf(stderr, "ALTNAME: %s, DIM: %s FIELD: %s\n",  */
				/* 	a->name, d->name, field->field_name); */
				memcpy(&fm->buffer[a->field->field_offset],
				       data,
				       a->field->field_size);
			    }
			}
		    }
		}
	    } else {
		log_error("adios_flexpath_write: error with variable creation: %s\n", f->name);
	    }
	} else {
	    //vector quantity
	    if (data)
	    {
		//we just need to copy the pointer stored in f->data
                // calculate size
                memcpy(&fm->buffer[field->field_offset], &data, sizeof(void *));

	    } else {
		log_error("adios_flexpath_write: no array data found for var: %s. Bad.\n", f->name);
	    }
	}
    }
}

extern void
adios_flexpath_close(struct adios_file_struct *fd, struct adios_method_struct *method)
{
    Flexpath_writer_file *file_data = find_open_file(method->group->name);

    fp_writer_log("DETAILS", 
                  "w_rank: %d entering adios_flexpath_close for reader_id: %d reader_step: %d\n",
                  file_data->rank,
                  file_data->current->replica_id,
                  file_data->current->reader_step);

    void *buffer = malloc(file_data->fm->size);

    struct adios_group_struct *g2 = fd->group;
    struct adios_var_struct *fields = g2->vars;
    while (fields) {
        if (fields->dimensions) {
            struct adios_dimension_struct *dims = fields->dimensions;
            int total_size = 1;
            //for each dimension
            while (dims) {
                int size = adios_get_dim_value (&dims->dimension);
                total_size *= size;
                dims = dims->next;
            }
            FMFieldList flist = file_data->fm->format->field_list;
            FMField *field = NULL;
	    char *fullname = resolve_path_name(fields->path, fields->name);
            field = internal_find_field(fullname, flist);

            total_size*=field->field_size;
            // malloc size
            void *pointer_data_copy = malloc(total_size);
            // while null
            while (pointer_data_copy==NULL) {
                sleep(1);
                void *pointer_data_copy = malloc(total_size);
                //block
            }

            void *temp = get_FMPtrField_by_name(flist, fields->name, file_data->fm->buffer, 0);
            memcpy(pointer_data_copy, temp, total_size);
            set_FMPtrField_by_name(flist, fields->name, file_data->fm->buffer, pointer_data_copy);
        }
        fields = fields->next;
    }

    memcpy(buffer, file_data->fm->buffer, file_data->fm->size);
    Flexpath_writer_replica *current = file_data->current;

    MPI_Barrier(file_data->mpiComm);

    double enqstart, enqtotal;
    enqstart = MPI_Wtime();
    fp_writer_log("DETAILS", "w_rank: %d enqueuing data for reader_id: %d.\n",
		  file_data->rank, file_data->current->replica_id);
    DataStore_obj *obj = malloc(sizeof(*obj));
    if (!obj) {
        fprintf(stderr, "error: cannot create DataStore_obj\n");
        exit(1);
    }
    obj->data = buffer;
    obj->tag = DATA_BUFFER;
    enqtotal = MPI_Wtime() - enqstart;
//    int c = 0;

    // now gather offsets and send them via MPI to root
    struct adios_group_struct *g = fd->group;
    struct adios_var_struct *list = g->vars;

    evgroup *gp = malloc(sizeof(*gp));
    gp->replica_id = file_data->replica_id;
    gp->process_id = file_data->rank;
    gp->step = file_data->current->writer_step;
    if (file_data->globalCount == 0) { // no global vars.
	gp->num_vars = 0;
	gp->vars = NULL;
    }

    else {
        // process local offsets here
	int num_gbl_vars = 0;
        global_var *gbl_vars = NULL;
	int num_vars = 0;
	int myrank = file_data->rank;
	int commsize = file_data->size;

	double offset_start = dgettimeofday();
	while (list) {
	    char *fullname = resolve_path_name(list->path, list->name);
	    //int num_local_offsets = 0;
	    uint64_t *local_offsets = NULL;
	    uint64_t *local_dimensions = NULL;
	    uint64_t *global_dimensions = NULL; // same at each rank.
	    int num_local_offsets = get_var_offsets(list, g,
						    &local_offsets,
						    &local_dimensions,
						    &global_dimensions);

	    if (file_data->host_language == FP_FORTRAN_MODE) {
		reverse_dims(local_offsets, num_local_offsets);
		reverse_dims(local_dimensions, num_local_offsets);
		reverse_dims(global_dimensions, num_local_offsets);
	    }
	    
	    if (num_local_offsets > 0) {
		uint64_t *all_offsets = NULL;
		uint64_t *all_local_dims = NULL;

		int buf_size = num_local_offsets * commsize * sizeof(uint64_t);
		all_offsets = malloc(buf_size);
		all_local_dims = malloc(buf_size);

		int arr_size = num_local_offsets * sizeof(uint64_t);
		MPI_Allgather(local_offsets, arr_size, MPI_BYTE,
			      all_offsets, arr_size, MPI_BYTE,
			      file_data->mpiComm);

		MPI_Allgather(local_dimensions, arr_size, MPI_BYTE,
			      all_local_dims, arr_size, MPI_BYTE,
			      file_data->mpiComm);

		num_gbl_vars++;
		offset_struct *ostruct = malloc(sizeof(offset_struct));
		ostruct->offsets_per_rank = num_local_offsets;
		ostruct->total_offsets = num_local_offsets * commsize;
		ostruct->local_offsets = all_offsets;
		ostruct->local_dimensions = all_local_dims;
		ostruct->global_dimensions = global_dimensions;

		gbl_vars = realloc(gbl_vars, sizeof(global_var) * num_gbl_vars);
		gbl_vars[num_gbl_vars - 1].name = fullname;
		gbl_vars[num_gbl_vars - 1].noffset_structs = 1;
		gbl_vars[num_gbl_vars - 1].offsets = ostruct;

	    }
	    list=list->next;
	}

	gp->num_vars = num_gbl_vars;
	gp->step = file_data->current->writer_step;
	gp->vars = gbl_vars;
	//file_data->gp = gp;
    }

    file_data->attrs = set_size_atom(file_data->attrs, file_data->size);

    MPI_Barrier(file_data->mpiComm);
    fp_writer_log("DETAILS", 
                  "w_rank: %d enqueuing datads for reader_id: %d with id: %d and tag: %d and count: %d\n",
                  file_data->rank, file_data->current->replica_id,
                  file_data->current->writer_step, obj->tag, DSget_count(file_data->current->datads));
    obj->metadata = gp;
    DSadd_obj(file_data->current->datads, obj, file_data->current->writer_step, DATA_BUFFER);
    fp_writer_log("DETAILS", 
                  "w_rank: %d after enqueuing datads for reader_id: %d.\n",
                  file_data->rank, file_data->current->replica_id);
    MPI_Barrier(file_data->mpiComm);
    fprintf(stderr, 
	    "BARRIER:%d:global:%d to reader_id:%d:step:%d\n", 
	    file_data->rank, global, 
	    file_data->current->replica_id, file_data->current->writer_step);
//#ifndef _NOMPI
    int count = DSget_count(file_data->current->datads);
    if (file_data->rank == 0) {
	report_queue_length(local.ctx,
			    file_data->current->container_name,
			    count,
			    enqtotal,
			    file_data->current->replica_id);
    }
    MPI_Barrier(file_data->mpiComm);
//#endif
    next_replica_rr(file_data);
    fp_writer_log("DETAILS", 
                  "w_rank: %d leaving adios_flexpath_close\n",
                  file_data->rank);
    // CNT: move ahead to next replica to which we are writing.
}

// wait until all open files have finished sending data to shutdown
extern void
adios_flexpath_finalize(int mype, struct adios_method_struct *method)
{
    // first tell manager you are finalizing.
    Flexpath_writer_file *file_data = local.openFiles;
    file_data->finalized = 1;

    /* printf("writer_id: %d rank: %d finalizing.\n", */
    /* 	   file_data->replica_id, file_data->rank); */

    Flexpath_writer_replica *list = file_data->replicas;
    // loop until all of them have queue_count of 0;
    while (list) {
        pthread_mutex_lock(&list->data_mutex);
        int a;
        if ((a = DSget_count(list->datads))) {
	    printf("finalizing: rank: %d queue_count: %d\n", file_data->rank, a);
            pthread_mutex_unlock(&list->data_mutex);
            list = list;
            CMsleep(local.cm, 1);
        } else {
            pthread_mutex_unlock(&list->data_mutex);
            list = list->next;
        }
    }
    fp_writer_log("DETAILS",
            "w_rank: %d LEAVING finalize.\n", file_data->rank);
    //CManager_close(local.cm);
}

// provides unknown functionality
extern enum ADIOS_FLAG
adios_flexpath_should_buffer (struct adios_file_struct * fd,struct adios_method_struct * method)
{
    return adios_flag_no;
}

// provides unknown functionality
extern void
adios_flexpath_end_iteration(struct adios_method_struct *method)
{
}

// provides unknown functionality
extern void
adios_flexpath_start_calculation(struct adios_method_struct *method)
{
}

// provides unknown functionality
extern void
adios_flexpath_stop_calculation(struct adios_method_struct *method)
{
}

// provides unknown functionality
extern void
adios_flexpath_get_write_buffer(struct adios_file_struct *fd,
				struct adios_var_struct *v,
				uint64_t *size,
				void **buffer,
				struct adios_method_struct *method)
{
    uint64_t mem_allowed;

    if (*size == 0) {
        *buffer = 0;
        return;
    }

    if (v->data && v->free_data == adios_flag_yes) {
        adios_method_buffer_free (v->data_size);
        free (v->data);
        v->data = NULL;
    }

    mem_allowed = adios_method_buffer_alloc (*size);
    if (mem_allowed == *size) {
        *buffer = malloc (*size);
        if (!*buffer) {
            adios_method_buffer_free (mem_allowed);
            log_error ("ERROR: Out of memory allocating %llu bytes for %s in %s:%s()\n"
                    ,*size, v->name, __FILE__, __func__
                    );
            v->got_buffer = adios_flag_no;
            v->free_data = adios_flag_no;
            v->data_size = 0;
            v->data = 0;
            *size = 0;
            *buffer = 0;
        }
        else{
            v->got_buffer = adios_flag_yes;
            v->free_data = adios_flag_yes;
            v->data_size = mem_allowed;
            v->data = *buffer;
        }
    }
    else{
        adios_method_buffer_free (mem_allowed);
        log_error ("OVERFLOW: error: cannot allocate requested buffer of %llu "
                         "bytes for %s in %s:%s()\n"
                ,*size
                ,v->name
                ,__FILE__, __func__
                );
        *size = 0;
        *buffer = 0;
    }
}

// should not be called from write, reason for inclusion here unknown
void
adios_flexpath_read(struct adios_file_struct *fd,
		    struct adios_var_struct *f,
		    void *buffer,
		    uint64_t buffer_size,
		    struct adios_method_struct *method)
{
}

#else // print empty version of all functions (if HAVE_FLEXPATH == 0)

void
adios_flexpath_read(struct adios_file_struct *fd,
		    struct adios_var_struct *f,
		    void *buffer,
		    struct adios_method_struct *method)
{
}

extern void
adios_flexpath_get_write_buffer(struct adios_file_struct *fd,
				struct adios_var_struct *f,
				unsigned long long *size,
				void **buffer,
				struct adios_method_struct *method)
{
}

extern void
adios_flexpath_stop_calculation(struct adios_method_struct *method)
{
}

extern void
adios_flexpath_start_calculation(struct adios_method_struct *method)
{
}

extern void
adios_flexpath_end_iteration(struct adios_method_struct *method)
{
}


#endif
