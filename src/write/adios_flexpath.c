
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

#include <pthread.h>

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

#if HAVE_FLEXPATH==1

// // evpath libraries
#include <evpath.h>
#include <cod.h>
#include "core/flexpath.h"
#include <sys/queue.h>

/************************* Structure and Type Definitions ***********************/
// used for messages in the control queue
typedef enum {VAR=0, DATA_FLUSH, OPEN, CLOSE, INIT, EVGROUP_FLUSH, DATA_BUFFER} FlexpathMessageType;

// maintains connection information
typedef struct _flexpath_stone {
    int myNum;
    int theirNum;
    int step;
    int opened;
    int created;
    int condition;
    char* contact;
} FlexpathStone;

// maintains variable dimension information
typedef struct _flexpath_var_node {
    char* varName;
    struct _flexpath_var_node* dimensions;
    struct _flexpath_var_node* next;
    int rank;
} FlexpathVarNode;

// used to construct message queues
typedef struct _flexpath_queue_node {
    void* data;
    FlexpathMessageType type;
    struct _flexpath_queue_node* next;
} FlexpathQueueNode;

// used to sanitize names
typedef struct _flexpath_name_table {
    char *originalName;
    char *mangledName;
    LIST_ENTRY(_flexpath_name_table) entries;
} FlexpathNameTable;

// used to sanitize names
typedef struct _flexpath_alt_name {
    char *name;
    FMField *field;
    LIST_ENTRY(_flexpath_alt_name) entries;
} FlexpathAltName;

// used to sanitize names
typedef struct _flexpath_dim_names {
    char *name;
    LIST_HEAD(alts, _flexpath_alt_name) altList;
    LIST_ENTRY(_flexpath_dim_names) entries;
} FlexpathDimNames;

// structure for file data (metadata and buffer)
typedef struct _flexpath_fm_structure {
    FMStructDescRec *format;
    int size;
    unsigned char *buffer;
    FMFormat ioFormat;
    attr_list attrList;	
    LIST_HEAD(tableHead, _flexpath_name_table) nameList;
    LIST_HEAD(dims, _flexpath_dim_names) dimList;
} FlexpathFMStructure;

// information used per each flexpath file
typedef struct _flexpath_write_file_data {
    // MPI stuff
    MPI_Comm mpiComm;
    int rank;
    int size;
    int host_language;
    // EVPath stuff
    EVstone multiStone;
    EVstone sinkStone;

    EVsource dataSource;
    EVsource offsetSource;
    EVsource dropSource;
    EVsource opSource;
    EVsource stepSource;

    EVaction multi_action;
    FlexpathStone* bridges;
    int numBridges;
    attr_list attrs;

    // server state
    int maxQueueSize;
    int openCount;
    int readerStep;
    int writerStep; // how many times has the writer called closed?
    int finalized; // have we finalized?
    int use_ctrl_thread;

    FlexpathFMStructure* fm;
    FlexpathVarNode* askedVars;
    FlexpathVarNode* writtenVars;
    FlexpathVarNode* formatVars;
    FlexpathQueueNode* controlQueue;
    FlexpathQueueNode* dataQueue;   
    pthread_mutex_t openMutex;
    pthread_mutex_t controlMutex;
    pthread_mutex_t dataMutex;
    pthread_cond_t controlCondition;
    pthread_cond_t dataCondition; //fill
    pthread_t ctrl_thr_id;    

    // global array distribution data
    int globalCount;
    evgroup *gp;

    // for maintaining open file list
    struct _flexpath_write_file_data* next;
    char* name;
} FlexpathWriteFileData;

typedef struct _flexpath_write_data {
    int rank;
    FlexpathWriteFileData* openFiles;
    CManager cm;
} FlexpathWriteData;

/************************* Global Variable Declarations *************************/
// used for sanitizing names
#define OPLEN 7
static char opList[OPLEN] = { '+', '-', '*', '/', '.', '>', '<' };
static char *opRepList[OPLEN] = { "_plus_", "_minus_", "_mult_", "_div_", "_dot_", "_greater_", "_less_" };

// used for global communication and data structures
FlexpathWriteData flexpathWriteData;

/**************************** Function Definitions *********************************/

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

char*
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

// add an attr for each dimension to an attr_list
void 
set_attr_dimensions(char* varName, char* altName, int numDims, attr_list attrs) 
{
    char atomName[200] = "";
    char dimNum[10];
    strcat(atomName, FP_DIM_ATTR_NAME);
    strcat(atomName, "_");
    strcat(atomName, varName);
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


attr_list 
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
attr_list 
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
attr_list 
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

void
evgroup_msg_free(void *eventData, void *clientData)
{
    
    /* evgroup *msg = (evgroup*)eventData; */
    /* int num_vars = msg->num_vars; */
    /* int i; */
    /* for (i=0; i<num_vars; i++) { */
    /* 	free(msg->vars[i].offsets); */
    /* } */
    /* free(msg); */
}

void
drop_evgroup_msg_free(void *eventData, void *clientData)
{
    free(eventData);
}

void
update_step_msg_free(void *eventData, void *clientData)
{
        update_step_msg *msg = (update_step_msg*)eventData;
        free(msg);
} 

// free data packets once EVPath is finished with them
void 
data_free(void* eventData, void* clientData) 
{
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*)clientData;
    FMfree_var_rec_elements(fileData->fm->ioFormat, eventData);
    free(eventData);
}

// free op packets once EVPath is finished with them
void 
op_free(void* eventData, void* clientData) 
{
    fp_write_log("OP", "freeing an op message\n");
    op_msg* op = (op_msg*) eventData;
    if (op->file_name) {
        free(op->file_name);
    }
    free(op);
}

// message queue count
int 
queue_count(FlexpathQueueNode** queue) 
{
    if (*queue==NULL) {
        return 0;
    }
    int count = 1;
    FlexpathQueueNode* current = *queue;
    while (current && current->next) {
        count++;
        current = current->next;
    }
    return count;
}

// message queue add to head
void 
threaded_enqueue(
    FlexpathQueueNode **queue, 
    void* item, 
    FlexpathMessageType type, 
    pthread_mutex_t *mutex, 
    pthread_cond_t *condition,
    int max_size) 
{
    pthread_mutex_lock(mutex);
    if (max_size > 0) {
	while (queue_count(queue) > max_size) {
	    pthread_cond_wait(condition, mutex);
	}
    }
    FlexpathQueueNode* newNode = malloc(sizeof(FlexpathQueueNode));
    newNode->data = item;
    newNode->type = type;
    newNode->next = *queue;
    *queue = newNode;
    pthread_cond_broadcast(condition);
    pthread_mutex_unlock(mutex);
}

// remove from tail of a message queue
FlexpathQueueNode* 
threaded_dequeue(
    FlexpathQueueNode **queue, 
    pthread_mutex_t *mutex, 
    pthread_cond_t *condition, 
    int signal_dequeue) 
{
    pthread_mutex_lock(mutex);
    while (queue_count(queue) == 0) {
        pthread_cond_wait(condition, mutex);
    }
    FlexpathQueueNode *tail;
    FlexpathQueueNode *prev = NULL;
    tail = *queue;
    while (tail && tail->next) {
        prev=tail;
        tail=tail->next;
    }
    if (prev) {
        prev->next = NULL;
    } else {
        *queue = NULL;
    }
    pthread_mutex_unlock(mutex);
    if (signal_dequeue==1) {
        pthread_cond_broadcast(condition);
    }
    return tail;
}

// peek at tail of message queue
FlexpathQueueNode* 
threaded_peek(FlexpathQueueNode** queue, 
	      pthread_mutex_t *mutex, 
	      pthread_cond_t *condition) 
{
    pthread_mutex_lock(mutex);
    int q = queue_count(queue);
    if (q == 0) {	
	pthread_cond_wait(condition, mutex);
    }
    FlexpathQueueNode* tail;
    tail = *queue;
    while (tail && tail->next) {
        tail=tail->next;
    }
    pthread_mutex_unlock(mutex);
    return tail;
}

// add new var to a var list
FlexpathVarNode* 
add_var(FlexpathVarNode* queue, char* varName, FlexpathVarNode* dims, int rank)
{
    if (queue) {
        queue->next=add_var(queue->next, varName, dims, rank);
        return queue;
    } else {
        queue = (FlexpathVarNode*) malloc(sizeof(FlexpathVarNode));
        queue->varName = strdup(varName);
        queue->dimensions = dims;
        queue->next = NULL;
        queue->rank = rank;
        return queue;
    }
}

// free a var list
void 
free_vars(FlexpathVarNode* queue)
{
    if (queue) {
        free_vars(queue->next);
        free(queue->varName);
        free(queue);
    }
}

// search a var list
FlexpathVarNode* 
queue_contains(FlexpathVarNode* queue, const char* name, int rank) 
{
    int compare_rank = 0;
    if (rank >= 0 ) {
        compare_rank = 1;
    }
    FlexpathVarNode* tmp = queue;
    while (tmp) {
        if (strcmp(tmp->varName, name)==0) {
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
get_alt_name(char *name, char *dimName) 
{
    int len = strlen(name) + strlen(dimName) + strlen("FPDIM_") + 2;
    char *newName = malloc(sizeof(char) * len);
    strcpy(newName, "FPDIM_");
    strcat(newName, dimName);
    strcat(newName, "_");
    strcat(newName, name);
    return newName;
}

// lookup a dimensions real name
static FlexpathAltName*
find_alt_name(FlexpathFMStructure *currentFm, char *dimName, char *varName) 
{
    char *altName = get_alt_name(varName, dimName);
    FlexpathDimNames *d = NULL;

    // move to dim name in fm dim name list
    for (d = currentFm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
        if (!strcmp(d->name, dimName)) {
	    break;
	}
    }

    // if reached end of list - create list with current dim name at head
    if (d == NULL) {
        d = (FlexpathDimNames *) malloc(sizeof(FlexpathDimNames));
        d->name = dimName;
        LIST_INIT(&d->altList);
        LIST_INSERT_HEAD(&currentFm->dimList, d, entries);
    }

    // create FlexpathAltName structure and field with alternative name in it 
    FlexpathAltName *a = (FlexpathAltName *) malloc(sizeof(FlexpathAltName));
    a->name = altName;
    FMField *field = (FMField *) malloc(sizeof(FMField));
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
int 
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
	uint64_t global = adios_get_dim_value(&dim_list->global_dimension);
	if (global == 0) { return 0;}
        uint64_t *local_offsets = (uint64_t*)malloc(sizeof(uint64_t) * ndims);
        uint64_t *local_sizes = (uint64_t*)malloc(sizeof(uint64_t) * ndims);
	uint64_t *global_sizes = (uint64_t*)malloc(sizeof(uint64_t) * ndims);
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

// creates multiqueue function to handle ctrl messages for given bridge stones 
char *multiqueue_action = "{\n\
    int found = 0;\n\
    int flush_data_count = 0; \n\
    int my_rank = -1;\n\
    attr_list mine;\n\
    if (EVcount_varMsg()>0) {\n\
        EVdiscard_and_submit_varMsg(0, 0);\n\
    }\n\
    if (EVcount_update_step_msg() > 1) {\n\
        EVdiscard_update_step_msg(0);\n\
    }\n\
    if (EVcount_drop_evgroup_msg()>0) {\n\
       if (EVcount_evgroup()>0) {\n\
          EVdiscard_evgroup(0);\n\
       }\n\
       EVdiscard_and_submit_drop_evgroup_msg(0,0);\n\
    }\n\
    if (EVcount_op_msg()>0) {\n\
        op_msg *msg = EVdata_op_msg(0);\n\
        mine = EVget_attrs_op_msg(0);\n\
        found = attr_ivalue(mine, \"fp_dst_rank\");\n\
        if (found > 0) {\n\
            EVdiscard_and_submit_op_msg(found, 0);\n\
        } else {\n\
            EVdiscard_and_submit_op_msg(0,0);\n\
        }\n\
    }\n\
    if (EVcount_flush()>0) {\n\
        flush *c = EVdata_flush(0);\n\
         if (c->type == 2) { \n\
             if (EVcount_evgroup()>0) {\n\
               evgroup *g = EVdata_evgroup(0); \n\
               g->condition = c->condition;\n\
               EVsubmit(c->process_id+1, g);\n\
               EVdiscard_flush(0);\n\
             }\n\
         }\n\
         else if (c->type == 3) {\n\
            if (EVcount_update_step_msg()>0) {\n\
               update_step_msg *stepmsg = EVdata_update_step_msg(0);\n\
               stepmsg->condition = c->condition;\n\
               EVsubmit(c->process_id+1, stepmsg);\n\
               EVdiscard_flush(0);\n\
            }\n\
          }\n\
         else {\n\
            EVdiscard_and_submit_flush(0,0);\n\
            flush_data_count++;\n\
         }\n\
    }\n\
    if (EVcount_anonymous()>0) {\n\
        mine = EVget_attrs_anonymous(0);\n\
        found = attr_ivalue(mine, \"fp_dst_rank\");\n\
        EVdiscard_and_submit_anonymous(found+1,0);\n\
    }\n\
 }";

// sets a field based on data type
void 
set_field(int type, FMFieldList* field_list_ptr, int fieldNo, int* size)
{
    FMFieldList field_list = *field_list_ptr;
    switch (type) {
    case adios_unknown:
	fprintf(stderr, "set_field: Bad Type Error %d\n", type);
	break;

    case adios_unsigned_integer:
	field_list[fieldNo].field_type = strdup("unsigned integer");
	field_list[fieldNo].field_size = sizeof(unsigned int);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(unsigned int);
	break;
	
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
	field_list[fieldNo].field_type = strdup("double");
	field_list[fieldNo].field_size = sizeof(double);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(double);
	break;

    case adios_long_double:
	field_list[fieldNo].field_type = strdup("double");
	field_list[fieldNo].field_size = sizeof(long double);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(long double);
	break;

    case adios_byte:
	field_list[fieldNo].field_type = strdup("char");
	field_list[fieldNo].field_size = sizeof(char);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(char);
	break;

    case adios_long:
	field_list[fieldNo].field_type = strdup("integer");
	field_list[fieldNo].field_size = sizeof(long);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(long);
	break;
	
    case adios_unsigned_long:
	field_list[fieldNo].field_type = strdup("unsigned integer");
	field_list[fieldNo].field_size = sizeof(unsigned long);
	field_list[fieldNo].field_offset = *size;
	*size += sizeof(unsigned long);
	break;

    case adios_complex:
        field_list[fieldNo].field_type = strdup("complex");
        field_list[fieldNo].field_size = sizeof(complex_dummy);
        field_list[fieldNo].field_offset = *size;
        *size += sizeof(complex_dummy);
        break;
        
    case adios_double_complex:
        field_list[fieldNo].field_type = strdup("double_complex");
        field_list[fieldNo].field_size = sizeof(double_complex_dummy);
        field_list[fieldNo].field_offset = *size;
        *size += sizeof(double_complex_dummy);
        break;

    default:
	fprintf(stderr, "set_field: Unknown Type Error %d\n", type);
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

// generic memory check for after mallocs
void 
mem_check(void* ptr, const char* str) 
{
    if (!ptr) {
        adios_error(err_no_memory, "Cannot allocate memory for flexpath %s.", str);
    }
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

// construct an fm structure based off the group xml file
FlexpathFMStructure* 
set_format(struct adios_group_struct *t, 
	   struct adios_var_struct *fields, 
	   FlexpathWriteFileData *fileData)
{
    FMStructDescRec *format = calloc(4, sizeof(FMStructDescRec));
    FlexpathFMStructure *currentFm = calloc(1, sizeof(FlexpathFMStructure));
    LIST_INIT(&currentFm->nameList);
    LIST_INIT(&currentFm->dimList);
    currentFm->format = format;
    format[0].format_name = strdup(t->name);

    if (t->hashtbl_vars->size(t->hashtbl_vars) == 0) {
	adios_error(err_invalid_group, "set_format: No Variables In Group\n");
	return NULL;
    }

    FMFieldList field_list = malloc(sizeof(FMField) * ((int)t->hashtbl_vars->size(t->hashtbl_vars) + 1));
    if (field_list == NULL) {
	adios_error(err_invalid_group, 
		    "set_format: Field List Memory Allocation Failed. t->hashtbl_vars->size: %d\n", 
		    t->hashtbl_vars->size(t->hashtbl_vars));
	return NULL;
    }

    int fieldNo = 0;
    int altvarcount = 0;

    // for each type look through all the fields
    struct adios_var_struct *adios_var;
    for (adios_var = t->vars; adios_var != NULL; adios_var = adios_var->next, fieldNo++) {
	char *fullname = resolve_path_name(adios_var->path, adios_var->name);

	// use the mangled name for the field.
	field_list[fieldNo].field_name = fullname;
        if (fullname!=NULL) {
            int num_dims = 0;
            char atom_name[200] = "";
            FlexpathVarNode *dims=NULL;
            if (adios_var->dimensions) {
                struct adios_dimension_struct *adim = adios_var->dimensions;  
		
                // attach appropriate attrs for dimensions	
                for (; adim != NULL; adim = adim->next) {
                    num_dims++;		    
		    // have to change get_alt_name to append FPVAR at the start of each varname.
                    char *vname = get_dim_name(&adim->dimension);
                    if (vname) {
			//printf("vname: %s\n", vname);
			//char *name = find_fixed_name(currentFm, vname);
			char *aname = get_alt_name(fullname, vname);
			//printf("aname: %s\n", aname);
			dims=add_var(dims, strdup(aname), NULL, 0);
			set_attr_dimensions(fullname, aname, num_dims, fileData->attrs);
		    }
                    char *gname = get_dim_name(&adim->global_dimension);
		    if (gname) {
			fileData->globalCount++;
			//char *name = find_fixed_name(currentFm, gname);
			char *aname = get_alt_name(fullname, gname);
			dims=add_var(dims, strdup(aname), NULL, 0);
			set_attr_dimensions(fullname, aname, num_dims, fileData->attrs);
		    }
		    if (adim->global_dimension.rank > 0) {
			fileData->globalCount++;
		    }
                }
            }
            // attach ndims attr
            strcat(atom_name, FP_NDIMS_ATTR_NAME);
            strcat(atom_name, "_");
            strcat(atom_name, fullname);
            atom_t ndims_atom = attr_atom_from_string(strdup(atom_name));
            add_int_attr(fileData->attrs, ndims_atom, num_dims);
            fileData->formatVars = add_var(fileData->formatVars, fullname, dims, 0);
        }
	// if its a single field
	if (!adios_var->dimensions) {
	    // set the field type size and offset approrpriately
	    set_field(adios_var->type, &field_list, fieldNo, &currentFm->size);
	} else {
	    //it's a vector!
	    struct adios_dimension_struct *d = adios_var->dimensions;
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
		    FlexpathAltName *a = find_alt_name(currentFm, vname, (char*)field_list[fieldNo].field_name);
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
		  
	    switch (adios_var->type) {
	    case adios_unknown:
		fprintf(stderr, "set_format: Bad Type Error\n");
		fieldNo--;
		break;
		      
	    case adios_integer:
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
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
		    (char *) malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "unsigned integer%s", dims);
		field_list[fieldNo].field_size = sizeof(unsigned int);
		      
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  } 
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(unsigned int));  } 
		break;

	    case adios_real:
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
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
		    (char *) malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "double%s", dims);
		field_list[fieldNo].field_size = sizeof(double);
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  } 
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(double));  } 
		break;
		
	    case adios_long_double:
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "double%s", dims);
		field_list[fieldNo].field_size = sizeof(long double);
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  } 
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(long double));  } 
		break;

	    case adios_byte:
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255, "char%s",
			 dims);
		field_list[fieldNo].field_size = sizeof(char);
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  } 
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(char));  } 
		break;

	    case adios_long: // needs to be unsigned integer in ffs
                // to distinguish on reader_side, I have to look at the size also
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "integer%s", dims);
		field_list[fieldNo].field_size = sizeof(long);
		      
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  } 
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(long));  } 
		break;
		
	    case adios_unsigned_long: // needs to be unsigned integer in ffs
                // to distinguish on reader_side, I have to look at the size also
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255,
			 "unsigned integer%s", dims);
		field_list[fieldNo].field_size = sizeof(unsigned long);
		      
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  } 
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(unsigned long));  } 
		break;

	    case adios_complex:
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255, "complex%s",
			 dims);
		field_list[fieldNo].field_size = sizeof(complex_dummy);
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  } 
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(complex_dummy));  } 
		break;

            case adios_double_complex:
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
		snprintf((char *) field_list[fieldNo].field_type, 255, "double_complex%s",
			 dims);
		field_list[fieldNo].field_size = sizeof(double_complex_dummy);
		field_list[fieldNo].field_offset = currentFm->size;
		if (v_offset == 0 ) // pointer to variably sized array
		{ currentFm->size += sizeof(void *);  } 
		else // statically sized array allocated inline
		{  currentFm->size += (v_offset * sizeof(double_complex_dummy)); } 
		break;

	    default:
		adios_error(err_invalid_group, 
			    "set_format: Unknown Type Error %d: name: %s\n", 
			    adios_var->type, field_list[fieldNo].field_name);
		fieldNo--;	      
		return NULL;
		//break;
	    }
	}

	fp_write_log("FORMAT","field: %s, %s, %d, %d\n", 
		     field_list[fieldNo].field_name, 
		     field_list[fieldNo].field_type,
		     field_list[fieldNo].field_size,
		     field_list[fieldNo].field_offset); 
    }

    FlexpathDimNames *d = NULL;
    field_list = (FMFieldList) realloc(field_list, sizeof(FMField) * (altvarcount + (int)t->hashtbl_vars->size(t->hashtbl_vars) + 1));

    for (d = currentFm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
	FlexpathAltName *a = NULL;
	for (a = d->altList.lh_first; a != NULL; a = a->entries.le_next) {
	    a->field->field_offset = currentFm->size;
	    currentFm->size += sizeof(int);
	    memcpy(&field_list[fieldNo], a->field, sizeof(FMField));
	    fieldNo++;
	}
    }

    for (; fieldNo < (t->hashtbl_vars->size(t->hashtbl_vars) + 1+altvarcount); fieldNo++) {
	field_list[fieldNo].field_type = NULL;
	field_list[fieldNo].field_name = NULL;
	field_list[fieldNo].field_offset = 0;
	field_list[fieldNo].field_size = 0;
    }

    format[0].field_list = field_list;
    format[1].format_name = strdup("complex");
    format[1].field_list = complex_dummy_field_list;
    format[2].format_name = strdup("double_complex");
    format[2].field_list = double_complex_dummy_field_list;

    format[0].struct_size = currentFm->size;
    format[1].struct_size = sizeof(complex_dummy);
    format[2].struct_size = sizeof(double_complex_dummy);
    currentFm->buffer = calloc(1, currentFm->size);

    return currentFm;
}

// copies buffer zeroing out arrays that havent been asked for
void* copy_buffer(void* buffer, int rank, FlexpathWriteFileData* fileData)
{
    char* temp = (char*)malloc(fileData->fm->size);
    memcpy(temp, buffer, fileData->fm->size);
    FMField *f = fileData->fm->format[0].field_list;
    while (f->field_name != NULL)
    {
        FlexpathVarNode* a;
        if (!queue_contains(fileData->askedVars, f->field_name, rank)) {
            if ((a=queue_contains(fileData->formatVars, f->field_name, -1)) 
	       && 
	       (a->dimensions != NULL)) {
                FlexpathVarNode* dim = a->dimensions;
                while (dim) {
                    FMField *f2 = fileData->fm->format[0].field_list;
                    while (f2->field_name != NULL) {
                        if (strcmp(f2->field_name, dim->varName)==0) {
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

void 
process_data_flush(FlexpathWriteFileData *fileData, 
		   Flush_msg *flushMsg, 
		   FlexpathQueueNode *dataNode)
{
    //fprintf(stderr, "writer:%d:processing flush for reader:%d:reader_step:%d:writer_step:%d\n", fileData->rank, flushMsg->process_id, fileData->readerStep, fileData->writerStep);
    void* temp = copy_buffer(dataNode->data, flushMsg->process_id, fileData);
   
    fileData->attrs = set_dst_rank_atom(fileData->attrs, flushMsg->process_id);
    fileData->attrs = set_flush_id_atom(fileData->attrs, flushMsg->id);
    
    if (!fileData->bridges[flushMsg->process_id].opened) {
	fileData->bridges[flushMsg->process_id].opened=1;
	fileData->openCount++;
    }
    //EVsubmit_general(fileData->dataSource, temp, data_free, fileData->attrs);
    EVsubmit_general(fileData->dataSource, temp, NULL, fileData->attrs);
    //fprintf(stderr, "writer:%d:processed flush for reader:%d:reader_step:%d:writer_step:%d\n", fileData->rank, flushMsg->process_id, fileData->readerStep, fileData->writerStep);
}

void
process_var_msg(FlexpathWriteFileData *fileData, Var_msg *varMsg)
{
    fileData->askedVars = add_var(fileData->askedVars, 
				  strdup(varMsg->var_name), 
				  NULL, 
				  varMsg->process_id);
}

void
process_open_msg(FlexpathWriteFileData *fileData, op_msg *open)
{
    fileData->bridges[open->process_id].step = open->step;
    fileData->bridges[open->process_id].condition = open->condition;
    if (!fileData->bridges[open->process_id].created) {
	fileData->bridges[open->process_id].created = 1;
	fileData->bridges[open->process_id].myNum = 
	    EVcreate_bridge_action(
		flexpathWriteData.cm, 
		attr_list_from_string(fileData->bridges[open->process_id].contact), 
		fileData->bridges[open->process_id].theirNum);		    
	
	EVaction_set_output(flexpathWriteData.cm, 
			    fileData->multiStone, 
			    fileData->multi_action, 
			    open->process_id+1, 
			    fileData->bridges[open->process_id].myNum);
    }	
	
    if (open->step == fileData->readerStep) {
	pthread_mutex_lock(&fileData->openMutex);
	fileData->openCount++;  
	fileData->bridges[open->process_id].opened = 1;
	pthread_mutex_unlock(&fileData->openMutex);

	op_msg *ack = malloc(sizeof(op_msg));
	ack->file_name = strdup(fileData->name);
	ack->process_id = fileData->rank;
	ack->step = fileData->readerStep;
	ack->type = 2;
	ack->condition = open->condition;
	fileData->attrs = set_dst_rank_atom(fileData->attrs, open->process_id+1);
	EVsubmit_general(fileData->opSource, ack, op_free, fileData->attrs);
    } 
    else if (open->step < fileData->readerStep) {
	log_error("Flexpath method control_thread: Received Past Step Open\n");
    } 
    else {
	fp_write_log("STEP", "recieved op with future step\n");
    }
}

void
process_close_msg(FlexpathWriteFileData *fileData, op_msg *close)
{

    pthread_mutex_lock(&fileData->openMutex);
    fileData->openCount--;
    fileData->bridges[close->process_id].opened=0;
    fileData->bridges[close->process_id].condition = close->condition;
    pthread_mutex_unlock(&fileData->openMutex);

    if (fileData->openCount==0) {
	FlexpathQueueNode* node = threaded_dequeue(&fileData->dataQueue, 
						   &fileData->dataMutex, 
						   &fileData->dataCondition, 1);
	FMfree_var_rec_elements(fileData->fm->ioFormat, node->data);

	drop_evgroup_msg *dropMsg = malloc(sizeof(drop_evgroup_msg));
	dropMsg->step = fileData->readerStep;
	int wait = CMCondition_get(flexpathWriteData.cm, NULL);
	dropMsg->condition = wait;
	EVsubmit_general(fileData->dropSource, dropMsg, drop_evgroup_msg_free, fileData->attrs);
	//EVsubmit_general(fileData->dropSource, dropMsg, NULL, fileData->attrs);
	// Will have to change when not using ctrl thread.
	CMCondition_wait(flexpathWriteData.cm,  wait); 		    
		     
	fileData->readerStep++;
    }
		
    op_msg *ack = malloc(sizeof(op_msg));
    ack->file_name = strdup(fileData->name);
    ack->process_id = fileData->rank;
    ack->step = fileData->readerStep;
    ack->type = 2;
    ack->condition = close->condition;
    fileData->attrs = set_dst_rank_atom(fileData->attrs, close->process_id + 1);
    EVsubmit_general(fileData->opSource, 
		     ack, 
		     op_free, 
		     fileData->attrs);		
		
}


static int 
var_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*) client_data;
    Var_msg* msg = (Var_msg*) vevent;
    EVtake_event_buffer(cm, vevent);
    threaded_enqueue(&fileData->controlQueue, msg, VAR, 
		     &fileData->controlMutex, &fileData->controlCondition, -1);
    return 0;
}

static int 
flush_handler(CManager cm, void* vevent, void* client_data, attr_list attrs) 
{
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*) client_data;
    Flush_msg* msg = (Flush_msg*) vevent;
    int err = EVtake_event_buffer(cm, vevent);
    //fprintf(stderr, "writer:%d:got_flush for reader:%d:reader_step:%d:writer_step:%d\n", fileData->rank, msg->process_id, fileData->readerStep, fileData->writerStep);
    threaded_enqueue(&fileData->controlQueue, msg, DATA_FLUSH, 
		     &fileData->controlMutex, &fileData->controlCondition,
		     -1);
    //fprintf(stderr, "writer:%d:enqueued flush for reader:%d:reader_step:%d:writer_step:%d\n", fileData->rank, msg->process_id, fileData->readerStep, fileData->writerStep);
    return 0;
}

static int
drop_evgroup_handler(CManager cm, void *vevent, void *client_data, attr_list attrs) {
    drop_evgroup_msg *msg = vevent;
    // will have to change when not using control thread.
    CMCondition_signal(cm, msg->condition);    
    return 0;
}

static int 
op_handler(CManager cm, void* vevent, void* client_data, attr_list attrs) 
{
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*) client_data;
    op_msg* msg = (op_msg*) vevent;
    EVtake_event_buffer(cm, vevent);
    if (msg->type == OPEN_MSG) {
        threaded_enqueue(&fileData->controlQueue, msg, OPEN, 
			 &fileData->controlMutex, &fileData->controlCondition, -1);
    } else if (msg->type == CLOSE_MSG) {
        threaded_enqueue(&fileData->controlQueue, msg, CLOSE, 
			 &fileData->controlMutex, &fileData->controlCondition, -1);  			
    }
    return 0;
}

// processes messages from control queue
void 
control_thread(void *arg) 
{
    FlexpathWriteFileData *fileData = (FlexpathWriteFileData*)arg;
    int rank = fileData->rank;
    FlexpathQueueNode *controlMsg;
    FlexpathQueueNode *dataNode;
    while (1) {
	if ((controlMsg = threaded_dequeue(&fileData->controlQueue, 
	    &fileData->controlMutex, &fileData->controlCondition, 0))) {
	    if (controlMsg->type==VAR) {
		Var_msg *varMsg = (Var_msg*) controlMsg->data;
		process_var_msg(fileData, varMsg);
		EVreturn_event_buffer(flexpathWriteData.cm,controlMsg->data);
	    }
	    else if (controlMsg->type==DATA_FLUSH) {
		Flush_msg *flushMsg = (Flush_msg*)controlMsg->data;
		dataNode = threaded_peek(&fileData->dataQueue, 
					 &fileData->dataMutex, 
					 &fileData->dataCondition);
		process_data_flush(fileData, flushMsg, dataNode);

	    }
	    else if (controlMsg->type==OPEN) {
                op_msg *open = (op_msg*) controlMsg->data;
		process_open_msg(fileData, open);                
		EVreturn_event_buffer(flexpathWriteData.cm, open);
            }
	    else if (controlMsg->type==CLOSE) {
                op_msg* close = (op_msg*) controlMsg->data;
		process_close_msg(fileData, close);
		EVreturn_event_buffer(flexpathWriteData.cm, close);
	    }
	    else {
		log_error("control_thread: Unrecognized Control Message\n");
	    }
	}
    }
    return;
}

// adds an open file handle to global open file list
void 
add_open_file(FlexpathWriteFileData* newFile) 
{
    FlexpathWriteFileData* last = flexpathWriteData.openFiles;
    while (last && last->next) {
        last = last->next;
    }
    if (last) {
        last->next = newFile;
    } else {
        flexpathWriteData.openFiles = newFile;
    }
}

// searches for an open file handle
FlexpathWriteFileData* 
find_open_file(char* name) 
{
    FlexpathWriteFileData* file = flexpathWriteData.openFiles;
    while (file && strcmp(file->name, name)) {
        file = file->next;
    }
    return file;
}


// Initializes flexpath write local data structures
extern void 
adios_flexpath_init(const PairStruct *params, struct adios_method_struct *method) 
{
    setenv("CMSelfFormats", "1", 1);
    // global data structure creation
    flexpathWriteData.rank = -1;
    flexpathWriteData.openFiles = NULL;
    
    // setup CM
    flexpathWriteData.cm = CManager_create();
    atom_t CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    char * transport = getenv("CMTransport");
    if (transport == NULL) {
	int listened = 0;
	while (listened == 0) {
	    if (CMlisten(flexpathWriteData.cm) == 0) {
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		fprintf(stderr, "error: writer %d:pid:%d unable to initialize connection manager. Trying again.\n", rank, (int)getpid());
	    } else {
		listened = 1;
	    }
	}
    } else {
	attr_list listen_list = create_attr_list();
	add_attr(listen_list, CM_TRANSPORT, Attr_String, (attr_value)strdup(transport));
	CMlisten_specific(flexpathWriteData.cm, listen_list);
    }
    
    // configuration setup
    //gen_pthread_init();
    setenv("CMSelfFormats", "1", 1);
    
    // fork communications thread
    int forked = CMfork_comm_thread(flexpathWriteData.cm);   
    if (!forked) {
	fprintf(stderr, "Wrtier error forking comm thread\n");
    }
}

extern int 
adios_flexpath_open(struct adios_file_struct *fd, 
		    struct adios_method_struct *method, 
		    MPI_Comm comm) 
{    
    if ( fd == NULL || method == NULL) {
        perr("open: Bad input parameters\n");
        return -1;
    }

    // file creation
    if (find_open_file(method->group->name)) {
        // stream already open
        return 0;
    }

    FlexpathWriteFileData *fileData = malloc(sizeof(FlexpathWriteFileData));
    mem_check(fileData, "fileData");
    memset(fileData, 0, sizeof(FlexpathWriteFileData));
    
    fileData->maxQueueSize = 1;
    fileData->use_ctrl_thread = 1;
    if (method->parameters) {
        sscanf(method->parameters,"QUEUE_SIZE=%d;", &fileData->maxQueueSize);
    }
   
    // setup step state
    fileData->attrs = create_attr_list();
    fileData->openCount = 0;
    //fileData->readerStep = 0;

    pthread_mutex_init(&fileData->controlMutex, NULL);
    pthread_mutex_init(&fileData->dataMutex, NULL);
    pthread_mutex_init(&fileData->openMutex, NULL);
     
    pthread_cond_init(&fileData->controlCondition, NULL);
    pthread_cond_init(&fileData->dataCondition, NULL);

    // communication channel setup
    char writer_info_filename[200];
    char writer_ready_filename[200];
    char reader_info_filename[200];
    char reader_ready_filename[200];   

    int i=0;
    flexpathWriteData.rank = fileData->rank;
    fileData->globalCount = 0;

    // mpi setup
    MPI_Comm_dup(comm, &fileData->mpiComm);

    MPI_Comm_rank((fileData->mpiComm), &fileData->rank);
    MPI_Comm_size((fileData->mpiComm), &fileData->size);
    char *recv_buff = NULL;
    char sendmsg[CONTACT_LENGTH];
    if (fileData->rank == 0) {
        recv_buff = (char *) malloc(fileData->size*CONTACT_LENGTH*sizeof(char));
    }
        
    // send out contact string
    char * contact = attr_list_to_string(CMget_contact_list(flexpathWriteData.cm));
    fileData->multiStone = EValloc_stone(flexpathWriteData.cm);
    fileData->sinkStone = EValloc_stone(flexpathWriteData.cm);
    sprintf(&sendmsg[0], "%d:%s", fileData->multiStone, contact);
    MPI_Gather(sendmsg, CONTACT_LENGTH, MPI_CHAR, recv_buff, 
        CONTACT_LENGTH, MPI_CHAR, 0, (fileData->mpiComm));

    // rank 0 prints contact info to file
    if (fileData->rank == 0) {
        sprintf(writer_info_filename, "%s_%s", fd->name, "writer_info.txt");
        FILE* writer_info = fopen(writer_info_filename,"w");
        for (i=0; i<fileData->size; i++) {
            fprintf(writer_info, "%s\n",&recv_buff[i*CONTACT_LENGTH]); 
        }
        fclose(writer_info);
    }

    // poll file - race condition issues
    FILE* reader_ready = NULL;
    sprintf(reader_ready_filename, "%s_%s", fd->name, "reader_ready.txt");
    while (!reader_ready) {
	reader_ready = fopen(reader_ready_filename,"r");
    }
    fclose(reader_ready);

    // read contact list
    sprintf(reader_info_filename, "%s_%s", fd->name, "reader_info.txt");
    FILE* reader_info = fopen(reader_info_filename, "r");
    if (!reader_info) {
	reader_info = fopen(reader_info_filename, "r");
    }
    char in_contact[CONTACT_LENGTH] = "";
    int numBridges = 0;
    int stone_num;
    // build a bridge per line
    while (fscanf(reader_info, "%d:%s",&stone_num, in_contact)!=EOF) {
	//fprintf(stderr, "reader contact: %d:%s\n", stone_num, in_contact);
        fileData->bridges = realloc(fileData->bridges, sizeof(FlexpathStone) * (numBridges + 1));
        attr_list contact_list = attr_list_from_string(in_contact);
        fileData->bridges[numBridges].opened = 0;
	fileData->bridges[numBridges].created = 0;
        fileData->bridges[numBridges].step = 0;
        fileData->bridges[numBridges].theirNum = stone_num;
        fileData->bridges[numBridges].contact = strdup(in_contact);
        numBridges += 1;
    }
    fileData->numBridges = numBridges;
    fclose(reader_info);

    MPI_Barrier((fileData->mpiComm));
    
    // cleanup of reader files (writer is done with it).
    if (fileData->rank == 0) {
	unlink(reader_info_filename);
	unlink(reader_ready_filename);
    }

	
    //process group format
    struct adios_group_struct *t = method->group;
    if (t == NULL) {
	adios_error(err_invalid_group, "Invalid group.\n");
	return err_invalid_group;
    }
    fileData->host_language = t->adios_host_language_fortran;
    struct adios_var_struct *fields = t->vars;
	
    if (fields == NULL) {
	adios_error(err_invalid_group, "Group has no variables.\n");
	return err_invalid_group;
    }	

    fileData->fm = set_format(t, fields, fileData);


    // attach rank attr and add file to open list
    fileData->name = strdup(method->group->name); 
    add_open_file(fileData);
    atom_t rank_atom = attr_atom_from_string(FP_RANK_ATTR_NAME);
    add_int_attr(fileData->attrs, rank_atom, fileData->rank);   

    //generate multiqueue function that sends formats or all data based on flush msg

    FMStructDescList queue_list[] = {flush_format_list, 
				     var_format_list, 
				     op_format_list, 
				     evgroup_format_list,
				     drop_evgroup_msg_format_list,
				     data_format_list,
				     update_step_msg_format_list,
				     NULL};
    char* q_action_spec = create_multityped_action_spec(queue_list, 
							multiqueue_action); 
    fileData->multi_action = EVassoc_multi_action(flexpathWriteData.cm, 
						  fileData->multiStone, 
						  q_action_spec, 
						  NULL);
    fileData->dataSource = EVcreate_submit_handle(flexpathWriteData.cm, 
						  fileData->multiStone, 
						  fileData->fm->format);						 

    fileData->opSource = EVcreate_submit_handle(flexpathWriteData.cm, 
						fileData->multiStone, 
						op_format_list); 
    
    fileData->offsetSource = EVcreate_submit_handle(flexpathWriteData.cm, 
						    fileData->multiStone, 
						    evgroup_format_list);
    fileData->dropSource = EVcreate_submit_handle(flexpathWriteData.cm, 
						  fileData->multiStone, 
						  drop_evgroup_msg_format_list);
    
    fileData->stepSource = EVcreate_submit_handle(flexpathWriteData.cm,
						  fileData->multiStone, 
						  update_step_msg_format_list);

    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
			    var_format_list, var_handler, fileData);
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
			    op_format_list, op_handler, fileData);
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
			    drop_evgroup_msg_format_list, drop_evgroup_handler, fileData);
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
	flush_format_list, flush_handler, fileData);

    //link multiqueue to sink
    EVaction_set_output(flexpathWriteData.cm, fileData->multiStone, 
        fileData->multi_action, 0, fileData->sinkStone);
	
    FMContext my_context = create_local_FMcontext();
    fileData->fm->ioFormat = register_data_format(my_context, fileData->fm->format);
    
    sprintf(writer_ready_filename, "%s_%s", fd->name, "writer_ready.txt");
    if (fileData->rank == 0) {
        FILE* writer_info = fopen(writer_ready_filename, "w");
        fprintf(writer_info, "ready");
        fclose(writer_info);
    }       
    
    pthread_create(&fileData->ctrl_thr_id, NULL, (void*)&control_thread, fileData);   
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
    FlexpathWriteFileData* fileData = find_open_file(method->group->name);
    FlexpathFMStructure* fm = fileData->fm;

    if (fm == NULL)
    {
	log_error("adios_flexpath_write: something has gone wrong with format registration: %s\n", 
		  f->name);
	return;
    }
    
    FMFieldList flist = fm->format[0].field_list;
    FMField *field = NULL;
    char *fullname = resolve_path_name(f->path, f->name);
    field = internal_find_field(fullname, flist);

    if (field != NULL) {
	//scalar quantity
	if (!f->dimensions) {
	    if (data) {
		//why wouldn't it have data?
		if (f->type == adios_string) {
		    char *tmpstr = strdup((char*)data);
		    set_FMPtrField_by_name(flist, fullname, fm->buffer, tmpstr);
		} else {
		    memcpy(&fm->buffer[field->field_offset], data, field->field_size);
		}

		//scalar quantities can have FlexpathAltNames also so assign those
		if (field->field_name != NULL) {
					
		    FlexpathDimNames *d = NULL;
		    for (d = fm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
			if (!strcmp(d->name, field->field_name)) {
			    //matches
			    //check if there are FlexpathAltNames
			    FlexpathAltName *a = NULL;
			    for (a = d->altList.lh_first; a != NULL; a = a->entries.le_next) {
				if (f->type == adios_string) {
				    char *tmpstr = strdup((char*)data);
				    set_FMPtrField_by_name(flist, fullname, fm->buffer, tmpstr);
				    //(strcpy(&fm->buffer[a->field->field_offset], (char*)data));
				} else {
				    memcpy(&fm->buffer[a->field->field_offset], 
					   data, 
					   a->field->field_size);
				}
			    }
			}
		    }
		}
	    } else {
		log_error("adios_flexpath_write: error with variable creation: %s\n", f->name);
	    }
	} else {
	    //vector quantity
	    if (data) {	    
                struct adios_dimension_struct *dims = f->dimensions;
                int arraysize = 0;
                while (dims) {
                    int size = adios_get_dim_value(&dims->dimension);
                    if (arraysize) {
                        arraysize *= size;
                    }
                    else
                        arraysize = size;
                    dims = dims->next;
                }
                arraysize *= field->field_size;
                void *datacpy = malloc(arraysize);
                //void *temp = get_FMPtrField_by_name(flist, fullname, fm->buffer, 0);
                memcpy(datacpy, data, arraysize);
                set_FMPtrField_by_name(flist, fullname, fm->buffer, datacpy);
	    } else {
		log_error("adios_flexpath_write: no array data found for var: %s. Bad.\n", f->name);	
	    }
	}
    }
}

extern void 
adios_flexpath_close(struct adios_file_struct *fd, struct adios_method_struct *method) 
{
    FlexpathWriteFileData *fileData = find_open_file(method->group->name);
    void *buffer = malloc(fileData->fm->size);    
    memcpy(buffer, fileData->fm->buffer, fileData->fm->size);

    threaded_enqueue(&fileData->dataQueue, buffer, 
		     DATA_BUFFER,
		     &fileData->dataMutex, 
		     &fileData->dataCondition,
		     fileData->maxQueueSize);
    
    // now gather offsets and send them via MPI to root
    struct adios_group_struct * g = fd->group;
    struct adios_pg_struct * pg = fd->pgs_written;
    evgroup *gp = malloc(sizeof(evgroup));    
    gp->group_name = strdup(method->group->name);
    gp->process_id = fileData->rank;
    if (fileData->globalCount == 0 ) {

	gp->num_vars = 0;
	gp->step = fileData->writerStep;
	gp->vars = NULL;
	EVsubmit_general(fileData->offsetSource, gp, evgroup_msg_free, fileData->attrs);
    }

    else {	
        // process local offsets here	
	int num_gbl_vars = 0;
        global_var * gbl_vars = NULL;
	int num_vars = 0;
	int myrank = fileData->rank;
	int commsize = fileData->size;

        while (pg) {
            struct adios_var_struct * list = pg->vars_written;
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
                // flip for fortran here.
                if (fileData->host_language == FP_FORTRAN_MODE) {
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
                            fileData->mpiComm);

                    MPI_Allgather(local_dimensions, arr_size, MPI_BYTE, 
                            all_local_dims, arr_size, MPI_BYTE,
                            fileData->mpiComm);

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

            pg = pg->next;
        }

	gp->num_vars = num_gbl_vars;
	gp->step = fileData->writerStep;
	gp->vars = gbl_vars;
	//fileData->gp = gp;       
    }
    
    update_step_msg *stepmsg = malloc(sizeof(update_step_msg));
    stepmsg->finalized = 0;
    stepmsg->step = fileData->writerStep;
    stepmsg->condition = -1;
    EVsubmit_general(fileData->stepSource, stepmsg, update_step_msg_free, fileData->attrs);
    
    fileData->attrs = set_size_atom(fileData->attrs, fileData->size);
    EVsubmit_general(fileData->offsetSource, gp, evgroup_msg_free, fileData->attrs);

    fileData->writerStep++;
}

// wait until all open files have finished sending data to shutdown
extern void 
adios_flexpath_finalize(int mype, struct adios_method_struct *method) 
{
    FlexpathWriteFileData* fileData = flexpathWriteData.openFiles;
    while (fileData) {

	update_step_msg *stepmsg = malloc(sizeof(update_step_msg));
	stepmsg->finalized = 1;
	stepmsg->step = fileData->writerStep - 1;
	stepmsg->condition = -1;
	EVsubmit_general(fileData->stepSource, stepmsg, update_step_msg_free, fileData->attrs);

        pthread_mutex_lock(&fileData->dataMutex);
        while (fileData->dataQueue != NULL) {
	    pthread_cond_wait(&fileData->dataCondition, &fileData->dataMutex);
	}
	pthread_mutex_unlock(&fileData->dataMutex);

	fileData->finalized = 1;
	fileData = fileData->next;	    
    }
}

// provides unknown functionality
extern enum BUFFERING_STRATEGY 
adios_flexpath_should_buffer (struct adios_file_struct * fd,struct adios_method_struct * method) 
{
    return no_buffering;  
}

extern void 
adios_flexpath_buffer_overflow (struct adios_file_struct * fd, 
                                struct adios_method_struct * method)
{
    // this call never happens without shared buffering
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

    if (v->adata && v->free_data == adios_flag_yes) {   
        adios_method_buffer_free (v->data_size);
        free (v->adata);
        v->data = v->adata = NULL;
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
        else {        
            v->got_buffer = adios_flag_yes;
            v->free_data = adios_flag_yes;
            v->data_size = mem_allowed;
            v->data = *buffer;
        }
    }
    else {    
        adios_method_buffer_free (mem_allowed);
        log_error ("OVERFLOW: Cannot allocate requested buffer of %llu "
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
