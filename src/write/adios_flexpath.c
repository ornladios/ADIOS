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

#include "public/adios_mpi.h"
#include "public/adios_error.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"
#include "core/adios_logger.h"

// // system libraries
// #include <stdio.h>
// #include <stdlib.h>
#if HAVE_FLEXPATH==1

// // evpath libraries
#include <evpath.h>
//#include <gen_thread.h>

// // local libraries
// #include "config.h"
// #include "core/adios_internals.h"
// #include "core/adios_transport_hooks.h"
// #include "core/util.h"
// #include "public/adios.h"
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

    // EVPath stuff
    EVstone multiStone;
    EVstone sinkStone;
    EVsource formatSource;
    EVsource dataSource;
    EVsource offsetSource;
    EVsource dropSource;
    EVsource opSource;
    EVsource stepSource;
    EVaction multi_action;
    FlexpathStone* bridges;
    int numBridges;
    attr_list attrs;

    int num_reader_coordinators;
    int *reader_coordinators;

    // server state
    int maxQueueSize;
    int openCount;
    int currentStep;
    int writerStep; // how many times has the writer called closed?
    int finalized; // have we finalized?

    pthread_mutex_t openMutex;
    FlexpathFMStructure* fm;
    FlexpathVarNode* askedVars;
    FlexpathVarNode* writtenVars;
    FlexpathVarNode* formatVars;
    FlexpathQueueNode* controlQueue;
    FlexpathQueueNode* dataQueue;   
    FlexpathQueueNode* evgroupQueue;
    pthread_mutex_t controlMutex;
    pthread_mutex_t dataMutex;
    pthread_mutex_t dataMutex2;
    pthread_mutex_t evgroupMutex;
    pthread_cond_t controlCondition;
    pthread_cond_t dataCondition; //fill
    pthread_cond_t dataCondition2; //empty
    pthread_cond_t evgroupCondition;
    pthread_t ctrl_thr_id;    

    // global array distribution data
    int globalCount;
    int sentGlobalOffsets;
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

// add an attr for each dimension to an attr_list
void set_attr_dimensions(char* varName, char* altName, int numDims, attr_list attrs) {
    fp_write_log("ATTR", "adding dim attr %s and ndim attr %d\n", varName, numDims);
    char atomName[200] = "";
    char dimNum[10];
    strcat(atomName, varName);
    strcat(atomName, "_");
    strcat(atomName, FP_DIM_ATTR_NAME);
    strcat(atomName, "_");
    sprintf(dimNum, "%d", numDims);
    strcat(atomName, dimNum);
    atom_t dimAtom = attr_atom_from_string(atomName);
    add_string_attr(attrs, dimAtom, altName);
    atomName[0] = '\0';
    strcat(atomName, altName);
    strcat(atomName, "_");
    strcat(atomName, FP_NDIMS_ATTR_NAME);
    atom_t ndimsAtom = attr_atom_from_string(atomName);
    add_int_attr(attrs, ndimsAtom, 0);
}

// free format packets once EVPath is finished with them
void format_free(void* eventData, void* clientData) {
    fp_write_log("FORMAT", "freeing a format message\n");
    Format_msg* format = (Format_msg*) eventData;
    free(format);
}

// free data packets once EVPath is finished with them
void data_free(void* eventData, void* clientData) {
    fp_write_log("DATA", "freeing a data message\n");
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*)clientData;
    FMfree_var_rec_elements(fileData->fm->ioFormat, eventData);
    free(eventData);
}

// free op packets once EVPath is finished with them
void op_free(void* eventData, void* clientData) {
    fp_write_log("OP", "freeing an op message\n");
    op_msg* op = (op_msg*) eventData;
    if(op->file_name) {
        free(op->file_name);
    }
    free(op);
}

void step_free(void* eventData, void* clientData) {
    fp_write_log("OP", "freeing an update_step message\n");
    free(eventData);
}


// message queue add to head
void threaded_enqueue(FlexpathQueueNode** queue, void* item, FlexpathMessageType type, pthread_mutex_t *mutex, pthread_cond_t *condition) {
    fp_write_log("QUEUE", "enqueing a message\n");
    pthread_mutex_lock(mutex);
    fp_write_log("MUTEX","lock 2\n");
    FlexpathQueueNode* newNode = (FlexpathQueueNode*) malloc(sizeof(FlexpathQueueNode));
    newNode->data = item;
    newNode->type = type;
    newNode->next = *queue;
    *queue = newNode;
    pthread_cond_broadcast(condition);
    fp_write_log("MUTEX","unlock 2\n");
    pthread_mutex_unlock(mutex);
}

// message queue count
int queue_count(FlexpathQueueNode** queue) {
    fp_write_log("QUEUE", "counting a queue\n");
    if(*queue==NULL) {
        return 0;
    }
    int count = 1;
    FlexpathQueueNode* current = *queue;
    while(current && current->next) {
        count++;
        current = current->next;
    }
    fp_write_log("QUEUE", "returning count\n");
    return count;
}

// remove from tail of a message queue
FlexpathQueueNode* threaded_dequeue(FlexpathQueueNode** queue, 
				    pthread_mutex_t *mutex, 
				    pthread_cond_t *condition, 
				    pthread_cond_t *condition2, 
				    int signal_dequeue) {
    fp_write_log("QUEUE", "dequeue\n");
    pthread_mutex_lock(mutex);
    fp_write_log("MUTEX","lock 4\n");
    while(*queue==NULL) {
        fp_write_log("QUEUE", "queue is null\n");
        pthread_cond_wait(condition, mutex);
    }
    FlexpathQueueNode* tail;
    FlexpathQueueNode* prev = NULL;
    tail = *queue;
    while(tail && tail->next) {
        prev=tail;
        tail=tail->next;
    }
    if(prev) {
        prev->next = NULL;
    } else {
        *queue = NULL;
    }
    fp_write_log("MUTEX","unlock 4\n");
    pthread_mutex_unlock(mutex);
    fp_write_log("QUEUE", "exiting dequeue queue:%p ret:%p\n", *queue, tail);
    if(signal_dequeue==1) {
        pthread_cond_broadcast(condition2);
    }
    return tail;
}

// peek at tail of message queue
FlexpathQueueNode* threaded_peek(FlexpathQueueNode** queue, pthread_mutex_t *mutex, pthread_cond_t *condition) {
    int q = queue_count(queue);
    fp_write_log("QUEUE", "peeking at a queue\n");
    fp_write_log("QUEUE", "queue count %d\n", q);
    pthread_mutex_lock(mutex);
    fp_write_log("QUEUE", "recieved lock\n");
    fp_write_log("MUTEX","lock 5\n");
    if(*queue==NULL) {
        fp_write_log("QUEUE", "null about to wait\n");
        pthread_cond_wait(condition, mutex);
        fp_write_log("QUEUE", "signaled with queue %p\n", *queue);
    }
    FlexpathQueueNode* tail;
    tail = *queue;
    while(tail && tail->next) {
        tail=tail->next;
    }
    fp_write_log("MUTEX","unlock 5\n");
    pthread_mutex_unlock(mutex);
    fp_write_log("QUEUE", "returning %p\n", tail);
    return tail;
}

// add new var to a var list
FlexpathVarNode* add_var(FlexpathVarNode* queue, char* varName, FlexpathVarNode* dims, int rank){
    if(queue) {
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
void free_vars(FlexpathVarNode* queue){
    if(queue) {
        free_vars(queue->next);
        free(queue->varName);
        free(queue);
    }
}

// search a var list
FlexpathVarNode* queue_contains(FlexpathVarNode* queue, const char* name, int rank) {
    int compare_rank = 0;
    if(rank >= 0 ) {
        compare_rank = 1;
    }
    FlexpathVarNode* tmp = queue;
    while(tmp) {
        if(strcmp(tmp->varName, name)==0){
            if(compare_rank) {
                if(tmp->rank == rank) {
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

// sanitize a name
char* get_fixed_name(char* name) {
    char* oldName = strdup(name);
    char* newName = (char*) malloc(sizeof(char) * 255);
    int i;
    for (i=0; i< OPLEN; i++){
        char op[] = {opList[i], '\0'};
        char* opRep=opRepList[i];
        char* token = strtok(oldName, op);
        char* lastTok=NULL;
	strcpy(newName, "");
	while(token != NULL){
	    strcat(newName, token);
            if((token = strtok(NULL, op))) {
	        strcat(newName, opRep);
	        lastTok = token;
	    }
        }
        if(lastTok!=NULL && (strlen(newName)-strlen(lastTok)-1>0)) {
            newName[strlen(newName)-strlen(lastTok)-1]='\0';
	}
        free(oldName);
	oldName = strdup(newName);
    }
    free(oldName);
    return newName;
}

// return name with operators removed by using the lookup list
static char* find_fixed_name(FlexpathFMStructure *fm, char *name) {
    FlexpathNameTable *node;
    for (node = fm->nameList.lh_first; node != NULL; node = node->entries.le_next) {
        if (!strcmp(node->originalName, name)) {
	    return node->mangledName;
        }
    }
    return name;
}

// returns a name with the dimension prepended
static char *get_alt_name(char *name, char *dimName) {
    int len = strlen(name) + strlen(dimName) + 2;
    char *newName = (char *) malloc(sizeof(char) * len);
    strcpy(newName, dimName);
    strcat(newName, "_");
    strcat(newName, name);
    return newName;
}

// lookup a dimensions real name
static FlexpathAltName *find_alt_name(FlexpathFMStructure *currentFm, char *dimName, char *varName) {
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
int get_local_offsets(struct adios_var_struct * v, struct adios_group_struct * g, int** offsets, int** dimensions)
{
    struct adios_dimension_struct * dim_list = v->dimensions;	    

    int ndims = 0;
    while (dim_list) {
        ndims++;
        dim_list = dim_list->next;
    }
    dim_list = v->dimensions;	    

    if(ndims){		
        int * local_offsets = (int*) malloc(sizeof(int) * (ndims));
        int * local_dimensions = (int*) malloc(sizeof(int) * (ndims));
        int n = 0; 
        while(dim_list) {		
            local_dimensions[n] = (int) adios_get_dim_value (&dim_list->dimension);
            local_offsets[n] = (int) adios_get_dim_value (&dim_list->local_offset);
            dim_list=dim_list->next;
            n++;
        }
        *offsets = local_offsets;	   
        *dimensions = local_dimensions;
    } else {
        *offsets = NULL;
        *dimensions = NULL;
    }
    return ndims;
}

// creates multiqueue function to handle ctrl messages for given bridge stones 
char *multiqueue_action = "{\n\
    int found = 0;\n\
    int flush_data_count = 0; \n\
    int my_rank = -1;\n\
    attr_list mine;\n\
    if(EVcount_varMsg()>0) {\n\
        EVdiscard_and_submit_varMsg(0, 0);\n\
    }\n\
    if(EVcount_drop_evgroup_msg()>0) {\n\
       if(EVcount_evgroup()>0){\n\
          EVdiscard_evgroup(0);\n\
       }\n\
       EVdiscard_and_submit_drop_evgroup_msg(0,0);\n\
    }\n\
    if(EVcount_update_step_msg()>0) {\n\
        mine = EVget_attrs_update_step_msg(0);\n\
        found = attr_ivalue(mine, \"fp_dst_rank\");\n\
        if(found > 0) {\n\
            EVdiscard_and_submit_update_step_msg(found, 0);\n\
        }\n\
    }\n\
    if(EVcount_op_msg()>0) {\n\
        mine = EVget_attrs_op_msg(0);\n\
        found = attr_ivalue(mine, \"fp_dst_rank\");\n\
        if(found > 0) {\n\
            EVdiscard_and_submit_op_msg(found, 0);\n\
        } else {\n\
            EVdiscard_and_submit_op_msg(0,0);\n\
        }\n\
    }\n\
    if(EVcount_flush()>0) {\n\
        flush* c = EVdata_flush(0);\n\
         if(c->type == 2){ \n\
             if(EVcount_evgroup()>0){\n\
               evgroup *g = EVdata_evgroup(0); \n\
               g->condition = c->condition;\n\
               EVsubmit(c->rank+1, g);\n\
               EVdiscard_flush(0);\n\
             }\n\
        } else {\n\
            EVdiscard_and_submit_flush(0,0);\n\
            flush_data_count++;\n\
        }\n\
    }\n\
    if(EVcount_anonymous()>0){\n\
        mine = EVget_attrs_anonymous(0);\n\
        found = attr_ivalue(mine, \"fp_dst_rank\");\n\
        EVdiscard_and_submit_anonymous(found+1,0);\n\
    }\n\
 }";

// sets a field based on data type
void set_field(int type, FMFieldList* field_list_ptr, int fieldNo, int* size){
    FMFieldList field_list = *field_list_ptr;
    switch (type) {
    case adios_unknown:
	perr("set_field: Bad Type Error\n");
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
	perr("set_field: Unknown Type Error\n");
	break;
    }
    *field_list_ptr = field_list;
}

// find a field in a given field list
static FMField *internal_find_field(char *name, FMFieldList flist) {
    FMField *f = flist;
    while (f->field_name != NULL && strcmp(f->field_name, name)) {
	f++;
    }
    return f;
}

// generic memory check for after mallocs
void mem_check(void* ptr, const char* str) {
    if(!ptr) {
        adios_error(err_no_memory, "Cannot allocate memory for flexpath %s.", str);
    }
}


static char * get_dim_name (struct adios_dimension_item_struct *d)
{
    char *vname = NULL;
    if (d->var) {
        vname = d->var->name;
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
FlexpathFMStructure* set_format(struct adios_group_struct* t, struct adios_var_struct* fields, FlexpathWriteFileData* fileData)
{
    FMStructDescRec *format = malloc(sizeof(FMStructDescRec)*2);
    mem_check(format, "format");
    memset(format, 0, sizeof(FMStructDescRec)*2);
    
    FlexpathFMStructure *currentFm = malloc(sizeof(FlexpathFMStructure));
    mem_check(currentFm, "currentFm");
    memset(currentFm, 0, sizeof(FlexpathFMStructure));

    LIST_INIT(&currentFm->nameList);
    LIST_INIT(&currentFm->dimList);
    currentFm->format = format;
    format->format_name = strdup(t->name);

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
    struct adios_var_struct *f;
    for (f = t->vars; f != NULL; f = f->next, fieldNo++) {
	char *tempName = get_fixed_name(f->name);
	if (strcmp(tempName, f->name)) {
	    FlexpathNameTable *nameNode = (FlexpathNameTable *) malloc(sizeof(FlexpathNameTable));
	    nameNode->originalName = strdup(f->name);
	    nameNode->mangledName = strdup(tempName);
	    LIST_INSERT_HEAD(&currentFm->nameList, nameNode, entries);
	}

	// use the mangled name for the field.
	field_list[fieldNo].field_name = tempName;
        if(tempName!=NULL) {
            int num_dims = 0;
            char atom_name[200] = "";
            FlexpathVarNode* dims=NULL;
            if(f->dimensions) {
                struct adios_dimension_struct* adim = f->dimensions;  
	
                // attach appropriate attrs for dimensions	
                for(; adim != NULL; adim = adim->next) {
                    num_dims++;		    
                    
                    char *vname = get_dim_name(&adim->dimension);
                    if (vname) {
			char *name = find_fixed_name(currentFm, vname);
			char *aname = get_alt_name(tempName,  name);
			dims=add_var(dims, strdup(aname), NULL, 0);
			set_attr_dimensions(tempName, aname, num_dims, fileData->attrs);
		    }
                    char *gname = get_dim_name(&adim->global_dimension);
		    if(gname) {
			fileData->globalCount++;
			char *name = find_fixed_name(currentFm, gname);
			char *aname = get_alt_name(tempName, name);
			dims=add_var(dims, strdup(aname), NULL, 0);
			set_attr_dimensions(tempName, aname, num_dims, fileData->attrs);
		    }
                }
            }
            // attach ndims attr
            strcat(atom_name, tempName);
            strcat(atom_name, "_");
            strcat(atom_name, FP_NDIMS_ATTR_NAME);
            atom_t ndims_atom = attr_atom_from_string(strdup(atom_name));
            add_int_attr(fileData->attrs, ndims_atom, num_dims);
            fileData->formatVars = add_var(fileData->formatVars, tempName, dims, 0);
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
		    char *name = find_fixed_name(currentFm, vname);
		    FlexpathAltName *a = find_alt_name(currentFm, name, (char*)field_list[fieldNo].field_name);
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
		  
	    while(currentFm->size % 8 != 0) {
		currentFm->size ++;					
	    }
		  
	    switch (f->type) {
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

	    case adios_unsigned_long:
		field_list[fieldNo].field_type =
		    (char *) malloc(sizeof(char) * 255);
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

	    default:
		adios_error(err_invalid_group, "set_format: Unknown Type Error %d: name: %s\n", f->type, field_list[fieldNo].field_name);
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

    format->field_list = field_list;
    currentFm->format->struct_size = currentFm->size;

    currentFm->buffer = malloc(currentFm->size);
    memset(currentFm->buffer, 0, currentFm->size);

    return currentFm;
}

// copies buffer zeroing out arrays that havent been asked for
void* copy_buffer(void* buffer, int rank, FlexpathWriteFileData* fileData){
    char* temp = (char*)malloc(fileData->fm->size);
    memcpy(temp, buffer, fileData->fm->size);
    FMField *f = fileData->fm->format->field_list;
    while (f->field_name != NULL)
    {
        FlexpathVarNode* a;
        if(!queue_contains(fileData->askedVars, f->field_name, rank)) {
            if((a=queue_contains(fileData->formatVars, f->field_name, -1)) 
	       && 
	       (a->dimensions != NULL)) {
                FlexpathVarNode* dim = a->dimensions;
                while(dim) {
                    FMField *f2 = fileData->fm->format->field_list;
                    while(f2->field_name != NULL) {
                        if(strcmp(f2->field_name, dim->varName)==0) {
                            break;
                        }
                        f2++;
                    }
                    if(f2->field_name != NULL) {
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

// terminal action for var messages: enqueues
static int var_handler(CManager cm, void *vevent, void *client_data, attr_list attrs){
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*) client_data;
    Var_msg* msg = (Var_msg*) vevent;
    EVtake_event_buffer(cm, msg);
    fp_write_log("MSG", "recieved var_msg : rank %d\n", msg->rank);
    threaded_enqueue(&fileData->controlQueue, msg, VAR, 
        &fileData->controlMutex, &fileData->controlCondition);
    return 0;
}

// terminal action for flush messages: enqueues
static int 
flush_handler(CManager cm, void* vevent, void* client_data, attr_list attrs) {
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*) client_data;
    Flush_msg* msg = (Flush_msg*) vevent;
    EVtake_event_buffer(cm, msg);
    fp_write_log("MSG", "recieved flush : rank %d type data\n", msg->rank);
    threaded_enqueue(&fileData->controlQueue, msg, DATA_FLUSH, 
        &fileData->controlMutex, &fileData->controlCondition);
    return 0;
}

static int
drop_evgroup_handler(CManager cm, void *vevent, void *client_data, attr_list attrs){
    drop_evgroup_msg *msg = vevent;
    CMCondition_signal(cm, msg->condition);    
    return 0;
}

// terminal action for op messages: enqueues
static int 
op_handler(CManager cm, void* vevent, void* client_data, attr_list attrs) {
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*) client_data;
    op_msg* msg = (op_msg*) vevent;
    EVtake_event_buffer(cm, msg);
    fp_write_log("MSG", "recieved op_msg : rank %d type %d: condition: %d step: %d\n", 
		 msg->process_id, msg->type, msg->condition, msg->step);
    if(msg->type == OPEN_MSG) {
        threaded_enqueue(&fileData->controlQueue, msg, OPEN, 
            &fileData->controlMutex, &fileData->controlCondition);
    } else if(msg->type == CLOSE_MSG) {
        threaded_enqueue(&fileData->controlQueue, msg, CLOSE, 
			 &fileData->controlMutex, &fileData->controlCondition);
    } else if(msg->type == INIT_MSG) {
	threaded_enqueue(&fileData->controlQueue, msg, INIT,
			 &fileData->controlMutex, &fileData->controlCondition);
			
    }
    return 0;
}

// sets a size atom
attr_list set_size_atom(attr_list attrs, int value) {
    atom_t dst_atom = attr_atom_from_string("fp_size");
    int size;
    if(!get_int_attr(attrs, dst_atom, &size)) {
        add_int_attr(attrs, dst_atom, value);
    }
    set_int_attr(attrs, dst_atom, value);
    return attrs;
}

// sets a dst rank atom
attr_list 
set_dst_rank_atom(attr_list attrs, int value) {
    atom_t dst_atom = attr_atom_from_string("fp_dst_rank");
    int dst;
    if(!get_int_attr(attrs, dst_atom, &dst)) {
        add_int_attr(attrs, dst_atom, value);
    }
    set_int_attr(attrs, dst_atom, value);
    return attrs;
}

// sets a dst condition atom
attr_list 
set_dst_condition_atom(attr_list attrs, int condition){
    atom_t dst_atom = attr_atom_from_string("fp_dst_condition");
    int dst;
    if(!get_int_attr(attrs, dst_atom, &dst)){
	add_int_attr(attrs, dst_atom, condition);
    }
    set_int_attr(attrs, dst_atom, condition);
    return attrs;
}

void
send_update_step_msgs(FlexpathWriteFileData *fileData, int step)
{
    int i;
    for(i = 0; i<fileData->num_reader_coordinators; i++){
	update_step_msg msg;
	msg.process_id = fileData->rank;
	msg.step = step;
	msg.finalized = fileData->finalized;
	int dest_rank = fileData->reader_coordinators[i];
	fileData->attrs = set_dst_rank_atom(fileData->attrs, dest_rank+1);
	EVsubmit(fileData->stepSource, &msg, fileData->attrs);
    }
}

// processes messages from control queue
void 
control_thread(void* arg) 
{
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*)arg;
    int rank = fileData->rank;
    FlexpathQueueNode* controlMsg;
    FlexpathQueueNode* dataNode;
    while(1) {
        fp_write_log("CONTROL", "control message attempts dequeue\n");
	if((controlMsg = threaded_dequeue(&fileData->controlQueue, 
	    &fileData->controlMutex, &fileData->controlCondition, NULL, 0))) {
            fp_write_log("CONTROL", "control message dequeued\n");
	    if(controlMsg->type==VAR) {
		Var_msg* varMsg = (Var_msg*) controlMsg->data;
		fileData->askedVars = add_var(fileData->askedVars, 
		    strdup(varMsg->var_name), NULL, varMsg->rank);
		EVreturn_event_buffer(flexpathWriteData.cm,controlMsg->data);
	    } else if(controlMsg->type==DATA_FLUSH) {
                fp_write_log("DATAMUTEX", "in use 1\n"); 
		dataNode = threaded_peek(&fileData->dataQueue, 
		    &fileData->dataMutex, &fileData->dataCondition);
                fp_write_log("DATAMUTEX", "no use 1\n"); 
		Flush_msg* flushMsg = (Flush_msg*) controlMsg->data;
		fp_write_log("QUEUE", "dataNode:%p, flushMsg:%p\n", dataNode, flushMsg);
                void* temp = copy_buffer(dataNode->data, flushMsg->rank, fileData);
		fileData->attrs = set_dst_rank_atom(fileData->attrs, flushMsg->rank);
		fileData->attrs = set_dst_condition_atom(fileData->attrs, flushMsg->condition);
		if(!fileData->bridges[flushMsg->rank].opened) {
                  fileData->bridges[flushMsg->rank].opened=1;
                  fileData->openCount++;
                }
		EVsubmit_general(fileData->dataSource, temp, data_free, fileData->attrs);
	    } else if(controlMsg->type==OPEN) {
                op_msg* open = (op_msg*) controlMsg->data;
                fileData->bridges[open->process_id].step = open->step;
                fileData->bridges[open->process_id].condition = open->condition;
		if(!fileData->bridges[open->process_id].created){
		    fileData->bridges[open->process_id].myNum = 
			EVcreate_bridge_action(flexpathWriteData.cm, 
					       attr_list_from_string(fileData->bridges[open->process_id].contact), 
					       fileData->bridges[open->process_id].theirNum);
		    
		    EVaction_set_output(flexpathWriteData.cm, 
					fileData->multiStone, 
					fileData->multi_action, 
					open->process_id+1, 
					fileData->bridges[open->process_id].myNum);				    
		}		
		if(open->step < fileData->currentStep) {
		    log_error("Flexpath method control_thread: Received Past Step Open\n");
                } else if (open->step == fileData->currentStep){
                    pthread_mutex_lock(&fileData->openMutex);
                    fileData->openCount++;  
                    fileData->bridges[open->process_id].opened = 1;
		    pthread_mutex_unlock(&fileData->openMutex);
                    op_msg* ack = malloc(sizeof(op_msg));
                    ack->file_name = strdup(fileData->name);
                    ack->process_id = fileData->rank;
                    ack->step = fileData->currentStep;
                    ack->type = 2;
		    ack->condition = open->condition;
                    fileData->attrs = set_dst_rank_atom(fileData->attrs, open->process_id+1);
                    EVsubmit_general(fileData->opSource, ack, op_free, fileData->attrs);
                } else {
                    fp_write_log("STEP", "recieved op with future step\n");
                }
		EVreturn_event_buffer(flexpathWriteData.cm, open);
            } else if(controlMsg->type==CLOSE) {
                op_msg* close = (op_msg*) controlMsg->data;
		pthread_mutex_lock(&fileData->openMutex);
                fp_write_log("MUTEX","lock 7\n");
		fileData->openCount--;
                fileData->bridges[close->process_id].opened=0;
                fp_write_log("MUTEX","unlock 7\n");
		pthread_mutex_unlock(&fileData->openMutex);
                 if(fileData->openCount==0) {
                    fp_write_log("STEP", "advancing\n");
                    fp_write_log("DATAMUTEX", "in use 2\n"); 
		    FlexpathQueueNode* node = threaded_dequeue(&fileData->dataQueue, 
		        &fileData->dataMutex, &fileData->dataCondition, &fileData->dataCondition2, 1);
                    fp_write_log("DATAMUTEX", "no use 2\n"); 
                    int q = queue_count(&fileData->dataQueue);
                    fp_write_log("QUEUE", "after step queue count now %d\n", q);
                    FMfree_var_rec_elements(fileData->fm->ioFormat, node->data);

		    drop_evgroup_msg dropMsg;
		    dropMsg.step = fileData->currentStep;
		    dropMsg.condition = CMCondition_get(flexpathWriteData.cm, NULL);
		    EVsubmit(fileData->dropSource, &dropMsg, fileData->attrs);
		    CMCondition_wait(flexpathWriteData.cm,  dropMsg.condition); 
                    fileData->currentStep++;                    
                    int i;
                    //for all bridges if step == currentstep send ack
		    // this block gets repeated in finalize.  gets repeated
		    // only AFTER sending finalize messages.  
		    // do it for everyone that has opened.
		    
                    for(i=0; i<fileData->numBridges; i++) {
                      if(fileData->bridges[i].step==fileData->currentStep) {
                        fileData->openCount++;
                        fileData->bridges[i].opened = 1;
                        op_msg* ack = malloc(sizeof(op_msg));
                        ack->file_name = strdup(fileData->name);
                        ack->process_id = fileData->rank;
                        ack->step = fileData->currentStep;
                        ack->type = 2;
			ack->condition = fileData->bridges[i].condition;
                        fileData->attrs = set_dst_rank_atom(fileData->attrs, i+1);
                        EVsubmit_general(fileData->opSource, 
					 ack, 
					 op_free, 
					 fileData->attrs);
                      }
                    }
		 }
		 EVreturn_event_buffer(flexpathWriteData.cm, close);
	    }else if(controlMsg->type == INIT){ 
		fp_write_log("DATAMUTEX", "in use 1\n"); 
		dataNode = threaded_peek(&fileData->dataQueue, 
		    &fileData->dataMutex, &fileData->dataCondition);
		op_msg* initMsg = (op_msg*) controlMsg->data;
		fileData->num_reader_coordinators++;
		fileData->reader_coordinators = realloc(fileData->reader_coordinators,
							sizeof(int)*fileData->num_reader_coordinators);
		fileData->reader_coordinators[fileData->num_reader_coordinators -1] = initMsg->process_id;
		void* temp = copy_buffer(dataNode->data, 
					 initMsg->process_id, fileData);
		fileData->attrs = set_dst_rank_atom(fileData->attrs, 
						    initMsg->process_id);
		fileData->attrs = set_dst_condition_atom(fileData->attrs, 
							 initMsg->condition);	
		fileData->bridges[initMsg->process_id].created = 1;
		fileData->bridges[initMsg->process_id].myNum = 
		    EVcreate_bridge_action(flexpathWriteData.cm, 
					   attr_list_from_string(fileData->bridges[initMsg->process_id].contact), 
					   fileData->bridges[initMsg->process_id].theirNum);
		    
		EVaction_set_output(flexpathWriteData.cm, 
				    fileData->multiStone, 
				    fileData->multi_action, 
				    initMsg->process_id+1, 
				    fileData->bridges[initMsg->process_id].myNum);

		EVsubmit_general(fileData->dataSource, 
				 temp, data_free, fileData->attrs);
		EVreturn_event_buffer(flexpathWriteData.cm, initMsg);
	    }
	    else{
		log_error("control_thread: Unrecognized Control Message\n");
	    }
	}
    }
    return;
}

// adds an open file handle to global open file list
void add_open_file(FlexpathWriteFileData* newFile) {
    FlexpathWriteFileData* last = flexpathWriteData.openFiles;
    while(last && last->next) {
        last = last->next;
    }
    if(last) {
        last->next = newFile;
    } else {
        flexpathWriteData.openFiles = newFile;
    }
}

// searches for an open file handle
FlexpathWriteFileData* find_open_file(char* name) {
    FlexpathWriteFileData* file = flexpathWriteData.openFiles;
    while(file && strcmp(file->name, name)) {
        file = file->next;
    }
    return file;
}


// Initializes flexpath write local data structures
extern void adios_flexpath_init(const PairStruct *params, struct adios_method_struct *method) 
{
    setenv("CMSelfFormats", "1", 1);
    // global data structure creation
    flexpathWriteData.rank = -1;
    flexpathWriteData.openFiles = NULL;
    
    // setup CM
    flexpathWriteData.cm = CManager_create();
    atom_t CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    char * transport = getenv("CMTransport");
    if(transport == NULL){
	fp_write_log("SETUP","transport is null\n");
	if(CMlisten(flexpathWriteData.cm) == 0) {
	    perr( "error: unable to initialize connection manager.\n");
	    exit(1);
	}
    } else {
	fp_write_log("SETUP", "writer transport: %s\n", transport);
	attr_list listen_list = create_attr_list();
	add_attr(listen_list, CM_TRANSPORT, Attr_String, (attr_value)strdup(transport));
	CMlisten_specific(flexpathWriteData.cm, listen_list);
    }
    
    // configuration setup
    //gen_pthread_init();
    setenv("CMSelfFormats", "1", 1);
    
    // fork communications thread
    int forked = CMfork_comm_thread(flexpathWriteData.cm);   
    if(!forked) {
         perr( "error forking comm thread\n");
    }
}

// opens a new adios file for writes
extern int 
adios_flexpath_open(struct adios_file_struct *fd, struct adios_method_struct *method, MPI_Comm comm) 
{    
    if( fd == NULL || method == NULL) {
        perr("open: Bad input parameters\n");
        return -1;
    }

    // file creation
    if(find_open_file(method->group->name)) {
        // stream already open
        return 0;
    }

    FlexpathWriteFileData *fileData = malloc(sizeof(FlexpathWriteFileData));
    mem_check(fileData, "fileData");
    memset(fileData, 0, sizeof(FlexpathWriteFileData));
    
    fileData->maxQueueSize=0;
    if(method->parameters) {
        sscanf(method->parameters,"QUEUE_SIZE=%d;",&fileData->maxQueueSize);
    }
    
    // setup step state
    fileData->attrs = create_attr_list();
    fileData->openCount = 0;
    //fileData->currentStep = 0;

    // setup mutexs
    pthread_mutex_t ctrlm = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t dm =  PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t dm2 =  PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t om =  PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t evm = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t cc =  PTHREAD_COND_INITIALIZER;
    pthread_cond_t dc =  PTHREAD_COND_INITIALIZER;
    pthread_cond_t dc2 =  PTHREAD_COND_INITIALIZER;
    pthread_cond_t evc = PTHREAD_COND_INITIALIZER;

    fileData->controlMutex = ctrlm;
    fileData->dataMutex = dm;
    fileData->dataMutex2 = dm2;
    fileData->openMutex = om;
    fileData->evgroupMutex = evm;
     
    // setup conditions
    fileData->controlCondition = cc;
    fileData->dataCondition = dc;
    fileData->dataCondition2 = dc2;
    fileData->evgroupCondition = evc;
    // communication channel setup
    char writer_info_filename[200];
    char writer_ready_filename[200];
    char reader_info_filename[200];
    char reader_ready_filename[200];
    
    /*
    // Titan filesystem specific
    char * filebase = "/tmp/work/jdayal3/titan/";
    sprintf(writer_info_filename, "%s", filebase);
    sprintf(writer_ready_filename, "%s", filebase);
    sprintf(reader_info_filename, "%s", filebase);
    sprintf(reader_ready_filename, "%s", filebase);
    */

    int i=0;
    flexpathWriteData.rank = fileData->rank;
    fileData->globalCount = 0;
    fileData->sentGlobalOffsets = 0;

    // mpi setup
    MPI_Comm_dup(comm, &fileData->mpiComm);

    MPI_Comm_rank((fileData->mpiComm), &fileData->rank);
    MPI_Comm_size((fileData->mpiComm), &fileData->size);
    char *recv_buff = NULL;
    char sendmsg[CONTACT_STR_LEN];
    if(fileData->rank == 0) {
        recv_buff = (char *) malloc(fileData->size*CONTACT_STR_LEN*sizeof(char));
    }
        
    // send out contact string
    char * contact = attr_list_to_string(CMget_contact_list(flexpathWriteData.cm));
    fileData->multiStone = EValloc_stone(flexpathWriteData.cm);
    fileData->sinkStone = EValloc_stone(flexpathWriteData.cm);
    sprintf(&sendmsg[0], "%d:%s", fileData->multiStone, contact);
    MPI_Gather(sendmsg, CONTACT_STR_LEN, MPI_CHAR, recv_buff, 
        CONTACT_STR_LEN, MPI_CHAR, 0, (fileData->mpiComm));

    // rank 0 prints contact info to file
    if(fileData->rank == 0) {
        sprintf(writer_info_filename, "%s_%s", fd->name, "writer_info.txt");
        FILE* writer_info = fopen(writer_info_filename,"w");
        for(i=0; i<fileData->size; i++) {
            fprintf(writer_info, "%s\n",&recv_buff[i*CONTACT_STR_LEN]); 
        }
        fclose(writer_info);
    }

    // poll file - race condition issues
    FILE* reader_ready = NULL;
    sprintf(reader_ready_filename, "%s_%s", fd->name, "reader_ready.txt");
    while(!reader_ready){
	reader_ready = fopen(reader_ready_filename,"r");
    }
    fclose(reader_ready);

    // read contact list
    sprintf(reader_info_filename, "%s_%s", fd->name, "reader_info.txt");
    FILE* reader_info = fopen(reader_info_filename, "r");
    if(!reader_info){
	reader_info = fopen(reader_info_filename, "r");
    }
    char in_contact[CONTACT_STR_LEN] = "";
    int numBridges = 0;
    int stone_num;
    // build a bridge per line
    while(fscanf(reader_info, "%d:%s",&stone_num, in_contact)!=EOF){
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
    if(fileData->rank == 0){
	unlink(reader_info_filename);
	unlink(reader_ready_filename);
    }
	
    //process group format
    struct adios_group_struct *t = method->group;

    if(t == NULL){
	adios_error(err_invalid_group, "Invalid group.\n");
	return err_invalid_group;
    }
    struct adios_var_struct *fields = t->vars;
	
    if(fields == NULL){
	adios_error(err_invalid_group, "Group has no variables.\n");
	return err_invalid_group;
    }	

    fileData->fm = set_format(t, fields, fileData);
    fp_write_log("SETUP", "set format complete\n");

    // attach rank attr and add file to open list
    fp_write_log("FILE", "opening file %s\n", method->group->name);
    fileData->name = strdup(method->group->name); 
    add_open_file(fileData);
    atom_t rank_atom = attr_atom_from_string(FP_RANK_ATTR_NAME);
    add_int_attr(fileData->attrs, rank_atom, fileData->rank);

    //generate multiqueue function that sends formats or all data based on flush msg
    fp_write_log("SETUP", "setup graph\n");
    FMStructDescList queue_list[] = {flush_format_list, 
				     format_format_list, 
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
						  fileData->multiStone, q_action_spec, NULL);
    fileData->formatSource = EVcreate_submit_handle(flexpathWriteData.cm, 
						    fileData->multiStone, format_format_list);
    fileData->dataSource = EVcreate_submit_handle_free(flexpathWriteData.cm, 
						       fileData->multiStone, fileData->fm->format, data_free,  NULL); 
    fileData->opSource = EVcreate_submit_handle_free(flexpathWriteData.cm, 
						     fileData->multiStone, op_format_list, op_free,  NULL); 
    fileData->offsetSource = EVcreate_submit_handle(flexpathWriteData.cm, 
						    fileData->multiStone, evgroup_format_list);
    fileData->dropSource = EVcreate_submit_handle(flexpathWriteData.cm, 
						  fileData->multiStone, drop_evgroup_msg_format_list);
    fileData->stepSource = EVcreate_submit_handle(flexpathWriteData.cm, 
						    fileData->multiStone, update_step_msg_format_list);
    
    
    fp_write_log("SETUP", "setup terminal actions\n");
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
			    var_format_list, var_handler, fileData);
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
			    op_format_list, op_handler, fileData);
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
			    drop_evgroup_msg_format_list, drop_evgroup_handler, fileData);
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
	flush_format_list, flush_handler, fileData);

    //link multiqueue to sink
    fp_write_log("SETUP", "linking stones\n");
    EVaction_set_output(flexpathWriteData.cm, fileData->multiStone, 
        fileData->multi_action, 0, fileData->sinkStone);
   
    fp_write_log("SETUP", "arranged evpath graph\n");
	
    FMContext my_context = create_local_FMcontext();
    fileData->fm->ioFormat = register_data_format(my_context, fileData->fm->format);
    
    fp_write_log("SETUP", "indicating to reader that ready\n");
    sprintf(writer_ready_filename, "%s_%s", fd->name, "writer_ready.txt");
    if(fileData->rank == 0) {
        FILE* writer_info = fopen(writer_ready_filename, "w");
        fprintf(writer_info, "ready");
        fclose(writer_info);
    }
        
    fp_write_log("SETUP", "fork control thread\n");
    
    pthread_create(&fileData->ctrl_thr_id, NULL, (void*)&control_thread, fileData);   
    return 0;	
}




//  writes data to multiqueue
extern void adios_flexpath_write(
    struct adios_file_struct *fd, 
    struct adios_var_struct *f, 
    void *data, 
    struct adios_method_struct *method) 
{
    fp_write_log("FILE", "entering flexpath file %s write\n", method->group->name);
    FlexpathWriteFileData* fileData = find_open_file(method->group->name);
    FlexpathFMStructure* fm = fileData->fm;

    if (fm == NULL)
    {
	log_error("adios_flexpath_write: something has gone wrong with format registration: %s\n", f->name);
	return;
    }
    
    FMFieldList flist = fm->format->field_list;
    FMField *field = NULL;
    char *fixedname = find_fixed_name(fm, f->name);
    field = internal_find_field(fixedname, flist);
    if (field != NULL) {
	//scalar quantity
	if (!f->dimensions) {
	    if (data) {
		//why wouldn't it have data?
		memcpy(&fm->buffer[field->field_offset], data, field->field_size);

		//scalar quantities can have FlexpathAltNames also so assign those
		if(field->field_name != NULL) {
					
		    FlexpathDimNames *d = NULL;
		    for (d = fm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
			if (!strcmp(d->name, field->field_name)) {
			    //matches
			    //check if there are FlexpathAltNames
			    FlexpathAltName *a = NULL;
			    for (a = d->altList.lh_first; a != NULL; a = a->entries.le_next) {
				//use the FlexpathAltName field to get the data into the buffer
				memcpy(&fm->buffer[a->field->field_offset], 
				       data, 
				       a->field->field_size);
                		//int *testingint = (int*)&fm->buffer[a->field->field_offset];
		        	//perr( "writing %s to %s at %d %d\n", f->name, a->name, a->field->field_offset, (int)*testingint);
			    }
			}
		    }
		}
	    } else {
		log_error("adios_flexpath_write: something has gone wrong with variable creation: %s\n", f->name);
	    }
	} else {
	    //vector quantity
	    if (data)
	    {
                //perr( "copying vector pointer\n");
		//we just need to copy the pointer stored in f->data
                // calculate size
                memcpy(&fm->buffer[field->field_offset], &data, sizeof(void *));

	    } else {
		log_error("adios_flexpath_write: no array data found for var: %s. Bad.\n", f->name);
		//perr( "no data for vector %s\n", f->name);
	    }
	}
    }
    //perr( "successfully copied data to buffer\n");
}

extern void 
adios_flexpath_close(struct adios_file_struct *fd, struct adios_method_struct *method) 
{
    fp_write_log("FILE", "file close %s\n", method->group->name);
    FlexpathWriteFileData* fileData = find_open_file(method->group->name);
    void* buffer = malloc(fileData->fm->size);

    struct adios_group_struct * g2 = fd->group;
    struct adios_var_struct * fields = g2->vars;
    while(fields) {       
        if(fields->dimensions) {
            struct adios_dimension_struct* dims = fields->dimensions;
            int total_size = 1;
            //for each dimension
            while(dims) {    
                int size = adios_get_dim_value (&dims->dimension);
                total_size *= size;
                dims = dims->next;
            }		
            FMFieldList flist = fileData->fm->format->field_list;
            FMField *field = NULL;
            char *fixedname = find_fixed_name(fileData->fm, fields->name);
            field = internal_find_field(fixedname, flist);
            //perr( "field offset %d size %d\n", field->field_offset, field->field_size);

            total_size*=field->field_size;
            // malloc size
            void* pointer_data_copy = malloc(total_size);
            // while null
            while(pointer_data_copy==NULL) { 
                perr("mallocing space for user buffer failed, trying again soon\n");
                sleep(1);
                void* pointer_data_copy = malloc(total_size);
                //block
            }
                
            fp_write_log("DATA","Attempting to get pointer to user data\n");
            void* temp = get_FMPtrField_by_name(flist, fields->name, fileData->fm->buffer, 0);
            fp_write_log("DATA","Copying user data to new space\n");
            memcpy(pointer_data_copy, temp, total_size);
            fp_write_log("DATA","Setting pointer to new space\n");
            set_FMPtrField_by_name(flist, fields->name, fileData->fm->buffer, pointer_data_copy);
        }    
        fields = fields->next;
    }

    
    memcpy(buffer, fileData->fm->buffer, fileData->fm->size);

    fp_write_log("DATAMUTEX", "in use 3\n"); 
    threaded_enqueue(&fileData->dataQueue, buffer, 
		     DATA_BUFFER,
		     &fileData->dataMutex, 
		     &fileData->dataCondition);
    
    int c = 0;
 
    // now gather offsets and send them via MPI to root
    struct adios_group_struct * g = fd->group;
    struct adios_var_struct * list = g->vars;

    if(fileData->globalCount > 0){	
	fp_write_log("BOUNDING", "check offsets\n");
        // process local offsets here	
	int num_gbl_vars = 0;
        global_var * gbl_vars = NULL;
	int num_vars = 0;
	int myrank = fileData->rank;
	int commsize = fileData->size;

	while(list){
	    //int num_local_offsets = 0;
	    int * local_offsets = NULL;
	    int * local_dimensions = NULL;
	    int num_local_offsets = get_local_offsets(list, g, &local_offsets, &local_dimensions);
	    
	    if(num_local_offsets > 0){
		int * all_offsets = NULL;
		int * all_local_dims = NULL;
		
		int buf_size = num_local_offsets * commsize * sizeof(int);		    
		all_offsets = (int*)malloc(buf_size);		
		all_local_dims = (int*)malloc(buf_size);
		

		MPI_Allgather(local_offsets, num_local_offsets, MPI_INT, 
			      all_offsets, num_local_offsets, MPI_INT,
			      fileData->mpiComm);

		MPI_Allgather(local_dimensions, num_local_offsets, MPI_INT, 
			      all_local_dims, num_local_offsets, MPI_INT,
			      fileData->mpiComm);

		
		num_gbl_vars++;
		offset_struct * ostruct = (offset_struct*)malloc(sizeof(offset_struct));
		ostruct->offsets_per_rank = num_local_offsets;
		ostruct->total_offsets = num_local_offsets * commsize;
		ostruct->local_offsets = all_offsets;
		ostruct->local_dimensions = all_local_dims;
		gbl_vars = realloc(gbl_vars, sizeof(global_var) * num_gbl_vars);
		gbl_vars[num_gbl_vars - 1].name = strdup(list->name);
		gbl_vars[num_gbl_vars - 1].noffset_structs = 1;
		gbl_vars[num_gbl_vars - 1].offsets = ostruct;

	    }
	    list=list->next;
            free (local_offsets);
            free (local_dimensions);
	}
	    
	evgroup * gp = malloc(sizeof(evgroup));
	gp->num_vars = num_gbl_vars;
	gp->step = fileData->writerStep;
	gp->vars = gbl_vars;
	fileData->gp = gp;       
	fileData->attrs = set_size_atom(fileData->attrs, fileData->size);
	EVsubmit(fileData->offsetSource, gp, fileData->attrs);
	fileData->sentGlobalOffsets = 1;
    }
    //send_update_step_msgs(fileData, fileData->writerStep);
    fileData->writerStep++;
    while((c=queue_count(&fileData->dataQueue))>fileData->maxQueueSize) {
        fp_write_log("QUEUE", "waiting for queue to be below max size\n");
        pthread_cond_wait(&fileData->dataCondition2, &fileData->dataMutex2);
        fp_write_log("QUEUE", "wakeup on queue size\n");
    }
    fp_write_log("FILE", "file close %s exiting\n", method->group->name);
}

// wait until all open files have finished sending data to shutdown
extern void adios_flexpath_finalize(int mype, struct adios_method_struct *method) 
{
    FlexpathWriteFileData* fileData = flexpathWriteData.openFiles;
    log_info("Flexpath method entered finalize: %d\n", fileData->rank);
    fp_write_log("FILE", "Entered finalize\n");
    while(fileData) {
        //fp_write_log("DATAMUTEX", "in use 4\n"); 
        //pthread_mutex_lock(fileData->dataMutex2);
        //fp_write_log("MUTEX","lock 1\n");
        while(fileData->dataQueue!=NULL) {
            fp_write_log("FILE", "waiting on %s to empty data\n", fileData->name);
	    pthread_cond_wait(&fileData->dataCondition2, &fileData->dataMutex2);
	}
	//fp_write_log("MUTEX","unlock 1\n");
	//pthread_mutex_unlock(fileData->dataMutex2);
	//fp_write_log("DATAMUTEX", "no use 4\n"); 
	fileData->finalized = 1;
	//send_update_step_msgs(fileData, fileData->writerStep);
	fileData = fileData->next;
	    
    }
}

// provides unknown functionality
extern enum ADIOS_FLAG adios_flexpath_should_buffer (struct adios_file_struct * fd,struct adios_method_struct * method) {
    fp_write_log("UNIMPLEMENTED", "adios_flexpath_should_buffer\n");
    return adios_flag_unknown;
}

// provides unknown functionality
extern void adios_flexpath_end_iteration(struct adios_method_struct *method) {
    fp_write_log("UNIMPLEMENTED", "adios_flexpath_end_iteration\n");
}

// provides unknown functionality
extern void adios_flexpath_start_calculation(struct adios_method_struct *method) {
    fp_write_log("UNIMPLEMENTED", "adios_flexpath_start_calculation\n");
}

// provides unknown functionality
extern void adios_flexpath_stop_calculation(struct adios_method_struct *method) {
    fp_write_log("UNIMPLEMENTED", "adios_flexpath_stop_calculation\n");
}

// provides unknown functionality
extern void adios_flexpath_get_write_buffer(struct adios_file_struct *fd,struct adios_var_struct *f, uint64_t *size, void **buffer, struct adios_method_struct *method) {
    fp_write_log("UNIMPLEMENTED", "adios_flexpath_get_write_buffer\n");
}

// should not be called from write, reason for inclusion here unknown
void adios_flexpath_read(struct adios_file_struct *fd, struct adios_var_struct *f, void *buffer, uint64_t buffer_size, struct adios_method_struct *method) {
    fp_write_log("UNIMPLEMENTED", "adios_flexpath_read\n");
}

#else // print empty version of all functions (if HAVE_FLEXPATH == 0)

void adios_flexpath_read(struct adios_file_struct *fd, struct adios_var_struct *f, void *buffer, struct adios_method_struct *method) {
}

extern void adios_flexpath_get_write_buffer(struct adios_file_struct *fd, struct adios_var_struct *f, unsigned long long *size, void **buffer, struct adios_method_struct *method) {
}

extern void adios_flexpath_stop_calculation(struct adios_method_struct *method) {
}

extern void adios_flexpath_start_calculation(struct adios_method_struct *method) {
}

extern void adios_flexpath_end_iteration(struct adios_method_struct *method) {
}

extern void adios_flexpath_finalize(int mype, struct adios_method_struct *method) {
}

extern void adios_flexpath_close(struct adios_file_struct *fd, struct adios_method_struct *method) {
}

extern void adios_flexpath_write(struct adios_file_struct *fd, struct adios_var_struct *f, void *data, struct adios_method_struct *method) {
}

extern void adios_flexpath_open(struct adios_file_struct *fd, struct adios_method_struct *method) {
}

extern void adios_flexpath_init(const PairStruct *params, struct adios_method_struct *method) {
}

enum ADIOS_FLAG adios_flexpath_should_buffer (struct adios_file_struct * fd, struct adios_method_struct * method) {
}

#endif
