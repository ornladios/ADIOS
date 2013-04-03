/*
    adios_flexpath.c
    uses evpath for io in conjunction with read/read_flexpath.c
*/


#if NO_FLEXPATH == 0

// system libraries
#include <stdio.h>
#include <stdlib.h>
#include <sys/queue.h>

// evpath libraries
#include <evpath.h>
#include <gen_thread.h>

// local libraries
#include "config.h"
#include "core/adios_internals.h"
#include "core/adios_transport_hooks.h"
#include "core/util.h"
#include "public/adios.h"
#include "public/flexpath.h"

/************************* Structure and Type Definitions ***********************/
// used for messages in the control queue
typedef enum {VAR=0, DATA_FLUSH, OPEN, CLOSE, DATA_BUFFER, OFFSET_MSG} FlexpathMessageType;

// maintains connection information
typedef struct _flexpath_stone {
    int myNum;
    int theirNum;
    int step;
    int opened;
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
    MPI_Comm * mpiComm;
    int rank;
    int size;

    // EVPath stuff
    EVstone multiStone;
    EVstone sinkStone;
    EVsource formatSource;
    EVsource dataSource;
    EVsource offsetSource;
    EVsource opSource;
    FlexpathStone* bridges;
    int numBridges;
    attr_list attrs;

    // server state
    int maxQueueSize;
    int openCount;
    int currentStep;
    thr_mutex_t openMutex;
    FlexpathFMStructure* fm;
    FlexpathVarNode* askedVars;
    FlexpathVarNode* writtenVars;
    FlexpathVarNode* formatVars;
    FlexpathQueueNode* controlQueue;
    FlexpathQueueNode* dataQueue;    
    thr_mutex_t controlMutex;
    thr_mutex_t dataMutex;
    thr_condition_t controlCondition;
    thr_condition_t dataCondition;

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

// checks for a valid mpi communicator 
static void adios_var_to_comm(const char* commName, enum ADIOS_FLAG hostLanguageFortran, void* data, MPI_Comm* comm) {
    fp_write_log("SETUP", "running var to comm function\n");
    if(data) {    
        int t = *(int*)data;
        if(!commName) {        
            if(!t) {            
                fprintf(stderr, "communicator not provided and none "
			 "listed in XML.  Defaulting to MPI_COMM_SELF\n");
                *comm = MPI_COMM_SELF;
            } else {            
                if(hostLanguageFortran == adios_flag_yes) {                
                    *comm = MPI_Comm_f2c(t);
                } else {                
                    *comm = *(MPI_Comm*)data;
                }
            }
        } else {        
            if(!strcmp(commName, "")) {            
                if(!t) {                
                    fprintf(stderr, "communicator not provided and none "
			     "listed in XML.  Defaulting to MPI_COMM_SELF\n");
                    *comm = MPI_COMM_SELF;
                } else {                
                    if(hostLanguageFortran == adios_flag_yes) {                    
                        *comm = MPI_Comm_f2c(t);
                    } else {                    
                        *comm = *(MPI_Comm*)data;
                    }
                }
            } else {            
                if(!t) {                
                    fprintf(stderr, "communicator not provided but one "
			     "listed in XML.  Defaulting to MPI_COMM_WORLD\n");
                    *comm = MPI_COMM_WORLD;
                } else {                
		    if(hostLanguageFortran == adios_flag_yes){                    
                        *comm = MPI_Comm_f2c(t);
                    } else {                   
                        *comm = *(MPI_Comm*)data;
                    }
                }
            }
        }
    } else {
        fprintf(stderr, "coordination-communication not provided. "
		 "Using MPI_COMM_WORLD instead\n");
        *comm = MPI_COMM_WORLD;
    }
}

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
    if(op->
    free(op);
}

// message queue add to head
void threaded_enqueue(FlexpathQueueNode** queue, void* item, FlexpathMessageType type, thr_mutex_t mutex, thr_condition_t condition) {
    fp_write_log("QUEUE", "enqueing a message\n");
    thr_mutex_lock(mutex);
    FlexpathQueueNode* newNode = (FlexpathQueueNode*) malloc(sizeof(FlexpathQueueNode));
    newNode->data = item;
    newNode->type = type;
    newNode->next = *queue;
    *queue = newNode;
    thr_condition_signal(condition);
    thr_mutex_unlock(mutex);
}

// message queue count
int queue_count(FlexpathQueueNode** queue, thr_mutex_t mutex) {
    fp_write_log("QUEUE", "counting a queue\n");
    thr_mutex_lock(mutex);
    if(*queue==NULL) {
        thr_mutex_unlock(mutex);
        return 0;
    }
    int count = 1;
    FlexpathQueueNode* current = *queue;
    while(current && current->next) {
        count++;
        current = current->next;
    }
    thr_mutex_unlock(mutex);
    return count;
}

// remove from tail of a message queue
FlexpathQueueNode* threaded_dequeue(FlexpathQueueNode** queue, thr_mutex_t mutex, thr_condition_t condition) {
    fp_write_log("QUEUE", "dequeue\n");
    thr_mutex_lock(mutex);
    fp_write_log("QUEUE", "should wait %d\n", gen_thr_initialized());
    while(*queue==NULL) {
        thr_condition_wait(condition, mutex);
    }
    fp_write_log("QUEUE", "but does not\n");
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
    thr_condition_signal(condition);
    thr_mutex_unlock(mutex);
    return tail;
}

// peek at tail of message queue
FlexpathQueueNode* threaded_peek(FlexpathQueueNode** queue, thr_mutex_t mutex, thr_condition_t condition) {
    fp_write_log("QUEUE", "peeking at a queue\n");
    thr_mutex_lock(mutex);
    if(*queue==NULL) {
        thr_mutex_unlock(mutex);
        thr_condition_wait(condition, mutex);
        thr_mutex_lock(mutex);
    }
    FlexpathQueueNode* tail;
    tail = *queue;
    while(tail && tail->next) {
        tail=tail->next;
    }
    thr_mutex_unlock(mutex);
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
int get_local_offsets(struct adios_var_struct * list, struct adios_group_struct * g, int** offsets, int** dimensions)
{
    //perr("\t\t\toffsets for var: %s\n", list->name);	    
    struct adios_dimension_struct * dim_list = list->dimensions;	    
    if(dim_list){		
	// if this var has a global dimension, then by default, it has local_offset
	uint16_t gdim_id = dim_list->global_dimension.id;
	uint16_t ldim_id = dim_list->dimension.id;
	if(gdim_id > 0) {	   
	    int num_local_offsets = 0;
	    int * local_offsets = NULL;		
	    int * local_dimensions = NULL;
	    int curr_offset = 0;
	    while(dim_list) {		
		uint16_t offset_id = dim_list->local_offset.id;
		uint16_t ldim_id = dim_list->dimension.id;
		if(offset_id > 0) {							       
		    struct adios_var_struct * tmp_var = adios_find_var_by_id(g->vars, 
									     offset_id);
		    local_offsets = realloc(local_offsets, sizeof(int) * (num_local_offsets+1));
		    memcpy(&local_offsets[curr_offset], tmp_var->data, sizeof(int));
		    // no id, so it must be a literal in the xml doc
		} else {
		    local_offsets = realloc(local_offsets, sizeof(int) * (num_local_offsets+1));
		    local_offsets[curr_offset] = (int)dim_list->local_offset.rank;
		}
		if(ldim_id > 0) {
		    struct adios_var_struct * tmp_var = adios_find_var_by_id(g->vars, ldim_id);
		    local_dimensions = realloc(local_dimensions, sizeof(int) * (num_local_offsets+1));
		    memcpy(&local_dimensions[curr_offset], tmp_var->data, sizeof(int));
		} else {
		    local_dimensions = realloc(local_dimensions, sizeof(int) * (num_local_offsets+1));
		    local_dimensions[curr_offset] = (int)dim_list->dimension.rank;
		}
		dim_list=dim_list->next;
		curr_offset++;
		num_local_offsets++;
	    }
	    offsets[0] = local_offsets;	   
	    dimensions[0] = local_dimensions;
	    return num_local_offsets;
	}	
    }
    offsets = NULL;
    return 0;
}

// creates multiqueue function to handle ctrl messages for given bridge stones 
char* multiqueue_action = "{\n\
    int found = 0;\n\
    int my_rank = -1;\n\
    attr_list mine;\n\
    if(EVcount_varMsg()>0) {\n\
        EVdiscard_and_submit_varMsg(0, 0);\n\
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
    if(EVcount_formatMsg()>0) {\n\
        formatMsg* msg = EVdata_formatMsg(0);\n\
        mine=EVget_attrs_formatMsg(0);\n\
        my_rank= attr_ivalue(mine, \"fp_rank_num\");\n\
    }\n\
    if(EVcount_flush()>0) {\n\
        flush* c = EVdata_flush(0);\n\
        if(c->type == 0) {\n\
            if(EVcount_formatMsg()>0) {\n\
                formatMsg* msg = EVdata_formatMsg(0);\n\
                msg->condition = c->condition;\n\
                printf(\"condition in format msg: \\%d \\n\", msg->condition);\n\
                EVdiscard_flush(0);\n\
                EVsubmit(c->rank+1, msg);\n\
            }\n\
        } else {\n\
            EVdiscard_and_submit_flush(0,0);\n\
        }\n\
    }\n\
    if(EVcount_evgroup()>0){\n\
        evgroup* g = EVdata_evgroup(0);\n\
        mine = EVget_attrs_evgroup(0);\n\
        found = attr_ivalue(mine, \"fp_size\");\n\
        int k;\n\
        for(k=0; k<found; k++){\n\
            EVsubmit(k+1,g);\n\
        }\n\
        EVdiscard_evgroup(0);\n\
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

// construct an fm structure based off the group xml file
FlexpathFMStructure* set_format(struct adios_group_struct* t,struct adios_var_struct* fields, FlexpathWriteFileData* fileData){
    FMStructDescRec *format = (FMStructDescRec*) malloc(sizeof(FMStructDescRec)*2);
    mem_check(format, "format");
    memset(format, 0, sizeof(FMStructDescRec)*2);
    
    FlexpathFMStructure *currentFm = (FlexpathFMStructure *) malloc(sizeof(FlexpathFMStructure));
    mem_check(currentFm, "currentFm");
    memset(currentFm, 0, sizeof(FlexpathFMStructure));

    LIST_INIT(&currentFm->nameList);
    LIST_INIT(&currentFm->dimList);
    currentFm->format = format;
    format->format_name = strdup(t->name);

    if (t->var_count == 0) {
	perr("set_format: No Variables In Group\n");
	return NULL;
    }

    FMFieldList field_list = (FMFieldList) malloc(sizeof(FMField) * (t->var_count + 1));
    if (field_list == NULL) {
	perr("set_format: Field List Memory Allocation Failed");
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
		    uint16_t dim_id = adim->dimension.id;
		    uint16_t gdim_id = adim->global_dimension.id;
		    uint16_t local_id = adim->local_offset.id;
		    if(dim_id > 0) {		    
			struct adios_var_struct *tmp_var = adios_find_var_by_id(t->vars, dim_id);
			char *name = find_fixed_name(currentFm, tmp_var->name);
			char *aname = get_alt_name(tempName,  name);
			dims=add_var(dims, strdup(aname), NULL, 0);
			set_attr_dimensions(tempName, aname, num_dims, fileData->attrs);
		    }
		    if(gdim_id> 0) {
			fileData->globalCount++;
			struct adios_var_struct *tmp_var = adios_find_var_by_id(t->vars, gdim_id);
			char *name = find_fixed_name(currentFm, tmp_var->name);
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
		if (d->dimension.id) {
		    struct adios_var_struct *tmp_var = adios_find_var_by_id(t->vars, d->dimension.id);
		    char *name = find_fixed_name(currentFm, tmp_var->name);
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
		perr( "set_format: Bad Type Error\n");
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
		perr("set_format: Unknown Type Error %d\n", f->type);
		fieldNo--;
		break;
	    }
	}

	fp_write_log("FORMAT","field: %s, %s, %d, %d\n", field_list[fieldNo].field_name, field_list[fieldNo].field_type,field_list[fieldNo].field_size,field_list[fieldNo].field_offset); 
    }

    FlexpathDimNames *d = NULL;
    field_list = (FMFieldList) realloc(field_list, sizeof(FMField) * (altvarcount + t->var_count + 1));

    for (d = currentFm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
	FlexpathAltName *a = NULL;
	for (a = d->altList.lh_first; a != NULL; a = a->entries.le_next) {
	    a->field->field_offset = currentFm->size;
	    currentFm->size += sizeof(int);
	    memcpy(&field_list[fieldNo], a->field, sizeof(FMField));
	    fieldNo++;
	}
    }

    for (; fieldNo < (t->var_count + 1+altvarcount); fieldNo++) {
	field_list[fieldNo].field_type = NULL;
	field_list[fieldNo].field_name = NULL;
	field_list[fieldNo].field_offset = 0;
	field_list[fieldNo].field_size = 0;
    }

    format->field_list = field_list;
    currentFm->format->struct_size = currentFm->size;

    currentFm->buffer = (unsigned char *) malloc(currentFm->size);
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
        fileData->controlMutex, fileData->controlCondition);
    return 0;
}

// terminal action for flush messages: enqueues
static int flush_handler(CManager cm, void* vevent, void* client_data, attr_list attrs) {
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*) client_data;
    Flush_msg* msg = (Flush_msg*) vevent;
    EVtake_event_buffer(cm, msg);
    fp_write_log("MSG", "recieved flush : rank %d type data\n", msg->rank);
    threaded_enqueue(&fileData->controlQueue, msg, DATA_FLUSH, 
        fileData->controlMutex, fileData->controlCondition);
    return 0;
}

// terminal action for op messages: enqueues
static int op_handler(CManager cm, void* vevent, void* client_data, attr_list attrs) {
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*) client_data;
    op_msg* msg = (op_msg*) vevent;
    EVtake_event_buffer(cm, msg);
    fp_write_log("MSG", "recieved op_msg : rank %d type %d: condition: %d\n", 
        msg->process_id, msg->type, msg->condition);
    if(msg->type == 1) {
        threaded_enqueue(&fileData->controlQueue, msg, OPEN, 
            fileData->controlMutex, fileData->controlCondition);
    } else {
        threaded_enqueue(&fileData->controlQueue, msg, CLOSE, 
            fileData->controlMutex, fileData->controlCondition);
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
attr_list set_dst_rank_atom(attr_list attrs, int value) {
    atom_t dst_atom = attr_atom_from_string("fp_dst_rank");
    int dst;
    if(!get_int_attr(attrs, dst_atom, &dst)) {
        add_int_attr(attrs, dst_atom, value);
    }
    set_int_attr(attrs, dst_atom, value);
    return attrs;
}

// sets a dst condition atom
attr_list set_dst_condition_atom(attr_list attrs, int condition){
    atom_t dst_atom = attr_atom_from_string("fp_dst_condition");
    int dst;
    if(!get_int_attr(attrs, dst_atom, &dst)){
	add_int_attr(attrs, dst_atom, condition);
    }
    set_int_attr(attrs, dst_atom, condition);
    return attrs;
}

// processes messages from control queue
int control_thread(void* arg) {
    FlexpathWriteFileData* fileData = (FlexpathWriteFileData*)arg;
    int rank = fileData->rank;
    FlexpathQueueNode* controlMsg;
    FlexpathQueueNode* dataNode;
    while(1) {
        fp_write_log("CONTROL", "control message attempts dequeue\n");
	if((controlMsg = threaded_dequeue(&fileData->controlQueue, 
	    fileData->controlMutex, fileData->controlCondition))) {
            fp_write_log("CONTROL", "control message dequeued\n");
	    if(controlMsg->type==VAR) {
		Var_msg* varMsg = (Var_msg*) controlMsg->data;
		fileData->askedVars = add_var(fileData->askedVars, 
		    strdup(varMsg->var_name), NULL, varMsg->rank);
		EVreturn_event_buffer(flexpathWriteData.cm,controlMsg->data);
	    } else if(controlMsg->type==DATA_FLUSH) { 
		dataNode = threaded_peek(&fileData->dataQueue, 
		    fileData->dataMutex, &fileData->dataCondition);
		Flush_msg* flushMsg = (Flush_msg*) controlMsg->data;
		void* temp = copy_buffer(dataNode->data, flushMsg->rank, fileData);
		fileData->attrs = set_dst_rank_atom(fileData->attrs, flushMsg->rank);
		fileData->attrs = set_dst_condition_atom(fileData->attrs, flushMsg->condition);
		if(!fileData->bridges[flushMsg->rank].opened) {
                  fileData->bridges[flushMsg->rank].opened=1;
                  fileData->openCount++;
                }
		fp_write_log("MSG", " sending data_msg : rank %d step %d\n", 
                    flushMsg->rank, fileData->currentStep);
		EVsubmit_general(fileData->dataSource, temp, data_free, fileData->attrs);
	    } else if(controlMsg->type==OPEN) {
                op_msg* open = (op_msg*) controlMsg->data;
                fileData->bridges[open->process_id].step = open->step;
                fileData->bridges[open->process_id].condition = open->condition;
                if(open->step < fileData->currentStep) {
                    perr("control_thread: Recieved Past Step Open\n");
                } else if (open->step == fileData->currentStep){
                    fp_write_log("STEP", "recieved op with current step\n");
                    thr_mutex_lock(fileData->openMutex);
                    fileData->openCount++;  
                    fileData->bridges[open->process_id].opened = 1;
		    thr_mutex_unlock(fileData->openMutex);
                    op_msg* ack = (op_msg*) malloc(sizeof(op_msg));
                    ack->file_name = strdup(method->group->name);
                    ack->process_id = fileData->rank;
                    ack->step = fileData->currentStep;
                    ack->type = 2;
		    ack->condition = open->condition;
                    fileData->attrs = set_dst_rank_atom(fileData->attrs, open->process_id+1);
		    fp_write_log("MSG", " sending op_msg : dst %d step %d type ack\n",
                        open->process_id, fileData->currentStep); 
                    EVsubmit_general(fileData->opSource, ack, op_free, fileData->attrs);
                } else {
                    fp_write_log("STEP", "recieved op with future step\n");
                }
            } else if(controlMsg->type==CLOSE) {
                op_msg* close = (op_msg*) controlMsg->data;
		thr_mutex_lock(fileData->openMutex);
		fileData->openCount--;
                fileData->bridges[close->process_id].opened=0;
		thr_mutex_unlock(fileData->openMutex);
                 if(fileData->openCount==0) {
                    fp_write_log("STEP", "advancing\n");
		    FlexpathQueueNode* node = threaded_dequeue(&fileData->dataQueue, 
		        fileData->dataMutex, fileData->dataCondition);
                    FMfree_var_rec_elements(fileData->fm->ioFormat, node->data);
                    fileData->currentStep++;
                    
                    int i;
                    //for all bridges if step == currentstep send ack
                    for(i=0; i<fileData->numBridges; i++) {
                      if(fileData->bridges[i].step==fileData->currentStep) {
                        fileData->openCount++;
                        fileData->bridges[i].opened = 1;
                        op_msg* ack = (op_msg*) malloc(sizeof(op_msg));
                        ack->file_name = strdup(method->group->name);
                        ack->process_id = fileData->rank;
                        ack->step = fileData->currentStep;
                        ack->type = 2;
			ack->condition = fileData->bridges[i].condition;
                        fileData->attrs = set_dst_rank_atom(fileData->attrs, i+1);
		        fp_write_log("MSG", " sending op_msg : dst %d step %d type ack\n",
                            i, fileData->currentStep);
                        EVsubmit_general(fileData->opSource, ack, op_free, fileData->attrs);
                      }
                    }
		}
	    } else {
		perr("control_thread: Unrecognized Control Message\n");
	    }
	}
    }
    return 0;
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
    gen_pthread_init();
    setenv("CMSelfFormats", "1", 1);
    
    // fork communications thread
    int forked = CMfork_comm_thread(flexpathWriteData.cm);   
    if(!forked) {
         perr( "error forking comm thread\n");
    }
}

// opens a new adios file for writes
extern int adios_flexpath_open(struct adios_file_struct *fd, struct adios_method_struct *method, void*comm) 
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

    FlexpathWriteFileData* fileData = (FlexpathWriteFileData *) malloc(sizeof(FlexpathWriteFileData));
    mem_check(fileData, "fileData");
    memset(fileData, 0, sizeof(FlexpathWriteFileData));
    fileData->maxQueueSize=0;
    if(method->parameters) {
        sscanf(method->parameters,"QUEUE_SIZE=%d;",&fileData->maxQueueSize);
        fp_write_log("SETUP", "setting max queue size to %d\n", fileData->maxQueueSize);
    }
    
    // setup step state
    fileData->attrs = create_attr_list();
    fileData->openCount = 0;
    fileData->currentStep = 0;

    // setup mutexs
    fileData->controlMutex = thr_mutex_alloc();
    fileData->dataMutex = thr_mutex_alloc();
    fileData->openMutex = thr_mutex_alloc();
    
    // setup conditions
    fileData->controlCondition = thr_condition_alloc();
    fileData->dataCondition = thr_condition_alloc();

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
    MPI_Comm * group_comm = (MPI_Comm*)malloc(sizeof(MPI_Comm));
    adios_var_to_comm(fd->group->group_comm,
	fd->group->adios_host_language_fortran, comm, group_comm);
    fileData->mpiComm = (MPI_Comm*)group_comm;
    MPI_Comm_rank(*(fileData->mpiComm), &fileData->rank);
    MPI_Comm_size(*(fileData->mpiComm), &fileData->size);
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
        CONTACT_STR_LEN, MPI_CHAR, 0, *(fileData->mpiComm));

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
    char in_contact[CONTACT_STR_LEN] = "";
    int numBridges = 0;
    int stone_num;
    // build a bridge per line
    while(fscanf(reader_info, "%d:%s",&stone_num, in_contact)!=EOF){
        fileData->bridges = realloc(fileData->bridges, sizeof(FlexpathStone) * (numBridges + 1));
        attr_list contact_list = attr_list_from_string(in_contact);
        fileData->bridges[numBridges].myNum = EVcreate_bridge_action(flexpathWriteData.cm, contact_list, stone_num);
        fileData->bridges[numBridges].opened = 0;
        fileData->bridges[numBridges].step = 0;
        fileData->bridges[numBridges].theirNum = stone_num;
        fileData->bridges[numBridges].contact = strdup(in_contact);
        numBridges += 1;
    }
    fileData->numBridges = numBridges;
    fclose(reader_info);

    MPI_Barrier(*(fileData->mpiComm));
    
    // cleanup of reader files (writer is done with it).
    if(fileData->rank == 0){
	unlink(reader_info_filename);
	unlink(reader_ready_filename);
    }
	
    //process group format
    struct adios_group_struct *t = method->group;
    struct adios_var_struct *fields = t->vars;
    if(t == NULL)
	perr("t is null\n");
    if(fields == NULL)
	perr("t is null\n");

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
    FMStructDescList queue_list[] = {flush_format_list, format_format_list, 
	var_format_list, op_format_list, evgroup_format_list, 
        data_format_list, NULL};
    char* q_action_spec = create_multityped_action_spec(queue_list, 
        multiqueue_action); 
    EVaction multi_action = EVassoc_multi_action(flexpathWriteData.cm, 
	fileData->multiStone, q_action_spec, NULL);
    fileData->formatSource = EVcreate_submit_handle(flexpathWriteData.cm, 
        fileData->multiStone, format_format_list);
    fileData->dataSource = EVcreate_submit_handle_free(flexpathWriteData.cm, 
        fileData->multiStone, fileData->fm->format, data_free,  NULL); 
    fileData->opSource = EVcreate_submit_handle_free(flexpathWriteData.cm, 
        fileData->multiStone, op_format_list, op_free,  NULL); 
    fileData->offsetSource = EVcreate_submit_handle(flexpathWriteData.cm, 
	fileData->multiStone, evgroup_format_list);
    
    fp_write_log("SETUP", "setup terminal actions\n");
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
	var_format_list, var_handler, fileData);
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
	op_format_list, op_handler, fileData);
    EVassoc_terminal_action(flexpathWriteData.cm, fileData->sinkStone, 
	flush_format_list, flush_handler, fileData);

    //link multiqueue to sink
    fp_write_log("SETUP", "linking stones\n");
    EVaction_set_output(flexpathWriteData.cm, fileData->multiStone, 
        multi_action, 0, fileData->sinkStone);

    //link up multiqueue ports to bridge stones
    for(i=0; i<numBridges; i++) {
        EVaction_set_output(flexpathWriteData.cm, 
            fileData->multiStone, multi_action, i+1, fileData->bridges[i].myNum);
    }
    
    fp_write_log("SETUP", "arranged evpath graph\n");
	
    //store format id in multiqueue
    Format_msg *initial_format_msg = malloc(sizeof(Format_msg));
    FMContext my_context = create_local_FMcontext();	
    fileData->fm->ioFormat = register_data_format(my_context, fileData->fm->format);
    int id_len;
    char* temp = get_server_ID_FMformat(fileData->fm->ioFormat, &id_len);
    initial_format_msg->format_id = temp;
    initial_format_msg->id_len = id_len;
    int rep_len;
    char *temp2 = get_server_rep_FMformat(fileData->fm->ioFormat, &rep_len);
    initial_format_msg->rep_id = temp2;
    initial_format_msg->rep_id_len = rep_len;
    
    fp_write_log("SETUP", "submitting format stuff\n");
    EVsubmit_general(fileData->formatSource, initial_format_msg, format_free, fileData->attrs);
    
    fp_write_log("SETUP", "indicating to reader that ready\n");
    sprintf(writer_ready_filename, "%s_%s", fd->name, "writer_ready.txt");
    if(fileData->rank == 0) {
        FILE* writer_info = fopen(writer_ready_filename, "w");
        fprintf(writer_info, "ready");
        fclose(writer_info);
    }
        
    fp_write_log("SETUP", "fork control thread\n");
    thr_thread_t forked_thread = thr_fork(control_thread, fileData);
    if(!forked_thread) {
        perr("on open ERROR forking control thread");
    }
   
    return 0;	
}




//  writes data to multiqueue
extern void adios_flexpath_write(struct adios_file_struct *fd, struct adios_var_struct *f, void *data, struct adios_method_struct *method) {
    fp_write_log("FILE", "entering flexpath file %s write\n", method->group->name);
    FlexpathWriteFileData* fileData = find_open_file(method->group->name);
    FlexpathFMStructure* fm = fileData->fm;

    if (fm == NULL)
    {
	return;

    }
    
    FMFieldList flist = fm->format->field_list;
    FMField *field = NULL;
    char *fixedname = find_fixed_name(fm, f->name);
    field = internal_find_field(fixedname, flist);
    //perr( "found field %s\n", field->field_name);
    if (field != NULL) {
	if (!f->dimensions) {
	    //scalar quantity
            //perr( "copying scalar value\n");
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
		//perr( "no data for  scalar %s\n", f->name);
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
        
        //perr( "field name: %s\n", fields->name);
        if(fields->dimensions) {
            //perr( "field is an array\n");
            struct adios_dimension_struct* dims = fields->dimensions;
            //perr( "field dims: %p\n", dims);
    
            int total_size = 1;
            //for each dimension
            while(dims) {    
                struct adios_var_struct* temp = adios_find_var_by_id(g2->vars, dims->dimension.id);            
                int size = *(int*)temp->data;
                //perr( "dim %s size %d\n", temp->name, *(int*)temp->data);
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
            //perr( "field %s field size %d\n", fields->name, total_size);
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

    threaded_enqueue(&fileData->dataQueue, buffer, 
        DATA_BUFFER, fileData->dataMutex, fileData->dataCondition);
    
    int c = 0;
 
    // now gather offsets and send them via MPI to root
    struct adios_group_struct * g = fd->group;
    struct adios_var_struct * list = g->vars;

    if(fileData->globalCount > 0 && !fileData->sentGlobalOffsets){	
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
		if(myrank == 0){
		    int buf_size = num_local_offsets * commsize * sizeof(int);		    
		    all_offsets = (int*)malloc(buf_size);		
		    all_local_dims = (int*)malloc(buf_size);
		}

		MPI_Gather(local_offsets, num_local_offsets, MPI_INT, 
			   all_offsets, num_local_offsets, MPI_INT,
			   0, *fileData->mpiComm);

		MPI_Gather(local_dimensions, num_local_offsets, MPI_INT, 
			   all_local_dims, num_local_offsets, MPI_INT,
			   0, *fileData->mpiComm);

		if(myrank == 0){
		    num_gbl_vars++;
		    offset_struct * ostruct = (offset_struct*)malloc(sizeof(offset_struct));
		    ostruct->offsets_per_rank = num_local_offsets;
		    ostruct->total_offsets = num_local_offsets * commsize;
		    ostruct->local_offsets = all_offsets;
		    ostruct->local_dimensions = all_local_dims;
		    gbl_vars = realloc(gbl_vars, sizeof(global_var) * num_gbl_vars);
		    gbl_vars[num_gbl_vars - 1].name = strdup(list->name);
		    gbl_vars[num_gbl_vars - 1].noffset_structs = 1;
		    perr("\n\n\n\t\tnoffset_structs: %d\n", 
		     gbl_vars[num_gbl_vars - 1].noffset_structs);
		    gbl_vars[num_gbl_vars - 1].offsets = ostruct;
		    int i;			   
		    i = 0;
		    perr("\t\t\tall offsets for var: %s\n", list->name);
		    while(i<commsize * num_local_offsets){
			int j;
			
			for(j=0; j<num_local_offsets;j++){
			    perr("\t\t\t%d ", all_offsets[i]);  
			    i++;
			}			
			perr("\n");				
			
		    }
		    perr("\n");
		    perr("\t\t\tall local_dims for var: %s\n", list->name);
		    i = 0;
		    
		    while(i<commsize * num_local_offsets){
			int j;				
			for(j=0; j<num_local_offsets;j++){
			    perr("\t\t\t%d ", all_local_dims[i]);  
			    i++;
			}			
			perr("\n");
			
		    }
		    perr("\n");
		    
		}		
	    }
	    list=list->next;
	}
	if(myrank == 0){
	    int i;
	    
	    for(i=0; i<num_gbl_vars; i++){
		perr("global_var: %s has local offsets\n", gbl_vars[i].name);
	    }
	    
	    evgroup * gp = (evgroup*)malloc(sizeof(evgroup));
	    gp->num_vars = num_gbl_vars;
	    perr("num global vars %d\n", num_gbl_vars);
            gp->vars = gbl_vars;
	    fileData->gp = gp;
	    fileData->attrs = set_size_atom(fileData->attrs, fileData->size);
	    perr("size:%d\n\n\n", fileData->size);
            EVsubmit(fileData->offsetSource, gp, fileData->attrs);
	}
	fileData->sentGlobalOffsets = 1;
    }

    while((c=queue_count(&fileData->dataQueue, fileData->dataMutex))>fileData->maxQueueSize) {
        fp_write_log("QUEUE", "waiting for queue to be below max size\n");
        thr_condition_wait(fileData->dataCondition, fileData->dataMutex);
    }
    fp_write_log("FILE", "file close %s exiting\n", method->group->name);
}

// wait until all open files have finished sending data to shutdown
extern void adios_flexpath_finalize(int mype, struct adios_method_struct *method) {
    FlexpathWriteFileData* fileData = flexpathWriteData.openFiles;
    while(fileData) {
        thr_mutex_lock(fileData->dataMutex);
        while(fileData->dataQueue!=NULL) {
            thr_condition_wait(fileData->dataCondition, fileData->dataMutex);
        }
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

#else // print empty version of all functions (if NO_FLEXPATH == 1)

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
