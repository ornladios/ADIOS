/*
    adios_flexpath.c
    Originally copied from adios_datatap.c
    Goal: to use evpath for io in conjunction with read/read_flexpath.c
*/


#if NO_FLEXPATH == 0

// system libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ffs.h>
#include <atl.h>
#include <evpath.h>
#include <fm.h>
#include <gen_thread.h>
#include <sys/queue.h>
#include <pthread.h>

// local libraries
#include "config.h"
#include "public/adios.h"
#include "core/adios_internals.h"
#include "core/adios_transport_hooks.h"
#include "core/util.h"
#include "public/flexpath.h"

//static LIST_HEAD(listhead, _fm_structure) globallist;


/************************* Global Variable Declarations *************************/
#define STARTINGSIZE 16
#define OPLEN 7
static char opList[OPLEN] = { '+', '-', '*', '/', '.', '>', '<' };
static char *opRepList[OPLEN] = { "_plus_", "_minus_", "_mult_", "_div_", "_dot_", "_greater_", "_less_" };
typedef enum {VAR=0, DATA_FLUSH, OPEN, CLOSE, DATA_BUFFER, OFFSET_MSG} MsgType;

/************************* EVPath/FFS Definitions *******************************/
/* NOTE: may need to move to a shared reader/writer header at some point for versioning control */

typedef struct _stone {
    int myNum;
    int theirNum;
    int step;
    char* contact;
} Stone;

typedef struct _var_node {
    char* varName;
    struct _var_node* dims;
    struct _var_node* next;
    int rank;
} VarNode;

typedef struct _QueueNode {
    void* data;
    MsgType type;
    struct _QueueNode* next;
} QueueNode;

/************************* Structure and Type Definitions ***********************/

typedef struct _name_table {
    char *originalName;
    char *mangledName;
    LIST_ENTRY(_name_table) entries;
} NameTable;

typedef struct _alt_name {
    char *name;
    FMField *field;                 //contains the field list
    LIST_ENTRY(_alt_name) entries;
} AltName;

typedef struct _dim_names {
    char *name;
    LIST_HEAD(alts, _alt_name) altList;
    LIST_ENTRY(_dim_names) entries;
} DimNames;

struct _fm_structure {
    FMStructDescRec *format;
    int size;	            //in bytes - no padding
    unsigned char *buffer;    //big enough to hold top-level data.
    int sndCount;
    FMFormat ioFormat;
    attr_list attrList;	
    LIST_HEAD(tableHead, _name_table) nameList;
    LIST_HEAD(dims, _dim_names) dimList;
};

typedef struct _local_write_data {
    // MPI stuff
    MPI_Comm * mpiComm;
    int rank;
    int size;

    // EVPath stuff
    CManager cm;
    EVstone multiStone;
    EVstone sinkStone;
    EVsource formatSource;
    EVsource dataSource;
    EVsource offsetSource;
    EVsource opSource;
    Stone* bridges;
    int numBridges;
    attr_list attrs;
    atom_t CM_TRANSPORT;

    // server state
    int openCount;
    thr_mutex_t openMutex;
    int setupCorrect;
    int cycleId;
    char *pFile;
    struct _fm_structure *fm;
    VarNode* askedVars;
    VarNode* writtenVars;
    VarNode* formatVars;
    QueueNode* controlQueue;
    thr_mutex_t controlMutex;
    int controlCondition;
    QueueNode* dataQueue;    
    thr_mutex_t dataMutex;
    int dataCondition;

    //global array distribution data;
    int global_count; //field to keep track if this file is handling global arrays.
    int sent_global_offsets;
    evgroup *gp;
} FlexpathWriteData;

// local storage pointer
FlexpathWriteData* localWriteData = NULL;
int currentStep = 0;

/**************************** Function Definitions *********************************/

static void adios_var_to_comm (
    const char * comm_name
    ,enum ADIOS_FLAG host_language_fortran
    ,void * data
    ,MPI_Comm * comm
    )
{
    if (data){    
        int t = *(int *) data;
        if (!comm_name){        
            if (!t){            
                fprintf (stderr, "communicator not provided and none "
			 "listed in XML.  Defaulting to "
			 "MPI_COMM_SELF\n"
		    );
                *comm = MPI_COMM_SELF;
            }
            else{            
                if (host_language_fortran == adios_flag_yes){                
                    *comm = MPI_Comm_f2c (t);
                }
                else{                
                    *comm = *(MPI_Comm *) data;
                }
            }
        }
	else{        
            if (!strcmp (comm_name, "")){            
                if (!t){                
                    fprintf (stderr, "communicator not provided and none "
			     "listed in XML.  Defaulting to "
			     "MPI_COMM_SELF\n"
			);
                    *comm = MPI_COMM_SELF;
                }
                else{                
                    if (host_language_fortran == adios_flag_yes){                    
                        *comm = MPI_Comm_f2c (t);
                    }
                    else{                    
                        *comm = *(MPI_Comm *) data;
                    }
                }
            }
            else{            
                if (!t){                
                    fprintf (stderr, "communicator not provided but one "
			     "listed in XML.  Defaulting to "
			     "MPI_COMM_WORLD\n"
			);
                    *comm = MPI_COMM_WORLD;
                }
                else{                
		    if(host_language_fortran == adios_flag_yes){                    
                        *comm = MPI_Comm_f2c (t);
                    }
                    else{                   
                        *comm = *(MPI_Comm *) data;
                    }
                }
            }
        }
    }
    else
    {
        fprintf (stderr, "coordination-communication not provided. "
		 "Using MPI_COMM_WORLD instead\n"
	    );
        *comm = MPI_COMM_WORLD;
    }
}

void
set_attr_dimensions(char* tempName, char * aname, int num_dims)
{
    char atom_name[200];
    atom_name[0] = '\0';
    strcat(atom_name, tempName);
    strcat(atom_name, "_");
    strcat(atom_name, FP_DIM_ATTR_NAME);
    strcat(atom_name, "_");
    char dim_num[10] = "";
    sprintf(dim_num, "%d", num_dims);
    strcat(atom_name, dim_num);
    atom_t dim_atom = attr_atom_from_string(atom_name);
    add_string_attr(localWriteData->attrs, dim_atom, aname);
    //fprintf(stderr, "added attr %s with value %s\n\n\n", atom_name, aname);
    // attach 0 ndims for that dimension alt name
    atom_name[0] = '\0';
    strcat(atom_name, aname);
    strcat(atom_name, "_");
    strcat(atom_name, FP_NDIMS_ATTR_NAME);
    atom_t ndims_atom = attr_atom_from_string(atom_name);
    add_int_attr(localWriteData->attrs, ndims_atom, 0);
}

void data_free(void* temp, void* a) {
    //fprintf(stderr, "\n\n\n\n\nfree\n\n\n\n\n\n");
    //free(temp);
}

void op_free(void* temp, void* a) {
    op_msg* op = (op_msg*) temp;
    //fprintf(stderr, "\n\nop free\n\n\n");
    free(op);
}


//add to head
void threaded_enqueue(QueueNode** queue, void* item, MsgType type, thr_mutex_t mutex, int* cond) {
    //if(localWriteData->rank>0) fprintf(stderr, "rank %d enter enqueue\n", localWriteData->rank);
    //fprintf(stderr, "attempting enqueue\n");
    thr_mutex_lock(mutex);
    //fprintf(stderr, "got lock\n");
    QueueNode* newNode = (QueueNode*) malloc(sizeof(QueueNode));
    newNode->data = item;
    newNode->type = type;
    newNode->next = *queue;
    *queue = newNode;
    //fprintf(stderr, "made node with data %p\n", newNode->data);
    //if(newNode->next) {
    //    fprintf(stderr, "next node has data %p\n", newNode->next->data);
    //}
    //fprintf(stderr, "released lock\n");
    //if(localWriteData->rank>0) fprintf(stderr, "rank %d signalling %d\n", localWriteData->rank, *cond);
    //fprintf(stderr, "start signal1\n");
    CMCondition_signal(localWriteData->cm, *cond);
    *cond = CMCondition_get(localWriteData->cm, NULL);
    //fprintf(stderr, "done signal1\n");

    //if(localWriteData->rank>0) fprintf(stderr, "rank %d exit enqueue\n", localWriteData->rank);
    //reset condition
    thr_mutex_unlock(mutex);
     //fprintf(stderr, "signaled exiting\n");
}

//remove from tail
QueueNode* threaded_dequeue(QueueNode** queue, thr_mutex_t mutex, int* cond) {
    //if(localWriteData->rank>0) fprintf(stderr, "rank %d enter dequeue\n", localWriteData->rank);
    thr_mutex_lock(mutex);
    if(*queue==NULL) {
        //empty queue, wait for item
        thr_mutex_unlock(mutex);
        //if(localWriteData->rank>0) fprintf(stderr, "rank %d dequeue waiting on %d\n", localWriteData->rank, *cond);
        //fprintf(stderr, "start wait1\n");
        int res = CMCondition_wait(localWriteData->cm, *cond);
        //fprintf(stderr, "done wait1 %d\n", res);
        thr_mutex_lock(mutex);
        //if(localWriteData->rank>0) fprintf(stderr, "rank %d dequeue new cond %d\n", localWriteData->rank, *cond);
    }
    QueueNode* tail;
    QueueNode* prev = NULL;
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
    thr_mutex_unlock(mutex);
    //if(localWriteData->rank>0) fprintf(stderr, "rank %d exit dequeue\n", localWriteData->rank);
    return tail;  // match holds either the matching item or NULL
}

//peek at tail
QueueNode* threaded_peek(QueueNode** queue, thr_mutex_t mutex, int* cond) {
    //fprintf(stderr, "rank %d enter peek\n", localWriteData->rank);
    thr_mutex_lock(mutex);
    if(*queue==NULL) {
        //if(localWriteData->rank>0) fprintf(stderr, "rank %d peek waiting\n", localWriteData->rank);
        thr_mutex_unlock(mutex);
        //fprintf(stderr, "start wait2\n");
        int res = CMCondition_wait(localWriteData->cm, *cond);
        //fprintf(stderr, "done wait2 %d\n", res);
        thr_mutex_lock(mutex);
    }
    QueueNode* tail;
    tail = *queue;
    while(tail && tail->next) {
        tail=tail->next;
    }
    thr_mutex_unlock(mutex);
    //fprintf(stderr, "rank %d exit peek\n", localWriteData->rank);
    return tail;
}


VarNode* add_var(VarNode* queue, char* varName, VarNode* dims, int rank){
    if(queue) {
        queue->next=add_var(queue->next, varName, dims, rank);
        return queue;
    } else {
        queue = (VarNode*) malloc(sizeof(VarNode));
        queue->varName = strdup(varName);
        queue->dims = dims;
        queue->next = NULL;
        queue->rank = rank;
        return queue;
    }
}

void free_vars(VarNode* queue){
    if(queue) {
        free_vars(queue->next);
        free(queue->varName);
        free(queue);
    }
}

VarNode* queue_contains(VarNode* queue, const char* name, int rank) {
    int compare_rank = 0;
    if(rank >= 0 ) {
        compare_rank = 1;
    }
    VarNode* tmp = queue;
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


void print_format(FMStructDescRec *format) {
    FMField *f =  format->field_list;
    int x = 0;
    for (x=0; 1; x++) {
      f = &format->field_list[x];
      if (f == NULL || f->field_name == NULL || f->field_size == 0)
	break;
    }
}

// return name with operators removed
static char *get_fixed_name(char *name) {
    char *oldName = strdup(name);
    char *newName = (char *) malloc(sizeof(char) * 255);
    int i;
    for (i=0; i< OPLEN; i++){
        char op[]={opList[i],'\0'};
        char *opRep=opRepList[i];
        char *token = strtok(oldName, op);
        char *lastTok=NULL;
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

// return name with operators removed from lookup list
static char *find_fixed_name(struct _fm_structure *fm, char *name) {
    NameTable *node;
    for (node = fm->nameList.lh_first; node != NULL; node = node->entries.le_next) {
        if (!strcmp(node->originalName, name)) {
	    return node->mangledName;
        }
    }
    return name;
}

// return a list of all the names associated with the variable
static char *get_alt_name(char *name, char *dimName) {
    //fprintf(stderr, "debug: get_alt_name\n");
    int len = strlen(name) + strlen(dimName) + 2;
    char *newName = (char *) malloc(sizeof(char) * len);
    strcpy(newName, dimName);
    strcat(newName, "_");
    strcat(newName, name);
    return newName;
}

// returns a name with the appropriate dimension appended
static AltName *find_alt_name(struct _fm_structure *currentFm, char *dimName, char *varName) {
    //fprintf(stderr, "debug: findAltName\n");
    char *altName = get_alt_name(varName, dimName);
    DimNames *d = NULL;

    // move to dim name in fm dim name list
    for (d = currentFm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
        if (!strcmp(d->name, dimName)) {
	    break;
	}
    }

    // if reached end of list - create list with current dim name at head
    if (d == NULL) {
        d = (DimNames *) malloc(sizeof(DimNames));
        d->name = dimName;
        LIST_INIT(&d->altList);
        LIST_INSERT_HEAD(&currentFm->dimList, d, entries);
    }

    // d now points to an entry in fm dimList with name == dimname
    
    // create AltName structure and field with alternative name in it 
    AltName *a = (AltName *) malloc(sizeof(AltName));
    a->name = altName;
    FMField *field = (FMField *) malloc(sizeof(FMField));
    a->field = field;
    field->field_name = strdup(altName);
    // TO FIX: Should really check datatype (another paramater?)
    field->field_type = strdup("integer");
    field->field_size = sizeof(int);
    field->field_offset = -1;
    
    // insert AltName structure into dimname list
    LIST_INSERT_HEAD(&d->altList, a, entries);
    return a;
}

int
get_local_offsets(struct adios_var_struct * list, 
		  struct adios_group_struct * g, 
		  int** offsets,
		  int** dimensions)
{
    //perr("\t\t\toffsets for var: %s\n", list->name);	    
    struct adios_dimension_struct * dim_list = list->dimensions;	    
    if(dim_list){		
	// if this var has a global dimension, then by default, it has local_offset
	uint16_t gdim_id = dim_list->global_dimension.id;
	uint16_t ldim_id = dim_list->dimension.id;
	if(gdim_id > 0){	   
	    int num_local_offsets = 0;
	    int * local_offsets = NULL;		
	    int * local_dimensions = NULL;
	    int curr_offset = 0;
	    while(dim_list){		
		uint16_t offset_id = dim_list->local_offset.id;
		uint16_t ldim_id = dim_list->dimension.id;
		if(offset_id > 0){							       
		    struct adios_var_struct * tmp_var = adios_find_var_by_id(g->vars, 
									     offset_id);
		    local_offsets = realloc(local_offsets, sizeof(int) * (num_local_offsets+1));
		    memcpy(&local_offsets[curr_offset], tmp_var->data, sizeof(int));
		    // no id, so it must be a literal in the xml doc
		}else{
		    local_offsets = realloc(local_offsets, sizeof(int) * (num_local_offsets+1));
		    local_offsets[curr_offset] = (int)dim_list->local_offset.rank;
		}
		if(ldim_id > 0){
		    struct adios_var_struct * tmp_var = adios_find_var_by_id(g->vars, ldim_id);
		    local_dimensions = realloc(local_dimensions, sizeof(int) * (num_local_offsets+1));
		    memcpy(&local_dimensions[curr_offset], tmp_var->data, sizeof(int));
		}else{
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

void set_field(int type, FMFieldList* field_list_ptr, int fieldNo, int* size){
    //fprintf(stderr, "debug: set_field\n");
    FMFieldList field_list = *field_list_ptr;
    switch (type) {
	case adios_unknown:
	  //fprintf(stderr, "bad type error\n");
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
	  //fprintf(stderr, "unknown type error\n");
	  break;
    }
    *field_list_ptr = field_list;
}


static FMField *internal_find_field(char *name, FMFieldList flist) {
    //fprintf(stderr, "debug: internal_find_field\n");
    FMField *f = flist;
    while (f->field_name != NULL && strcmp(f->field_name, name))
    {
	f++;
    }

    return f;
}

struct _fm_structure* setFormat(struct adios_group_struct* t,struct adios_var_struct* fields){
 //iterate through all the types
    //create a format rec -- size 2 OK because only single-level structs allowed through ADIOS.
    FMStructDescRec *format = (FMStructDescRec *) malloc(sizeof(FMStructDescRec) * 2);
    //fprintf(stderr, "entering setFormat\n");

    // attach rank attr
    atom_t rank_atom = attr_atom_from_string(FP_RANK_ATTR_NAME);
    add_int_attr(localWriteData->attrs, rank_atom, localWriteData->rank);
    //fprintf(stderr, "added attr %s with value %d\n", FP_RANK_ATTR_NAME, localWriteData->rank);

    if (format == NULL) {
	perror("memory allocation failed");
	return NULL;
    }
    memset(format, 0, sizeof(FMStructDescRec) * 2);
    struct _fm_structure *currentFm = (struct _fm_structure *) malloc(sizeof(struct _fm_structure));
    if (currentFm == NULL) {
        perror("memory allocation failed");
	return NULL;
    }

    memset(currentFm, 0, sizeof(struct _fm_structure));

    //store writer's communicator size in attribute list.
    currentFm->attrList = create_attr_list();
    set_int_attr(currentFm->attrList, attr_atom_from_string("mpisize"), localWriteData->size);
	  
    LIST_INIT(&currentFm->nameList);
    LIST_INIT(&currentFm->dimList);
  
    //associate the FMStructDescRec with the _fm_structure
    currentFm->format = format;
    /*TO FIX: Should be sure that the PG name t->name doesn't need name mangling to be a valid format name string.*/
    format->format_name = strdup(t->name);

    //allocate field list
    if (t->var_count == 0) {
	//fprintf(stderr, "no variables in this group - possibly an error\n");
	return NULL;
    }

    int altvarcount = 0;

    FMFieldList field_list = (FMFieldList) malloc(sizeof(FMField) * (t->var_count + 1));
    if (field_list == NULL) {
	perror("memory allocation failed");
	return NULL;
    }

    //keep count of the total number of fields -- will count up to t->var_count
    int fieldNo = 0;

    //for each type look through all the fields
    struct adios_var_struct *f;
    for (f = t->vars; f != NULL; f = f->next, fieldNo++) {
	//make the field list
	//check name for + - * / (operators) and replace them
	char *tempName = get_fixed_name(f->name);
	if (strcmp(tempName, f->name)) {
	    //strings don't match
	    //add to name list
	    NameTable *nameNode = (NameTable *) malloc(sizeof(NameTable));
	    nameNode->originalName = strdup(f->name);
	    nameNode->mangledName = strdup(tempName);
	    LIST_INSERT_HEAD(&currentFm->nameList, nameNode, entries);
	}

	//use the mangled name for the field.
	field_list[fieldNo].field_name = tempName;
        if(tempName!=NULL) {
            //fprintf(stderr, "looking at var %s\n", tempName);
            int num_dims = 0;
            char atom_name[200];
            VarNode* dims=NULL;
            if(f->dimensions) {
                struct adios_dimension_struct* adim = f->dimensions;  
		
                for(; adim != NULL; adim = adim->next){
                    num_dims++;		    
		    // this part is for local arrays only.
		    uint16_t dim_id = adim->dimension.id;
		    uint16_t gdim_id = adim->global_dimension.id;
		    uint16_t local_id = adim->local_offset.id;
		    if(dim_id > 0){		    
			struct adios_var_struct *tmp_var = adios_find_var_by_id(t->vars, dim_id);
			char *name = find_fixed_name(currentFm, tmp_var->name);
			char *aname = get_alt_name(tempName,  name);
			//perr("\t\t\ttmp_var->name: %s name: %s, aname: %s\n", tmp_var->name, name, aname);
			dims=add_var(dims, strdup(aname), NULL, 0);
			// attach a dimension attr
			set_attr_dimensions(tempName, aname, num_dims);
		    }
		    if(gdim_id> 0){
			localWriteData->global_count++;
			struct adios_var_struct *tmp_var = adios_find_var_by_id(t->vars, gdim_id);
			char *name = find_fixed_name(currentFm, tmp_var->name);
			char *aname = get_alt_name(tempName, name);
			//perr("\t\t\tgbl tmp_var->name: %s name: %s, aname: %s\n", tmp_var->name, name, aname);
			dims=add_var(dims, strdup(aname), NULL, 0);
			set_attr_dimensions(tempName, aname, num_dims);			
		    }
                }
            }
            // attach ndims attr
            atom_name[0] = '\0';
            strcat(atom_name, tempName);
            strcat(atom_name, "_");
            strcat(atom_name, FP_NDIMS_ATTR_NAME);
            atom_t ndims_atom = attr_atom_from_string(strdup(atom_name));
            add_int_attr(localWriteData->attrs, ndims_atom, num_dims);
            //fprintf(stderr, "added attr %s with value %d\n\n\n", atom_name, num_dims);
            //fprintf(stderr, "adding %s %d to %p\n", tempName, num_dims, localWriteData->formatVars);
            localWriteData->formatVars = add_var(localWriteData->formatVars, tempName, dims, 0);
        }
	// if its a single field
	if (!f->dimensions) {
	    // set the field type size and offset approrpriately
	    set_field(f->type, &field_list, fieldNo, &currentFm->size);
	}
	else
	{
	    //it's a vector!
	    //find out the dimensions by walking the dimension list

#define DIMSIZE 10240
#define ELSIZE 256
	    struct adios_dimension_struct *d = f->dimensions;
	    char dims[DIMSIZE] = { 0 };
	    char el[ELSIZE] = { 0 };
	    int v_offset=-1;
		  
		  
	    //create the textual representation of the dimensions
	    for (; d != NULL; d = d->next)
	    {
		//for each dimension just take the upper_bound
		if (d->dimension.id)
		{
		    //find_fixed_name returns the mangled name from the original name
		    struct adios_var_struct *tmp_var = adios_find_var_by_id(t->vars, d->dimension.id);
		    char *name =
			find_fixed_name(currentFm, 
				      tmp_var->name);
		    //create the alternate name for this variable and the array its defining
		    AltName *a = find_alt_name(currentFm, name,
					     (char*)field_list[fieldNo].field_name);
		    //AltName is a new variable that we need to add to the field list
		    altvarcount++;
			  
		    snprintf(el, ELSIZE, "[%s]", a->name);
		    //If any of the dimentions of the array are variable, we'll have a pointer in the data struct.
		    v_offset = 0;
			  
		}
		else			//it's a number
		{
		    snprintf(el, ELSIZE, "[%llu]", d->dimension.rank);
		    //if it's a number the offset will be the size of the variable*rank
		    //if it's multidimensional with any variable, v_offset=0.  If they're all integer, v_offset = -1*allocation
		    v_offset *= d->dimension.rank;
		}
		strncat(dims, el, DIMSIZE);
	    }
	    v_offset *= -1; // make it positive.  Actual offset = v_offset*sizeof(var).
		  
	    while(currentFm->size % 8 != 0)
	    {
		currentFm->size ++;					
	    }
		  
	    switch (f->type)
	    {
	    case adios_unknown:
		//fprintf(stderr, "bad type error\n");
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
		//To Fix:  Do we have to worry about char[30]?  I think so... did the v_offset logic work in that case?
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
		//fprintf(stderr, "unknown type error %d\n", f->type);
		fieldNo--;
		break;
	    } //end switch

	} //end vector handling case

	//      fprintf(formatfile, "%s, %s, %d, %d\n", field_list[fieldNo].field_name, field_list[fieldNo].field_type,field_list[fieldNo].field_size,field_list[fieldNo].field_offset); 


    } //end handling of basic variable list

    DimNames *d = NULL;
    field_list = (FMFieldList) realloc(field_list, sizeof(FMField) * (altvarcount + t->var_count + 1));

    for (d = currentFm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
	AltName *a = NULL;
	for (a = d->altList.lh_first; a != NULL; a = a->entries.le_next)
	{
	    a->field->field_offset = currentFm->size;
	    currentFm->size += sizeof(int);
	    memcpy(&field_list[fieldNo], a->field, sizeof(FMField));
	    fieldNo++;

	}

    }

    for (; fieldNo < (t->var_count + 1+altvarcount); fieldNo++)
    {
	field_list[fieldNo].field_type = NULL;
	field_list[fieldNo].field_name = NULL;
	field_list[fieldNo].field_offset = 0;
	field_list[fieldNo].field_size = 0;
    }

    format->field_list = field_list;


    currentFm->format->struct_size = currentFm->size;

    // create buffer to hold top-level (non-variable array) data.

    currentFm->buffer = (unsigned char *) malloc(currentFm->size);
    memset(currentFm->buffer, 0, currentFm->size);

    currentFm->sndCount = 0;

    return currentFm;
}




// copies buffer zeroing out arrays that havent been asked for
void* copy_buffer(void* buffer, int rank){
    char* temp = (char*)malloc(localWriteData->fm->size);
    memcpy(temp, buffer, localWriteData->fm->size);
    FMField *f = localWriteData->fm->format->field_list;
    while (f->field_name != NULL)
    {
        //if we wrote an array and didnt get asked for it, zero it out
        VarNode* a;
        if(!queue_contains(localWriteData->askedVars, f->field_name, rank)) {
            if((a=queue_contains(localWriteData->formatVars, f->field_name, -1)) 
	       && 
	       (a->dims != NULL)) {
                //fprintf(stderr, "rank %d won't send array %s\n", localWriteData->rank, f->field_name);
                VarNode* dim = a->dims;
                while(dim) {
                    //fprintf(stderr, "looking at %s\n", dim->varName);
                    FMField *f2 = localWriteData->fm->format->field_list;
                    while(f2->field_name != NULL) {
                        //fprintf(stderr, "compared to %s\n", f2->field_name);
                        if(strcmp(f2->field_name, dim->varName)==0) {
                            break;
                        }
                        f2++;
                    }
                    if(f2->field_name != NULL) {
                        //fprintf(stderr, "zero out dim %s\n", f2->field_name);
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




static int var_handler(CManager cm, void *vevent, void *client_data, attr_list attrs){
    Var_msg* msg = (Var_msg*) vevent;
    EVtake_event_buffer(cm, msg);
    fprintf(stderr, "rank %d <- var_msg : rank %d\n", localWriteData->rank, msg->rank);
    threaded_enqueue(&localWriteData->controlQueue, 
		     msg, 
		     VAR, 
		     localWriteData->controlMutex, 
		     &localWriteData->controlCondition);
    return 0;
}

static int flush_handler(CManager cm, void* vevent, void* client_data, attr_list attrs) {
    Flush_msg* msg = (Flush_msg*) vevent;
    EVtake_event_buffer(cm, msg);
    // Don't have to take buffer since we're only copying the rank onto the stack.
    fprintf(stderr, 
	    "rank %d <- flush data : rank %d\n", 
	    localWriteData->rank, msg->rank);
    threaded_enqueue(&localWriteData->controlQueue, 
		     msg, 
		     DATA_FLUSH, 
		     localWriteData->controlMutex, 
		     &localWriteData->controlCondition);
    return 0;
}

static int op_handler(CManager cm, void* vevent, void* client_data, attr_list attrs) {
    op_msg* msg = (op_msg*) vevent;
    EVtake_event_buffer(cm, msg);
    fprintf(stderr, "rank %d <- op_msg : rank %d type %d\n", localWriteData->rank, msg->process_id, msg->type);
    if(msg->type == 1) {
        //fprintf(stderr, "recieved open\n");
        threaded_enqueue(&localWriteData->controlQueue, msg, OPEN, localWriteData->controlMutex, &localWriteData->controlCondition);
    } else {
        //fprintf(stderr, "recieved close\n");
        threaded_enqueue(&localWriteData->controlQueue, NULL, CLOSE, localWriteData->controlMutex, &localWriteData->controlCondition);
    }
    return 0;
}

attr_list set_size_atom(attr_list attrs, int value) {
    atom_t dst_atom = attr_atom_from_string("fp_size");
    int size;
    if(!get_int_attr(attrs, dst_atom, &size)) {
        add_int_attr(attrs, dst_atom, value);
    }
    set_int_attr(attrs, dst_atom, value);
    return attrs;
}

attr_list set_dst_rank_atom(attr_list attrs, int value) {
    atom_t dst_atom = attr_atom_from_string("fp_dst_rank");
    int dst;
    if(!get_int_attr(attrs, dst_atom, &dst)) {
        add_int_attr(attrs, dst_atom, value);
    }
    set_int_attr(attrs, dst_atom, value);
    return attrs;
}


int control_thread(void* arg) {
    int rank = *((int*) arg);
    //fprintf(stderr, "entered control thread on rank %d\n", rank);
    QueueNode* controlMsg;
    QueueNode* dataNode;
    //fprintf(stderr, "waiting for control message\n");
    while(1) {
	if((controlMsg = threaded_dequeue(&localWriteData->controlQueue, 
					  localWriteData->controlMutex, 
					  &localWriteData->controlCondition))) {
	    if(controlMsg->type==VAR) {
		//add to localWriteData->askedVars
		Var_msg* varMsg = (Var_msg*) controlMsg->data;
		//fprintf(stderr, "recieved var message---------------------++++++%d\n",varMsg->rank);
		localWriteData->askedVars = add_var(localWriteData->askedVars, 
						    strdup(varMsg->var_name), 
						    NULL, 
						    varMsg->rank);
		EVreturn_event_buffer(localWriteData->cm,controlMsg->data);
	    } else if(controlMsg->type==DATA_FLUSH) { 
		//fprintf(stderr, "recieved data flush message\n");
		//make copy of buffer
		dataNode = threaded_peek(&localWriteData->dataQueue, 
					 localWriteData->dataMutex, 
					 &localWriteData->dataCondition);
		//fprintf(stderr, "peeked at dataNode %p\n", dataNode->data);
		Flush_msg* flushMsg = (Flush_msg*) controlMsg->data;
		//fprintf(stderr, "looking at flush msg\n");
		void* temp = copy_buffer(dataNode->data, flushMsg->rank);
		//add dst attr 
		//fprintf(stderr, "adding dst attr\n");
		localWriteData->attrs = set_dst_rank_atom(localWriteData->attrs, flushMsg->rank);
		//send data on multiqueue stone
		//fprintf(stderr, "submitting data\n");
		EVsubmit_general(localWriteData->dataSource, temp, data_free, localWriteData->attrs);
	    } else if(controlMsg->type==OPEN) {
		fprintf(stderr, "recieved open message\n");
                op_msg* open = (op_msg*) controlMsg->data;
                fprintf(stderr, "rank %d has step %d\n", localWriteData->rank, open->step);
                localWriteData->bridges[open->process_id].step = open->step;
                if(open->step < currentStep) {
                    fprintf(stderr, "error! recieved open for past step...\n");
                } else if (open->step == currentStep){
                    fprintf(stderr, "equal to step\n");
                    thr_mutex_lock(localWriteData->openMutex);
		    if(localWriteData->openCount==-1) localWriteData->openCount=0;
                    localWriteData->openCount++;  
                    fprintf(stderr, "opencount %d\n", localWriteData->openCount);      
		    thr_mutex_unlock(localWriteData->openMutex);
                    fprintf(stderr, "send ack\n");
                    op_msg* ack = (op_msg*) malloc(sizeof(op_msg));
                    ack->file_name = "hey";
                    ack->process_id = localWriteData->rank;
                    ack->step = currentStep;
                    ack->type = 2;
                    localWriteData->attrs = set_dst_rank_atom(localWriteData->attrs, open->process_id+1);
                    EVsubmit_general(localWriteData->opSource, ack, op_free, localWriteData->attrs);
                    fprintf(stderr, "continue\n");
                } else {
                    fprintf(stderr, "future step\n");
                }
            } else if(controlMsg->type==CLOSE) {
		//fprintf(stderr, "recieved close message\n");
		thr_mutex_lock(localWriteData->openMutex);
		localWriteData->openCount--;
                fprintf(stderr, "opencount %d\n", localWriteData->openCount);
		thr_mutex_unlock(localWriteData->openMutex);
		fprintf(stderr, "unlocked...\n");
                 if(localWriteData->openCount==0) {
		    threaded_dequeue(&localWriteData->dataQueue, 
                                    localWriteData->dataMutex, 
				     &localWriteData->dataCondition);
                    fprintf(stderr, "end of step %d\n", currentStep);
                    currentStep++;
                    //for all bridges if step == currentstep send ack
                    int i;
                    for(i=0; i<localWriteData->numBridges; i++) {
                      if(localWriteData->bridges[i].step==currentStep) {
                        localWriteData->openCount++;
                        op_msg* ack = (op_msg*) malloc(sizeof(op_msg));
                        ack->file_name = "hey";
                        ack->process_id = localWriteData->rank;
                        ack->step = currentStep;
                        ack->type = 2;
                        localWriteData->attrs = set_dst_rank_atom(localWriteData->attrs, i+1);
                        EVsubmit_general(localWriteData->opSource, ack, op_free, localWriteData->attrs);
                      }
                    }
		}
	    } else {
		fprintf(stderr, "unrecognized control message in control thread\n");
	    }
	}
    }
    //fprintf(stderr, "exiting control thread\n");
    return 0;
}

// Flexpath Functions

/*
 * Initializes flexpath write local data structures for client writes
 * - malloc space for global values
 * - store reference to new connection manager instance
 */

extern void 
adios_flexpath_init(const PairStruct *params, struct adios_method_struct *method) 
{
    setenv("CMSelfFormats", "1", 1);
    localWriteData = (FlexpathWriteData *) malloc(sizeof(FlexpathWriteData));
    if(!localWriteData) {
        adios_error(err_no_memory, "Cannot allocate memory for flexpath.");
    }
    memset(localWriteData, 0, sizeof(FlexpathWriteData));
    localWriteData->CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");

    attr_list listen_list = NULL;
    char * transport = NULL;
    transport = getenv("CMTransport");

    if(method!=NULL) {
        //fprintf(stderr, "todo: recieved non null method struct\n");
    }
    
    localWriteData->attrs = create_attr_list();
    gen_pthread_init();
    localWriteData->cm = CManager_create();

    if(transport == NULL){
	perr("transport is null\n");
	if(CMlisten(localWriteData->cm) == 0) {
	    fprintf(stderr, "error: unable to initialize connection manager.\n");
	    exit(1);
	}
    }else{
	perr("writer transport: %s\n", transport);
	listen_list = create_attr_list();
	add_attr(listen_list, localWriteData->CM_TRANSPORT, Attr_String, (attr_value)strdup(transport));
	CMlisten_specific(localWriteData->cm, listen_list);
    }

    // setup mutexs
    localWriteData->controlMutex = thr_mutex_alloc();
    localWriteData->dataMutex = thr_mutex_alloc();
    localWriteData->openMutex = thr_mutex_alloc();
    localWriteData->openCount = -1;
    // setup conditions
    localWriteData->controlCondition = CMCondition_get(localWriteData->cm, NULL);
    localWriteData->dataCondition = CMCondition_get(localWriteData->cm, NULL);

    // fork communications thread
    int forked = CMfork_comm_thread(localWriteData->cm);   
    if(!forked) {
         fprintf(stderr, "error forking comm thread\n");
    }


    //fprintf(stderr, "debug: flexpath write structures initialized\n");
    return;
}

/*
 * Sets up local data structure for a series of writes to an adios file.
 * - store reference to MPI_Comm and get MPI information
 */

void 
print_array(char* arr, int len) 
{
    //fprintf(stderr, "array:");
    int i;
    /*
    for (i = 0; i < len; i++) {
        fprintf(stderr, " 0x%x", arr[i]);
    }
    fprintf(stderr, "\n");
    */
}

extern int 
adios_flexpath_open(struct adios_file_struct *fd, struct adios_method_struct *method, void*comm) 
{ 
    //perr("opening\n");
    //char * filebase = "/tmp/work/jdayal3/titan/";
    char writer_info_filename[200];
    char writer_ready_filename[200];
    char reader_info_filename[200];
    char reader_ready_filename[200];
    
    /*
    sprintf(writer_info_filename, "%s", filebase);
    sprintf(writer_ready_filename, "%s", filebase);
    sprintf(reader_info_filename, "%s", filebase);
    sprintf(reader_ready_filename, "%s", filebase);
    */

    char *recv_buff = NULL;
    int i;
    char sendmsg[CONTACT_STR_LEN];
    //fprintf(stderr, "debug: entering adios_flexpath_open\n");


    if( fd == NULL || method == NULL) {
        //fprintf(stderr, "Bad input parameters\n");
        return -1;
    }
    /*
    MPI_Comm * group_comm = (MPI_Comm*)malloc(sizeof(MPI_Comm));
    adios_var_to_comm(fd->group->group_comm,
		      fd->group->adios_host_language_fortran,
		      comm,		      
		      group_comm);
    */
    localWriteData->mpiComm = (MPI_Comm*)comm;
    //localWriteData->mpiComm = group_comm;
    localWriteData->global_count = 0;
    localWriteData->sent_global_offsets = 0;
    // if evpath graph is not created and correct, exchange contact list and setup
    if(localWriteData->setupCorrect == 0) {
        // get some mpi information
        MPI_Comm_rank(*(localWriteData->mpiComm), &localWriteData->rank);
        MPI_Comm_size(*(localWriteData->mpiComm), &localWriteData->size);

        //if rank 0 allocate buffer
        if(localWriteData->rank == 0) {
            
            recv_buff = (char *) malloc(localWriteData->size*CONTACT_STR_LEN*sizeof(char));
        }
        
        // get contact string
        char * contact = attr_list_to_string(CMget_contact_list(localWriteData->cm));

	// allocate our stones
	localWriteData->multiStone = EValloc_stone(localWriteData->cm);
        localWriteData->sinkStone = EValloc_stone(localWriteData->cm);

        // send out contact string
	sprintf(&sendmsg[0], "%d:%s", localWriteData->multiStone, contact);
        MPI_Gather(sendmsg, 
		   CONTACT_STR_LEN, 
		   MPI_CHAR, recv_buff, 
		   CONTACT_STR_LEN, 
		   MPI_CHAR, 
		   0, 
		   *(localWriteData->mpiComm));

        //if rank 0 write contact info to file
        if(localWriteData->rank == 0) {
	    //fprintf(stderr, "rank 0 writting info\n");
	    sprintf(writer_info_filename, "%s_%s", fd->name, "writer_info.txt");
            FILE* writer_info = fopen(writer_info_filename,"w");
            for(i=0; i<localWriteData->size; i++) {
                fprintf(writer_info, "%s\n",&recv_buff[i*CONTACT_STR_LEN]); 
            }
            fclose(writer_info);
        }

        //poll file - race condition issues
	FILE* reader_ready = NULL;
	sprintf(reader_ready_filename, "%s_%s", fd->name, "reader_ready.txt");
        //fprintf(stderr, "polling reader_ready.txt file\n");
	while(!reader_ready){
	    reader_ready = fopen(reader_ready_filename,"r");
	}
	fclose(reader_ready);

        //read contact list
	//fprintf(stderr, "processing file\n");
	sprintf(reader_info_filename, "%s_%s", fd->name, "reader_info.txt");
        FILE* reader_info = fopen(reader_info_filename, "r");
	while(!reader_info){
	  reader_info = fopen(reader_info_filename, "r");
	}
        char in_contact[CONTACT_STR_LEN] = "";
       
        int numBridges = 0;
        int stone_num;
        //fprintf(stderr, "rank %d enters loop\n",localWriteData->rank);
        while(fscanf(reader_info, "%d:%s",&stone_num, in_contact)!=EOF){
            //fprintf(stderr, "rank %d recieved contact %s for stone %d\n", localWriteData->rank, in_contact, stone_num);
            //for each line create bridge stone
            localWriteData->bridges = realloc(localWriteData->bridges, sizeof(Stone) * (numBridges + 1));
	    //fprintf(stderr, "resized bridge array\n");
            attr_list contact_list = attr_list_from_string(in_contact);
            //fprintf(stderr, "generated contact list\n");
            localWriteData->bridges[numBridges].myNum = EVcreate_bridge_action(localWriteData->cm, contact_list, stone_num);
            //fprintf(stderr, "created bridge action\n");
            localWriteData->bridges[numBridges].step = 0;
            localWriteData->bridges[numBridges].theirNum = stone_num;
            localWriteData->bridges[numBridges].contact = strdup(in_contact);
            numBridges += 1;
        }
        //fprintf(stderr, "rank %d exits loop\n", localWriteData->rank);
        localWriteData->numBridges = numBridges;
        fclose(reader_info);
	// cleanup of reader file (writer is done with it).
	if(localWriteData->rank == 0){
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

	localWriteData->fm = setFormat(t, fields);
	print_format(localWriteData->fm->format);

        //generate multiqueue function that sends formats or all data based on flush msg
        FMStructDescList queue_list[] = {flush_format_list, 
					 format_format_list, 
					 var_format_list, 
					 op_format_list, 
					 evgroup_format_list,
					 data_format_list,					 
					 NULL};
        char* q_action_spec = create_multityped_action_spec(queue_list, multiqueue_action); 
        EVaction multi_action = EVassoc_multi_action(localWriteData->cm, 
						     localWriteData->multiStone, 
						     q_action_spec, 
						     NULL);
	localWriteData->formatSource = EVcreate_submit_handle(localWriteData->cm, 
							      localWriteData->multiStone, 
							      format_format_list);
        localWriteData->dataSource = EVcreate_submit_handle_free(localWriteData->cm, localWriteData->multiStone, localWriteData->fm->format, data_free,  NULL); 
        localWriteData->opSource = EVcreate_submit_handle_free(localWriteData->cm, localWriteData->multiStone, op_format_list, op_free,  NULL); 
	localWriteData->offsetSource = EVcreate_submit_handle(localWriteData->cm, 
							      localWriteData->multiStone, 
							      evgroup_format_list);
        EVassoc_terminal_action(localWriteData->cm, 
				localWriteData->sinkStone, 
				var_format_list, 
				var_handler, 
				NULL);
        EVassoc_terminal_action(localWriteData->cm, 
				localWriteData->sinkStone, 
				op_format_list, 
				op_handler, 
				NULL);
        EVassoc_terminal_action(localWriteData->cm, 
				localWriteData->sinkStone, 
				flush_format_list, 
				flush_handler, 
				NULL);

        //link multiqueue to sink
        EVaction_set_output(localWriteData->cm, 
			    localWriteData->multiStone, 
			    multi_action, 
			    0, 
			    localWriteData->sinkStone);

        //link up multiqueue ports to bridge stones
        int i=0;
        for(i=0; i<numBridges; i++) {
            //fprintf(stderr, "linking port %d to rank %d\n", i+1, i);
            EVaction_set_output(localWriteData->cm, 
				localWriteData->multiStone, 
				multi_action, 
				i+1, 
				localWriteData->bridges[i].myNum);
	}
	
	//store format id in multiqueue
	Format_msg *initial_format_msg = malloc(sizeof(Format_msg));
        FMContext my_context = create_local_FMcontext();	
        FMFormat my_format = register_data_format(my_context, localWriteData->fm->format);

        int id_len;
	char* temp = get_server_ID_FMformat(my_format, &id_len);
        
	/*for(i=0; i<id_len; i++) {
            temp[i]=temp[i]+1;
	    }*/
        initial_format_msg->format_id = temp;
        initial_format_msg->id_len = id_len;

	int rep_len;
	char* temp2 = get_server_rep_FMformat(my_format, &rep_len);
	//for(i=0; i<id_len; i++)
	//    temp2[i]=temp2[i]+1;
	initial_format_msg->rep_id = temp2;
	initial_format_msg->rep_id_len = rep_len;
        print_array(temp, id_len);
        CMsleep(localWriteData->cm, 1); 
	/* 
	 * this is a memory leak of initial_format_msg (and my_context, etc).
	 * Fix by changing to EVsubmit_general and using a free call.
	 */
	EVsubmit(localWriteData->formatSource, initial_format_msg, localWriteData->attrs);
        //fprintf(stderr, "notifing reader its ok to send msgs\n");
        CMsleep(localWriteData->cm, 1); 
	sprintf(writer_ready_filename, "%s_%s", fd->name, "writer_ready.txt");
        if(localWriteData->rank == 0) {
            FILE* writer_info = fopen(writer_ready_filename, "w");
            fprintf(writer_info, "ready");
            fclose(writer_info);
        }
        
        CMsleep(localWriteData->cm, 5); 
	//indicate evpath is setup correctly
	localWriteData->setupCorrect = 1;
        // fork control thread
        thr_thread_t forked_thread = thr_fork(control_thread, &localWriteData->rank);
        if(forked_thread) {
            //fprintf(stderr, "successfully forked control thread\n");
        } else {
            //fprintf(stderr, "error forking control thread\n");
        }
    //fprintf(stderr, "continuing\n");
    }

    return 0;
	
}


extern enum ADIOS_FLAG adios_flexpath_should_buffer (struct adios_file_struct * fd,struct adios_method_struct * method) {
    //fprintf(stderr, "debug: adios_flexpath_should_buffer\n");
  return adios_flag_unknown;
}


//  writes data to multiqueue
extern void adios_flexpath_write(struct adios_file_struct *fd, struct adios_var_struct *f, void *data, struct adios_method_struct *method) {
    //fprintf(stderr, "debug: adios_flexpath_write\n");
    
    //fprintf(stderr, "submitting var %s\n", f->name);
    
    struct _fm_structure* fm = localWriteData->fm;

    if (fm == NULL)
    {
	//fprintf(stderr, "group or fm is null - improperly initialized\n");
	return;

    }
    
    FMFieldList flist = fm->format->field_list;
    FMField *field = NULL;
    char *fixedname = find_fixed_name(fm, f->name);
    field = internal_find_field(fixedname, flist);
    //fprintf(stderr, "found field %s\n", field->field_name);
    if (field != NULL) {
	if (!f->dimensions) {
	    //scalar quantity
            //fprintf(stderr, "copying scalar value\n");
	    if (data) {
		//why wouldn't it have data?
		memcpy(&fm->buffer[field->field_offset], data, field->field_size);

		//scalar quantities can have AltNames also so assign those
		if(field->field_name != NULL) {
					
		    DimNames *d = NULL;
		    for (d = fm->dimList.lh_first; d != NULL; d = d->entries.le_next) {
			if (!strcmp(d->name, field->field_name)) {
			    //matches
			    //check if there are AltNames
			    AltName *a = NULL;
			    for (a = d->altList.lh_first; a != NULL; a = a->entries.le_next) {
				//use the AltName field to get the data into the buffer
				memcpy(&fm->buffer[a->field->field_offset], 
				       data, 
				       a->field->field_size);
                		//int *testingint = (int*)&fm->buffer[a->field->field_offset];
		        	//fprintf(stderr, "writing %s to %s at %d %d\n", f->name, a->name, a->field->field_offset, (int)*testingint);
			    }
			}
		    }
		}
	    } else {
		//fprintf(stderr, "no data for  scalar %s\n", f->name);
	    }
	} else {
	    //vector quantity
	    if (data)
	    {
                //fprintf(stderr, "copying vector pointer\n");
		//we just need to copy the pointer stored in f->data
                // calculate size
                memcpy(&fm->buffer[field->field_offset], &data, sizeof(void *));

	    } else {
		//fprintf(stderr, "no data for vector %s\n", f->name);
	    }
	}
    }
    //fprintf(stderr, "successfully copied data to buffer\n");
}

extern void 
adios_flexpath_close(struct adios_file_struct *fd, struct adios_method_struct *method) 
{
    //fprintf(stderr, "debug: adios_flexpath_close\n");
    //fprintf(stderr, "enqueueing data\n");   
    
    //no copy
    //void* buffer =  localWriteData->fm->buffer;
    
    //copy
    void* buffer = malloc(localWriteData->fm->size);

    struct adios_group_struct * g2 = fd->group;
    struct adios_var_struct * fields = g2->vars;
    while(fields) {
        
        //fprintf(stderr, "field name: %s\n", fields->name);
        if(fields->dimensions) {
            //fprintf(stderr, "field is an array\n");
            struct adios_dimension_struct* dims = fields->dimensions;
            //fprintf(stderr, "field dims: %p\n", dims);
    
            int total_size = 1;
            //for each dimension
            while(dims) {    
                struct adios_var_struct* temp = adios_find_var_by_id(g2->vars, dims->dimension.id);            
                int size = *(int*)temp->data;
                //fprintf(stderr, "dim %s size %d\n", temp->name, *(int*)temp->data);
                total_size *= size;
                dims = dims->next;
            }		
            FMFieldList flist = localWriteData->fm->format->field_list;
            FMField *field = NULL;
            char *fixedname = find_fixed_name(localWriteData->fm, fields->name);
            field = internal_find_field(fixedname, flist);
            //fprintf(stderr, "field offset %d size %d\n", field->field_offset, field->field_size);

            total_size*=field->field_size;
            // malloc size
            //fprintf(stderr, "field %s field size %d\n", fields->name, total_size);
            void* pointer_data_copy = malloc(total_size);
            // while null
            while(pointer_data_copy==NULL) { 
                sleep(1);
                void* pointer_data_copy = malloc(total_size);
                //block
            }
                
            // memcpy data
            char* cur_offset = (char*)&localWriteData->fm->buffer[field->field_offset];
            //fprintf(stderr, "pointer %p\n", cur_offset); 
            void* aptr8 = (void*)(*((unsigned long*)cur_offset));
            //fprintf(stderr, "copying data to %p from %p of size %d\n", pointer_data_copy, aptr8, total_size);
            memcpy(pointer_data_copy, aptr8, total_size);
            // memcpy pointer
            memcpy(&localWriteData->fm->buffer[field->field_offset], &pointer_data_copy, sizeof(void *));
        }    
        fields = fields->next;
    }

    
    memcpy(buffer, localWriteData->fm->buffer, localWriteData->fm->size);

    //fprintf(stderr, "got buffer %p\n", buffer);
    threaded_enqueue(&localWriteData->dataQueue, 
		     buffer, 
		     DATA_BUFFER, 
		     localWriteData->dataMutex, 
		     &localWriteData->dataCondition);
    //fprintf(stderr, "sucessfully enqueued\n");
    localWriteData->fm->sndCount++;
    // now gather offsets and send them via MPI to root
    struct adios_group_struct * g = fd->group;
    struct adios_var_struct * list = g->vars;

    if(localWriteData->global_count > 0 && !localWriteData->sent_global_offsets){	
	// process local offsets here	
	int num_gbl_vars = 0;
        global_var * gbl_vars = NULL;
	int num_vars = 0;
	int myrank = localWriteData->rank;
	int commsize = localWriteData->size;

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
			   0, *localWriteData->mpiComm);

		MPI_Gather(local_dimensions, num_local_offsets, MPI_INT, 
			   all_local_dims, num_local_offsets, MPI_INT,
			   0, *localWriteData->mpiComm);

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
		    //perr("\n\n\n\t\tnoffset_structs: %d\n", 
		    // gbl_vars[num_gbl_vars - 1].noffset_structs);
		    gbl_vars[num_gbl_vars - 1].offsets = ostruct;
		    int i;			   
		    i = 0;
		    //perr("\t\t\tall offsets for var: %s\n", list->name);
		    while(i<commsize * num_local_offsets){
			int j;
			/*
			for(j=0; j<num_local_offsets;j++){
			    perr("\t\t\t%d ", all_offsets[i]);  
			    i++;
			}			
			perr("\n");				
			*/
		    }
		    //perr("\n");
		    //perr("\t\t\tall local_dims for var: %s\n", list->name);
		    i = 0;
		    /*
		    while(i<commsize * num_local_offsets){
			int j;				
			for(j=0; j<num_local_offsets;j++){
			    perr("\t\t\t%d ", all_local_dims[i]);  
			    i++;
			}			
			perr("\n");
			
		    }
		    perr("\n");
		    */
		}		
	    }
	    list=list->next;
	    //okay, no more dims, now send it to root.
	}
	if(myrank == 0){
	    int i;
	    /*
	    for(i=0; i<num_gbl_vars; i++){
		perr("global_var: %s has local offsets\n", gbl_vars[i].name);
	    }
	    */
	    evgroup * gp = (evgroup*)malloc(sizeof(evgroup));
	    gp->num_vars = num_gbl_vars;
	    //perr("num global vars %d\n", num_gbl_vars);
            gp->vars = gbl_vars;
	    localWriteData->gp = gp;
	    localWriteData->attrs = set_size_atom(localWriteData->attrs, localWriteData->size);
	    //perr("size:%d\n\n\n", localWriteData->size);
            EVsubmit(localWriteData->offsetSource, gp, localWriteData->attrs);
	    /*
	    threaded_enqueue(&localWriteData->offsetQueue, 
			     gp, 
			     OFFSET_MSG, 
			     localWriteData->offsetMutex, 
			     &localWriteData->offsetCondition);
	    */
	}
	localWriteData->sent_global_offsets = 1;
    }
    //fprintf(stderr, "exiting close\n");
}

extern void adios_flexpath_finalize(int mype, struct adios_method_struct *method) {
    //fprintf(stderr, "debug: adios_flexpath_finalize\n");
    while(localWriteData->dataQueue!=NULL) {
        //fprintf(stderr, "rank %d data still in data queue\n", localWriteData->rank);
        sleep(1);
    }
}

extern void adios_flexpath_end_iteration(struct adios_method_struct *method) {
    //fprintf(stderr, "debug: adios_flexpath_end_iteration\n");
    struct _fm_structure *fm;
    FlexpathWriteData *mdata = (FlexpathWriteData *) method->method_data;
    fm = mdata->fm;

    if (fm == NULL)
	return;

    mdata->cycleId = 0;

}

extern void adios_flexpath_start_calculation(struct adios_method_struct *method) {
    //fprintf(stderr, "debug: adios_flexpath_start_calculation\n");
    struct _fm_structure *fm;
    FlexpathWriteData *mdata = (FlexpathWriteData *) method->method_data;
    fm = mdata->fm;
  
    if (fm == NULL)
	return;
  
}

extern void adios_flexpath_stop_calculation(struct adios_method_struct *method) {
    //fprintf(stderr, "debug: adios_flexpath_stop_calculation\n");
    struct _fm_structure *fm;
    FlexpathWriteData *mdata = (FlexpathWriteData *) method->method_data;
    fm = mdata->fm;
  
    if (fm == NULL)
	return;

    mdata->cycleId++;
}

extern void adios_flexpath_get_write_buffer(struct adios_file_struct *fd,struct adios_var_struct *f, uint64_t *size, void **buffer, struct adios_method_struct *method) {
    //fprintf(stderr, "debug: adios_flexpath_get_write_buffer\n");
  //MDW: What is this function for?????
  //fprintf(stderr, "adios_flexpath_write_get_buffer: flexpath disabled, no portals support\n");
}

void adios_flexpath_read(struct adios_file_struct *fd, struct adios_var_struct *f, void *buffer, uint64_t buffer_size, struct adios_method_struct *method) {
    // not required by write api
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
