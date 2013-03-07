/*
    read_flexpath.c
    
    Originally copied from read_datatap.c
    Goal: to create evpath io connection layer in conjunction with 
    write/adios_flexpath.c

*/





// system libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/queue.h>
#include <sys/socket.h>
#include <sys/times.h>
#include <netinet/in.h>
#include <sys/time.h>
#include <sys/uio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>

// evpath libraries
#include <ffs.h>
#include <atl.h>
#include <evpath.h>

// local libraries
#include "config.h"
#include "public/adios.h"
#include "public/adios_read_v2.h"
#include "core/adios_read_hooks.h"
#include "public/adios_error.h"
#include "core/globals.h"
#include "public/flexpath.h"

// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif


#define MAX(a, b) ((a<b)?(b):(a))
#define MIN(a, b) ((a>b)?(b):(a))

/*
 * Contains start & counts for each dimension for each writer.
 */
typedef struct _array_displ
{
    int writer_rank;
    int ndims;
    int * start;
    int * count;    
}array_displacements;

typedef struct _bridge_info
{
    EVstone bridge_stone;
    EVsource flush_source;
    EVsource var_source;
    EVsource op_source;
    int their_num;
    char * contact;
    int created;
    int step;
}bridge_info;

typedef struct _flexpath_var_chunk
{
    int has_data;
    int rank;
    void *data;
    uint64_t *local_bounds; // nodims
    uint64_t *global_bounds; // ndims
    uint64_t *global_offsets; // ndims
    struct _flexpath_var_chunk *next;
} flexpath_var_chunk, *flexpath_var_chunk_p;

typedef struct _flexpath_var_info
{
    int id;
    char *varname;
    char *varpath;
    enum ADIOS_DATATYPES type;
    //add selector typ
    uint64_t data_size; // type size, not arrays size
    int time_dim; // -1 means no time dimension
    int ndims;
    int * dims; // ndims size (if ndims>0)
    uint64_t array_size; // not relevant for scalars
    int num_chunks;
    int was_scheduled;
    flexpath_var_chunk *chunks;
    int num_displ;
    array_displacements * displ;
    const ADIOS_SELECTION * sel;
    struct _flexpath_var_info *next;
} flexpath_var_info, *flexpath_var_info_p;

typedef struct _flexpath_file_data
{
    char * file_name;
    char * group_name; // assuming one group per file right now.

    EVstone ctrl_stone;
    EVstone split_stone; // to be added to ctrl_stone
    EVaction split_action;
    EVstone data_stone;

    MPI_Comm comm;
    int rank;
    int size;
    int valid;
    struct _flexpath_file_data * next;
    int num_bridges;
    bridge_info *bridges;
    FMFormat current_format;
    int polling;
    int num_vars;
    flexpath_var_info * var_list;
    int num_gp; // for array distribution.
    int valid_evgroup;
    evgroup * gp;

    int* sendees;
    int num_sendees;

} flexpath_file_data;


flexpath_file_data *
new_flexpath_file_data(const char * fname);

flexpath_var_info*
new_flexpath_var_info(const char * varname, int id, uint64_t data_size);

flexpath_var_info*
new_flexpath_var_info(const char * varname, int id, uint64_t data_size)
{
    flexpath_var_info * var = malloc(sizeof(flexpath_var_info));
    if(var == NULL){
	perr("Error creating new var: %s\n", varname);
	exit(1);
    }
    
    var->varname = strdup(varname);
    var->id = id;
    var->data_size = data_size;
    var->chunks = NULL;
    var->sel = NULL;
    var->dims = NULL;
    var->displ = NULL;
    var->was_scheduled = 0;    
    var->time_dim = 0;
    var->ndims = 0;
    var->next = NULL;
    return var;
}


flexpath_file_data*
new_flexpath_file_data(const char * fname)
{
    flexpath_file_data * fp = malloc(sizeof(flexpath_file_data));
    if(fp == NULL){
	perr("Cannot create data for new file.\n");
	exit(1);
    }
    fp->file_name = strdup(fname);
    fp->group_name = NULL;    
    fp->next = NULL;
    fp->var_list = NULL;
    fp->gp = NULL;
    fp->bridges = NULL;
    fp->current_format = NULL;
    
    fp->valid = 0;
    fp->num_bridges = 0;
    fp->num_gp = 0;
    fp->valid_evgroup = 0;
    fp->polling = 1;
    fp->num_vars = 0;
    fp->sendees = NULL;
    fp->num_sendees = 0;    
    return fp;        
}

flexpath_file_data * file_data_list = NULL;

typedef struct _local_read_data
{
    // MPI stuff
    MPI_Comm fp_comm;
    int fp_comm_rank;
    int fp_comm_size;

    // EVPath stuff
    CManager fp_cm;
    EVstone ctrl_stone;
    EVstone data_stone;
    atom_t CM_TRANSPORT;

    // server state
    int fp_server_ready;
    int num_io_dumps;
    // TODO: timestep

} flexpath_read_data, *flexpath_read_data_p;


static int compare_var_name(const char* varname, const flexpath_var_info *v);
// this sructure holds all global data for flexpath read  methods
flexpath_read_data* fp_read_data = NULL;
int ackCondition;

#define VAR_BITMAP_SIZE 16

ADIOS_VARINFO*
convert_file_info(flexpath_var_info * current_var,
		  ADIOS_VARINFO * v,
		  const char* varname,
		  const ADIOS_FILE* gp);

flexpath_var_info *
find_fp_var(flexpath_var_info * var_list, const char * varname)
{
    while(var_list){
	if(!compare_var_name(varname, var_list)){
	    return var_list;
	}
	else
	    var_list = var_list->next;
    }
    return NULL;
}

// compare used-providd varname with the full path name of variable v
// return zero if matches and non-zero otherwise
static int
compare_var_name (const char *varname, const flexpath_var_info *v)
{
    if (varname[0] == '/') { // varname is full path
        char fullpath[256];
        if(!strcmp(v->varpath, "/")) {
            sprintf(fullpath, "/%s", v->varname);
        }
        else {
            sprintf(fullpath, "%s/%s", v->varpath, v->varname);
        }
        return strcmp(fullpath, varname);
    }
    else { // varname doesn't include path
        return strcmp(v->varname, varname);
    }
}

global_var* 
find_gbl_var(global_var * vars, char * name, int num_vars)
{
    global_var * retvar = NULL;
    int i;
    for(i=0; i<num_vars; i++){
	if(!strcmp(vars[i].name, name))
	    return &vars[i];
    }
    return retvar;
}

static FMField
*find_field (const char *name, const FMFieldList flist)
{
    FMField *f = flist;
    while (f->field_name != NULL)
    {
        if(!strcmp(name, f->field_name))
            return f;
        else
            f++;
    }
    return NULL;
}



static int op_msg_handler(CManager cm, void *vevent, void *client_data, attr_list attrs) {
    op_msg* msg = (op_msg*)vevent;
    //fprintf(stderr, "recieved op_msg type %d step %d\n", msg->type, msg->step);
    if(msg->type==2) {
        //fprintf(stderr, "signal ackCondition\n");
        CMCondition_signal(fp_read_data->fp_cm, ackCondition);
        ackCondition = CMCondition_get(fp_read_data->fp_cm, NULL);
    }
    return 0;
}

/*
 * Should only be invoked from rank 0.  might need a better way to go about this.
 */
static int
group_msg_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    //perr("debug: group_msg_handler called %d.\n", file_data_list->rank);
    EVtake_event_buffer(fp_read_data->fp_cm, vevent);
    evgroup * msg = (evgroup*)vevent;
    flexpath_file_data * fp = (flexpath_file_data*)client_data;
    fp->gp = msg;
    fp->valid_evgroup = 1;
    global_var * vars = msg->vars;
    int num_vars = msg->num_vars;


    /*
    int i;
    
    for(i=0; i<num_vars; i++){
	int noffset_structs = vars[i].noffset_structs;
	perr("noffset_structs: %d\n", noffset_structs);
	int j;
	for(j=0; j<noffset_structs; j++){
	    //perr("i: %d, j: %d\n", i, j);
	    int total_offsets = vars[i].offsets[j].total_offsets;
	    int k;
	    perr("offsets for var: %s\n", vars[i].name);
	    for(k=0; k<total_offsets; k++){
		perr("%d\t", vars[i].offsets[j].local_offsets[k]);
	    }
	    perr("\n");
	    perr("local_dims for var: %s\n", vars[i].name);
	    for(k=0; k<total_offsets; k++){
		perr("%d\t", vars[i].offsets[j].local_dimensions[k]);
	    }
	}
    }  
    */  

    return 0;

}


void 
print_int_arr(char * tag, int * arr, int count)
{
    int i;
    perr("%s: ", tag);
    for(i=0; i<count; i++){
	perr("%d ", arr[i]);
    }
    perr("\n");
}

array_displacements*
find_displacement(array_displacements* list, int rank, int num_displ){
    int i;
    for(i=0; i<num_displ; i++){
	if(list[i].writer_rank == rank)
	    return &list[i];	
    }
    return NULL;
}

int
linearize_displ(int * offset, int * sizes, int ndim, int data_size)
{
    //print_int_arr("linearize_displ offsets", offset, ndim);
    //print_int_arr("linearize_displ sizes", sizes, ndim);
    int i;
    int retval = 0;
    for(i = 0; i<ndim - 1; i++){
	retval += (offset[i] * sizes[i+1])*data_size;       
    }
    retval+=offset[ndim-1]*data_size;
    //perr("\t\t\tretval: %d\n\n\n", retval);
    return retval;
}


void 
copyoffsets(int dim, // dimension index
	    int ndims, // number of dimensions
	    int data_size, // data size
	    int* disp_start, // start array from array_displacements struct
	    int* disp_count, // count array from array_displacements struct
	    int* writer_count, // local dimensions from all writers; from offset_struct
	    uint64_t* reader_count, // the count field from reader's selector
	    char* from, 
	    char* to) 
{
    /*
    perr("\n\n\ncopyoffsets call:\n");
    print_int_arr("disp_start", disp_start, ndims);
    print_int_arr("disp_count", disp_count, ndims);
    print_int_arr("writer_count", writer_count, ndims);
    print_int_arr("reader_count", (int*)reader_count, ndims);
    */
    if(dim==ndims-1) {
        int* reader_count_copy = (int*)malloc(sizeof(int)*ndims);
        int i=0;
        for(i=0; i<ndims; i++){
            reader_count_copy[i]=reader_count[i];
        }
        int s = linearize_displ(disp_start, writer_count, ndims, data_size);
        int e = linearize_displ(disp_start, reader_count_copy, ndims, data_size);
        //perr("copying %d from %d to %d\n", disp_count[dim]*data_size, s, e);
	free(reader_count_copy);
        memcpy(to, from+s,  data_size*disp_count[ndims-1]);
    } else {
        int i;
        for(i=0; i<disp_count[dim]; i++) {
            int* disp_startcpy = malloc(sizeof(int)*ndims);
            memcpy(disp_startcpy, disp_start, sizeof(int)*ndims);
            disp_startcpy[dim] += i;
            copyoffsets(dim+1, 
			ndims, 
			data_size, 
			disp_start,
			disp_count, 
			writer_count, 
			reader_count, 
			from, 
			to);
        }
    }
}

/*
 * gets the data and puts it in the appropriate flexpath_var_struct.  Need to do a bit different
 * things if it's a scalar vs. an array.
 * Format for this is gathered from the file_data_list->my_format field
 */
static int
data_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    //perr("Data handler called.\n");
    int rank;
    flexpath_file_data * file_data = (flexpath_file_data*)client_data;
    char * buffer = (char*)vevent;

    get_int_attr(attrs, attr_atom_from_string(FP_RANK_ATTR_NAME), &rank);
    FMFormat format = file_data->current_format;

    if(!file_data_list->current_format){
	adios_error(err_file_open_error, "file not opened correctly.  Format not specified.\n");
    }

    FMStructDescList struct_list = FMcopy_struct_list(format_list_of_FMFormat(format));
    FMField *f = struct_list[0].field_list;
    char * curr_offset = NULL;
    int i = 0, l =0;
    int j=0;
    while(f->field_name != NULL){		
	curr_offset = &buffer[f->field_offset];
        char atom_name[200] = "";
	flexpath_var_info * var = find_fp_var(file_data->var_list, strdup(f->field_name));
	
	if(!var){
	    adios_error(err_file_open_error,
			"file not opened correctly.  var does not match format.\n");
	    return -1;
	}
        strcat(atom_name, f->field_name);
        strcat(atom_name, "_");
        strcat(atom_name, FP_NDIMS_ATTR_NAME);
        int num_dims;
        int i;
        get_int_attr(attrs, attr_atom_from_string(strdup(atom_name)), &num_dims);
	// scalar
	if(num_dims == 0){
	    flexpath_var_chunk * curr_chunk = &var->chunks[0];
	    curr_chunk->global_offsets = NULL;
	    curr_chunk->global_bounds = NULL;
	    curr_chunk->local_bounds = NULL;
	    curr_chunk->data = malloc(var->data_size);
	    memcpy(curr_chunk->data, curr_offset, var->data_size);
	    curr_chunk->has_data = 1;
	    // else it's an array
	}else{
            if(var->sel == NULL)
	    {// var hasn't been scheduled yet.  		
	    }
	    else if(var->sel->type == ADIOS_SELECTION_WRITEBLOCK){
		var->ndims = num_dims;
		var->dims = (int*)malloc(sizeof(int)*num_dims);
		if(var->was_scheduled == 1){
		    var->array_size = var->data_size;
		    for(i=0; i<num_dims; i++){
			char* dim = malloc(200*sizeof(char));
			atom_name[0] ='\0';
			strcat(atom_name, f->field_name);
			strcat(atom_name, "_");
			strcat(atom_name, FP_DIM_ATTR_NAME);
			strcat(atom_name, "_");
			char dim_num[10] = "";
			sprintf(dim_num, "%d", i+1);
			strcat(atom_name, dim_num);
			get_string_attr(attrs, attr_atom_from_string(atom_name), &dim);

			FMField * temp_f = find_field(dim, f);
			if(!temp_f){
			    adios_error(err_invalid_varname,
					"Could not find fieldname: %s\n",
					dim);
			}
			else{
			    // since it's a dimension, field size should be int.
			    char * temp_offset = &buffer[temp_f->field_offset];
			    int * temp_data = (int*)malloc(f->field_size);
			    memcpy(temp_data, temp_offset, f->field_size);
			    var->dims[i] = *temp_data;
			    var->array_size = var->array_size * var->dims[i];
			}
		    }
		    void* aptr8 = (void*)(*((unsigned long*)curr_offset));
		    memcpy(var->chunks[0].data, aptr8, var->array_size);
		}

	    }
	    else if(var->sel->type == ADIOS_SELECTION_BOUNDINGBOX){
		//perr("\t\t\tvar type: %d\n\n\n", (int)var->type);
		int i;
                global_var* gv = find_gbl_var(file_data_list->gp->vars, 
					      var->varname, 
					      file_data_list->gp->num_vars);
                int * writer_count = gv->offsets[0].local_dimensions;
                uint64_t * reader_count = var->sel->u.bb.count;
		array_displacements * disp = find_displacement(var->displ, 
							       rank, 
							       var->num_displ);
		void* aptr8 = (void*)(*((unsigned long*)curr_offset));
		double * temp = (double*)curr_offset;
		/*
		perr("%s, %d first value: %f\n", var->varname, (int)var->data_size, temp[2]);
		perr("copying offsets.\n");
		perr("\t\t\tndims: %d\n", disp->ndims);
		print_int_arr("\t\t\tdisp->start", disp->start, disp->ndims);
		print_int_arr("\t\t\tdisp->count", disp->count, disp->ndims);
		print_int_arr("\t\t\twriter_count", writer_count, disp->ndims);
		print_int_arr("\t\t\treader_count", (int*)reader_count, disp->ndims);		
                perr("\t\t\tcurent_offset %p", var->chunks[0].data);
		*/
                copyoffsets(0, 
			    disp->ndims,
			    f->field_size,
			    disp->start,
			    disp->count, 
			    writer_count, 
			    reader_count, 
			    (char*)aptr8, 
			    (char*)var->chunks[0].data);	
	    }
	}
        j++;
        f++;
    }
    file_data_list->polling = 0;
    return 0;
}

static int
format_handler(CManager cm, void *vevent, void *client_data, attr_list attrs) {
    Format_msg* msg = (Format_msg*)vevent;
    flexpath_file_data * fp = (flexpath_file_data*)client_data;
    char* arr = msg->format_id;
    char* rep_id = msg->rep_id;
    int rep_id_len = msg->rep_id_len;
    int len = msg->id_len;
    int i;

    FMContext my_context = create_local_FMcontext();
    if(my_context!=NULL) {
        //FMFormat my_format = FMformat_from_ID(my_context, arr);
	FMFormat my_format = load_external_format_FMcontext(my_context, arr, len, rep_id);
	if(!my_format)
	{
	    adios_error(err_file_open_error,
			"Could not get FMFormat from format server.");
	    return err_file_open_error;
	}
	fp->current_format = my_format;
    }
    fp->polling = 0;
    return 0;
}

/*
 * Initializes flexpath read structures for a client read
 * - malloc space for global values
 * - store reference to MPI_Comm and get MPI information
 * - store reference to a new evpath connection manager instance
 */
int
adios_read_flexpath_init_method (MPI_Comm comm, PairStruct* params)
{
    setenv("CMSelfFormats", "1", 1);
    fp_read_data = (flexpath_read_data *) malloc(sizeof(flexpath_read_data));     
    if(!fp_read_data) {
        adios_error(err_no_memory, "Cannot allocate memory for flexpath.");
        return -1;
    }
    memset(fp_read_data, 0, sizeof(flexpath_read_data));
    
    fp_read_data->CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    attr_list listen_list = NULL;
    char * transport = NULL;
    transport = getenv("CMTransport");

    // setup MPI stuffs
    fp_read_data->fp_comm = comm;
    MPI_Comm_size(fp_read_data->fp_comm, &(fp_read_data->fp_comm_size));
    MPI_Comm_rank(fp_read_data->fp_comm, &(fp_read_data->fp_comm_rank));

    // setup connection manager
    gen_pthread_init();
    fp_read_data->fp_cm = CManager_create();
    //perr("after cmanager_create\n");
    if(transport == NULL){
	if(CMlisten(fp_read_data->fp_cm) == 0) {
	    perr( "error: unable to initialize connection manager.\n");
	}
    }else{
	perr("reader transport: %s\n", transport);
	listen_list = create_attr_list();
	add_attr(listen_list, fp_read_data->CM_TRANSPORT, Attr_String, 
		 (attr_value)strdup(transport));
	CMlisten_specific(fp_read_data->fp_cm, listen_list);
    }
    if(CMfork_comm_thread(fp_read_data->fp_cm)) {/*perr( "forked\n");*/}
    return 0;
}


void build_bridge(bridge_info* bridge) {
    attr_list contact_list = attr_list_from_string(bridge->contact);

    bridge->bridge_stone =
        EVcreate_bridge_action(fp_read_data->fp_cm,
            contact_list,
            (EVstone)bridge->their_num);

    bridge->flush_source =
        EVcreate_submit_handle(fp_read_data->fp_cm,
            bridge->bridge_stone,
            flush_format_list);

    bridge->var_source =
	EVcreate_submit_handle(fp_read_data->fp_cm,
	    bridge->bridge_stone,
	    var_format_list);

    bridge->op_source =
	EVcreate_submit_handle(fp_read_data->fp_cm,
	    bridge->bridge_stone,
	    op_format_list);

    bridge->created = 1;

}

/*
 * Still have work to do here.  will have to move code from open_file, and
 * change it so that we can support the timeouts and lock_modes.
 */
ADIOS_FILE*
adios_read_flexpath_open_stream(const char * fname,
				MPI_Comm comm,
                                enum ADIOS_LOCKMODE lock_mode,
				float timeout_sec)
{
    return adios_read_flexpath_open_file(fname, comm);
}


/*
 * Sets up local data structure for series of reads on an adios file
 * - create evpath graph and structures
 * -- create evpath control stone (outgoing)
 * -- create evpath data stone (incoming)
 * -- rank 0 dumps contact info to file
 * -- create connections using contact info from file
 */
ADIOS_FILE*
adios_read_flexpath_open_file(const char * fname, MPI_Comm comm)
{
    //perr( "debug: entering adios_read_flexpath_fopen\n");
    //search linked list for connection, if not there add connection send open, otherwise use connection
    // connection information for this file does not yet exist.
    // establish graph for this file.
    ADIOS_FILE *fp = malloc(sizeof(ADIOS_FILE));    
    if(!fp){
	adios_error (err_no_memory, "Cannot allocate memory for file info.\n");
	return NULL;
    }    
    if(file_data_list == NULL){	
	ackCondition = CMCondition_get(fp_read_data->fp_cm, NULL);
	
	adios_errno = 0;
 	file_data_list = new_flexpath_file_data(fname);
	file_data_list->data_stone = EValloc_stone(fp_read_data->fp_cm);	
	file_data_list->comm = comm;

#ifndef _NOMPI      
	MPI_Comm_size(file_data_list->comm, &(file_data_list->size));
	MPI_Comm_rank(file_data_list->comm, &(file_data_list->rank));
#endif

	EVassoc_terminal_action(fp_read_data->fp_cm,
				file_data_list->data_stone,
				op_format_list,
				op_msg_handler,
				file_data_list);
	
        EVassoc_terminal_action(fp_read_data->fp_cm,
				file_data_list->data_stone,
				format_format_list,
				format_handler,
				file_data_list);

	EVassoc_terminal_action(fp_read_data->fp_cm,
				file_data_list->data_stone,
				evgroup_format_list,
				group_msg_handler,
				file_data_list);

	char writer_ready_filename[200];
	char writer_info_filename[200];
	char reader_ready_filename[200];
	char reader_info_filename[200];
	
	sprintf(reader_ready_filename, "%s_%s", fname, READER_READY_FILE);
	sprintf(reader_info_filename, "%s_%s", fname, READER_CONTACT_FILE);
	sprintf(writer_ready_filename, "%s_%s", fname, WRITER_READY_FILE);
	sprintf(writer_info_filename, "%s_%s", fname, WRITER_CONTACT_FILE);
	
	char * string_list;
	char data_contact_info[50];
	string_list = attr_list_to_string(CMget_contact_list(fp_read_data->fp_cm));
	sprintf(&data_contact_info[0], "%d:%s", file_data_list->data_stone, string_list);
	char * recvbuf;
	if(file_data_list->rank == 0){	
	    recvbuf = (char*)malloc(sizeof(char)*50*(file_data_list->size));
	}

#ifndef _NOMPI
	MPI_Gather(data_contact_info, 50, MPI_CHAR, recvbuf,
		   50, MPI_CHAR, 0, file_data_list->comm);
#endif

	if(file_data_list->rank == 0)
	{
	    // print our own contact information

	    FILE * fp_out = fopen(reader_info_filename, "w");
            int i;
	    if(!fp_out)
	    {
		adios_error(err_file_open_error,
			    "File for contact info could not be opened for writing.\n");
		exit(1);
	    }
            for(i=0; i<file_data_list->size; i++) {
	        fprintf(fp_out,"%s\n", &recvbuf[i*50]);
            }
	    fclose(fp_out);
            fp_out = fopen(reader_ready_filename, "w");
            fprintf(fp_out, "ready");
            fclose(fp_out);
        }

	FILE * read_ready = fopen(reader_ready_filename, "w");
	fprintf(read_ready, "ready");
	fclose(read_ready);
	//may need to switch to rank 0 and mpi broadcast
        FILE * fp_in = fopen(writer_ready_filename,"r");
	while(!fp_in) {
            CMsleep(fp_read_data->fp_cm, 1);
	    fp_in = fopen(writer_ready_filename, "r");
        }
        fclose(fp_in);
        fp_in = fopen(writer_info_filename, "r");
        while(!fp_in)
            fp_in = fopen(writer_info_filename, "r");

	char in_contact[50] = "";
	file_data_list->bridges = malloc(sizeof(bridge_info));
	int num_bridges = 0;
	int their_stone;

        // change to read all numbers, dont create stones, turn bridge array into linked list
        while(fscanf(fp_in, "%d:%s", &their_stone, in_contact) != EOF){	
	    file_data_list->bridges = realloc(file_data_list->bridges,
					      sizeof(bridge_info) * (num_bridges+1));
	    file_data_list->bridges[num_bridges].their_num = their_stone;
            file_data_list->bridges[num_bridges].contact = strdup(in_contact);
            file_data_list->bridges[num_bridges].created = 0;
            file_data_list->bridges[num_bridges].step = 0;
	    num_bridges++;
	}
	fclose(fp_in);
        build_bridge(&file_data_list->bridges[0]);

        //fprintf(stderr, "waiting on ackCondition\n");        

	file_data_list->num_bridges = num_bridges;
	// clean up of writer's files
	if(fp_read_data->fp_comm_rank == 0){
	    unlink(writer_info_filename);
	    unlink(writer_ready_filename);
	}	
    }
    fp->fh = (uint64_t)file_data_list;
    //fp->current_step = file_data_list->bridges[0].step;
    op_msg open_msg;    
    open_msg.process_id = file_data_list->rank;
    open_msg.file_name = strdup(file_data_list->file_name);
    open_msg.type = 1;
    open_msg.step = fp->current_step;
    EVsubmit(file_data_list->bridges[0].op_source, &open_msg, NULL);

    CMCondition_wait(fp_read_data->fp_cm, ackCondition);

    Flush_msg msg;
    msg.rank = fp_read_data->fp_comm_rank;
    msg.type = FORMAT;
    file_data_list->polling = 1;
    // telling writer to flush format.
    EVsubmit(file_data_list->bridges[0].flush_source, &msg, NULL);
    while(file_data_list->polling) {
        CMsleep(fp_read_data->fp_cm, 1);
    }

    int var_count = 0;
    FMStructDescList struct_list = FMcopy_struct_list(format_list_of_FMFormat(file_data_list->current_format));
    FMField *f = struct_list[0].field_list;
    // need to construct var lists here.
    while(f->field_name != NULL)
    {
	//flexpath_var_info * curr_var = (flexpath_var_info*)malloc(sizeof(flexpath_var_info));
	flexpath_var_info * curr_var = new_flexpath_var_info(f->field_name, var_count, f->field_size);
	//curr_var->num_chunks = file_data_list->num_bridges;
	curr_var->num_chunks = 1;
	curr_var->chunks = (flexpath_var_chunk*)malloc(sizeof(flexpath_var_chunk)*curr_var->num_chunks);
	int i;
	memset(curr_var->chunks, 0, sizeof(flexpath_var_chunk)*curr_var->num_chunks);
	for(i = 0; i<curr_var->num_chunks; i++)
	{
	    flexpath_var_chunk * c = &curr_var->chunks[i];
	    c->has_data = 0;
	    c->data = NULL;
	}
	flexpath_var_info * temp = file_data_list->var_list;
	
	curr_var->next = temp;
	file_data_list->var_list = curr_var;
	curr_var->type = adios_integer;
	// because we're only doing scalars here, we know the dims is 0	
        var_count++;
        f++;
    }

    fp->nvars = var_count;
    fp->var_namelist = (char **) malloc(var_count * sizeof(char *));
    f = struct_list[0].field_list;  // f is top-level field list 
    int i=0;
    while(f->field_name != NULL) {
        fp->var_namelist[i++] = strdup(f->field_name);
        //perr("added var name %s\n", f->field_name);
        f++;
    }
    // setting up terminal action for data
    EVassoc_terminal_action(fp_read_data->fp_cm,
			    file_data_list->data_stone,
			    struct_list, data_handler,
			    (void*)file_data_list);

    // telling the writer to flush data
    //TODO: send to split stone to all writers
    msg.type = DATA;
    file_data_list->polling = 1;
    EVsubmit(file_data_list->bridges[0].flush_source, &msg, NULL);

    while(file_data_list->polling) {
        CMsleep(fp_read_data->fp_cm, 1);
    }

    return fp;
}

int adios_read_flexpath_finalize_method ()
{
    return 0;
}

void adios_read_flexpath_release_step(ADIOS_FILE *fp) {
    //perr( "debug:entering release_step\n");

    int i;
    for(i=0; i<file_data_list->num_bridges; i++) {
        if(file_data_list->bridges[i].created) {
            op_msg close;
            close.step = fp->current_step;
            close.type = 0;
            close.process_id = file_data_list->rank;
            close.file_name = fp->path;
            EVsubmit(file_data_list->bridges[i].op_source, &close, NULL);
        }
    }
}
int adios_read_flexpath_advance_step(ADIOS_FILE *fp, int last, float timeout_sec) {
    //fprintf(stderr, "debug:entering advance step\n");
    int i=0;
    for(i=0; i<file_data_list->num_bridges; i++) {
        //fprintf(stderr, "close bridge %d\n", i);
        if(file_data_list->bridges[i].created) {
            op_msg close;
            close.step = fp->current_step;
            close.type = 0;
            close.process_id = file_data_list->rank;
            close.file_name = "test";
            //fprintf(stderr, "submitting close\n");
            EVsubmit(file_data_list->bridges[i].op_source, &close, NULL);
            ///fprintf(stderr, "continuing\n");
        }
        //fprintf(stderr, "reopen bridge %d and wait for ack\n", i);
        if(file_data_list->bridges[i].created) {
            op_msg open;
            open.step = fp->current_step +1;
            open.type = 1;
            open.process_id = file_data_list->rank;
            open.file_name = "test";
            EVsubmit(file_data_list->bridges[i].op_source, &open, NULL);

            //fprintf(stderr, "waiting on ack from bridge %d\n", i);
            CMCondition_wait(fp_read_data->fp_cm, ackCondition);
            //fprintf(stderr, "resuming\n");
        }
    }
    fp->current_step++;
   return fp->current_step;
}

int adios_read_flexpath_close(ADIOS_FILE * fp)
{
    flexpath_file_data * file = (flexpath_file_data*)fp->fh;
    int i;
    op_msg msg;
    msg.type=0;
    msg.file_name = strdup(file->file_name);
    msg.process_id = file->rank;

    //send to each opened link
    for(i = 0; i<file->num_bridges; i++){
        if(file->bridges[i].created) {
            msg.step = file->bridges[i].step;
	    EVsubmit(file->bridges[i].op_source, &msg, NULL);
        }
    }
    // start to cleanup.  Clean up var_lists for now, as the
    // data has already been COPIED over to ADIOS_VARINFO structs
    // that the user maintains a copy of.  Leave stone info in place
    // will handle "invalid" graph setups later
    /* flexpath_var_info * v = file->var_list; */
    /* while(v) */
    /* { */
    /* 	flexpath_var_info * v_tmp = v; */
    /* 	//free(v->varpath); */
    /* 	//free(v->varname); */
    /* 	if(v->chunks == NULL) */
    /* 	{ */
    /* 	    perr( "NULL\n"); */
    /* 	    break; */
    /* 	} */
    /* 	//flexpath_var_chunk * c = &v->chunks[0]; */
    /* 	//flexpath_var_chunk * c_tmp = c; */
    /* 	// free chunks; data has already been copied to user */
    /* 	int i; */
    /* 	for(i = 0; i<v->num_chunks; i++) */
    /* 	{ */
    /* 	    flexpath_var_chunk * c = &v->chunks[i]; */
    /* 	    // have to do this, because for arrays data is actually pre-alloced by user, */
    /* 	    // and not flexpath, so user still has a copy of the data.  for scalars, this isn't so. */
    /* 	    if( (c->data) &&(v->ndims == 0)) */
    /* 	    { */
		
    /* 		//free(c->data); */
    /* 		if(c->global_bounds) */
    /* 		    free(c->global_bounds); */
    /* 		if(c->global_offsets) */
    /* 		    free(c->global_offsets); */
    /* 		if(c->local_bounds) */
    /* 		    free(c->local_bounds); */
    /* 	    } */
    /* 	} */
    /* 	v=v->next; */
    /* 	//free(v_tmp); */
    /* 	if(v == NULL) */
    /* 	    break; */
    /* } */
    return 0;
}

ADIOS_FILE *adios_read_flexpath_fopen(const char *fname, MPI_Comm comm) {
   return 0;
}

int adios_read_flexpath_is_var_timed(const ADIOS_FILE* fp, int varid) { return 0; }

void adios_read_flexpath_get_groupinfo(const ADIOS_FILE *fp, int *ngroups, char ***group_namelist, int **nvars_per_group, int **nattrs_per_group) {}

int adios_read_flexpath_check_reads(const ADIOS_FILE* fp, ADIOS_VARCHUNK** chunk) { perr( "debug:adios function check reads\n"); return 0; }

int adios_read_flexpath_perform_reads(const ADIOS_FILE* fp, int blocking)
{
    //perr( "debug: rank=%d adios function perform reads with blocking: %d\n",
//	  fp_read_data->fp_comm_rank, blocking);
    flexpath_file_data * fd = (flexpath_file_data*)fp->fh;
    Flush_msg msg;
    msg.rank = fp_read_data->fp_comm_rank;
    msg.type = DATA;
    int i;
    int num_sendees = fd->num_sendees;
    for(i = 0; i<num_sendees; i++)
    {
	int sendee = fd->sendees[i];
	/*
        if(file_data_list->bridges[sendee].created) {
            perr("rank %d sending flush to %d\n", fp_read_data->fp_comm_rank, sendee);
        } else {
            perr("rank %d bridge %d not built!!!\n", fp_read_data->fp_comm_rank, sendee);
        }
	*/
	EVsubmit(file_data_list->bridges[sendee].flush_source, &msg, NULL);
    }
    if(blocking)
    {
	fd->polling = 1;
	while(fd->polling)
	{
            CMsleep(fp_read_data->fp_cm, 1);
	}
    }
    return 0;
}
int
adios_read_flexpath_inq_var_blockinfo(const ADIOS_FILE* fp,
				      ADIOS_VARINFO* varinfo)
{ /*perr( "debug:adios function inq var block info\n");*/ return 0; }
int
adios_read_flexpath_inq_var_stat(const ADIOS_FILE* fp,
				 ADIOS_VARINFO* varinfo,
				 int per_step_stat,
				 int per_block_stat)
{ /*perr( "debug:adios function inq var stat\n");*/ return 0; }
void adiosread_flexpath_release_step (ADIOS_FILE *fp);



array_displacements*
get_writer_displacements(int rank, const ADIOS_SELECTION * sel, global_var* gvar){
    int ndims = sel->u.bb.ndim;
    array_displacements * displ = (array_displacements*)malloc(sizeof(array_displacements));
    displ->writer_rank = rank;
    //displ->strides = (stride*)malloc(sizeof(stride)*ndims);
    displ->start = (int*)malloc(sizeof(int) * ndims);
    displ->count = (int*)malloc(sizeof(int) * ndims);    
    displ->ndims = ndims;
    int * offsets = gvar->offsets[0].local_offsets;
    int * local_dims = gvar->offsets[0].local_dimensions;
    int pos = rank * gvar->offsets[0].offsets_per_rank;
    // malloc of ndims size;
    //for each dim
    int i;
    for(i=0; i<ndims; i++){	
	//perr("\t\t\t%d selector start: %llu, selector count: %llu\n", i, sel->u.bb.start[i], sel->u.bb.count[i]);
	//perr("\t\t\t%d offsets[%d]: %d local_dims[%d]: %d\n", i, pos+i, offsets[pos+i], pos+i, local_dims[pos+i]);
	if(sel->u.bb.start[i] >= offsets[pos+i]){
	    int start = sel->u.bb.start[i] - offsets[pos+i];
	    displ->start[i] = start;
	}
	if((sel->u.bb.start[i] + sel->u.bb.count[i] - 1) <= (offsets[pos+i] + local_dims[pos+i] - 1)){	   
	    int count = ((sel->u.bb.start[i] + sel->u.bb.count[i] - 1) - offsets[pos+i]) - displ->start[i] + 1;
	    displ->count[i] = count;
	    
	}else{
	    int count = (local_dims[pos+i] - 1) - displ->start[i] + 1;
	    displ->count[i] = count;
	}
	
	//perr("\t\t\tdispl->start[%d] = %d, count = %d\n", i, displ->start[i], displ->count[i]);
    }
    //print_int_arr("displ->start: ", displ->start, displ->ndims);
    //print_int_arr("displ->count: ", displ->count, displ->ndims);
    return displ;
}

int
need_writer(int j, const ADIOS_SELECTION* sel, evgroup_ptr gp, char* varname) {
    //perr("\n\n\n\n Checking rank %d against a selector\n", j);

    while(file_data_list->gp==NULL) {
        //perr("rank %d waiting for group info\n", file_data_list->rank);
        CMsleep(fp_read_data->fp_cm,1);
    }

    //select var from group
    global_var * gvar = find_gbl_var(gp->vars, varname, gp->num_vars);

    //for each dimension
    int i=0;
    offset_struct var_offsets = gvar->offsets[0];
    for(i=0; i< var_offsets.offsets_per_rank; i++){
        //select sel offsets
        int sel_offset = sel->u.bb.start[i];
        //grab sel dimensions(size)
        int sel_size = sel->u.bb.count[i];
        //perr("sel offset %d with val %d and size %d\n", i, sel_offset, sel_size);


        //select rank offsets
        int rank_offset = var_offsets.local_offsets[j*var_offsets.offsets_per_rank+i];
        //grab rank dimencsions(size)
        int rank_size =var_offsets.local_dimensions[j*var_offsets.offsets_per_rank+i];
        //perr("rank offset %d with val %d and size %d\n", i, rank_offset, rank_size);

        //if rank offset < selector offset and rank offset +size-1 > selector offset
	
        if((rank_offset <= sel_offset) && (rank_offset + rank_size - 1 >=sel_offset)) {
	    // perr("matched overlap type 1\n");
        }
        //if rank offset < selector offset + selector size -1 and rank offset+size-1 > selector offset +selector size -1
        else if((rank_offset <= sel_offset + sel_size - 1) && (rank_offset+rank_size-1>=sel_offset+sel_size-1)) {
            //perr("matched overlap type 2\n");
        } else {
            //perr("overlap not present\n\n");
            return 0;
        }
    }
    //perr("overlap detected\n\n");
    return 1;
}


int adios_read_flexpath_schedule_read_byid(const ADIOS_FILE * fp,
					   const ADIOS_SELECTION * sel,
					   int varid,
					   int from_steps,
					   int nsteps,
					   void * data)
{
    //perr( "debug:schedule_read_byid\n");
    flexpath_file_data * fd = (flexpath_file_data*)(fp->fh);
    flexpath_var_info * v = fd->var_list;
    while(v){
        if(v->id == varid)
        	break;
        else
    	v=v->next;
    }
    if(!v){
        adios_error(err_invalid_varid,
    		"Invalid variable id: %d\n",
    		varid);
        return err_invalid_varid;
    }
    //store the user allocated buffer.
    flexpath_var_chunk * chunk = &v->chunks[0];  
    chunk->data = data;
    v->was_scheduled = 1;
    if(nsteps != 1){
	adios_error (err_invalid_timestep,
                     "Only one step can be read from a stream at a time. "
                     "You requested %d steps in adios_schedule_read()\n", nsteps);
        return err_invalid_timestep;
    }
    v->sel = sel;
    switch(sel->type)
    {
    case ADIOS_SELECTION_WRITEBLOCK:
    {
	int writer_index = sel->u.block.index;
	if(writer_index > fd->num_bridges){
	    adios_error(err_out_of_bound,
			"No process exists on the writer side matching the index.\n");
	}


	//perr( "rank %d sending var message\n", fd->rank);
	//send to what is specified in selector
        int i = 0;
        int found = 0;
        for(i=0; i<fd->num_sendees; i++) {
            if(fd->sendees[i]==writer_index) {
                found=1;
                break;
            }
        }
        if(!found) {
            fd->num_sendees+=1;
            fd->sendees=realloc(fd->sendees, fd->num_sendees*sizeof(int));
            fd->sendees[fd->num_sendees-1] = writer_index;
        }
        if(!fd->bridges[writer_index].created) {
            //perr("rank %d building bridge to %d\n", fp_read_data->fp_comm_rank, writer_index);
            build_bridge(&(fd->bridges[writer_index]));
            op_msg open_msg;
            open_msg.process_id = file_data_list->rank;
            open_msg.file_name = "hey";
            open_msg.type = 1;
            open_msg.step = 0;
            EVsubmit(fd->bridges[writer_index].op_source, &open_msg, NULL);
            CMCondition_wait(fp_read_data->fp_cm, ackCondition);
            //fprintf(stderr, "resuming\n");
	    Var_msg var;
	    var.rank = fd->rank;
            var.var_name = strdup(v->varname);
            /*perr("1 sending %s from %d to %p aka %d\n",
		 var.var_name,
		 var.rank,
		 fd->bridges[writer_index].var_source,
		 writer_index);
	    */
	    EVsubmit(fd->bridges[writer_index].var_source, &var, NULL);
        } else {
	    Var_msg var;
	    var.rank = fd->rank;
            var.var_name = strdup(v->varname);
	    /*
            perr("2 sending %s from %d to %p aka %d\n",
		 var.var_name,
		 var.rank,
		 fd->bridges[writer_index].var_source,
		 writer_index);
	    */
	    EVsubmit(fd->bridges[writer_index].var_source, &var, NULL);
        }
    //perr("rank %d sent var msg to %d\n", fp_read_data->fp_comm_rank, writer_index);
	break;
    }
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        int j=0;
	int need_count = 0;
	array_displacements * all_disp = NULL;	
        for(j=0; j<fd->size; j++) {
            int reader=0;	    	    
            if(need_writer(j, sel, fd->gp, v->varname)==1){		
		need_count++;
                reader = j;
		global_var * gvar = find_gbl_var(fd->gp->vars, v->varname, fd->gp->num_vars);
		array_displacements * displ = get_writer_displacements(j, sel, gvar);
		all_disp = realloc(all_disp, sizeof(array_displacements)*need_count);
		all_disp[need_count-1] = *displ;
            } else {
                continue;
            }
            int i = 0;
            int found = 0;
            for(i=0; i<fd->num_sendees; i++) {
                if(fd->sendees[i]==reader) {
                    found=1;
                    break;
                }
            }
            if(!found) {
                fd->num_sendees+=1;
                fd->sendees=realloc(fd->sendees, fd->num_sendees*sizeof(int));
                fd->sendees[fd->num_sendees-1] = reader;
            }
            if(!fd->bridges[reader].created) {
                //perr("rank %d building bridge to %d\n", 
		//    fp_read_data->fp_comm_rank, reader);
                build_bridge(&(fd->bridges[reader]));
                op_msg open_msg;
                open_msg.process_id = file_data_list->rank;
                open_msg.file_name = "hey";
                open_msg.type = 1;
                EVsubmit(fd->bridges[reader].op_source, &open_msg, NULL);
	        Var_msg var;
	        var.rank = fd->rank;
                var.var_name = strdup(v->varname);
                //perr("1 sending %s from %d to %p aka %d\n", 
		//   var.var_name, var.rank, fd->bridges[reader].var_source, reader);
	        EVsubmit(fd->bridges[reader].var_source, &var, NULL);
            } else {
	        Var_msg var;
	        var.rank = fd->rank;
                var.var_name = strdup(v->varname);
                //perr("2 sending %s from %d to %p aka %d\n", 
		//   var.var_name, var.rank, fd->bridges[reader].var_source, reader);
	        EVsubmit(fd->bridges[reader].var_source, &var, NULL);
            }
            //perr("rank %d sent var msg to %d\n", 
	    // fp_read_data->fp_comm_rank,reader);
	}
	v->displ = all_disp;
	v->num_displ = need_count;
        break;
    }
    case ADIOS_SELECTION_AUTO:
    {
	adios_error(err_operation_not_supported,
		    "ADIOS_SELECTION_AUTO not yet supported by flexpath.");
	break;
    }
    case ADIOS_SELECTION_POINTS:
    {
	adios_error(err_operation_not_supported,
		    "ADIOS_SELECTION_POINTS not yet supported by flexpath.");
	break;
    }
    }
    return 0;
}

int adios_read_flexpath_schedule_read(const ADIOS_FILE *fp,
			const ADIOS_SELECTION * sel,
			const char * varname,
			int from_steps,
			int nsteps,
			void * data)
{
    //perr( "in schedule_read\n");
    return 0;
}

int adios_read_flexpath_fclose(ADIOS_FILE *fp)
{
    //perr( "debug: adios_read_flexpath_fclose\n");
    return 0;
}

int * adios_read_flexpath_gopen (ADIOS_FILE *fp, const char *grpname)
{
    //perr( "debug: adios_read_flexpath_gopen\n");
    return NULL;
}

int * adios_read_flexpath_gopen_byid (ADIOS_FILE *fp, int grpid)
{
    //perr( "debug: adios_read_flexpath_gopen_byid\n");
    return NULL;
}

int adios_read_flexpath_gclose (int *gp)
{
    perr( "debug: adios_read_flexpath_gclose\n");
//    flexpath_read_file_data *ds = (flexpath_read_file_data *) gp->fp->fh;
    adios_errno = 0;
    int i;
//    free_namelist ((gp->var_namelist),gp->vars_count);
    free(gp);
    //perr( "im here %s:%d\n",__FILE__,__LINE__);
    return 0;

}

int adios_read_flexpath_get_attr (int *gp, const char *attrname,
                                 enum ADIOS_DATATYPES *type,
                                 int *size, void **data)
{
    //perr( "debug: adios_read_flexpath_get_attr\n");
    // TODO: borrowed from dimes
    adios_error(err_invalid_read_method, "adios_read_flexpath_get_attr is not implemented.");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}

int adios_read_flexpath_get_attr_byid (const ADIOS_FILE *fp, int attrid,
                                      enum ADIOS_DATATYPES *type,
                                      int *size, void **data)
{
//    perr( "debug: adios_read_flexpath_get_attr_byid\n");
    // TODO: borrowed from dimes
    adios_error(err_invalid_read_method, "adios_read_flexpath_get_attr_byid is not implemented.");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}

ADIOS_VARINFO* adios_read_flexpath_inq_var(const ADIOS_FILE * fp, const char* varname)
{
    //perr( "debug: adios_read_flexpath_inq_var\n");

    ADIOS_VARINFO* v = malloc(sizeof(ADIOS_VARINFO));
    if(!v) {
        adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
        return NULL;
    }
    memset(v, 0, sizeof(ADIOS_VARINFO));

    flexpath_file_data *ds = (flexpath_file_data *) fp->fh;
    flexpath_var_info *current_var = find_fp_var(ds->var_list, varname);
    if(current_var) {
	v = convert_file_info(current_var, v, varname, fp);
	return v;
    }
    else {
        adios_error(err_invalid_varname, "Cannot find var %s\n", varname);
        return NULL;
    }
}

ADIOS_VARINFO * adios_read_flexpath_inq_var_byid (const ADIOS_FILE * fp, int varid)
{
    //perr( "debug: inq_var_byid\n");
    if(varid >= 0 && varid < fp->nvars) {
        return adios_read_flexpath_inq_var(fp, fp->var_namelist[varid]);
    }
    else {
        adios_error(err_invalid_varid, "Cannot find var %d\n", varid);
        return NULL;
    }
}

ADIOS_VARINFO* 
convert_file_info(flexpath_var_info * current_var,
				 ADIOS_VARINFO * v, const char* varname,
				 const ADIOS_FILE * fp)
{
    int i;
    current_var->type = v->type;
    for(i = 0; i < fp->nvars; i ++) {
	if(!strcmp(fp->var_namelist[i], varname)) {
	    v->varid = i; // TODO: this may not be cmpatible with BP
	    break;
	}
    }
    v->type = current_var->type;
    v->ndim = current_var->ndims;
    //v->timedim = current_var->time_dim;
    if(v->ndim == 0){    
	//int value_size = common_read_type_size(v->type, current_var->chunks->data);
	int value_size = current_var->data_size;
	v->value = malloc(value_size);
	if(!v->value) {
	    adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
	    return NULL;
	}
	flexpath_var_chunk * chunk = &current_var->chunks[0];
	memcpy(v->value, chunk->data, value_size);
    }else{ // arrays
	perr( "ARRAYS\n");
	v->dims = (uint64_t *) malloc(v->ndim * sizeof(uint64_t));
	if(!v->dims) {
	    adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
	    return NULL;
	}
	int k;
	for(k = 0; k < v->ndim; k ++) {
	    //v->dims[k] = ds->pgs[i].vars[j].global_bounds[k];
	    v->dims[k] = current_var->chunks->global_bounds[k];
	}
    }
    return v;
}


void adios_read_flexpath_free_varinfo (ADIOS_VARINFO *vp)
{
    //perr( "debug: adios_read_flexpath_free_varinfo\n");
    return;
}

int64_t adios_read_flexpath_read_var (int *gp, const char *varname,
                                     const uint64_t *start, const uint64_t *count,
                                     void *data)
{
    //perr( "debug: adios_read_flexpath_read_var\n");
    return (int64_t)0;
}

int64_t adios_read_flexpath_read_var_byid (int *gp, int varid,
                                          const uint64_t *start,
                                          const uint64_t *count,
                                          void *data)
{
    //perr( "debug: adios_read_flexpath_read_var_byid\n");
    return (int64_t)0;
}

void adios_read_flexpath_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{
    //perr( "debug: adios_read_flexpath_reset_dimension_order\n");
    adios_error(err_invalid_read_method, "adios_read_flexpath_reset_dimension_order is not implemented.");
}
