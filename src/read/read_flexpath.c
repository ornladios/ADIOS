/*
    read_flexpath.c
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
//#include <gen_thread.h>
#include <evpath.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

// local libraries
#include "config.h"
#include "public/adios.h"
#include "public/adios_types.h"
#include "public/adios_read_v2.h"
#include "core/adios_read_hooks.h"
#include "core/adios_logger.h"
#include "core/common_read.h"
#include "public/adios_error.h"
#include "core/flexpath.h"

#include "core/transforms/adios_transforms_common.h" // NCSU ALACRITY-ADIOS

// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define FP_BATCH_SIZE 4

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef struct _bridge_info
{
    EVstone bridge_stone;
    EVsource flush_source;
    EVsource var_source;
    EVsource op_source;
    int their_num;
    char *contact;
    int created;
    int opened;
    int step;
    int scheduled;
}bridge_info;

typedef struct _flexpath_read_request
{
    long pad1, pad2, pad3, pad4, pad5;
    int num_pending;
    int num_completed;
    int condition;
}flexpath_read_request;

typedef struct _flexpath_var_chunk
{
    int has_data;
    int rank; // the writer's rank this chunk represents. not used or needed right now
    void *data;
    void *user_buf;
} flexpath_var_chunk;

/*
 * Contains start & counts for each dimension for a writer_rank.
 */
typedef struct _array_displ
{
    int writer_rank;
    int ndims;
    uint64_t pos;
    uint64_t *start;
    uint64_t *count;
}array_displacements;

typedef struct _flexpath_var
{
    int id;
    char *varname;
    char *varpath;

    enum ADIOS_DATATYPES type;
    uint64_t type_size; // type size, not arrays size

    int time_dim; // -1 means it is not a time dimension
    int ndims;
    uint64_t *global_dims; // ndims size (if ndims>0)
    uint64_t *local_dims; // for local arrays
    uint64_t array_size; // not relevant for scalars

    int num_chunks;
    flexpath_var_chunk *chunks;

    int num_displ;
    array_displacements *displ;

    ADIOS_SELECTION *sel;
    uint64_t start_position;

    struct _flexpath_var *next;
} flexpath_var;

typedef struct _flexpath_reader_file
{
    char *file_name;
    char *group_name; // assuming one group per file right now.
    int host_language;

    EVstone stone;

    MPI_Comm comm;
    int rank;
    int size;

    int num_bridges;
    bridge_info *bridges;
    int writer_coordinator;

    int num_vars;
    flexpath_var *var_list;
    evgroup *gp;

    int writer_finalized;
    int last_writer_step;
    int mystep;
    int num_sendees;
    int *sendees;

    uint64_t data_read; // for perf measurements.
    double time_in; // for perf measurements.
    flexpath_read_request req;
    pthread_mutex_t data_mutex;
    pthread_cond_t data_condition;
} flexpath_reader_file;

typedef struct _local_read_data
{
    // MPI stuff
    MPI_Comm comm;
    int rank;
    int size;

    // EVPath stuff
    CManager cm;
    atom_t CM_TRANSPORT;
} flexpath_read_data;

flexpath_read_data* fp_read_data = NULL;

/********** Helper functions. **********/

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

void build_bridge(bridge_info* bridge)
{
    attr_list contact_list = attr_list_from_string(bridge->contact);
    if (bridge->created == 0) {
	bridge->bridge_stone =
	    EVcreate_bridge_action(fp_read_data->cm,
				   contact_list,
				   (EVstone)bridge->their_num);

	bridge->flush_source =
	    EVcreate_submit_handle(fp_read_data->cm,
				   bridge->bridge_stone,
				   flush_format_list);

	bridge->var_source =
	    EVcreate_submit_handle(fp_read_data->cm,
				   bridge->bridge_stone,
				   var_format_list);

	bridge->op_source =
	    EVcreate_submit_handle(fp_read_data->cm,
				   bridge->bridge_stone,
				   op_format_list);

	bridge->created = 1;
    }
}

void
free_displacements(array_displacements *displ, int num)
{
    if (displ) {
	int i;
	for (i=0; i<num; i++) {
	    free(displ[i].start);
	    free(displ[i].count);
	}
	free(displ);
    }
}

void
free_evgroup(evgroup *gp)
{
    EVreturn_event_buffer(fp_read_data->cm, gp);
}

flexpath_var*
new_flexpath_var(const char *varname, int id, uint64_t type_size)
{
    flexpath_var *var = malloc(sizeof(flexpath_var));
    if (var == NULL) {
	log_error("Error creating new var: %s\n", varname);
	return NULL;
    }

    memset(var, 0, sizeof(flexpath_var));
    // free this when freeing vars.
    var->varname = strdup(varname);
    var->time_dim = -1;
    var->id = id;
    var->type_size = type_size;
    var->displ = NULL;
    return var;
}

flexpath_reader_file*
new_flexpath_reader_file(const char *fname)
{
    flexpath_reader_file *fp = calloc(1, sizeof(*fp));
    if (fp == NULL) {
	log_error("Cannot create data for new file.\n");
	exit(1);
    }
    fp->file_name = strdup(fname);
    fp->writer_coordinator = -1;
    fp->last_writer_step = -1;
    fp->req = (flexpath_read_request) {.pad1 = 0, .pad2 = 0, .pad3 = 0, .pad4 = 0, .pad5=0,.num_pending = 0, .num_completed = 0, .condition = -1};
    pthread_mutex_init(&fp->data_mutex, NULL);
    pthread_cond_init(&fp->data_condition, NULL);
    return fp;
}

enum ADIOS_DATATYPES
ffs_type_to_adios_type(const char *ffs_type, int size)

{
    char *bracket = "[";
    size_t posfound = strcspn(ffs_type, bracket);
    char *filtered_type = NULL;
    if (strlen(ffs_type) == strlen(bracket)) {
        filtered_type = strdup(ffs_type);
    }
    else {
        filtered_type = malloc(posfound+1);
        memset(filtered_type, '\0', posfound+1);
        filtered_type = strncpy(filtered_type, ffs_type, posfound);
    }

    if (!strcmp("integer", filtered_type)) {
	if (size == sizeof(int)) {
	    return adios_integer;
	} else if (size == sizeof(long)) {
	    return adios_long;
	}
    }
    else if (!strcmp("unsigned integer", filtered_type)) {
	    if (size == sizeof(unsigned int)) {
		return adios_unsigned_integer;
	    } else if (size == sizeof(unsigned long)) {
		return adios_unsigned_long;
	    }
	}	   
    else if (!strcmp("float", filtered_type))
	return adios_real;
    else if (!strcmp("string", filtered_type))
	return adios_string;
    else if (!strcmp("double", filtered_type)) {
	if (size == sizeof(double)) {
	    return adios_double;
	} else if (size == sizeof(long double)) {
	    return adios_long_double;
	}
    }
    else if (!strcmp("char", filtered_type))
	return adios_byte;
    else if (!strcmp("complex", filtered_type))
	return adios_complex;
    else if (!strcmp("double_complex", filtered_type))
        return adios_double_complex;
    else {
	fprintf(stderr, "returning unknown for: ffs_type: %s\n", ffs_type);
	return adios_unknown;
    }
}

ADIOS_VARINFO*
convert_var_info(flexpath_var * fpvar,
		 ADIOS_VARINFO * v,
		 const char* varname,
		 const ADIOS_FILE *adiosfile)
{
    int i;
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    v->type = fpvar->type;
    v->ndim = fpvar->ndims;
    // needs to change. Has to get information from write.
    v->nsteps = 1;
    v->nblocks = malloc(sizeof(int)*v->nsteps);
    v->sum_nblocks = 1;
    v->nblocks[0] = 1;
    v->statistics = NULL;
    v->blockinfo = NULL;

    if (v->ndim == 0) {
	int value_size = fpvar->type_size;
	v->value = malloc(value_size);
	if(!v->value) {
	    adios_error(err_no_memory,
			"Cannot allocate buffer in adios_read_flexpath_inq_var()");
	    return NULL;
	}
	flexpath_var_chunk * chunk = &fpvar->chunks[0];
	memcpy(v->value, chunk->data, value_size);
	v->global = 0;
    } else { // arrays
	v->dims = (uint64_t*)malloc(v->ndim * sizeof(uint64_t));
	if (!v->dims) {
	    adios_error(err_no_memory,
			"Cannot allocate buffer in adios_read_flexpath_inq_var()");
	    return NULL;
	}
	// broken.  fix. -- why did I put this comment here?
	int cpysize = fpvar->ndims*sizeof(uint64_t);
	if (fpvar->global_dims) {
	    v->global = 1;
	    memcpy(v->dims, fpvar->global_dims, cpysize);
	}
	else {
	    v->global = 0;
	}
    }
    return v;
}

flexpath_var *
find_fp_var(flexpath_var * var_list, const char * varname)
{
    while (var_list) {
	if (!strcmp(varname, var_list->varname)) {
	    return var_list;
	}
	var_list = var_list->next;
    }
    return NULL;
}

global_var*
find_gbl_var(global_var *vars, const char *name, int num_vars)
{
    int i;
    for (i=0; i<num_vars; i++) {
	if (!strcmp(vars[i].name, name))
	    return &vars[i];
    }
    return NULL;
}

static FMField*
find_field_by_name (const char *name, const FMFieldList flist)
{
    FMField *f = flist;
    while (f->field_name != NULL)
    {
        if (!strcmp(name, f->field_name))
            return f;
        else
            f++;
    }
    return NULL;
}

static uint64_t
calc_ffspacket_size(FMField *f, attr_list attrs, void *buffer)
{
    uint64_t size = 0;
    while (f->field_name) {
        char atom_name[200] = "";

        strcat(atom_name, FP_NDIMS_ATTR_NAME);
        strcat(atom_name, "_");
        strcat(atom_name, f->field_name);
        int num_dims = 0;
        get_int_attr(attrs, attr_atom_from_string(atom_name), &num_dims);
	if (num_dims == 0) {
	    size += (uint64_t)f->field_size;
	}
	else {
	    int i;
	    uint64_t tmpsize = (uint64_t)f->field_size;
	    for (i=0; i<num_dims; i++) {
		char *dim;
		char atom_name[200] ="";
		strcat(atom_name, FP_DIM_ATTR_NAME);
		strcat(atom_name, "_");
		strcat(atom_name, f->field_name);
		strcat(atom_name, "_");
		char dim_num[10] = "";
		sprintf(dim_num, "%d", i+1);
		strcat(atom_name, dim_num);
		get_string_attr(attrs, attr_atom_from_string(atom_name), &dim);

		FMField *temp_field = find_field_by_name(dim, f);
		uint64_t *temp_val = get_FMfieldAddr_by_name(temp_field,
							     temp_field->field_name,
							     buffer);
		uint64_t dimval = *temp_val;
		tmpsize *= dimval;
	    }
	    size += tmpsize;
	}
	f++;
    }
    return size;
}

void
print_displacement(array_displacements *disp, int myrank)
{
    printf("rank: %d, writer displacements for writer: %d\n", myrank, disp->writer_rank);
    printf("\tndims: %d, pos: %d\n", disp->ndims, (int)disp->pos);

    int i;
    printf("\tstarts: ");
    for (i = 0; i<disp->ndims; i++) {
	printf("%d ", (int)disp->start[i]);
    }
    printf("\n");
    printf("\tcounts: ");
    for (i = 0; i<disp->ndims; i++) {
	printf("%d ", (int)disp->count[i]);
    }
    printf("\n");
}

/*
 * Finds the array displacements for a writer identified by its rank.
 */
array_displacements*
find_displacement(array_displacements *list, int rank, int num_displ)
{
    int i;
    for (i=0; i<num_displ; i++) {
	if (list[i].writer_rank == rank)
	    return &list[i];
    }
    return NULL;
}


array_displacements*
get_writer_displacements(
    int writer_rank,
    const ADIOS_SELECTION * sel,
    global_var* gvar,
    uint64_t *size)
{
    int ndims = sel->u.bb.ndim;
    //where do I free these?
    array_displacements *displ = malloc(sizeof(array_displacements));
    displ->writer_rank = writer_rank;

    displ->start = malloc(sizeof(uint64_t) * ndims);
    displ->count = malloc(sizeof(uint64_t) * ndims);
    memset(displ->start, 0, sizeof(uint64_t) * ndims);
    memset(displ->count, 0, sizeof(uint64_t) * ndims);

    displ->ndims = ndims;
    uint64_t *offsets = gvar->offsets[0].local_offsets;
    uint64_t *local_dims = gvar->offsets[0].local_dimensions;
    uint64_t pos = writer_rank * gvar->offsets[0].offsets_per_rank;

    int i;
    int _size = 1;
    for (i=0; i<ndims; i++) {
	if (sel->u.bb.start[i] >= offsets[pos+i]) {
	    int start = sel->u.bb.start[i] - offsets[pos+i];
	    displ->start[i] = start;
	}
	if ((sel->u.bb.start[i] + sel->u.bb.count[i] - 1) <=
	   (offsets[pos+i] + local_dims[pos+i] - 1)) {
	    int count = ((sel->u.bb.start[i] + sel->u.bb.count[i] - 1) -
			 offsets[pos+i]) - displ->start[i] + 1;
	    displ->count[i] = count;


	} else {
	    int count = (local_dims[pos+i] - 1) - displ->start[i] + 1;
	    displ->count[i] = count;
	}
	_size *= displ->count[i];
    }
    *size = _size;
    return displ;
}

int
need_writer(
    flexpath_reader_file *fp,
    int writer,
    const ADIOS_SELECTION *sel,
    evgroup_ptr gp,
    char *varname)
{
    //select var from group
    global_var *gvar = find_gbl_var(gp->vars, varname, gp->num_vars);

    //for each dimension
    int i=0;
    offset_struct var_offsets = gvar->offsets[0];
    for (i=0; i< var_offsets.offsets_per_rank; i++) {
	int pos = writer*(var_offsets.offsets_per_rank) + i;

        uint64_t sel_offset = sel->u.bb.start[i];
        uint64_t sel_size = sel->u.bb.count[i];

        uint64_t rank_offset = var_offsets.local_offsets[pos];
        uint64_t rank_size = var_offsets.local_dimensions[pos];

	/* printf("need writer rank: %d writer: %d sel_start: %d sel_count: %d rank_offset: %d rank_size: %d\n", */
	/*        fp->rank, writer, (int)sel_offset, (int)sel_size, (int)rank_offset, (int)rank_size); */

        if ((rank_offset <= sel_offset) && (rank_offset + rank_size - 1 >=sel_offset)) {
	     log_debug("matched overlap type 1\n");
        }

        else if ((rank_offset <= sel_offset + sel_size - 1) && \
		(rank_offset+rank_size>=sel_offset+sel_size-1)) {
        } else if ((sel_offset <= rank_offset) && (rank_offset+rank_size<= sel_offset+sel_size-1)) {
        } else {
            return 0;
        }
    }
    return 1;
}

void
free_fmstructdesclist(FMStructDescList struct_list)
{
    FMField *f = struct_list[0].field_list;

    //cant free field_name because it's const.
    /* FMField *temp = f; */
    /* while (temp->field_name) { */
    /* 	free(temp->field_name); */
    /* 	temp++; */
    /* } */
    free(f);
    free(struct_list[0].opt_info);
    free(struct_list);
}

int
get_ndims_attr(const char *field_name, attr_list attrs)
{
    char atom_name[200] = "";
    strcat(atom_name, FP_NDIMS_ATTR_NAME);
    strcat(atom_name, "_");
    strcat(atom_name, field_name);
    int num_dims;
    atom_t atm = attr_atom_from_string(atom_name);
    get_int_attr(attrs, atm, &num_dims);
    return num_dims;
}

flexpath_var*
setup_flexpath_vars(FMField *f, int *num)
{
    flexpath_var *vars = NULL;
    int var_count = 0;

    while (f->field_name != NULL) {
	flexpath_var *curr_var = new_flexpath_var(f->field_name,
						  var_count,
						  f->field_size);
	curr_var->num_chunks = 1;
	curr_var->chunks =  malloc(sizeof(flexpath_var_chunk)*curr_var->num_chunks);
	memset(curr_var->chunks, 0, sizeof(flexpath_var_chunk)*curr_var->num_chunks);
	curr_var->sel = NULL;
	curr_var->type = ffs_type_to_adios_type(f->field_type, f->field_size);
	flexpath_var *temp = vars;
	curr_var->next = temp;
	vars = curr_var;
	if (strncmp(f->field_name, "FPDIM", 5)) {
	    var_count++;
	}
	f++;
    }
    *num = var_count;
    return vars;
}

/*****************Messages to writer procs**********************/

void
send_open_msg(flexpath_reader_file *fp, int destination)
{
    if (!fp->bridges[destination].created) {
	build_bridge(&(fp->bridges[destination]));
    }
    op_msg msg;
    msg.process_id = fp->rank;
    msg.file_name = fp->file_name;
    msg.step = fp->mystep;
    msg.type = OPEN_MSG;
    int cond = CMCondition_get(fp_read_data->cm, NULL);
    msg.condition = cond;

    EVsubmit(fp->bridges[destination].op_source, &msg, NULL);
    CMCondition_wait(fp_read_data->cm, cond);
    fp->bridges[destination].opened = 1;
}

void
send_close_msg(flexpath_reader_file *fp, int destination)
{
    if (!fp->bridges[destination].created) {
	build_bridge(&(fp->bridges[destination]));
    }
    op_msg msg;
    msg.process_id = fp->rank;
    msg.file_name = fp->file_name;
    msg.step = fp->mystep;
    msg.type = CLOSE_MSG;
    //msg.condition = -1;
    int cond = CMCondition_get(fp_read_data->cm, NULL);
    msg.condition = cond;
    EVsubmit(fp->bridges[destination].op_source, &msg, NULL);
    CMCondition_wait(fp_read_data->cm, cond);
    fp->bridges[destination].opened = 0;
}

void
send_flush_msg(flexpath_reader_file *fp, int destination, Flush_type type, int use_condition)
{
    Flush_msg msg;
    msg.type = type;
    msg.process_id = fp->rank;
    msg.id = fp->mystep;

    if (use_condition)
	msg.condition = CMCondition_get(fp_read_data->cm, NULL);
    else
	msg.condition = -1;
    // maybe check to see if the bridge is create first.
    EVsubmit(fp->bridges[destination].flush_source, &msg, NULL);
    if (use_condition) {
	CMCondition_wait(fp_read_data->cm, msg.condition);
    }
}

void
send_var_message(flexpath_reader_file *fp, int destination, char *varname)
{
        int i = 0;
        int found = 0;
        for (i=0; i<fp->num_sendees; i++) {
            if (fp->sendees[i]==destination) {
                found=1;
                break;
            }
        }
        if (!found) {
            fp->num_sendees+=1;
            fp->sendees=realloc(fp->sendees, fp->num_sendees*sizeof(int));
            fp->sendees[fp->num_sendees-1] = destination;
        }
        if (!fp->bridges[destination].created) {
            build_bridge(&(fp->bridges[destination]));
	}
	if (!fp->bridges[destination].opened) {
	    fp->bridges[destination].opened = 1;
	    send_open_msg(fp, destination);
	}
	Var_msg var;
	var.process_id = fp->rank;
	var.var_name = varname;
	EVsubmit(fp->bridges[destination].var_source, &var, NULL);
}

/********** EVPath Handlers **********/

static int
update_step_msg_handler(
    CManager cm,
    void *vevent,
    void *client_data,
    attr_list attrs)
{
    update_step_msg *msg = (update_step_msg*)vevent;
    ADIOS_FILE *adiosfile = (ADIOS_FILE*)client_data;
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;

    fp->last_writer_step = msg->step;
    fp->writer_finalized = msg->finalized;
    adiosfile->last_step = msg->step;
    CMCondition_signal(fp_read_data->cm, msg->condition);
    return 0;
}

static int
op_msg_handler(CManager cm, void *vevent, void *client_data, attr_list attrs) {
    op_msg* msg = (op_msg*)vevent;
    ADIOS_FILE *adiosfile = (ADIOS_FILE*)client_data;
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    if (msg->type==ACK_MSG) {
	if (msg->condition != -1) {
	    CMCondition_signal(fp_read_data->cm, msg->condition);
	}
        //ackCondition = CMCondition_get(fp_read_data->cm, NULL);
    }
    if (msg->type == EOS_MSG) {
	adios_errno = err_end_of_stream;
	CMCondition_signal(fp_read_data->cm, msg->condition);
    }
    return 0;
}

static int
group_msg_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    EVtake_event_buffer(fp_read_data->cm, vevent);
    evgroup *msg = (evgroup*)vevent;
    ADIOS_FILE *adiosfile = client_data;
    flexpath_reader_file * fp = (flexpath_reader_file*)adiosfile->fh;
    if (fp->group_name == NULL) {
        fp->group_name = strdup(msg->group_name);
    }
    fp->gp = msg;
    int i;
    for (i = 0; i<msg->num_vars; i++) {
	global_var *gblvar = &msg->vars[i];
	flexpath_var *fpvar = find_fp_var(fp->var_list, gblvar->name);
	if (fpvar) {
	    offset_struct *offset = &gblvar->offsets[0];
	    uint64_t *local_dimensions = offset->local_dimensions;
	    uint64_t *local_offsets = offset->local_offsets;
	    uint64_t *global_dimensions = offset->global_dimensions;

	    fpvar->ndims = offset->offsets_per_rank;
	    fpvar->global_dims = malloc(sizeof(uint64_t)*fpvar->ndims);
	    memcpy(fpvar->global_dims, global_dimensions, sizeof(uint64_t)*fpvar->ndims);
	} else {
	    adios_error(err_corrupted_variable,
			"Mismatch between global variables and variables specified %s.",
			gblvar->name);
	    return err_corrupted_variable;
	}
    }
    CMCondition_signal(fp_read_data->cm, msg->condition);
    return 0;
}


int
increment_index(int64_t ndim, uint64_t *dimen_array, uint64_t *index_array)
{
    ndim--;
    while (ndim >= 0) {
	index_array[ndim]++;
	if (index_array[ndim] < dimen_array[ndim]) {
	    return 1;
	}
	index_array[ndim] = 0;
	ndim--;
    }
    return 0;
}

void
map_local_to_global_index(uint64_t ndim, uint64_t *local_index, uint64_t *local_offsets, uint64_t *global_index)
{
    int i;
    for (i=0; i < ndim; i++) {
	global_index[i] = local_index[i] + local_offsets[i];
    }
}

void
map_global_to_local_index(uint64_t ndim, uint64_t *global_index, uint64_t *local_offsets, uint64_t *local_index)
{
    int i;
    for (i=0; i < ndim; i++) {
	local_index[i] = global_index[i] - local_offsets[i];
    }
}

int
find_offset(uint64_t ndim, uint64_t *size, uint64_t *index)
{
    int offset = 0;
    int i;
    for (i=0; i< ndim; i++) {
	offset = index[i] + (size[i] * offset);
    }
    return offset;
}

/*
 *  element_size is the byte size of the array elements
 *  dims is the number of dimensions in the variable
 *  global_dimens is an array, dims long, giving the size of each dimension
 *  partial_offsets is an array, dims long, giving the starting offsets per dimension
 *      of this data block in the global array
 *  partial_counts is an array, dims long, giving the size per dimension
 *      of this data block in the global array
 *  selection_offsets is an array, dims long, giving the starting offsets in the global array
 *      of the output selection.
 *  selection_counts is an array, dims long, giving the size per dimension
 *      of the output selection.
 *  data is the input, a slab of the global array
 *  selection is the output, to be filled with the selection array.
 */
void
extract_selection_from_partial(int element_size, uint64_t dims, uint64_t *global_dimens,
			       uint64_t *partial_offsets, uint64_t *partial_counts,
			       uint64_t *selection_offsets, uint64_t *selection_counts,
			       char *data, char *selection)
{
    int block_size;
    int source_block_stride;
    int dest_block_stride;
    int source_block_start_offset;
    int dest_block_start_offset;
    int block_count;
    int dim;
    int operant_dims;
    int operant_element_size;

    block_size = 1;
    operant_dims = dims;
    operant_element_size = element_size;
    for (dim = dims-1; dim >= 0; dim--) {
	if ((global_dimens[dim] == partial_counts[dim]) &&
	    (selection_counts[dim] == partial_counts[dim])) {
	    block_size *= global_dimens[dim];
	    operant_dims--;   /* last dimension doesn't matter, we got all and we want all */
	    operant_element_size *= global_dimens[dim];
	} else {
	    int left = MAX(partial_offsets[dim], selection_offsets[dim]);
	    int right = MIN(partial_offsets[dim] + partial_counts[dim],
			    selection_offsets[dim] + selection_counts[dim]);
	    block_size *= (right - left);
	    break;
	}
    }
    source_block_stride = partial_counts[operant_dims-1] * operant_element_size;
    dest_block_stride = selection_counts[operant_dims-1] * operant_element_size;

    /* calculate first selected element and count */
    block_count = 1;
    uint64_t *first_index = malloc(dims * sizeof(first_index[0]));
    for (dim = 0; dim < dims; dim++) {
	int left = MAX(partial_offsets[dim], selection_offsets[dim]);
	int right = MIN(partial_offsets[dim] + partial_counts[dim],
			selection_offsets[dim] + selection_counts[dim]);
	if (dim < operant_dims-1) {
	    block_count *= (right-left);
	}
	first_index[dim] = left;
    }
    uint64_t *selection_index = malloc(dims * sizeof(selection_index[0]));
    map_global_to_local_index(dims, first_index, selection_offsets, selection_index);
    dest_block_start_offset = find_offset(dims, selection_counts, selection_index);
    free(selection_index);
    dest_block_start_offset *= element_size;

    uint64_t *partial_index = malloc(dims * sizeof(selection_index[0]));
    map_global_to_local_index(dims, first_index, partial_offsets, partial_index);
    source_block_start_offset = find_offset(dims, partial_counts, partial_index);
    free(partial_index);
    source_block_start_offset *= element_size;

    data += source_block_start_offset;
    selection += dest_block_start_offset;
    int i;
    for (i=0 ; i < block_count; i++) {
	memcpy(selection, data, block_size * element_size);
	data += source_block_stride;
	selection += dest_block_stride;
    }
}

    static int
raw_handler(CManager cm, void *vevent, int len, void *client_data, attr_list attrs)
{
    ADIOS_FILE *adiosfile = client_data;
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    int writer_rank;
    int flush_id;
    double data_start;
    get_double_attr(attrs, attr_atom_from_string("fp_starttime"), &data_start);
    get_int_attr(attrs, attr_atom_from_string(FP_RANK_ATTR_NAME), &writer_rank);
    get_int_attr(attrs, attr_atom_from_string("fp_flush_id"), &flush_id);
    /* fprintf(stderr, "\treader rank:%d:got data from writer:%d:step:%d\n", */
    /* 	    fp->rank, writer_rank, fp->mystep); */
    FMContext context = CMget_FMcontext(cm);
    void *base_data = FMheader_skip(context, vevent);
    FMFormat format = FMformat_from_ID(context, vevent);

    // copy //FMfree_struct_desc_list call
    FMStructDescList struct_list =
	FMcopy_struct_list(format_list_of_FMFormat(format));
    FMField *f = struct_list[0].field_list;

#if 0
    uint64_t packet_size = calc_ffspacket_size(f, attrs, base_data);
    fp->data_read += packet_size;
#endif
    /* setting up initial vars from the format list that comes along with the
       message. Message contains both an FFS description and the data. */
    if (fp->num_vars == 0) {
	int var_count = 0;
	fp->var_list = setup_flexpath_vars(f, &var_count);

	adiosfile->var_namelist = malloc(var_count * sizeof(char *));
	int i = 0;
	while (f->field_name != NULL) {
	    if (strncmp(f->field_name, "FPDIM", 5)) {
		adiosfile->var_namelist[i++] = strdup(f->field_name);
	    }
	    f++;
	}
	adiosfile->nvars = var_count;
	fp->num_vars = var_count;
    }

    f = struct_list[0].field_list;
    char *curr_offset = NULL;

    while (f->field_name) {
        char atom_name[200] = "";
    	flexpath_var *var = find_fp_var(fp->var_list, f->field_name);

    	if (!var) {
    	    adios_error(err_file_open_error,
    			"file not opened correctly.  var does not match format.\n");
    	    return err_file_open_error;
    	}

	int num_dims = get_ndims_attr(f->field_name, attrs);
    	var->ndims = num_dims;

	flexpath_var_chunk *curr_chunk = &var->chunks[0];

	// has the var been scheduled?
	if (var->sel) {
	    if (var->sel->type == ADIOS_SELECTION_WRITEBLOCK) {
		if (num_dims == 0) { // writeblock selection for scalar
		    if (var->sel->u.block.index == writer_rank) {
			void *tmp_data = get_FMfieldAddr_by_name(f,
								 f->field_name,
								 base_data);
			memcpy(var->chunks[0].user_buf, tmp_data, f->field_size);
		    }
		}
		else { // writeblock selection for arrays
		    if (var->sel->u.block.index == writer_rank) {
			var->array_size = var->type_size;
			int i;
			for (i=0; i<num_dims; i++) {
			    char *dim;
			    atom_name[0] ='\0';
			    strcat(atom_name, FP_DIM_ATTR_NAME);
			    strcat(atom_name, "_");
			    strcat(atom_name, f->field_name);
			    strcat(atom_name, "_");
			    char dim_num[10] = "";
			    sprintf(dim_num, "%d", i+1);
			    strcat(atom_name, dim_num);
			    get_string_attr(attrs, attr_atom_from_string(atom_name), &dim);

			    FMField *temp_field = find_field_by_name(dim, f);
			    if (!temp_field) {
				adios_error(err_corrupted_variable,
					    "Could not find fieldname: %s\n",
					    dim);
			    } else {
				int *temp_data = get_FMfieldAddr_by_name
				    (
					temp_field,
					temp_field->field_name,
					base_data
				    );

				uint64_t dim = (uint64_t)(*temp_data);
				var->array_size = var->array_size * dim;
			    }
			}
			void *arrays_data  = get_FMPtrField_by_name(f,
								    f->field_name,
								    base_data,
								    1);
			memcpy(var->chunks[0].user_buf, arrays_data, var->array_size);
		    }
		}
	    }
	    else if (var->sel->type == ADIOS_SELECTION_BOUNDINGBOX) {

		if (var->ndims > 0) { // arrays
		    int i;
		    global_var *gv = find_gbl_var(fp->gp->vars,
						  var->varname,
						  fp->gp->num_vars);

		    array_displacements *disp = find_displacement(var->displ,
								  writer_rank,
								  var->num_displ);
		    if (disp) { // does this writer hold a chunk we've asked for

			//print_displacement(disp, fp->rank);

			uint64_t *global_sel_start = var->sel->u.bb.start;
			uint64_t *global_sel_count = var->sel->u.bb.count;
			uint64_t *temp = gv->offsets[0].local_dimensions;
			uint64_t *temp2 = gv->offsets[0].local_offsets;
			uint64_t *global_dimensions = gv->offsets[0].global_dimensions;
			int offsets_per_rank = gv->offsets[0].offsets_per_rank;
			uint64_t *writer_sizes = &temp[offsets_per_rank * writer_rank];
			uint64_t *writer_offsets = &temp2[offsets_per_rank * writer_rank];

			char *writer_array = (char*)get_FMPtrField_by_name(f,
									   f->field_name,
									   base_data, 1);
			char *reader_array = (char*)var->chunks[0].user_buf;

			extract_selection_from_partial(f->field_size, disp->ndims, global_dimensions,
						       writer_offsets, writer_sizes,
						       global_sel_start, global_sel_count,
						       writer_array, reader_array);
		    }
		}
	    }
	}
	else { //var has not been scheduled;
	    if (num_dims == 0) { // only worry about scalars
		flexpath_var_chunk *chunk = &var->chunks[0];
		if (!chunk->has_data) {
		    void *tmp_data = get_FMfieldAddr_by_name(f, f->field_name, base_data);
		    chunk->data = malloc(f->field_size);
		    memcpy(chunk->data, tmp_data, f->field_size);
		    chunk->has_data = 1;
		}
	    }
	}
        f++;
    }


    fp->req.num_completed++;
    /* fprintf(stderr, "\t\treader rank:%d:step:%d:num_completed:%d:num_pending:%d\n", */
    /* 	    fp->rank, fp->mystep, fp->req.num_completed, fp->req.num_pending); */
    if (fp->req.num_completed == fp->req.num_pending) {
	/* fprintf(stderr, "\t\treader rank:%d:step:%d:signalling_on:%d\n", */
	/* 	fp->rank, fp->mystep, fp->req.condition); */
        CMCondition_signal(fp_read_data->cm, fp->req.condition);
    }


    free_fmstructdesclist(struct_list);
    return 0;
}

/********** Core ADIOS Read functions. **********/

/*
 * Gathers basic MPI information; sets up reader CM.
 */
int
adios_read_flexpath_init_method (MPI_Comm comm, PairStruct* params)
{
    setenv("CMSelfFormats", "1", 1);
    fp_read_data = malloc(sizeof(flexpath_read_data));
    if (!fp_read_data) {
        adios_error(err_no_memory, "Cannot allocate memory for flexpath.");
        return -1;
    }
    memset(fp_read_data, 0, sizeof(flexpath_read_data));

    fp_read_data->CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    attr_list listen_list = NULL;
    char * transport = NULL;
    transport = getenv("CMTransport");

    // setup MPI stuffs
    fp_read_data->comm = comm;
    MPI_Comm_size(fp_read_data->comm, &(fp_read_data->size));
    MPI_Comm_rank(fp_read_data->comm, &(fp_read_data->rank));

    fp_read_data->cm = CManager_create();
    if (transport == NULL) {
      int listened = 0;
      while (listened == 0) {
	  if (CMlisten(fp_read_data->cm) == 0) {
	      fprintf(stderr, "Flexpath ERROR: reader %d:pid:%d unable to initialize connection manager. Trying again.\n",
		      fp_read_data->rank, (int)getpid());
	  } else {
	      listened = 1;
	  }
      }
    } else {
	listen_list = create_attr_list();
	add_attr(listen_list, fp_read_data->CM_TRANSPORT, Attr_String,
		 (attr_value)strdup(transport));
	CMlisten_specific(fp_read_data->cm, listen_list);
    }
    int forked = CMfork_comm_thread(fp_read_data->cm);
    if (!forked) {
	fprintf(stderr, "reader %d failed to fork comm_thread.\n", fp_read_data->rank);
	/*log_debug( "forked\n");*/
    }
    return 0;
}

ADIOS_FILE*
adios_read_flexpath_open_file(const char * fname, MPI_Comm comm)
{
    adios_error (err_operation_not_supported,
                 "FLEXPATH staging method does not support file mode for reading. "
                 "Use adios_read_open() to open a staged dataset.\n");
    return NULL;
}

/*
 * Still have work to do here.
 * Change it so that we can support the timeouts and lock_modes.
 */
/*
 * Sets up local data structure for series of reads on an adios file
 * - create evpath graph and structures
 * -- create evpath control stone (outgoing)
 * -- create evpath data stone (incoming)
 * -- rank 0 dumps contact info to file
 * -- create connections using contact info from file
 */
ADIOS_FILE*
adios_read_flexpath_open(const char * fname,
			 MPI_Comm comm,
			 enum ADIOS_LOCKMODE lock_mode,
			 float timeout_sec)
{
    fp_log("FUNC", "entering flexpath_open\n");
    ADIOS_FILE *adiosfile = malloc(sizeof(ADIOS_FILE));
    if (!adiosfile) {
	adios_error (err_no_memory,
		     "Cannot allocate memory for file info.\n");
	return NULL;
    }

    flexpath_reader_file *fp = new_flexpath_reader_file(fname);
    fp->host_language = futils_is_called_from_fortran();
    adios_errno = 0;
    fp->stone = EValloc_stone(fp_read_data->cm);
    fp->comm = comm;

    MPI_Comm_size(fp->comm, &(fp->size));
    MPI_Comm_rank(fp->comm, &(fp->rank));

    EVassoc_terminal_action(fp_read_data->cm,
			    fp->stone,
			    op_format_list,
			    op_msg_handler,
			    adiosfile);

    EVassoc_terminal_action(fp_read_data->cm,
			    fp->stone,
			    update_step_msg_format_list,
			    update_step_msg_handler,
			    adiosfile);


    EVassoc_terminal_action(fp_read_data->cm,
			    fp->stone,
			    evgroup_format_list,
			    group_msg_handler,
			    adiosfile);

    EVassoc_raw_terminal_action(fp_read_data->cm,
				fp->stone,
				raw_handler,
				adiosfile);

    /* Gather the contact info from the other readers
       and write it to a file. Create a ready file so
       that the writer knows it can parse this file. */
    char writer_ready_filename[200];
    char writer_info_filename[200];
    char reader_ready_filename[200];
    char reader_info_filename[200];

    sprintf(reader_ready_filename, "%s_%s", fname, READER_READY_FILE);
    sprintf(reader_info_filename, "%s_%s", fname, READER_CONTACT_FILE);
    sprintf(writer_ready_filename, "%s_%s", fname, WRITER_READY_FILE);
    sprintf(writer_info_filename, "%s_%s", fname, WRITER_CONTACT_FILE);

    char *string_list;
    char data_contact_info[CONTACT_LENGTH];
    string_list = attr_list_to_string(CMget_contact_list(fp_read_data->cm));
    sprintf(&data_contact_info[0], "%d:%s", fp->stone, string_list);
    free(string_list);

    char * recvbuf;
    if (fp->rank == 0) {
	recvbuf = (char*)malloc(sizeof(char)*CONTACT_LENGTH*(fp->size));
    }

    MPI_Gather(data_contact_info, CONTACT_LENGTH, MPI_CHAR, recvbuf,
	       CONTACT_LENGTH, MPI_CHAR, 0, fp->comm);

    if (fp->rank == 0) {
	// print our own contact information
	FILE * fp_out = fopen(reader_info_filename, "w");
	int i;
	if (!fp_out) {
	    adios_error(err_file_open_error,
			"File for contact info could not be opened for writing.\n");
	    exit(1);
	}
	for (i=0; i<fp->size; i++) {
	    fprintf(fp_out,"%s\n", &recvbuf[i*CONTACT_LENGTH]);
	}
	fclose(fp_out);
	free(recvbuf);

	FILE * read_ready = fopen(reader_ready_filename, "w");
	fprintf(read_ready, "ready");
	fclose(read_ready);
    }
    MPI_Barrier(fp->comm);

    FILE * fp_in = fopen(writer_ready_filename,"r");
    while (!fp_in) {
	//CMsleep(fp_read_data->cm, 1);
	fp_in = fopen(writer_ready_filename, "r");
    }
    fclose(fp_in);

    fp_in = fopen(writer_info_filename, "r");
    while (!fp_in) {
	//CMsleep(fp_read_data->cm, 1);
	fp_in = fopen(writer_info_filename, "r");
    }

    char in_contact[CONTACT_LENGTH] = "";
    //fp->bridges = malloc(sizeof(bridge_info));
    int num_bridges = 0;
    int their_stone;

    // change to read all numbers, dont create stones, turn bridge array into linked list
    while (fscanf(fp_in, "%d:%s", &their_stone, in_contact) != EOF) {
	//fprintf(stderr, "writer contact: %d:%s\n", their_stone, in_contact);
	fp->bridges = realloc(fp->bridges,
					  sizeof(bridge_info) * (num_bridges+1));
	fp->bridges[num_bridges].their_num = their_stone;
	fp->bridges[num_bridges].contact = strdup(in_contact);
	fp->bridges[num_bridges].created = 0;
	fp->bridges[num_bridges].step = 0;
	fp->bridges[num_bridges].opened = 0;
	fp->bridges[num_bridges].scheduled = 0;
	num_bridges++;
    }
    fclose(fp_in);
    fp->num_bridges = num_bridges;
    // clean up of writer's files
    MPI_Barrier(fp->comm);
    if (fp->rank == 0) {
	unlink(writer_info_filename);
	unlink(writer_ready_filename);
    }

    adiosfile->fh = (uint64_t)fp;
    adiosfile->current_step = 0;

    /* Init with a writer to get initial scalar
       data so we can handle inq_var calls and
       also populate the ADIOS_FILE struct. */
    if (fp->size < num_bridges) {
    	int mystart = (num_bridges/fp->size) * fp->rank;
    	int myend = (num_bridges/fp->size) * (fp->rank+1);
    	fp->writer_coordinator = mystart;
    	int z;
    	for (z=mystart; z<myend; z++) {
    	    build_bridge(&fp->bridges[z]);
    	}
    }
    else {
	int writer_rank = fp->rank % num_bridges;
	build_bridge(&fp->bridges[writer_rank]);
	fp->writer_coordinator = writer_rank;
    }

    // requesting initial data.
    send_open_msg(fp, fp->writer_coordinator);


    fp->data_read = 0;
    fp->req.condition = CMCondition_get(fp_read_data->cm, NULL);
    fp->req.num_pending = 1;
    fp->req.num_completed = 0;
    send_flush_msg(fp, fp->writer_coordinator, DATA, 0);
    CMCondition_wait(fp_read_data->cm, fp->req.condition);

    send_flush_msg(fp, fp->writer_coordinator, EVGROUP, 1);
    fp->data_read = 0;
    // this has to change. Writer needs to have some way of
    // taking the attributes out of the xml document
    // and sending them over ffs encoded. Not yet implemented.
    // the rest of this info for adiosfile gets filled in raw_handler.
    adiosfile->nattrs = 0;
    adiosfile->attr_namelist = NULL;
    // first step is at least one, otherwise raw_handler will not execute.
    // in reality, writer might be further along, so we might have to make
    // the writer explitly send across messages each time it calls close, to
    // indicate which timesteps are available.
    adiosfile->last_step = 1;
    adiosfile->path = strdup(fname);
    // verifies these two fields. It's not BP, so no BP version.
    // It's a stream, so how can the file size be known?
    adiosfile->version = -1;
    adiosfile->file_size = 0;
    adios_errno = err_no_error;
    fp_log("FUNC", "leaving flexpath_open\n");
    return adiosfile;
}

int adios_read_flexpath_finalize_method ()
{
    return 0;
}

void adios_read_flexpath_release_step(ADIOS_FILE *adiosfile) {
    int i;
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    MPI_Barrier(fp->comm);
    for (i=0; i<fp->num_bridges; i++) {
        if (fp->bridges[i].created && !fp->bridges[i].opened) {
	    send_open_msg(fp, i);
        }
    }
    free_evgroup(fp->gp);
    fp->gp = NULL;

    flexpath_var *tmpvars = fp->var_list;
    while (tmpvars) {

	if (tmpvars->ndims > 0) {
	    free(tmpvars->global_dims);
	    tmpvars->ndims = 0;
	}
	if (tmpvars->displ) {
	    free_displacements(tmpvars->displ, tmpvars->num_displ);
	    tmpvars->displ = NULL;
	}

	if (tmpvars->sel) {
	    common_read_selection_delete(tmpvars->sel);
	    tmpvars->sel = NULL;
	}

	tmpvars->sel = NULL;
	for (i=0; i<tmpvars->num_chunks; i++) {
	    flexpath_var_chunk *chunk = &tmpvars->chunks[i];
	    if (chunk->has_data) {
		free(chunk->data);
		chunk->data = NULL;
		chunk->has_data = 0;
	    }
	    chunk->rank = 0;
	}
	tmpvars = tmpvars->next;
    }
}

int
adios_read_flexpath_advance_step(ADIOS_FILE *adiosfile, int last, float timeout_sec)
{    
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    //fprintf(stderr, "reader_rank:%d:calling advance_step:step:%d\n", fp->rank, fp->mystep);
    MPI_Barrier(fp->comm);
    int count = 0; // for perf measurements
    //fprintf(stderr, "reader_rank:%d:advance_step:sending flush step to coordinator:%d:step:%d\n", fp->rank, fp->writer_coordinator,fp->mystep);
    send_flush_msg(fp, fp->writer_coordinator, STEP, 1);
    //fprintf(stderr, "reader_rank:%d:advance_step:advance_step:after sending flush step to:%d:step:%d\n", fp->rank, fp->writer_coordinator,fp->mystep);
    //put this on a timer, so to speak, for timeout_sec
    while (fp->mystep == fp->last_writer_step) {
	if (fp->writer_finalized) {
	    adios_errno = err_end_of_stream;
	    return err_end_of_stream;
	}
	CMsleep(fp_read_data->cm, 1);
	//fprintf(stderr, "reader_rank:%d:advance_step:while loop sending flush step to coordinator:%d:step:%d\n", fp->rank, fp->writer_coordinator,fp->mystep);
	send_flush_msg(fp, fp->writer_coordinator, STEP, 1);
	//fprintf(stderr, "reader_rank:%d:advance_step:after while loop sending flush step to coordinator:%d:step:%d\n", fp->rank, fp->writer_coordinator,fp->mystep);
    }

    int i=0;
    for (i=0; i<fp->num_bridges; i++) {
        if (fp->bridges[i].created && fp->bridges[i].opened) {
	    count++;
	    send_close_msg(fp, i);
	}
    }
    MPI_Barrier(fp->comm);

    count = 0;
    adiosfile->current_step++;
    fp->mystep = adiosfile->current_step;

    for (i=0; i<fp->num_bridges; i++) {
	if (fp->bridges[i].created && !fp->bridges[i].opened) {
	    send_open_msg(fp, i);
	    count++;
        }
    }
    // need to remove selectors from each var now.
    fp->req.num_completed = 0;
    fp->req.num_pending = 1;
    fp->req.condition = CMCondition_get(fp_read_data->cm, NULL);
    //fprintf(stderr, "reader_rank:%d:advance_step:while loop sending flush DATA to coordinator:%d:step:%d\n", fp->rank, fp->writer_coordinator,fp->mystep);
    send_flush_msg(fp, fp->writer_coordinator, DATA, 0);
    CMCondition_wait(fp_read_data->cm, fp->req.condition);
    //fprintf(stderr, "reader_rank:%d:advance_step:while loop after flush DATA to coordinator:%d:step:%d\n", fp->rank, fp->writer_coordinator,fp->mystep);
    // should only happen if there are more steps available.
    // writer should have advanced.
    //fprintf(stderr, "reader_rank:%d:advance_step:while loop sending flush EVGROUP to coordinator:%d:step:%d\n", fp->rank, fp->writer_coordinator,fp->mystep);
    send_flush_msg(fp, fp->writer_coordinator, EVGROUP, 1);
    //fprintf(stderr, "reader_rank:%d:advance_step:while loop after flush EVGROUP to coordinator:%d:step:%d\n", fp->rank, fp->writer_coordinator,fp->mystep);
    return 0;
}

int adios_read_flexpath_close(ADIOS_FILE * fp)
{
    flexpath_reader_file *file = (flexpath_reader_file*)fp->fh;
    //send to each opened link
    int i;
    for (i = 0; i<file->num_bridges; i++) {
        if (file->bridges[i].created && file->bridges[i].opened) {
	    send_close_msg(file, i);
        }
    }
    /*
    start to cleanup.  Clean up var_lists for now, as the
    data has already been copied over to ADIOS_VARINFO structs
    that the user maintains a copy of.
    */
    flexpath_var *v = file->var_list;
    while (v) {
    	// free chunks; data has already been copied to user
    	int i;
    	for (i = 0; i<v->num_chunks; i++) {
    	    flexpath_var_chunk *c = &v->chunks[i];
	    if (!c)
		log_error("FLEXPATH: %s This should not happen! line %d\n",__func__,__LINE__);
	    //free(c->data);
	    c->data = NULL;
	    free(c);
	}
	flexpath_var *tmp = v->next;
	free(v);
	v = tmp;
    	//v=v->next;
    }
    return 0;
}

ADIOS_FILE *adios_read_flexpath_fopen(const char *fname, MPI_Comm comm) {
   return 0;
}

int adios_read_flexpath_is_var_timed(const ADIOS_FILE* fp, int varid) { return 0; }

void adios_read_flexpath_get_groupinfo(
    const ADIOS_FILE *adiosfile,
    int *ngroups,
    char ***group_namelist,
    uint32_t **nvars_per_group,
    uint32_t **nattrs_per_group)
{
    flexpath_reader_file *fp;
    if (adiosfile) {
        fp = (flexpath_reader_file *) adiosfile->fh;
        *ngroups = 1;
        *group_namelist = (char **) malloc (sizeof (char*));
        *group_namelist[0] = strdup (fp->group_name);
    }

}

int adios_read_flexpath_check_reads(const ADIOS_FILE* fp, ADIOS_VARCHUNK** chunk) { log_debug( "flexpath:adios function check reads\n"); return 0; }

int adios_read_flexpath_perform_reads(const ADIOS_FILE *adiosfile, int blocking)
{
    fp_log("FUNC", "entering perform_reads.\n");
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    fp->data_read = 0;
    int i;
    int batchcount = 0;
    int num_sendees = fp->num_sendees;
    int total_sent = 0;
    fp->time_in = 0.00;
    fp->req.num_completed = 0;
    fp->req.num_pending = FP_BATCH_SIZE;
    fp->req.condition = CMCondition_get(fp_read_data->cm, NULL);

    for (i = 0; i<num_sendees; i++) {
	int sendee = fp->sendees[i];
	batchcount++;
	total_sent++;
		/* fprintf(stderr, "reader rank:%d:flush_data to writer:%d:of:%d:step:%d:batch:%d:total_sent:%d\n", */
		/* 	fp->rank, sendee, num_sendees, fp->mystep, batchcount, total_sent); */
	send_flush_msg(fp, sendee, DATA, 0);

	if ((total_sent % FP_BATCH_SIZE == 0) || (total_sent == num_sendees)) {
            fp->req.num_pending = batchcount;

	    /* fprintf(stderr, "\t\treader rank:%d:blocking on:%d:step:%d\n", */
	    /* 	    fp->rank, fp->req.condition, fp->mystep); */
            CMCondition_wait(fp_read_data->cm, fp->req.condition);
	    /* fprintf(stderr, "\t\treader rank:%d:after blocking:%d:step:%d\n", */
	    /* 	    fp->rank, fp->req.condition, fp->mystep); */
	    fp->req.num_completed = 0;
	    //fp->req.num_pending = 0;
	    //total_sent = 0;
            batchcount = 0;
            fp->req.condition = CMCondition_get(fp_read_data->cm, NULL);
	}

    }

    free(fp->sendees);
    fp->sendees = NULL;
    fp->num_sendees = 0;

    flexpath_var *tmpvars = fp->var_list;
    while (tmpvars) {
	if (tmpvars->displ) {
	    free_displacements(tmpvars->displ, tmpvars->num_displ);
	    tmpvars->displ = NULL;
	}

	if (tmpvars->sel) {
	    common_read_selection_delete(tmpvars->sel);
	    tmpvars->sel = NULL;
	}

	tmpvars = tmpvars->next;
    }
    fp_log("FUNC", "leaving perform_reads.\n");
    return 0;
}

int
adios_read_flexpath_inq_var_blockinfo(const ADIOS_FILE* fp,
				      ADIOS_VARINFO* varinfo)
{ /*log_debug( "flexpath:adios function inq var block info\n");*/ return 0; }

int
adios_read_flexpath_inq_var_stat(const ADIOS_FILE* fp,
				 ADIOS_VARINFO* varinfo,
				 int per_step_stat,
				 int per_block_stat)
{ /*log_debug( "flexpath:adios function inq var stat\n");*/ return 0; }


int
adios_read_flexpath_schedule_read_byid(const ADIOS_FILE *adiosfile,
				       const ADIOS_SELECTION *sel,
				       int varid,
				       int from_steps,
				       int nsteps,
				       void *data)
{
    fp_log("FUNC", "entering schedule_read_byid\n");
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    flexpath_var *fpvar = fp->var_list;

    while (fpvar) {
        if (fpvar->id == varid)
        	break;
        else
	    fpvar=fpvar->next;
    }
    if (!fpvar) {
        adios_error(err_invalid_varid,
		    "Invalid variable id: %d\n",
		    varid);
        return err_invalid_varid;
    }
    //store the user allocated buffer.
    flexpath_var_chunk *chunk = &fpvar->chunks[0];
    if (nsteps != 1) {
	adios_error (err_invalid_timestep,
                     "Only one step can be read from a stream at a time. "
                     "You requested % steps in adios_schedule_read()\n",
		     nsteps);
        return err_invalid_timestep;
    }
    // this is done so that the user can do multiple schedule_read/perform_reads
    // within before doing release/advance step. Might need a better way to
    // manage the ADIOS selections.
    if (fpvar->sel) {
	common_read_selection_delete(fpvar->sel);
	fpvar->sel = NULL;
    }
    if (!sel) { // null selection; read whole variable
	//TODO: This will have to be fixed for local arrays,
	// but dataspaces doesn't have local arrays so there
	// are no use cases for it.
	uint64_t *starts = calloc(fpvar->ndims, sizeof(uint64_t));
	uint64_t *counts = calloc(fpvar->ndims, sizeof(uint64_t));
	memcpy(counts, fpvar->global_dims, fpvar->ndims*sizeof(uint64_t));
	fpvar->sel = common_read_selection_boundingbox(fpvar->ndims, starts, counts);
    } else {
	fpvar->sel = copy_selection(sel);
    }

    switch(fpvar->sel->type)
    {
    case ADIOS_SELECTION_WRITEBLOCK:
    {
        chunk->user_buf = data;
        fpvar->start_position = 0;
	int writer_index = fpvar->sel->u.block.index;
	if (writer_index > fp->num_bridges) {
	    adios_error(err_out_of_bound,
			"No process exists on the writer side matching the index.\n");
	    return err_out_of_bound;
	}
	send_var_message(fp, writer_index, fpvar->varname);
	break;
    }
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        // boundingbox for a scalar; handle it as we do with inq_var
        if (fpvar->ndims == 0) {
            memcpy(data, chunk->data, fpvar->type_size);
        } else {
            if (fp->host_language == FP_FORTRAN_MODE) {
                reverse_dims(sel->u.bb.start, sel->u.bb.ndim);
                reverse_dims(sel->u.bb.count, sel->u.bb.ndim);
            }
            chunk->user_buf = data;
            fpvar->start_position = 0;
            free_displacements(fpvar->displ, fpvar->num_displ);
            fpvar->displ = NULL;
            int j=0;
            int need_count = 0;
            array_displacements *all_disp = NULL;
            uint64_t pos = 0;
            for (j=0; j<fp->num_bridges; j++) {
                int destination=0;
                if (need_writer(fp, j, fpvar->sel, fp->gp, fpvar->varname)==1) {
		    //printf("\t\trank: %d need_writer: %d\n", fp->rank, j);
                    uint64_t _pos = 0;
                    need_count++;
                    destination = j;
                    global_var *gvar = find_gbl_var(fp->gp->vars, fpvar->varname, fp->gp->num_vars);
                    // displ is freed in release_step.
                    array_displacements *displ = get_writer_displacements(j, fpvar->sel, gvar, &_pos);
                    displ->pos = pos;
                    _pos *= (uint64_t)fpvar->type_size;
                    pos += _pos;

                    all_disp = realloc(all_disp, sizeof(array_displacements)*need_count);
                    all_disp[need_count-1] = *displ;
                    send_var_message(fp, j, fpvar->varname);
                }
            }
            fpvar->displ = all_disp;
            fpvar->num_displ = need_count;
        }
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
    fp_log("FUNC", "entering schedule_read_byid\n");
    return 0;
}

int
adios_read_flexpath_schedule_read(const ADIOS_FILE *adiosfile,
			const ADIOS_SELECTION * sel,
			const char * varname,
			int from_steps,
			int nsteps,
			void * data)
{
    fprintf(stderr, "schedule_read is called\n");
    return 0;
}

int
adios_read_flexpath_get_attr (int *gp, const char *attrname,
                                 enum ADIOS_DATATYPES *type,
                                 int *size, void **data)
{
    //log_debug( "debug: adios_read_flexpath_get_attr\n");
    // TODO: borrowed from dimes
    adios_error(err_invalid_read_method,
		"adios_read_flexpath_get_attr is not implemented.");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}

int
adios_read_flexpath_get_attr_byid (const ADIOS_FILE *adiosfile, int attrid,
				   enum ADIOS_DATATYPES *type,
				   int *size, void **data)
{
//    log_debug( "debug: adios_read_flexpath_get_attr_byid\n");
    // TODO: borrowed from dimes
    adios_error(err_invalid_read_method,
		"adios_read_flexpath_get_attr_byid is not implemented.");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}

ADIOS_VARINFO*
adios_read_flexpath_inq_var(const ADIOS_FILE * adiosfile, const char* varname)
{
    fp_log("FUNC", "entering flexpath_inq_var\n");
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    ADIOS_VARINFO *v = NULL;

    flexpath_var *fpvar = find_fp_var(fp->var_list, varname);
    if (fpvar) {
        v = calloc(1, sizeof(ADIOS_VARINFO));

        if (!v) {
            adios_error(err_no_memory,
                        "Cannot allocate buffer in adios_read_flexpath_inq_var()");
            return NULL;
        }

	v = convert_var_info(fpvar, v, varname, adiosfile);
	fp_log("FUNC", "leaving flexpath_inq_var\n");
    }
    else {
        adios_error(err_invalid_varname, "Cannot find var %s\n", varname);
    }
    return v;
}

ADIOS_VARINFO*
adios_read_flexpath_inq_var_byid (const ADIOS_FILE * adiosfile, int varid)
{
    fp_log("FUNC", "entering flexpath_inq_var_byid\n");
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    if (varid >= 0 && varid < adiosfile->nvars) {
	ADIOS_VARINFO *v = adios_read_flexpath_inq_var(adiosfile, adiosfile->var_namelist[varid]);
	fp_log("FUNC", "leaving flexpath_inq_var_byid\n");
	return v;
    }
    else {
        adios_error(err_invalid_varid, "FLEXPATH method: Cannot find var %d\n", varid);
        return NULL;
    }
}

void
adios_read_flexpath_free_varinfo (ADIOS_VARINFO *adiosvar)
{
    //log_debug( "debug: adios_read_flexpath_free_varinfo\n");
    fprintf(stderr, "adios_read_flexpath_free_varinfo called\n");
    return;
}


ADIOS_TRANSINFO*
adios_read_flexpath_inq_var_transinfo(const ADIOS_FILE *gp, const ADIOS_VARINFO *vi)
{
    //adios_error(err_operation_not_supported, "Flexpath does not yet support transforms: var_transinfo.\n");
    ADIOS_TRANSINFO *trans = malloc(sizeof(ADIOS_TRANSINFO));
    memset(trans, 0, sizeof(ADIOS_TRANSINFO));
    trans->transform_type = adios_transform_none;
    return trans;
}


int
adios_read_flexpath_inq_var_trans_blockinfo(const ADIOS_FILE *gp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti)
{
    adios_error(err_operation_not_supported, "Flexpath does not yet support transforms: trans_blockinfo.\n");
    return (int64_t)0;
}

int  adios_read_flexpath_get_dimension_order (const ADIOS_FILE *adiosfile)
{
    flexpath_reader_file *fp = (flexpath_reader_file*)adiosfile->fh;
    return (fp->host_language == FP_FORTRAN_MODE);
}

void
adios_read_flexpath_reset_dimension_order (const ADIOS_FILE *adiosfile, int is_fortran)
{
    //log_debug( "debug: adios_read_flexpath_reset_dimension_order\n");
    adios_error(err_invalid_read_method, "adios_read_flexpath_reset_dimension_order is not implemented.");
}
