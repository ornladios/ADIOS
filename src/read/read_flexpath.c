/*
    Read_flexpath.c
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
#include <stdarg.h>

// evpath libraries
#include <ffs.h>
#include <atl.h>
//#include <gen_thread.h>
#include <evpath.h>

//#ifndef _NOMPI
#include <data_store.h>
#include <data_store_globals.h>
#include <container_replica.h>
//#include "core/globals.h"
//#endif

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

// local libraries
#include "config.h"
#include "public/adios.h"
#include "public/adios_types.h"
#include "core/globals.h"

#include "public/adios_read_v2.h"
#include "core/adios_read_hooks.h"
#include "core/adios_logger.h"
#include "public/adios_error.h"
#include "core/flexpath.h"

#include "core/transforms/adios_transforms_common.h" // NCSU ALACRITY-ADIOS

#ifdef _NOMPI
#include <mpirelay_client.h>
/* #else */
/* #include "public/adios_mpi.h" */
#endif


// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define FP_BATCH_SIZE 4

#define _MAX(x, y) (((x) > (y)) ? (x) : (y))
#define _MIN(x, y) (((x) < (y)) ? (x) : (y))

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

typedef struct _Bridge_info
{
    EVstone bridge_stone;
    EVsource flush_source;
    EVsource var_source;
    EVsource op_source;
    int their_num;
    char * contact;
    int created;
    int opened;
    int step;
    int scheduled;
}Bridge_info;

typedef struct _Flexpath_reader_replica
{
    char *container_name;
    int replica_id;
    int writer_coordinator;
    int num_bridges; // number of corresponding writers.
    int active;
    
    int mystep;
    int count;
    int writer_finalized;
    int last_writer_step;

    Bridge_info *bridges;
    struct _Flexpath_reader_replica *next;
}Flexpath_reader_replica;

typedef struct _Flexpath_var_chunk
{
    int has_data;
    int rank;
    void *data;
    void *user_buf;
    uint64_t *local_bounds; // nodims
    uint64_t *global_bounds; // ndims
    uint64_t *global_offsets; // ndims
} Flexpath_var_chunk;

typedef struct _Flexpath_var
{
    int id;
    char *varname;
    char *varpath;

    enum ADIOS_DATATYPES type;
    uint64_t type_size; // type size, not arrays size

    int num_dims;
    uint64_t *global_dims; // ndims size (if ndims>0)
    uint64_t *local_dims; // for local arrays
    uint64_t array_size; // not relevant for scalars

    int num_chunks;
    Flexpath_var_chunk *chunks;

    int num_displ;
    array_displacements *displ;

    ADIOS_SELECTION *sel;
    uint64_t start_position;

    struct _Flexpath_var *next;
} Flexpath_var;

typedef struct _Flexpath_reader_file
{
    char *file_name;
    char *group_name; // assuming one group per file right now.
    int host_language;
    
    char *container_name;
    int replica_id;
    EVstone stone;

#ifndef _NOMPI
    MPI_Comm comm;
#else
    MPIRelay_client *client;
#endif
    int rank;
    int size;
    int valid;

    pthread_mutex_t repmutex;
    Flexpath_reader_replica *replicas;
    Flexpath_reader_replica *current;

    int num_vars;
    Flexpath_var * var_list;
    int z;
    int num_gp;
    evgroup * gp;

    int num_sendees;
    int *sendees;
    int ackCondition;

    int pending_requests;
    int completed_requests;
    uint64_t data_read; // for perf measurements.
    double time_in; // for perf measurements.
    pthread_mutex_t data_mutex;
    pthread_cond_t data_condition;
} Flexpath_reader_file;

typedef struct _local_read_data
{
    // MPI stuff
#ifndef _NOMPI
    MPI_Comm comm;
#else
    MPIRelay_client *client;
#endif
    int rank;
    int size;

    int port;
    int tunnel;
    // EVPath stuff
    CManager cm;
    EVstone stone;
    atom_t CM_TRANSPORT;
//#ifndef _NOMPI
    Replica_context *ctx;
//#endif
} Flexpath_reader_local;

static Flexpath_reader_local *local = NULL;
static int global_step = 0;
/********** Helper functions. **********/

static double dgettimeofday( void )
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

void
fp_reader_log(const char *log, const char *format, ...)
{
    //if (local->rank != 15) return;

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

        sprintf(header, "READER:%s", log);
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


static Update_contact_msg*
exchange_update_msg(Flexpath_reader_file *file_data, Update_contact_msg *op)
{
#ifndef _NOMPI
    MPI_Barrier(file_data->comm);
#else
    MPIRelay_barrier(file_data->client);
#endif
    
    double start = MPI_Wtime();
    if (file_data->rank != 0) {
	op = malloc(sizeof(Update_contact_msg));
    }
    
#ifndef _NOMPI
    MPI_Bcast(&op->num_replicas, 1, MPI_INT, 0, file_data->comm);
#else
    MPIRelay_bcast(file_data->client, &op->num_replicas, 1, MPIInt, 0);
#endif
    
    fp_reader_log("BIGBIG", "reader got num_replicas: %d\n", op->num_replicas);
    if (file_data->rank != 0) {
	op->replicas = malloc(sizeof(Replica_info_msg)*op->num_replicas);
	memset(op->replicas, 0, sizeof(Replica_info_msg)*op->num_replicas);
    }
    int i;
    for (i = 0; i < op->num_replicas; i++) {
#ifndef _NOMPI
	MPI_Bcast(&op->replicas[i].replica_id, 1, MPI_INT, 0, file_data->comm);
	MPI_Bcast(&op->replicas[i].num_writers, 1, MPI_INT, 0, file_data->comm);
#else
	MPIRelay_bcast(file_data->client, &op->replicas[i].replica_id, 1, MPIInt, 0);
	MPIRelay_bcast(file_data->client, &op->replicas[i].num_writers, 1, MPIInt, 0);
#endif
	//printf ("num_writers: %d\n", op->replicas[i].num_writers);
	int num_writers = op->replicas[i].num_writers;
	int bufsize = sizeof(char) * num_writers * CONTACT_LENGTH;
	char *upstream_writers = malloc(bufsize);
	memset(upstream_writers, 0, bufsize);

	if(file_data->rank == 0) {
	    int pos = 0;
	    int j;
	    for (j = 0; j<num_writers; j++) {
		strncpy(&upstream_writers[pos], 
                        op->replicas[i].writers[j].endpoint, 
                        CONTACT_LENGTH);
		pos += CONTACT_LENGTH;
	    }
	}
#ifndef _NOMPI
	MPI_Bcast(upstream_writers, bufsize, MPI_BYTE, 0, file_data->comm);
#else
	MPIRelay_bcast(file_data->client, upstream_writers, bufsize, MPIByte, 0);	   
#endif
	
	if (file_data->rank != 0) {
	    op->replicas[i].writers = malloc(sizeof(Contact_info) * num_writers);
	    int j;
	    int pos = 0;

	    for (j = 0; j<num_writers; j++) {
		op->replicas[i].writers[j].endpoint = strndup(&upstream_writers[pos], CONTACT_LENGTH);
		pos += CONTACT_LENGTH;
	    }
	}
    }
    if (file_data->rank == 0) {
        double end = MPI_Wtime();
        printf("protocol:%lf:total:%lf:replica_id:%d:reader_dispersed\n", end, end-start, file_data->replica_id);
    }
    return op;
}
//#endif

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

void
build_bridge(Bridge_info* bridge) {
    attr_list contact_list = attr_list_from_string(bridge->contact);
    //printf("\tbridge to: %d:%s\n", bridge->their_num, bridge->contact);
    if (bridge->created == 0) {
	bridge->bridge_stone =
	    EVcreate_bridge_action(local->cm,
				   contact_list,
				   (EVstone)bridge->their_num);

	bridge->flush_source =
	    EVcreate_submit_handle(local->cm,
				   bridge->bridge_stone,
				   flush_format_list);

	bridge->var_source =
	    EVcreate_submit_handle(local->cm,
				   bridge->bridge_stone,
				   var_format_list);

	bridge->op_source =
	    EVcreate_submit_handle(local->cm,
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
    EVreturn_event_buffer(local->cm, gp);
}

Flexpath_var*
new_flexpath_var(const char *varname, int id, uint64_t type_size)
{
    Flexpath_var *var = malloc(sizeof(Flexpath_var));
    if (var == NULL) {
	log_error("Error creating new var: %s\n", varname);
	return NULL;
    }

    memset(var, 0, sizeof(Flexpath_var));
    // free this when freeing vars.
    var->varname = strdup(varname);
    var->id = id;
    var->type_size = type_size;
    var->displ = NULL;
    return var;
}

Flexpath_reader_file*
new_Flexpath_reader_file(const char *fname)
{
    Flexpath_reader_file * fp = malloc(sizeof(Flexpath_reader_file));
    if (fp == NULL) {
	log_error("Cannot create data for new file.\n");
	exit(1);
    }
    memset(fp, 0, sizeof(Flexpath_reader_file));
    fp->file_name = strdup(fname);
    //fp->last_writer_step = -1;

    pthread_mutex_init(&fp->repmutex, NULL);
    pthread_mutex_init(&fp->data_mutex, NULL);
    pthread_cond_init(&fp->data_condition, NULL);
    return fp;
}

enum ADIOS_DATATYPES
ffs_type_to_adios_type(const char *ffs_type)
{
    if (!strcmp("integer", ffs_type))
	return adios_integer;
    else if (!strcmp("float", ffs_type))
	return adios_real;
    else if (!strcmp("string", ffs_type))
	return adios_string;
    else if (!strcmp("double", ffs_type))
	return adios_double;
    else if (!strcmp("char", ffs_type))
	return adios_byte;
    else
	return adios_unknown;
}

ADIOS_VARINFO*
convert_var_info(Flexpath_var * fpvar,
		 ADIOS_VARINFO * v,
		 const char* varname,
		 const ADIOS_FILE *adiosfile)
{
    int i;
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    v->type = fpvar->type;
    v->ndim = fpvar->num_dims;
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
	if (!v->value) {
	    adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
	    return NULL;
	}
	Flexpath_var_chunk * chunk = &fpvar->chunks[0];
	memcpy(v->value, chunk->data, value_size);
	v->global = 0;
    } else { // arrays
	v->dims = (uint64_t*)malloc(v->ndim * sizeof(uint64_t));
	if (!v->dims) {
	    adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
	    return NULL;
	}
	// broken.  fix.
	int cpysize = fpvar->num_dims*sizeof(uint64_t);
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

Flexpath_var *
find_fp_var(Flexpath_var * var_list, const char * varname)
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
	else{
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

/*
 * Finds the array displacements for a writer identified by its rank.
 */
array_displacements*
find_displacement(array_displacements* list, int rank, int num_displ) {
    int i;
    for (i=0; i<num_displ; i++) {
	if (list[i].writer_rank == rank)
	    return &list[i];
    }
    return NULL;
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
map_local_to_global_index(uint64_t ndim,
                          uint64_t *local_index,
                          uint64_t *local_offsets,
                          uint64_t *global_index)
{
    int i;
    for (i=0; i < ndim; i++) {
        global_index[i] = local_index[i] + local_offsets[i];
    }
}
void
map_global_to_local_index(uint64_t ndim,
                          uint64_t *global_index,
                          uint64_t *local_offsets,
                          uint64_t *local_index)
{
    int i;
    for (i=0; i < ndim; i++) {
        local_index[i] = global_index[i] - local_offsets[i];
    }
}
int
index_in_selection(uint64_t ndim,
                   uint64_t *global_index,
                   uint64_t *selection_offsets,
                   uint64_t *selection_counts)
{
    int i;
    for (i=0; i < ndim; i++) {
        if ((global_index[i] < selection_offsets[i]) ||
            (global_index[i] >= selection_offsets[i] + selection_counts[i])) {
            return 0;
        }
    }
    return 1;
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
	    int left = _MAX(partial_offsets[dim], selection_offsets[dim]);
	    int right = _MIN(partial_offsets[dim] + partial_counts[dim],
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
	int left = _MAX(partial_offsets[dim], selection_offsets[dim]);
	int right = _MIN(partial_offsets[dim] + partial_counts[dim],
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


	}else{
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
    Flexpath_reader_file *fp,
    int writer,
    const ADIOS_SELECTION *sel,
    evgroup_ptr gp,
    char* varname)
{
    //select var from group
    global_var *gvar = find_gbl_var(gp->vars, varname, gp->num_vars);
    if (!gvar)
	printf("GVAR IS NULL!!!\n");
    //for each dimension
    int i=0;
    offset_struct var_offsets = gvar->offsets[0];
    for (i=0; i< var_offsets.offsets_per_rank; i++) {
	int pos = writer*(var_offsets.offsets_per_rank) + i;

        uint64_t sel_offset = sel->u.bb.start[i];
        uint64_t sel_size = sel->u.bb.count[i];

        uint64_t rank_offset = var_offsets.local_offsets[pos];
        uint64_t rank_size = var_offsets.local_dimensions[pos];

        if ((rank_offset <= sel_offset) && (rank_offset + rank_size - 1 >=sel_offset)) {
	     log_debug("matched overlap type 1\n");
        }

        else if ((rank_offset <= sel_offset + sel_size - 1) && \
		(rank_offset+rank_size>=sel_offset+sel_size-1)) {
        }else if ((sel_offset <= rank_offset) && (rank_offset+rank_size<= sel_offset+sel_size-1)) {
        }else {
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
    FMField *temp = f;
    while (temp->field_name) {
    	//free(temp->field_name);
        //free(temp->field_type);
    	temp++;
    }
    free(f);
    free(struct_list[0].opt_info);
    free(struct_list->format_name);
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

Flexpath_var*
setup_flexpath_vars(FMField *f, int *num)
{
    Flexpath_var *vars = NULL;
    int var_count = 0;

    while (f->field_name != NULL) {
	Flexpath_var *curr_var = new_flexpath_var(f->field_name,
						  var_count,
						  f->field_size);
	curr_var->num_chunks = 1;
	curr_var->chunks =  malloc(sizeof(Flexpath_var_chunk)*curr_var->num_chunks);
	memset(curr_var->chunks, 0, sizeof(Flexpath_var_chunk)*curr_var->num_chunks);
	curr_var->sel = NULL;
	curr_var->type = ffs_type_to_adios_type(f->field_type);
	Flexpath_var *temp = vars;
	curr_var->next = temp;
	vars = curr_var;
	var_count++;
	f++;
    }
    *num = var_count;
    return vars;
}

/*****************Messages to writer procs**********************/

void
send_open_msg(Flexpath_reader_file *fp, Flexpath_reader_replica *rep, int destination)
{
    fp_reader_log("DETAILS", 
		  "replica_id: %d rank: %d entering send_open step: %d to writer_replica: %d writer_rank: %d\n",
		  fp->replica_id, fp->rank, rep->mystep, rep->replica_id, destination);
    if (!rep->bridges[destination].created) {
	build_bridge(&(rep->bridges[destination]));
    }
    op_msg msg;
    msg.process_id = fp->rank;
    msg.replica_id = fp->replica_id;
    msg.file_name = fp->file_name;
    msg.step = rep->mystep;
    msg.type = OPEN_MSG;
    int cond = CMCondition_get(local->cm, NULL);
    msg.condition = cond;

    EVsubmit(rep->bridges[destination].op_source, &msg, NULL);
    CMCondition_wait(local->cm, cond);
    rep->bridges[destination].opened = 1;
    fp_reader_log("DETAILS", 
		  "replica_id: %d rank: %d after open step: %d to writer_replica: %d writer_rank: %d\n",
		  fp->replica_id, fp->rank, rep->mystep, rep->replica_id, destination);
}

void
send_close_msg(Flexpath_reader_file *fp,
	       Flexpath_reader_replica *rep,
	       int destination)
{
    if (!rep->bridges[destination].created) {
	build_bridge(&(rep->bridges[destination]));
    }
    op_msg msg;
    msg.process_id = fp->rank;
    msg.replica_id = fp->replica_id;
    msg.file_name = fp->file_name;
    msg.step = rep->mystep;
    msg.type = CLOSE_MSG;
    //msg.condition = -1;
    int cond = CMCondition_get(local->cm, NULL);
    msg.condition = cond;
    fp_reader_log(
	"DETAILS",
	"replica_id: %d rank: %d sending close step: %d with condition: %d to writer_replica: %d writer_rank: %d\n",
	fp->replica_id, fp->rank, rep->mystep, cond, rep->replica_id, destination);

    EVsubmit(rep->bridges[destination].op_source, &msg, NULL);
    CMCondition_wait(local->cm, cond);
    rep->bridges[destination].opened = 0;
    fp_reader_log("DETAILS", 
		  "replica_id: %d rank: %d after close %d to writer_replica: %d writer_rank: %d\n",
		  fp->replica_id, fp->rank, rep->mystep, rep->replica_id, destination);
}

void
send_flush_msg(Flexpath_reader_file *fp,
	       Flexpath_reader_replica *rep,
	       int destination,
	       Flush_type type,
	       int use_condition)
{
    Flush_msg msg;
    msg.type = type;
    msg.process_id = fp->rank;
    msg.replica_id = fp->replica_id;
    msg.step = rep->mystep;
    msg.id = rep->mystep;

    char *ftype = NULL;
    if (type == DATA) {
	ftype = "DATA";
    } else if (type == EVGROUP) {
	ftype = "EVGROUP";
    } else if (type == STEP) {
	ftype = "STEP";
    } else {
	ftype = "ERROR";
    }

    int condition = -1;
    if (use_condition)
	condition = CMCondition_get(local->cm, NULL);
    msg.condition = condition;

    if (!rep->bridges[destination].created) {
	build_bridge(&(rep->bridges[destination]));
    }
    // maybe check to see if the bridge is create first.
    fp_reader_log(
	"DETAILS", 
	"replica_id: %d rank: %d sending flush of type: %s to writer_replica: %d writer_rank: %d with step: %d finalized: %d count: %d\n",
	fp->replica_id, 
        fp->rank, 
        ftype, 
        rep->replica_id, 
        destination, 
        rep->mystep, 
        rep->writer_finalized, 
        rep->count);

    EVsubmit(rep->bridges[destination].flush_source, &msg, NULL);
    if (use_condition) {
	CMCondition_wait(local->cm, msg.condition);
    }

    fp_reader_log(
	"DETAILS",
	"\t\treplica_id: %d rank: %d after flush of type %s to writer_replica: %d writer_rank: %d with step: %d\n",
	fp->replica_id, fp->rank, ftype, rep->replica_id, destination, rep->mystep);
}

void
send_var_message(Flexpath_reader_file *fp, Flexpath_reader_replica *rep, int destination, char *varname)
{
    fp_reader_log("DETAILS", 
		  "replica_id: %d rank: %d sending var to writer_replica: %d writer_rank: %d with step: %d\n",
		  fp->replica_id,
		  fp->rank,
		  rep->replica_id,
		  destination,
		  rep->mystep);
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
    if (!rep->bridges[destination].created) {
	build_bridge(&(rep->bridges[destination]));
    }
    if (!rep->bridges[destination].opened) {
	rep->bridges[destination].opened = 1;
	send_open_msg(fp, rep, destination);
    }

    Var_msg var;
    var.process_id = fp->rank;
    var.replica_id = fp->replica_id;
    var.var_name = varname;
    var.step = rep->mystep;
    EVsubmit(rep->bridges[destination].var_source, &var, NULL);
    fp_reader_log("DETAILS", 
		  "replica_id: %d rank: %d after send var to writer_replica: %d writer_rank: %d with step: %d\n",
		  fp->replica_id, fp->rank, rep->replica_id, destination, rep->mystep);
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
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    fp_reader_log(
	"DETAILS",
        "\treplica_id: %d rank: %d got step from writer_replica: %d writer_rank: %d count: %d writer_step: %d mystep: %d finalized: %d\n",
        fp->replica_id,
        fp->rank, 
	msg->replica_id,
        msg->process_id,
        msg->step_count,
        msg->step,
        fp->current->mystep,
        msg->finalized);
    int condition = -1;
    get_int_attr(attrs, attr_atom_from_string("fp_dst_condition"), &condition);

    fp->current->last_writer_step = msg->step;
    fp->current->count = msg->step_count;
    fp->current->writer_finalized = msg->finalized;
    CMCondition_signal(local->cm, condition);
    return 0;
}

static int
op_msg_handler(CManager cm,
	       void *vevent,
	       void *client_data,
	       attr_list attrs)
{
    int condition = -1;
    get_int_attr(attrs, attr_atom_from_string("fp_dst_condition"), &condition);

    op_msg* msg = (op_msg*)vevent;
    ADIOS_FILE *adiosfile = (ADIOS_FILE*)client_data;
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    fp_reader_log(
	"HANDLER",
	"\treplica_id: %d rank: %d got op_msg type: %d with step: %d from writer_replica: %d writer_rank: %d\n",
	fp->replica_id, fp->rank, msg->type, msg->step, msg->replica_id, msg->process_id);
    if (msg->type==ACK_MSG) {
	if (condition != -1) {
	    /* printf("\t\t\t\treplica_id: %d rank: %d op_msg_handler type: %d actually signaling: %d\n", */
	    /* 	   fp->replica_id, fp->rank, msg->type, condition); */
	    CMCondition_signal(local->cm, condition);
	}
        //ackCondition = CMCondition_get(local->cm, NULL);
    }
    if (msg->type == EOS_MSG) {
	adios_errno = err_end_of_stream;
	CMCondition_signal(local->cm, condition);
    }
    fp_reader_log(
	"HANDLER",
	"\treplica_id: %d rank: %d leaving op_msg_handler with step: %dfor type: %d from writer_replica: %d writer_rank: %d\n",
	fp->replica_id, fp->rank, msg->step, msg->type, msg->replica_id, msg->process_id);
    return 0;
}

static int
evgroup_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    EVtake_event_buffer(local->cm, vevent);
    ADIOS_FILE *adiosfile = client_data;
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    evgroup *msg = (evgroup*)vevent;
    evgroup *gp = msg;

    fp_reader_log("HANDLER", 
		  "\treplica_id: %d rank: %d got evgp_msg from writer_replica: %d writer_rank: %d for step: %d\n",
		  fp->replica_id, fp->rank, msg->replica_id, msg->process_id, msg->step);
    int condition = -1;
    get_int_attr(attrs, attr_atom_from_string("fp_dst_condition"), &condition);

    fp->gp = msg;
    int i;
    for (i = 0; i<msg->num_vars; i++) {
	global_var *gblvar = &msg->vars[i];
	Flexpath_var *fpvar = find_fp_var(fp->var_list, gblvar->name);
	if (fpvar) {
	    offset_struct *offset = &gblvar->offsets[0];
	    uint64_t *local_dimensions = offset->local_dimensions;
	    uint64_t *local_offsets = offset->local_offsets;
	    uint64_t *global_dimensions = offset->global_dimensions;

	    fpvar->num_dims = offset->offsets_per_rank;
	    fpvar->global_dims = malloc(sizeof(uint64_t)*fpvar->num_dims);
	    memcpy(fpvar->global_dims, global_dimensions, sizeof(uint64_t)*fpvar->num_dims);
	} else {
	    adios_error(err_corrupted_variable,
			"Mismatch between global variables and variables specified %s.",
			gblvar->name);
	    return err_corrupted_variable;
	}
    }
    CMCondition_signal(local->cm, condition);
    fp_reader_log("HANDLER",
		  "\treplica_id: %d rank: %d leaving evgp_handler from writer_replica: %d writer_rank: %d for step: %d\n",
		  fp->replica_id, fp->rank, msg->replica_id, msg->process_id, msg->step);
    return 0;
}


static int
raw_handler(CManager cm, void *vevent, int len, void *client_data, attr_list attrs)
{
    ADIOS_FILE *adiosfile = client_data;
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    //printf("\t\t\trank: %d raw_handler\n", fp->rank);
    double data_end = dgettimeofday();
    if (fp->time_in == 0.00)
	fp->time_in = data_end; // used for perf measurements only

    int condition;
    int writer_rank;
    int flush_id;
    double data_start;
    get_double_attr(attrs, attr_atom_from_string("fp_starttime"), &data_start);
    get_int_attr(attrs, attr_atom_from_string("fp_dst_condition"), &condition);
    get_int_attr(attrs, attr_atom_from_string(FP_SENDER_RANK), &writer_rank);
    get_int_attr(attrs, attr_atom_from_string("fp_flush_id"), &flush_id);

    fp_reader_log("HANDLER",
		  "\treplica_id: %d rank: %d entering raw_handler from writer_replica: %d writer_rank: %d\n",
		  fp->replica_id, fp->rank, fp->current->replica_id, writer_rank);

    double format_start = dgettimeofday();

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
	    adiosfile->var_namelist[i++] = strdup(f->field_name);
	    f++;
	}
	adiosfile->nvars = var_count;
	fp->num_vars = var_count;
    }

    f = struct_list[0].field_list;
    char *curr_offset = NULL;

    while (f->field_name) {
        char atom_name[200] = "";
    	Flexpath_var *var = find_fp_var(fp->var_list, f->field_name);

    	if (!var) {
    	    adios_error(err_file_open_error,
    			"file not opened correctly.  var does not match format.\n");
    	    return err_file_open_error;
    	}

	int num_dims = get_ndims_attr(f->field_name, attrs);
    	var->num_dims = num_dims;

	Flexpath_var_chunk *curr_chunk = &var->chunks[0];

	// has the var been scheduled?
	if (var->sel) {
	    if (var->sel->type == ADIOS_SELECTION_WRITEBLOCK) {
		if (num_dims == 0) { // writeblock selection for scalar
		    if (var->sel->u.block.index == writer_rank) {
			void *tmp_data = get_FMfieldAddr_by_name(f, f->field_name, base_data);
			memcpy(var->chunks[0].user_buf, tmp_data, f->field_size);
		    }
		}
		else { // writeblock selection for arrays
		    /* if (var->num_dims == 0) { */
		    /* 	var->global_dims = malloc(sizeof(uint64_t)*num_dims); */
		    /* } */
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
			    }
			    else{
				int *temp_data = get_FMfieldAddr_by_name(temp_field,
									 temp_field->field_name,
									 base_data);

				uint64_t dim = (uint64_t)(*temp_data);
				var->array_size = var->array_size * dim;
			    }
			}
			void *arrays_data  = get_FMPtrField_by_name(f, f->field_name, base_data, 1);
			memcpy(var->chunks[0].user_buf, arrays_data, var->array_size);
		    }
		}
	    }
	    else if (var->sel->type == ADIOS_SELECTION_BOUNDINGBOX) {
		if (num_dims == 0) { // scalars; throw error
		    adios_error(err_offset_required,
				"Only scalars can be scheduled with write_block selection.\n");
		}
		else{ // arrays
		    int i;
		    global_var *gv = find_gbl_var(fp->gp->vars,
						  var->varname,
						  fp->gp->num_vars);
		    array_displacements * disp = find_displacement(var->displ,
								   writer_rank,
								   var->num_displ);
		    if (disp) { // does this writer hold a chunk we've asked for, for this var?
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


			extract_selection_from_partial(f->field_size,
                                                       disp->ndims,
                                                       global_dimensions,
                                                       writer_offsets, writer_sizes,
                                                       global_sel_start,
                                                       global_sel_count,
                                                       writer_array,
                                                       reader_array);
		    }
		}
	    }
	}
	else { //var has not been scheduled;
	    if (num_dims == 0) { // only worry about scalars
		Flexpath_var_chunk *chunk = &var->chunks[0];
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

    if (condition == -1) {
	fp->completed_requests++;
	if (fp->completed_requests == fp->pending_requests) {
	    pthread_mutex_lock(&fp->data_mutex);
	    pthread_cond_signal(&fp->data_condition);
	    pthread_mutex_unlock(&fp->data_mutex);
	}
    }
    else{
	CMCondition_signal(local->cm, condition);
    }
        
    free_fmstructdesclist(struct_list);
    fp_reader_log("HANDLER",
		  "\treplica_id: %d rank: %d leaving raw_handler from writer_replica: %d writer_rank: %d\n",
		  fp->replica_id, fp->rank, fp->current->replica_id, writer_rank);
    //printf("rank: %d leaving raw_handler\n", fp->rank);
    return 0;
}

static int
replica_finalized(Flexpath_reader_replica *rep)
{
    if (!rep->writer_finalized)
        return 0;
    if (rep->count > 0)
        return 0;
    return 1;
}

static int
all_finalized(Flexpath_reader_replica *list, pthread_mutex_t *lock)
{
    pthread_mutex_lock(lock);
    while (list) {
        if (!replica_finalized(list)) {
            pthread_mutex_unlock(lock);
            return 0;
        }
        list = list->next;                        
    }
    return 1;
}


static void
next_replica_rr(Flexpath_reader_file *fp)
{
  pthread_mutex_lock(&fp->repmutex);
  if (!fp->current->next) {
    fp->current = fp->replicas;
  }
  else {
    fp->current = fp->current->next;
  }
  pthread_mutex_unlock(&fp->repmutex);
}

static void
add_replica_tolist(Flexpath_reader_replica **list, Flexpath_reader_replica *rep)
{
    if (!(*list)) {
	*list = rep;
    } else {
	Flexpath_reader_replica *tmp = *list;
	while (tmp && tmp->next) {
	    tmp = tmp->next;
	}
	tmp->next = rep;
    }
}

static void
process_update_contact(Flexpath_reader_file *fp, Update_contact_msg *msg)
{
    int i;
    for (i = 0; i<msg->num_replicas; i++) {
	Replica_info_msg *rep = &msg->replicas[i];
	int num_writers = rep->num_writers;
        
	Flexpath_reader_replica *fprep = calloc(1, sizeof(Flexpath_reader_replica));
        fprep->last_writer_step = -1;
        fprep->active = 1;         
	fprep->mystep = 0;
	//fprep->container_name = strdup(msg->container_name);
	fprep->replica_id = rep->replica_id;
	fprep->num_bridges = num_writers;
	fprep->bridges = malloc(sizeof(Bridge_info) * fprep->num_bridges);
	int j;
	for (j = 0; j<num_writers; j++) {
	    char in_contact[CONTACT_LENGTH];
	    int their_stone;
	    //printf("%s ", rep->writers[j].endpoint);
	    sscanf(rep->writers[j].endpoint, "%d:%s", &their_stone, in_contact);
	    //fprintf(stderr, "writer contact: %d:%s\n", their_stone, in_contact);
	    fprep->bridges[j].their_num = their_stone;
	    fprep->bridges[j].contact = strdup(&in_contact[0]);
	    fprep->bridges[j].created = 0;
	    fprep->bridges[j].step = 0;
	    fprep->bridges[j].opened = 0;
	    fprep->bridges[j].scheduled = 0;
	}

	if (fp->size < rep->num_writers) {
	    int mystart = (num_writers/fp->size) * fp->rank;
	    int myend = (num_writers/fp->size) * (fp->rank+1);
	    fprep->writer_coordinator = mystart;
	    int z;
	    for (z=mystart; z<myend; z++) {
		build_bridge(&fprep->bridges[z]);
	    }
	}
	else {
	    int writer_rank = fp->rank % num_writers;
	    build_bridge(&fprep->bridges[writer_rank]);
	    fprep->writer_coordinator = writer_rank;
	}

	pthread_mutex_lock(&fp->repmutex);
	add_replica_tolist(&fp->replicas, fprep);
	pthread_mutex_unlock(&fp->repmutex);
    }
}
//#endif
/********** Core ADIOS Read functions. **********/

/*
 * Gathers basic MPI information; sets up reader CM.
 */
int
adios_read_flexpath_init_method (MPI_Comm comm, PairStruct* params)
{
    fprintf(stderr, "\t\tflexpath init called.\n");
    //fp_reader_log("FUNC", "reader init.\n");
    setenv("CMSelfFormats", "1", 1);
    local = malloc(sizeof(Flexpath_reader_local));
    fprintf(stderr, "\t\tflexpath init after malloc.\n");
    if (!local) {
        adios_error(err_no_memory, "Cannot allocate memory for flexpath.");
        return -1;
    }
    memset(local, 0, sizeof(Flexpath_reader_local));
    fprintf(stderr, "\t\tflexpath init after memset.\n");
    local->CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    attr_list listen_list = NULL;
    char * transport = NULL;
    transport = getenv("CMTransport");
    fprintf(stderr, "\t\tflexpath init after getenv.\n");
    // setup MPI stuffs
#ifndef _NOMPI
    fprintf(stderr, "\t\tflexpath init this should NOT be called..\n");
    local->comm = comm;
    MPI_Comm_size(local->comm, &(local->size));
    MPI_Comm_rank(local->comm, &(local->rank));
#else
    int set = 0;
    fprintf(stderr, "\t\tflexpath init before get_relay_client.\n");
    local->client = adios_read_get_relay_client();
    fprintf(stderr, "\t\tflexpath init after get_relay_client.\n");
    if (!local->client) {
	fprintf(stderr, "\t\tflexpath init client is null. problem.\n");
    }
    local->size = MPIRelay_client_size(local->client);
    local->rank = MPIRelay_client_rank(local->client);
#endif
    
    local->cm = CManager_create();
    if (transport == NULL) {
	while (CMlisten(local->cm) == 0) {
	    fprintf(stderr,
		    "Flexpath ERROR: reader %d unable to initconnection manager.\n",
		    local->rank);
	}
    }else{
	listen_list = create_attr_list();
	add_attr(listen_list, local->CM_TRANSPORT, Attr_String,
		 (attr_value)strdup(transport));
	CMlisten_specific(local->cm, listen_list);
    }
    int forked = CMfork_comm_thread(local->cm);
    if (!forked) {
	fprintf(stderr, "reader %d failed to fork comm_thread.\n", local->rank);
	/*log_debug( "forked\n");*/
    }

    local->stone = EValloc_stone(local->cm);

    int name_set;
    int id_set;
    char *cname = globals_adios_get_container_name(&name_set);
    int replica_id = globals_adios_get_application_id(&id_set);

    fprintf(stderr, "\t\t NAME: %s ID: %d\n", cname, replica_id);
    
    char *string_list;
    char data_contact_info[CONTACT_LENGTH];
    if (!local->tunnel) {
	string_list = attr_list_to_string(CMget_contact_list(local->cm));
	sprintf(&data_contact_info[0], "%d:%s", local->stone, string_list);
	free(string_list);
    } else {
	attr_list string_attr = create_attr_list();
	add_int_attr(string_attr, attr_atom_from_string("IP_PORT"), local->port-1000);
	add_string_attr(string_attr, attr_atom_from_string("IP_HOST"), "localhost");
	add_int_attr(string_attr, attr_atom_from_string("IP_ADDR"), 0x7f000001);
	char *tunnel_string = attr_list_to_string(string_attr);
	sprintf(&data_contact_info[0], "%d:%s", local->stone, tunnel_string);
    }
    
    local->ctx = replica_init(local->cm,
			      local->size,
			      local->rank,
			      replica_id,
			      &data_contact_info[0],
			      cname,
			      READER_REPLICA);

    // JAI: FIX THIS!@!!!!
#ifndef _NOMPI    
    char *allcontacts = gather_contacts(local->comm, &data_contact_info[0], 0, local->rank);
    char *allhosts = gather_hostnames(local->comm, 0, local->rank);
#else
    char *allcontacts = gather_contacts2(local->client, &data_contact_info[0], 0, local->rank);
    char *allhosts = gather_hostnames2(local->client, 0, local->rank);
#endif
    
    if (local->rank == 0) {
	report_new_replica(local->ctx, allcontacts, allhosts);
    }

    //fp_reader_log("FUNC", "reader init leaving\n");
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
    ADIOS_FILE *adiosfile = malloc(sizeof(ADIOS_FILE));
    if (!adiosfile) {
	adios_error (err_no_memory,
		     "Cannot allocate memory for file info.\n");
	return NULL;
    }

    Flexpath_reader_file *fp = new_Flexpath_reader_file(fname);
    fp->host_language = futils_is_called_from_fortran();
    adios_errno = 0;
    fp->stone = local->stone;
    
#ifndef _NOMPI    
    fp->comm = comm;
    MPI_Comm_size(fp->comm, &(fp->size));
    MPI_Comm_rank(fp->comm, &(fp->rank));
#else
    fp->client = local->client;
    fp->size = local->size;
    fp->rank = local->rank;
#endif    
    //printf("rank: %d opening.\n", fp->rank);

    fp->replica_id = local->ctx->replica_id;
    fp_reader_log("FUNC", "=====replica_id: %d, rank: %d entering reader flexpath_open.=====\n",
		  fp->replica_id, fp->rank);

    EVassoc_terminal_action(local->cm,
			    fp->stone,
			    op_format_list,
			    op_msg_handler,
			    adiosfile);

    EVassoc_terminal_action(local->cm,
			    fp->stone,
			    update_step_msg_format_list,
			    update_step_msg_handler,
			    adiosfile);


    EVassoc_terminal_action(local->cm,
			    fp->stone,
			    evgroup_format_list,
			    evgroup_handler,
			    adiosfile);

    EVassoc_raw_terminal_action(local->cm,
				fp->stone,
				raw_handler,
				adiosfile);

    /* Gather the contact info from the other readers
       and write it to a file. Create a ready file so
       that the writer knows it can parse this file. */
    Replica_context *ctx = local->ctx;

    int optype = REPLICA_NO_OP;
    void *op = NULL;

    if (fp->rank == 0) {
	while (!op) {
	    //printf("why here?\n");
	    op = replica_get_operation(ctx, &optype);
	    CMsleep(local->cm, 1);
	}
    }
    
#ifndef _NOMPI
    MPI_Bcast(&optype, 1, MPI_INT, 0, fp->comm);
#else
    MPIRelay_bcast(fp->client, &optype, 1, MPIInt, 0);
#endif
    
    //printf("got here.\n");
    if (optype == REPLICA_ADD_CONTACTS) {
	Update_contact_msg *upstream = exchange_update_msg(fp, op);
	process_update_contact(fp, upstream);
    }
    else {
	// not valid control message at this point. what to do?
    }
    fp->current = fp->replicas;

#ifndef _NOMPI
    MPI_Barrier(fp->comm);
#else
    MPIRelay_barrier(fp->client);
#endif
    
    adiosfile->fh = (uint64_t)fp;
    adiosfile->current_step = 0;

    send_flush_msg(fp, fp->current, fp->current->writer_coordinator, STEP, 1);
    while (fp->current->count == 0) {	
	if (all_finalized(fp->replicas, &fp->repmutex)) {
	    adios_errno = err_end_of_stream;
	    return NULL;
	}
	next_replica_rr(fp);
	send_flush_msg(fp, fp->current, fp->current->writer_coordinator, STEP, 1);
	CMsleep(local->cm, 1);
    }


    /* Init with a writer to get initial scalar
       data so we can handle inq_var calls and
       also populate the ADIOS_FILE struct. */

    // requesting initial data.
    send_open_msg(fp, fp->current, fp->current->writer_coordinator);

    fp->data_read = 0;
    send_flush_msg(fp, fp->current, fp->current->writer_coordinator, DATA, 1);
    send_flush_msg(fp, fp->current, fp->current->writer_coordinator, EVGROUP, 1);
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

    fp_reader_log("FUNC", "=====replica_id: %d rank: %d leaving flexpath_open.=====\n",
		  fp->replica_id, fp->rank);
    return adiosfile;
}

int
adios_read_flexpath_finalize_method ()
{
    return 0;
}

void
adios_read_flexpath_release_step(ADIOS_FILE *adiosfile)
{
    int i;
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;

    fp_reader_log("FUNC",
                  "=====replica_id :%d rank: %d entering release_step=====\n",
                  fp->replica_id, fp->rank);
    /*TODO: optimize this. The whole open/close/step stuff needs to be optimized.*/
    for (i=0; i<fp->current->num_bridges; i++) {
        if ( (!fp->current->bridges[i].opened) && (fp->current->bridges[i].created)) {
            send_open_msg(fp, fp->current, i);
        }        
    }
    
#ifndef _NOMPI
    MPI_Barrier(fp->comm);
#else
    MPIRelay_barrier(fp->client);    
#endif
    
    for (i=0; i<fp->current->num_bridges; i++) {
	if (fp->current->bridges[i].opened) {
	    send_close_msg(fp, fp->current, i);
	}
    }
    fp->current->count--;
    free_evgroup(fp->gp);
    fp->gp = NULL;

    Flexpath_var *tmpvars = fp->var_list;
    while (tmpvars) {
	if (tmpvars->num_dims > 0) {
	    free(tmpvars->global_dims);
	    tmpvars->num_dims = 0;
	}
	free_displacements(tmpvars->displ, tmpvars->num_displ);
	tmpvars->displ = NULL;
	if (tmpvars->sel) {
	    free_selection(tmpvars->sel);
	}
	tmpvars->sel = NULL;
	for (i=0; i<tmpvars->num_chunks; i++) {
	    Flexpath_var_chunk *chunk = &tmpvars->chunks[i];
	    if (chunk->has_data) {
		free(chunk->data);
		chunk->data = NULL;
		chunk->has_data = 0;
	    }
	    free(chunk->local_bounds);
	    free(chunk->global_offsets);
	    free(chunk->global_bounds);

	    chunk->global_offsets = NULL;
	    chunk->global_bounds = NULL;
	    chunk->local_bounds = NULL;
	    chunk->rank = 0;
	}
	tmpvars = tmpvars->next;
    }

#ifndef _NOMPI
    MPI_Barrier(fp->comm);
#else
    MPIRelay_barrier(fp->client);    
#endif
    fp_reader_log("FUNC",
                  "=====replica_id: %d rank: %d leaving release_step=====\n",
                  fp->replica_id, fp->rank);
}

int
adios_read_flexpath_advance_step(ADIOS_FILE *adiosfile, int last, float timeout_sec)
{
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    fp_reader_log("FUNC", "=====replica_id: %d rank: %d entering advance_step.=====\n",
	   fp->replica_id, fp->rank);

    fp->current->mystep++;
    /* if (all_finalized(fp->replicas, &fp->repmutex)) { */
    /* 	adios_errno = err_end_of_stream; */
    /* 	return err_end_of_stream; */
    /* } */


    int optype = REPLICA_NO_OP;
    void *op = NULL;

    if (fp->rank == 0) {
	op = replica_get_operation(local->ctx, &optype);
    }
#ifndef _NOMPI
    MPI_Bcast(&optype, 1, MPI_INT, 0, fp->comm);
#else
    MPIRelay_bcast(fp->client, &optype, 1, MPIInt, 0);
#endif
    
    if (optype == REPLICA_ADD_CONTACTS) {
	Update_contact_msg *upstream = exchange_update_msg(fp, (Update_contact_msg*)op);
	process_update_contact(fp, upstream);
#ifndef _NOMPI
	MPI_Barrier(fp->comm);
#else
	MPIRelay_barrier(fp->client);
#endif
	
    }
    else if (optype == REPLICA_REMOVE_CONTACTS) {

    }
    else if (optype == REPLICA_PAUSE) {

    }
    else {
	// not valid control message at this point. what to do?
    }

    next_replica_rr(fp); // HAVE COUNT PIGGY_BACK CLOSE ACK TO AVOID THIS PROBLEM
    if (!replica_finalized(fp->current)) {
        send_flush_msg(fp, fp->current, fp->current->writer_coordinator, STEP, 1);
    }

    //send_flush_msg(fp, fp->current, fp->current->writer_coordinator, STEP, 1);
    while (fp->current->count == 0) {
	if (all_finalized(fp->replicas, &fp->repmutex)) {
	    adios_errno = err_end_of_stream;
	    return err_end_of_stream;
	}
	next_replica_rr(fp);
        if (!replica_finalized(fp->current)) {
            send_flush_msg(fp, fp->current, fp->current->writer_coordinator, STEP, 1);
        }
	CMsleep(local->cm, 1);
    }
#ifndef _NOMPI
    MPI_Barrier(fp->comm);
#else
    MPIRelay_barrier(fp->client);
#endif
    /* if (all_finalized(fp->replicas, &fp->repmutex)) { */
    /* 	adios_errno = err_end_of_stream; */
    /* 	return err_end_of_stream; */
    /* } */
    /* printf("replica_id: %d rank: %d global_step: %d to writer: %d coord: %d\n", */
    /* 	   fp->replica_id, fp->rank, global_step, fp->current->replica_id, */
    /* 	   fp->current->writer_coordinator); */
    send_open_msg(fp, fp->current, fp->current->writer_coordinator);

    send_flush_msg(fp, fp->current, fp->current->writer_coordinator, DATA, 1);
    send_flush_msg(fp, fp->current, fp->current->writer_coordinator, EVGROUP, 1);

#ifndef _NOMPI
    MPI_Barrier(fp->comm);
#else
    MPIRelay_barrier(fp->client);
#endif
    
    global_step++;
    //printf("replica_id: %d rank: %d leaving advance_step\n", fp->replica_id, fp->rank);
    fp_reader_log("FUNC", "=====replica_id: %d rank: %d leaving advance_step.=====\n",
	   fp->replica_id, fp->rank);
    return 0;
}

int adios_read_flexpath_close(ADIOS_FILE *file)
{
    Flexpath_reader_file *fp = (Flexpath_reader_file*)file->fh;
    fp_reader_log("FUNC", "=====replica_id: %d rank: %d entering close.=====\n", fp->replica_id, fp->rank);
    //send to each opened link
    // has to be repeated for each replica
    int i;
    for (i = 0; i<fp->current->num_bridges; i++) {
        if (fp->current->bridges[i].created && fp->current->bridges[i].opened) {
	    send_close_msg(fp, fp->current, i);
        }
    }
    /*
    start to cleanup.  Clean up var_lists for now, as the
    data has already been copied over to ADIOS_VARINFO structs
    that the user maintains a copy of.
    */
    Flexpath_var *v = fp->var_list;
    while (v) {
    	// free chunks; data has already been copied to user
    	int i;
    	for (i = 0; i<v->num_chunks; i++) {
    	    Flexpath_var_chunk *c = &v->chunks[i];
	    if (!c)
		log_error("FLEXPATH: %s This should not happen! line %d\n",__func__,__LINE__);
	    //free(c->data);
	    free(c->global_bounds);
	    free(c->global_offsets);
	    free(c->local_bounds);
	    c->data = NULL;
	    c->global_bounds = NULL;
	    c->global_offsets = NULL;
	    c->local_bounds = NULL;
	    free(c);
	}
	Flexpath_var *tmp = v->next;
	free(v);
	v = tmp;
    	//v=v->next;
    }
    CManager_close(local->cm);
    fp_reader_log("FUNC", "=====replica_id: %d rank: %d leaving close.=====\n", fp->replica_id, fp->rank);
    return 0;
}

ADIOS_FILE *adios_read_flexpath_fopen(const char *fname, MPI_Comm comm) {
   return 0;
}

int adios_read_flexpath_is_var_timed(const ADIOS_FILE* fp, int varid) { return 0; }

void adios_read_flexpath_get_groupinfo(const ADIOS_FILE *fp, int *ngroups, char ***group_namelist, uint32_t **nvars_per_group, uint32_t **nattrs_per_group) {}

int adios_read_flexpath_check_reads(const ADIOS_FILE* fp, ADIOS_VARCHUNK** chunk) { log_debug( "flexpath:adios function check reads\n"); return 0; }

int adios_read_flexpath_perform_reads(const ADIOS_FILE *adiosfile, int blocking)
{
    Flexpath_reader_file * fp = (Flexpath_reader_file*)adiosfile->fh;
    fp_reader_log("FUNC", "=====replica_id: %d rank: %d entering perform_reads.=====\n", fp->replica_id, fp->rank);
    fp->data_read = 0;
    int i,j;
    int num_sendees = fp->num_sendees;
    int total_sent = 0;
    fp->time_in = 0.00;

    //printf("before send flush perform rank: %d, mystep: %d\n", fp->rank, fp->mystep);
    for (i = 0; i<num_sendees; i++) {
	pthread_mutex_lock(&fp->data_mutex);
	int sendee = fp->sendees[i];
	fp->pending_requests++;
	total_sent++;

	send_flush_msg(fp, fp->current, sendee, DATA, 0);
	if ((total_sent % FP_BATCH_SIZE == 0) || (total_sent = num_sendees)) {
	    pthread_cond_wait(&fp->data_condition, &fp->data_mutex);
	    pthread_mutex_unlock(&fp->data_mutex);
	    fp->completed_requests = 0;
	    fp->pending_requests = 0;
	    total_sent = 0;
	}

    }
    //printf("after send flush perform rank: %d, mystep: %d num: %d\n", fp->rank, fp->mystep, num_sendees);

    free(fp->sendees);
    fp->sendees = NULL;
    fp->num_sendees = 0;

    fp_reader_log("FUNC", "=====replica_id: %d rank: %d leaving perform_reads.=====\n", fp->replica_id, fp->rank);
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
    Flexpath_reader_file * fp = (Flexpath_reader_file*)adiosfile->fh;
    fp_reader_log("FUNC", "=====replica_id: %d, rank: %d entering schedule_read_byid.=====\n",
		  fp->replica_id, fp->rank);
    //printf("replica_id: %d rank: %d entering schedule read\n", fp->replica_id, fp->rank);
    Flexpath_var *var = fp->var_list;
    while (var) {
        if (var->id == varid)
        	break;
        else
	    var=var->next;
    }
    if (!var) {
        adios_error(err_invalid_varid,
		    "Invalid variable id: %d\n",
		    varid);
        return err_invalid_varid;
    }

    //store the user allocated buffer.
    Flexpath_var_chunk *chunk = &var->chunks[0];
    chunk->user_buf = data;
    var->start_position = 0;
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
    if (var->sel) {
	free_selection(var->sel);
    }
    var->sel = copy_selection(sel);

    switch(var->sel->type)
    {
    case ADIOS_SELECTION_WRITEBLOCK:
    {
	int writer_index = var->sel->u.block.index;
	if (writer_index > fp->current->num_bridges) {
	    adios_error(err_out_of_bound,
			"No process exists on the writer side matching the index.\n");
	    return err_out_of_bound;
	}
	send_var_message(fp, fp->current, writer_index, var->varname);
	break;
    }
    case ADIOS_SELECTION_BOUNDINGBOX:
    {
        if (var->num_dims == 0) {
            memcpy(data, chunk->data, var->type_size);
        } else {
            if (fp->host_language == FP_FORTRAN_MODE) {
                reverse_dims(sel->u.bb.start, sel->u.bb.ndim);
                reverse_dims(sel->u.bb.count, sel->u.bb.ndim);
            }
            free_displacements(var->displ, var->num_displ);
            var->displ = NULL;
            int j=0;
            int need_count = 0;
            array_displacements * all_disp = NULL;
            uint64_t pos = 0;
            double sched_start = MPI_Wtime();
            for (j=0; j<fp->current->num_bridges; j++) {
                int destination=0;
                if (need_writer(fp, j, var->sel, fp->gp, var->varname)==1) {
                    uint64_t _pos = 0;
                    need_count++;
                    destination = j;
                    global_var *gvar = find_gbl_var(fp->gp->vars, var->varname, fp->gp->num_vars);
                    // TODO: memory leak here. have to free these at some point.
                    array_displacements *displ = get_writer_displacements(j, var->sel, gvar, &_pos);
                    displ->pos = pos;
                    _pos *= (uint64_t)var->type_size;
                    pos += _pos;

                    all_disp = realloc(all_disp, sizeof(array_displacements)*need_count);
                    all_disp[need_count-1] = *displ;
                    send_var_message(fp, fp->current, j, var->varname);
                }
            }
            double sched_end = MPI_Wtime();
            var->displ = all_disp;
            var->num_displ = need_count;
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
    fp_reader_log("FUNC", "=====replica_id: %d, rank: %d leaving schedule_read_byid.=====\n",
		  fp->replica_id, fp->rank);
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
    //fp_reader_log("FUNC", "entering flexpath_inq_var\n");
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    ADIOS_VARINFO *v = malloc(sizeof(ADIOS_VARINFO));

    if (!v) {
        adios_error(err_no_memory,
		    "Cannot allocate buffer in adios_read_datatap_inq_var()");
        return NULL;
    }
    memset(v, 0, sizeof(ADIOS_VARINFO));

    Flexpath_var *fpvar = find_fp_var(fp->var_list, varname);
    if (fpvar) {
	v = convert_var_info(fpvar, v, varname, adiosfile);
	//fp_reader_log("FUNC", "leaving flexpath_inq_var\n");
	return v;
    }
    else {
        adios_error(err_invalid_varname, "Cannot find var %s\n", varname);
        return NULL;
    }
}

ADIOS_VARINFO*
adios_read_flexpath_inq_var_byid (const ADIOS_FILE * adiosfile, int varid)
{
    //fp_reader_log("FUNC", "entering flexpath_inq_var_byid\n");
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    if (varid >= 0 && varid < adiosfile->nvars) {
	ADIOS_VARINFO *v = adios_read_flexpath_inq_var(adiosfile, adiosfile->var_namelist[varid]);
	//fp_reader_log("FUNC", "leaving flexpath_inq_var_byid\n");
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
    adios_error(err_operation_not_supported,
		"Flexpath does not yet support transforms: trans_blockinfo.\n");
    return (int64_t)0;
}

int  adios_read_flexpath_get_dimension_order (const ADIOS_FILE *adiosfile)
{
    Flexpath_reader_file *fp = (Flexpath_reader_file*)adiosfile->fh;
    return (fp->host_language == FP_FORTRAN_MODE);
}

void
adios_read_flexpath_reset_dimension_order (const ADIOS_FILE *adiosfile, int is_fortran)
{
    //log_debug( "debug: adios_read_flexpath_reset_dimension_order\n");
    adios_error(err_invalid_read_method,
		"adios_read_flexpath_reset_dimension_order is not implemented.");
}
