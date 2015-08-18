#include "config.h"

#if NO_DATATAP == 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ffs.h>
#include <atl.h>
#include <evpath.h>
#include "public/adios_mpi.h"

#include <pthread.h>
#include "adios.h"
#include "adios_read.h"
#include "adios_read_hooks.h"
#include "adios_error.h"
#include "globals.h"

#include <sys/queue.h>
#if HAVE_PORTALS == 1
#include <thin_portal.h>
#elif HAVE_INFINIBAND == 1
#include <thin_ib.h>
#endif

#include <sys/socket.h>
#include <sys/times.h>
#include <netinet/in.h>
#include <sys/time.h>
#include <sys/uio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>

#include <gen_thread.h>
//#include "queue.h"
//#include "get_clock.h"
//#include "attributes.h"
//#include "memwatch.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif


#define DT_MAX_QUEUE_LENGTH 512 

typedef struct _Queue_Item
{
    uint32_t length;
    char * data;
    FMStructDescList var_list;
    int32_t rank;     // rank of client from which the chunk comes
} QueueItem;

typedef struct _datatap_var_chunk
{
    int rank;
    void *data;
    uint64_t *local_bounds; // ndims
    uint64_t *global_bounds; // ndims
    uint64_t *global_offsets; // ndims
    struct _datata_var_chunk *next;
} datatap_var_chunk, *datatap_var_chunk_p;

typedef struct _datatap_var_info
{
    int id;
    char *varname;
    char *varpath;
    enum ADIOS_DATATYPES type;
    uint64_t data_size;   
    int time_dim; // -1 means no time dimension
    int ndims;
    int num_chunks;
    datatap_var_chunk *chunks;
    struct _datatap_var_info *next;
} datatap_var_info, *datatap_var_info_p;


//typedef struct _datatap_pg_info
//{
//    int rank;    
//    int num_vars;    
//    datatap_var_info *vars;    
//} datatap_pg_info, *datatap_pg_info_p;

#define VAR_BITMAP_SIZE 16

// TODO
typedef struct _datatap_read_file_data
{
    char *file_name;
    char *group_name; // TODO: assume one group in file
    file_info *f_info;
//    int timestep;    // TODO: it's already in file_info
    int num_vars;
    datatap_var_info *vars;
//    FMStructDescList var_list; // TODO

    // TODO: replicated meta-data info from peer readers
    // it is a list of vars which is not present locally
    int num_vars_peer;
    datatap_var_info *vars_peer;

    MPI_Comm comm;
    int my_rank;
    int comm_size;

    int host_language_fortran; // 1 for Fortran; 0 for C

    char var_bitmap[VAR_BITMAP_SIZE]; // 128 bit for 128 var

//    int num_vars_read;    
//    datatap_var_info *vars_read;    
//    datatap_var_info *vars_read_tail;    
} datatap_read_file_data, *datatap_read_file_data_p;

typedef struct _datatap_read_method_data
{
    pthread_t dt_server_thread;
//    Queue *dt_queue;
//    uint32_t dt_queue_max_length;
//    pthread_mutex_t mutex;
//    pthread_cond_t cond1;
//    pthread_cond_t cond2;
    MPI_Comm dt_comm;
    int dt_comm_rank;
    int dt_comm_size;
    CManager dt_cm;
    IOhandle *dt_handle;
    int dt_server_ready; 
//    int dt_server_stop; 
    int num_io_dumps;    // TODO: timestep
} datatap_read_method_data, *datatap_read_method_data_p;


// this sructure holds all global data for Datatap Read method
datatap_read_method_data *dt_read_data = NULL;

// compare used-providd varname with the full path name of variable v
// return zero if matches and non-zero otherwise
static int compare_var_name (char *varname, datatap_var_info *v) 
{
    if (varname[0] == '/') { // varname is full path
        char fullpath[256];
        if(!strcmp(v->varpath, "/")) {
            sprintf(fullpath, "/%s\0", v->varname);  
        } 
        else {   
            sprintf(fullpath, "%s/%s\0", v->varpath, v->varname);  
        }
        return strcmp(fullpath, varname);
    }
    else { // varname doesn't include path
        return strcmp(v->varname, varname);
    }
}

static FMField *find_field (char *name, char *path, FMFieldList flist)
{
    char *temp_name;
    char *full_path_name = get_full_path_name(name, path);
    temp_name = getFixedName(full_path_name);
    free(full_path_name);
    FMField *f = flist;
    while (f->field_name != NULL) {
        if(!strcmp(temp_name, f->field_name)) {
            free(temp_name);
            return f;
        }
        else {
            f++;
        }
#if 0
        char *name_pos = f->field_name + (strlen(f->field_name) - strlen(name));
        if(!strncmp(path, f->field_name, strlen(path)) &&
           !strcmp(name, name_pos)) {
            return f;
        }
        else {
            f++;
        }
#endif
    }
    free(temp_name);
    return f;
}

int64_t read_array(datatap_read_file_data *ds, datatap_var_info *var, 
                    uint64_t *start, uint64_t *count, void *data)
{
    int type_size = common_read_type_size(var->type, NULL);
fprintf(stderr, "im here type_size %d %s:%d\n", type_size, __FILE__,__LINE__);
    int64_t total_size = 0;

    // go over the whole list of data chunks and reorganize arrays
    Queue *data_q = ds->f_info->dt_queue;

    // TODO: we should be able to access the queue since the dt server will not
    // touch it any more
    ListElmt *current_chunk = list_head(data_q);
    while(current_chunk) {
        QueueItem *qi = (QueueItem *)current_chunk->data;  
        FMFieldList filed_list = qi->var_list->field_list; // TODO
        
        // first, find the var   
        FMField *f = find_field(var->varname, var->varpath, filed_list);

        if(!f) { 
            // actually this will not happen because the filed_list will
            // contain all vars even though some may not be written
            current_chunk = current_chunk->next;
            continue;
        }

        char *source_addr = (char *)qi->data + f->field_offset; // this is the address
        void *source_addr2;
  
        // here we need to distinguish static and dynamic arrays
        int dim_are_nums = 1;

        {
            char *dim_start = strchr(f->field_type, '[');
            char *dim_end = strrchr(f->field_type, ']');
            while(dim_start < dim_end) {
                if(*dim_start == '[' || *dim_start == ']') {
                    dim_start ++;
                    continue;
                } else if(!isdigit(*dim_start)) {
                    dim_are_nums = 0;
                    break;
                }
                else {
                    dim_start ++;
                    continue;
                }
            }
        }

        if(dim_are_nums) {
            source_addr2 = source_addr;
        } else {
            switch(var->type) {
                case adios_byte:
                case adios_unsigned_byte:
                    source_addr2 = *((char **)source_addr);
                    break;

                case adios_string: // TODO
                    return strlen((char *)source_addr) + 1;
                    //source_addr2 = *((char **)source_addr);
                    //break;

                case adios_short:
                case adios_unsigned_short:
                    source_addr2 = *((short int **)source_addr);
                    break;

                case adios_integer:
                case adios_unsigned_integer:
                    source_addr2 = *((int **)source_addr);
                    break;

                case adios_long:
                case adios_unsigned_long:
                    source_addr2 = *((long int **)source_addr);
                    break;

                case adios_real:
                    source_addr2 = *((float **)source_addr);
                    break;

                case adios_double:
                    source_addr2 = *((double **)source_addr);
                    break;

                case adios_long_double:
                    source_addr2 = *((long double **)source_addr);
                    break;

                case adios_complex:
                case adios_double_complex:
                default:
                    adios_error(err_invalid_read_method, "complex data type is not supported.");
                    return -1;
            }
        }

        // find the array dimension info 
        int j;
        datatap_var_chunk *chunk = var->chunks; 
        while(chunk) {
            if(chunk->rank == qi->rank) {
                break;
            } 
            else {
                chunk = chunk->next;
            }
        }

        if(!chunk) { // this chunk doesn't contain the var so skip it
            current_chunk = current_chunk->next;
            continue;
        }

        // now we copy data into user's read buffer
        uint64_t global_start_r = 0;
        for(j = var->ndims-1; j >= 0; j --) {
            global_start_r *= chunk->global_bounds[j];
            global_start_r += start[j];
        }


        total_size = type_size;
        for(j = var->ndims-1; j >= 0; j --) {
            total_size *= count[j];
        }

        for(j = 0; j < var->ndims; j ++) {
            int lower = chunk->global_offsets[j] - start[j];
            int higher = (start[j]+count[j]) - (chunk->global_offsets[j]+chunk->local_bounds[j]);
            if(lower < 0 || higher < 0) {
                // this means this chunk has some data which is not needed by this reader
                // TODO: for pixie3d case, this will not happen
                adios_error(err_invalid_read_method, "Datatap cannot support this decomposition.\n");
                return -1;
            }
        }

        // find the smallest slice
        int min_dim = 0;
        uint64_t stride_size = type_size;
        for(j = 0; j < var->ndims; j ++) {
            int lower = chunk->global_offsets[j] - start[j];
            int higher = (start[j]+count[j]) - (chunk->global_offsets[j]+chunk->local_bounds[j]);
            if(lower == 0 && higher == 0) {
                // this means this dimnesion fits, let's move on to the next dimension
                stride_size *= chunk->local_bounds[j];
                min_dim ++; 
                continue;
            }
            else if((lower > 0 && higher >= 0) ||
                    (lower >= 0 && higher > 0)) {
                // this means this dimension is strided. this dimension is the highest possible 
                // dimension
                break;
            }
        }
 
        uint64_t num_strides = 1;
        for(j = min_dim; j < var->ndims; j ++) {
            num_strides *= chunk->local_bounds[j];
        }

        uint64_t current_pos[var->ndims]; // track current stride's starting global offset
        for(j = 0; j < var->ndims; j ++) {
            current_pos[j] = chunk->global_offsets[j];
        }
        uint64_t k = 0;       
        char *stride_start_addr = source_addr2; 
        int should_stop = 0;
        for(; k < num_strides; k ++) {
            // calculate the coordinates within the bounding box
            uint64_t my_start = 0;
            for(j = var->ndims-1; j >= 0; j --) {
                my_start *= count[j];
                my_start += (current_pos[j] - chunk->global_offsets[j]);
            } 
            char *position = (char *)data + my_start * type_size;

            memcpy(position, stride_start_addr, stride_size);

            // now advance to copy next stride
            stride_start_addr += stride_size; 
            should_stop = 0;
            for(j = min_dim; j < var->ndims; j ++) {
                if(!should_stop) {
                    if(current_pos[j] == chunk->global_offsets[j]+chunk->local_bounds[j]-1) {
                        // don't set should_stop so we move on to next dimension
                        current_pos[j] = chunk->global_offsets[j];
                    }
                    else {
                        current_pos[j] ++;
                        should_stop = 1;
                    }
                } 
                else {
                    break;
                }
            }
        }  

#if 0
        for(j = var->ndims-1; j >= 0; j --) {
            // for each sub-chunk, test if it falls into the reading region
            int lower = (chunk->global_offsets[j] >= start[j]) ? 1: 0;
            int higher = ((chunk->global_offsets[j]+chunk->local_bounds[j]) <=
                 (start[j]+count[j])) ? 1: 0;
            if(lower && higher) {
                if(j == 0) { // copy the whole chunks
                    uint64_t my_start = 0;
                    uint64_t my_size = type_size;
                    int k;
                    for(k = var->ndims-1; k >= 0; k --) {
                        my_start *= chunk->global_bounds[k];
                        my_start += chunk->global_offsets[k];
                        my_size *= chunk->local_bounds[k];
                    }
                    
                    int y, z;
                    for(z = 0; z < ; z ++) {
                        uint64_t z_offset = pix_record->zoffset + k;
                        for(y = 0; y < pix_record->ysize; y ++) {
                            unsigned long long y_offset = pix_record->yoffset + j;
                            unsigned long long position = z_offset * pix_record->nyd_plus_2 * pix_record->nxd_plus_2 +
                               y_offset * pix_record->nxd_plus_2 + pix_record->xoffset;
                            position -= p_data->start_pos;
                            memcpy(&(p_data->buffers[m].buffer[position]), v, pix_record->xsize * sizeof(double));
                        }
                    }


                    char *position = (char *)data + (my_start - global_start_r) * type_size;
                    memcpy(position, source_addr2, my_size);


fprintf(stderr, "%d rank=%d name = %s position=%p data=%p mystart = %lu global_start_r = %lu mysize = %lu %s:%d\n", ds->my_rank,chunk->rank,var->varname,position,data,my_start, global_start_r, my_size, __FILE__,__LINE__);
                    total_size += my_size;
fprintf(stderr, "im here %d %s:%d\n", ds->my_rank,__FILE__,__LINE__);
                    break;
                }
                else {
fprintf(stderr, "im here %s:%d\n", __FILE__,__LINE__);
                    // TODO
                    continue;
                }
            }
            else {
                // TODO
                adios_error(err_invalid_read_method, "Datatap cannot support this decomposition.\n");
                return -1;
            }
        }
        // we are done with this chunk, now extract data from next chunk
#endif
        

        current_chunk = current_chunk->next;
    }
fprintf(stderr, "im here read rank %d var %s total size %ld %s:%d\n", ds->my_rank, var->varname, total_size,__FILE__,__LINE__);

    // TODO: the data size is messed up because of ghost zone
    total_size = type_size;
    int j;
    for(j = var->ndims-1; j >= 0; j --) {
        total_size *= count[j]; 
    }
    return total_size;       
}

#if 0
int reorganize_array (QueueItem *qi, datatap_read_file_data *ds)
{
fprintf(stderr, "im here rank %d %s:%d\n",dt_read_data->dt_comm_rank,__FILE__,__LINE__);
    // go through all variables user wants to read
    datatap_pg_info *pg = NULL;
    int i;
    for(i = 0; i < ds->num_pgs; i ++) {
        if(ds->pgs[i].rank == qi->rank) {
            pg = &(ds->pgs[i]);                   
            break;
        }
    }
    if(!pg) {
        adios_error(err_unspecified, "cannot find pg.\n");
        return -1;        
    }
    
fprintf(stderr, "im here rank %d %s:%d\n",dt_read_data->dt_comm_rank,__FILE__,__LINE__);
    datatap_var_info *var = ds->vars_read;
    FMFieldList filed_list = qi->var_list->field_list; // TODO
    while(var) {
        // first search within the pg for this var        
fprintf(stderr, "im here rank %d var name %s %s:%d\n",dt_read_data->dt_comm_rank,var->varname,__FILE__,__LINE__);
        datatap_var_info *v = NULL;
        for(i = 0; i < pg->num_vars; i ++) {
fprintf(stderr, "im here rank %d var name %s i %d %s:%d\n",dt_read_data->dt_comm_rank,var->varname,i,__FILE__,__LINE__);
            //if(!strcmp(pg->vars[i].varname, var->varname)) {
            if(!compare_var_name(var->varname, &(pg->vars[i]))) {
fprintf(stderr, "im here rank %d var name %s i %d %s:%d\n",dt_read_data->dt_comm_rank,var->varname,i,__FILE__,__LINE__);
                // got it
                v = &(pg->vars[i]);
fprintf(stderr, "im here rank %d var name %s i %d %s:%d\n",dt_read_data->dt_comm_rank,var->varname,i,__FILE__,__LINE__);
                 
                int type_size = common_read_type_size(v->type, NULL); 
fprintf(stderr, "im here rank %d var name %s i %d %s:%d\n",dt_read_data->dt_comm_rank,var->varname,i,__FILE__,__LINE__);
                FMField *f = find_field(v->varname, v->varpath, filed_list);
fprintf(stderr, "im here rank %d var name %s i %d %s:%d\n",dt_read_data->dt_comm_rank,var->varname,i,__FILE__,__LINE__);
                char *source_addr = (char *)qi->data + f->field_offset; // this is the address
                void *source_addr2;
fprintf(stderr, "im here rank %d var name %s i %d %s:%d\n",dt_read_data->dt_comm_rank,var->varname,i,__FILE__,__LINE__);
                switch(v->type) {
                    case adios_byte:
                    case adios_unsigned_byte:
                        source_addr2 = *((char **)source_addr);
                        break;
          
                    case adios_string: // TODO
                        return strlen((char *)source_addr) + 1;
                        source_addr2 = *((char **)source_addr);
                        break;
         
                    case adios_short:
                    case adios_unsigned_short:
                        source_addr2 = *((short int **)source_addr);
                        break;
        
                    case adios_integer:
                    case adios_unsigned_integer:
                        source_addr2 = *((int **)source_addr);
                        break;
        
                    case adios_long:
                    case adios_unsigned_long:
                        source_addr2 = *((long int **)source_addr);
                        break;
        
                    case adios_real:
                        source_addr2 = *((float **)source_addr);
                        break;
        
                    case adios_double:
                        source_addr2 = *((double **)source_addr);
                        break;
        
                    case adios_long_double:
                        source_addr2 = *((long double **)source_addr);
                        break;
        
                    case adios_complex:
                    case adios_double_complex:
                    default:
                        adios_error(err_invalid_read_method, "complex data type is not supported.");
                        return -1;
                }
                //memcopy(var->read_buffer, source_addr2, var->read_buffer_size);
                                      
                // now we need to calculate how to map chunks into read buffer
                // according to dimension specification
                // also we need to figure out which part should be shuffled to
                // other peer processes and which part should be moved in from
                // other processes
                // search along the slowest changing dimension and determine if 
                // the sub-chunk falls into the reading region. If so, copy the 
                // whole sub-chunk, otherwise search along the second slowest
                // changing dimension within that subchunk, identify which part
                // should go where and copy the remaining part to local buffer
fprintf(stderr, "im here rank %d var name %s i %d %s:%d\n",dt_read_data->dt_comm_rank,var->varname,i,__FILE__,__LINE__);
                int j;
                uint64_t global_start_r = 0; 
                for(j = var->ndims-1; j >= 0; j --) {
                    if(var->global_bounds[j]) {
                        global_start_r *= var->global_bounds[j];
                        global_start_r += var->global_offsets[j];
                    }
                }
                
if(!strcmp(v->varname,"bconds")){
    int *bconds = (int *)source_addr2;
    int t;
    for(t=0;t<48;t++) {
        fprintf(stderr, "bconds[%d]=%d\n",t,bconds[t]);
    }

}

if(!strcmp(v->varname,"v1")){
    double *v1 = (double *)source_addr2;
    int t;
    for(t=0;t<100;t++) {
        fprintf(stderr, "v1[%d]=%f\n",t,v1[t]);
    }

}

fprintf(stderr, "im here rank %d var name %s i %d %s:%d\n",dt_read_data->dt_comm_rank,var->varname,i,__FILE__,__LINE__);
                for(j = v->ndims-1; j >= 0; j --) {
                    // for each sub-chunk, test if it falls into the reading region
                    int lower = (v->global_offsets[j] >= var->global_offsets[j]) ? 1: 0;
                    int higher = ((v->global_offsets[j]+v->local_bounds[j]) <= 
                         (var->global_offsets[j]+var->local_bounds[j])) ? 1: 0; 
fprintf(stderr, "im here rank=%d j %d name %s go %lu lb %lu go %lu lb %lu %s:%d\n", dt_read_data->dt_comm_rank,j,var->varname,v->global_offsets[j],v->local_bounds[j],var->global_offsets[j],var->local_bounds[j],__FILE__,__LINE__);
                    if(lower && higher) { 
fprintf(stderr, "im here j %d %s:%d\n", j,__FILE__,__LINE__);
                        if(j == 0) { // copy the whole chunks
                            uint64_t my_start = 0;
                            uint64_t my_size = type_size;
                            int k;
                            for(k = v->ndims-1; k >= 0; k --) {
                                my_start *= v->global_bounds[k];
                                my_start += v->global_offsets[k];
                                my_size *= v->local_bounds[k];  
                            }
                            char *position = (char *)var->data + (my_start - global_start_r) * type_size;                            
fprintf(stderr, "rank=%d name = %s position=%lu data=%lu mystart = %lu global_start_r = %lu mysize = %lu %s:%d\n", dt_read_data->dt_comm_rank,var->varname,position,var->data,my_start, global_start_r, my_size, __FILE__,__LINE__);
                            memcpy(position, source_addr2, my_size);      
fprintf(stderr, "im here j %d %s:%d\n", j,__FILE__,__LINE__);
                            break;
                        }
                        else {
                            // TODO 
                            continue; 
                        }
                    }
                    else {
fprintf(stderr, "im here rank=%d j %d name %s go %lu lb %lu go %lu lb %lu %s:%d\n", dt_read_data->dt_comm_rank,j,var->varname,v->global_offsets[j],v->local_bounds[j],var->global_offsets[j],var->local_bounds[j],__FILE__,__LINE__);
fprintf(stderr, "im here j %d %s:%d\n", j,__FILE__,__LINE__);
                        // TODO
                        adios_error(err_invalid_read_method, "Datatap cannot support this decomposition.\n");
                        return -1;
                    } 
fprintf(stderr, "im here j %d %s:%d\n", j,__FILE__,__LINE__);
                }
                break;
                // end of data re-organization for this var
            } 
        }             
        
        // now process the next var to read        
        var = var->next;
    }    
fprintf(stderr, "im here rank %d %s:%d\n",dt_read_data->dt_comm_rank,__FILE__,__LINE__);

    return 0;
}
#endif

void data_handler (void *data, int length, void *user_data, attr_list attr, int rank, void *timing, file_info *f)
{
    recvtime *r = (recvtime*)timing;
    IOhandle *h = (IOhandle*)user_data;
    elapsedtime *e = updateTimes(h, r, length);
    
    // decode the data and insert the data into the queue
    int decoded_length = FFS_est_decode_length(h->iocontext, data, length);

    // TODO: make sure we free this later (after writing to hdf5 file)
    char *decoded_data = (char *)malloc(decoded_length);

    if(!decoded_data) {
        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
        exit(-1);
    }

    FFSTypeHandle ffshandle = FFSTypeHandle_from_encode(h->iocontext, data);
    FMFormat form = FMFormat_of_original(ffshandle);

    // TODO
    FMStructDescList var_list = get_localized_formats(form);
    establish_conversion(h->iocontext, ffshandle, var_list);
    FFSdecode_to_buffer(h->iocontext, data, decoded_data);

    // The encoded data can be recycled now
    returnbuffer(h, data, length);

    QueueItem * qi =(QueueItem *) malloc(sizeof(QueueItem));
    if(!qi) {
        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
        exit(-1);
    }
    qi->data = decoded_data;
    qi->length = decoded_length;
    qi->var_list = var_list;
    qi->rank = rank; // TODO: hack it!   

    // put message in queue
    pthread_mutex_lock(&(f->mutex));
//    while(queue_size(f->dt_queue) >= f->dt_queue_max_length) {
//        // wait until queue is not full
//        pthread_cond_signal(&(dt_read_data->cond2));
//        pthread_cond_wait(&(dt_read_data->cond1), &(dt_read_data->mutex));
//    }
    queue_enqueue(f->dt_queue, qi);
    if(queue_size(f->dt_queue) == f->num_chunks) {
        pthread_cond_signal(&(f->cond2));
    }
    pthread_mutex_unlock(&(f->mutex));


fprintf(stderr, "im here rank= %d client rank=%d data handler done %s:%d\n",h->rank,rank,__FILE__,__LINE__);
    
    // TODO: check if it's time to stop
}

void * dt_server_thread_func (void *arg)
{
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    MPI_Comm orig_comm = (MPI_Comm) arg;

#if 0 // open-mpi   
    // duplicate a MPI communicator for synchronization between dt servers
    int rc = MPI_Comm_dup(orig_comm, &(dt_read_data->dt_comm));
    if(rc != MPI_SUCCESS) {
        error(err_unspecified, "Cannot duplicate communicator for Datatap.");
        pthread_exit(NULL);
    }
#else 
    dt_read_data->dt_comm = orig_comm;
#endif

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);

#ifdef _NOMPI
    dt_read_data->dt_comm_size = 1;
    dt_read_data->dt_comm_rank = 0;
#else 
    MPI_Comm_size(dt_read_data->dt_comm, &(dt_read_data->dt_comm_size));
    MPI_Comm_rank(dt_read_data->dt_comm, &(dt_read_data->dt_comm_rank));
#endif
  
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    // initialize ptlpbio interface
    dt_read_data->dt_cm = CManager_create();
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    CMlisten_specific(dt_read_data->dt_cm, NULL);

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    lrand48();

    dt_read_data->dt_handle = EVthin_portals_listen(dt_read_data->dt_cm, 120,
                              0, data_handler, dt_read_data->dt_comm);

    char param_file[30];
    int appid, was_set;
    appid = globals_adios_get_application_id(&was_set);
    if(!was_set) {
        adios_error(err_unspecified, "Application ID was not set.");
        sprintf(param_file, "datatap_param\0");
    }
    else {
        sprintf(param_file, "datatap_param%d\0", appid);
    }

    // dt server(rank 0) gather contact info from other servers and write to
    // a file so upstream writers can connect to this application
    outputConnectionData(param_file, dt_read_data->dt_handle);

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    dt_read_data->dt_server_ready = 1;

    // serve the network
    CMrun_network(dt_read_data->dt_cm);

    // TODO: cleanup and exit
    CManager_close(dt_read_data->dt_cm);
    return NULL;
}

int adios_read_datatap_init (MPI_Comm comm)
{
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    setenv("CMSelfFormats", "1", 1);

    // initialize Datatap read method structure
    dt_read_data = (datatap_read_method_data *) malloc(sizeof(datatap_read_method_data));
    if(!dt_read_data) {
        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
        return -1;
    }    
    memset(dt_read_data, 0, sizeof(datatap_read_method_data));

    // enable threading for EVPath
    gen_pthread_init();

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
#if 0
    // set up queue for incoming data chunks
    dt_read_data->dt_queue =(Queue *) calloc(1, sizeof(Queue));
    if(!dt_read_data->dt_queue) {
        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
        return -1;
    }
    queue_init(dt_read_data->dt_queue, free);
    dt_read_data->dt_queue_max_length = DT_MAX_QUEUE_LENGTH;

    pthread_mutex_init(&(dt_read_data->mutex), NULL);
    pthread_cond_init(&(dt_read_data->cond1), NULL);
    pthread_cond_init(&(dt_read_data->cond2), NULL);
#endif    

    // fork the thread to poll network 
    int rc = pthread_create(&(dt_read_data->dt_server_thread), NULL, 
                            dt_server_thread_func, (void *)comm);
    if(rc) {
        adios_error(err_unspecified, "Failed to create Datatap server thread.");
        return -1; 
    }

    // TODO: wait until the dt server thread is ready
    while(!dt_read_data->dt_server_ready) { }
    
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    return 0;
}

int adios_read_datatap_finalize ()
{
    // notify and wait for dt server thread to exit
    datatap_stop_server(dt_read_data->dt_handle);
    pthread_join(dt_read_data->dt_server_thread, NULL);

#if 0
    // TODO: we need a datatap cleanup function
    pthread_mutex_destroy(&(dt_read_data->mutex));
    pthread_cond_destroy(&(dt_read_data->cond1));
    pthread_cond_destroy(&(dt_read_data->cond2));    
    free(dt_read_data->dt_queue);
#endif

    free(dt_read_data);
    return 0;
}

ADIOS_FILE *adios_read_datatap_fopen(const char *fname, MPI_Comm comm)
{
    ADIOS_FILE *fp;
    adios_errno = 0;

    datatap_read_file_data *ds = (datatap_read_file_data *) malloc(sizeof(datatap_read_file_data));
    if(!ds) {
        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
        return NULL;                
    }
    ds->file_name = strdup(fname);    
     
    // here we need to syncrhonize with other reader processes to see what they find
    // 1: file available locally
    // 0: file not found
    int my_status, global_status;
    int total_num_readers, my_rank;

#ifdef _NOMPI
    total_num_readers = 1;
    my_rank = 0;
#else
    MPI_Comm_size(comm, &total_num_readers);
    MPI_Comm_rank(comm, &my_rank);
#endif

    ds->comm = comm;
    ds->my_rank = my_rank;
    ds->comm_size = total_num_readers;

    // first we need to make sure the file has been 'written'
    // TODO: always start from timstep 0
    int is_EOF = 0;
    ds->f_info = datatap_get_file_info(dt_read_data->dt_handle, fname, dt_read_data->num_io_dumps, &is_EOF);
    if(ds->f_info) {
        // now I should tell my peer readers that I have seen this file available locally
        my_status = 1;
        if(total_num_readers > 1) {
            int rc = MPI_Allreduce(&my_status, &global_status, 1, MPI_INT, MPI_SUM, comm);
            if(rc != MPI_SUCCESS) {
                fprintf(stderr, "something bad happened somewhere.\n");
                free(ds->file_name);
                free(ds);
                return NULL;
            } 
        }
        // now I can move on to process meta-data
fprintf(stderr, "im here rank %d move on %s:%d\n", my_rank, __FILE__,__LINE__);
    }
    else {
        if(is_EOF) {
            adios_errno = err_end_of_file;
            adios_error(err_end_of_file, "Reach the end of file (%s).", fname);
            free(ds->file_name);
            free(ds);
    
            if(total_num_readers > 1) {               
                // now I should wait for my peer readers to see EOF in their local context
                int rc = MPI_Barrier(comm);
                if(rc != MPI_SUCCESS) {
                    fprintf(stderr, "something bad happened somewhere.\n");
                }
                return NULL; 
            }
        }
        else {
            if(total_num_readers > 1) {

                // I didn't see the file available in my local context, but my peer readers may
                // have seen it, so I should ask them to figure out
fprintf(stderr, "im here rank %d ask %s:%d\n", my_rank, __FILE__,__LINE__);
                my_status = 0;
                MPI_Allreduce(&my_status, &global_status, 1, MPI_INT, MPI_SUM, comm);
                if(global_status == 0) { // no one see file available
fprintf(stderr, "im here rank %d no file %s:%d\n", my_rank, __FILE__,__LINE__);
                    adios_errno = err_file_not_found_error;
                    adios_error(err_file_not_found_error, "Cannot find file (%s).", fname);
fprintf(stderr, "im here rank %d no file %s:%d\n", my_rank, __FILE__,__LINE__);
                    free(ds->file_name);
                    free(ds);
                    return NULL;
                }
                else { 
                    // some one has seen this file available, so I should keep polling instead of return
                    do {
fprintf(stderr, "im here rank %d %s:%d\n", my_rank, __FILE__,__LINE__);
                        is_EOF = 0;
                        ds->f_info = datatap_get_file_info(dt_read_data->dt_handle, 
                            fname, dt_read_data->num_io_dumps, &is_EOF);
                        if(ds->f_info) { 
                            // now we get it!
fprintf(stderr, "im here rank %d %s:%d\n", my_rank, __FILE__,__LINE__);
                            break;
                        }
                        else {
                            if(is_EOF) {
                                adios_errno = err_end_of_file;
                                adios_error(err_end_of_file, "Reach the end of file (%s).", fname);
                                free(ds->file_name);
                                free(ds);

                                int rc = MPI_Barrier(comm);
                                if(rc != MPI_SUCCESS) {
                                    fprintf(stderr, "something bad happened somewhere.\n");
                                }
                                return NULL;
                            }
                            else {
                                usleep(10000); // TODO: we need to set a proper value
fprintf(stderr, "im here rank %d %s:%d\n", my_rank, __FILE__,__LINE__);
                                continue;
                            }
                        }  
                    } 
                    while (1);
                }
            }
            else {
                adios_errno = err_file_not_found_error;
                adios_error(err_file_not_found_error, "Cannot find file (%s).", fname);
fprintf(stderr, "im here rank %d no file %s:%d\n", my_rank, __FILE__,__LINE__);
                free(ds->file_name);
                free(ds);
                return NULL;
            }  
        }
    }       
fprintf(stderr, "im here rank %d see file %s:%d\n", my_rank, __FILE__,__LINE__);
    
    // we don't know what to read yet
    ds->num_vars = 0;
    ds->vars = NULL;
//    ds->num_vars_read = 0;
//    ds->vars_read = NULL;
//    ds->vars_read_tail = NULL;

    // TODO: add a loop over chunks
    int i;
    for(i = 0; i < ds->f_info->num_chunks; i ++) {  
        chunk_info *current_chunk = &(ds->f_info->chunks[i]);  
        
fprintf(stderr, "im here %s:%d\n", __FILE__,__LINE__);
        // parse the var info
        char *current_pos = current_chunk->var_info_buffer;
        int total_var_info_size = *(int *)current_pos; // total size
        current_pos += 4;

        if(i == 0) {
            ds->host_language_fortran = *(enum ADIOS_FLAG *)current_pos; // host language
        }
        current_pos += 4;

        int group_name_len = *(int *)current_pos; // size of group name
        current_pos += 4; 
  
        if(i == 0) { // TODO: let's assume one group name
            ds->group_name = strdup(current_pos); // group name
        }
        current_pos += group_name_len;

        int num_vars_in_pg = *(int *)current_pos; // total num of vars  
        current_pos += 4;
fprintf(stderr, "im here info size %d num vars %d %s:%d\n", total_var_info_size,num_vars_in_pg,__FILE__,__LINE__);

        datatap_var_info *current_var;
        char *end = current_chunk->var_info_buffer + total_var_info_size;
        while(current_pos < end) {
            // size of this var info
            int var_info_size = *(int *) current_pos;             
            current_pos += sizeof(int);        
            int var_id = *(int *) current_pos;
            current_pos += sizeof(int);        
            int varname_len = *(int *) current_pos;        
            current_pos += sizeof(int);        
            char *varname = current_pos;  
            current_pos += varname_len;
            int varpath_len = *(int *) current_pos;
            current_pos += sizeof(int);
            char *varpath = current_pos;
            current_pos += varpath_len;

            // now we go through the current list of vars to see if this is new
            current_var = ds->vars;
            while(current_var != NULL) {
                // TODO: compare var id
                //if(!strcmp(current_var->varname, varname) && 
                //    !strcmp(current_var->varpath, varpath)) {
                if(var_id == current_var->id) {
                    break;
                }
                else {
                    current_var = current_var->next;
                }
            }

            if(!current_var) { // this is a new var
                current_var = (datatap_var_info *) malloc(sizeof(datatap_var_info));
                if(!current_var) {
                    adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                    return NULL;
                }
                current_var->id = var_id;
                current_var->varname = strdup(varname);
                current_var->varpath = strdup(varpath);
                current_var->type = *(enum ADIOS_DATATYPES *) current_pos;
                current_pos += sizeof(enum ADIOS_DATATYPES);
                current_var->time_dim = *(int *) current_pos;
                current_pos += sizeof(int);
                current_var->ndims = *(int *) current_pos;
                current_pos += sizeof(int);
                current_var->num_chunks = 0;
                current_var->chunks = NULL;
                if(!current_var->ndims) { // scalars and strings
                    current_var->data_size = common_read_type_size(current_var->type, current_pos);
                }
fprintf(stderr, "im here %s %s %d %s:%d\n", current_var->varname, current_var->varpath, current_var->ndims, __FILE__,__LINE__);
                current_var->next = ds->vars;
                ds->vars = current_var;
                ds->num_vars ++;
            }
            else {
                current_pos += sizeof(enum ADIOS_DATATYPES);
                //current_var->time_dim = *(int *) current_pos;
                current_pos += sizeof(int);
                //current_var->ndims = *(int *) current_pos;
                current_pos += sizeof(int);
            }            
            
            datatap_var_chunk *new_chunk = (datatap_var_chunk *) malloc(sizeof(datatap_var_chunk));             
            if(!new_chunk) {
                    adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                    return NULL;
            }
            new_chunk->next = current_var->chunks;
            current_var->chunks = new_chunk;
            current_var->num_chunks ++; 

            new_chunk->rank = current_chunk->rank;
            if(!current_var->ndims) { // scalars and strings
                // copy data value             
                new_chunk->data = malloc(current_var->data_size);
                if(!new_chunk->data) {
                    adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                    return NULL;                                  
                }
                memcpy(new_chunk->data, current_pos, current_var->data_size);         
                current_pos += current_var->data_size;
            }
            else { // arrays
                new_chunk->local_bounds = (uint64_t *) malloc(current_var->ndims * sizeof(uint64_t));  
                if(!new_chunk->local_bounds) {
                    adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                    return NULL;                                  
                }
                new_chunk->global_bounds = (uint64_t *) malloc(current_var->ndims * sizeof(uint64_t));  
                if(!new_chunk->global_bounds) {
                    adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                    return NULL;                                  
                }
                new_chunk->global_offsets = (uint64_t *) malloc(current_var->ndims * sizeof(uint64_t));  
                if(!new_chunk->global_offsets) {
                    adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                    return NULL;                                  
                }
                int i;
                for(i = 0; i < current_var->ndims; i ++) {
                    new_chunk->local_bounds[i] = *(uint64_t *)current_pos;
                    current_pos += sizeof(uint64_t);  
                    new_chunk->global_bounds[i] = *(uint64_t *)current_pos;
                    current_pos += sizeof(uint64_t);  
                    new_chunk->global_offsets[i] = *(uint64_t *)current_pos;
                    current_pos += sizeof(uint64_t);  
                }                          
            }
        }    
        // TODO: at this point, we no longer need var_info_buffer so free it
        // free(current_chunk->var_info_buffer);
    }
        

    // TODO: replicate meta-data among peer reader
    ds->num_vars_peer = 0;
    ds->vars_peer = NULL;

    if(total_num_readers > 1)
    {
        char *send_buf = NULL;
        char *recv_buf = NULL; 
        int send_count = VAR_INFO_SIZE, recv_count = VAR_INFO_SIZE;
        if(my_rank == 0) {
            for(i = 0; i < ds->f_info->num_chunks; i ++) {
                chunk_info *current_chunk = &(ds->f_info->chunks[i]);
fprintf(stderr, "im here rank %d %d %s:%d\n", ds->my_rank, current_chunk->rank,__FILE__,__LINE__);
                if(current_chunk->rank == 0) {
                    send_buf = ds->f_info->chunks[i].var_info_buffer;
                    break;
                } 
            }
            
            // TODO: before we do this, we need to make sure this one has sth special
            // the only case we need to deal with is local array with only one chunk
            // we will check this when calling read_var

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
            int rc = MPI_Bcast(send_buf, send_count, MPI_BYTE, 0, comm);
            if(rc != MPI_SUCCESS) {
                fprintf(stderr, "rank %d: MPI_Scatter returns error (%d). %s:%d\n",
                    my_rank, rc, __FILE__, __LINE__);
                return NULL;
            }
        }
        else {
            recv_buf = (char *) malloc(VAR_INFO_SIZE); 
            if(!recv_buf) {
                adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                return NULL;
            }

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
            int rc = MPI_Bcast(recv_buf, recv_count, MPI_BYTE, 0, comm);
            if(rc != MPI_SUCCESS) {
                fprintf(stderr, "rank %d: MPI_Scatter returns error (%d). %s:%d\n",
                    my_rank, rc, __FILE__, __LINE__);
                return NULL;
            }
        }

        if(my_rank != 0) {
            // parse the var_info buffer
            char *current_pos = recv_buf;
            int total_var_info_size = *(int *)current_pos; // total size
            current_pos += 4;

            current_pos += 4; // host language

            int group_name_len = *(int *)current_pos; // size of group name
            current_pos += 4;
            current_pos += group_name_len;

            int num_vars = *(int *)current_pos; // total num of vars
            current_pos += 4;

            datatap_var_info *current_var;
            char *end = recv_buf + total_var_info_size;
            while(current_pos < end) {
                // size of this var info
                int var_info_size = *(int *) current_pos;
                char *start_of_next_var = current_pos + var_info_size;
                current_pos += sizeof(int);
                int var_id = *(int *) current_pos;
                current_pos += sizeof(int);
                int varname_len = *(int *) current_pos;
                current_pos += sizeof(int);
                char *varname = current_pos;
                current_pos += varname_len;
                int varpath_len = *(int *) current_pos;
                current_pos += sizeof(int);
                char *varpath = current_pos;
                current_pos += varpath_len;

                // now we go through the current list of vars to see if this is new
                current_var = ds->vars;
                while(current_var != NULL) {
                    //if(!strcmp(current_var->varname, varname) &&
                    //    !strcmp(current_var->varpath, varpath)) {
                    if(var_id == current_var->id) {
                        // this var is locally available, so skip it
                        break;
                    }
                    else {
                        current_var = current_var->next;
                    }
                }

                if(current_var) { // this var is locally available
                    current_pos = start_of_next_var;
                    continue;
                }

                current_var = (datatap_var_info *) malloc(sizeof(datatap_var_info));
                if(!current_var) {
                    adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                    return NULL;
                }
                current_var->id = var_id;
                current_var->varname = strdup(varname);
                current_var->varpath = strdup(varpath);
                current_var->type = *(enum ADIOS_DATATYPES *) current_pos;
                current_pos += sizeof(enum ADIOS_DATATYPES);
                current_var->time_dim = *(int *) current_pos;
                current_pos += sizeof(int);
                current_var->ndims = *(int *) current_pos;
                current_pos += sizeof(int);
                current_var->num_chunks = 0;
                current_var->chunks = NULL;
                if(!current_var->ndims) { // scalars and strings
                    current_var->data_size = common_read_type_size(current_var->type, current_pos);
                }
                current_var->next = ds->vars_peer;
                ds->vars_peer = current_var;
                ds->num_vars_peer ++;

                datatap_var_chunk *new_chunk = (datatap_var_chunk *) malloc(sizeof(datatap_var_chunk));
                if(!new_chunk) {
                    adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                    return NULL;
                }
                new_chunk->next = current_var->chunks;
                current_var->chunks = new_chunk;
                current_var->num_chunks ++;

                new_chunk->rank = 0;
                if(!current_var->ndims) { // scalars and strings
                    // copy data value
                    new_chunk->data = malloc(current_var->data_size);
                    if(!new_chunk->data) {
                        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                        return NULL;
                    }
                    memcpy(new_chunk->data, current_pos, current_var->data_size);
                    current_pos += current_var->data_size;
                }
                else { // arrays
                    new_chunk->local_bounds = (uint64_t *) malloc(current_var->ndims * sizeof(uint64_t));
                    if(!new_chunk->local_bounds) {
                        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                        return NULL;
                    }
                    new_chunk->global_bounds = (uint64_t *) malloc(current_var->ndims * sizeof(uint64_t));
                    if(!new_chunk->global_bounds) {
                        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                        return NULL;
                    }
                    new_chunk->global_offsets = (uint64_t *) malloc(current_var->ndims * sizeof(uint64_t));
                    if(!new_chunk->global_offsets) {
                        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                        return NULL;
                    }
                    int i;
                    for(i = 0; i < current_var->ndims; i ++) {
                        new_chunk->local_bounds[i] = *(uint64_t *)current_pos;
                        current_pos += sizeof(uint64_t);
                        new_chunk->global_bounds[i] = *(uint64_t *)current_pos;
                        current_pos += sizeof(uint64_t);
                        new_chunk->global_offsets[i] = *(uint64_t *)current_pos;
                        current_pos += sizeof(uint64_t);
                    }
                }
            }
            
            free(recv_buf);   
        }
    }

fprintf(stderr, "im here %d %d %d %s:%d\n", ds->host_language_fortran, adios_flag_yes, futils_is_called_from_fortran(), __FILE__, __LINE__);
    // TODO: here we need to ajust array dimension if the reader is in C 
    // but writer is in Fortran or vice versa
    if(ds->host_language_fortran == adios_flag_yes && !futils_is_called_from_fortran()) {
        // reader is in C but writer is in Fortran, there are several things to adjust:
        // array dimension index starts from 1 in Fortran --> start from 0 in C
        // change array dimension order (including time dimension)
        datatap_var_info *v = ds->vars;        
        while(v) {
            datatap_var_chunk *chunk = v->chunks; 
            while(chunk) {
                int i;
                uint64_t temp;
                for(i = 0; i < v->ndims/2; i ++) {
                    temp = chunk->local_bounds[v->ndims-i-1];
                    chunk->local_bounds[v->ndims-i-1] = chunk->local_bounds[i];
                    chunk->local_bounds[i] = temp;
                    temp = chunk->global_bounds[v->ndims-i-1];
                    chunk->global_bounds[v->ndims-i-1] = chunk->global_bounds[i];
                    chunk->global_bounds[i] = temp;
                    temp = chunk->global_offsets[v->ndims-i-1];
                    chunk->global_offsets[v->ndims-i-1] = chunk->global_offsets[i];
                    chunk->global_offsets[i] = temp;
                }
                chunk = chunk->next;
            }
            if(v->time_dim > 0) { // -1 means no time dimension
                v->time_dim = v->ndims - v->time_dim;  
            }
            v = v->next;
        }

        v = ds->vars_peer;
        while(v) {
            datatap_var_chunk *chunk = v->chunks;
            while(chunk) {
                int i;
                for(i = 0; i < v->ndims; i ++) {
                    chunk->local_bounds[i] --;
                    //chunk->global_offsets[i] --;
                }
                uint64_t temp;
                for(i = 0; i < v->ndims/2; i ++) {
                    temp = chunk->local_bounds[v->ndims-i-1];
                    chunk->local_bounds[v->ndims-i-1] = chunk->local_bounds[i];
                    chunk->local_bounds[i] = temp;
                    temp = chunk->global_bounds[v->ndims-i-1];
                    chunk->global_bounds[v->ndims-i-1] = chunk->global_bounds[i];
                    chunk->global_bounds[i] = temp;
                    temp = chunk->global_offsets[v->ndims-i-1];
                    chunk->global_offsets[v->ndims-i-1] = chunk->global_offsets[i];
                    chunk->global_offsets[i] = temp;
                }
                chunk = chunk->next;
            }
            v->time_dim = v->ndims - v->time_dim - 1;
            v = v->next;
        }
    }
    else if(ds->host_language_fortran == adios_flag_no && futils_is_called_from_fortran()) {
        // adjuct dimension C --> Fortran  
        // TODO: for the demo, this will not happen so leave it here
    }

    fp = (ADIOS_FILE *) malloc(sizeof(ADIOS_FILE));
    if(!fp) {
        adios_error(err_no_memory, "Cannot allocate memory for file info.");
        free(ds);
        return NULL;
    }
 
    fp->fh = (uint64_t) ds;
    fp->groups_count = 1; // TODO: assume one group per file
    fp->vars_count = 0; // TODO: just do not use this filed
    fp->attrs_count = 0; // TODO: do not support attributes yet

    // TODO: since we require fopen/fclose for each timestep, 
    // so there is always only 1 timestep in file
    // TODO: set ntimesteps to be max to make pixie3d working
    fp->tidx_start = ds->f_info->timestep;  
    fp->ntimesteps = INT32_MAX;
fprintf(stderr, "im here timestep %d %s:%d\n", fp->tidx_start = ds->f_info->timestep, __FILE__,__LINE__);

    fp->file_size = 0; 
    fp->version = 1;
    fp->endianness = 0; 
    alloc_namelist(&fp->group_namelist, fp->groups_count); 
    for (i = 0; i < fp->groups_count; i++) {
        if (!fp->group_namelist[i])  {
            adios_error(err_no_memory, "Cannot allocate buffer in adios_fopen()");
            return NULL;
        }
        else {
            strcpy(fp->group_namelist[i], ds->group_name); 
        }
    }
fprintf(stderr, "im here %s:%d\n", __FILE__,__LINE__);

    // TODO: return code of adios_errno for fopen:
    // 0: file metadata is available and everything is success
    // we need a value telling that the file is not available yet
    // we also need a value telling that writer is finishing so don't
    // read this file any more

    return fp; 
}

int adios_read_datatap_fclose(ADIOS_FILE *fp)
{
fprintf(stderr, "im here rank %d %s:%d\n",dt_read_data->dt_comm_rank,__FILE__,__LINE__);
    datatap_read_file_data *ds = (datatap_read_file_data *) fp->fh;

    adios_errno = 0;

    dt_read_data->num_io_dumps ++;

    // recycle the queue
    QueueItem *qi;

    while((!queue_dequeue(ds->f_info->dt_queue, &qi)) && qi) {
        free(qi->data);
        free(qi);  
    }

    // now we 'read' the data into user-provided buffers
    // note that the data is actually moved by dt servr thread
    // here we just re-organize the data into the distribution
    // user wants

#if 0
    // process the incoming data chunks 
    void *d = NULL;
    QueueItem *qi;
    uint64_t num_chunks_processed = 0;
    int file_done = 0;

    while(1) {
        pthread_mutex_lock(&(ds->f_info->mutex));
                while(queue_size(ds->f_info->dt_queue) == 0) {
        fprintf(stderr, "im here num_chunks_processed %d %d %s:%d\n",num_chunks_processed, ds->num_pgs,__FILE__,__LINE__);
            // check if it's time to finish this file
            if(num_chunks_processed == ds->num_pgs) { // TODO
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
                dt_read_data->num_io_dumps ++;
                pthread_mutex_unlock(&(ds->f_info->mutex));
                file_done = 1;
                break;
            }
            else { 
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
                pthread_cond_signal(&(ds->f_info->cond1));
                pthread_cond_wait(&(ds->f_info->cond2), &(ds->f_info->mutex));
            }
        }
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);

        if(file_done) {
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
            break;
        }

        // start of busy time
        while((!queue_dequeue(ds->f_info->dt_queue, &d)) && d) {
            //pthread_cond_signal(&(ds->f_info->cond1));
            pthread_mutex_unlock(&(ds->f_info->mutex));

            qi =(QueueItem *) d;

fprintf(stderr, "im here rank %d %s:%d\n",dt_read_data->dt_comm_rank,__FILE__,__LINE__);
            int rc = reorganize_array (qi, ds);
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
            if(rc) {
                adios_error(err_unspecified, "Error in reorganize_array() function.");
                exit(-1);
            }

            // recycle decoded data buffer
            free(qi->data);
            free(qi);

            num_chunks_processed ++;

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
            pthread_mutex_unlock(&(ds->f_info->mutex));
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
            pthread_mutex_lock(&(ds->f_info->mutex));
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
        }
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);

        pthread_mutex_unlock(&(ds->f_info->mutex));
    }
#endif    

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    // notify dt server thread that we are done with this file 
    int rc = datatap_release_file(dt_read_data->dt_handle, ds->f_info);
    if(rc != 0) {
        adios_error(err_unspecified, "Could not close file.");
        return -1;
    }

    free_namelist((fp->group_namelist), fp->groups_count);
    if (ds->file_name) { 
        free(ds->file_name); 
        ds->file_name = NULL; 
    }
    
    if (ds->group_name) {
        free(ds->group_name);
        ds->group_name = NULL;
    }

    // TODO release pg_info and var_info 
    datatap_var_info *v = ds->vars;
    datatap_var_info *tmp_v;
    while(v) {
        tmp_v = v;
        free(v->varname);
        free(v->varpath);
        datatap_var_chunk *chunk = v->chunks;
        datatap_var_chunk *tmp_chunk;
        while(chunk) {
            tmp_chunk = chunk;
            if(v->ndims) {
                free(chunk->local_bounds);
                free(chunk->global_bounds);
                free(chunk->global_offsets);
            }
            else {
                free(chunk->data);     
            }
            chunk = chunk->next;
            free(tmp_chunk);
        }
        v = v->next;
        free(tmp_v);
    }

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    v = ds->vars_peer;
    while(v) {
        tmp_v = v;
        free(v->varname);
        free(v->varpath);
        datatap_var_chunk *chunk = v->chunks;
        datatap_var_chunk *tmp_chunk;
        while(chunk) {
            tmp_chunk = chunk;
            if(v->ndims) {
                free(chunk->local_bounds);
                free(chunk->global_bounds);
                free(chunk->global_offsets);
            }
            else {
                free(chunk->data);
            }
            chunk = chunk->next;
            free(tmp_chunk);
        }
        v = v->next;
        free(tmp_v);
    }
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);

    free(ds);
    free(fp);
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    return 0;
}

ADIOS_GROUP * adios_read_datatap_gopen (ADIOS_FILE *fp, const char *grpname)
{
fprintf(stderr, "im here %s:%d\n", __FILE__,__LINE__);
    if(!grpname) {
        adios_error(err_invalid_group, "Group name is not valid");
        return NULL;
    }
    int i;
    for(i = 0; i < fp->groups_count; i ++) {   
        if(!strcmp(grpname, fp->group_namelist[i])) {
            return adios_read_datatap_gopen_byid(fp, i);
        }
    }
    adios_error(err_invalid_group, "Group %s is not valid", grpname);
    return NULL;
}

ADIOS_GROUP * adios_read_datatap_gopen_byid (ADIOS_FILE *fp, int grpid)
{
fprintf(stderr, "im here %s:%d\n", __FILE__,__LINE__);
    datatap_read_file_data *ds = (datatap_read_file_data *) fp->fh;
    ADIOS_GROUP * gp;

    adios_errno = 0;

    gp = (ADIOS_GROUP *) malloc(sizeof(ADIOS_GROUP));
    if (!gp) {
        adios_error(err_no_memory, "Could not allocate memory for group info");
        return NULL;
    }

    // TODO: again, assume one group per file
    gp->grpid = grpid;
    gp->gh = (uint64_t) 0; // TODO: we should re-organize the metadata
    gp->fp = fp;
    gp->attrs_count = 0; // attributes are not supported yet
    
fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
    // generate a list of variables with distinct names among all pgs
    gp->vars_count = ds->num_vars + ds->num_vars_peer;
fprintf(stderr, "im here count %d %s:%d\n", gp->vars_count,__FILE__,__LINE__);
    
    // to return a globally consistently ordered var list, we sort the list by var id
    datatap_var_info **vars = (datatap_var_info *) malloc(sizeof(datatap_var_info *) * gp->vars_count);
    if(!vars) {
        adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_gopen_byid()");
        return NULL;
    }
    memset(vars, 0, sizeof(datatap_var_info *) * gp->vars_count);

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
    datatap_var_info *current_var = ds->vars_peer;
    datatap_var_info *var_to_sort = current_var;
    while(var_to_sort) {
        int j;
        for(j = 0; j < gp->vars_count; j ++) {
            if(vars[j] == NULL) {
                vars[j] = var_to_sort;
                current_var = current_var->next;
                var_to_sort = current_var;
                break; 
            }
            else if(vars[j]->id < var_to_sort->id) {
                // move vars[j] 
                datatap_var_info *tmp; 
                tmp = vars[j];
                vars[j] = var_to_sort;
                var_to_sort = tmp;
                break;
            }
        }
    }
fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);

    current_var = ds->vars;
    var_to_sort = current_var;
    while(var_to_sort) {
        int j;
        for(j = 0; j < gp->vars_count; j ++) {
            if(vars[j] == NULL) {
                vars[j] = var_to_sort;
                current_var = current_var->next;
                var_to_sort = current_var;
                break;
            }
            else if(vars[j]->id < var_to_sort->id) {
                // move vars[j]
                datatap_var_info *tmp;
                tmp = vars[j];
                vars[j] = var_to_sort;
                var_to_sort = tmp;
                break;
            }
        }
    }
    // now the vars list is in descendent order
    gp->var_namelist = (char **) malloc(gp->vars_count * sizeof(char *));
    if(!gp->var_namelist) {
        adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_gopen_byid()");
        return NULL;
    }

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
    int i;
    for(i = 0; i < gp->vars_count; i ++) {
        int index = gp->vars_count - i - 1;
        gp->var_namelist[i] = (char *) malloc(strlen(vars[index]->varname) +
            strlen(vars[index]->varpath) + 2);
        if(!gp->var_namelist[i]) {
            adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_gopen_byid()");
            return NULL;
        }
        // return full path name
        // TODO: make sure the size of var_namelist[j] is enough
        if(!strcmp(vars[index]->varpath, "/")) {
            sprintf(gp->var_namelist[i], "/%s\0\0", vars[index]->varname);
        }
        else {
            sprintf(gp->var_namelist[i], "%s/%s\0", vars[index]->varpath,
                vars[index]->varname);
        }
    }
fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);

    // here we construct a bitmap of var list and send it to rank 0 so rank 0 knows 
    // which vars are only avaialble in its local memory
    memset(ds->var_bitmap, 0, VAR_BITMAP_SIZE);
    current_var = ds->vars_peer;
    while(current_var) {           
        int byte_pos = current_var->id / 8;
        int bit_pos = current_var->id % 8;
        unsigned char mask = 0x01 << bit_pos;
        ds->var_bitmap[byte_pos] = ds->var_bitmap[byte_pos] | mask;
        current_var = current_var->next;
    }

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
    // send back to rank 0 so it knows which var is missing on other processes
    unsigned char *combined_bitmap = NULL;
    int combined_size = ds->comm_size * VAR_BITMAP_SIZE;
    if(ds->my_rank == 0) {
        combined_bitmap = (unsigned char *) malloc(combined_size);
        if(!combined_size) {
            adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_gopen_byid()");
            return NULL;
        }
    }
    MPI_Gather(ds->var_bitmap, VAR_BITMAP_SIZE, MPI_BYTE, combined_bitmap, VAR_BITMAP_SIZE, MPI_BYTE, 0, ds->comm);

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
    if(ds->my_rank == 0) {
        // now determine which var is missing on other processes 
        for(i = 0; i < ds->comm_size; i ++) {
            int j;
            for(j = 0; j < VAR_BITMAP_SIZE; j ++) {
                ds->var_bitmap[j] = ds->var_bitmap[j] | combined_bitmap[i*VAR_BITMAP_SIZE+j];
            }
        }
        free(combined_bitmap);
        // if the corresponding bit is set in ds->var_bitmap, then the var is missing on other processes
    }

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
#if 0
    do {
        if(i == 0) { // the first pg
            gp->vars_count = ds->pgs[0].num_vars;
            gp->var_namelist = (char **) malloc(gp->vars_count * sizeof(char *));
            if(!gp->var_namelist) {
                adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_gopen_byid()");
                return NULL;
            }
            for(j = 0; j < gp->vars_count; j ++) {
                gp->var_namelist[j] = (char *) malloc(strlen(ds->pgs[i].vars[j].varname) + 
                    strlen(ds->pgs[i].vars[j].varpath) + 2);
                if(!gp->var_namelist[j]) {
                    adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_gopen_byid()");
                    return NULL;
                }   
                // return full path name   
                // TODO: make sure the size of var_namelist[j] is enough
                if(!strcmp(ds->pgs[i].vars[j].varpath, "/")) {
                    sprintf(gp->var_namelist[j], "/%s\0\0", ds->pgs[i].vars[j].varname);
                }
                else {
                    sprintf(gp->var_namelist[j], "%s/%s\0", ds->pgs[i].vars[j].varpath, 
                        ds->pgs[i].vars[j].varname);
                }    
            }
            k = gp->vars_count;
            i ++;
            continue;
        } 

        // go over the vars in the ith pg 
        for(j = 0; j < ds->pgs[i].num_vars; j ++) {
            // first, we need to make sure the var is not seen before
            int t;
            int is_new = 1;
            char fullname[256];
            if(!strcmp(ds->pgs[i].vars[j].varpath, "/")) {
                sprintf(fullname, "/%s\0", ds->pgs[i].vars[j].varname);
            }
            else { 
                sprintf(fullname, "%s/%s\0", ds->pgs[i].vars[j].varpath, ds->pgs[i].vars[j].varname);  
            } 
            for(t = 0; t < k; t++) {
                if(!strcmp(fullname, gp->var_namelist[t])) { 
                    is_new = 0;
                    break;  
                } 
            }

            if(!is_new) continue;

            // now add this var to list
            char **temp = gp->var_namelist;  
//            temp = (char **) realloc(gp->var_namelist, (k+1) * sizeof(char*));
            temp = (char **) malloc((k+1) * sizeof(char*));
            if(!temp) {
                adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_gopen_byid()");
                return NULL;
            }
            else {
                memcpy(temp, gp->var_namelist, k*sizeof(char *));
                free(gp->var_namelist);
                gp->var_namelist = temp;

            }
             
            // return full path name
            gp->var_namelist[k] = strdup(fullname);
            gp->vars_count ++; 
            k ++; 
        }
 
        i ++;
    }
    while(i < ds->num_pgs);
#endif 

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);

    gp->attr_namelist = 0;

    gp->timestep = ds->f_info->timestep;
    gp->lasttimestep = ds->f_info->timestep;

    // now we need to wait here until we see all data are ready
    pthread_mutex_lock(&(ds->f_info->mutex));
    while(queue_size(ds->f_info->dt_queue) != ds->f_info->num_chunks) {
        pthread_cond_wait(&(ds->f_info->cond2), &(ds->f_info->mutex));
    }
    pthread_mutex_unlock(&(ds->f_info->mutex));

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
    return gp;
}

int adios_read_datatap_gclose (ADIOS_GROUP *gp)
{
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
//    datatap_read_file_data *ds = (datatap_read_file_data *) gp->fp->fh;
    adios_errno = 0;
    free_namelist ((gp->attr_namelist),gp->attrs_count);
    int i;
    for(i = 0; i < gp->vars_count; i ++) {
        free(gp->var_namelist[i]); 
    }
    free(gp->var_namelist);
//    free_namelist ((gp->var_namelist),gp->vars_count);
    free(gp);
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    return 0;

}

int adios_read_datatap_get_attr (ADIOS_GROUP *gp, const char *attrname, 
                                 enum ADIOS_DATATYPES *type,
                                 int *size, void **data)
{
    // TODO: borrowed from dimes
    adios_error(err_invalid_read_method, "adios_read_datatap_get_attr is not implemented.");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}

int adios_read_datatap_get_attr_byid (ADIOS_GROUP *gp, int attrid, 
                                      enum ADIOS_DATATYPES *type, 
                                      int *size, void **data)
{
    // TODO: borrowed from dimes
    adios_error(err_invalid_read_method, "adios_read_datatap_get_attr_byid is not implemented.");
    *size = 0;
    *type = adios_unknown;
    *data = 0;
    return adios_errno;
}

ADIOS_VARINFO * adios_read_datatap_inq_var (ADIOS_GROUP *gp, const char *varname) 
{
    // TODO: usually user will read those variables reperesenting dimensions directly
//    error(err_invalid_read_method, "adios_read_datatap_inq_var is not implemented.");

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
    // find the var among all pgs
    ADIOS_VARINFO *v = (ADIOS_VARINFO *) malloc(sizeof(ADIOS_VARINFO));
    if(!v) {
        adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
        return NULL;
    }
    memset(v, 0, sizeof(ADIOS_VARINFO));

    datatap_read_file_data *ds = (datatap_read_file_data *) gp->fp->fh;
    
    int found = 0;
    datatap_var_info *current_var = ds->vars;
    while(current_var) {
        if(!compare_var_name(varname, current_var)) { 
            found = 1;
            break;
        }
        else {
            current_var = current_var->next;
        } 
    } 

    if(!found) {
        current_var = ds->vars_peer;
        while(current_var) {
            if(!compare_var_name(varname, current_var)) {
                found = 1;
                break;
            }
            else {
                current_var = current_var->next;
            }
        }
    }
  
    if(found) {
        v->grpid = gp->grpid;
        int i;
        for(i = 0; i < gp->vars_count; i ++) {
            if(!strcmp(gp->var_namelist[i], varname)) {
                v->varid = i; // TODO: this may not be cmpatible with BP
                break;
            }
        }
        v->type = current_var->type;
        v->ndim = current_var->ndims;
        v->timedim = current_var->time_dim;
        if(!v->ndim) { // scalar and string
            if(v->timedim != -1) { // scalar with time dimension
                v->ndim = 1;
                v->dims = (uint64_t *) malloc(sizeof(uint64_t));
                if(!v->dims) {
                    adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
                    return NULL;
                }
                v->dims[0] = 1; // TODO: only one timestep in the file
            }
            //int value_size = common_read_type_size(v->type, current_var->chunks->data);
            int value_size = current_var->data_size;
            v->value = malloc(value_size);
            if(!v->value) {
                adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
                return NULL;
            }
            memcpy(v->value, current_var->chunks->data, value_size);
        }
        else { // arrays
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
    else {
        adios_error(err_invalid_varname, "Cannot find var %s\n", varname);
        return NULL;
    }

    
#if 0
    int i, j;
    for(i = 0; i < ds->num_pgs; i ++) {
        for(j = 0; j < ds->pgs[i].num_vars; j ++) {
            // the parameter varname can be full path or just var name, so we
            // need to first find it by matching the name
            if(!compare_var_name(varname, &(ds->pgs[i].vars[j]))) {     
                v->grpid = gp->grpid;
                v->varid = j; // TODO: this may not be cmpatible with BP 
                v->type = ds->pgs[i].vars[j].type;
                v->ndim = ds->pgs[i].vars[j].ndims; 
                v->timedim = ds->pgs[i].vars[j].time_dim;
                if(!v->ndim) { // scalar
                    if(v->timedim != -1) { // scalar with time dimension
                        v->ndim = 1;
                        v->dims = (uint64_t *) malloc(sizeof(uint64_t));
                        if(!v->dims) {
                            adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
                            return NULL;
                        }
                        v->dims[0] = 1; // TODO: only one timestep in the file
                    }
                    int value_size = common_read_type_size(v->type, ds->pgs[i].vars[j].data);
                    v->value = malloc(value_size); 
                    if(!v->value) {
                        adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
                        return NULL;
                    }
                    memcpy(v->value, ds->pgs[i].vars[j].data, value_size);
                }
                else { // arrays  
                    v->dims = (uint64_t *) malloc(v->ndim * sizeof(uint64_t));   
                    if(!v->dims) {
                        adios_error(err_no_memory, "Cannot allocate buffer in adios_read_datatap_inq_var()");
                        return NULL;
                    }
                    int k;
                    for(k = 0; k < v->ndim; k ++) {
                        v->dims[k] = ds->pgs[i].vars[j].global_bounds[k]; 
                    }
                }
                return v;
            }
        }
    }

    return NULL;    
#endif

fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
}

ADIOS_VARINFO * adios_read_datatap_inq_var_byid (ADIOS_GROUP *gp, int varid)
{
    if(varid >= 0 && varid < gp->vars_count) {
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
        return adios_read_datatap_inq_var(gp, gp->var_namelist[varid]);
    }
    else {
        adios_error(err_invalid_varid, "Cannot find var %d\n", varid);
        return NULL;
    }
}

void adios_read_datatap_free_varinfo (ADIOS_VARINFO *vp)
{
    if(!vp) return;

    if(!vp->ndim) { // scalar
        if(vp->timedim != -1) { // scalar with time dimension
            free(vp->dims);
        }
        free(vp->value);
    }
    else { // arrays
        free(vp->dims);
    }
}

int64_t adios_read_datatap_read_var (ADIOS_GROUP *gp, const char *varname,
                                     const uint64_t *start, const uint64_t *count,
                                     void *data)
{
fprintf(stderr, "im here read var %s addr %p %s:%d\n", varname, data,__FILE__,__LINE__);
    int64_t total_size;
    datatap_read_file_data *ds = (datatap_read_file_data *) gp->fp->fh;
    int found = 0;

fprintf(stderr, "im here rank %d %s %s:%d\n", ds->my_rank,varname, __FILE__,__LINE__);
    datatap_var_info *current_var = ds->vars;
    while(current_var) {
fprintf(stderr, "im here rank %d %s %s %s %s:%d\n", ds->my_rank, varname, current_var->varname,current_var->varpath,__FILE__,__LINE__);
        if(!compare_var_name(varname, current_var)) {
fprintf(stderr, "im here rank %d %s %s %s %s:%d\n", ds->my_rank, varname, current_var->varname,current_var->varpath,__FILE__,__LINE__);
            // found it locally
            found = 1;
            if(!current_var->ndims) { // scalar
fprintf(stderr, "im here rank %d %s %s %s %s:%d\n", ds->my_rank, varname, current_var->varname,current_var->varpath,__FILE__,__LINE__);
                // TODO: check time dimension if there is any
                if(current_var->time_dim != -1 &&
                    (gp->fp->tidx_start != start[0] || count[0] > 1)) {
                // TODO: check time dimension if there is any
                    adios_error(err_no_data_at_timestep, "Specified time step is not available.");
                    return -1;
                }

                total_size = current_var->data_size;
                memcpy(data, current_var->chunks->data, total_size);
                return total_size;
            }
            else { // arrays
fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
                // TODO: check time dimension if there is any
                if(current_var->time_dim != -1) {
                    uint64_t ti = current_var->time_dim;
                    if(futils_is_called_from_fortran()) {
                        ti --;
                    }
                    // TODO: in Fortran index starts from 1 but in C index starts from 0
                    if(count[ti] > 1 || start[ti] != current_var->chunks->global_offsets[ti]) {
                        adios_error(err_no_data_at_timestep, "Specified time step is not available.");
                        return -1;
                    }
                }

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
                total_size = read_array(ds, current_var, start, count, data);
fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);

                // rank 0 should be responsible for send data out to peer readers
                if(ds->comm_size > 1) {
                    if(ds->my_rank == 0) {
                        // check if this array is missing on peer readers  
                        //if(current_var->num_chunks == 1) {
                        int byte_pos = current_var->id / 8;
                        int bit_pos = current_var->id % 8;
                        unsigned char mask = 0x01 << bit_pos;
                        //if(current_var->num_chunks == 1) {
                        if((ds->var_bitmap[byte_pos] & mask) != 0x00) {
fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
                            int rc = MPI_Bcast(data, total_size, MPI_BYTE, 0, ds->comm); 
                            if(rc != MPI_SUCCESS) {
                                fprintf(stderr, "rank %d: MPI_Bcast() returns error (%d). %s:%d\n",
                                    ds->my_rank, rc, __FILE__, __LINE__);
                                return -1;
                            }
                        }
                    } 
                }
                return total_size;
            }
        }
        else {
            current_var = current_var->next;
fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
        }
    }
fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);

    if(!found) {
        current_var = ds->vars_peer;
        while(current_var) {
            if(!compare_var_name(varname, current_var)) {
                // found it on remote peer reader
                if(!current_var->ndims) { // scalar
                    // TODO: check time dimension if there is any
                    if(current_var->time_dim != -1 &&
                        (gp->fp->tidx_start != start[0] || count[0] > 1)) {
                        adios_error(err_no_data_at_timestep, "Specified time step is not available.");
                        return -1;
                    }

                    total_size = current_var->data_size;
                    memcpy(data, current_var->chunks->data, total_size);
                    return total_size;
                }
                else { // arrays
                    // TODO: check time dimension if there is any
                    if(current_var->time_dim != -1) {
                        uint64_t ti = current_var->time_dim;
                        // TODO: in Fortran index starts from 1 but in C index starts from 0
                        if(futils_is_called_from_fortran()) {
                            ti --;
                        }

                        if(count[ti] > 1 || start[ti] != current_var->chunks->global_offsets[ti]) {
                            adios_error(err_no_data_at_timestep, "Specified time step is not available.");
                            return -1;
                        }
                    }

                    if(ds->my_rank != 0) {
                        int i = 0;
                        total_size = common_read_type_size(current_var->type, NULL); 
                        for(; i < current_var->ndims; i ++) {
                             total_size *= count[i];
                        }
                                     
                        if(ds->comm_size > 1) {
fprintf(stderr, "im here rank %d %s %ld %s:%d\n", ds->my_rank, current_var->varname,total_size, __FILE__,__LINE__);
                            // check if this array is missing on peer readers
                            int rc = MPI_Bcast(data, total_size, MPI_BYTE, 0, ds->comm);
fprintf(stderr, "im here rank %d %s %ld %s:%d\n", ds->my_rank, current_var->varname,total_size, __FILE__,__LINE__);
                            if(rc != MPI_SUCCESS) {
                                fprintf(stderr, "rank %d: MPI_Bcast() returns error (%d). %s:%d\n",
                                    ds->my_rank, rc, __FILE__, __LINE__);
                                return -1;
                            }
fprintf(stderr, "im here rank %d %s %ld %s:%d\n", ds->my_rank, current_var->varname,total_size, __FILE__,__LINE__);
                        }
                    }

fprintf(stderr, "im here read rank %d var %s addr %p total size %ld %s:%d\n", ds->my_rank, varname, data,total_size,__FILE__,__LINE__);
                    return total_size;
                }
            }
            else {
                current_var = current_var->next;
            }
        }
    }

fprintf(stderr, "im here rank %d %s:%d\n", ds->my_rank, __FILE__,__LINE__);
    adios_error(err_invalid_varname, "Cannot find var %s\n", varname);
    return -1;



#if 0
    // TODO: search through all pgs to find this var
    int p;
    for(p = 0; p < ds->num_pgs; p ++) {
        datatap_pg_info *current_pg = &(ds->pgs[p]);    
        datatap_var_info *var_info = NULL;
        int v;
        for(v = 0; v < current_pg->num_vars; v ++) {
            if(!compare_var_name(varname, &(current_pg->vars[v]))) {     
                // found it
                found = 1;

                if(!current_pg->vars[v].ndims) { // scalar
                    // TODO: check time dimension if there is any
                    if(current_pg->vars[v].time_dim != -1 && 
                        (gp->fp->tidx_start != start[0] || count[0] > 1)) {
                        adios_error(err_no_data_at_timestep, "Specified time step is not available.");
                        return -1;
                    }
  
                    total_size = common_read_type_size(current_pg->vars[v].type, current_pg->vars[v].data);                
                    memcpy(data, current_pg->vars[v].data, total_size);
                    return total_size;
                }
                else { // arrays
                    // TODO: check time dimension if there is any
                    if(current_pg->vars[v].time_dim != -1) {
                        // TODO: in Fortran index starts from 1 but in C index starts from 0 
                        uint64_t ti = current_pg->vars[v].time_dim;
                        if(count[ti] > 1 || start[ti] != current_pg->vars[v].global_offsets[ti]) {
                            adios_error(err_no_data_at_timestep, "Specified time step is not available.");
                            return -1;
                        } 
                    }

                    // Datatap batches per-variable reads, so here we only record the buffer address
                    datatap_var_info *var_info = (datatap_var_info *) malloc(sizeof(datatap_var_info));
                    if(!var_info) {
                        adios_error(err_no_memory, "Could not allocate memory for group info");
                        return -1;
                    }                
                    memcpy(var_info, &(current_pg->vars[v]), sizeof(datatap_var_info));
                    var_info->local_bounds = (uint64_t *) malloc(var_info->ndims * sizeof(uint64_t));
                    if(!var_info->local_bounds) {
                        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                        return NULL;
                    }
                    var_info->global_bounds = (uint64_t *) malloc(var_info->ndims * sizeof(uint64_t));
                    if(!var_info->global_bounds) {
                        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                        return NULL;
                    }
                    var_info->global_offsets = (uint64_t *) malloc(var_info->ndims * sizeof(uint64_t));
                    if(!var_info->global_offsets) {
                        adios_error(err_no_memory, "Cannot allocate memory for Datatap.");
                        return NULL;
                    }

                    total_size = 1;
                    int i;
                    for(i = 0; i < var_info->ndims; i ++) {
                        var_info->local_bounds[i] = count[i];
                        // TODO: it seems that user won't read from N-dimensional array into M-dimension
                        var_info->global_bounds[i] = current_pg->vars[v].global_bounds[i];  
                        var_info->global_offsets[i] = start[i];
                        total_size = total_size * count[i];
                    }
                    total_size *= common_read_type_size(var_info->type, NULL);
                    var_info->data = data;
                    var_info->data_size = total_size;
                    
                    // now we add this to the list of vars to read
                    if(!ds->vars_read) {
                        ds->vars_read = var_info;                
                        ds->vars_read_tail = var_info;
                        var_info->next = NULL;                   
                    }
                    else {
                        ds->vars_read_tail->next = var_info;     
                        ds->vars_read_tail = var_info;
                        var_info->next = NULL; 
                    }
                    ds->num_vars_read ++;                  
                    return total_size;
                }            
            }   
        }    
    }
fprintf(stderr, "im here read var %s addr %p %s:%d\n", varname, data,__FILE__,__LINE__);

    if(found) {
        return total_size;
    }

    adios_error(err_invalid_varname, "Cannot find var %s\n", varname);
    return -1;    
#endif 

}

int64_t adios_read_datatap_read_var_byid (ADIOS_GROUP *gp, int varid,
                                          const uint64_t *start,
                                          const uint64_t *count,
                                          void *data)
{
    if(varid >= 0 && varid < gp->vars_count) {
fprintf(stderr, "im here %s:%d\n",__FILE__,__LINE__);
        return adios_read_datatap_read_var (gp, gp->var_namelist[varid], start, count, data);
    }
    else {
        adios_error(err_invalid_varid, "Cannot find var %d\n", varid);
        return -1;
    }

}

int adios_read_datatap_get_dimension_order (const ADIOS_FILE *fp)
{
    return 0;
}

void adios_read_datatap_reset_dimension_order (ADIOS_FILE *fp, int is_fortran)
{
    // TODO
    adios_error(err_invalid_read_method, "adios_read_datatap_reset_dimension_order is not implemented.");
}
#endif

