#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>

#include "config.h"

#include "bpindexing_thread.h"
#include "globals.h"
#include "precedence.h"

#include "../adios_internals.h"


/* bpindexing thread is started from bpworker thread
   In a non-threaded environment, bpindexing_doindex() is called from bpworker.
*/

struct indexing_data_struct indexing_data;
static bool indexing_inited = false;

#if HAVE_PTHREAD
void bpindexing_doindex(void);
void * bpindexing_thread_main (void *arg)
{   
    bool should_exit = false;
    pthread_mutex_init (&indexing_mutex, NULL);
    pthread_cond_init  (&indexing_cv, NULL);

    while (!gd.terminate && !should_exit) {
        //log_debug("rank %d: indexing thread: block on condition variable.\n", gd.mpi_rank);
        pthread_mutex_lock (&indexing_mutex);
        pthread_cond_wait(&indexing_cv, &indexing_mutex);
        //log_debug("rank %d: indexing thread: woke up from blocking, request=%d.\n", 
        //          gd.mpi_rank, indexing_data.indexing_request);
        switch (indexing_data.indexing_request) {
            case INDEXING_REQUEST_INDEX:
                log_debug("rank %d: indexing thread: worker asked to perform indexing.\n", gd.mpi_rank);
                bpindexing_doindex();
                log_debug("rank %d: indexing thread: indexing completed.\n", gd.mpi_rank);
                break;
            case INDEXING_REQUEST_EXIT:
                //log_debug("rank %d: indexing thread: worker asked to exit.\n", gd.mpi_rank);
                should_exit = true;
                break;
            default:
                log_error("rank %d: indexing thread: Unknown request from worker.\n", gd.mpi_rank);
                break;
        }
        pthread_mutex_unlock(&indexing_mutex);
    }

    if (gd.terminate) {
        log_debug("rank %d: indexing thread: exit due to failure.\n", gd.mpi_rank);
    } else {
        log_debug("rank %d: indexing thread: exit due to request.\n", gd.mpi_rank);
    }

    pthread_mutex_destroy (&indexing_mutex);
    pthread_cond_destroy  (&indexing_cv);
    pthread_exit(NULL);
    return NULL; // just to avoid compiler warning
}
#endif

static void adjust_index (uint64_t offset_to_add, int subfile_index,
                          struct adios_index_var_struct_v1 * vars_root,
                          struct adios_index_attribute_struct_v1 * attrs_root
                         )
{
    while (vars_root) {   
        vars_root->characteristics [0].offset += offset_to_add;
        vars_root->characteristics [0].payload_offset += offset_to_add;
        vars_root->characteristics [0].file_index = subfile_index;
        vars_root = vars_root->next;
    }
    while (attrs_root) {
        attrs_root->characteristics [0].offset += offset_to_add;
        attrs_root->characteristics [0].payload_offset += offset_to_add;
        attrs_root->characteristics [0].file_index = subfile_index;
        attrs_root = attrs_root->next;
    }
}   


void bpindexing_doindex (void)
{
    struct adios_bp_buffer_struct_v1 b;
    struct adios_index_process_group_struct_v1 *pg_root,    *single_pg_root;
    struct adios_index_var_struct_v1           *vars_root,  *single_vars_root;
    struct adios_index_attribute_struct_v1     *attrs_root, *single_attrs_root;

    int i;
    uint64_t file_offset = 0;
    struct globals_client_data * cd;

    pg_root = NULL;
    vars_root = NULL;
    attrs_root = NULL;

    /* 
        STEP 1: Parse each client's index and merge into one large index 
    */
    for (i = 0; i < gd.nc; i++)
    {
        cd = &gd.clientdata[gd.order[i]];
        b.buff = gd.rdma_buffer + cd->idx_offset;
        b.length = cd->idx_size;
        b.offset = 0;    
        single_pg_root = NULL;
        single_vars_root = NULL;
        single_attrs_root = NULL;

        /* build lists from idx buffer received from client */
        adios_parse_process_group_index_v1 (&b, &single_pg_root);
        adios_parse_vars_index_v1 (&b, &single_vars_root);
        adios_parse_attributes_index_v1 (&b, &single_attrs_root);

        /* adjust the offset from 0 the client assumed to the
           actual file offset + the subfile index to server rank */
        adjust_index (file_offset, gd.mpi_rank, single_vars_root, single_attrs_root);

        // Note: merge also frees elements in single_* lists
        adios_merge_index_v1 (&pg_root,&vars_root,&attrs_root,
                single_pg_root, single_vars_root, single_attrs_root);

        // we calculate next file offset here too, not wait on worker to
        // calculate it later in the pull cycle
        file_offset += cd->pg_size; 
    }

    /* Write index into a buffer for output */
    char * buffer = 0;
    uint64_t buffer_size = gd.index_size+4;
    uint64_t buffer_offset = 0;
    uint32_t flag = 0;
    buffer = malloc (buffer_size);     // index size ~ total index size of all clients
                                       // +4: BP version at the end
                                       // Note: it is freed after writing in writer thread
    
    if (buffer) {
        // Note: this function increases the buffer if necessary
        // Note: buffer_size >= filled area, buffer_offset = filled area
        adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset,
                              file_offset, pg_root, vars_root, attrs_root);
        log_debug("rank %d: index buffer size=%d, pg_idx_offs=%lld var_idx_offs=%lld "
                  "attr_idx_offs=%lld\n", gd.mpi_rank, buffer_offset, 
                  *(uint64_t*)&(buffer[buffer_offset-24]),
                  *(uint64_t*)&(buffer[buffer_offset-16]),
                  *(uint64_t*)&(buffer[buffer_offset-8])
                  );

        adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);
    } else {
        log_error("rank %d: Indexing ERROR: Cannot allocate memory for output "
                  "index buffer of size %lld\n", gd.mpi_rank, gd.index_size);
    }

    // easy half done: local index built
    indexing_data.laddr = (uint64_t)buffer;
    indexing_data.lsize = (uint64_t)buffer_offset;
    log_debug("rank %d: Indexing step 1: local index built\n", gd.mpi_rank);



    /* 
        STEP 2: Gather and build global index on rank 0 
        Assumption: total size of indexes is still <2GB because of
                    int arrays in MPI_Gatherv
    */
    precedence_get(gd.prec_writer_indexing, 1);
    MPI_Comm newcomm;
    MPI_Comm_dup(gd.mpi_comm, &newcomm);
    
    if (gd.mpi_rank > 0) {
        // just send local index
        int bufsize32 = (int) buffer_offset; // just deal with 32bit (<2GB)
        // if index is >2GB, we are already screwed, so assumption is feasible
        MPI_Gather (&bufsize32, 1, MPI_INT, 0, 0, MPI_INT,0, newcomm);
        MPI_Gatherv (buffer, bufsize32, MPI_BYTE,0, 0, 0, MPI_BYTE,0, newcomm);
        log_debug("rank %d: Indexing: done with mpi gatherv\n", gd.mpi_rank);

    } else {

        /* gather index sizes first, then all indices */
        int * index_sizes = malloc (sizeof(int) * gd.mpi_size);
        int * index_offsets = malloc (sizeof(int) * gd.mpi_size);
        char * recv_buffer = 0;
        uint32_t size = 0; 
        uint32_t total_size = 0; // sum of size of index from rank 1..gd.mpi_size-1

        MPI_Gather (&size, 1, MPI_INT, index_sizes, 1, MPI_INT, 0, newcomm);

        for (i = 0; i < gd.mpi_size; i++) {                
            index_offsets [i] = total_size;
            total_size += index_sizes [i];
        }                
        // Note, index_offsets[0] = index_offsets[1] = 0 since index_sizes[0] = 0
        // rank 1's index will start from 0 offset in gathered buffer

        recv_buffer = malloc (total_size);
        if (recv_buffer) {
            log_debug("rank %d: Indexing step 2: gather local indices, total size=%d\n", 
                      gd.mpi_rank, total_size);

            MPI_Gatherv (&size, 0, MPI_BYTE, 
                         recv_buffer, index_sizes, index_offsets,
                         MPI_BYTE, 0, newcomm
                        );

            /* parse each index and merge it into rank 0's index */
            for (i = 1; i < gd.mpi_size; i++) {                
                b.buff = recv_buffer + index_offsets [i];
                b.length = index_sizes[i];
                b.offset = 0;    
                single_pg_root = NULL;
                single_vars_root = NULL;
                single_attrs_root = NULL;

                /* build lists from idx buffer received from server i */
                adios_parse_process_group_index_v1 (&b, &single_pg_root);
                adios_parse_vars_index_v1 (&b, &single_vars_root);
                adios_parse_attributes_index_v1 (&b, &single_attrs_root);

                // Note: merge also frees elements in single_* lists
                adios_merge_index_v1 (&pg_root,&vars_root,&attrs_root,
                        single_pg_root, single_vars_root, single_attrs_root);
            } 
            free (recv_buffer);

            /* Write global index into a buffer for output */
            char * gbuffer = 0;
            uint64_t gbuffer_size = buffer_offset + total_size;
            uint64_t gbuffer_offset = 0;
            gbuffer = malloc (gbuffer_size); // Note: it is freed after writing in writer thread

            if (gbuffer) {
                // Note: this function increases the buffer if necessary
                // Note: buffer_size >= filled area, buffer_offset = filled area
                adios_write_index_v1 (&gbuffer, &gbuffer_size, &gbuffer_offset,
                                      0, pg_root, vars_root, attrs_root);
                log_debug("rank %d: global index buffer size=%d, pg_idx_offs=%lld "
                          "var_idx_offs=%lld attr_idx_offs=%lld\n", 
                          gd.mpi_rank, gbuffer_offset, 
                          *(uint64_t*)&(gbuffer[gbuffer_offset-24]),
                          *(uint64_t*)&(gbuffer[gbuffer_offset-16]),
                          *(uint64_t*)&(gbuffer[gbuffer_offset-8])
                        );

                flag |= ADIOS_VERSION_HAVE_SUBFILE;
                adios_write_version_flag_v1 (&gbuffer, &gbuffer_size, &gbuffer_offset, flag);
            } else {
                log_error("rank %d: Indexing ERROR: Cannot allocate memory for output "
                        "global index buffer of size %lld\n", gd.mpi_rank, gbuffer_size);
            }

            // second half done: global index built
            indexing_data.gaddr = (uint64_t)gbuffer;
            indexing_data.gsize = (uint64_t)gbuffer_offset;
            log_debug("rank %d: Indexing step 2: global index built\n", gd.mpi_rank);

        } else {
            log_error("rank %d: Indexing ERROR: Could not allocate buffer with "
                      "size %d to receive metadata from all staging servers\n", 
                      gd.mpi_rank, total_size);
            indexing_data.gaddr = 0;
            indexing_data.gsize = 0;
        }

        free (index_sizes);
        free (index_offsets);
    }
    precedence_release(gd.prec_writer_indexing);
    
    log_debug("rank %d: Indexing: clear and complete doindex\n", gd.mpi_rank);
    adios_clear_index_v1 (pg_root, vars_root, attrs_root);
    indexing_data.indexing_completed = true;
}

