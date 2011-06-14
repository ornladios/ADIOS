/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>   // nanosleep

// see if we have MPI or other tools
#include "config.h"

// xml parser
#include <mxml.h>

#include "adios.h"
//#include "adios_transport_hooks.h"
//#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "buffer.h"
#include "adios_error.h"

// include RDMA client functions
#include "ardma/ardma_client.h"

/*
#define adios_logger(verbose_level, ...) if (adios_staging_verbose >= verbose_level) fprintf (stderr, __VA_ARGS__); 

#define log_error(...) adios_logger(0, __VA_ARGS__)
#define log_warn(...) adios_logger(1, __VA_ARGS__)
#define log_info(...) adios_logger(2, __VA_ARGS__)
#define log_debug(...) adios_logger(3, __VA_ARGS__)
*/
static unsigned int adios_staging_verbose = 3;



static int adios_staging_initialized = 0;

/* we pack data into our private buffer (no shared buffer) */
struct adios_staging_buffer_struct 
{
    char * buf; 
    uint64_t size;
    uint64_t bytes_written;
    uint64_t offset;
};

// struct passed on between function calls from init to finalize
struct adios_staging_data_struct
{
    //uint64_t vars_start;
    //uint64_t vars_header_size;
    MPI_Comm group_comm;
    int rank;
    int size;

    // this method needs to have its own buffers
    struct adios_staging_buffer_struct pg; // holds the PG buffer to be sent out
    struct adios_staging_buffer_struct index; // holds the metadata index

    int memory_registered; // true if buffers have been already registered to DMA layer
    uint64_t memory_buffer_size; // the size of (to be) allocated memory. Should not use more than that

    struct ardma_client_connection *acc; // returned by ardma_client_connect()
    int timestep; // each close() will increase it to separate different requests
    char *path; // output file path
};


/* Functions and variables to save and set fd-> buffer to our private buffer and 
   then restore it */
    static struct adios_file_struct fd_save;  
    static void SWAP_FD_BUF(struct adios_file_struct *fd, struct adios_staging_buffer_struct *pg) 
    {
        fd_save.buffer = fd->buffer;
        fd_save.offset = fd->offset;
        fd_save.bytes_written = fd->bytes_written;
        fd_save.buffer_size = fd->buffer_size;
    
        fd->buffer = pg->buf;
        fd->offset = pg->offset;
        fd->bytes_written = pg->bytes_written;
        fd->buffer_size = pg->size;
    }

    static void UNSWAP_FD_BUF(struct adios_file_struct *fd, struct adios_staging_buffer_struct *pg) 
    {
        // sanity check: we must not allow realloc in buffer management
        if (pg->buf != fd->buffer) {
            log_error("ERROR in ADIOS STAGING method: buffer address has been changed.\n"
                "This is a bug. Please report it to the ADIOS team. UNSWAP_FD_BUF():\n"
                " fd (b=%x, size=%lld, offset=%lld, bytes_written=%lld\n"
                " pg (b=%x, size=%lld, offset=%lld, bytes_written=%lld\n",
                fd->buffer, fd->buffer_size, fd->offset, fd->bytes_written,
                pg->buf, pg->size, pg->offset, pg->bytes_written);
        }
        pg->offset = fd->offset;
        pg->bytes_written = fd->bytes_written;
        pg->size = fd->buffer_size;
    
        fd->buffer = fd_save.buffer;
        fd->offset = fd_save.offset;
        fd->bytes_written = fd_save.bytes_written;
        fd->buffer_size = fd_save.buffer_size;
    }


#ifdef HAVE_MPI
static void adios_var_to_comm (const char * comm_name
                              ,enum ADIOS_FLAG host_language_fortran
                              ,void * data
                              ,MPI_Comm * comm
                              )
{
    if (data)
    {
        int t = *(int *) data;

        if (!comm_name)
        {
            if (!t)
            {
                log_error ("communicator not provided and none "
                                 "listed in XML.  Defaulting to "
                                 "MPI_COMM_SELF\n"
                        );

                *comm = MPI_COMM_SELF;
            }
            else
            {
                if (host_language_fortran == adios_flag_yes)
                {
                    *comm = MPI_Comm_f2c (t);
                }
                else
                {
                    *comm = *(MPI_Comm *) data;
                }
            }
        }
        else
        {
            if (!strcmp (comm_name, ""))
            {
                if (!t)
                {
                    log_error ("communicator not provided and none "
                                     "listed in XML.  Defaulting to "
                                     "MPI_COMM_SELF\n"
                            );

                    *comm = MPI_COMM_SELF;
                }
                else
                {
                    if (host_language_fortran == adios_flag_yes)
                    {
                        *comm = MPI_Comm_f2c (t);
                    }
                    else
                    {
                        *comm = *(MPI_Comm *) data;
                    }
                }
            }
            else
            {
                if (!t)
                {
                    log_error ("communicator not provided but one "
                                     "listed in XML.  Defaulting to "
                                     "MPI_COMM_WORLD\n"
                            );

                    *comm = MPI_COMM_WORLD;
                }
                else
                {
                    if (host_language_fortran == adios_flag_yes)
                    {
                        *comm = MPI_Comm_f2c (t);
                    }
                    else
                    {
                        *comm = *(MPI_Comm *) data;
                    }
                }
            }
        }
    }
    else
    {
        log_error ("coordination-communication not provided. "
                         "Using MPI_COMM_SELF instead\n"
                );

        *comm = MPI_COMM_SELF;
    }
}
#endif

void adios_staging_init (const char * parameters ,struct adios_method_struct * method
                      )
{
    struct adios_staging_data_struct * md = 0;

    if (adios_staging_initialized)
    {
        fprintf (stderr, "ADIOS STAGING method has been already initialized\n");
        return;
    }

    method->method_data = malloc (sizeof (struct adios_staging_data_struct));
    md = (struct adios_staging_data_struct *) method->method_data;
    memset (md, 0, sizeof(struct adios_staging_data_struct));
    md->group_comm = MPI_COMM_NULL;
    /*
    md->rank = 0;
    md->size = 0;
    md->memory_buffer_size = 0;
    md->memory_registered = 0;
    md->pg.buf = NULL;
    md->pg.size = 0;
    md->pg.offset = 0;
    md->pg.bytes_written = 0;
    md->acc = NULL;
    md->timestep = 0;
    */

    // parse the parameters into key=value segments for optional settings
    if (parameters && strlen(parameters)>0)
    {
        char * p = strdup (parameters);
        char * token = strtok (p, ";");

        while (token)
        {
            char * equal_pos = strchr (token, '=');
            if (!equal_pos) {
                fprintf(stderr, "STAGING parameter '%s' is not recognized\n", token);
                continue;
            }
            int equal = equal_pos - token + 1;
            int len = strlen (token);
            char * key = malloc (len + 1);
            char * value = malloc (len + 1);
            strncpy (key, token, equal);
            key [equal - 1] = 0;
            strncpy (value, equal_pos + 1, len - equal);
            value [len - equal] = 0;

            if (key && value) {
                if (strcasecmp (key, "memory_buffer_size") == 0) {
                    long long v = atoll (value);
                    md->memory_buffer_size = (uint64_t)v;
                } else if (strcasecmp (key, "verbose") == 0) {
                    int v = atoi (value);
                    adios_staging_verbose = (unsigned int)v;
                }
                else {
                    fprintf (stderr, "STAGING parameter: key: {%s} "
                            "value: {%s} not recognized. Ignored\n"
                            ,key, value);
                }
            }

            if (key) { free (key); key = 0; }
            if (value) { free (value); value = 0; }
            token = strtok (NULL, ";");
        } 

        if (p) { free (p); p = 0; }

    } /* if parameters */

    if (ardma_client_init(adios_staging_verbose)) {
        log_error ("ADIOS STAGING method failed to initialize the RDMA layer\n");
        return;
    }

    adios_staging_initialized = 1;
}


static enum ADIOS_ERRCODES previous_staging_status (struct adios_staging_data_struct * md) {
    
    // this function works and returns ARDMA_STAGING_READY if connection has failed above
    enum ARDMA_STAGING_STATUS s = ardma_client_staging_status(md->acc);
    int c = (int)s;
    int allcompleted;
    MPI_Allreduce (&c, &allcompleted, 1, MPI_INT, MPI_MAX, md->group_comm);
    log_debug ("%s, rank=%d, acc=%lld, c=%d c_all=%d\n", __func__, md->rank, md->acc, c, allcompleted);
    if ((enum ARDMA_STAGING_STATUS)allcompleted == ARDMA_STAGING_INPROGRESS) 
        return err_staging_in_progress;
    else if ((enum ARDMA_STAGING_STATUS)allcompleted == ARDMA_STAGING_READY) 
        return err_no_error;
    return err_staging_failed;
}

/*enum ADIOS_ERRCODES should be for all method's open() */
int adios_staging_open (struct adios_file_struct * fd
        ,struct adios_method_struct * method, void * comm
        )
{
    struct adios_staging_data_struct * md = 
        (struct adios_staging_data_struct *) method->method_data;

    if (!adios_staging_initialized)
    {
        log_error ("%s, method was not initialized\n", __func__);
        return err_staging_uninitialized;
    }

#ifdef HAVE_MPI
    // Need to figure out how many processes will send data to the staging area
    adios_var_to_comm (fd->group->group_comm
                      ,fd->group->adios_host_language_fortran
                      ,comm
                      ,&md->group_comm
                      );

    if (md->group_comm == MPI_COMM_NULL)
    {
        md->group_comm = MPI_COMM_SELF;
    }

    MPI_Comm_rank (md->group_comm, &md->rank);
    MPI_Comm_size (md->group_comm, &md->size);
#endif

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            log_error ("ADIOS STAGING method does not support reading\n");
            return err_invalid_method_mode;
        }

        case adios_mode_write:
        case adios_mode_append:
        {
            
            /* Connect to staging server, in the very first open */
            if (!md->acc) {
                log_debug ("%s, rank=%d, try to connect to server\n", __func__, md->rank);
                md->acc = ardma_client_connect (md->group_comm);
                MPI_Barrier (md->group_comm);
                /* Note: the ARDMA layer assumes we are not sending messages until everybody
                   has connected */
            }
            // we do not return error here if connection failed, because below we have
            // an allreduce, which needs every process

            /* Return error if previous staging has not yet completed for all processes */
            // this function works and returns err_no_error if connection has failed above
            enum ADIOS_ERRCODES aerr = previous_staging_status(md);
            if (aerr != err_no_error)
                return aerr;

            // if connection failed on any process disconnect and return error now
            int c = (md->acc != NULL);
            int allconnected;
            MPI_Allreduce (&c, &allconnected, 1, MPI_INT, MPI_MIN, md->group_comm);
            if (!allconnected) {
                if (!md->acc) {
                    ardma_client_disconnect(&md->acc);
                }
                return err_connection_failed;
            }

            fd->base_offset = 0;
            fd->pg_start_in_file = 0;
            break;
        }

        default:
        {
            log_error ("rank %d Error: Unknown file mode: %d\n", md->rank, fd->mode);
            return err_invalid_method_mode;
        }
    }

    /* determine file path now */
    md->path = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (md->path, "%s%s", method->base_path, fd->name);
    log_debug("rank %d: output file len = %d, path = %s\n", md->rank, strlen (method->base_path) + strlen (fd->name), md->path);

    return err_no_error;
}

enum ADIOS_FLAG adios_staging_should_buffer (struct adios_file_struct * fd
                                          ,struct adios_method_struct * method
                                          )
{
    struct adios_staging_data_struct * md = 
           (struct adios_staging_data_struct *) method->method_data;

    if (!adios_staging_initialized) 
        return fd->shared_buffer;
    if (!md->acc) 
        return fd->shared_buffer;

    /* 
      This method needs to handle memory itself, because it needs to lock and register
      the buffer to the DMA layer, which should be done infrequently. 
      ADIOS shared buffer use malloc/free in each open/close cycle. 
    */

    // initialize for this round of open/write/close
    md->pg.size = 0; 

    if (fd->mode == adios_mode_read) 
    {
        log_error ("ERROR: ADIOS STAGING method does not support reading\n");
        return fd->shared_buffer;
    }

    /* we allocate the buffer for ourselves and lock/register to DMA layer */
    /* NOTE: above in common_adios_group_size, the size will be given back to ADIOS
       since we return adios_flag_no. So we use extra memory not what allowed 
       for ADIOS in the XML file */
    if (md->memory_buffer_size == 0) {
        /* no XML option and first open: use PG size as max memory buffer */
        md->memory_buffer_size = fd->write_size_bytes;
    }

    /* we allocate the metadata index buffer for ourselves too.
       Ask adios if we are allowed to use more memory.*/
    uint64_t index_size = 1024*1024; // FIXME: lets allocate 1MB for index

    if (!md->memory_registered || md->memory_buffer_size < fd->write_size_bytes) {

        if (md->memory_buffer_size < fd->write_size_bytes) {
            log_warn (
                    "WARNING: ADIOS STAGING method: your requested buffer size (%lld) is more "
                    "than allocated (%lld).\n Reallocating memory for RDMA is costly, so use "
                    "larger parameter 'memory_buffer_size' for this method.\n"
                    ,fd->write_size_bytes, md->memory_buffer_size);
            return adios_flag_no;  
        }

        /* allocate and register memory to RDMA layer at 
           first time and if current buffer is not big enough */
        md->pg.buf = (char *) ardma_client_register_pg_memory (md->acc, 
                                                               md->memory_buffer_size, 
                                                               md->pg.buf);

        if (!md->pg.buf) {
            log_error (
                    "ERROR: ADIOS STAGING method: could not allocate RDMA buffer size %lld.\n"
                    "Will not output anything\n"
                    ,md->memory_buffer_size);
            return adios_flag_no;  
        }

        if (!md->memory_registered) { 
            /* register index memory to RDMA layer only at first time*/
            md->index.buf = (char *) ardma_client_register_index_memory (md->acc, 
                                                                      index_size, 
                                                                      NULL);

            if (!md->index.buf) {
                log_error (
                        "ERROR: ADIOS STAGING method: could not allocate buffer size %lld\n"
                        "for metadata. Will not output anything\n"
                        ,md->memory_buffer_size);
                ardma_client_deregister_pg_memory (md->acc, md->pg.buf);
                md->pg.buf = NULL;
                md->pg.size = 0;
                return adios_flag_no;  
            }
        }
        md->memory_registered = 1;

    }

    /* initialize PG buffer */
    md->pg.size = fd->write_size_bytes;
    memset(md->pg.buf, 0, md->pg.size);
    md->pg.offset = 0;
    md->pg.bytes_written = 0;


    md->index.size = index_size;
    md->index.offset = 0;
    md->index.bytes_written = 0;

    /* FIXME: lock and register md->pg.buf and md->pg.index in DMA */

    md->memory_registered = 1;

    /* We build the PG and metadata index the same way as common_adios does in 
       case of a shared buffer. We need to put our buffer into fd for each 
       operation and then restore the original fd */

    // write the process group header
    SWAP_FD_BUF(fd, &md->pg);
    adios_write_process_group_header_v1 (fd, md->pg.size);
    // setup for writing vars
    adios_write_open_vars_v1 (fd);
    UNSWAP_FD_BUF(fd, &md->pg); 


    /* Prevent ADIOS from allocating a shared buffer */
    return adios_flag_no;  
}

void adios_staging_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,void * data
                       ,struct adios_method_struct * method
                       )
{
    struct adios_staging_data_struct * md = 
           (struct adios_staging_data_struct *) method->method_data;
    if (!adios_staging_initialized) return;
    if (!md->acc) return; // no connection, so do nothing
    if (md->pg.size == 0) return; // no memory, so do nothing

    // Encode this variable to PG
    SWAP_FD_BUF(fd, &md->pg);
    // var payload sent for sizing information
    adios_write_var_header_v1 (fd, v);
    // write payload
    adios_write_var_payload_v1 (fd, v);
    UNSWAP_FD_BUF(fd, &md->pg); 

    log_debug("%s: rank=%d var=%s, b=%x, size=%lld, offset=%lld, bytes_written=%lld\n",
            __func__, md->rank, v->name, md->pg.buf, md->pg.size, md->pg.offset, md->pg.bytes_written);

}


void adios_staging_close (struct adios_file_struct * fd
                       ,struct adios_method_struct * method
                       )
{
    struct adios_staging_data_struct * md = 
            (struct adios_staging_data_struct *) method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_var_struct_v1 * vars_root = 0;
    struct adios_index_attribute_struct_v1 * attrs_root = 0;

    uint64_t index_start = 0; //fd->base_offset + fd->offset;

    if (!adios_staging_initialized) return;
    if (!md->acc) return; // no connection, so do nothing
    if (md->pg.size == 0) return; // no memory, so do nothing

    
    switch (fd->mode)
    {
        case adios_mode_write:
        case adios_mode_append:
        {
            // Finish the PG with the attributes
            SWAP_FD_BUF(fd, &md->pg);
            adios_write_close_vars_v1 (fd);
            adios_write_open_attributes_v1 (fd);

            while (a)
            { 
                adios_write_attribute_v1 (fd, a);
                a = a->next;
            }
            adios_write_close_attributes_v1 (fd);
            UNSWAP_FD_BUF(fd, &md->pg);

            log_debug("%s: attributes added, b=%x, size=%lld, offset=%lld, bytes_written=%lld\n",
                    __func__, md->pg.buf, md->pg.size, md->pg.offset, md->pg.bytes_written);

            // build index
            adios_build_index_v1 (fd, &pg_root, &vars_root,&attrs_root);
            // write index into buffer
            char * buf_save = md->index.buf; // check that it does not change
            adios_write_index_v1 (&md->index.buf, &md->index.size, &md->index.offset
                                 ,index_start, pg_root, vars_root,attrs_root);

            if (md->index.buf != buf_save) {
                log_error (
                        "ERROR in ADIOS STAGING method: index buffer address has been changed.\n"
                        "This is a bug. Please report it to the ADIOS team. adios_staging_close():\n"
                        " buf_save=%x, new buf=%x size=%lld, offset=%lld\n",
                        buf_save, md->index.buf, md->index.size, md->index.offset);
            }

            log_debug("%s: index buffer b=%x, size=%lld, offset=%lld, bytes_written=%lld\n",
                    __func__, md->index.buf, md->index.size, md->index.offset, md->index.bytes_written);

            // FIXME: Here we should send the request to staging server with DMA
            // send out md->pg.buf and md->index.buf
            ardma_client_send_request (md->acc, 
                                       md->pg.offset /*actual size*/, 
                                       md->index.offset, 
                                       md->path,
                                       md->timestep);


            adios_clear_index_v1 (pg_root, vars_root, attrs_root);
            break;
        }


        case adios_mode_read:
        {
            log_error ("ERROR: ADIOS STAGING does not support read mode\n");
            return;
        }

        default:
        {
            log_error ("Unknown file mode: %d\n", fd->mode);
            return;
        }
    }
    
    free(md->path);
    md->timestep++;

}

void adios_staging_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_staging_data_struct * md = 
            (struct adios_staging_data_struct *) method->method_data;

    struct timespec treq = {.tv_sec=1, .tv_nsec=0}; // 1s sleeptime
    struct timespec trem;
    /* Block until last staging is completed, otherwise the server
       may not pull off all data before we exit */
    if (md->acc) {
        enum ADIOS_ERRCODES aerr = previous_staging_status(md);
        if (aerr == err_staging_in_progress && md->rank == 0 ) 
            log_warn ("adios_finalize: Staging is still in progress. Wait until finished...\n");
        while (aerr == err_staging_in_progress) {
            nanosleep(&treq, &trem);
            aerr = previous_staging_status(md);
        }
    }

    /* Deregister from DMA and free the private buffer now */
    if (md->pg.size > 0) {
        ardma_client_deregister_pg_memory (md->acc, md->pg.buf);
        md->pg.size = 0;
        md->pg.offset = 0;
    }
    if (md->index.size > 0) {
        ardma_client_deregister_index_memory (md->acc, md->index.buf);
        md->index.size = 0;
        md->index.offset = 0;
    }
    if (adios_staging_initialized)
        adios_staging_initialized = 0;

    /* FIXME: Disconnect from staging server */
    if (md->acc) {
        ardma_client_disconnect(&md->acc);
        md->acc = NULL;
    }
}

void adios_staging_end_iteration (struct adios_method_struct * method)
{
}

void adios_staging_start_calculation (struct adios_method_struct * method)
{
}

void adios_staging_stop_calculation (struct adios_method_struct * method)
{
}

void adios_staging_get_write_buffer (struct adios_file_struct * fd
                                  ,struct adios_var_struct * v
                                  ,uint64_t * size
                                  ,void ** buffer
                                  ,struct adios_method_struct * method
                                  )
{
    uint64_t mem_allowed;

    log_error ("ERROR: ADIOS STAGING method does not support adios_get_write_buffer()\n");
    *buffer = 0;
}

void adios_staging_read (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * buffer
                      ,uint64_t buffer_size
                      ,struct adios_method_struct * method
                      )
{
    log_error ("ERROR: ADIOS STAGING method does not support adios_read()\n");
}

