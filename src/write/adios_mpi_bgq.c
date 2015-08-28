/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <spi/include/kernel/location.h>
#include <spi/include/kernel/process.h>
#include <firmware/include/personality.h>

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

#include "public/adios_error.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"
#include "core/adios_logger.h"

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#   include "core/adios_timing.h"
#endif


enum ADIOS_MPI_BGQ_IO_TYPE
{
    ADIOS_MPI_BGQ_IO_NONE   = 0,
    ADIOS_MPI_BGQ_IO_AG     = 1, // simple all to one aggregation
    ADIOS_MPI_BGQ_IO_BG     = 2, // Brigade aggregation
    ADIOS_MPI_BGQ_IO_SIMPLE = 3, // simple write out multiple files and each group
                                  // write out a shared file.
};

static int adios_mpi_bgq_initialized = 0;

#define is_aggregator(rank)  md->g_is_aggregator[rank]
#define FREE(v) \
  if (v)        \
  {             \
      free (v); \
      v = 0;    \
  }             \

#define MAX_AGG_BUF      704643072
// Temporarily change BLOCK_UNIT to 1, which disable write aligment.
// GPFS doesn't seem to support holes in a file
// in the same way as Lustre. Setting BLOCK_UNIT to 8 MiB
// which is GPFS block size, will inflate the file size and
// that concerns large runs. Q. Liu, 10/19/2014
#define BLOCK_UNIT     1

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#   define START_TIMER(t) adios_timing_go (fd->group->timing_obj, (t) ) 
#else
#   define START_TIMER(t) ; 
#endif

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#   define STOP_TIMER(t) adios_timing_stop (fd->group->timing_obj, (t) )
#else
#   define STOP_TIMER(t) ;
#endif

struct adios_MPI_data_struct
{
    MPI_File fh;
    MPI_File mfh;
    char * subfile_name;
    MPI_Request req;
    MPI_Status status;
    MPI_Comm group_comm;
    int rank;
    int size;
    MPI_Comm file_comm;
    int file_comm_rank;
    int file_comm_size;

    struct adios_bp_buffer_struct_v1 b;
    struct adios_index_struct_v1 * index;

    int * g_is_aggregator;
    int g_num_aggregators;
    int g_color2;
    int partition_id;
    int n_partitions;
    MPI_Offset * g_offsets;
    struct adios_MPI_thread_data_open * open_thread_data;
    enum ADIOS_MPI_BGQ_IO_TYPE g_io_type;
    int g_have_mdf;
};

struct adios_MPI_thread_data_open
{
    struct adios_MPI_data_struct * md;
    char * parameters;
};

struct adios_MPI_thread_data_write
{
    MPI_File * fh;
    uint64_t * base_offset;
    void * aggr_buff;
    int * total_data_size;
};

#if defined(__APPLE__)
#       include <sys/param.h>
#       include <sys/mount.h>
#else
#       include <sys/statfs.h>
#endif

static void trim_spaces (char * str)
{
    char * t = str, * p = NULL;
    while (*t != '\0')
    {
        if (*t == ' ')
        {
            p = t + 1;
            strcpy (t, p);
        }
        else
            t++;
    }

}

void adios_mpi_bgq_add_offset (uint64_t var_offset_to_add
                              ,uint64_t attr_offset_to_add
                              ,struct adios_index_var_struct_v1 * vars_root
                              ,struct adios_index_attribute_struct_v1 * attrs_root
                              )
{
    while (vars_root)
    {
        vars_root->characteristics [0].offset += var_offset_to_add;
        vars_root->characteristics [0].payload_offset += var_offset_to_add;
        vars_root = vars_root->next;
    }

    while (attrs_root)
    {
        attrs_root->characteristics [0].offset += attr_offset_to_add;
        attrs_root->characteristics [0].payload_offset += attr_offset_to_add;
        attrs_root = attrs_root->next;
    }
}

static void
adios_mpi_bgq_set_aggregation_parameters(char * parameters, struct adios_MPI_data_struct * md)
{
    int err = 0, flag, i, aggr_group_size, remain, index;
    int nproc = md->size, rank = md->rank;
    char value[64], *temp_string, *p_count,*p_size;

    temp_string = (char *) malloc (strlen (parameters) + 1);

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "aggregation_type"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");

        if (!q)
            md->g_io_type = atoi (q + 1);
        else
            md->g_io_type = atoi (p + 1);
    }
    else
    {
        // by default, use BG. Currently only simple type
        // is supported. Q. Liu, 11-29-2013.
        md->g_io_type = ADIOS_MPI_BGQ_IO_SIMPLE;
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "have_metadata_file"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");

        if (!q)
            md->g_have_mdf = atoi (q + 1);
        else
            md->g_have_mdf = atoi (p + 1);
    }
    else
    {
        // by default, write metadata file. 
        md->g_have_mdf = 1;
    }

    free (temp_string);
/*
    md->g_is_aggregator = (int *) malloc (nproc * sizeof(int));
    assert (md->g_is_aggregator);
    memset (md->g_is_aggregator, 0, nproc * sizeof(int));

    aggr_group_size = nproc / md->g_num_aggregators;
    remain = nproc - (int) aggr_group_size * md->g_num_aggregators;

    index = 0;
    for (i = 0; i < md->g_num_aggregators; i++)
    {
        md->g_is_aggregator[index] = 1;

        if (i < remain)
        {
            index += aggr_group_size + 1;
        }
        else
        {
            index += aggr_group_size;
        }
    }
*/
/*
    if (remain == 0)
    {
        md->g_color2 = rank % aggr_group_size;
    }
    else
    {
        if (rank < (aggr_group_size + 1) * remain)
        {
            md->g_color1 = rank / (aggr_group_size + 1);
            md->g_color2 = rank % (aggr_group_size + 1);
        }
        else
        {
            md->g_color1 = remain + (rank - (aggr_group_size + 1) * remain) / aggr_group_size;
            md->g_color2 = (rank - (aggr_group_size + 1) * remain)% aggr_group_size;
        }
    }
*/
}

static void adios_mpi_bgq_buffer_write (char ** buffer, uint64_t * buffer_size
                                       ,uint64_t * buffer_offset
                                       ,const void * data, uint64_t size
                                       )
{
    if (*buffer_offset + size > *buffer_size || *buffer == 0)
    {
        char * b = realloc (*buffer, *buffer_offset + size + 1000);
        if (b)
        {
            *buffer = b;
            *buffer_size = (*buffer_offset + size + 1000);
        }
        else
        {
            adios_error (err_no_memory, "Cannot allocate memory in adios_mpi_bgq_buffer_write.  "
                             "Requested: %llu\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}

static int
adios_mpi_bgq_get_striping_unit(MPI_File fh, char *filename)
{
    struct statfs fsbuf;
    int err, flag;
    uint64_t striping_unit = 1048576;
    char     value[64];
    MPI_Info info_used;

    // get striping_unit from MPI hint if it has
    MPI_File_get_info(fh, &info_used);
    MPI_Info_get(info_used, "striping_unit", 63, value, &flag);
    MPI_Info_free(&info_used);

    if (flag) return atoi(value);

    // if striping_unit is not set in MPI file info, get it from system
    err = statfs(filename, &fsbuf);
    if (err == -1) {
        log_warn ("Warning: statfs failed %s %s.\n",filename,strerror(errno));
        return striping_unit;
    }

    // set the file striping size
    return striping_unit;
}

static uint64_t
adios_mpi_bgq_striping_unit_write(MPI_File   fh
                                 ,MPI_Offset offset
                                 ,void       *buf
                                 ,uint64_t   len
                                 )
{
    uint64_t err = -1;
    MPI_Status status;
    uint64_t total_written = 0;
    uint64_t to_write = len;
    int write_len = 0;
    int count;
    char * buf_ptr = buf;

    if (len == 0)
    {
        log_debug ("shouldn't be here\n");
        return 0;
    }

    if (offset == -1) // use current position
        MPI_File_get_position(fh, &offset);
    else
        MPI_File_seek (fh, offset, MPI_SEEK_SET);

    while (total_written < len)
    {
        write_len = (to_write > MAX_MPIWRITE_SIZE) ? MAX_MPIWRITE_SIZE : to_write;
        MPI_File_write (fh, buf_ptr, write_len, MPI_BYTE, &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
        if (count != write_len)
        {
            err = count;
            break;
        }
        total_written += count;
        buf_ptr += count;
        to_write -= count;
        err = total_written;
    }

    return err;
}

struct adios_var_struct * adios_mpi_bgq_copy_var (struct adios_var_struct * v)
{
    struct adios_var_struct * v_new = (struct adios_var_struct *) 
                            malloc (sizeof (struct adios_var_struct));
    if (v_new == 0)
    {
        adios_error (err_no_memory, "MPI_BGQ method: Cannot allocate %d bytes "
            "to duplicate variable structure in adios_mpi_bgq_copy_var()\n",
            sizeof (struct adios_var_struct));
        return 0;
    }

    v_new->name = strdup (v->name);
    v_new->path = strdup (v->path);
    v_new->type = v->type;
    v_new->got_buffer = v->got_buffer;
    v_new->is_dim = v->is_dim;
    v_new->write_offset = v->write_offset;
    v_new->stats = 0;
    v_new->free_data = v->free_data;
    v_new->data = 0;
    v_new->adata = 0;
    v_new->data_size = v->data_size;
    v_new->next = 0;

    //struct adios_dimension_struct * dimensions;

    return v_new;
}

void adios_mpi_bgq_append_var (struct adios_file_struct * fd, struct adios_var_struct * v)
{
    struct adios_var_struct * root = fd->group->vars;
    if (!root)
    {
        root = v;
        return;
    }

    while (root->next)
    {
        root = root->next;
    }

    root->next = v;
}

void adios_mpi_bgq_build_global_index_v1 (char * fname
                                          ,struct adios_index_process_group_struct_v1 * pg_root
                                          ,struct adios_index_var_struct_v1 * vars_root
                                          ,struct adios_index_attribute_struct_v1 * attrs_root
                                          )
{
    int len;
    char * s;

    while (vars_root)
    {
        // Add, e.g., "/restart.bp.0/" in the beginning
        len = 1 + strlen (fname) + 1 + strlen (vars_root->var_path) + 1;
        s = (char *)malloc (len);

        sprintf (s, "%s%s%s%s", "/", fname, "/", vars_root->var_path);
        FREE (vars_root->var_path);
        vars_root->var_path = s;

        vars_root = vars_root->next;
    }

    while (attrs_root)
    {
        len = 1 + strlen (fname) + 1 + strlen (attrs_root->attr_path) + 1;
        s = (char *)malloc (len);

        sprintf (s, "%s%s%s%s", "/", fname, "/", attrs_root->attr_path);
        FREE (attrs_root->attr_path);
        attrs_root->attr_path = s;

        attrs_root = attrs_root->next;
    }

}


void * adios_mpi_bgq_do_mkdir (void * param)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) param;
    // 4 bytes for ".dir" 
    char * dir_name = malloc (strlen (fd->name) + 4 + 1);
    sprintf (dir_name, "%s%s", fd->name, ".dir");
    
    mkdir (dir_name, S_IRWXU | S_IRWXG);
  
    free (dir_name);

    return NULL;
}

void * adios_mpi_bgq_do_open_thread (void * param)
{
    struct adios_MPI_thread_data_open * td = (struct adios_MPI_thread_data_open *) param;
    int ret;

    if (td->md->file_comm_rank == 0)
    {
        ret = MPI_File_open (MPI_COMM_SELF, td->md->subfile_name
                            ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                            ,MPI_INFO_NULL
                            ,&td->md->fh
                            );

        assert (ret == MPI_SUCCESS);
        MPI_File_close (&td->md->fh);  
    }

    MPI_Barrier (td->md->file_comm);

    ret = MPI_File_open (td->md->file_comm, td->md->subfile_name
                        ,MPI_MODE_WRONLY
                        ,MPI_INFO_NULL
                        ,&td->md->fh
                        );
    return NULL;
}


void * adios_mpi_bgq_do_write_thread (void * param)
{
    struct adios_MPI_thread_data_write * td = (struct adios_MPI_thread_data_write *) param;

    uint64_t count = adios_mpi_bgq_striping_unit_write(
                               *(td->fh)
                              ,*(td->base_offset)
                              ,td->aggr_buff
                              ,*(td->total_data_size)
                              );

    if (count != *(td->total_data_size))
    {
        adios_error (err_unspecified, "Error in adios_mpi_bgq_striping_unit_write(). "
            "count = %llu != thread's total_data_size = %llu\n",
            count, *(td->total_data_size));
    }

    return NULL;
}


void adios_mpi_bgq_init (const PairStruct * parameters
                         ,struct adios_method_struct * method
                         )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    if (!adios_mpi_bgq_initialized)
    {
        adios_mpi_bgq_initialized = 1;
    }

    method->method_data = malloc (sizeof (struct adios_MPI_data_struct));
    md = (struct adios_MPI_data_struct *) method->method_data;
    md->fh = 0;
    md->mfh = 0;
    md->subfile_name = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->group_comm = MPI_COMM_NULL;
    md->rank = 0;
    md->size = 0;
    md->file_comm = MPI_COMM_NULL;
    md->file_comm_rank = 0;
    md->file_comm_size = 0;
    md->index = adios_alloc_index_v1(1); // with hashtables;

    md->g_is_aggregator = 0;
    md->g_num_aggregators = 0;
    md->g_color2 = 0;
    md->partition_id = 0;
    md->n_partitions = 0;
    md->g_offsets = 0;
    md->open_thread_data = 0;
    md->g_io_type = ADIOS_MPI_BGQ_IO_SIMPLE;
    md->g_have_mdf = 1;

    adios_buffer_struct_init (&md->b);
}

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
// Indices for the timer object
static int ADIOS_TIMER_COMM         = ADIOS_TIMING_MAX_USER_TIMERS + 0;
static int ADIOS_TIMER_IO           = ADIOS_TIMING_MAX_USER_TIMERS + 1;
static int ADIOS_TIMER_LOCALMD      = ADIOS_TIMING_MAX_USER_TIMERS + 2;
static int ADIOS_TIMER_GLOBALMD     = ADIOS_TIMING_MAX_USER_TIMERS + 3;
static int ADIOS_TIMER_AD_OPEN      = ADIOS_TIMING_MAX_USER_TIMERS + 4;
static int ADIOS_TIMER_AD_WRITE     = ADIOS_TIMING_MAX_USER_TIMERS + 5;
static int ADIOS_TIMER_OVERFLOW     = ADIOS_TIMING_MAX_USER_TIMERS + 6;
static int ADIOS_TIMER_AD_CLOSE     = ADIOS_TIMING_MAX_USER_TIMERS + 7;
#endif

int adios_mpi_bgq_open (struct adios_file_struct * fd
                       ,struct adios_method_struct * method, MPI_Comm comm
                       )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

    md->group_comm = comm;
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
    int timer_count = 8;
    char ** timer_names = (char**) malloc (timer_count * sizeof (char*) );
    timer_names [0] = "Communication   ";
    timer_names [1] = "I/O             ";
    timer_names [2] = "Local metadata  ";
    timer_names [3] = "Global metadata ";
    timer_names [4] = "adios_open()    ";
    timer_names [5] = "adios_write()   ";
    timer_names [6] = "adios_overflow()";
    timer_names [7] = "adios_close()   ";

    // Ensure both timing objects exist
    // timing_obj should get created at every open
    // prev_timing_obj should only be created at the first open
    if (fd->group)
    {
        if (!fd->group->timing_obj)
            fd->group->timing_obj = adios_timing_create (timer_count, timer_names);

        if (!fd->group->prev_timing_obj)
            fd->group->prev_timing_obj = adios_timing_create (timer_count, timer_names);
    }
#endif

    START_TIMER (ADIOS_TIMER_AD_OPEN);

    // The following code is BGQ specific
#define UNIT_A 2
#define UNIT_B 2
#define UNIT_C 4
#define UNIT_D 4
#define UNIT_E 2
    Personality_t pers;
    int Anodes, Bnodes, Cnodes, Dnodes, Enodes;
    int Acoord, Bcoord, Ccoord, Dcoord, Ecoord;
    int A_color, B_color, C_color, D_color, E_color;
    int A_blocks, B_blocks, C_blocks, D_blocks, E_blocks;
    int color;

    Kernel_GetPersonality(&pers, sizeof(pers));

    Anodes = pers.Network_Config.Anodes;
    Acoord = pers.Network_Config.Acoord;
    Bnodes = pers.Network_Config.Bnodes;
    Bcoord = pers.Network_Config.Bcoord;
    Cnodes = pers.Network_Config.Cnodes;
    Ccoord = pers.Network_Config.Ccoord;
    Dnodes = pers.Network_Config.Dnodes;
    Dcoord = pers.Network_Config.Dcoord;
    Enodes = pers.Network_Config.Enodes;
    Ecoord = pers.Network_Config.Ecoord;

    A_color  = Acoord / UNIT_A;
    B_color  = Bcoord / UNIT_B;
    C_color  = Ccoord / UNIT_C;
    D_color  = Dcoord / UNIT_D;
    E_color  = Ecoord / UNIT_E;

    // Number of blocks
    A_blocks = Anodes / UNIT_A;
    B_blocks = Bnodes / UNIT_B;
    C_blocks = Cnodes / UNIT_C;
    D_blocks = Dnodes / UNIT_D;
    E_blocks = Enodes / UNIT_E;

    color = (A_color * (B_blocks * C_blocks * D_blocks * E_blocks))
            + (B_color * (C_blocks * D_blocks * E_blocks))
            + (C_color * ( D_blocks * E_blocks))
            + (D_color * ( E_blocks))
            + E_color;

    md->partition_id = color;
    md->n_partitions = A_blocks * B_blocks * C_blocks * D_blocks * E_blocks;

    log_debug ("color = %d, n_partitions = %d, (%d,%d,%d,%d,%d)\n", 
               color, md->n_partitions, A_blocks, B_blocks, C_blocks, D_blocks, E_blocks);
    MPI_Comm_split (md->group_comm, color, md->rank, &md->file_comm);
    MPI_Comm_rank (md->file_comm, &md->file_comm_rank);
    MPI_Comm_size (md->file_comm, &md->file_comm_size);

    fd->group->process_id = md->rank;
    // we have to wait for the group_size (should_buffer)
    // to calculate stripe sizes from output sizes of the processes
    // before we can do an open for any of the modes

    int i;
    char * name, * name_no_path, * ch;
    char * d_name;
    int err;
    int sig;    // used for coordinating the MPI_File_open
    uint16_t flag;


    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);

    switch (fd->mode)
    {
        case adios_mode_read:
        case adios_mode_append:
        case adios_mode_update:
        {
            adios_error (err_invalid_file_mode, "MPI_BGQ method: specified mode is not supported.\n");
            break;
        }

        case adios_mode_write:
        {
            adios_mpi_bgq_set_aggregation_parameters (method->parameters, md);

            if (md->partition_id == 0 && md->file_comm_rank == 0)
            {
                // open metadata file
                unlink (fd->name);

                if (md->g_have_mdf)
                {
                    MPI_File_open (MPI_COMM_SELF, fd->name
                                  ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                  ,MPI_INFO_NULL
                                  ,&md->mfh
                                  );
                }

                adios_mpi_bgq_do_mkdir (fd);
            }

            MPI_Barrier (md->group_comm);

            // Check if fd->name contains path
            if (ch = strrchr (fd->name, '/'))
            {
                name_no_path = malloc (strlen (ch + 1) + 1); 
                strcpy (name_no_path, ch + 1); 
            }
            else
            {
                name_no_path = malloc (strlen (fd->name) + 1);
                strcpy (name_no_path, fd->name);
            }

            name = realloc (name, strlen (fd->name) + 5 + strlen (method->base_path) + strlen (name_no_path) + 1 + 10 + 1);
            // create the subfile name, e.g. restart.bp.1
            // 1 for '.' + 10 for subfile index + 1 for '\0'
            sprintf (name, "%s%s%s%s.%d", fd->name, ".dir/", method->base_path, name_no_path, md->partition_id);
            md->subfile_name = strdup (name);
            fd->subfile_index = (uint32_t)md->partition_id;

            free (name_no_path);

            // open subfiles. Everyone should do that, which is different
            // from AMR method.
            md->open_thread_data = (struct adios_MPI_thread_data_open *) malloc (sizeof (struct adios_MPI_thread_data_open));
            md->open_thread_data->md = md;
            md->open_thread_data->parameters = method->parameters;

            adios_mpi_bgq_do_open_thread ((void *) md->open_thread_data);
            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode, "MPI_BGQ method: Unknown file mode requested: %d\n", fd->mode);
            free (name);
            return 0;
        }
    }

    free (name);

    STOP_TIMER (ADIOS_TIMER_AD_OPEN);
    return 1;
}

static void build_file_offsets (struct adios_MPI_data_struct *md,
                                       struct adios_file_struct *fd)
{
    int i;
    if (md->group_comm != MPI_COMM_NULL)
    {
        if (md->file_comm_rank == 0)
        {
            MPI_Offset * offsets = malloc (sizeof (MPI_Offset) * md->file_comm_size);

            // round up to GPFS block size (8 MiB)
            if (fd->bytes_written % BLOCK_UNIT)
                offsets [0] =  (fd->bytes_written / BLOCK_UNIT + 1) * BLOCK_UNIT;
            else
                offsets [0] = fd->bytes_written;

            MPI_Gather (MPI_IN_PLACE, 1, MPI_LONG_LONG
                    ,offsets, 1, MPI_LONG_LONG
                    ,0, md->file_comm
                    );

            uint64_t last_pgsize = offsets [0];
            offsets [0] = md->b.end_of_pgs; // = 0 or where to append to existing data (if append was supported)
            for (i = 1; i < md->file_comm_size; i++)
            {
                uint64_t this_offset = offsets [i];
                offsets [i] = offsets [i - 1] + last_pgsize;
                last_pgsize = this_offset;
            }

            md->b.pg_index_offset =   offsets [md->file_comm_size - 1]
                                    + last_pgsize;
            MPI_Scatter (offsets, 1, MPI_LONG_LONG
                    ,MPI_IN_PLACE, 1, MPI_LONG_LONG
                    ,0, md->file_comm
                    );

            fd->current_pg->pg_start_in_file = offsets[0];
            free (offsets);
        }
        else
        {
            MPI_Offset offset[1];
            if (fd->bytes_written % BLOCK_UNIT)
            {
                offset[0] =  (fd->bytes_written / BLOCK_UNIT + 1)
                    * BLOCK_UNIT;
            }
            else
            {
                offset[0] = fd->bytes_written;
            }

            MPI_Gather (offset, 1, MPI_LONG_LONG
                    ,0, 1, MPI_LONG_LONG
                    ,0, md->file_comm
                    );

            MPI_Scatter (0, 1, MPI_LONG_LONG
                    ,offset, 1, MPI_LONG_LONG
                    ,0, md->file_comm
                    );
            fd->current_pg->pg_start_in_file = offset[0];
        }
    }
    else
    {
        md->b.pg_index_offset = fd->bytes_written;     
        fd->current_pg->pg_start_in_file = md->b.end_of_pgs; // 0 or where to append to existing data
    }
}


enum BUFFERING_STRATEGY adios_mpi_bgq_should_buffer (struct adios_file_struct * fd
                                                    ,struct adios_method_struct * method
                                                    )
{
    return stop_on_overflow;
}

void adios_mpi_bgq_write (struct adios_file_struct * fd
                         ,struct adios_var_struct * v
                         ,const void * data
                         ,struct adios_method_struct * method
                         )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;
    START_TIMER (ADIOS_TIMER_AD_WRITE);
    if (v->got_buffer == adios_flag_yes)
    {
        if (data != v->data)  // if the user didn't give back the same thing
        {
            if (v->free_data == adios_flag_yes)
            {
                free (v->adata);
                adios_method_buffer_free (v->data_size);
            }
        }
        else
        {
            // we already saved all of the info, so we're ok.
            return;
        }
    }

    STOP_TIMER (ADIOS_TIMER_AD_WRITE);
}

void adios_mpi_bgq_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
    uint64_t mem_allowed;

    if (*size == 0)
    {
        *buffer = 0;

        return;
    }

    if (v->adata && v->free_data)
    {
        adios_method_buffer_free (v->data_size);
        free (v->adata);
    }

    mem_allowed = adios_method_buffer_alloc (*size);
    if (mem_allowed == *size)
    {
        *buffer = malloc (*size);
        if (!*buffer)
        {
            adios_method_buffer_free (mem_allowed);
            adios_error (err_no_memory, 
                    "MPI_BGQ method: Out of memory allocating %llu bytes for variable %s\n",
                    *size ,v->name);
            v->got_buffer = adios_flag_no;
            v->free_data = adios_flag_no;
            v->data_size = 0;
            v->data = 0;
            *size = 0;
            *buffer = 0;
        }
        else
        {
            v->got_buffer = adios_flag_yes;
            v->free_data = adios_flag_yes;
            v->data_size = mem_allowed;
            v->data = *buffer;
        }
    }
    else
    {
        adios_method_buffer_free (mem_allowed);
        adios_error (err_buffer_overflow, 
                "MPI_BGQ method: OVERFLOW: Cannot allocate requested buffer of %llu "
                "bytes for %s. Allowed max size is %llu\n",
                *size, v->name, mem_allowed);
        *size = 0;
        *buffer = 0;
    }
}

void adios_mpi_bgq_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )
{
    v->data = v->adata = buffer;
    v->data_size = buffer_size;
}

static void adios_mpi_bgq_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;
    struct adios_var_struct * v = fd->group->vars;

    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer = 0;
    data.buffer_len = 0;

    uint32_t version = md->b.version & ADIOS_VERSION_NUM_MASK;
    switch (version)
    {
        case 1:
        case 2:
        case 3:
        {
            // the three section headers
            struct adios_process_group_header_struct_v1 pg_header;
            struct adios_vars_header_struct_v1 vars_header;
            struct adios_attributes_header_struct_v1 attrs_header;

            struct adios_var_header_struct_v1 var_header;
            struct adios_var_payload_struct_v1 var_payload;
            struct adios_attribute_struct_v1 attribute;

            uint64_t i;

            adios_init_buffer_read_process_group (&md->b);
            MPI_File_seek (md->fh, md->b.read_pg_offset
                          ,MPI_SEEK_SET
                          );
            MPI_File_read (md->fh, md->b.buff, md->b.read_pg_size, MPI_BYTE
                          ,&md->status
                          );
            adios_parse_process_group_header_v1 (&md->b, &pg_header);

            adios_parse_vars_header_v1 (&md->b, &vars_header);

            for (i = 0; i < vars_header.count; i++)
            {
                memset (&var_payload, 0
                       ,sizeof (struct adios_var_payload_struct_v1)
                       );
                adios_parse_var_data_header_v1 (&md->b, &var_header);

                struct adios_var_struct * v1 = v;
                while (v1)
                {
                    if (   strcasecmp (var_header.name, v1->name)
                        || strcasecmp (var_header.path, v1->path)
                       )
                    {
                        v1 = v1->next;
                    }
                    else
                    {
                        break;
                    }
                }

                if (v1)
                {
                    var_payload.payload = v1->adata;
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,&var_payload
                                                    ,v1->data_size
                                                    );
                }
                else
                {
                    log_warn ("MPI AMR method read: skipping name: %s path: %s\n",
                           var_header.name, var_header.path);
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,NULL, 0
                                                    );
                }

                adios_clear_var_header_v1 (&var_header);
            }

#if 1
            adios_parse_attributes_header_v1 (&md->b, &attrs_header);

            for (i = 0; i < attrs_header.count; i++)
            {
                adios_parse_attribute_v1 (&md->b, &attribute);
                adios_clear_attribute_v1 (&attribute);
            }
#endif
            adios_clear_process_group_header_v1 (&pg_header);

            break;
        }

        default:
            adios_error (err_invalid_file_version, 
                    "MPI_BGQ method read: file version unknown: %u\n", version);
            return;
    }

    adios_buffer_struct_clear (&md->b);
}

static
uint32_t adios_mpi_bgq_calculate_attributes_size (struct adios_file_struct * fd)
{
    uint32_t overhead = 0;
    struct adios_attribute_struct * a = fd->group->attributes;

    overhead += 2; // attributes count
    overhead += 8; // attributes length

    while (a)
    {
        overhead += adios_calc_attribute_overhead_v1 (a);

        a = a->next;
    }

    return overhead;
}


void adios_mpi_bgq_simple_close (struct adios_file_struct * fd
                                 ,struct adios_method_struct * method
                                 )
{

    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    START_TIMER (ADIOS_TIMER_AD_CLOSE);
    switch (fd->mode)
    {
        case adios_mode_read:
        case adios_mode_append:
        case adios_mode_update:
        {
            adios_error (err_invalid_file_mode, 
                    "Only \"w\" mode is supported by MPI_BGQ\n");
            break;
        }
        case adios_mode_write:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            int * disp = 0, * sendbuf = 0, * recvbuf = 0, * attr_sizes = 0;
            void * * recv_buff = 0;
            struct adios_MPI_thread_data_write write_thread_data;
            int i, max_data_size = 0, total_data_size = 0, total_data_size1 = 0;
            MPI_Comm new_comm2;

            START_TIMER (ADIOS_TIMER_COMM);
            MPI_Comm_split (md->group_comm, md->file_comm_rank, md->rank, &new_comm2);
            STOP_TIMER (ADIOS_TIMER_COMM);

            /* Write data block (PG) */
            // if we need to write > 2 GB, need to do it in parts
            // since count is limited to MAX_MPIWRITE_SIZE (signed 32-bit max).
            uint64_t bytes_written = 0; 
            int32_t to_write = 0; 
            if (fd->bytes_written > MAX_MPIWRITE_SIZE)
            {    
                to_write = MAX_MPIWRITE_SIZE;
            }    
            else 
            {    
                to_write = (int32_t) fd->bytes_written;
            }    

            // figure out the offsets
            // before writing out the buffer and build the index based on target offsets
            build_file_offsets (md, fd);

            START_TIMER (ADIOS_TIMER_IO);
            while (bytes_written < fd->bytes_written)
            {
                // everyone writes their data
                MPI_File_seek (md->fh, fd->current_pg->pg_start_in_file + bytes_written
                        ,MPI_SEEK_SET
                        );   
                int err = MPI_File_write (md->fh, fd->buffer + bytes_written
                        ,to_write, MPI_BYTE, &md->status
                        );

                if (err != MPI_SUCCESS)
                {
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    adios_error (err_write_error,
                            "MPI_BGQ method, rank %d: adios_close(): writing of buffered data "
                            "[%llu..%llu] to file %s failed: '%s'\n",
                            md->rank, bytes_written, bytes_written+to_write-1,
                            fd->name, e);
                }
                bytes_written += to_write;
                if (fd->bytes_written > bytes_written)
                {
                    if (fd->bytes_written - bytes_written > MAX_MPIWRITE_SIZE)
                    {
                        to_write = MAX_MPIWRITE_SIZE;
                    }
                    else
                    {
                        to_write = fd->bytes_written - bytes_written;
                    }
                }
            }
            STOP_TIMER (ADIOS_TIMER_IO);

            /* Build the local index */
            START_TIMER (ADIOS_TIMER_LOCALMD);
            // build index appending to any existing index
            adios_build_index_v1 (fd, md->index);
            STOP_TIMER (ADIOS_TIMER_LOCALMD);

            // if collective, gather the indexes from the rest and call
            if (md->group_comm != MPI_COMM_NULL)
            {
                // Collect index from the processors
                // that belong to the same file comm.
                if (md->file_comm_rank == 0)
                {
                    int * index_sizes = malloc (4 * md->file_comm_size);
                    int * index_offsets = malloc (4 * md->file_comm_size);
                    char * recv_buffer = 0;
                    uint32_t size = 0, total_size = 0;
                    int i;

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, md->file_comm
                               );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    for (i = 0; i < md->file_comm_size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);
                    assert (recv_buffer);

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, md->file_comm
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    START_TIMER (ADIOS_TIMER_LOCALMD);
                    for (i = 1; i < md->file_comm_size; i++)
                    {
                        md->b.buff = recv_buffer + index_offsets [i];
                        md->b.length = index_sizes [i];
                        md->b.offset = 0;

                        adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (md->index, new_pg_root,
                                              new_vars_root ,new_attrs_root,0);
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }

                    md->b.buff = buffer_save;
                    md->b.length = buffer_size_save;
                    md->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);
                    STOP_TIMER (ADIOS_TIMER_LOCALMD);
                }
                else
                {
                    START_TIMER (ADIOS_TIMER_LOCALMD);
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,0, md->index);
                    STOP_TIMER (ADIOS_TIMER_LOCALMD);

                    uint32_t temp_buffer_size = buffer_offset;

                    START_TIMER (ADIOS_TIMER_COMM);
/*
                    MPI_Gather ((uint32_t *)&buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->file_comm
                               );
*/
                    MPI_Gather (&temp_buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->file_comm
                               );
/*
                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->file_comm
                                );
*/
                    MPI_Gatherv (buffer, temp_buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->file_comm
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);
                }
            }

            // Rank 0 within each file comm writes
            // out indexes in each subfile
            if (md->file_comm_rank == 0)
            {
                int err;

                START_TIMER (ADIOS_TIMER_LOCALMD);
                uint64_t index_start = md->b.pg_index_offset;

                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,index_start, md->index);
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                MPI_File_seek (md->fh, md->b.pg_index_offset, MPI_SEEK_SET);
                {
                    uint64_t total_written = 0;
                    uint64_t to_write = buffer_offset;
                    int write_len = 0;
                    int count;
                    char * buf_ptr = buffer;
      
                    while (total_written < buffer_offset)
                    {
                        write_len = (to_write > MAX_MPIWRITE_SIZE) ? MAX_MPIWRITE_SIZE : to_write;
                        err = MPI_File_write (md->fh, buf_ptr, write_len, MPI_BYTE, &md->status);

                        MPI_Get_count(&md->status, MPI_BYTE, &count);
                        if (count != write_len)
                        {
                            log_error ("MPI_BGQ method, rank %d: Need to do multi-write 1 (tried: "
                                    "%d wrote: %d) errno %d\n",
                                    md->rank, write_len, count, errno);
                            err = count;
                            break;
                        }
                        total_written += count;
                        buf_ptr += count;
                        to_write -= count;
                        //err = total_written;
                    }
                }
                if (err != MPI_SUCCESS)
                {
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    adios_error (err_write_error,
                            "MPI_BGQ method, rank %d: adios_close(): writing of index data "
                            "of %llu bytes to file %s failed: '%s'\n",
                            md->rank, buffer_offset, fd->name, e);
                }
                STOP_TIMER (ADIOS_TIMER_LOCALMD);
            }

            if (md->g_have_mdf)
            {
                // collect index among aggregators
                if (md->file_comm_rank == 0)
                {
                    if (md->partition_id == 0)
                    {
                        int * index_sizes = malloc (4 * md->n_partitions);
                        int * index_offsets = malloc (4 * md->n_partitions);
                        char * recv_buffer = 0;
                        uint32_t size = 0, total_size = 0;

                        START_TIMER (ADIOS_TIMER_COMM);
                        MPI_Gather (&size, 1, MPI_INT
                                   ,index_sizes, 1, MPI_INT
                                   ,0, new_comm2
                                   );
                        STOP_TIMER (ADIOS_TIMER_COMM);

                        for (i = 0; i < md->n_partitions; i++)
                        {
                            index_offsets [i] = total_size;
                            total_size += index_sizes [i];
                        }

                        recv_buffer = malloc (total_size);
                        assert (recv_buffer);

                        START_TIMER (ADIOS_TIMER_COMM);
                        MPI_Gatherv (&size, 0, MPI_BYTE
                                    ,recv_buffer, index_sizes, index_offsets
                                    ,MPI_BYTE, 0, new_comm2
                                    );
                        STOP_TIMER (ADIOS_TIMER_COMM);

                        char * buffer_save = md->b.buff;
                        uint64_t buffer_size_save = md->b.length;
                        uint64_t offset_save = md->b.offset;

                        START_TIMER (ADIOS_TIMER_GLOBALMD);
                        for (i = 1; i < md->n_partitions; i++)
                        {
                            md->b.buff = recv_buffer + index_offsets [i];
                            md->b.length = index_sizes [i];
                            md->b.offset = 0;

                            adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                            adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                            adios_parse_attributes_index_v1 (&md->b
                                                            ,&new_attrs_root
                                                            );

                            adios_merge_index_v1 (md->index
                                                 ,new_pg_root, new_vars_root
                                                 ,new_attrs_root, 0
                                                 );
                            new_pg_root = 0;
                            new_vars_root = 0;
                            new_attrs_root = 0;
                        }

                        md->b.buff = buffer_save;
                        md->b.length = buffer_size_save;
                        md->b.offset = offset_save;

                        free (recv_buffer);
                        free (index_sizes);
                        free (index_offsets);
                        STOP_TIMER (ADIOS_TIMER_GLOBALMD);
                    }
                    else
                    {
                        char * buffer2 = 0;
                        uint64_t buffer_size2 = 0;
                        uint64_t buffer_offset2 = 0;

                        START_TIMER (ADIOS_TIMER_GLOBALMD);
                        adios_write_index_v1 (&buffer2, &buffer_size2, &buffer_offset2
                                             ,0, md->index
                                             );
                        STOP_TIMER (ADIOS_TIMER_GLOBALMD);
                        uint32_t temp_buffer_size2 = buffer_offset2;

                        START_TIMER (ADIOS_TIMER_COMM);
/*
                        MPI_Gather (&buffer_size2, 1, MPI_INT
                                   ,0, 0, MPI_INT
                                   ,0, new_comm2
                                   );
                        MPI_Gatherv (buffer2, buffer_size2, MPI_BYTE
                                    ,0, 0, 0, MPI_BYTE
                                    ,0, new_comm2
                                    );
*/
                        MPI_Gather (&temp_buffer_size2, 1, MPI_INT
                                   ,0, 0, MPI_INT
                                   ,0, new_comm2
                                   );
                        MPI_Gatherv (buffer2, temp_buffer_size2, MPI_BYTE
                                    ,0, 0, 0, MPI_BYTE
                                    ,0, new_comm2
                                    );
                        STOP_TIMER (ADIOS_TIMER_COMM);

                        if (buffer2)
                        {
                            free (buffer2);
                            buffer2 = 0;
                            buffer_size2 = 0;
                            buffer_offset2 = 0;
                        }
                    }
                }

                // write out the metadata file from rank 0
                if (md->partition_id == 0 && md->file_comm_rank == 0)
                {
                    MPI_File m_file;
                    char * global_index_buffer = 0;
                    uint64_t global_index_buffer_size = 0;
                    uint64_t global_index_buffer_offset = 0;
                    uint64_t global_index_start = 0;
                    uint16_t flag = 0;

                    START_TIMER (ADIOS_TIMER_GLOBALMD);
                    adios_write_index_v1 (&global_index_buffer, &global_index_buffer_size
                                         ,&global_index_buffer_offset, global_index_start
                                         ,md->index
                                         );

                    flag |= ADIOS_VERSION_HAVE_SUBFILE;

                    adios_write_version_flag_v1 (&global_index_buffer
                                                ,&global_index_buffer_size
                                                ,&global_index_buffer_offset
                                                ,flag
                                                );

                    adios_mpi_bgq_striping_unit_write(
                                      md->mfh,
                                      -1,
                                      global_index_buffer,
                                      global_index_buffer_offset
                                      );

                    if (global_index_buffer)
                    {
                        free (global_index_buffer);
                        global_index_buffer = 0;
                        global_index_buffer_size = 0;
                        global_index_buffer_offset = 0;
                    }
                    STOP_TIMER (ADIOS_TIMER_GLOBALMD);
                }
            }

            FREE (buffer);
            buffer_size = 0;
            buffer_offset = 0;

            adios_clear_index_v1 (md->index);
            //md->index = 0;
            md->g_num_aggregators = 0;
            md->g_color2 = 0;

            FREE (md->subfile_name);
            FREE (md->g_is_aggregator);
            FREE (md->g_offsets);
            FREE (md->open_thread_data);
        }

        break;
    }

    if (md && md->fh)
        MPI_File_close (&md->fh);

    if (md && md->g_have_mdf && md->mfh)
        MPI_File_close (&md->mfh);

    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
    {   
        md->group_comm = MPI_COMM_NULL;
    }

    md->fh = 0;
    md->mfh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));

    STOP_TIMER (ADIOS_TIMER_AD_CLOSE);

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS

    //Finished timing this cycle, swap the timing buffers
    adios_timing_destroy(fd->group->prev_timing_obj);
    fd->group->prev_timing_obj = fd->group->timing_obj;
    fd->group->timing_obj = 0;

    // prev_timing_obj points to unwritten timing info, timing_obj is
    // ready to allocate at the next open

#endif

    return;
}

void adios_mpi_bgq_bg_close (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            )
{
#if 0
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        case adios_mode_append:
        case adios_mode_update:
        {
            adios_error (err_invalid_file_mode, 
                    "Only \"w\" mode is supported by MPI_BGQ\n");
            break;
        }
        case adios_mode_write:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = md->b.pg_index_offset, index_start1;
            int * pg_sizes = 0, * disp = 0, * sendbuf = 0, * recvbuf = 0, * attr_sizes = 0;
            void * aggr_buff = 0, * recv_buff = 0;
            struct adios_MPI_thread_data_write write_thread_data;
            int i, new_rank, new_group_size, new_rank2, new_group_size2, max_data_size = 0, total_data_size = 0, total_data_size1 = 0;
            MPI_Comm new_comm2;

            MPI_Comm_split (md->group_comm, md->file_comm_rank, md->rank, &new_comm2);
            MPI_Comm_rank (new_comm2, &new_rank2);
            MPI_Comm_size (new_comm2, &new_group_size2);

            // if not merge PG's on the aggregator side
            if (fd->shared_buffer == adios_flag_no)
            {
                log_warn ("The ADIOS buffer in the XML is not large enough for buffering.\n");
            }
            else
            {
                struct adios_bp_buffer_struct_v1 b;
                struct adios_process_group_header_struct_v1 pg_header;
                struct adios_vars_header_struct_v1 vars_header;
                int pg_size, header_size;
                uint64_t vars_count_offset;
                MPI_Request request;
                MPI_Status status;

                pg_size = fd->bytes_written;

                pg_sizes = (int *) malloc (new_group_size * 4);
                disp = (int *) malloc (new_group_size * 4);
                if (pg_sizes == 0 || disp == 0)
                {
                    adios_error (err_no_memory, "MPI_BGQ method: Cannot allocate memory "
                                "for merging process blocks (mpi_amr_bg_close)\n");
                    return;
                }

                MPI_Allgather (&pg_size, 1, MPI_INT
                              ,pg_sizes, 1, MPI_INT
                              ,md->file_comm);

                disp[0] = 0;
                max_data_size = pg_size;

                for (i = 1; i < new_group_size; i++)
                {
                    disp[i] = disp[i - 1] + pg_sizes[i - 1];
                    max_data_size = (pg_sizes[i] > max_data_size) ? pg_sizes[i] : max_data_size;
                }

                if (is_aggregator (md->rank))
                {
                    if (2 * max_data_size > MAX_AGG_BUF)
                    {
                        log_warn ("MPI_BGQ method (BG): The max allowed aggregation "
                                "buffer is %llu bytes.\n"
                                "But this ADIOS method needs %llu bytes for aggregation\n",
                                MAX_AGG_BUF, 2 * max_data_size);
                    }

                    aggr_buff = malloc (max_data_size);
                    recv_buff = malloc (max_data_size);
                    if (aggr_buff == 0 || recv_buff == 0)
                    {
                        adios_error (err_no_memory, "MPI_BGQ method (BG): Cannot allocate "
                                    "2 x %llu bytes for aggregation buffers.\n", 
                                    max_data_size);
                        return;
                    }
                }
                else
                {
                    if (max_data_size > MAX_AGG_BUF)
                    {
                        log_warn ("MPI_BGQ method (BG): The max allowed aggregation "
                                  "buffer is %llu bytes.\n",
                                  MAX_AGG_BUF);
                    }

                    recv_buff = malloc (max_data_size);
                    if (recv_buff == 0)
                    {
                        adios_error (err_no_memory, "MPI_BGQ method (BG): Cannot allocate "
                                    "%llu bytes for receive buffer.\n", 
                                    max_data_size);
                        return;
                    }
                }

                total_data_size = disp[new_group_size - 1]
                                + pg_sizes[new_group_size - 1];

                if (is_aggregator (md->rank))
                {
                    uint64_t index_start1 = 0;
                    for (i = 0; i < new_group_size; i++)
                    {
                        if (i + 1 < new_group_size)
                        {
                            MPI_Irecv (recv_buff, pg_sizes[i + 1], MPI_BYTE, new_rank + 1
                                      ,0, md->file_comm, &request);
                        }

                        write_thread_data.fh = &md->fh;
                        write_thread_data.base_offset = &index_start1;
                        write_thread_data.aggr_buff = (i == 0) ? fd->buffer : aggr_buff;
                        write_thread_data.total_data_size = &pg_sizes[i];

                        // This write call is not threaded
                        adios_mpi_bgq_do_write_thread ((void *) &write_thread_data);

                        index_start1 += pg_sizes[i];

                        if (i + 1 < new_group_size)
                        {
                            MPI_Wait (&request, &status);

                            memcpy (aggr_buff, recv_buff, pg_sizes[i + 1]);
                        }
                    }
                }
                else
                {
                    if (new_rank == new_group_size - 1)
                    {
                        MPI_Send (fd->buffer, pg_size, MPI_BYTE, new_rank - 1
                                 ,0, md->file_comm);
                    }
                    else
                    {
                        for (i = new_rank + 1; i < new_group_size; i++)
                        {
                            // Recv data from upstream rank
                            MPI_Irecv (recv_buff, pg_sizes[i], MPI_BYTE, new_rank + 1
                                      ,0, md->file_comm, &request);

                            if (i == new_rank + 1)
                                // Send my data to downstream rank
                                MPI_Send (fd->buffer, pg_size, MPI_BYTE, new_rank - 1
                                         ,0, md->file_comm);

                            MPI_Wait (&request, &status);
                            // Send it to downstream rank
                            MPI_Send (recv_buff, pg_sizes[i], MPI_BYTE, new_rank - 1
                                     ,0, md->file_comm);
                        }
                    }
                }

                FREE (aggr_buff);
                FREE (recv_buff);
            }

            // build index appending to any existing index
            adios_build_index_v1 (fd, md->index);

            if (fd->shared_buffer == adios_flag_yes)
            {
                if (!is_aggregator(md->rank))
                {
                    uint64_t var_offset_to_add = 0, attr_offset_to_add = 0;
                    uint64_t var_base_offset = 0, attr_base_offset = 0;

                    // change to relative offset
                    if (md->old_vars_root)
                    {
                        var_base_offset = md->old_vars_root->characteristics [0].offset;
                    }

                    if (md->old_attrs_root)
                    {
                        attr_base_offset = md->old_attrs_root->characteristics [0].offset;
                    }

                    for (i = 0; i < new_rank; i++)
                    {
                        attr_offset_to_add += pg_sizes[i];
                        var_offset_to_add += pg_sizes[i];
                    }

                    adios_mpi_bgq_add_offset (var_offset_to_add
                                             ,attr_offset_to_add
                                             ,md->old_vars_root
                                             ,md->old_attrs_root
                                             );
                }

                // pg_sizes, disp are no longer needed from this point on.
                FREE (pg_sizes);
                FREE (disp);
            }

            // if collective, gather the indexes from the rest and call
            if (md->group_comm != MPI_COMM_NULL)
            {
                // Collect index from all MPI processors
                if (is_aggregator (md->rank))
                {
                    int * index_sizes = malloc (4 * new_group_size);
                    int * index_offsets = malloc (4 * new_group_size);
                    char * recv_buffer = 0;
                    uint32_t size = 0, total_size = 0;
                    int i;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, md->file_comm
                               );

                    for (i = 0; i < new_group_size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);

                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, md->file_comm
                                );

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    for (i = 1; i < new_group_size; i++)
                    {
                        md->b.buff = recv_buffer + index_offsets [i];
                        md->b.length = index_sizes [i];
                        md->b.offset = 0;

                        adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (md->index
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root, 0
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }

                    md->b.buff = buffer_save;
                    md->b.length = buffer_size_save;
                    md->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);
                }
                else
                {
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,0, md->index
                                         );

                    uint32_t temp_buffer_size = buffer_size;

/*
                    MPI_Gather (&buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, new_comm
                               );
*/
                    MPI_Gather (&temp_buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->file_comm
                               );
/*
                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, new_comm
                                );
*/
                    MPI_Gatherv (buffer, temp_buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->file_comm
                                );
                }
            }

            // write out indexes in each subfile
            if (is_aggregator (md->rank))
            {
                uint32_t flag = 0;
                index_start = total_data_size;

                adios_write_index_v1 (&buffer, &buffer_size
                                     ,&buffer_offset, index_start
                                     ,md->index
                                     );

                adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);

                if (fd->shared_buffer == adios_flag_yes)
                {
                    index_start = -1;
                    total_data_size1 = buffer_offset;

                    write_thread_data.fh = &md->fh;
                    write_thread_data.base_offset = &index_start;
                    write_thread_data.aggr_buff = buffer;
                    write_thread_data.total_data_size = &total_data_size1;
                    adios_mpi_bgq_do_write_thread ((void *) &write_thread_data); 
                }
            }

            // collect index among aggregators
            if (is_aggregator (md->rank))
            {
                if (md->rank == 0)
                {
                    int * index_sizes = malloc (4 * new_group_size2);
                    int * index_offsets = malloc (4 * new_group_size2);
                    char * recv_buffer = 0;
                    uint32_t size = 0, total_size = 0;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, new_comm2
                               );

                    for (i = 0; i < new_group_size2; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);

                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, new_comm2
                                );

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    for (i = 1; i < new_group_size2; i++)
                    {
                        md->b.buff = recv_buffer + index_offsets [i];
                        md->b.length = index_sizes [i];
                        md->b.offset = 0;

                        adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (md->index
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root, 0
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }

                    md->b.buff = buffer_save;
                    md->b.length = buffer_size_save;
                    md->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);
                }
                else
                {
                    char * buffer2 = 0;
                    uint64_t buffer_size2 = 0;
                    uint64_t buffer_offset2 = 0;

                    adios_write_index_v1 (&buffer2, &buffer_size2, &buffer_offset2
                                         ,0, md->index
                                         );
                    uint32_t temp_buffer_size2 = buffer_size2;

/*
                    MPI_Gather (&buffer_size2, 1, MPI_INT
                               ,0, 0, MPI_INT
                               ,0, new_comm2
                               );
                    MPI_Gatherv (buffer2, buffer_size2, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, new_comm2
                                );
*/
                    MPI_Gather (&temp_buffer_size2, 1, MPI_INT
                               ,0, 0, MPI_INT
                               ,0, new_comm2
                               );
                    MPI_Gatherv (buffer2, temp_buffer_size2, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, new_comm2
                                );


                    if (buffer2)
                    {
                        free (buffer2);
                        buffer2 = 0;
                        buffer_size2 = 0;
                        buffer_offset2 = 0;
                    }
                }
            }

            // write out the metadata file from rank 0
            if (md->rank == 0)
            {
                MPI_File m_file;
                char * global_index_buffer = 0;
                uint64_t global_index_buffer_size = 0;
                uint64_t global_index_buffer_offset = 0;
                uint64_t global_index_start = 0;
                uint16_t flag = 0;

                adios_write_index_v1 (&global_index_buffer, &global_index_buffer_size
                                     ,&global_index_buffer_offset, global_index_start
                                     ,md->index
                                     );

                flag |= ADIOS_VERSION_HAVE_SUBFILE;

                adios_write_version_flag_v1 (&global_index_buffer
                                            ,&global_index_buffer_size
                                            ,&global_index_buffer_offset
                                            ,flag
                                            );

                adios_mpi_bgq_striping_unit_write(
                                  md->mfh,
                                  -1,
                                  global_index_buffer,
                                  global_index_buffer_offset
                                  );

                if (global_index_buffer)
                {
                    free (global_index_buffer);
                    global_index_buffer = 0;
                    global_index_buffer_size = 0;
                    global_index_buffer_offset = 0;
                }
            }

            if (is_aggregator (md->rank))
            {
                FREE (aggr_buff);
            }

            FREE (buffer);
            buffer_size = 0;
            buffer_offset = 0;

            adios_clear_index_v1 (md->index);
            md->index = 0;

            md->g_num_aggregators = 0;
            md->g_color2 = 0;

            FREE (md->subfile_name);
            FREE (md->g_is_aggregator);
            FREE (md->g_offsets);
            FREE (md->open_thread_data);
        }

        break;
    }

    if (md && md->fh)
        MPI_File_close (&md->fh);

    if (md && md->mfh)
        MPI_File_close (&md->mfh);

    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
    {   
        md->group_comm = MPI_COMM_NULL;
    }

    md->fh = 0;
    md->mfh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
#endif
    return;
}

void adios_mpi_bgq_ag_close (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            )
{
#if 0
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        case adios_mode_append:
        {
            adios_error (err_invalid_file_mode, 
                        "Only \"w\" mode is supported by MPI_BGQ Aggregation IO\n");
            break;
        }
        case adios_mode_write:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = md->b.pg_index_offset, index_start1;
            int * pg_sizes = 0, * disp = 0, * sendbuf = 0, * recvbuf = 0, * attr_sizes = 0;
            void * aggr_buff = 0;
            struct adios_MPI_thread_data_write write_thread_data;
            int i, new_rank, new_group_size, new_rank2, new_group_size2, total_data_size = 0, total_data_size1 = 0;;
            MPI_Comm new_comm, new_comm2;

/*
            MPI_Comm_split (md->group_comm, md->g_color1, md->rank, &new_comm);
            MPI_Comm_rank (new_comm, &new_rank);
            MPI_Comm_size (new_comm, &new_group_size);
*/
            MPI_Comm_split (md->group_comm, md->g_color2, md->rank, &new_comm2);
            MPI_Comm_rank (new_comm2, &new_rank2);
            MPI_Comm_size (new_comm2, &new_group_size2);

            if (fd->shared_buffer == adios_flag_no)
            {
                MPI_Offset new_off;
                // set it up so that it will start at 0, but have correct sizes
                MPI_File_get_position (md->fh, &new_off);
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_vars_v1 (fd);
                // fd->vars_start gets updated with the size written

                uint64_t count;
                if (is_aggregator(md->rank))
                {
                    count = adios_mpi_bgq_striping_unit_write(
                                   md->fh
                                  ,md->vars_start
                                  ,fd->buffer
                                  ,md->vars_header_size
                                  );

                    if (count != md->vars_header_size)
                    {
                        log_warn ("d:MPI_BGQ method tried to write %llu, only wrote %d\n",
                                md->vars_header_size, count);
                    }
                }
                fd->offset = 0;
                fd->bytes_written = 0;
                adios_shared_buffer_free (&md->b);

                adios_write_open_attributes_v1 (fd);
                md->vars_start = new_off;
                md->vars_header_size = fd->offset;

                MPI_File_seek (md->fh, new_off + md->vars_header_size
                              ,MPI_SEEK_SET
                              ); // go back to end, but after attr header

                fd->base_offset += fd->offset;  // add size of header
                fd->offset = 0;
                fd->bytes_written = 0;

                if (!fd->group->process_id) { // from ADIOS 1.4, only rank 0 writes attributes
                    while (a)
                    {
                        adios_write_attribute_v1 (fd, a);

                        int bytes_written[new_group_size];
                        int disp[new_group_size];
                        int total_size = 0;
                        void * aggr_buff;

                        MPI_Gather (&fd->bytes_written, 1, MPI_INT
                                ,bytes_written, 1, MPI_INT
                                ,0, new_comm
                                );

                        disp[0] = 0;
                        for (i = 1; i < new_group_size; i++)
                        {
                            disp[i] = disp[i - 1] + bytes_written[i - 1];
                        }
                        total_size += disp[new_group_size - 1]
                            + bytes_written[new_group_size - 1];

                        if (is_aggregator(md->rank))
                        {
                            aggr_buff = malloc (total_size);
                            if (aggr_buff == 0)
                            {
                                adios_error (err_no_memory, 
                                        "MPI_BGQ method (AG): Cannot allocate aggregation buffer.\n"
                                        "Need to increase the number of aggregators.\n"
                                        );
                                return;
                            }
                        }

                        MPI_Gatherv (fd->buffer, fd->bytes_written, MPI_BYTE
                                ,aggr_buff, bytes_written, disp, MPI_BYTE
                                ,0, new_comm);

                        if (is_aggregator (md->rank))
                        {
                            count = adios_mpi_bgq_striping_unit_write(
                                    md->fh,
                                    -1,
                                    aggr_buff, //fd->buffer,
                                    total_size //fd->bytes_written,
                                    );

                            if (count != total_size)
                            {
                                log_warn ("e:MPI_BGQ method tried to write %llu, only wrote %llu\n",
                                        fd->bytes_written, count);
                            }
                        }

                        // Broadcast new offsets to all processors in the communicator.
                        uint64_t new_offsets[new_group_size];

                        if (is_aggregator (md->rank))
                        {
                            new_offsets[0] = a->write_offset;
                            for (i = 1; i < new_group_size; i++)
                            {
                                new_offsets[i] = new_offsets[i - 1] + bytes_written[i - 1];
                            }
                        }

                        MPI_Bcast (new_offsets, new_group_size, MPI_LONG_LONG, 0, new_comm);

                        a->write_offset = new_offsets[new_rank];

                        fd->base_offset += count;
                        fd->offset = 0;
                        fd->bytes_written = 0;
                        adios_shared_buffer_free (&md->b);

                        a = a->next;
                    }
                }

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);

                // fd->vars_start gets updated with the size written
                if (is_aggregator(md->rank))
                {
                    *(uint16_t *)fd->buffer = *(uint16_t *)fd->buffer * new_group_size;

                    count = adios_mpi_bgq_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size
                                  );

                    if (count != md->vars_header_size)
                    {
                        log_warn ("f:MPI_BGQ method tried to write %llu, only wrote %llu\n",
                                  md->vars_header_size, count);
                    }
                }

                fd->offset = 0;
                fd->bytes_written = 0;

                MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
            }

            // if not merge PG's on the aggregator side
            if (fd->shared_buffer == adios_flag_yes)
            {
                struct adios_bp_buffer_struct_v1 b;
                struct adios_process_group_header_struct_v1 pg_header;
                struct adios_vars_header_struct_v1 vars_header;
                int pg_size, header_size;
                uint64_t vars_count_offset;

                pg_size = fd->bytes_written;
                pg_sizes = (int *) malloc (new_group_size * 4);
                disp = (int *) malloc (new_group_size * 4);
                if (pg_sizes == 0 || disp == 0)
                {
                    adios_error (err_no_memory, 
                            "MPI_BGQ method (AG): Cannot allocate buffers (%d bytes) "
                            "for merging process blocks.\n",
                            2*4*new_group_size
                            );
                    return;
                }

                MPI_Allgather (&pg_size, 1, MPI_INT
                              ,pg_sizes, 1, MPI_INT
                              ,new_comm);

                disp[0] = 0;
                for (i = 1; i < new_group_size; i++)
                {
                    disp[i] = disp[i - 1] + pg_sizes[i - 1];
                }
                total_data_size = disp[new_group_size - 1]
                                + pg_sizes[new_group_size - 1];

                if (is_aggregator (md->rank))
                {
                    if (total_data_size > MAX_AGG_BUF)
                    {
                        log_warn ("The max allowed aggregation buffer is %llu. Requested %llu.\n"
                                "Need to increase the number of aggregators.\n",
                                MAX_AGG_BUF, total_data_size);
                    }
                    aggr_buff = malloc (total_data_size);
                    if (aggr_buff == 0)
                    {
                        adios_error (err_no_memory, 
                                "MPI_BGQ method (AG): Cannot allocate %llu bytes "
                                "for aggregation buffer.\n"
                                "Need to increase the number of aggregators.\n",
                                total_data_size);
                        return;
                    }
                }
                else
                {
                }

                MPI_Gatherv (fd->buffer, pg_size, MPI_BYTE
                            ,aggr_buff, pg_sizes, disp, MPI_BYTE
                            ,0, new_comm);
            }

            // build index appending to any existing index
            adios_build_index_v1 (fd, md->index);

            if (fd->shared_buffer == adios_flag_yes)
            {
                if (!is_aggregator(md->rank))
                {
                    uint64_t var_offset_to_add = 0, attr_offset_to_add = 0;
                    uint64_t var_base_offset = 0, attr_base_offset = 0;

                    // change to relative offset
                    if (md->old_vars_root)
                    {
                        var_base_offset = md->old_vars_root->characteristics [0].offset;
                    }

                    if (md->old_attrs_root)
                    {
                        attr_base_offset = md->old_attrs_root->characteristics [0].offset;
                    }
/*
                    adios_mpi_amr_subtract_offset (var_base_offset
                                                   ,var_base_offset
                                                   ,md->old_vars_root
                                                   ,md->old_attrs_root
                                                   );
*/

                    for (i = 0; i < new_rank; i++)
                    {
                        attr_offset_to_add += pg_sizes[i];
                        var_offset_to_add += pg_sizes[i];
                    }

                    adios_mpi_bgq_add_offset (var_offset_to_add
                                              ,attr_offset_to_add
                                              ,md->old_vars_root
                                              ,md->old_attrs_root
                                              );
                }

                // pg_sizes, disp are no longer needed from this point on.
                free (pg_sizes);
                free (disp);
            }

            // if collective, gather the indexes from the rest and call
            if (md->group_comm != MPI_COMM_NULL)
            {
                // Collect index from all MPI processors
                if (is_aggregator (md->rank))
                {
                    int * index_sizes = malloc (4 * new_group_size);
                    int * index_offsets = malloc (4 * new_group_size);
                    char * recv_buffer = 0;
                    uint32_t size = 0, total_size = 0;
                    int i;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, new_comm
                               );

                    for (i = 0; i < new_group_size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    } 

                    recv_buffer = malloc (total_size);

                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, new_comm
                                );

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    for (i = 1; i < new_group_size; i++)
                    {
                        md->b.buff = recv_buffer + index_offsets [i];
                        md->b.length = index_sizes [i];
                        md->b.offset = 0;

                        adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (md->index
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root, 0
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }

                    md->b.buff = buffer_save;
                    md->b.length = buffer_size_save;
                    md->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);
                }
                else
                {
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,0, md->index
                                         );

                    MPI_Gather (&buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, new_comm
                               );
                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, new_comm
                                );
                }
            }

            // write out indexes in each subfile
            if (is_aggregator (md->rank))
            {
                uint32_t flag = 0;
#if 0
                if (fd->shared_buffer == adios_flag_yes)
                {
                    pthread_join (t, NULL);
                    FREE (aggr_buff);
                }

                MPI_File_get_position (md->fh, (MPI_Offset *)&index_start);
#endif
#if 1
                index_start = total_data_size;
#endif
                adios_write_index_v1 (&buffer, &buffer_size
                                     ,&buffer_offset, index_start
                                     ,md->index
                                     );
//FIXME
                //adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset, flag);
                adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);

                if (fd->shared_buffer == adios_flag_yes)
                {
                    aggr_buff = realloc (aggr_buff, total_data_size + buffer_offset);
                    memcpy (aggr_buff + total_data_size, buffer, buffer_offset); 

                    // Waiting for the subfile to open if pthread is enabled
                    index_start1 = 0;
                    total_data_size1 = total_data_size + buffer_offset;

                    write_thread_data.fh = &md->fh;
                    write_thread_data.base_offset = &index_start1;
                    write_thread_data.aggr_buff = aggr_buff;
                    write_thread_data.total_data_size = &total_data_size1;

                    adios_mpi_bgq_do_write_thread ((void *) &write_thread_data);
                }
#if 0
                adios_mpi_amr_striping_unit_write(
                                  md->fh,
                                  -1,
                                  buffer,
                                  buffer_offset,
                                  md->block_unit);
#endif
            }

            // collect index among aggregators
            if (is_aggregator (md->rank))
            {
                if (md->rank == 0)
                {
                    int * index_sizes = malloc (4 * new_group_size2);
                    int * index_offsets = malloc (4 * new_group_size2);
                    char * recv_buffer = 0;
                    uint32_t size = 0, total_size = 0;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, new_comm2
                               );

                    for (i = 0; i < new_group_size2; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);

                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, new_comm2
                                );

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    for (i = 1; i < new_group_size2; i++)
                    {
                        md->b.buff = recv_buffer + index_offsets [i];
                        md->b.length = index_sizes [i];
                        md->b.offset = 0;

                        adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (md->index
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root, 0
                                             );
                        new_pg_root = 0;
                        new_vars_root = 0;
                        new_attrs_root = 0;
                    }

                    md->b.buff = buffer_save;
                    md->b.length = buffer_size_save;
                    md->b.offset = offset_save;

                    free (recv_buffer);
                    free (index_sizes);
                    free (index_offsets);
                }
                else
                {
                    char * buffer2 = 0;
                    uint64_t buffer_size2 = 0;
                    uint64_t buffer_offset2 = 0;

                    adios_write_index_v1 (&buffer2, &buffer_size2, &buffer_offset2
                                         ,0, md->index
                                         );

                    MPI_Gather (&buffer_size2, 1, MPI_INT
                               ,0, 0, MPI_INT
                               ,0, new_comm2
                               );
                    MPI_Gatherv (buffer2, buffer_size2, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, new_comm2
                                );

                    if (buffer2)
                    {
                        free (buffer2);
                        buffer2 = 0;
                        buffer_size2 = 0;
                        buffer_offset2 = 0;
                    }
                }
            }

            // write out the metadata file from rank 0
            if (md->rank == 0)
            {
                MPI_File m_file;
                char * global_index_buffer = 0;
                uint64_t global_index_buffer_size = 0;
                uint64_t global_index_buffer_offset = 0;
                uint64_t global_index_start = 0;
                uint16_t flag = 0;

                adios_write_index_v1 (&global_index_buffer, &global_index_buffer_size
                                     ,&global_index_buffer_offset, global_index_start
                                     ,md->index
                                     );

                flag |= ADIOS_VERSION_HAVE_SUBFILE;

                adios_write_version_flag_v1 (&global_index_buffer
                                            ,&global_index_buffer_size
                                            ,&global_index_buffer_offset
                                            ,flag
                                            );
/*
                adios_write_version_v1 (&global_index_buffer
                                       ,&global_index_buffer_size
                                       ,&global_index_buffer_offset
                                       );
*/
#if 0
                unlink (fd->name);

                MPI_File_open (MPI_COMM_SELF, fd->name
                              ,MPI_MODE_RDWR | MPI_MODE_CREATE
                              ,MPI_INFO_NULL, &m_file
                              );
#endif
                adios_mpi_bgq_striping_unit_write(
                                  md->mfh,
                                  -1,
                                  global_index_buffer,
                                  global_index_buffer_offset
                                  );

                if (global_index_buffer)
                {
                    free (global_index_buffer);
                    global_index_buffer = 0;
                    global_index_buffer_size = 0;
                    global_index_buffer_offset = 0;
                }
            }

            if (is_aggregator (md->rank))
            {
                FREE (aggr_buff);
            }
            FREE (buffer);
            buffer_size = 0;
            buffer_offset = 0;

            adios_clear_index_v1 (md->index);
            md->index = 0;
            md->g_num_aggregators = 0;
            md->g_color2 = 0;

            FREE (md->subfile_name);
            FREE (md->g_is_aggregator);
            FREE (md->g_offsets);
            FREE (md->open_thread_data);
            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode, 
                    "MPI_BGQ method (AG): Unknown file mode (%d) at close time\n", 
                    fd->mode);
        }
    }

    if (md && md->fh)
        MPI_File_close (&md->fh);

    if (md && md->mfh)
        MPI_File_close (&md->mfh);

    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
    {
        md->group_comm = MPI_COMM_NULL;
    }

    md->fh = 0;
    md->mfh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
#endif
    return;
}

void adios_mpi_bgq_buffer_overflow (struct adios_file_struct * fd, 
                                    struct adios_method_struct * method)
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    adios_error (err_buffer_overflow, 
            "rank %d: MPI_BGQ method only works with complete buffering of data between adios_open() "
            "and adios_close(). Portions of global arrays, that do not fit into the "
            "buffer on some processors will not be written by this method to %s\n", 
            md->rank, fd->name);
}

void adios_mpi_bgq_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    if (md->g_io_type == ADIOS_MPI_BGQ_IO_AG)
    {
        adios_mpi_bgq_ag_close (fd, method);
    }
    else if (md->g_io_type == ADIOS_MPI_BGQ_IO_BG)
    {
        adios_mpi_bgq_bg_close (fd, method);
    }
    else if (md->g_io_type == ADIOS_MPI_BGQ_IO_SIMPLE)
    {
        adios_mpi_bgq_simple_close (fd, method);
    }
    else
    {
        adios_error (err_invalid_write_method, "MPI_BGQ method: unknown I/O type",
                md->g_io_type);
        return;
    }
}

void adios_mpi_bgq_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    adios_free_index_v1 (md->index);

    if (adios_mpi_bgq_initialized)
        adios_mpi_bgq_initialized = 0;
}

void adios_mpi_bgq_end_iteration (struct adios_method_struct * method)
{
}

void adios_mpi_bgq_start_calculation (struct adios_method_struct * method)
{
}

void adios_mpi_bgq_stop_calculation (struct adios_method_struct * method)
{
}
