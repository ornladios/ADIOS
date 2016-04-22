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

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#include "core/adios_timing.h"
#endif


enum ADIOS_MPI_AMR_IO_TYPE
{
    ADIOS_MPI_AMR_IO_NONE = 0,
    ADIOS_MPI_AMR_IO_AG   = 1, // simple all to one aggregation
    ADIOS_MPI_AMR_IO_BG   = 2, // Brigade aggregation
};

static int adios_mpi_amr_initialized = 0;

//#define is_aggregator(rank)  md->g_is_aggregator[rank]
#define is_aggregator(rank)  (md->g_color2 == 0)
#define FREE(v) \
  if (v)        \
  {             \
      free (v); \
      v = 0;    \
  }             \

#define SHIM_FOOTER_SIZE 4
#define ATTR_COUNT_SIZE  2
#define ATTR_LEN_SIZE    8
#define MAX_AGG_BUF      704643072
#define DEFAULT_NUM_OST  672
#define DEFAULT_STRIPE_COUNT 1
#define DEFAULT_STRIPE_SIZE  1024*1024



#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#define START_TIMER(t) adios_timing_go (fd->group->timing_obj, (t) ) 
#else
#define START_TIMER(t) ; 
#endif

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
#define STOP_TIMER(t) adios_timing_stop (fd->group->timing_obj, (t) )
#else
#define STOP_TIMER(t) ;
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

    struct adios_bp_buffer_struct_v1 b;

    struct adios_index_struct_v1 * index;

    uint64_t vars_start;
    uint64_t vars_header_size;

    int * g_is_aggregator;
    int g_num_aggregators;
    int g_have_mdf;
    int g_merging_pgs;
    int g_num_ost;
    int is_local_fs;
    int g_threading;
    int is_color_set; // whether 'color' is set from XML.
    int g_color1;
    int g_color2;
    MPI_Comm g_comm1;
    MPI_Comm g_comm2;
    MPI_Offset * g_offsets;
    int * g_ost_skipping_list;
    pthread_t g_sot;
    pthread_t g_swt; // subfile open thread, metadata file open thread, subfile write thread
    struct adios_MPI_thread_data_open * open_thread_data;
    struct adios_MPI_thread_data_reopen * reopen_thread_data;
    enum ADIOS_MPI_AMR_IO_TYPE g_io_type;
};

struct adios_MPI_thread_data_open
{
    struct adios_MPI_data_struct * md;
    char * parameters;
};

struct adios_MPI_thread_data_reopen
{
    struct adios_MPI_data_struct * md;
    struct adios_file_struct * fd;
};

struct adios_MPI_thread_data_write
{
    MPI_File * fh;
    uint64_t * base_offset;
    void * aggr_buff;
    uint64_t * total_data_size;
};

#if defined(__APPLE__)
#       include <sys/param.h>
#       include <sys/mount.h>
#else
#       include <sys/statfs.h>
#endif

// this should be determined at configure time
//#define ADIOS_LUSTRE

//#ifdef ADIOS_LUSTRE
#include <sys/ioctl.h>
//#include <lustre/lustre_user.h>
//#endif
// from /usr/include/lustre/lustre_user.h
#define LUSTRE_SUPER_MAGIC 0x0BD00BD0
#  define LOV_USER_MAGIC 0x0BD10BD0
#  define LL_IOC_LOV_SETSTRIPE  _IOW ('f', 154, long)
#  define LL_IOC_LOV_GETSTRIPE  _IOW ('f', 155, long)
#define O_LOV_DELAY_CREATE 0100000000

struct lov_user_ost_data {           // per-stripe data structure
        uint64_t l_object_id;        // OST object ID
        uint64_t l_object_gr;        // OST object group (creating MDS number)
        uint32_t l_ost_gen;          // generation of this OST index
        uint32_t l_ost_idx;          // OST index in LOV
} __attribute__((packed));
struct lov_user_md {                 // LOV EA user data (host-endian)
        uint32_t lmm_magic;          // magic number = LOV_USER_MAGIC_V1
        uint32_t lmm_pattern;        // LOV_PATTERN_RAID0, LOV_PATTERN_RAID1
        uint64_t lmm_object_id;      // LOV object ID
        uint64_t lmm_object_gr;      // LOV object group
        uint32_t lmm_stripe_size;    // size of stripe in bytes
        uint16_t lmm_stripe_count;   // num stripes in use for this object
        uint16_t lmm_stripe_offset;  // starting stripe offset in lmm_objects
        struct lov_user_ost_data  lmm_objects[0]; // per-stripe data
} __attribute__((packed));
struct obd_uuid {
        char uuid[40];
};

#ifdef HAVE_FGR
#include "fgr.h"
#endif


int * allocOSTList (int n_ost)
{
    int * ost_list = (int *) malloc (n_ost * sizeof (int));

    if (ost_list == 0)
    {   
        adios_error (err_no_memory, 
                "Can not malloc %d bytes in allocOSTList() in MPI_AMR method\n",
                n_ost*sizeof(int));
        return 0;
    }
    memset (ost_list, 0, n_ost * sizeof (int));

    return ost_list;
}

/* Parse the XML transport parameter to get the list of OST's to skip */
int * parseOSTSkipping (int * ost_list, char * str, int n_ost)
{
    char * p = 0, * dash = 0;
    char n[16];
    int ost_id1, ost_id2, i;

    if (ost_list == 0)
    {
        log_warn ("MPI_AMR method: Pointer ost_list is null.\n");
        return 0;
    }

    p = strtok (str, ",");
    while (p)
    {
        dash = strchr (p, '-');
        if (!dash)
        {
            ost_id1 = atoi (p);
            ost_id2 = ost_id1;
        }
        else
        {
            strncpy (n, p, dash - p);
            n[dash - p] = '\0';
            ost_id1 = atoi (n);
  
            strncpy (n, dash + 1, strlen (dash + 1));
            n[strlen (dash + 1)] = '\0';
            ost_id2 = atoi (n);
        }

        for (i = ost_id1; i <= ost_id2; i++)
        {
            ost_list[i] = 1;
        }
          
        p = strtok (NULL, ",");
    }

    return ost_list;
}

#ifdef HAVE_FGR
int find_myost (MPI_Comm comm)
{
    uint32_t * nids, * osts, myid, ost_id;
    int i, nnids = get_unique_nids (comm, &nids);

    osts = (uint32_t *) malloc (nnids * 4);

    if (fgr_nid2ost (nids, osts, nnids, ATLAS) == true)
    {
/*
int rank;
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
if (rank == 0)
{         
        printf ("nids:");
        for (i = 0; i < nnids; i++)
        {
            printf ("%d:%d ", nids[i], osts[i]);
        }
        printf ("\n");
}
*/
        myid = nid_atoi();
        for (i = 0; i < nnids; i++)
        {
            if (nids[i] == myid)
            {
                ost_id = osts[i];
                break;
            }
        }

        if (i == nnids)
        {
printf ("something is wrong\n");
            // something is wrong
        }

        free (nids);
        free (osts);
#define ATLAS_OFFSET 1008
        return (ost_id >= ATLAS_OFFSET ? ost_id - ATLAS_OFFSET : ost_id);
    }
    else
    {
printf ("something is wrong with FGR\n");
        free (nids);
        free (osts);
        return -1;
    }
}
#endif

static void
//adios_mpi_amr_set_striping_unit(MPI_File fh, char *filename, char *parameters)
adios_mpi_amr_set_striping_unit(struct adios_MPI_data_struct * md, char *parameters)
{
    char * filename = md->subfile_name;
    int err = 0;
    uint64_t striping_unit = 0;
    uint16_t striping_count = 0;
    char     *temp_string, *p_count,*p_size;
    int fd, old_mask, perm, n_ost_skipping, n_ost, n, i, should_striping;
    int random_offset_flag;

    temp_string = (char *) malloc (strlen (parameters) + 1);
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_count = strstr (temp_string, "striping")) )
    {
        char * p = strchr (p_count, '=');
        char * q = strtok (p, ";");
        if (!q)
            should_striping = atoi (q + 1);
        else
            should_striping = atoi (p + 1);
    }
    else
    {
        should_striping = 1;
    }

    if (should_striping == 0)
    {
        return;
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_count = strstr (temp_string, "stripe_count")) )
    {
        char * p = strchr (p_count, '=');
        char * q = strtok (p, ";");
        if (!q)
            striping_count = atoi (q + 1);
        else
            striping_count = atoi (p + 1);
    }
    else
    {
        // By default, set stripe count to 1 to maximize concurrency.
        striping_count = DEFAULT_STRIPE_COUNT;
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_count = strstr (temp_string, "random_offset")) )
    {
        char * p = strchr (p_count, '=');
        char * q = strtok (p, ";");
        if (!q)
            random_offset_flag = atoi (q + 1);
        else
            random_offset_flag = atoi (p + 1);
    }
    else
    {
        // By default, set stripe count to 1 to maximize concurrency.
        random_offset_flag = 0;
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_size = strstr (temp_string, "stripe_size")) )
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            striping_unit = atoi(q + 1);
        else
            striping_unit = atoi(p + 1);
    }
    else
    {
        //stripe size shouldn't matter here. Simply set it to 1 MB here. 
        striping_unit = DEFAULT_STRIPE_SIZE;
    }

    free (temp_string);

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ 0666;

    fd =  open(filename, O_RDONLY | O_CREAT | O_LOV_DELAY_CREATE, perm);
    if (fd != -1) {
        struct lov_user_md lum;
        lum.lmm_magic = LOV_USER_MAGIC;
        lum.lmm_pattern = 0;
        lum.lmm_stripe_size = striping_unit;
        lum.lmm_stripe_count = striping_count;

        // calculate the # of ost's to skip
        n_ost_skipping = 0;
        for (i = 0; i < md->g_num_ost; i++)
        {
            if (md->g_ost_skipping_list[i] == 1)
            {
                n_ost_skipping++;
            }
        }

        // the actual # of ost that can be used
        n_ost = md->g_num_ost - n_ost_skipping;
        if (n_ost <= 0)
        {
            log_warn ("MPI_AMR method: No OST to use. Set num_ost=NNN in the adios config xml file.\n");
            return;
        }

        i = 0;
        n = 0;
        while (i < md->g_num_ost)
        {
            if (md->g_ost_skipping_list[i] == 0)
            {
                n++;
                if (n - 1 == md->g_color1 % n_ost)
                    break;
            }
            
            i++;
        }

#ifdef HAVE_FGR
       int ost_id = find_myost (md->g_comm2);
       if (ost_id >= 0)
       {
           lum.lmm_stripe_offset = ost_id;
       }
       else
       {
printf ("why is here\n");
           lum.lmm_stripe_offset = (random_offset_flag ? -1 : i);
       }
#else
        lum.lmm_stripe_offset = (random_offset_flag ? -1 : i);
#endif
        ioctl (fd, LL_IOC_LOV_SETSTRIPE
              ,(void *) &lum
              );

        if (err == 0 && lum.lmm_stripe_size > 0) {
            striping_unit = lum.lmm_stripe_size;
        }

        close(fd);
    }
    else
    {
        log_warn("MPI_AMR method: open to set lustre striping failed on file %s %s rank = %d.\n",filename,strerror(errno), md->rank);
    }
}

static void
adios_mpi_amr_set_have_mdf (char * parameters, struct adios_MPI_data_struct * md)
{
    char *temp_string, *p_size;

    temp_string = (char *) malloc (strlen (parameters) + 1);
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_size = strstr (temp_string, "have_metadata_file")) )
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
}

static void
adios_mpi_amr_set_aggregation_parameters(char * parameters, struct adios_MPI_data_struct * md)
{
    int i, aggr_group_size, remain, index;
    int nproc = md->size, rank = md->rank;
    char *temp_string, *p_size;

    temp_string = (char *) malloc (strlen (parameters) + 1);

    // set up the number of OST to use
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_size = strstr (temp_string, "num_ost")) )
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            md->g_num_ost = atoi(q + 1);
        else
            md->g_num_ost = atoi(p + 1);
    }
    else
    {
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_size = strstr (temp_string, "local-fs")) )
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            md->is_local_fs = atoi(q + 1);
        else
            md->is_local_fs = atoi(p + 1);
    }
    else
    {
        md->is_local_fs = 0;
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    // set up # of aggregators
    if ( (p_size = strstr (temp_string, "num_aggregators")) )
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            md->g_num_aggregators = atoi(q + 1);
        else
            md->g_num_aggregators = atoi(p + 1);
    }
    else
    {
        if (nproc <= md->g_num_ost)
        {
            md->g_num_aggregators = nproc;
        }
        else
        {
            md->g_num_aggregators = md->g_num_ost;
        }
    }

    // Get 'color' parameter. If 'color' is set,
    // the num_aggregators will be disregarded.
    // The actual # of aggregators will be caculated
    // according to color. 
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_size = strstr (temp_string, "color")) )
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");

        md->is_color_set = 1;
        if (!q)
            md->g_color1 = atoi (q + 1);
        else
            md->g_color1 = atoi (p + 1);
    }
    else
    {
        // by default, use BG
        md->g_io_type = ADIOS_MPI_AMR_IO_BG;
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_size = strstr (temp_string, "have_metadata_file")) )
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

    // set up whether to thread IO ops
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_size = strstr (temp_string, "threading")) )
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            md->g_threading = atoi(q + 1);
        else
            md->g_threading = atoi(p + 1);
    }
    else
    {
        // by default, threading is disabled.
        md->g_threading = 0;
    }

    // set up which ost's to skip
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    md->g_ost_skipping_list = allocOSTList (md->g_num_ost);

    if ( (p_size = strstr (temp_string, "osts_to_skip")) )
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");

        if (!q)
            md->g_ost_skipping_list = parseOSTSkipping (md->g_ost_skipping_list, q + 1, md->g_num_ost);
        else
            md->g_ost_skipping_list = parseOSTSkipping (md->g_ost_skipping_list, p + 1, md->g_num_ost);
    }

    // set up which ost's to skip
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if ( (p_size = strstr (temp_string, "aggregation_type")) )
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
        // by default, use BG
        md->g_io_type = ADIOS_MPI_AMR_IO_BG;
    }

    free (temp_string);

    if (md->g_num_aggregators > nproc || md->g_num_aggregators <= 0)
    {
        md->g_num_aggregators = nproc;  //no aggregation
    }

    md->g_is_aggregator = (int *) malloc (nproc * sizeof(int));
    if (md->g_is_aggregator == 0)
    {
        adios_error (err_no_memory, "Can not malloc %d bytes in MPI_AMR method, adios_mpi_amr_set_aggregation_parameters()\n", nproc*sizeof(int));
        return;
    }
    memset (md->g_is_aggregator, 0, nproc * sizeof(int));

    if (!md->is_color_set)
    {
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

        if (remain == 0)
        {
            md->g_color1 = rank / aggr_group_size;
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

        MPI_Comm_split (md->group_comm, md->g_color1, md->rank, &md->g_comm1);
        MPI_Comm_split (md->group_comm, md->g_color2, md->rank, &md->g_comm2);
    }
    else // if color is set
    {
        MPI_Comm_split (md->group_comm, md->g_color1, md->rank, &md->g_comm1);
        MPI_Comm_rank (md->g_comm1, &md->g_color2);
    }
}

/*
static void adios_mpi_amr_buffer_write (char ** buffer, uint64_t * buffer_size
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
            adios_error (err_no_memory, "Cannot allocate memory in adios_mpi_amr_buffer_write.  "
                             "Requested: %llu\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}
*/

static uint64_t
adios_mpi_amr_striping_unit_write(MPI_File   fh
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
        return 0;

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

struct adios_var_struct * adios_mpi_amr_copy_var (struct adios_var_struct * v)
{
    struct adios_var_struct * v_new = (struct adios_var_struct *) 
                            malloc (sizeof (struct adios_var_struct));
    if (v_new == 0)
    {
        adios_error (err_no_memory, "MPI_AMR method: Cannot allocate %d bytes "
            "to duplicate variable structure in adios_mpi_amr_copy_var()\n",
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

void adios_mpi_amr_append_var (struct adios_file_struct * fd, struct adios_var_struct * v)
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

void adios_mpi_amr_add_offset (uint64_t pg_offset_to_add, 
                               uint64_t var_offset_to_add,
                               uint64_t attr_offset_to_add,
                               struct adios_index_struct_v1 * index
                              )
{
    struct adios_index_process_group_struct_v1 *pg_root = index->pg_root;
    struct adios_index_var_struct_v1 * vars_root = index->vars_root;
    struct adios_index_attribute_struct_v1 * attrs_root = index->attrs_root;
    while (pg_root)
    {
        pg_root->offset_in_file += pg_offset_to_add;
        pg_root = pg_root->next;
    }

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

void adios_mpi_amr_subtract_offset (uint64_t var_offset_to_subtract
                                    ,uint64_t attr_offset_to_subtract
                                    ,struct adios_index_struct_v1 * index
                                    )
{
    struct adios_index_var_struct_v1 * vars_root = index->vars_root;
    struct adios_index_attribute_struct_v1 * attrs_root = index->attrs_root;
    while (vars_root)
    {
        vars_root->characteristics [0].offset -= var_offset_to_subtract;
        vars_root->characteristics [0].payload_offset -= var_offset_to_subtract;
        vars_root = vars_root->next;
    }

    while (attrs_root)
    {
        attrs_root->characteristics [0].offset -= attr_offset_to_subtract;
        attrs_root->characteristics [0].payload_offset -= attr_offset_to_subtract;
        attrs_root = attrs_root->next;
    }
}


void adios_mpi_amr_build_global_index_v1 (char * fname
                                         ,struct adios_index_struct_v1 * index
                                          )
{
    //struct adios_index_process_group_struct_v1 * pg_root = index->pg_root;
    struct adios_index_var_struct_v1 * vars_root = index->vars_root;
    struct adios_index_attribute_struct_v1 * attrs_root = index->attrs_root;
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


void * adios_mpi_amr_do_mkdir (char * path)
{
    // 4 bytes for ".dir" 
    char * dir_name = malloc (strlen (path) + 4 + 1);
    sprintf (dir_name, "%s%s", path, ".dir");
    
    mkdir (dir_name, S_IRWXU | S_IRWXG);
  
    free (dir_name);

    return NULL;
}

void * adios_mpi_amr_do_open_thread (void * param)
{
    struct adios_MPI_thread_data_open * td = (struct adios_MPI_thread_data_open *) param;
    int err;

    unlink (td->md->subfile_name);
    if (td->parameters)
    {
        adios_mpi_amr_set_striping_unit (td->md, td->parameters);

    }

    err = MPI_File_open (MPI_COMM_SELF, td->md->subfile_name
                        ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                        ,MPI_INFO_NULL
                        ,&td->md->fh
                        );

    if (err != MPI_SUCCESS)
    {
        char e [MPI_MAX_ERROR_STRING];
        int len = 0;
        memset (e, 0, MPI_MAX_ERROR_STRING);
        MPI_Error_string (err, e, &len);
        adios_error (err_file_open_error,
                     "MPI_AMR method: MPI open failed for %s: '%s'\n",
                     td->md->subfile_name, e);
    }

    return NULL;
}

// reopen a subfile for append and read/build the existing index
void * adios_mpi_amr_do_reopen_thread (void * param)
{
    struct adios_MPI_thread_data_reopen * td = (struct adios_MPI_thread_data_reopen *) param;
    struct adios_MPI_data_struct * md = td->md;
    struct adios_file_struct * fd = td->fd;
    int err;

    // open for read/write because we keep it open for all writing later
    err = MPI_File_open (MPI_COMM_SELF, td->md->subfile_name, 
                         MPI_MODE_RDWR, MPI_INFO_NULL, &td->md->fh);

    if (err == MPI_SUCCESS) 
    {
        // read in the index
        MPI_Offset file_size;
        MPI_File_get_size (md->fh, &file_size);
        md->b.file_size = file_size;

        adios_init_buffer_read_version (&md->b);
        MPI_File_seek (md->fh, md->b.file_size - md->b.length, MPI_SEEK_SET);
        MPI_File_read (md->fh, md->b.buff, md->b.length, MPI_BYTE, &md->status);
        adios_parse_version (&md->b, &md->b.version);

        adios_init_buffer_read_index_offsets (&md->b);
        // already in the buffer
        adios_parse_index_offsets_v1 (&md->b);

        adios_init_buffer_read_process_group_index (&md->b);
        MPI_File_seek (md->fh, md->b.pg_index_offset, MPI_SEEK_SET);
        MPI_File_read (md->fh, md->b.buff, md->b.pg_size, MPI_BYTE, &md->status);

        adios_parse_process_group_index_v1 (&md->b, &md->index->pg_root, &md->index->pg_tail);

        // find the largest time index so we can append properly
        struct adios_index_process_group_struct_v1 * p;
        uint32_t max_time_index = 0;
        p = md->index->pg_root;
        while (p)
        {
            if (p->time_index > max_time_index)
                max_time_index = p->time_index;
            p = p->next;
        }
        fd->group->time_index = ++max_time_index;

        adios_init_buffer_read_vars_index (&md->b);
        MPI_File_seek (md->fh, md->b.vars_index_offset, MPI_SEEK_SET);
        MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE, &md->status);
        adios_parse_vars_index_v1 (&md->b, &md->index->vars_root, 
                                    md->index->hashtbl_vars,
                                   &md->index->vars_tail);

        adios_init_buffer_read_attributes_index (&md->b);
        MPI_File_seek (md->fh, md->b.attrs_index_offset, MPI_SEEK_SET);
        MPI_File_read (md->fh, md->b.buff, md->b.attrs_size, MPI_BYTE, &md->status);
        adios_parse_attributes_index_v1 (&md->b, &md->index->attrs_root);

        // remember the end of the last PG in file. The new PG written now
        // will have an offset updated from here. 
        // md->b.end_of_pgs points to this position

       // printf ("rank %d: APPEND: end_of_pgs=%llu bytes_written=%llu  pg_index_offset=%llu\n", 
       //         md->rank, md->b.end_of_pgs, fd->bytes_written, md->b.pg_index_offset);
    } 
    else 
    {
        // create as write-only now
        err = MPI_File_open (MPI_COMM_SELF, td->md->subfile_name
                ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                ,MPI_INFO_NULL
                ,&td->md->fh
                );

        if (err != MPI_SUCCESS)
        {
            char e [MPI_MAX_ERROR_STRING];
            int len = 0;
            memset (e, 0, MPI_MAX_ERROR_STRING);
            MPI_Error_string (err, e, &len);
            adios_error (err_file_open_error,
                    "MPI_AMR method: MPI open failed for %s: '%s'\n",
                    td->md->subfile_name, e);
            td->md->fh = NULL;
        }
        /* FIXME: how does this error propagate back to the user and later adios calls? */

        md->b.file_size = 0;
    }

    return NULL;
}

void * adios_mpi_amr_do_write_thread (void * param)
{
    struct adios_MPI_thread_data_write * td = (struct adios_MPI_thread_data_write *) param;

    uint64_t count = adios_mpi_amr_striping_unit_write(
                               *(td->fh)
                              ,*(td->base_offset)
                              ,td->aggr_buff
                              ,*(td->total_data_size)
                              );

    if (count != *(td->total_data_size))
    {
        adios_error (err_unspecified, "Error in adios_mpi_amr_striping_unit_write(). "
            "count = %llu != thread's total_data_size = %llu\n",
            count, td->total_data_size);
    }

    return NULL;
}

void adios_mpi_amr_init (const PairStruct * parameters
                         ,struct adios_method_struct * method
                         )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    if (!adios_mpi_amr_initialized)
    {
        adios_mpi_amr_initialized = 1;
    }

    method->method_data = malloc (sizeof (struct adios_MPI_data_struct));
    md = (struct adios_MPI_data_struct *) method->method_data;
    md->fh = 0;
    md->mfh = 0;
    md->subfile_name = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->rank = 0;
    md->size = 0;
    md->group_comm = method->init_comm; //unused, adios_open sets current comm
    md->index = adios_alloc_index_v1(1); // with hashtables
    md->vars_start = 0;
    md->vars_header_size = 0;

    md->g_is_aggregator = 0;
    md->g_num_aggregators = 0;
    md->g_have_mdf = 1;
    md->g_merging_pgs = 0;
    md->g_num_ost = 0;
    md->is_local_fs = 0;
    md->g_threading = 0;
    md->is_color_set = 0;
    md->g_color1 = 0;
    md->g_color2 = 0;
    md->g_offsets = 0;
    md->g_ost_skipping_list = 0;
    md->open_thread_data = 0;
    md->reopen_thread_data = 0;
    md->g_io_type = ADIOS_MPI_AMR_IO_BG;

    adios_buffer_struct_init (&md->b);

#ifdef HAVE_FGR
    if (fgr_init (0) == false)
    {
        adios_error (err_fgr, "fgr_init() error\n");
    }
#endif
}


static char * get_subfile_name (char *base_path, char *filename, int color)
{
    char *ch, *name_no_path, *subfilename;
    // Check if fd->name contains path
    if ( (ch = strrchr (filename, '/')) )
    {
        name_no_path = malloc (strlen (ch + 1) + 1); 
        strcpy (name_no_path, ch + 1); 
    }
    else
    {
        name_no_path = malloc (strlen (filename) + 1);
        strcpy (name_no_path, filename);
    }

    subfilename = malloc (strlen (base_path) + strlen (filename) + 5 + strlen (name_no_path) + 1 + 10 + 1);
    // create the subfile name, e.g. restart.bp.1
    // 1 for '.' + 10 for subfile index + 1 for '\0'
    sprintf (subfilename, "%s%s%s%s.%d", base_path, filename, ".dir/", name_no_path, color);
    free (name_no_path);
    return subfilename;
}


#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS
// Indices for the timer object
int ADIOS_TIMER_COMM        = ADIOS_TIMING_MAX_USER_TIMERS + 0;
int ADIOS_TIMER_IO          = ADIOS_TIMING_MAX_USER_TIMERS + 1;
int ADIOS_TIMER_LOCALMD     = ADIOS_TIMING_MAX_USER_TIMERS + 2;
int ADIOS_TIMER_GLOBALMD    = ADIOS_TIMING_MAX_USER_TIMERS + 3;
int ADIOS_TIMER_AD_OPEN     = ADIOS_TIMING_MAX_USER_TIMERS + 4;
int ADIOS_TIMER_AD_WRITE    = ADIOS_TIMING_MAX_USER_TIMERS + 5;
int ADIOS_TIMER_AD_OVERFLOW = ADIOS_TIMING_MAX_USER_TIMERS + 6;
int ADIOS_TIMER_AD_CLOSE    = ADIOS_TIMING_MAX_USER_TIMERS + 7;
#endif


int adios_mpi_amr_open (struct adios_file_struct * fd
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

    fd->group->process_id = md->rank;
    
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

    // need to dealloc/realloc buffer because of supporting append mode
    // b.length and others won't be re-initialized otherwise.
    adios_buffer_struct_clear (&md->b);


    char * name;
    //int sig;    // used for coordinating the MPI_File_open

    START_TIMER (ADIOS_TIMER_AD_OPEN);

    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            adios_error (err_invalid_file_mode, "MPI_AMR method: Read mode is not supported.\n");
            break;
        }

        case adios_mode_write:
        {
            if (md->rank == 0)
            {
                struct lov_user_md lum;
                int f;

                // open metadata file
                unlink (name);

                adios_mpi_amr_set_have_mdf (method->parameters, md);
                if (md->g_have_mdf)
                {
                    f = open(name, O_CREAT | O_RDWR | O_LOV_DELAY_CREATE, 0644);
                    if (f == -1)
                    {
//                        adios_error (err_file_open_error,"MPI_AMR method: open() failed: %s\n", strerror(errno));
                        adios_error (err_file_open_error,"MPI_AMR method: open() failed: %s\n", name);
                        return -1;
                    }

                    lum.lmm_magic = LOV_USER_MAGIC;
                    lum.lmm_pattern = 0;
                    lum.lmm_stripe_size = DEFAULT_STRIPE_SIZE;
                    lum.lmm_stripe_count = 1;
                    lum.lmm_stripe_offset = -1;

                    ioctl (f, LL_IOC_LOV_SETSTRIPE ,(void *) &lum);
#ifdef HAVE_LUSTRE
                    struct obd_uuid uuids[1024];
                    int rc;
                    md->g_num_ost = 1024;
                    rc = llapi_lov_get_uuids(f, uuids, &md->g_num_ost);
                    if (rc != 0)
                    {
                        log_warn ("MPI_AMR method: Lustre get uuids failed after creating the file: %s\n" ,
                                  strerror(errno));
                    }

#endif 
                    close (f);

                    MPI_File_open (MPI_COMM_SELF, name
                                  ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                  ,MPI_INFO_NULL
                                  ,&md->mfh
                                  );
                }

//                adios_mpi_amr_do_mkdir (name);
            }

            MPI_Bcast (&md->g_num_ost, 1, MPI_INT, 0, md->group_comm);
          
            adios_mpi_amr_set_aggregation_parameters (method->parameters, md);

            if (is_aggregator (md->rank))
            {
                if (md->is_local_fs)
                {
                    adios_mpi_amr_do_mkdir (name);
                }
                else
                {
                    if (md->rank == 0)
                    {
                        adios_mpi_amr_do_mkdir (name);
                    }
                }
            
                MPI_Barrier (md->g_comm2);
            }


            md->subfile_name = get_subfile_name (method->base_path, fd->name, md->g_color1);
            fd->subfile_index = (uint32_t)md->g_color1;


            if (is_aggregator(md->rank))
            {
                // open subfiles
                md->open_thread_data = (struct adios_MPI_thread_data_open *) malloc (sizeof (struct adios_MPI_thread_data_open));
                md->open_thread_data->md = md;
                md->open_thread_data->parameters = method->parameters;

                if (md->g_threading)
                {
                    pthread_create (&md->g_sot, NULL
                            ,adios_mpi_amr_do_open_thread
                            ,(void *) md->open_thread_data
                            );
                }
                else
                {
                    adios_mpi_amr_do_open_thread ((void *) md->open_thread_data);
                }
            }

            /*log_debug ("rank %d: WRITE: end_of_pgs=%llu bytes_written=%llu pg_start_in_file=%llu"
                      "  write_size_bytes=%llu  pg_index_offset=%llu\n", md->rank,
                      md->b.end_of_pgs, fd->bytes_written, fd->pg_start_in_file, 
                      fd->write_size_bytes, md->b.pg_index_offset);*/
            break;
        }

        case adios_mode_append:
        case adios_mode_update:
        {
            // rank 0 opens the global metadata (if that exists)
            // we throw away the content because it's a replica of metadata present in subfiles
            // so the index merging procedure will result in a complete metadata again
            if (md->rank == 0)
            {
                int f;
                md->g_num_ost = 1024; // default num_ost, maybe updated below

                // open metadata file
                adios_mpi_amr_set_have_mdf (method->parameters, md);
                if (md->g_have_mdf)
                {
                    f = open(name, O_RDWR, 0644);
                    if (f == -1)
                    {
                        adios_error (err_file_open_error,"MPI_AMR method: open() failed at append: %s\n", name);
                        return -1;
                    }

#ifdef HAVE_LUSTRE
                    struct obd_uuid uuids[1024];
                    int rc;
                    rc = llapi_lov_get_uuids(f, uuids, &md->g_num_ost);
                    if (rc != 0)
                    {
                        log_warn ("MPI_AMR method: Lustre get uuids failed after opening the file: %s\n" ,
                                  strerror(errno));
                    }

#endif 
                    close (f);
                    MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_WRONLY ,MPI_INFO_NULL ,&md->mfh);
                }
            }

            MPI_Bcast (&md->g_num_ost, 1, MPI_INT, 0, md->group_comm);

            adios_mpi_amr_set_aggregation_parameters (method->parameters, md); // is_aggregator() works after this

            // Every aggregator opens a subfile and reads in the metadata
            md->subfile_name = get_subfile_name (method->base_path, fd->name, md->g_color1);
            fd->subfile_index = (uint32_t)md->g_color1;

            if (is_aggregator(md->rank))
            {
                // open subfiles
                md->reopen_thread_data = (struct adios_MPI_thread_data_reopen *) 
                    malloc (sizeof (struct adios_MPI_thread_data_reopen));
                md->reopen_thread_data->md = md;
                md->reopen_thread_data->fd = fd;

                // there is a MPI_Bcast in the function below, so this cannot be threaded as is now
                /*
                   if (md->g_threading)
                   {
                   pthread_create (&md->g_sot, NULL,
                   adios_mpi_amr_do_reopen_thread,
                   (void *) md->reopen_thread_data);
                   }
                   else
                 */
                {
                    //log_debug ("rank %d: APPEND: reopen subfile...\n", md->rank);
                    adios_mpi_amr_do_reopen_thread ((void *) md->reopen_thread_data);
                }

                MPI_Bcast (&fd->group->time_index, 1, MPI_INT, 0, md->g_comm1);
                MPI_Bcast (&md->b.pg_index_offset, 1, MPI_LONG_LONG, 0, md->g_comm1);
            } 
            else //non-aggregators
            {
                MPI_Bcast (&fd->group->time_index, 1, MPI_INT, 0, md->g_comm1);
                MPI_Bcast (&md->b.pg_index_offset, 1, MPI_LONG_LONG, 0, md->g_comm1);
            }

            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode, "MPI_AMR method: Unknown file mode requested: %d\n", fd->mode);
            free (name);
            return adios_flag_no;
        }
    }

    free (name);

    STOP_TIMER (ADIOS_TIMER_AD_OPEN);
    return 1;
}


enum BUFFERING_STRATEGY adios_mpi_amr_should_buffer (struct adios_file_struct * fd
                                                    ,struct adios_method_struct * method
                                                    )
{
    return stop_on_overflow;
}

void adios_mpi_amr_write (struct adios_file_struct * fd
                         ,struct adios_var_struct * v
                         ,const void * data
                         ,struct adios_method_struct * method
                         )
{
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

void adios_mpi_amr_get_write_buffer (struct adios_file_struct * fd
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
                    "MPI_AMR method: Out of memory allocating %llu bytes for variable %s\n",
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
                "MPI_AMR method: OVERFLOW: Cannot allocate requested buffer of %llu "
                "bytes for %s. Allowed max size is %llu\n",
                *size, v->name, mem_allowed);
        *size = 0;
        *buffer = 0;
    }
}

void adios_mpi_amr_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )
{
    v->data = buffer;
    v->data_size = buffer_size;
}


/*static
uint32_t adios_mpi_amr_calculate_attributes_size (struct adios_file_struct * fd)
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
}*/


/* Help routine to send data size greater than 2 GB */
int adios_MPI_Send(void *buf, uint64_t count, int dest, int tag,
             MPI_Comm comm)
{
    while (count > INT32_MAX)
    {
        MPI_Send (buf, INT32_MAX, MPI_BYTE, dest, tag, comm);
        count -= INT32_MAX;
        buf += INT32_MAX;
    }

    if (count)
    {
        int temp_count = (int) count;
        MPI_Send (buf, temp_count, MPI_BYTE, dest, tag, comm);
    }

    return 0;
}


/* Help routine to receive data size greater than 2 GB */
int adios_MPI_Recv(void *buf, uint64_t count, int source,
                    int tag, MPI_Comm comm, MPI_Status *status)
{
    while (count > INT32_MAX)
    {
        MPI_Recv (buf, INT32_MAX, MPI_BYTE, source, tag, comm, status);
        count -= INT32_MAX;
        buf += INT32_MAX;
    }

    if (count)
    {
        int temp_count = (int) count;
        MPI_Recv (buf, temp_count, MPI_BYTE, source, tag, comm, status);
    }

    return 0;

}

void adios_mpi_amr_bg_close (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            adios_error (err_invalid_file_mode, 
                    "Only \"w\" mode is supported by MPI_AMR Brigade IO\n");
            break;
        }
        case adios_mode_update:
        case adios_mode_append:
        case adios_mode_write:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start1;
            uint64_t * pg_sizes = 0, * disp = 0;
            void * aggr_buff = 0, * recv_buff = 0;
            struct adios_MPI_thread_data_write write_thread_data;
            int i, new_rank, new_group_size, new_rank2, new_group_size2;
            uint64_t max_data_size = 0, total_data_size = 0, total_data_size1 = 0;
            START_TIMER (ADIOS_TIMER_COMM);
            //MPI_Comm_split (md->group_comm, md->g_color1, md->rank, &md->g_comm1);
            MPI_Comm_rank (md->g_comm1, &new_rank);
            MPI_Comm_size (md->g_comm1, &new_group_size);

            //MPI_Comm_split (md->group_comm, md->g_color2, md->rank, &md->g_comm2);
            MPI_Comm_rank (md->g_comm2, &new_rank2);
            MPI_Comm_size (md->g_comm2, &new_group_size2);
            STOP_TIMER (ADIOS_TIMER_COMM);

            // if not merge PG's on the aggregator side
            if (!md->g_merging_pgs)
            {
                //printf ("do not merge pg\n");
                uint64_t pg_size;
                MPI_Status status;

                pg_size = fd->bytes_written;
                pg_sizes = (uint64_t *) malloc (new_group_size * 8);
                disp = (uint64_t *) malloc (new_group_size * 8);
                if (pg_sizes == 0 || disp == 0)
                {
                    adios_error (err_no_memory, "MPI_AMR method: Cannot allocate memory "
                                "for merging process blocks (mpi_amr_bg_close)\n");
                    return;
                }

                START_TIMER (ADIOS_TIMER_COMM);
                MPI_Allgather (&pg_size, 1, MPI_UNSIGNED_LONG_LONG
                              ,pg_sizes, 1, MPI_UNSIGNED_LONG_LONG
                              ,md->g_comm1);
                STOP_TIMER (ADIOS_TIMER_COMM);

                disp[0] = 0;
                max_data_size = pg_size;

                for (i = 1; i < new_group_size; i++)
                {
                    disp[i] = disp[i - 1] + pg_sizes[i - 1];
                    max_data_size = (pg_sizes[i] > max_data_size) ? pg_sizes[i] : max_data_size;
                }

                if (is_aggregator (md->rank))
                {
                    aggr_buff = malloc (max_data_size);
                    recv_buff = malloc (max_data_size);
                    if (aggr_buff == 0 || recv_buff == 0)
                    {
                        adios_error (err_no_memory, "MPI_AMR method (with brigade strategy): Cannot allocate "
                                    "2 x %lu bytes for aggregation buffers. "
                                    "An aggregator process needs a buffer to hold one process' output for writing, "
                                    "while it needs another buffer to concurrently receive another process' output for "
                                    "subsequent writing.\n", 
                                    max_data_size);
                        return;
                    }
                }
                else
                {
                    recv_buff = malloc (max_data_size);
                    if (recv_buff == 0)
                    {
                        adios_error (err_no_memory, "MPI_AMR method (with brigade strategy): Cannot allocate "
                                    "%lu bytes for receive buffer in a non-aggregator process. "
                                    "This method needs an extra buffer in every process to pass data along "
                                    "towards the aggregator.\n", 
                                    max_data_size);
                        return;
                    }
                }

                total_data_size = disp[new_group_size - 1]
                                + pg_sizes[new_group_size - 1];

                if (is_aggregator (md->rank))
                {
                    if (md->g_threading)
                    {
                        pthread_join (md->g_sot, NULL);
                    }

                    index_start1 = md->b.pg_index_offset; // starting point to write data at this moment
                    for (i = 0; i < new_group_size; i++)
                    {
                        if (i + 1 < new_group_size)
                        {
                            START_TIMER (ADIOS_TIMER_COMM);
                            adios_MPI_Recv (recv_buff, pg_sizes[i + 1], new_rank + 1
                                      ,0, md->g_comm1, &status);
                            STOP_TIMER (ADIOS_TIMER_COMM);
                        }

                        write_thread_data.fh = &md->fh;
                        write_thread_data.base_offset = &index_start1;
                        write_thread_data.aggr_buff = (i == 0) ? fd->buffer : aggr_buff;
                        write_thread_data.total_data_size = &pg_sizes[i];

                        //printf ("rank %d: Write PG to subfile %d, offset=%llu, size=%u\n", md->rank,
                        //       fd->subfile_index, *write_thread_data.base_offset, pg_sizes[i]);

                        // This write call is not threaded
                        START_TIMER (ADIOS_TIMER_IO);
                        adios_mpi_amr_do_write_thread ((void *) &write_thread_data);
                        STOP_TIMER (ADIOS_TIMER_IO);

                        index_start1 += pg_sizes[i];

                        if (i + 1 < new_group_size)
                        {
                            memcpy (aggr_buff, recv_buff, pg_sizes[i + 1]);
                        }
                    }
                }
                else
                {
                    if (new_rank == new_group_size - 1)
                    {
                        START_TIMER (ADIOS_TIMER_COMM);
                        adios_MPI_Send (fd->buffer, pg_size, new_rank - 1
                                 ,0, md->g_comm1);
                        STOP_TIMER (ADIOS_TIMER_COMM);
                    }
                    else
                    {
                        for (i = new_rank + 1; i < new_group_size; i++)
                        {
                            START_TIMER (ADIOS_TIMER_COMM);
                            // Recv data from upstream rank
                            adios_MPI_Recv (recv_buff, pg_sizes[i], new_rank + 1
                                      ,0, md->g_comm1, &status);

                            if (i == new_rank + 1)
                                // Send my data to downstream rank
                                adios_MPI_Send (fd->buffer, pg_size, new_rank - 1
                                         ,0, md->g_comm1);

                            //MPI_Wait (&request, &status);
                            // Send it to downstream rank
                            adios_MPI_Send (recv_buff, pg_sizes[i], new_rank - 1
                                     ,0, md->g_comm1);
                            STOP_TIMER (ADIOS_TIMER_COMM);
                        }
                    }
                }

                FREE (aggr_buff);
                FREE (recv_buff);
            }
            else 
            {
                // Merge PG's on the aggregator side
                log_warn ("MPI_AMR method (BG): Merging process blocks is not supported yet\n");
            }


            /* Before building index, determine the actual offset where this PG starts */
            fd->current_pg->pg_start_in_file = md->b.pg_index_offset; // aggregator's starting offset (!0 on append)
            if (!md->g_merging_pgs)
            {
                for (i = 0; i < new_rank; i++)
                {
                    fd->current_pg->pg_start_in_file += pg_sizes[i];
                }

                FREE (pg_sizes);
                FREE (disp);
            }

            // build index appending to any existing index
            //printf ("rank %d: Build index with PG offset %llu\n", md->rank, fd->current_pg->pg_start_in_file);
            adios_build_index_v1 (fd, md->index);

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

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, md->g_comm1
                               );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    for (i = 0; i < new_group_size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, md->g_comm1
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);

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
                        if (md->g_merging_pgs)
                            new_pg_root = 0;

                        adios_merge_index_v1 (md->index, new_pg_root, 
                                              new_vars_root, new_attrs_root, 0);
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
                                         ,0, md->index);

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->g_comm1
                               );
                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->g_comm1
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);
                }
            }

            // write out indexes in each subfile
            if (is_aggregator (md->rank))
            {

                uint32_t flag = 0;
                uint64_t index_start; // = end of PGs, where the index starts
                index_start  = md->b.pg_index_offset; // = end of PGs of previous timesteps
                index_start += total_data_size; //old index start before append + currently written PGs
                /*DEBUG*/
                /*log_debug ("rank %d: write index start=%llu  pg_index_offset=%llu  total_data_size=%u\n", 
                        md->rank, index_start, md->b.pg_index_offset, total_data_size);

                struct adios_index_process_group_struct_v1 *pg_root = md->index->pg_root;
                i=0;
                while (pg_root) {
                    log_debug ("rank %d: pg %d offset=%llu\n", md->rank, i, pg_root->offset_in_file);
                    pg_root = pg_root->next;
                    i++;
                }*/

                adios_write_index_v1 (&buffer, &buffer_size
                                     ,&buffer_offset, index_start
                                     ,md->index);

                adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);

                index_start = -1;
                total_data_size1 = buffer_offset;

                write_thread_data.fh = &md->fh;
                write_thread_data.base_offset = &index_start;
                write_thread_data.aggr_buff = buffer;
                write_thread_data.total_data_size = &total_data_size1;

                if (md->g_threading)
                {
                    pthread_create (&md->g_swt, NULL
                            ,adios_mpi_amr_do_write_thread
                            ,(void *) &write_thread_data
                            );
                }
                else
                {
                    START_TIMER (ADIOS_TIMER_IO);
                    adios_mpi_amr_do_write_thread ((void *) &write_thread_data); 
                    STOP_TIMER (ADIOS_TIMER_IO);
                }
            }

            // collect index among aggregators
            if (md->g_have_mdf)
            {
                if (is_aggregator (md->rank))
                {
                    if (md->rank == 0)
                    {
                        int * index_sizes = malloc (4 * new_group_size2);
                        int * index_offsets = malloc (4 * new_group_size2);
                        char * recv_buffer = 0;
                        uint32_t size = 0, total_size = 0;

                        START_TIMER (ADIOS_TIMER_COMM);
                        MPI_Gather (&size, 1, MPI_INT
                                   ,index_sizes, 1, MPI_INT
                                   ,0, md->g_comm2
                                   );
                        STOP_TIMER (ADIOS_TIMER_COMM);

                        for (i = 0; i < new_group_size2; i++)
                        {
                            index_offsets [i] = total_size;
                            total_size += index_sizes [i];
                        }

                        recv_buffer = malloc (total_size);

                        START_TIMER (ADIOS_TIMER_COMM);
                        MPI_Gatherv (&size, 0, MPI_BYTE
                                    ,recv_buffer, index_sizes, index_offsets
                                    ,MPI_BYTE, 0, md->g_comm2
                                    );
                        STOP_TIMER (ADIOS_TIMER_COMM);

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

                            // global index would become unsorted on main aggregator during merging 
                            // so sort timesteps if appending
                            adios_merge_index_v1 (md->index, new_pg_root, 
                                                  new_vars_root, new_attrs_root,
                                                  (fd->mode == adios_mode_append));
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
                                             ,0, md->index);
 
                        START_TIMER (ADIOS_TIMER_COMM);
                        MPI_Gather (&buffer_size2, 1, MPI_INT
                                   ,0, 0, MPI_INT
                                   ,0, md->g_comm2
                                   );
                        MPI_Gatherv (buffer2, buffer_size2, MPI_BYTE
                                    ,0, 0, 0, MPI_BYTE
                                    ,0, md->g_comm2
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
                if (md->rank == 0)
                {
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

                    START_TIMER (ADIOS_TIMER_GLOBALMD);
                    adios_mpi_amr_striping_unit_write(
                                      md->mfh,
                                      -1,
                                      global_index_buffer,
                                      global_index_buffer_offset
                                      );
                    STOP_TIMER (ADIOS_TIMER_GLOBALMD);

                    if (global_index_buffer)
                    {
                        free (global_index_buffer);
                        global_index_buffer = 0;
                        global_index_buffer_size = 0;
                        global_index_buffer_offset = 0;
                    }
                }
            }

            if (is_aggregator (md->rank))
            {
                if (md->g_threading)
                {
                    pthread_join (md->g_swt, NULL);
                }

                FREE (aggr_buff);
            }

            FREE (buffer);
            buffer_size = 0;
            buffer_offset = 0;


            md->g_num_aggregators = 0;
            md->g_have_mdf = 1;
            md->is_color_set = 0;
            md->g_color1 = 0;
            md->g_color2 = 0;

            FREE (md->subfile_name);
            FREE (md->g_is_aggregator);
            FREE (md->g_offsets);
            FREE (md->open_thread_data);
            FREE (md->reopen_thread_data);
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
    if (md->g_ost_skipping_list) {
        free (md->g_ost_skipping_list);
        md->g_ost_skipping_list = NULL;
    }

    adios_clear_index_v1 (md->index);
    return;
}

void adios_mpi_amr_ag_close (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            adios_error (err_invalid_file_mode, 
                        "Only \"w\" mode is supported by MPI_AMR Aggregation IO\n");
            break;
        }
        case adios_mode_write:
        case adios_mode_append:
        case adios_mode_update:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start1;
            uint64_t * pg_sizes = 0, * disp = 0;
            void * aggr_buff = 0;
            struct adios_MPI_thread_data_write write_thread_data;
            int i, new_rank, new_group_size, new_rank2, new_group_size2;
            uint64_t total_data_size = 0, total_data_size1 = 0;;

            START_TIMER (ADIOS_TIMER_COMM);
            //MPI_Comm_split (md->group_comm, md->g_color1, md->rank, &new_comm);
            MPI_Comm_rank (md->g_comm1, &new_rank);
            MPI_Comm_size (md->g_comm1, &new_group_size);

            //MPI_Comm_split (md->group_comm, md->g_color2, md->rank, &new_comm2);
            MPI_Comm_rank (md->g_comm2, &new_rank2);
            MPI_Comm_size (md->g_comm2, &new_group_size2);
            STOP_TIMER (ADIOS_TIMER_COMM);


            // if not merge PG's on the aggregator side
            if (!md->g_merging_pgs)
            {
                uint64_t pg_size;

                pg_size = fd->bytes_written;
                if (pg_size > INT32_MAX)
                {
                    log_warn ("Each processor writes out more than %d bytes, Not supported in aggregation mode.\n", 
                               INT32_MAX);
                }
                pg_sizes = (uint64_t *) malloc (new_group_size * 8);
                disp = (uint64_t *) malloc (new_group_size * 8);
                if (pg_sizes == 0 || disp == 0)
                {
                    adios_error (err_no_memory, 
                            "MPI_AMR method (AG): Cannot allocate buffers (%d bytes) "
                            "for merging process blocks.\n",
                            2*4*new_group_size
                            );
                    return;
                }

                START_TIMER (ADIOS_TIMER_COMM);
                MPI_Allgather (&pg_size, 1, MPI_UNSIGNED_LONG_LONG
                              ,pg_sizes, 1, MPI_UNSIGNED_LONG_LONG
                              ,md->g_comm1);
                STOP_TIMER (ADIOS_TIMER_COMM);

                disp[0] = 0;
                //if (md->rank==0) fprintf (stderr, "rank %d: pg_size[0]=%llu ", md->rank, pg_sizes[0]); 
                for (i = 1; i < new_group_size; i++)
                {
                    disp[i] = disp[i - 1] + pg_sizes[i - 1];
                    //if (md->rank==0) fprintf (stderr, "pg_size[%d]=%llu ", i, pg_sizes[i]); 
                }
                total_data_size = disp[new_group_size - 1]
                                + pg_sizes[new_group_size - 1];
                //if (md->rank==0) fprintf (stderr, "total=%llu\n", total_data_size); 

                if (is_aggregator (md->rank))
                {
                    aggr_buff = malloc (total_data_size);
                    if (aggr_buff == 0)
                    {
                        adios_error (err_no_memory, 
                                "MPI_AMR method (AG): Cannot allocate %lu bytes "
                                "for aggregation buffer.\n"
                                "Need to increase the number of aggregators.\n",
                                total_data_size);
                        return;
                    }
                }
                else
                {
                }

                START_TIMER (ADIOS_TIMER_COMM);
                // This needs to be changed in the future to support > 2 GB data.
                int * int_pg_sizes = (int*) malloc (new_group_size * sizeof(int));
                int * int_disp = (int*) malloc (new_group_size * sizeof(int));
                int int_total_data_size = (int) total_data_size;
                for (i = 0; i < new_group_size; i++)
                {
                    int_pg_sizes[i] = pg_sizes[i];
                    int_disp[i] = disp[i];
                }
                if (total_data_size != (uint64_t) int_total_data_size)
                {
                    adios_error (err_unspecified, 
                            "MPI_AGGRGATE with aggregate_type=1 does not handle >2GB buffers. "
                            "If each process writes less than 2GB, then increase the number of aggregators "
                            "so that the total data size on each aggregator is still < 2GB. "
                            "If some process has more than 2GB data, use aggregate_type=2 "
                            "(the default aggregation method)."
                            );
                }
                MPI_Gatherv (fd->buffer, (int)pg_size, MPI_BYTE
                            ,aggr_buff, int_pg_sizes, int_disp, MPI_BYTE
                            ,0, md->g_comm1);
                STOP_TIMER (ADIOS_TIMER_COMM);
            }
            else 
            {
                // Merge PG's on the aggregator side
                log_warn ("MPI_AMR method (AG): Merging process blocks is not supported yet\n");
            }

            /* Before building index, determine the actual offset where this PG starts */
            fd->current_pg->pg_start_in_file = md->b.pg_index_offset; // aggregator's starting offset (!0 on append)
            if (!md->g_merging_pgs)
            {
                for (i = 0; i < new_rank; i++)
                {
                    fd->current_pg->pg_start_in_file += pg_sizes[i];
                }

                FREE (pg_sizes);
                FREE (disp);
            }

            // build index appending to any existing index
            //fprintf (stderr,"rank %d: Build index with PG offset %llu\n", md->rank, fd->current_pg->pg_start_in_file);
            adios_build_index_v1 (fd, md->index);

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

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, md->g_comm1
                               );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    //if (md->rank==0) fprintf (stderr, "rank %d: ", md->rank);
                    for (i = 0; i < new_group_size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                        //if (md->rank==0) fprintf (stderr, "index_sizes[%d]=%d ", i, index_sizes[i]); 
                    } 
                    //if (md->rank==0) fprintf (stderr, " total index size=%u\n", total_size);

                    /*DEBUG if (md->rank==0) {
                        fprintf (stderr, "rank %d: ", md->rank);
                        for (i = 0; i < new_group_size; i++)
                        {
                            fprintf (stderr, "index_offsets[%d]=%d ", i, index_offsets[i]); 
                        } 
                        fprintf (stderr, " total index size=%u: ", total_size);
                    }*/

                    recv_buffer = malloc (total_size);

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, md->g_comm1
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);

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
                        if (md->g_merging_pgs)
                            new_pg_root = 0;

                        adios_merge_index_v1 (md->index, new_pg_root, 
                                              new_vars_root, new_attrs_root, 0);
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
                                         ,0, md->index);

                    //fprintf (stderr, "rank %d: buffer size = %llu buffer offset=%llu\n", md->rank, buffer_size, buffer_offset);
                    START_TIMER (ADIOS_TIMER_COMM);
                    uint32_t index_size = (uint32_t) buffer_offset;
                    MPI_Gather (&index_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->g_comm1
                               );
                    MPI_Gatherv (buffer, index_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->g_comm1
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);
                }
            }

            // write out indexes in each subfile
            if (is_aggregator (md->rank))
            {
                uint32_t flag = 0;
#if 0
                pthread_join (t, NULL);
                FREE (aggr_buff);

                MPI_File_get_position (md->fh, (MPI_Offset *)&index_start);
#endif
                uint64_t index_start; // = end of PGs, where the index starts
                index_start  = md->b.pg_index_offset; // = end of PGs of previous timesteps
                index_start += total_data_size; //old index start before append + currently written PGs
                /*DEBUG*/
                /*
                fprintf (stderr, "rank %d: write index start=%llu  pg_index_offset=%llu  total_data_size=%llu\n", 
                        md->rank, index_start, md->b.pg_index_offset, total_data_size);
                struct adios_index_process_group_struct_v1 *pg_root = md->index->pg_root;
                i=0;
                while (pg_root) {
                    log_warn ("rank %d: pg %d offset=%llu\n", md->rank, i, pg_root->offset_in_file);
                    pg_root = pg_root->next;
                    i++;
                }*/

                //fprintf (stderr, "rank %d: Before write_index: buffer size = %llu buffer offset=%llu\n", 
                //                  md->rank, buffer_size, buffer_offset);
                adios_write_index_v1 (&buffer, &buffer_size
                                     ,&buffer_offset, index_start
                                     ,md->index);
                //fprintf (stderr, "rank %d: After write_index: buffer size = %llu buffer offset=%llu\n", 
                //                  md->rank, buffer_size, buffer_offset);
//FIXME
                //adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset, flag);
                adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);

                aggr_buff = realloc (aggr_buff, total_data_size + buffer_offset);
                memcpy (aggr_buff + total_data_size, buffer, buffer_offset); 

                // Waiting for the subfile to open if pthread is enabled
                if (md->g_threading)
                {
                    pthread_join (md->g_sot, NULL);
                }

                index_start1 = md->b.pg_index_offset; // starting point to write data at this moment
                total_data_size1 = total_data_size + buffer_offset;
                //fprintf (stderr,"rank %d: Write index+data with PG offset %llu, sizes %llu + %llu = %llu bytes\n", 
                //        md->rank, index_start1, total_data_size, buffer_offset, total_data_size1);

                write_thread_data.fh = &md->fh;
                write_thread_data.base_offset = &index_start1;
                write_thread_data.aggr_buff = aggr_buff;
                write_thread_data.total_data_size = &total_data_size1;

                // Threading the write so that we can overlap write with index collection.
                if (md->g_threading)
                {
                    pthread_create (&md->g_swt, NULL
                            ,adios_mpi_amr_do_write_thread
                            ,(void *) &write_thread_data
                            );
                }
                else
                {
                    START_TIMER (ADIOS_TIMER_IO);
                    adios_mpi_amr_do_write_thread ((void *) &write_thread_data);
                    STOP_TIMER (ADIOS_TIMER_IO);
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

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, md->g_comm2
                               );
                    STOP_TIMER (ADIOS_TIMER_COMM);

                    for (i = 0; i < new_group_size2; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, md->g_comm2
                                );
                    STOP_TIMER (ADIOS_TIMER_COMM);

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

                        // global index would become unsorted on main aggregator during merging 
                        // so sort timesteps if appending
                        adios_merge_index_v1 (md->index, new_pg_root, 
                                              new_vars_root, new_attrs_root, 
                                              (fd->mode == adios_mode_append));
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
                                         ,0, md->index);

                    START_TIMER (ADIOS_TIMER_COMM);
                    MPI_Gather (&buffer_size2, 1, MPI_INT
                               ,0, 0, MPI_INT
                               ,0, md->g_comm2
                               );
                    MPI_Gatherv (buffer2, buffer_size2, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->g_comm2
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
            if (md->rank == 0)
            {
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
                START_TIMER (ADIOS_TIMER_GLOBALMD);
                adios_mpi_amr_striping_unit_write(
                                  md->mfh,
                                  -1,
                                  global_index_buffer,
                                  global_index_buffer_offset
                                  );
                STOP_TIMER (ADIOS_TIMER_GLOBALMD);

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
                if (md->g_threading)
                {
                    pthread_join (md->g_swt, NULL);
                }

                FREE (aggr_buff);
            }
            FREE (buffer);
            buffer_size = 0;
            buffer_offset = 0;

            md->g_num_aggregators = 0;
            md->g_have_mdf = 1;
            md->is_color_set = 0;
            md->g_color1 = 0;
            md->g_color2 = 0;

            FREE (md->subfile_name);
            FREE (md->g_is_aggregator);
            FREE (md->g_ost_skipping_list);
            FREE (md->g_offsets);
            FREE (md->open_thread_data);
            FREE (md->reopen_thread_data);
            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode, 
                    "MPI_AMR method (AG): Unknown file mode (%d) at close time\n", 
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
    if (md->g_ost_skipping_list) {
        free (md->g_ost_skipping_list);
        md->g_ost_skipping_list = NULL;
    }

    adios_clear_index_v1 (md->index);
}

void adios_mpi_amr_buffer_overflow (struct adios_file_struct * fd, 
                                    struct adios_method_struct * method)
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    adios_error (err_buffer_overflow, 
            "rank %d: MPI_AGGREGATE method only works with complete buffering of data between adios_open() "
            "and adios_close(). Portions of global arrays, that do not fit into the "
            "buffer on some processors will not be written by this method to %s\n", 
            md->rank, fd->name);
}

void adios_mpi_amr_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    START_TIMER (ADIOS_TIMER_AD_CLOSE);
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    if (md->g_io_type == ADIOS_MPI_AMR_IO_AG)
    {
        adios_mpi_amr_ag_close (fd, method);
    }
    else if (md->g_io_type == ADIOS_MPI_AMR_IO_BG)
    {
        adios_mpi_amr_bg_close (fd, method);
    }
    else
    {
        adios_error (err_invalid_write_method, "MPI_AMR method: unknown I/O type (%d). "
                "Only MPI_AMR_AGGREGATION and MPI_AMR_BRIGADE are supported\n",
                md->g_io_type);
        return;
    }
    STOP_TIMER (ADIOS_TIMER_AD_CLOSE);

#if defined ADIOS_TIMERS || defined ADIOS_TIMER_EVENTS

    //Finished timing this cycle, swap the timing buffers
    adios_timing_destroy(fd->group->prev_timing_obj);
    fd->group->prev_timing_obj = fd->group->timing_obj;
    fd->group->timing_obj = 0;

    // prev_timing_obj points to unwritten timing info, timing_obj is
    // ready to allocate at the next open

#endif

}

void adios_mpi_amr_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    adios_free_index_v1 (md->index);

#ifdef HAVE_FGR
    fgr_finalize ();
#endif
    if (adios_mpi_amr_initialized)
        adios_mpi_amr_initialized = 0;
}

void adios_mpi_amr_end_iteration (struct adios_method_struct * method)
{
}

void adios_mpi_amr_start_calculation (struct adios_method_struct * method)
{
}

void adios_mpi_amr_stop_calculation (struct adios_method_struct * method)
{
}
