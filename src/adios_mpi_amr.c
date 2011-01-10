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
// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "buffer.h"

enum ADIOS_MPI_AMR_IO_TYPE
{
    ADIOS_MPI_AMR_IO_NONE = 0,
    ADIOS_MPI_AMR_IO_AG   = 1, // simple all to one aggregation
    ADIOS_MPI_AMR_IO_BG   = 2, // Brigade aggregation
};

static int adios_mpi_amr_initialized = 0;
static int * g_is_aggregator = 0;
static int g_num_aggregators = 0;
static int g_merging_pgs = 0;
static int g_num_ost = 0;
static int g_threading = 0;
static int g_color1 = 0;
static int g_color2 = 0;
static MPI_Offset * g_offsets= 0;
static int * g_ost_skipping_list = 0;
static pthread_t g_sot, g_mot, g_swt; // subfile open thread, metadata file open thread, subfile write thread
static enum ADIOS_MPI_AMR_IO_TYPE g_io_type = ADIOS_MPI_AMR_IO_AG;

static struct adios_MPI_thread_data_open open_thread_data1;
static struct adios_MPI_thread_data_open open_thread_data2;

#define is_aggregator(rank)  g_is_aggregator[rank]
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
#define DEFAULT_STRIPE_COUNT 4

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

    void * comm; // temporary until moved from should_buffer to open

    struct adios_bp_buffer_struct_v1 b;

    struct adios_index_process_group_struct_v1 * old_pg_root;
    struct adios_index_var_struct_v1 * old_vars_root;
    struct adios_index_attribute_struct_v1 * old_attrs_root;

    uint64_t vars_start;
    uint64_t vars_header_size;

    uint64_t striping_unit;  // file system stripe size
    uint64_t block_unit;
};

struct adios_MPI_thread_data_open
{
    MPI_File * fh;
    char * name;
    uint64_t * striping_unit;
    char * parameters;
};

struct adios_MPI_thread_data_write
{
    MPI_File * fh;
    uint64_t * base_offset;
    void * aggr_buff;
    int * total_data_size;
    uint64_t * block_unit;
};

#if COLLECT_METRICS
// see adios_adaptive_finalize for what each represents
struct timeval t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25;

// Subtract the `struct timeval' values X and Y,
// storing the result in RESULT.
// Return 1 if the difference is negative, otherwise 0.
static int timeval_subtract (struct timeval * result
                            ,struct timeval * x, struct timeval * y
                            )
{
  // Perform the carry for the later subtraction by updating y.
  if (x->tv_usec < y->tv_usec)
  {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000)
  {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  // Compute the time remaining to wait.
  // tv_usec is certainly positive.
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  // Return 1 if result is negative.
  return x->tv_sec < y->tv_sec;
}

static
void print_metrics (struct adios_MPI_data_struct * md, int iteration)
{
    struct timeval diff;
    if (md->rank == 0)
    {
        timeval_subtract (&diff, &t2, &t1);
        printf ("cc\t%2d\tFile create (stripe setup):\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
        
        timeval_subtract (&diff, &t6, &t5);
        printf ("dd\t%2d\tMass file open:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
        
        timeval_subtract (&diff, &t17, &t16);
        printf ("ee\t%2d\tBuild file offsets:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
    }
    if (md->rank == md->size - 1)
    {   
        timeval_subtract (&diff, &t10, &t9);
        printf ("ff\t%2d\tGlobal index creation:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
        
        timeval_subtract (&diff, &t8, &t7);
        printf ("gg\t%2d\tAll writes complete (w/ local index):\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
        
        timeval_subtract (&diff, &t11, &t0);
        printf ("hh\t%2d\tTotal time:\t%02d.%06d\n"
               ,iteration, diff.tv_sec, diff.tv_usec);
    }
    
    timeval_subtract (&diff, &t13, &t12);
    printf ("ii\t%2d\tLocal index creation:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
    
    timeval_subtract (&diff, &t22, &t21);
    printf ("kk\t%2d\tshould buffer time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
    
    timeval_subtract (&diff, &t19, &t23);
    printf ("ll\t%2d\tclose startup time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
    
    timeval_subtract (&diff, &t19, &t0);
    printf ("mm\t%2d\tsetup time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
    
    timeval_subtract (&diff, &t14, &t20);
    printf ("nn\t%2d\tcleanup time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
    
    timeval_subtract (&diff, &t21, &t0);
    printf ("oo\t%2d\topen->should_buffer time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
    
    timeval_subtract (&diff, &t24, &t21);
    printf ("pp\t%2d\tshould_buffer->write1 time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
    
    timeval_subtract (&diff, &t25, &t24);
    printf ("qq1\t%2d\twrite1->write2 time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
    
    timeval_subtract (&diff, &t23, &t25);
    printf ("qq2\t%2d\twrite2->close start time:\t%6d\t%02d.%06d\n"
           ,iteration, md->rank, diff.tv_sec, diff.tv_usec);
}
#endif

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

int * allocOSTList (int n_ost)
{
    int * ost_list = (int *) malloc (n_ost * sizeof (int));

    if (ost_list == 0)
    {   
        fprintf (stderr, "can not malloc");
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
        fprintf (stderr, "Pointer ost_list is null.\n");
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

static void
adios_mpi_amr_set_striping_unit(MPI_File fh, char *filename, char *parameters)
{
    struct statfs fsbuf;
    int err = 0, flag;
    uint64_t striping_unit = 0;
    uint64_t block_unit = 0;
    uint16_t striping_count = 0;
    char     value[64], *temp_string, *p_count,*p_size;
    int fd, old_mask, perm, n_ost_skipping, n_ost, n, i;
    MPI_Info info_used;

//#ifndef ADIOS_LUSTRE
//    return 0;  // disable stripe-size I/O for non-Lustre file system
//#else
    temp_string = (char *) malloc (strlen (parameters) + 1);
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_count = strstr (temp_string, "stripe_count"))
    {
        char * p = strchr (p_count, '=');
        char * q = strtok (p, ";");
        if (!q)
            striping_count = atoi (q + 1);
        else
            striping_count = atoi (p + 1);
    }

    if (striping_count <= 0)
        striping_count = DEFAULT_STRIPE_COUNT;

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "stripe_size"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            striping_unit = atoi(q + 1);
        else
            striping_unit = atoi(p + 1);
    }

    if (striping_unit <= 0)
        striping_unit = 1048576;

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "block_size"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            block_unit = atoi(q + 1);
        else
            block_unit = atoi(p + 1);
    }
    else
    {
        block_unit = 0;
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
        for (i = 0; i < g_num_ost; i++)
        {
            if (g_ost_skipping_list[i] == 1)
            {
                n_ost_skipping++;
            }
        }

        // the actual # of ost that can be used
        n_ost = g_num_ost - n_ost_skipping;
        if (n_ost <= 0)
        {
            fprintf (stderr, "No OST to use.\n");
            return;
        }

        i = 0;
        while (i < g_num_ost)
        {
            if (g_ost_skipping_list[i] == 0)
            {
                n++;
                if (n - 1 == g_color1 % n_ost)
                    break;
            }
            
            i++;
        }

        lum.lmm_stripe_offset = i;
        ioctl (fd, LL_IOC_LOV_SETSTRIPE
              ,(void *) &lum
              );

        if (err == 0 && lum.lmm_stripe_size > 0) {
            striping_unit = lum.lmm_stripe_size;
        }
        close(fd);
    }
    else
        printf("Warning: open failed on file %s %s.\n",filename,strerror(errno));
//#endif
}

static void
adios_mpi_amr_set_block_unit(uint64_t *block_unit, char *parameters)
{
    char *temp_string, *p_count,*p_size;

    temp_string = (char *) malloc (strlen (parameters) + 1);
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "block_size"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            *block_unit = atoi(q + 1);
        else
            *block_unit = atoi(p + 1);
    }
    else
    {
        *block_unit = 0;
    }

    free (temp_string);
}

static void
adios_mpi_amr_set_aggregation_parameters(char * parameters, int nproc, int rank)
{
    int err = 0, flag, i, aggr_group_size, remain, index;
    char value[64], *temp_string, *p_count,*p_size;

    temp_string = (char *) malloc (strlen (parameters) + 1);
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    // set up # of aggregators
    if (p_size = strstr (temp_string, "num_aggregators"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            g_num_aggregators = atoi(q + 1);
        else
            g_num_aggregators = atoi(p + 1);
    }

    // set up whether to merge PG
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "merging_pgs"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            g_merging_pgs = atoi(q + 1);
        else
            g_merging_pgs = atoi(p + 1);
    }

    // set up whether to thread IO ops
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "threading"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            g_threading = atoi(q + 1);
        else
            g_threading = atoi(p + 1);
    }

    // set up the number of OST to use
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "num_ost"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");
        if (!q)
            g_num_ost = atoi(q + 1);
        else
            g_num_ost = atoi(p + 1);
    }

    // number of ost's is, by default, 672 (jaguar configuration).
    if (g_num_ost <= 0)
        g_num_ost = DEFAULT_NUM_OST;

    // set up which ost's to skip
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    g_ost_skipping_list = allocOSTList (g_num_ost);

    if (p_size = strstr (temp_string, "osts_to_skip"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");

        if (!q)
            g_ost_skipping_list = parseOSTSkipping (g_ost_skipping_list, q + 1, g_num_ost);
        else
            g_ost_skipping_list = parseOSTSkipping (g_ost_skipping_list, p + 1, g_num_ost);
    }

    // set up which ost's to skip
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "aggregation_type"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ";");

        if (!q)
            g_io_type = atoi (q + 1);
        else
            g_io_type = atoi (p + 1);
    }
    else
    {
        g_io_type = ADIOS_MPI_AMR_IO_AG;
    }

    free (temp_string);

    if (g_num_aggregators > nproc || g_num_aggregators <= 0)
    {
        g_num_aggregators = nproc;  //no aggregation
    }

    g_is_aggregator = (int *) malloc (nproc * sizeof(int));
    if (g_is_aggregator == 0)
    {
        fprintf (stderr, "can not malloc\n");
        return;
    }
    memset (g_is_aggregator, 0, nproc * sizeof(int));

    aggr_group_size = nproc / g_num_aggregators;
    remain = nproc - (int) aggr_group_size * g_num_aggregators;

    index = 0;
    for (i = 0; i < g_num_aggregators; i++)
    {
        g_is_aggregator[index] = 1;

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
        g_color1 = rank / aggr_group_size;
        g_color2 = rank % aggr_group_size;
    }
    else
    {
        if (rank < (aggr_group_size + 1) * remain)
        {
            g_color1 = rank / (aggr_group_size + 1);
            g_color2 = rank % (aggr_group_size + 1);
        }
        else
        {
            g_color1 = remain + (rank - (aggr_group_size + 1) * remain) / aggr_group_size;
            g_color2 = (rank - (aggr_group_size + 1) * remain)% aggr_group_size;
        }
    }
}

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
            fprintf (stderr, "Cannot allocate memory in adios_mpi_amr_buffer_write.  "
                             "Requested: %llu\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}

static int
adios_mpi_amr_get_striping_unit(MPI_File fh, char *filename)
{
    struct statfs fsbuf;
    int err, flag;
    uint64_t striping_unit = 1048576;
    char     value[64];
    MPI_Info info_used;

#if COLLECT_METRICS
    gettimeofday (&t1, NULL);
#endif

    // get striping_unit from MPI hint if it has
    MPI_File_get_info(fh, &info_used);
    MPI_Info_get(info_used, "striping_unit", 63, value, &flag);
    MPI_Info_free(&info_used);

    if (flag) return atoi(value);

    // if striping_unit is not set in MPI file info, get it from system
    err = statfs(filename, &fsbuf);
    if (err == -1) {
        printf("Warning: statfs failed %s %s.\n",filename,strerror(errno));
        return striping_unit;
    }

    if (!err && fsbuf.f_type == LUSTRE_SUPER_MAGIC) {
        int fd, old_mask, perm;

        old_mask = umask(022);
        umask(old_mask);
        perm = old_mask ^ 0666;

        fd =  open(filename, O_RDONLY, perm);
        if (fd != -1) {
            struct lov_user_md lum;
            lum.lmm_magic = LOV_USER_MAGIC;
            err = ioctl(fd, LL_IOC_LOV_GETSTRIPE, (void *) &lum);
            if (err == 0 && lum.lmm_stripe_size > 0) {
                striping_unit = lum.lmm_stripe_size;
            }
            close(fd);
        }
        else
            printf("Warning: open failed on file %s %s.\n",filename,strerror(errno));
    }

#if COLLECT_METRICS         
    gettimeofday (&t2, NULL);
#endif
    // set the file striping size
    return striping_unit;
}

static uint64_t
adios_mpi_amr_striping_unit_write(MPI_File    fh
                                  ,MPI_Offset  offset
                                  ,void       *buf
                                  ,uint64_t   len
                                  ,uint64_t   block_unit
                                  )
{
    uint64_t err = -1;
    MPI_Status status;

    if (len == 0)
        return 0;

    if (offset == -1) // use current position
        MPI_File_get_position(fh, &offset);
    else
        MPI_File_seek (fh, offset, MPI_SEEK_SET);

    if (block_unit > 0)
    {
        MPI_Offset  rem_off = offset;
        uint64_t    rem_size = len;
        char       *buf_ptr = buf;

        err = 0;
        while (rem_size > 0)
        {
            uint64_t rem_unit  = block_unit - rem_off % block_unit;
            int write_len = (rem_unit < rem_size) ? rem_unit : rem_size;
            int ret_len;

            MPI_File_write (fh, buf_ptr, write_len, MPI_BYTE, &status);
            MPI_Get_count(&status, MPI_BYTE, &ret_len);
            if (ret_len < 0) {err = ret_len; break;}
            err += ret_len;
            if (ret_len != write_len) break;
            buf_ptr  += write_len;
            rem_off  += write_len;
            rem_size -= write_len;
        }
    }
    else
    {
        uint64_t total_written = 0;
        uint64_t to_write = len;
        int write_len = 0;
        int count;
        char * buf_ptr = buf;
        while (total_written < len)
        {
            write_len = (to_write > INT32_MAX) ? INT32_MAX : to_write;
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
    }

    return err;
}

struct adios_var_struct * adios_mpi_amr_copy_var (struct adios_var_struct * v)
{
    struct adios_var_struct * v_new = (struct adios_var_struct *) 
                            malloc (sizeof (struct adios_var_struct));
    if (v_new == 0)
    {
        fprintf (stderr, "can not malloc\n");
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

void adios_mpi_amr_add_offset (uint64_t var_offset_to_add
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

void adios_mpi_amr_subtract_offset (uint64_t var_offset_to_subtract
                                    ,uint64_t attr_offset_to_subtract
                                    ,struct adios_index_var_struct_v1 * vars_root
                                    ,struct adios_index_attribute_struct_v1 * attrs_root
                                    )
{
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


int adios_mpi_amr_calc_aggregator_index (int rank)
{
    int j = rank - 1;

    if (is_aggregator (rank))
        return rank;

    while (j > -1)
    {
        if (is_aggregator (j))
            break;
        j--;
    }

    return j;
}

void * adios_mpi_amr_do_mkdir (void * param)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) param;
    // 4 bytes for ".dir" 
    char * dir_name = malloc (strlen (fd->name) + 4 + 1);
    sprintf (dir_name, "%s%s", fd->name, ".dir");
    
    mkdir (dir_name, S_IRWXU | S_IRWXG);
  
    free (dir_name);

    return NULL;
}

void * adios_mpi_amr_do_open_thread (void * param)
{
    struct adios_MPI_thread_data_open * td = (struct adios_MPI_thread_data_open *) param;

    unlink (td->name);
    if (td->parameters)
    {
        adios_mpi_amr_set_striping_unit (*td->fh
                                         ,td->name
                                         ,td->parameters
                                         );
    }

    MPI_File_open (MPI_COMM_SELF, td->name
                  ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                  ,MPI_INFO_NULL
                  ,td->fh
                  );

    *td->striping_unit = adios_mpi_amr_get_striping_unit (*td->fh, td->name);
 
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
                              ,*(td->block_unit));

    if (count != *(td->total_data_size))
    {
        fprintf (stderr, "Err in adios_mpi_amr_striping_unit_write()\n");
    }

    return NULL;
}

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
                fprintf (stderr, "communicator not provided and none "
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
                    fprintf (stderr, "communicator not provided and none "
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
                    fprintf (stderr, "communicator not provided but one "
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
        fprintf (stderr, "coordination-communication not provided. "
                         "Using MPI_COMM_WORLD instead\n"
                );

        *comm = MPI_COMM_WORLD;
    }
}

void adios_mpi_amr_init (const char * parameters
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
    md->group_comm = MPI_COMM_NULL;
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
    md->vars_start = 0;
    md->vars_header_size = 0;

    adios_buffer_struct_init (&md->b);
}

int adios_mpi_amr_open (struct adios_file_struct * fd
                        ,struct adios_method_struct * method, void * comm
                        )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

    // we have to wait for the group_size (should_buffer) to get the comm
    // before we can do an open for any of the modes
    md->comm = comm;

    return 1;
}

static
void build_offsets (struct adios_bp_buffer_struct_v1 * b
                   ,MPI_Offset * offsets, uint64_t size, char * group_name
                   ,struct adios_index_process_group_struct_v1 * pg_root
                   )
{
    while (pg_root)
    {
        if (!strcasecmp (pg_root->group_name, group_name))
        {
            MPI_Offset size = 0;

            if (pg_root->next)
            {
                size = pg_root->next->offset_in_file - pg_root->offset_in_file;
            }
            else
            {
                size = b->pg_index_offset - pg_root->offset_in_file;
            }

            offsets [pg_root->process_id * 3] = pg_root->offset_in_file;
            offsets [pg_root->process_id * 3 + 1] = size;
            offsets [pg_root->process_id * 3 + 2] = b->version;
        }

        pg_root = pg_root->next;
    }
}

enum ADIOS_FLAG adios_mpi_amr_should_buffer (struct adios_file_struct * fd
                                             ,struct adios_method_struct * method
                                             )
{
    int i;
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;
    char * name, * name_no_path, * ch;
    char * d_name;
    int err;
    int sig;    // used for coordinating the MPI_File_open

    int previous;
    int current;
    int next;
    uint16_t flag;


    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);

    adios_var_to_comm (fd->group->group_comm
                      ,fd->group->adios_host_language_fortran
                      ,md->comm
                      ,&md->group_comm
                      );

    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }

    fd->group->process_id = md->rank;

    if (md->rank == md->size - 1)
        next = -1;
    else
        next = md->rank + 1;
    previous = md->rank - 1;
    current = md->rank;

    fd->base_offset = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL, &md->fh
                                    );
                if (err != MPI_SUCCESS)
                {
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    fprintf (stderr, "MPI open read failed for %s: '%s'\n"
                            ,name, e
                            );
                    free (name);

                    return adios_flag_no;
                }

                MPI_Offset file_size;
                MPI_File_get_size (md->fh, &file_size);
                md->b.file_size = file_size;

                adios_init_buffer_read_version (&md->b);
                MPI_File_seek (md->fh, md->b.file_size - md->b.length
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.length, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_version (&md->b, &md->b.version);

                adios_init_buffer_read_index_offsets (&md->b);
                // already in the buffer
                adios_parse_index_offsets_v1 (&md->b);

                adios_init_buffer_read_process_group_index (&md->b);
                MPI_File_seek (md->fh, md->b.pg_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.pg_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_process_group_index_v1 (&md->b
                                                   ,&md->old_pg_root
                                                   );

#if 1
                adios_init_buffer_read_vars_index (&md->b);
                MPI_File_seek (md->fh, md->b.vars_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

                adios_init_buffer_read_attributes_index (&md->b);
                MPI_File_seek (md->fh, md->b.attrs_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.attrs_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_attributes_index_v1 (&md->b, &md->old_attrs_root);
#endif

                fd->base_offset = md->b.end_of_pgs;
            }

            if (   md->group_comm != MPI_COMM_NULL
                && md->group_comm != MPI_COMM_SELF
               )
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size * 3
                                                  );
                    memset (offsets, 0, sizeof (MPI_Offset) * md->size * 3);

                    // go through the pg index to build the offsets array
                    build_offsets (&md->b, offsets, md->size
                                  ,fd->group->name, md->old_pg_root
                                  );
                    MPI_Scatter (offsets, 3, MPI_LONG_LONG
                                ,offsets, 3, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
                    md->b.read_pg_offset = offsets [0];
                    md->b.read_pg_size = offsets [1];
                    free (offsets);
                }
                else
                {
                    MPI_Offset offset [3];
                    offset [0] = offset [1] = offset [2] = 0;

                    MPI_Scatter (offset, 3, MPI_LONG_LONG
                                ,offset, 3, MPI_LONG_LONG
                                ,0, md->group_comm
                                );

                    md->b.read_pg_offset = offset [0];
                    md->b.read_pg_size = offset [1];
                    md->b.version = offset [2];
                }
            }

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                // note rank 0 is already open
                // don't open it again here

                if (next != -1)
                {
                    MPI_Isend (&sig, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&sig, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&sig, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                        ,name, e
                        );
                free (name);

                return adios_flag_no;
            }

            break;
        }

        case adios_mode_write:
        {
            if (md->rank == 0)
            {
                adios_mpi_amr_do_mkdir (fd);
            }

            MPI_Barrier (md->group_comm);

            fd->base_offset = 0;
            fd->pg_start_in_file = 0;
            adios_mpi_amr_set_aggregation_parameters (method->parameters
                                                      ,md->size
                                                      ,md->rank
                                                      );
            //adios_mpi_amr_set_block_unit (&md->block_unit, method->parameters);

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
            sprintf (name, "%s%s%s%s.%d", fd->name, ".dir/", method->base_path, name_no_path, g_color1);
            md->subfile_name = strdup (name);
            fd->subfile_index = (uint32_t)g_color1;

            free (name_no_path);

            if (is_aggregator(md->rank))
            {
                if (fd->shared_buffer == adios_flag_yes)
                {
                    if (is_aggregator (md->rank))
                    {
                        // open subfiles
                        open_thread_data1.fh = &md->fh;
                        open_thread_data1.name = md->subfile_name;
                        open_thread_data1.striping_unit = &md->striping_unit;
                        open_thread_data1.parameters = method->parameters;

                        if (g_threading)
                        {
                            pthread_create (&g_sot, NULL
                                           ,adios_mpi_amr_do_open_thread
                                           ,(void *) &open_thread_data1
                                           );
                        }
                        else
                        {
                            adios_mpi_amr_do_open_thread ((void *) &open_thread_data1);
                        }

                        // open metadata file
                        if (md->rank == 0)
                        {
                            open_thread_data2.fh = &md->mfh;
                            open_thread_data2.name = fd->name;
                            open_thread_data2.striping_unit = &md->striping_unit;
                            open_thread_data2.parameters = method->parameters;

                            if (g_threading)
                            {
                                pthread_create (&g_mot, NULL
                                               ,adios_mpi_amr_do_open_thread
                                               ,(void *) &open_thread_data2
                                               );
                            }
                            else
                            {
                                adios_mpi_amr_do_open_thread ((void *) &open_thread_data2);
                            }
                        }
                    }
                }

                if (fd->shared_buffer == adios_flag_no)
                {
                    unlink (name);
                    if (md->rank == 0)
                    {
                        unlink (fd->name);
                    }

                    if (method->parameters)
                    {
                        adios_mpi_amr_set_striping_unit (md->fh
                                                         ,name
                                                         ,method->parameters
                                                         );
                    }

                    err = MPI_File_open (MPI_COMM_SELF, name
                                        ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                        ,MPI_INFO_NULL
                                        ,&md->fh
                                        );
                    md->striping_unit = adios_mpi_amr_get_striping_unit (md->fh, name);

                    if (err != MPI_SUCCESS)
                    {
                        char e [MPI_MAX_ERROR_STRING];
                        int len = 0;
                        memset (e, 0, MPI_MAX_ERROR_STRING);
                        MPI_Error_string (err, e, &len);
                        fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                                ,name, e
                                );
                        free (name);

                        return adios_flag_no;
                    }
                }
            }

            if (md->group_comm != MPI_COMM_NULL)
            {
                fd->base_offset = 0;
                fd->pg_start_in_file = fd->base_offset;
            }
            else
            {
                md->b.pg_index_offset = fd->write_size_bytes;
            }

            break;
        }

        case adios_mode_append:
        {
            int old_file = 1;
            adios_buffer_struct_clear (&md->b);

            err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                ,MPI_INFO_NULL, &md->fh
                                );

            if (err != MPI_SUCCESS)
            {
                old_file = 0;
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,MPI_INFO_NULL, &md->fh
                                    );

                if (err != MPI_SUCCESS)
                {
                    char e [MPI_MAX_ERROR_STRING];
                    int len = 0;
                    memset (e, 0, MPI_MAX_ERROR_STRING);
                    MPI_Error_string (err, e, &len);
                    fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                            ,name, e
                            );
                    free (name);

                    return adios_flag_no;
                }
                md->striping_unit = adios_mpi_amr_get_striping_unit(md->fh, name);
            }

            if (old_file)
            {
                if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
                {
                    if (err != MPI_SUCCESS)
                    {
                        md->b.file_size = 0;
                    }
                    else
                    {
                        MPI_Offset file_size;
                        MPI_File_get_size (md->fh, &file_size);
                        md->b.file_size = file_size;
                    }

                    adios_init_buffer_read_version (&md->b);
                    MPI_File_seek (md->fh, md->b.file_size - md->b.length
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.length, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_version (&md->b, &md->b.version);

                    adios_init_buffer_read_index_offsets (&md->b);
                    // already in the buffer
                    adios_parse_index_offsets_v1 (&md->b);

                    adios_init_buffer_read_process_group_index (&md->b);
                    MPI_File_seek (md->fh, md->b.pg_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.pg_size, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_process_group_index_v1 (&md->b
                                                       ,&md->old_pg_root
                                                       );

                    adios_init_buffer_read_vars_index (&md->b);
                    MPI_File_seek (md->fh, md->b.vars_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

                    adios_init_buffer_read_attributes_index (&md->b);
                    MPI_File_seek (md->fh, md->b.attrs_index_offset
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.attrs_size
                                  ,MPI_BYTE, &md->status
                                  );
                    adios_parse_attributes_index_v1 (&md->b
                                                    ,&md->old_attrs_root
                                                    );

                    fd->base_offset = md->b.end_of_pgs;
                    fd->pg_start_in_file = fd->base_offset;
                }
                else
                {
                    fd->base_offset = 0;
                    fd->pg_start_in_file = 0;
                }

                MPI_File_close (&md->fh);
            }
            else
            {
                fd->base_offset = 0;
                fd->pg_start_in_file = 0;
            }

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                // we know it exists, because we created it if it didn't
                // when reading the old file so can just open wronly
                // but adding the create for consistency with write mode
                // so it is easier to merge write/append later
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
                md->striping_unit = adios_mpi_amr_get_striping_unit(md->fh, name);
                if (next != -1)
                {
                    MPI_Isend (&sig, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&sig, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&sig, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
                md->striping_unit = adios_mpi_amr_get_striping_unit(md->fh, name);
            }

            if (err != MPI_SUCCESS)
            {
                char e [MPI_MAX_ERROR_STRING];
                int len = 0;
                memset (e, 0, MPI_MAX_ERROR_STRING);
                MPI_Error_string (err, e, &len);
                fprintf (stderr, "MPI open write failed for %s: '%s'\n"
                        ,name, e
                        );
                free (name);

                return adios_flag_no;
            }

            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size
                                                  );

                    if (fd->write_size_bytes % md->striping_unit)
                        offsets [0] =  (fd->write_size_bytes / md->striping_unit + 1)
                                     * md->striping_unit;
                    else
                        offsets [0] = fd->write_size_bytes;

                    MPI_Gather (offsets, 1, MPI_LONG_LONG
                               ,offsets, 1, MPI_LONG_LONG
                               ,0, md->group_comm
                               );

                    uint64_t last_offset = offsets [0];
                    offsets [0] = fd->base_offset;
                    for (i = 1; i < md->size; i++)
                    {
                        uint64_t this_offset = offsets [i];
                        offsets [i] = offsets [i - 1] + last_offset;
                        last_offset = this_offset;
                    }
                    md->b.pg_index_offset =   offsets [md->size - 1]
                                            + last_offset;
                    MPI_Scatter (offsets, 1, MPI_LONG_LONG
                                ,offsets, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
                    fd->base_offset = offsets [0];
                    fd->pg_start_in_file = fd->base_offset;
                    free (offsets);
                }
                else
                {
                    MPI_Offset offset;
                    if (fd->write_size_bytes % md->striping_unit)
                        offset =  (fd->write_size_bytes / md->striping_unit + 1)
                                     * md->striping_unit;
                    else
                        offset = fd->write_size_bytes;


                    MPI_Gather (&offset, 1, MPI_LONG_LONG
                               ,&offset, 1, MPI_LONG_LONG
                               ,0, md->group_comm
                               );

                    MPI_Scatter (&offset, 1, MPI_LONG_LONG
                                ,&offset, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
                    fd->base_offset = offset;
                    fd->pg_start_in_file = fd->base_offset;
                }
            }
            else
            {
                md->b.pg_index_offset = fd->write_size_bytes;
            }

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);

            free (name);

            return adios_flag_no;
        }
    }

    free (name);

    if (fd->shared_buffer == adios_flag_no && fd->mode != adios_mode_read)
    {
        uint64_t count;
        // write the process group header
        adios_write_process_group_header_v1 (fd, fd->write_size_bytes);

        if (is_aggregator (md->rank))
        {
            count = adios_mpi_amr_striping_unit_write(
                                  md->fh
                                 ,fd->base_offset
                                 ,fd->buffer
                                 ,fd->bytes_written
                                 ,md->block_unit
                                 );
            if (count != fd->bytes_written)
            {
                fprintf (stderr, "a:MPI method tried to write %llu, "
                                 "only wrote %llu\n"
                        ,fd->bytes_written
                        ,count
                        );
            }
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);

        // setup for writing vars
        adios_write_open_vars_v1 (fd);
        md->vars_start = fd->base_offset;
        md->vars_header_size = fd->offset;
        fd->base_offset += fd->offset;
        MPI_File_seek (md->fh, md->vars_header_size, MPI_SEEK_CUR);
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);
    }

    return fd->shared_buffer;
}

void adios_mpi_amr_write (struct adios_file_struct * fd
                         ,struct adios_var_struct * v
                         ,void * data
                         ,struct adios_method_struct * method
                         )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;
    if (v->got_buffer == adios_flag_yes)
    {
        if (data != v->data)  // if the user didn't give back the same thing
        {
            if (v->free_data == adios_flag_yes)
            {
                free (v->data);
                adios_method_buffer_free (v->data_size);
            }
        }
        else
        {
            // we already saved all of the info, so we're ok.
            return;
        }
    }

    if (fd->shared_buffer == adios_flag_no)
    {
        uint64_t total_size = 0;
        MPI_Comm new_comm;
        int i, new_rank, new_group_size;
        void * aggr_buff = 0;

        // var payload sent for sizing information
        adios_write_var_header_v1 (fd, v);
        adios_write_var_payload_v1 (fd, v);

        MPI_Comm_split (md->group_comm, g_color1, md->rank, &new_comm);
        MPI_Comm_rank (new_comm, &new_rank);
        MPI_Comm_size (new_comm, &new_group_size);

        int bytes_written[new_group_size];
        int disp[new_group_size];

        MPI_Gather (&fd->bytes_written, 1, MPI_INT
                   ,bytes_written, 1, MPI_INT
                   ,0, new_comm);

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
                fprintf (stderr, "Can not alloc aggregation buffer.\n"
                                 "Need to increase the number of aggregators.\n"
                        );
                return;
            }
        }
  
        MPI_Gatherv (fd->buffer, fd->bytes_written, MPI_BYTE
                    ,aggr_buff, bytes_written, disp, MPI_BYTE
                    ,0, new_comm);

        fd->vars_written += new_group_size - 1;

        uint64_t count = 0;
        if (is_aggregator(md->rank))
        {
            count = adios_mpi_amr_striping_unit_write(
                           md->fh
                          ,-1
                          ,aggr_buff
                          ,total_size
                          ,md->block_unit
                          );
            if (count != total_size)
            {
                fprintf (stderr, "b:MPI method tried to write %llu, "
                                 "only wrote %llu\n"
                        ,total_size
                        ,count
                        );
            }

            FREE (aggr_buff);
        }
        else
        {
            // Non-aggregators do nothing
        }

        // Broadcast new offsets to all processors in the communicator.
        uint64_t new_offsets[new_group_size];

        if (is_aggregator (md->rank))
        {
            new_offsets[0] = v->write_offset;
            for (i = 1; i < new_group_size; i++)
            {
                new_offsets[i] = new_offsets[i - 1] + bytes_written[i - 1];
            }
        }

        MPI_Bcast (new_offsets, new_group_size, MPI_LONG_LONG, 0, new_comm);
        v->write_offset = new_offsets[new_rank];

        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);
    }
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

    if (v->data && v->free_data)
    {
        adios_method_buffer_free (v->data_size);
        free (v->data);
    }

    mem_allowed = adios_method_buffer_alloc (*size);
    if (mem_allowed == *size)
    {
        *buffer = malloc (*size);
        if (!*buffer)
        {
            adios_method_buffer_free (mem_allowed);
            fprintf (stderr, "Out of memory allocating %llu bytes for %s\n"
                    ,*size, v->name
                    );
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
        fprintf (stderr, "OVERFLOW: Cannot allocate requested buffer of %llu "
                         "bytes for %s\n"
                ,*size
                ,v->name
                );
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

static void adios_mpi_amr_do_read (struct adios_file_struct * fd
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

    switch (md->b.version & ADIOS_VERSION_NUM_MASK)
    {
        case 1:
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
                    var_payload.payload = v1->data;
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,&var_payload
                                                    ,v1->data_size
                                                    );
                }
                else
                {
                    printf ("MPI read: skipping name: %s path: %s\n"
                           ,var_header.name, var_header.path
                           );
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
            fprintf (stderr, "MPI read: file version unknown: %u\n"
                    ,md->b.version
                    );
            return;
    }

    adios_buffer_struct_clear (&md->b);
}

static
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
}

void adios_mpi_amr_bg_close (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            )
{
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
            fprintf (stderr, "Only \"w\" is supported by MPI_AMR Brigade IO\n");
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
            MPI_Comm new_comm, new_comm2;

            MPI_Comm_split (md->group_comm, g_color1, md->rank, &new_comm);
            MPI_Comm_rank (new_comm, &new_rank);
            MPI_Comm_size (new_comm, &new_group_size);

            MPI_Comm_split (md->group_comm, g_color2, md->rank, &new_comm2);
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
                    count = adios_mpi_amr_striping_unit_write(
                                   md->fh
                                  ,md->vars_start
                                  ,fd->buffer
                                  ,md->vars_header_size
                                  ,md->block_unit
                                  );
                    if (count != md->vars_header_size)
                    {
                        fprintf (stderr, "d:MPI method tried to write %llu, "
                                         "only wrote %d\n"
                                ,md->vars_header_size
                                ,count
                                );
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
                            fprintf (stderr, "Can not alloc aggregation buffer.\n"
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
                        count = adios_mpi_amr_striping_unit_write(
                                          md->fh,
                                          -1,
                                          aggr_buff, //fd->buffer,
                                          total_size, //fd->bytes_written,
                                          md->block_unit);
                        if (count != total_size)
                        {
                            fprintf (stderr, "e:MPI method tried to write %llu, "
                                             "only wrote %llu\n"
                                     ,fd->bytes_written
                                     ,count
                                     );
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

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);

                // fd->vars_start gets updated with the size written
                if (is_aggregator(md->rank))
                {
                    *(uint16_t *)fd->buffer = *(uint16_t *)fd->buffer * new_group_size;
                    count = adios_mpi_amr_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size,
                                  md->block_unit);
                    if (count != md->vars_header_size)
                    {
                        fprintf (stderr, "f:MPI method tried to write %llu, "
                                         "only wrote %llu\n"
                                ,md->vars_header_size
                                ,count
                                );
                    }
                }

                fd->offset = 0;
                fd->bytes_written = 0;

                MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
            }

            // if not merge PG's on the aggregator side
            if (fd->shared_buffer == adios_flag_yes && !g_merging_pgs)
            {
                //printf ("do not merge pg\n");
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
                    fprintf (stderr, "can not malloc\n");
                    return;
                }

                MPI_Allgather (&pg_size, 1, MPI_INT
                              ,pg_sizes, 1, MPI_INT
                              ,new_comm);

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
                        fprintf (stderr, "Warning: The max allowed aggregation buffer is %llu bytes.\n"
                                "But this ADIOS method needs extra %llu bytes for aggregation\n"
                                ,MAX_AGG_BUF, 2 * max_data_size);
                        return;
                    }

                    aggr_buff = malloc (max_data_size);
                    recv_buff = malloc (max_data_size);
                    if (aggr_buff == 0 || recv_buff == 0)
                    {
                        fprintf (stderr, "Can not alloc %d bytes for aggregation buffer.\n"
                                ,max_data_size);
                        return;
                    }
                }
                else
                {
                    if (max_data_size > MAX_AGG_BUF)
                    {
                        fprintf (stderr, "The max allowed aggregation buffer is %llu.\n"
                                ,MAX_AGG_BUF);
                        return;
                    }

                    recv_buff = malloc (max_data_size);
                    if (recv_buff == 0)
                    {
                        fprintf (stderr, "Can not alloc %d bytes for aggregation buffer.\n"
                                ,max_data_size);
                        return;
                    }
                }

                total_data_size = disp[new_group_size - 1]
                                + pg_sizes[new_group_size - 1];

                if (is_aggregator (md->rank))
                {
                    if (g_threading)
                    {
                        pthread_join (g_sot, NULL);
                    }

                    index_start1 = 0;
                    for (i = 0; i < new_group_size; i++)
                    {
                        if (i + 1 < new_group_size)
                        {
                            MPI_Irecv (recv_buff, pg_sizes[i + 1], MPI_BYTE, new_rank + 1
                                      ,0, new_comm, &request);
                        }

                        write_thread_data.fh = &md->fh;
                        write_thread_data.base_offset = &index_start1;
                        write_thread_data.aggr_buff = (i == 0) ? fd->buffer : aggr_buff;
                        write_thread_data.total_data_size = &pg_sizes[i];
                        write_thread_data.block_unit = &md->block_unit;

                        // This write call is not threaded
                        adios_mpi_amr_do_write_thread ((void *) &write_thread_data);

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
                                 ,0, new_comm);
                    }
                    else
                    {
                        for (i = new_rank + 1; i < new_group_size; i++)
                        {
                            // Recv data from upstream rank
                            MPI_Irecv (recv_buff, pg_sizes[i], MPI_BYTE, new_rank + 1
                                      ,0, new_comm, &request);

                            if (i == new_rank + 1)
                                // Send my data to downstream rank
                                MPI_Send (fd->buffer, pg_size, MPI_BYTE, new_rank - 1
                                         ,0, new_comm);

                            MPI_Wait (&request, &status);
                            // Send it to downstream rank
                            MPI_Send (recv_buff, pg_sizes[i], MPI_BYTE, new_rank - 1
                                     ,0, new_comm);
                        }
                    }
                }

                FREE (aggr_buff);
                FREE (recv_buff);
            }

            // build index appending to any existing index
            adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root
                                 ,&md->old_attrs_root
                                 );

            if (fd->shared_buffer == adios_flag_yes && !g_merging_pgs)
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

                    adios_mpi_amr_add_offset (var_offset_to_add
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

                        adios_parse_process_group_index_v1 (&md->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );
                        if (g_merging_pgs)
                            new_pg_root = 0;

                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
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
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
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
                index_start = total_data_size;

                adios_write_index_v1 (&buffer, &buffer_size
                                     ,&buffer_offset, index_start
                                     ,md->old_pg_root
                                     ,md->old_vars_root
                                     ,md->old_attrs_root
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
                    write_thread_data.block_unit = &md->block_unit;

                    if (g_threading)
                    {
                        pthread_create (&g_swt, NULL
                                       ,adios_mpi_amr_do_write_thread
                                       ,(void *) &write_thread_data
                                       );
                    }
                    else
                    {
                        adios_mpi_amr_do_write_thread ((void *) &write_thread_data); 
                    }
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

                        adios_parse_process_group_index_v1 (&md->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
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
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
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
                                     ,md->old_pg_root, md->old_vars_root, md->old_attrs_root
                                     );

                flag |= ADIOS_VERSION_HAVE_SUBFILE;

                adios_write_version_flag_v1 (&global_index_buffer
                                            ,&global_index_buffer_size
                                            ,&global_index_buffer_offset
                                            ,flag
                                            );

                if (g_threading)
                {
                    pthread_join (g_mot, NULL);
                }

                adios_mpi_amr_striping_unit_write(
                                  md->mfh,
                                  -1,
                                  global_index_buffer,
                                  global_index_buffer_offset,
                                  md->block_unit);

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
                if (g_threading)
                {
                    pthread_join (g_swt, NULL);
                }

                FREE (aggr_buff);
            }

            FREE (buffer);
            buffer_size = 0;
            buffer_offset = 0;

            adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);
            adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                                 ,md->old_attrs_root
                                 );
            new_pg_root = 0;
            new_vars_root = 0;
            new_attrs_root = 0;
            md->old_pg_root = 0;
            md->old_vars_root = 0;
            md->old_attrs_root = 0;

            g_num_aggregators = 0;
            g_color1 = 0;
            g_color2 = 0;

            FREE (md->subfile_name);
            FREE (g_is_aggregator);
            FREE (g_offsets);
        }

        break;
    }

    return;
}

void adios_mpi_amr_ag_close (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            // read the index to find the place to start reading
            adios_mpi_amr_do_read (fd, method);
            struct adios_var_struct * v = fd->group->vars;
            while (v)
            {
                v->data = 0;
                v = v->next;
            }

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
            //pthread_t t, t1;
            //struct adios_MPI_thread_data_open open_thread_data1;
            //struct adios_MPI_thread_data_open open_thread_data2;
            struct adios_MPI_thread_data_write write_thread_data;
            int i, new_rank, new_group_size, new_rank2, new_group_size2, total_data_size = 0, total_data_size1 = 0;;
            MPI_Comm new_comm, new_comm2;

            MPI_Comm_split (md->group_comm, g_color1, md->rank, &new_comm);
            MPI_Comm_rank (new_comm, &new_rank);
            MPI_Comm_size (new_comm, &new_group_size);

            MPI_Comm_split (md->group_comm, g_color2, md->rank, &new_comm2);
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
                    count = adios_mpi_amr_striping_unit_write(
                                   md->fh
                                  ,md->vars_start
                                  ,fd->buffer
                                  ,md->vars_header_size
                                  ,md->block_unit
                                  );
                    if (count != md->vars_header_size)
                    {
                        fprintf (stderr, "d:MPI method tried to write %llu, "
                                         "only wrote %d\n"
                                ,md->vars_header_size
                                ,count
                                );
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
                            fprintf (stderr, "Can not alloc aggregation buffer.\n"
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
                        count = adios_mpi_amr_striping_unit_write(
                                          md->fh,
                                          -1,
                                          aggr_buff, //fd->buffer,
                                          total_size, //fd->bytes_written,
                                          md->block_unit);
                        if (count != total_size)
                        {
                            fprintf (stderr, "e:MPI method tried to write %llu, "
                                             "only wrote %llu\n"
                                     ,fd->bytes_written
                                     ,count
                                     );
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

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);

                // fd->vars_start gets updated with the size written
                if (is_aggregator(md->rank))
                {
                    *(uint16_t *)fd->buffer = *(uint16_t *)fd->buffer * new_group_size;
                    count = adios_mpi_amr_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size,
                                  md->block_unit);
                    if (count != md->vars_header_size)
                    {
                        fprintf (stderr, "f:MPI method tried to write %llu, "
                                         "only wrote %llu\n"
                                ,md->vars_header_size
                                ,count
                                );
                    }
                }

                fd->offset = 0;
                fd->bytes_written = 0;

                MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
            }

            // if not merge PG's on the aggregator side
            if (fd->shared_buffer == adios_flag_yes && !g_merging_pgs)
            {
                //printf ("do not merge pg\n");
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
                    fprintf (stderr, "can not malloc\n");
                    return;
                }
//moved to adios_group_size
#if 0
                if (is_aggregator (md->rank))
                {
                    // open subfiles
                    open_thread_data1.fh = &md->fh;
                    open_thread_data1.name = md->subfile_name;
                    open_thread_data1.striping_unit = &md->striping_unit;
                    open_thread_data1.parameters = method->parameters;

                    pthread_create (&t, NULL
                                   ,adios_mpi_amr_do_open_thread
                                   ,(void *) &open_thread_data1
                                   );

                    // open metadata file
                    if (md->rank == 0)
                    {
                        open_thread_data2.fh = &md->mfh;
                        open_thread_data2.name = fd->name;
                        open_thread_data2.striping_unit = &md->striping_unit;
                        open_thread_data2.parameters = method->parameters;

                        pthread_create (&t1, NULL
                                       ,adios_mpi_amr_do_open_thread
                                       ,(void *) &open_thread_data2
                                       );
                    }
                }
#endif
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
                        fprintf (stderr, "The max allowed aggregation buffer is %llu.\n"
                                         "Need to increase the number of aggregators.\n"
                                ,MAX_AGG_BUF);
                        return;
                    }
                    aggr_buff = malloc (total_data_size);
                    if (aggr_buff == 0)
                    {
                        fprintf (stderr, "Can not alloc %d bytes for aggregation buffer.\n"
                                         "Need to increase the number of aggregators.\n"
                                ,total_data_size);
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

            // Merge PG's on the aggregator side
            if (fd->shared_buffer == adios_flag_yes && g_merging_pgs)
            {
                printf ("do merge pg\n");
                // Merge PG's on the aggregator side
                struct adios_bp_buffer_struct_v1 b;
                struct adios_process_group_header_struct_v1 pg_header;
                struct adios_vars_header_struct_v1 vars_header;
                int pg_size, header_size;
                uint32_t attr_size;
                uint64_t vars_count_offset;

                // Unfortunately, we need to walk through the buffer to get vars count
                b.buff = fd->buffer;
                b.change_endianness = md->b.change_endianness;
                b.offset = 0;
                b.length = fd->bytes_written;

                adios_parse_process_group_header_v1 (&b, &pg_header);
                vars_count_offset = b.offset;
                adios_clear_process_group_header_v1 (&pg_header);

                adios_parse_vars_header_v1 (&b, &vars_header);
                header_size = b.offset;
                attr_size = adios_mpi_amr_calculate_attributes_size (fd);

                // reset vars count
                if (vars_header.count * new_group_size > UINT16_MAX)
                {
                    fprintf (stderr, "Vars count exceed UINT16_MAX");
                    return;
                }
                *(uint16_t *) (b.buff + vars_count_offset) = 
                                vars_header.count * new_group_size;

                // attributes size is save in the end
                adios_mpi_amr_buffer_write (&fd->buffer, &fd->buffer_size, &fd->offset, &attr_size, SHIM_FOOTER_SIZE);
                fd->bytes_written += SHIM_FOOTER_SIZE;
            
                // PG header, vars header, vars, attrs header, attrs, 4 bytes
                if (is_aggregator(md->rank))
                {
                    pg_size = fd->bytes_written;
                }
                else
                {
                    // Non-aggregator process doesn't need to send pg header + vars header
                    pg_size = fd->bytes_written - header_size;
                }

                pg_sizes = (int *) malloc (new_group_size * 4);
                attr_sizes = (int *) malloc (new_group_size * 4);
                disp = (int *) malloc (new_group_size * 4);
                sendbuf = (int *) malloc (2 * 4);
                recvbuf = (int *) malloc (new_group_size * 2 * 4);
                if (pg_sizes == 0 || attr_sizes == 0 || disp == 0
                 || sendbuf == 0 || recvbuf == 0)
                {
                    fprintf (stderr, "can not malloc\n");
                    return;
                }
  
                sendbuf[0] = pg_size;
                sendbuf[1] = attr_size + SHIM_FOOTER_SIZE;

//moved to adios_group_size
#if 0
                if (is_aggregator (md->rank))
                {
                    open_thread_data1.fh = &md->fh;
                    open_thread_data1.name = md->subfile_name;
                    open_thread_data1.striping_unit = &md->striping_unit;
                    open_thread_data1.parameters = method->parameters;

                    pthread_create (&t, NULL
                                   ,adios_mpi_amr_do_open_thread
                                   ,(void *) &open_thread_data1
                                   );

                    if (md->rank == 0)
                    {
                        open_thread_data2.fh = &md->mfh;
                        open_thread_data2.name = fd->name;
                        open_thread_data2.striping_unit = &md->striping_unit;
                        open_thread_data2.parameters = method->parameters;

                        pthread_create (&t1, NULL
                                       ,adios_mpi_amr_do_open_thread
                                       ,(void *) &open_thread_data2
                                       );
                    }
                }
#endif
                MPI_Allgather (sendbuf, 2, MPI_INT
                              ,recvbuf, 2, MPI_INT
                              ,new_comm);

                for (i = 0; i < new_group_size; i++)
                {
                    pg_sizes[i] = recvbuf[i * 2];
                    attr_sizes[i] = recvbuf[i * 2 + 1];
                }
   
                free (sendbuf);
                free (recvbuf);

                disp[0] = 0;
                for (i = 1; i < new_group_size; i++)
                {
                    disp[i] = disp[i - 1] + pg_sizes[i - 1];
                }
                total_data_size = disp[new_group_size - 1]
                                + pg_sizes[new_group_size - 1];

                if (is_aggregator (md->rank))
                {
                    aggr_buff = malloc (total_data_size);
                    if (aggr_buff == 0)
                    {
                        fprintf (stderr, "Can not alloc %d bytes for aggregation buffer.\n"
                                         "Need to increase the number of aggregators.\n"
                                ,total_data_size);
                        return;
                    }
                }
                else
                {
                }

                if (is_aggregator (md->rank))
                { 
                    uint32_t aggr_attr_size = 0;
                    void * aggr_attr_buff, * temp_aggr_buff, * temp_aggr_attr_buff;
                    uint16_t new_attr_count = 0, new_attr_len = 0;

                    MPI_Gatherv (fd->buffer, pg_size, MPI_BYTE
                                ,aggr_buff, pg_sizes, disp, MPI_BYTE
                                ,0, new_comm);

                    for (i= 0; i < new_group_size; i++)
                    {
                        aggr_attr_size += *(uint32_t *)(aggr_buff
                                                       + disp[i]
                                                       + pg_sizes[i]
                                                       - SHIM_FOOTER_SIZE
                                                       );
                    }

                    aggr_attr_buff = malloc (aggr_attr_size);
                    if (aggr_attr_buff == 0)
                    {
                        fprintf (stderr, "can not alloc %d bytes for aggregation buffer.\n"
                                         "Need to increase the number of aggregators.\n"
                                ,aggr_attr_size);
                        return;
                    }

                    temp_aggr_buff = aggr_buff
                                   + pg_sizes[0]
                                   - SHIM_FOOTER_SIZE
                                   - attr_size;
                    temp_aggr_attr_buff = aggr_attr_buff;

                    for (i= 0; i < new_group_size; i++)
                    {
                        uint32_t temp_attr_size = *(uint32_t *)(aggr_buff
                                                               + disp[i]
                                                               + pg_sizes[i]
                                                               - SHIM_FOOTER_SIZE
                                                               );
                        void * temp_attr_ptr = aggr_buff + disp[i] + pg_sizes[i] 
                                             - SHIM_FOOTER_SIZE - temp_attr_size;

                        new_attr_count += *(uint16_t *)temp_attr_ptr;

                        if (i == 0)
                        {
                            new_attr_len += *(uint64_t *)(temp_attr_ptr + ATTR_COUNT_SIZE);

                            memcpy (temp_aggr_attr_buff, temp_attr_ptr, temp_attr_size);
                            temp_aggr_attr_buff += temp_attr_size;
                        }
                        else
                        {
                            new_attr_len += *(uint64_t *)(temp_attr_ptr + ATTR_COUNT_SIZE)
                                          - ATTR_COUNT_SIZE - ATTR_LEN_SIZE;

                            memcpy (temp_aggr_attr_buff
                                   ,temp_attr_ptr + ATTR_COUNT_SIZE + ATTR_LEN_SIZE
                                   ,temp_attr_size - ATTR_COUNT_SIZE - ATTR_LEN_SIZE
                                   );
                            temp_aggr_attr_buff += temp_attr_size - ATTR_COUNT_SIZE - ATTR_LEN_SIZE;

                            memmove (temp_aggr_buff
                                    ,aggr_buff + disp[i]
                                    ,pg_sizes[i] - temp_attr_size - SHIM_FOOTER_SIZE
                                    );
                            temp_aggr_buff += pg_sizes[i] - temp_attr_size - SHIM_FOOTER_SIZE;
                        }
                    }

                    memcpy (temp_aggr_buff, aggr_attr_buff, aggr_attr_size);

                    *(uint16_t *)temp_aggr_buff = new_attr_count;  //attrs count
                    *(uint64_t *)(temp_aggr_buff + ATTR_COUNT_SIZE) = new_attr_len;  //attrs length

                    free (aggr_attr_buff);
                }
                else
                {
                    MPI_Gatherv (fd->buffer + header_size, pg_size, MPI_BYTE
                                ,aggr_buff, pg_sizes, disp, MPI_BYTE
                                ,0, new_comm);
                }

                uint64_t count = 0;
                if (is_aggregator(md->rank))
                {
#if 0
                    pthread_join (t, NULL);
          
                    write_thread_data.fh = &md->fh;
                    write_thread_data.base_offset = &fd->base_offset;
                    write_thread_data.aggr_buff = aggr_buff;
                    write_thread_data.total_data_size = &total_data_size;
                    write_thread_data.block_unit = &md->block_unit;
 
                    pthread_create (&t, NULL
                                   ,adios_mpi_amr_do_write_thread
                                   ,(void *) &write_thread_data
                                   );
#endif

#if 0 
                    count = adios_mpi_amr_striping_unit_write(
                               md->fh
                              ,fd->base_offset
                              ,aggr_buff
                              ,total_data_size
                              ,md->block_unit);

                    if (count != total_data_size)
                    {
                        fprintf (stderr, "Err in adios_mpi_amr_striping_unit_write()\n");
                        return;
                    }

                    free (aggr_buff);
#endif
                }
                else
                {
                }
            }

            // build index appending to any existing index
            adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root
                                 ,&md->old_attrs_root
                                 );

            if (fd->shared_buffer == adios_flag_yes && !g_merging_pgs)
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

                    adios_mpi_amr_add_offset (var_offset_to_add
                                              ,attr_offset_to_add
                                              ,md->old_vars_root
                                              ,md->old_attrs_root
                                              );
                }

                // pg_sizes, disp are no longer needed from this point on.
                free (pg_sizes);
                free (disp);
            }

            if (fd->shared_buffer == adios_flag_yes && g_merging_pgs)
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

                    adios_mpi_amr_subtract_offset (var_base_offset
                                                   ,attr_base_offset
                                                   ,md->old_vars_root
                                                   ,md->old_attrs_root
                                                   );

                    for (i = 0; i < new_group_size; i++)
                    {
                        attr_offset_to_add += pg_sizes[i] - attr_sizes[i];  
                    }

                    for (i = 0; i < new_rank; i++)
                    {
                        attr_offset_to_add += attr_sizes[i] - SHIM_FOOTER_SIZE;  
                        var_offset_to_add += pg_sizes[i] - attr_sizes[i]; 
                    }

                    adios_mpi_amr_add_offset (var_offset_to_add
                                              ,attr_offset_to_add
                                              ,md->old_vars_root
                                              ,md->old_attrs_root
                                              );
                }

                // pg_sizes, attr_sizs, disp are no longer needed from this point on.
                free (pg_sizes);
                free (attr_sizes);
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

                        adios_parse_process_group_index_v1 (&md->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );
                        if (g_merging_pgs)
                            new_pg_root = 0;

                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
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
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
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
                                     ,md->old_pg_root
                                     ,md->old_vars_root
                                     ,md->old_attrs_root
                                     );
//FIXME
                //adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset, flag);
                adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);

                if (fd->shared_buffer == adios_flag_yes)
                {
                    aggr_buff = realloc (aggr_buff, total_data_size + buffer_offset);
                    memcpy (aggr_buff + total_data_size, buffer, buffer_offset); 

                    // Waiting for the subfile to open if pthread is enabled
                    if (g_threading)
                    {
                        pthread_join (g_sot, NULL);
                    }

                    index_start1 = 0;
                    total_data_size1 = total_data_size + buffer_offset;

                    write_thread_data.fh = &md->fh;
                    write_thread_data.base_offset = &index_start1;
                    write_thread_data.aggr_buff = aggr_buff;
                    write_thread_data.total_data_size = &total_data_size1;
                    write_thread_data.block_unit = &md->block_unit;

                    // Threading the write so that we can overlap write with index collection.
                    if (g_threading)
                    {
                        pthread_create (&g_swt, NULL
                                       ,adios_mpi_amr_do_write_thread
                                       ,(void *) &write_thread_data
                                       );
                    }
                    else
                    {
                        adios_mpi_amr_do_write_thread ((void *) &write_thread_data);
                    }
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

                        adios_parse_process_group_index_v1 (&md->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );

                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
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
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
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
                                     ,md->old_pg_root, md->old_vars_root, md->old_attrs_root
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
                // Waiting for metadata file to open
                if (g_threading)
                {
                    pthread_join (g_mot, NULL);
                }

                adios_mpi_amr_striping_unit_write(
                                  md->mfh,
                                  -1,
                                  global_index_buffer,
                                  global_index_buffer_offset,
                                  md->block_unit);

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
                if (g_threading)
                {
                    pthread_join (g_swt, NULL);
                }

                FREE (aggr_buff);
            }
            FREE (buffer);
            buffer_size = 0;
            buffer_offset = 0;

            adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);
            adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                                 ,md->old_attrs_root
                                 );
            new_pg_root = 0;
            new_vars_root = 0;
            new_attrs_root = 0;
            md->old_pg_root = 0;
            md->old_vars_root = 0;
            md->old_attrs_root = 0;

            g_num_aggregators = 0;
            g_color1 = 0;
            g_color2 = 0;

            FREE (md->subfile_name);
            FREE (g_is_aggregator);
            FREE (g_ost_skipping_list);
            FREE (g_offsets);

            break;
        }

        case adios_mode_append:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = md->b.pg_index_offset;

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
                count = adios_mpi_amr_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size,
                                  md->block_unit);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "d:MPI method tried to write %llu, "
                                     "only wrote %llu\n"
                            ,md->vars_header_size
                            ,count
                            );
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

                while (a)
                {
                    adios_write_attribute_v1 (fd, a);
                    count = adios_mpi_amr_striping_unit_write(
                                  md->fh,
                                  -1,
                                  fd->buffer,
                                  fd->bytes_written,
                                  md->block_unit);
                    if (count != fd->bytes_written)
                    {
                        fprintf (stderr, "e:MPI method tried to write %llu, "
                                         "only wrote %llu\n"
                                ,fd->bytes_written
                                ,count
                                );
                    }
                    fd->base_offset += count;
                    fd->offset = 0;
                    fd->bytes_written = 0;
                    adios_shared_buffer_free (&md->b);

                    a = a->next;
                }

                // set it up so that it will start at 0, but have correct sizes
                fd->offset = fd->base_offset - md->vars_start;
                fd->vars_start = 0;
                fd->buffer_size = 0;
                adios_write_close_attributes_v1 (fd);
                // fd->vars_start gets updated with the size written
                count = adios_mpi_amr_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size,
                                  md->block_unit);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "f:MPI method tried to write %llu, "
                                     "only wrote %llu\n"
                            ,md->vars_header_size
                            ,count
                            );
                }
                fd->offset = 0;
                fd->bytes_written = 0;
            }

            // build index appending to any existing index
            adios_build_index_v1 (fd, &md->old_pg_root, &md->old_vars_root
                                 ,&md->old_attrs_root
                                 );
            // if collective, gather the indexes from the rest and call
            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    int * index_sizes = malloc (4 * md->size);
                    int * index_offsets = malloc (4 * md->size);
                    char * recv_buffer = 0;
                    uint32_t size = 0;
                    uint32_t total_size = 0;
                    int i;

                    MPI_Gather (&size, 1, MPI_INT
                               ,index_sizes, 1, MPI_INT
                               ,0, md->group_comm
                               );

                    for (i = 0; i < md->size; i++)
                    {
                        index_offsets [i] = total_size;
                        total_size += index_sizes [i];
                    }

                    recv_buffer = malloc (total_size);

                    MPI_Gatherv (&size, 0, MPI_BYTE
                                ,recv_buffer, index_sizes, index_offsets
                                ,MPI_BYTE, 0, md->group_comm
                                );

                    char * buffer_save = md->b.buff;
                    uint64_t buffer_size_save = md->b.length;
                    uint64_t offset_save = md->b.offset;

                    for (i = 1; i < md->size; i++)
                    {
                        md->b.buff = recv_buffer + index_offsets [i];
                        md->b.length = index_sizes [i];
                        md->b.offset = 0;

                        adios_parse_process_group_index_v1 (&md->b
                                                           ,&new_pg_root
                                                           );
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );
                        adios_merge_index_v1 (&md->old_pg_root
                                             ,&md->old_vars_root
                                             ,&md->old_attrs_root
                                             ,new_pg_root, new_vars_root
                                             ,new_attrs_root
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
                                         ,0, md->old_pg_root
                                         ,md->old_vars_root
                                         ,md->old_attrs_root
                                         );

                    MPI_Gather (&buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->group_comm
                               );
                    MPI_Gatherv (buffer, buffer_size, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->group_comm
                                );
                }
            }

            if (fd->shared_buffer == adios_flag_yes)
            {
                // everyone writes their data
                adios_mpi_amr_striping_unit_write(
                                  md->fh,
                                  fd->base_offset,
                                  fd->buffer,
                                  fd->bytes_written,
                                  md->block_unit);
            }

            if (md->rank == 0)
            {
                uint32_t flag = 0;

                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,index_start, md->old_pg_root
                                     ,md->old_vars_root
                                     ,md->old_attrs_root
                                     );

                flag |= ADIOS_VERSION_HAVE_SUBFILE;
                adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);
/*
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
*/

                adios_mpi_amr_striping_unit_write(
                                  md->fh,
                                  md->b.pg_index_offset,
                                  buffer,
                                  buffer_offset,
                                  md->block_unit);
            }

            free (buffer);

            adios_clear_index_v1 (new_pg_root, new_vars_root, new_attrs_root);
            adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                                 ,md->old_attrs_root
                                 );
            new_pg_root = 0;
            new_vars_root = 0;
            new_attrs_root = 0;
            md->old_pg_root = 0;
            md->old_vars_root = 0;
            md->old_attrs_root = 0;

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);
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
    md->group_comm = MPI_COMM_NULL;

    adios_clear_index_v1 (md->old_pg_root, md->old_vars_root
                         ,md->old_attrs_root
                         );
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
#if COLLECT_METRICS
    print_metrics (md, iteration++);
#endif
}

void adios_mpi_amr_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    if (g_io_type == ADIOS_MPI_AMR_IO_AG)
    {
        adios_mpi_amr_ag_close (fd, method);
    }
    else if (g_io_type == ADIOS_MPI_AMR_IO_BG)
    {
        adios_mpi_amr_bg_close (fd, method);
    }
    else
    {
        fprintf (stderr, "unknown I/O type. Only MPI_AMR_AGGREGATION and MPI_AMR_BRIGADE are supported\n");
        return;
    }

    //g_io_type = ADIOS_MPI_AMR_IO_NONE;
}

void adios_mpi_amr_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
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
