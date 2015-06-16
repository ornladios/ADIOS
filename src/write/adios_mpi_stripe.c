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

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/buffer.h"
#include "core/util.h"

static int adios_mpi_stripe_initialized = 0;

#define COLLECT_METRICS 0

struct adios_MPI_data_struct
{
    MPI_File fh;
    MPI_Request req;
    MPI_Status status;
    MPI_Comm group_comm;
    int rank;
    int size;

    struct adios_bp_buffer_struct_v1 b;

    struct adios_index_process_group_struct_v1 * old_pg_root;
    struct adios_index_var_struct_v1 * old_vars_root;
    struct adios_index_attribute_struct_v1 * old_attrs_root;

    uint64_t vars_start;
    uint64_t vars_header_size;

    uint64_t striping_unit;  // file system stripe size
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
#    include <sys/param.h>
#    include <sys/mount.h>
#else
#    include <sys/statfs.h>
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

static int 
adios_mpi_stripe_set_striping_unit(MPI_File fh, char *filename, char *parameters)
{
    struct statfs fsbuf;
    int err = 0, flag;
    uint64_t striping_unit = 0;
    uint16_t striping_count = 0;
    char     value[64], *temp_string, *p_count,*p_size;
    MPI_Info info_used;

#if COLLECT_METRICS
    gettimeofday (&t1, NULL);
#endif


//#ifndef ADIOS_LUSTRE
//    return 0;  // disable stripe-size I/O for non-Lustre file system
//#else

    temp_string = (char *) malloc (strlen (parameters) + 1);
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_count = strstr (temp_string, "stripe_count"))
    {
        char * p = strchr (p_count, '=');
        char * q = strtok (p, ",");
        if (!q)
            striping_count = atoi (q + 1);
        else
            striping_count = atoi (p + 1);
    }

    if (striping_count <= 0)
        striping_count = 4;

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);

    if (p_size = strstr (temp_string, "stripe_size"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ",");
        if (!q)
            striping_unit = atoi(q + 1);
        else
            striping_unit = atoi(p + 1);
    }

    if (striping_unit <= 0)
        striping_unit = 1048576;

    free (temp_string);

    int fd, old_mask, perm;

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
        lum.lmm_stripe_offset = -1;
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

#if COLLECT_METRICS         
    gettimeofday (&t2, NULL);
#endif
    // set the file striping size
    return striping_unit;
//#endif
}

static int
adios_mpi_stripe_get_striping_unit(MPI_File fh, char *filename)
{
    struct statfs fsbuf;
    int err, flag;
    uint64_t striping_unit = 1048576;
    char     value[64];
    MPI_Info info_used;

#if COLLECT_METRICS
    gettimeofday (&t1, NULL);
#endif

#ifndef ADIOS_LUSTRE
    return 0;  // disable stripe-size I/O for non-Lustre file system
#else
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
#endif
}

static uint64_t
adios_mpi_stripe_striping_unit_write(MPI_File    fh,
                              MPI_Offset  offset,
                              void       *buf,
                              uint64_t   len,
                              uint64_t   striping_unit)
{
    uint64_t err = -1;
    MPI_Status status;

    if (len == 0) return 0;

    if (offset == -1) // use current position
        MPI_File_get_position(fh, &offset);
    else
        MPI_File_seek (fh, offset, MPI_SEEK_SET);

    if (striping_unit > 0) {
        MPI_Offset  rem_off = offset;
        uint64_t    rem_size = len;
        char       *buf_ptr = buf;

        err = 0;
        while (rem_size > 0) {
            uint64_t rem_unit  = striping_unit - rem_off % striping_unit;
            int write_len = (rem_unit < rem_size) ? rem_unit : rem_size;
            int ret_len;

#ifdef _WKL_CHECK_STRIPE_IO
printf("adios_mpi_stripe_striping_unit_write offset=%12lld len=%12d\n",offset,write_len);offset+=write_len;
#endif
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
    else {
#ifdef _WKL_CHECK_STRIPE_IO
printf("adios_mpi_stripe_striping_unit_write offset=%12lld len=%12d\n",offset,len);
#endif
        uint64_t total_written = 0;
        uint64_t to_write = len;
        int write_len = 0;
        int count;
        char * buf_ptr = buf;
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
    }
    return err;
}

void adios_mpi_stripe_init (const PairStruct * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    if (!adios_mpi_stripe_initialized)
    {
        adios_mpi_stripe_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_MPI_data_struct));
    md = (struct adios_MPI_data_struct *) method->method_data;
    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->rank = 0;
    md->size = 0;
    md->group_comm = method->init_comm;//unused, adios_open sets current comm
    md->old_pg_root = 0;
    md->old_vars_root = 0;
    md->old_attrs_root = 0;
    md->vars_start = 0;
    md->vars_header_size = 0;

    adios_buffer_struct_init (&md->b);
}

int adios_mpi_stripe_open (struct adios_file_struct * fd
                   ,struct adios_method_struct * method, MPI_Comm comm
                   )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

#if COLLECT_METRICS
    gettimeofday (&t0, NULL); // only used on rank == size - 1, but we don't
                              // have the comm yet to get the rank/size
#endif
    // we have to wait for the group_size (should_buffer) to get the comm
    // before we can do an open for any of the modes
    md->group_comm = comm;
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }

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

enum ADIOS_FLAG adios_mpi_stripe_should_buffer (struct adios_file_struct * fd
                                        ,struct adios_method_struct * method
                                        )
{
    int i;
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data;
    char * name;
    int err;
    int flag;    // used for coordinating the MPI_File_open

    int previous;
    int current;
    int next;

#if COLLECT_METRICS
    gettimeofday (&t21, NULL);
#endif

    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);

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
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INT, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
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
            fd->base_offset = 0;
            fd->pg_start_in_file = 0;
#if COLLECT_METRICS                     
            gettimeofday (&t16, NULL);
#endif
            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size
                                                  );

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
                    MPI_Offset offset = fd->write_size_bytes;

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

#if COLLECT_METRICS
            gettimeofday (&t17, NULL);
#endif

#if COLLECT_METRICS   
            gettimeofday (&t5, NULL);
#endif
            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                unlink (name);  // make sure clean

                if (method->parameters)
                    md->striping_unit = adios_mpi_stripe_set_striping_unit (md->fh 
                                                                           ,name
                                                                           ,method->parameters);

                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );

                md->striping_unit = adios_mpi_stripe_get_striping_unit(md->fh, name);

                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INT, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
                md->striping_unit = adios_mpi_stripe_get_striping_unit(md->fh, name);
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
#if COLLECT_METRICS
            gettimeofday (&t6, NULL);
#endif

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
                md->striping_unit = adios_mpi_stripe_get_striping_unit(md->fh, name);
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

            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size
                                                  );

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
                    MPI_Offset offset = fd->write_size_bytes;

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
                md->striping_unit = adios_mpi_stripe_get_striping_unit(md->fh, name);
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&flag, 1, MPI_INT, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&flag, 1, MPI_INT, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
                md->striping_unit = adios_mpi_stripe_get_striping_unit(md->fh, name);
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
        // write the process group header
        adios_write_process_group_header_v1 (fd, fd->write_size_bytes);

        uint64_t count;
        count = adios_mpi_stripe_striping_unit_write(
                          md->fh,
                          fd->base_offset,
                          fd->buffer,
                          fd->bytes_written,
                          md->striping_unit);
        if (count != fd->bytes_written)
        {
            fprintf (stderr, "a:MPI method tried to write %llu, "
                             "only wrote %llu\n"
                    ,fd->bytes_written
                    ,count
                    );
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

#if COLLECT_METRICS
    gettimeofday (&t22, NULL);
#endif
    return fd->shared_buffer;
}

void adios_mpi_stripe_write (struct adios_file_struct * fd
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
        // var payload sent for sizing information
        adios_write_var_header_v1 (fd, v);

        uint64_t count;
        count = adios_mpi_stripe_striping_unit_write(
                          md->fh,
                          -1,
                          fd->buffer,
                          fd->bytes_written,
                          md->striping_unit);
        if (count != fd->bytes_written)
        {
            fprintf (stderr, "b:MPI method tried to write %llu, "
                             "only wrote %llu\n"
                    ,fd->bytes_written
                    ,count
                    );
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);

        // write payload
        // adios_write_var_payload_v1 (fd, v);
        uint64_t var_size = adios_get_var_size (v, v->data);
        if (fd->base_offset + var_size > fd->pg_start_in_file + fd->write_size_bytes) 
            fprintf (stderr, "adios_mpi_write exceeds pg bound. File is corrupted. "
                             "Need to enlarge group size. \n");
        count = adios_mpi_stripe_striping_unit_write(
                          md->fh,
                          -1,
                          v->data,
                          var_size,
                          md->striping_unit);
        if (count != var_size)
        {
            fprintf (stderr, "c:MPI method tried to write %llu, "
                             "only wrote %llu\n"
                    ,var_size
                    ,count
                    );
        }
        fd->base_offset += count;
        fd->offset = 0;
        fd->bytes_written = 0;
        adios_shared_buffer_free (&md->b);
    }
#if COLLECT_METRICS
    static int writes_seen = 0;

    if (writes_seen == 0) gettimeofday (&t24, NULL);
    else if (writes_seen == 1) gettimeofday (&t25, NULL);
    writes_seen++;
#endif
}

void adios_mpi_stripe_get_write_buffer (struct adios_file_struct * fd
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

void adios_mpi_stripe_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )
{
    v->data = buffer;
    v->data_size = buffer_size;
}

static void adios_mpi_stripe_do_read (struct adios_file_struct * fd
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
            fprintf (stderr, "MPI read: file version unknown: %u\n", version);
            return;
    }

    adios_buffer_struct_clear (&md->b);
}

void adios_mpi_stripe_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;
    struct adios_attribute_struct * a = fd->group->attributes;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;
#if COLLECT_METRICS
    gettimeofday (&t23, NULL);
    static int iteration = 0;
#endif

    switch (fd->mode)
    {
        case adios_mode_read:
        {
            // read the index to find the place to start reading
            adios_mpi_stripe_do_read (fd, method);
            struct adios_var_struct * v = fd->group->vars;
            while (v)
            {
                v->data = v->adata = 0;
                v = v->next;
            }

            break;
        }

        case adios_mode_write:
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
                int retlen;
                count = adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size,
                                  md->striping_unit);
                if (count != md->vars_header_size)
                {
                    fprintf (stderr, "d:MPI method tried to write %llu, "
                                     "only wrote %d\n"
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
                    if (fd->base_offset + fd->bytes_written > fd->pg_start_in_file + fd->write_size_bytes)
                        fprintf (stderr, "adios_mpi_write exceeds pg bound. File is corrupted. "
                                         "Need to enlarge group size. \n");
                    count = adios_mpi_stripe_striping_unit_write(
                                      md->fh,
                                      -1,
                                      fd->buffer,
                                      fd->bytes_written,
                                      md->striping_unit);
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
                count = adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size,
                                  md->striping_unit);
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

#if COLLECT_METRICS
            gettimeofday (&t19, NULL);
#endif
#if COLLECT_METRICS
            gettimeofday (&t7, NULL);
#endif
#if COLLECT_METRICS
            gettimeofday (&t12, NULL);
#endif
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

#if COLLECT_METRICS
            gettimeofday (&t13, NULL);
#endif
            if (fd->shared_buffer == adios_flag_yes)
            {
                // everyone writes their data
                if (fd->base_offset + fd->bytes_written > fd->pg_start_in_file + fd->write_size_bytes)
                    fprintf (stderr, "adios_mpi_write exceeds pg bound. File is corrupted. "
                             "Need to enlarge group size. \n");

                adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  fd->base_offset,
                                  fd->buffer,
                                  fd->bytes_written,
                                  md->striping_unit);
            }

            if (md->rank == 0)
            {
                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,index_start, md->old_pg_root
                                     ,md->old_vars_root
                                     ,md->old_attrs_root
                                     );
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  md->b.pg_index_offset,
                                  buffer,
                                  buffer_offset,
                                  md->striping_unit);
            }
#if COLLECT_METRICS
            gettimeofday (&t8, NULL);
#endif
#if COLLECT_METRICS
            gettimeofday (&t20, NULL);
#endif
#if COLLECT_METRICS
            gettimeofday (&t14, NULL);
#endif

            if (buffer)
            {
                free (buffer);
                buffer = 0;
                buffer_size = 0;
                buffer_offset = 0;
            }

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
#if COLLECT_METRICS
            gettimeofday (&t11, NULL);
            t15.tv_sec = t11.tv_sec;
            t15.tv_usec = t11.tv_usec;
#endif

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
                count = adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size,
                                  md->striping_unit);
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
                    count = adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  -1,
                                  fd->buffer,
                                  fd->bytes_written,
                                  md->striping_unit);
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
                count = adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  md->vars_start,
                                  fd->buffer,
                                  md->vars_header_size,
                                  md->striping_unit);
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
                adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  fd->base_offset,
                                  fd->buffer,
                                  fd->bytes_written,
                                  md->striping_unit);
            }

            if (md->rank == 0)
            {
                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,index_start, md->old_pg_root
                                     ,md->old_vars_root
                                     ,md->old_attrs_root
                                     );
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                adios_mpi_stripe_striping_unit_write(
                                  md->fh,
                                  md->b.pg_index_offset,
                                  buffer,
                                  buffer_offset,
                                  md->striping_unit);
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

    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
    {
        md->group_comm = MPI_COMM_NULL;
    }

    md->fh = 0;
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

void adios_mpi_stripe_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
    if (adios_mpi_stripe_initialized)
        adios_mpi_stripe_initialized = 0;
}

void adios_mpi_stripe_end_iteration (struct adios_method_struct * method)
{
}

void adios_mpi_stripe_start_calculation (struct adios_method_struct * method)
{
}

void adios_mpi_stripe_stop_calculation (struct adios_method_struct * method)
{
}
