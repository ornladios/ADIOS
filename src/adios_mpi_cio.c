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

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

#define STR_LEN 1000

MPI_Info adios_mpi_info = MPI_INFO_NULL;

static int adios_mpi_cio_initialized = 0;

struct adios_MPI_Collective_data_struct
{
    MPI_File fh;
    MPI_Request req;
    MPI_Status status;
    int rank;
    int size;
    void * buffer;
    uint64_t buffer_size;
    uint64_t start;
    int last_var_write_yes; // was the last item asked to write a write="yes"?
    uint64_t base_offset;

    void * comm;
};

static void adios_var_to_comm (const char * varname
                              ,struct adios_var_struct * vars
                              ,MPI_Comm * comm
                              ,enum ADIOS_FLAG unique_var_names
                              )
{
    if (varname)
    {
        struct adios_var_struct * var =
                    adios_find_var_by_name (vars, varname, unique_var_names);

        if (var)
        {
            if (var->data)
            {
                *comm = *(MPI_Comm *) var->data;
            }
            else
            {
                fprintf (stderr, "coordination-communication: %s not provided. "
                                 "Using MPI_COMM_WORLD instead\n"
                        ,varname
                        );
                *comm = MPI_COMM_WORLD;
            }
        }
        else
        {
            fprintf (stderr, "coordination-communicator: %s not found in "
                             "adios-group.  Using MPI_COMM_WORLD instead\n"
                    ,varname
                    );

            *comm = MPI_COMM_WORLD;
        }
    }
    else
    {
        *comm = MPI_COMM_NULL;
    }
}

void adios_mpi_cio_init (const char * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_MPI_Collective_data_struct * md;
    if (!adios_mpi_cio_initialized)
    {
        adios_mpi_cio_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_MPI_Collective_data_struct));
    md = (struct adios_MPI_Collective_data_struct *) method->method_data;
    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->rank = 0;
    md->size = 0;
    md->buffer = 0;
    md->buffer_size = 0;
    md->start = 0;
    md->last_var_write_yes = 0;
    md->base_offset = 0;
}

int adios_mpi_cio_open (struct adios_file_struct * fd
                       ,struct adios_method_struct * method, void * comm
                       )
{
    char name [STR_LEN];
    struct adios_MPI_Collective_data_struct * md = (struct adios_MPI_Collective_data_struct *)
                                                    method->method_data;

    md->comm = comm;
    sprintf (name, "%s%s", method->base_path, fd->name);
    MPI_Offset file_size_for_offset = 0;

    // if it is an append, add the base file size to the initial offset
    if (fd->mode == adios_mode_append)
    {
        MPI_Comm group_comm = MPI_COMM_NULL;
        int rank = -1;

        if (fd->group->group_comm)
        {
            adios_var_to_comm (fd->group->group_comm, fd->group->vars
                              ,&group_comm, fd->group->all_unique_var_names);
            if (group_comm != MPI_COMM_NULL)
                MPI_Comm_rank (group_comm, &rank);
        }

        if (group_comm == MPI_COMM_NULL || rank == 0)
        {
            int err;
            err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                ,adios_mpi_info, &md->fh
                                );
            if (err != MPI_SUCCESS)
            {
                file_size_for_offset = 0;
            }
            else
            {
                MPI_File_get_size (md->fh, &file_size_for_offset);
            }
printf ("file size: %llu\n", file_size_for_offset);
            MPI_File_close (&md->fh);  // just in case it opened
        }

        if (group_comm != MPI_COMM_NULL)
            MPI_Bcast (&file_size_for_offset, 1, MPI_LONG_LONG, 0, group_comm);

        md->base_offset += file_size_for_offset;
printf ("base offset: %llu\n", md->base_offset);
    }

    return 1;
}

enum ADIOS_FLAG adios_mpi_cio_should_buffer (struct adios_file_struct * fd
                                            ,struct adios_method_struct * method
                                            )
{
    return adios_flag_yes;   // as far as we care, buffer
}

void adios_mpi_cio_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,void * data
                     ,struct adios_method_struct * method
                     )
{
#if 0
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

    if (0) //v->copy_on_write == adios_flag_yes)
    {
        uint64_t var_size;
        uint64_t mem_allowed;

        var_size = adios_size_of_var (v, data);
        mem_allowed = adios_method_buffer_alloc (var_size);
        if (mem_allowed == var_size)
        {
            v->free_data = adios_flag_yes;
            v->data = adios_dupe_data_scalar (v->type, data);
            v->data_size = var_size;
        }
        else
        {
            adios_method_buffer_free (mem_allowed);
//TODO: need to handle overflow of allocatable space
            fprintf (stderr, "OVERFLOW in MPI when writing %s\n", v->name);
            v->free_data = adios_flag_no;
            v->data = 0;
            v->data_size = 0;
        }
    }
    else // otherwise, should be safe to hold pointer
    {
        v->data = data;
        v->free_data = adios_flag_no;
        v->data_size = 0;
    }
#endif
}

void adios_mpi_cio_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
#if 0
    uint64_t mem_allowed;

    if (*size == 0)
    {
        *size = adios_size_of_var (v, (void *) "");
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
#endif
}

void adios_mpi_cio_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )
{
    v->data = buffer;
    v->data_size = buffer_size;
}

static void adios_mpi_cio_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
#if 0
    char name [STR_LEN];
    struct adios_MPI_Collective_data_struct * md = (struct adios_MPI_Collective_data_struct *)
                                                      method->method_data; 
    struct adios_var_struct * v = 0;

    uint64_t my_data_len = adios_data_size (fd->group);

    uint64_t my_offset = 0;
    MPI_Status mstatus;
    MPI_Offset read_offset;
    int amode = MPI_MODE_RDONLY;
    int err;
    MPI_Comm group_comm;

    if (fd->group->group_comm)
    {
        adios_var_to_comm (fd->group->group_comm, fd->group->vars, &group_comm
                          ,fd->group->all_unique_var_names);
    }
    else
    {
        group_comm = MPI_COMM_NULL;
    }

    // this code section is important if we want to write from
    // several processes to a single file.  For one file per process
    // writes, this is unnecessary.

    // ******************************
    // Begin Coordinate file offet
    // ******************************
    //char buf1 [STR_LEN];
    //char buf2 [STR_LEN];

    sprintf (name, "%s%s", method->base_path, fd->name);
    // make a filename based on the size and current node
    // sprintf (buf1, "%%s.%%0%dd", (int) (log10 (size)));
    // sprintf (buf2, buf1, fd->name, fd->group->current);
    if (group_comm != MPI_COMM_NULL)
    {
        int do_free_info = 0;

        MPI_Scan(&my_data_len, &my_offset, 1, MPI_LONG_LONG, MPI_SUM, group_comm);
        my_offset -= my_data_len;
        read_offset = md->base_offset + my_offset;

        if (adios_mpi_info == MPI_INFO_NULL) {
            do_free_info = 1;
            MPI_Info_create(&adios_mpi_info);
        }
        MPI_Info_set(adios_mpi_info, "romio_no_indep_rw", "true");
        // this hint will enable aggregator I/O

        MPI_File_open (group_comm, name, MPI_MODE_RDONLY, adios_mpi_info, &md->fh);

        if (do_free_info) {
            MPI_Info_free(&adios_mpi_info);
            adios_mpi_info = MPI_INFO_NULL;
        }
    }
    else
    {
        MPI_File_open (MPI_COMM_SELF, name, amode, adios_mpi_info, &md->fh);
        read_offset = md->base_offset;
    }
    // *****************************
    // End Coordinate file offet
    // *****************************

    if (md->buffer_size < my_data_len)
    {
        uint64_t mem_needed = 0;
        uint64_t mem_allowed = 0;
        mem_needed = my_data_len - md->buffer_size;
        mem_allowed = adios_method_buffer_alloc (mem_needed);
        if (mem_needed != mem_allowed)
        {
            adios_method_buffer_free (mem_allowed);
            fprintf (stderr, "OVERFLOW!!\n");
        }
        
        md->buffer_size = md->buffer_size + mem_allowed;
        free (md->buffer);
        md->buffer = malloc (md->buffer_size);
    }
    
    err = MPI_File_read_at_all(md->fh, read_offset, md->buffer, my_data_len, MPI_BYTE, &mstatus);
    if (err != MPI_SUCCESS)
    {
        fprintf (stderr, "Could not read from file %s\n", name);
    }

    adios_parse_buffer (fd, md->buffer, my_data_len);

    v = fd->group->vars;
    while (v)
    {
        v->data = 0;
        v = v->next;
    }
#endif
}

static void adios_mpi_cio_do_write (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               )
{
#if 0
    struct adios_MPI_Collective_data_struct * md = (struct adios_MPI_Collective_data_struct *)
                                                       method->method_data;
    char name [STR_LEN];
    uint64_t my_data_len = adios_data_size (fd->group);
    uint64_t mem_needed = 0;
    uint64_t mem_allowed = 0;
    uint64_t my_offset = 0;
    MPI_Status mstatus;
    MPI_Offset write_offset = 0;
    int overflow = 0;
    uint64_t end = 0;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_attribute_struct * a = fd->group->attributes;
    int err;
    MPI_Comm group_comm;

    if (fd->group->group_comm)
    {
        adios_var_to_comm (fd->group->group_comm, fd->group->vars, &group_comm
                          ,fd->group->all_unique_var_names);
    }
    else
    {
        group_comm = MPI_COMM_NULL;
    }

    // this code section is important if we want to write from
    // several processes to a single file.  For one file per process
    // writes, this is unnecessary.

    // *****************************
    // Begin Coordinate file offet
    // *****************************
    //char buf1 [STR_LEN];
    //char buf2 [STR_LEN];

    sprintf (name, "%s%s", method->base_path, fd->name);
    // make a filename based on the size and current node
    // sprintf (buf1, "%%s.%%0%dd", (int) (log10 (size)));
    // sprintf (buf2, buf1, fd->name, fd->group->current);
    if (group_comm != MPI_COMM_NULL)
    {
        int do_free_info = 0;

        MPI_Scan(&my_data_len, &my_offset, 1, MPI_LONG_LONG, MPI_SUM, group_comm);
        my_offset -= my_data_len;
        write_offset = md->base_offset + my_offset;

        if (adios_mpi_info == MPI_INFO_NULL) {
            do_free_info = 1;
            MPI_Info_create(&adios_mpi_info);
        }
        MPI_Info_set(adios_mpi_info, "romio_no_indep_rw", "true");
        // this hint will enable aggregator I/O

        MPI_File_open (group_comm, name, MPI_MODE_CREATE|MPI_MODE_WRONLY, adios_mpi_info, &md->fh);

        if (do_free_info) {
            MPI_Info_free(&adios_mpi_info);
            adios_mpi_info = MPI_INFO_NULL;
        }
    }
    else
    {
        MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_WRONLY | MPI_MODE_CREATE
                      ,adios_mpi_info, &md->fh
                      );
        write_offset = md->base_offset;
    }

    // ******************************
    //  End Coordinate file offet
    // ******************************

    if (my_data_len > md->buffer_size)
    {
        mem_needed = my_data_len - md->buffer_size;
        mem_allowed = adios_method_buffer_alloc (mem_needed);
        if (mem_needed != mem_allowed)
        {
            adios_method_buffer_free (mem_allowed);
            fprintf (stderr, "OVERFLOW!!\n");
            overflow = 1;
        }

        md->buffer_size = md->buffer_size + mem_allowed;
        free (md->buffer);
        md->buffer = malloc (md->buffer_size);
    }

    while (v && !overflow)
    {
        overflow = adios_do_write_var (v
                                      ,md->buffer
                                      ,md->buffer_size
                                      ,md->start
                                      ,&end
                                      );

        md->start = end;
        v = v->next;
    }

    while (a && !overflow)
    {
        overflow = adios_do_write_attribute (a
                                            ,md->buffer
                                            ,md->buffer_size
                                            ,md->start
                                            ,&end
                                            );

        md->start = end;
        a = a->next;
    }

    if (overflow)
    {
        overflow = 0;
        int64_t file;
        char name [STR_LEN];

        sprintf (name, "%sOVERFLOW_%s", method->base_path, fd->name);
        bw_set_write (0); // set write to a file directly
        bw_fopen_ (name, &file);
        v = fd->group->vars;
        a = fd->group->attributes;

        while (v && !overflow)
        {
            overflow = adios_do_write_var (v
                                          ,&file
                                          ,9223372036854775807LL // LLONG_MAX
                                          ,0
                                          ,&end
                                          );

            v = v->next;
        }
        while (a && !overflow)
        {
            overflow = adios_do_write_attribute (a
                                                ,md->buffer
                                                ,md->buffer_size
                                                ,0
                                                ,&end
                                                );

            a = a->next;
        }

        bw_fclose_ (&file);
        bw_set_write (1); // set write to the buffer
    }
    else
    {
        MPI_File_write_at_all(md->fh, write_offset, md->buffer, end, MPI_BYTE, &mstatus);
    }

    // clear out the cached data for the caller to be able to continue to work
    v = fd->group->vars;
    while (v)
    {
        if (v->free_data == adios_flag_yes)
        {
            if (v->data)
            {
                free (v->data);
                adios_method_buffer_free (v->data_size);
            }
        }

        v->data = 0;
        v->data_size = 0;
        v->free_data = adios_flag_no;
        v->got_buffer = adios_flag_no;
        v = v->next;
    }
#endif
}

void adios_mpi_cio_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_Collective_data_struct * md = (struct adios_MPI_Collective_data_struct *)
                                                 method->method_data;

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        adios_mpi_cio_do_write (fd, method);
    }

    if (fd->mode == adios_mode_read)
    {
        adios_mpi_cio_do_read (fd, method);
    }

    if (md && md->fh)
        MPI_File_close (&md->fh);

    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->start = 0;
    md->last_var_write_yes = 0;
    md->base_offset = 0;
}

void adios_mpi_cio_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
    if (adios_mpi_cio_initialized)
        adios_mpi_cio_initialized = 0;
}

void adios_mpi_cio_end_iteration (struct adios_method_struct * method)
{
}

void adios_mpi_cio_start_calculation (struct adios_method_struct * method)
{
}

void adios_mpi_cio_stop_calculation (struct adios_method_struct * method)
{
}
