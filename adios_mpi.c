#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// mpi
#include "mpi.h"

// xml parser
#include <mxml.h>

#include "binpack-general.h"
#include "binpack-utils.h"
#include "br-utils.h"
#include "bw-utils.h"
#include "adios.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"

extern MPI_Comm adios_mpi_comm_world;
extern MPI_Comm adios_mpi_comm_self;
extern MPI_Info adios_mpi_info;

static int adios_mpi_initialized = 0;

struct adios_MPI_data_struct
{
    MPI_File fh;
    MPI_Request req;
    MPI_Status status;
    int rank;
    int size;
    void * buffer;
    unsigned long long buffer_size;
    unsigned long long start;
    int last_var_write_yes; // was the last item asked to write a write="yes"?
};

static void adios_var_to_comm (const char * varname
                              ,struct adios_var_struct * vars
                              ,MPI_Comm * comm
                              )
{
    if (varname)
    {
        struct adios_var_struct * var = adios_find_var_by_name (vars, varname);

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
                *comm = adios_mpi_comm_world;
            }
        }
        else
        {
            fprintf (stderr, "coordination-communicator: %s not found in "
                             "adios-group.  Using MPI_COMM_WORLD instead\n"
                    ,varname
                    );

            *comm = adios_mpi_comm_world;
        }
    }
    else
    {
        *comm = MPI_COMM_NULL;
    }
}

void adios_mpi_init (const char * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_MPI_data_struct * md;
    if (!adios_mpi_initialized)
    {
        adios_mpi_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_MPI_data_struct));
    md = (struct adios_MPI_data_struct *) method->method_data;
    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->rank = 0;
    md->size = 0;
    md->buffer = 0;
    md->buffer_size = 0;
    md->start = 0;
    md->last_var_write_yes = 0;
}

void adios_mpi_open (struct adios_file_struct * fd
                    ,struct adios_method_struct * method
                    )
{
    char name [STR_LEN];
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

    MPI_Offset file_size_for_offset = 0;

    // if it is an append, add the base file size to the initial offset
    if (fd->mode == adios_mode_append)
    {
        MPI_Comm group_comm = MPI_COMM_NULL;
        int rank = -1;

        if (fd->group->group_comm)
        {
            adios_var_to_comm (fd->group->group_comm, fd->group->vars, &group_comm);
            if (group_comm != MPI_COMM_NULL)
                MPI_Comm_rank (group_comm, &rank);
        }

        if (group_comm == MPI_COMM_NULL || rank == 0)
        {
            int err;
            err = MPI_File_open (adios_mpi_comm_self, name, MPI_MODE_RDONLY
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
            //printf ("group_comm %d, file size: %llu\n", group_comm,file_size_for_offset);
            MPI_File_close (&md->fh);  // just in case it opened
        }

        if (group_comm != MPI_COMM_NULL)
            MPI_Bcast (&file_size_for_offset, 1, MPI_LONG_LONG, 0, group_comm);

        fd->base_offset += file_size_for_offset;
        //printf ("base offset: %llu\n", fd->base_offset);
    }
}

void adios_mpi_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,void * data
                     ,struct adios_method_struct * method
                     )
{
    if (v->got_buffer)
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

    if (v->copy_on_write == adios_flag_yes)
    {
        unsigned long long var_size;
        unsigned long long mem_allowed;

        var_size = adios_size_of_var (v, data);
        mem_allowed = adios_method_buffer_alloc (var_size);
        if (mem_allowed == var_size)
        {
            v->free_data = adios_flag_yes;
            v->data = adios_dupe_data (v, data);
            v->data_size = var_size;
        }
        else
        {
            adios_method_buffer_free (mem_allowed);
// need to handle overflow of allocatable space
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
}

void adios_mpi_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,unsigned long long * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
    unsigned long long mem_allowed;

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
}

void adios_mpi_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,struct adios_method_struct * method
                    )
{
    v->data = buffer;
}

static void adios_mpi_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
    char name [STR_LEN];
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data; 
    struct adios_var_struct * v = 0;

    unsigned long long my_data_len = adios_data_size (fd->group);
    //printf("read: adios_data_size:%llu\n",my_data_len);
    unsigned long long my_offset = 0;
    MPI_Status mstatus;
    int amode = MPI_MODE_RDONLY;
    int err;
    MPI_Comm group_comm;
    int previous;
    int current;
    int next;
    int rank;
    int size;

    if (fd->group->group_comm)
    {
        adios_var_to_comm (fd->group->group_comm, fd->group->vars, &group_comm);

        MPI_Comm_rank (group_comm, &rank);
        MPI_Comm_size (group_comm, &size);
        if (rank == size - 1)
        {
            next = -1;
        }
        else
        {
            next = rank + 1;
        }
        previous = rank - 1;
        current = rank;
    }
    else
    {
        group_comm = MPI_COMM_NULL;
    }

    /* this code section is important if we want to write from
     * several processes to a single file.  For one file per process
     * writes, this is unnecessary.
     */
    /******************************
     * Begin Coordinate file offet
     ******************************/
    //char buf1 [STR_LEN];
    //char buf2 [STR_LEN];

    sprintf (name, "%s%s", method->base_path, fd->name);
    /* make a filename based on the size and current node */
    // sprintf (buf1, "%%s.%%0%dd", (int) (log10 (size)));
    // sprintf (buf2, buf1, fd->name, fd->group->current);
    if (group_comm != MPI_COMM_NULL)
    {
        if (previous == -1)
        {
            /* node 0 does an open */
            MPI_File_open (adios_mpi_comm_self, name, amode, adios_mpi_info, &md->fh);
            if (next != -1)
            {
                MPI_Isend (&my_data_len, 1, MPI_LONG_LONG, next
                          ,current, group_comm, &md->req
                          );
            }
        }
        else
        {
            MPI_Recv (&my_offset, 1, MPI_LONG_LONG, previous
                     ,previous, group_comm, &md->status
                     );
            if (next != -1)
            {
                unsigned long long new_offset = my_offset + my_data_len;
                MPI_Isend (&new_offset, 1, MPI_LONG_LONG, next
                          ,current, group_comm, &md->req
                          );
            }
            MPI_File_open (adios_mpi_comm_self, name, amode, adios_mpi_info
                          ,&md->fh
                          );
        }
        //printf ("previous=%d next=%d group_comm=%d fdoffset=%llu myoffset=%llu\n",previous,next,group_comm, fd->offset,my_offset);
        MPI_File_seek (md->fh, fd->base_offset + fd->offset + my_offset
                      ,MPI_SEEK_SET
                      );
    }
    else
    {
        MPI_File_open (adios_mpi_comm_self, name, amode, adios_mpi_info, &md->fh);
        MPI_File_seek (md->fh, fd->base_offset + fd->offset, MPI_SEEK_SET);
    }
    /******************************
     * End Coordinate file offet
     ******************************/

    if (md->buffer_size < my_data_len)
    {
        unsigned long long mem_needed = 0;
        unsigned long long mem_allowed = 0;
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
    err = MPI_File_read (md->fh, md->buffer, my_data_len, MPI_BYTE, &mstatus);
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

    if (next != -1)
        MPI_Wait (&md->req, &md->status);
}

static void adios_mpi_do_write (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                       method->method_data;
    char name [STR_LEN];
    unsigned long long my_data_len = adios_data_size (fd->group);
    unsigned long long mem_needed = 0;
    unsigned long long mem_allowed = 0;
    MPI_Offset my_offset = 0;
    MPI_Status mstatus;
    int overflow = 0;
    unsigned long long end = 0;
    struct adios_var_struct * v = fd->group->vars;
    struct adios_attribute_struct * a = fd->group->attributes;
    int err;
    MPI_Comm group_comm;
    int previous = -2;
    int current = -2;
    int next = -2;
    int rank = -2;
    int size = -2;

    if (fd->group->group_comm)
    {
        adios_var_to_comm (fd->group->group_comm, fd->group->vars, &group_comm);

        MPI_Comm_rank (group_comm, &rank);
        MPI_Comm_size (group_comm, &size);
        if (rank == size - 1)
        {
            next = -1;
        }
        else
        {
            next = rank + 1;
        }
        previous = rank - 1;
        current = rank;
    }
    else
    {
        group_comm = MPI_COMM_NULL;
    }

    /* this code section is important if we want to write from
     * several processes to a single file.  For one file per process
     * writes, this is unnecessary.
     */
    /******************************
     * Begin Coordinate file offet
     ******************************/
    //char buf1 [STR_LEN];
    //char buf2 [STR_LEN];

    sprintf (name, "%s%s", method->base_path, fd->name);
    /* make a filename based on the size and current node */
    // sprintf (buf1, "%%s.%%0%dd", (int) (log10 (size)));
    // sprintf (buf2, buf1, fd->name, fd->group->current);
    if (group_comm != MPI_COMM_NULL)
    {
        if (previous == -1)
        {
            int amode = MPI_MODE_WRONLY;
            if (fd->mode != adios_mode_append)
            {
                // truncate the file
                err = MPI_File_delete (name, adios_mpi_info);
                if (err != MPI_SUCCESS)
                {
                    fprintf (stderr, "Error truncating file %s for group %s\n"
                            ,name ,adios_file_mode_to_string (fd->mode)
                            );
                }
            }

            if (fd->mode != adios_mode_append)
                amode |= MPI_MODE_CREATE;
            /* node 0 does an open to create the file, if necessary */
            err = MPI_File_open (adios_mpi_comm_self, name, amode
                                ,adios_mpi_info, &md->fh
                                );
            if (err != MPI_SUCCESS)
            {
                amode |= MPI_MODE_CREATE;
                err = MPI_File_open (adios_mpi_comm_self, name, amode
                                    ,adios_mpi_info, &md->fh
                                    );
                if (err != MPI_SUCCESS)
                {
                    fprintf (stderr, "Error opening file %s for %s\n", name
                            ,adios_file_mode_to_string (fd->mode)
                            );
                }
            }
            if (next != -1)
            {
                MPI_Isend (&my_data_len, 1, MPI_LONG_LONG, next
                          ,current, group_comm, &md->req
                          );
            }
        }
        else
        {
            MPI_Recv (&my_offset, 1, MPI_LONG_LONG, previous
                     ,previous, group_comm
                     ,&md->status
                     );
            if (next != -1)
            {
                unsigned long long new_offset = my_offset + my_data_len;
                MPI_Isend (&new_offset, 1, MPI_LONG_LONG, next
                          ,current, group_comm, &md->req
                          );
            }
            MPI_File_open (adios_mpi_comm_self, name, MPI_MODE_WRONLY
                          ,adios_mpi_info, &md->fh
                          );
        }

        err = MPI_File_seek (md->fh, fd->base_offset + fd->offset + my_offset
                            ,MPI_SEEK_SET
                            );
        if (err != MPI_SUCCESS)
            fprintf (stderr, "cannot seek for some reason: %d\n", err);
    }
    else
    {
        int err_class;
        if (fd->mode != adios_mode_append)
        {
            // truncate the file
            err = MPI_File_delete (name, adios_mpi_info);
            MPI_Error_class (err, &err_class);
            if (err_class != MPI_SUCCESS)
            {
                fprintf (stderr, "Error %d truncating file %s for %s\n"
                        ,err_class
                        ,name, adios_file_mode_to_string (fd->mode)
                        );
            }
        }
        err = MPI_File_open (adios_mpi_comm_self, name
                            ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                            ,adios_mpi_info, &md->fh
                            );
        if (err != MPI_SUCCESS)
        {
            fprintf (stderr, "Error opening file %s for %s\n", name
                    ,adios_file_mode_to_string (fd->mode)
                    );
 
        }
        if (fd->mode != adios_mode_append)
        {
            err = MPI_File_get_size (md->fh, &my_offset);
            MPI_Error_class (err, &err_class);
            if (err_class != MPI_SUCCESS)
            {
                fprintf (stderr, "Error %d getting file size %s\n", err_class
                        ,name
                        );
            }
            //printf ("initial size: %llu\n", my_offset);
        }
        err = MPI_File_seek (md->fh, fd->base_offset + fd->offset, MPI_SEEK_SET);
        MPI_Error_class (err, &err_class);
        if (err_class != MPI_SUCCESS)
        {
            fprintf (stderr, "Error %d seeking file %s\n", err_class
                    ,name
                    );
        }
    }
    /******************************
     * End Coordinate file offet
     ******************************/

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
        long long file;
        char name [STR_LEN];

        sprintf (name, "%s%s_OVERFLOW%d", method->base_path, fd->name, rank);
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
//if (!fd->name [0] == 'T' || current == 0) 
        MPI_File_write (md->fh, md->buffer, end, MPI_BYTE, &mstatus);
    }

    if (group_comm != MPI_COMM_NULL && next != -1)
    {
        MPI_Wait (&md->req, &md->status);
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
}

void adios_mpi_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;

    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        adios_mpi_do_write (fd, method);
    }

    if (fd->mode == adios_mode_read)
    {
        adios_mpi_do_read (fd, method);
    }

    if (md && md->fh)
        MPI_File_close (&md->fh);

    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->start = 0;
    md->last_var_write_yes = 0;
}

void adios_mpi_finalize (int mype, struct adios_method_struct * method)
{
/* nothing to do here */
    if (adios_mpi_initialized)
        adios_mpi_initialized = 0;
}

void adios_mpi_end_iteration (struct adios_method_struct * method)
{
}

void adios_mpi_start_calculation (struct adios_method_struct * method)
{
}

void adios_mpi_stop_calculation (struct adios_method_struct * method)
{
}
