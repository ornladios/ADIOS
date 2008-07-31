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
#include "adios_bp_v1.h"
#include "adios_internals.h"

static int adios_mpi_initialized = 0;

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
};

static void adios_var_to_comm (enum ADIOS_FLAG host_language_fortran
                              ,void * data
                              ,MPI_Comm * comm
                              )
{
    if (data)
    {
        int t = *(int *) data;
        if (host_language_fortran == adios_flag_yes)
        {
            *comm = MPI_Comm_f2c (t);
        }
        else
        {
            *comm = *(MPI_Comm *) data;
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

void adios_mpi_init (const char * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
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
    md->group_comm = MPI_COMM_NULL;
    md->old_pg_root = 0;
    md->old_vars_root = 0;

    adios_buffer_struct_init (&md->b);
}

int adios_mpi_open (struct adios_file_struct * fd
                   ,struct adios_method_struct * method
                   )
{
    char name [STR_LEN];
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

    // we have to wait for the group_size (should_buffer) to get the comm
    // before we can do an open for any of the modes

    return 1;
}

int adios_mpi_should_buffer (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            ,void * comm
                            )
{
    uint32_t version = 0;
    int i;
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data; 
    char * name;
    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf (name, "%s%s", method->base_path, fd->name);

    adios_var_to_comm (fd->group->adios_host_language_fortran
                      ,comm
                      ,&md->group_comm
                      );
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }

#define DO_SCATTER 1
    switch (fd->mode)
    {
        case adios_mode_read:
        {
            int err;

            if (md->group_comm == MPI_COMM_NULL || md->rank == 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL, &md->fh
                                    );
                if (err != MPI_SUCCESS)
                {
                    free (name);

                    return 0;
                }

                MPI_Offset file_size;
                MPI_File_get_size (md->fh, &file_size);
                md->b.file_size = file_size;

                adios_init_buffer_read_version (&md->b);
                MPI_File_seek (md->fh, md->b.file_size - md->b.length
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.length
                              ,&md->status
                              );
                adios_parse_version (&md->b, &version);

                adios_init_buffer_read_index_offsets_v1 (&md->b);
                // already in the buffer
                adios_parse_index_offsets_v1 (&md->b);

                adios_init_buffer_read_process_group_index (&md->b);
                MPI_File_seek (md->fh, md->b.pg_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.pg_size
                              ,&md->status
                              );
                adios_parse_process_group_index_v1 (&md->b
                                                   ,&md->old_pg_root
                                                   );

                adios_init_buffer_read_vars_index (&md->b);
                MPI_File_seek (md->fh, md->b.vars_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.vars_size
                              ,&md->status
                              );
                adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

                fd->base_offset = md->b.end_of_pgs;
            }

            if (1) //md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size
                                                  );
                    memset (offsets, 0, md->size * sizeof (MPI_Offset));

                    // go through the pg index to build the offsets array
printf ("build the base offsets\n");

#if DO_SCATTER
                    MPI_Scatter (offsets, 1, MPI_LONG_LONG
                                ,offsets, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
#endif
                    fd->base_offset = offsets [0];
                }
                else
                {
                    MPI_Offset offset = 0;

#if DO_SCATTER
                    MPI_Scatter (&offset, 1, MPI_LONG_LONG
                                ,&offset, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
#endif

                    fd->base_offset = offset;
                }
            }
printf ("do the cascading open (MPI open read)\n");
            // note rank 0 is already open
            if (md->rank != 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL, &md->fh
                                    );
            }
            break;
        }

        case adios_mode_write:
        {
            int err;

            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size
                                                  );

                    MPI_Gather (offsets, 1, MPI_LONG_LONG
                               ,offsets, 1, MPI_LONG_LONG
                               ,0, md->group_comm
                               );

                    for (i = 0; i < md->size - 1; i++)
                    {
                        offsets [i + 1] += offsets [i];
                    }
#if DO_SCATTER
                    MPI_Scatter (offsets, 1, MPI_LONG_LONG
                                ,offsets, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
#endif
                    fd->base_offset = offsets [0];
                    free (offsets);
                }
                else
                {
                    MPI_Offset offset = fd->write_size_bytes;

                    MPI_Gather (&offset, 1, MPI_LONG_LONG
                               ,&offset, 1, MPI_LONG_LONG
                               ,0, md->group_comm
                               );

#if DO_SCATTER
                    MPI_Scatter (&offset, 1, MPI_LONG_LONG
                                ,&offset, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
#endif
                    fd->base_offset = offset;
                }
            }

printf ("do the cascading open thing (MPI open write)\n");
            err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_WRONLY
                                ,MPI_INFO_NULL, &md->fh
                                );
            if (err != MPI_SUCCESS)
            {
                free (name);

                return 0;
            }

            break;
        }

        case adios_mode_append:
        {
            adios_buffer_struct_clear (&md->b);

            int err;
            err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDWR
                                ,MPI_INFO_NULL, &md->fh
                                );

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
                MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.length
                              ,&md->status
                              );
                adios_parse_version (&md->b, &version);

                adios_init_buffer_read_index_offsets_v1 (&md->b);
                // already in the buffer
                adios_parse_index_offsets_v1 (&md->b);

                adios_init_buffer_read_process_group_index (&md->b);
                MPI_File_seek (md->fh, md->b.pg_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.pg_size
                              ,&md->status
                              );
                adios_parse_process_group_index_v1 (&md->b
                                                   ,&md->old_pg_root
                                                   );

                adios_init_buffer_read_vars_index (&md->b);
                MPI_File_seek (md->fh, md->b.vars_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.vars_size
                              ,&md->status
                              );
                adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

                fd->base_offset = md->b.end_of_pgs;
            }

            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size
                                                  );

                    MPI_Gather (offsets, 1, MPI_LONG_LONG
                               ,offsets, 1, MPI_LONG_LONG
                               ,0, md->group_comm
                               );

                    offsets [0] = fd->base_offset;
                    for (i = 0; i < md->size - 1; i++)
                    {
                        offsets [i + 1] += offsets [i];
                    }
#if DO_SCATTER
                    MPI_Scatter (offsets, 1, MPI_LONG_LONG
                                ,offsets, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
#endif
                }
                else
                {
                    MPI_Offset offset = fd->write_size_bytes;

                    MPI_Gather (&offset, 1, MPI_LONG_LONG
                               ,&offset, 1, MPI_LONG_LONG
                               ,0, md->group_comm
                               );

#if DO_SCATTER
                    MPI_Scatter (&offset, 1, MPI_LONG_LONG
                                ,&offset, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
#endif
                }
            }
printf ("do the cascading open thing (MPI open append)\n");
            if (md->rank != 0)
            {
                err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDWR
                                    ,MPI_INFO_NULL, &md->fh
                                    );
            }

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);

            free (name);

            return 0;
        }
    }

    free (name);

    return 1;   // as far as we care, buffer
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
}

void adios_mpi_get_write_buffer (struct adios_file_struct * fd
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
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data; 
    struct adios_var_struct * v = fd->group->vars;

    uint64_t element_size = 0;
    struct adios_bp_element_struct * element = NULL;
    struct adios_parse_buffer_struct data;

    data.vars = v;
    data.buffer = 0;
    data.buffer_len = 0;

    uint32_t version = 0;

    adios_init_buffer_read_version (&md->b);
    MPI_File_seek (md->fh, md->b.file_size - md->b.length
                  ,MPI_SEEK_SET
                  );
    MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.length
                  ,&md->status
                  );
    adios_parse_version (&md->b, &version);

    switch (version)
    {
        case 1:
        {
            struct adios_index_process_group_struct_v1 * pg_root = 0;
            struct adios_index_var_struct_v1 * vars_root = 0;

            adios_init_buffer_read_version (&md->b);
            MPI_File_seek (md->fh, md->b.file_size - md->b.length
                          ,MPI_SEEK_SET
                          );
            MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.length
                          ,&md->status
                          );
            adios_parse_version (&md->b, &version);
                               
            adios_init_buffer_read_index_offsets_v1 (&md->b);
            // already in the buffer
            adios_parse_index_offsets_v1 (&md->b);

            adios_init_buffer_read_process_group_index (&md->b);
            MPI_File_seek (md->fh, md->b.pg_index_offset
                          ,MPI_SEEK_SET
                          );
            MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.pg_size
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
            MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.vars_size
                          ,&md->status
                          );
            adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);
#endif

            // the three section headers
            struct adios_process_group_header_struct_v1 pg_header;
            struct adios_vars_header_struct_v1 vars_header;
            struct adios_attributes_header_struct_v1 attrs_header;

            struct adios_var_header_struct_v1 var_header;
            struct adios_var_payload_struct_v1 var_payload;
            struct adios_attribute_struct_v1 attribute;

            int i;

            adios_init_buffer_read_process_group (&md->b, 1, pg_root);
            MPI_File_seek (md->fh, md->b.read_pg_offset
                          ,MPI_SEEK_SET
                          );
            MPI_File_read (md->fh, md->b.buff, MPI_BYTE, md->b.read_pg_size
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
                        break;
                }

                if (v1)
                {
                    var_payload.payload = v1->data;
                    adios_parse_var_data_payload_v1 (&md->b, &var_header
                                                    ,&var_payload
                                                    );
                }
                else
                {
                    adios_parse_var_data_payload_v1 (&md->b, &var_header, NULL);
                }
            }

#if 1
            adios_parse_attributes_header_v1 (&md->b, &attrs_header);

            for (i = 0; i < attrs_header.count; i++)
            {
                adios_parse_attribute_v1 (&md->b, &attribute);
            }
#endif
            adios_clear_index_v1 (pg_root, vars_root);
            break;
        }

        default:
            fprintf (stderr, "MPI read: file version unknown: %u\n"
                    ,version
                    );
            return;
    }

    adios_buffer_struct_clear (&md->b);
}

static void adios_mpi_do_write (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               ,char * buffer
                               ,uint64_t buffer_size
                               )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;

    MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
    MPI_File_write (md->fh, fd->buffer, fd->bytes_written, MPI_BYTE, &md->status);
    MPI_File_write (md->fh, buffer, buffer_size, MPI_BYTE, &md->status);
}

void adios_mpi_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;

    switch (fd->mode)
    {
        case adios_mode_write:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = fd->base_offset + fd->offset;

            // build index
            adios_build_index_v1 (fd, &new_pg_root, &new_vars_root);
            // if collective, gather the indexes from the rest and call
            // adios_merge_index_v1 (&new_pg_root, &new_vars_root, pg, vars);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, new_pg_root, new_vars_root
                                 );
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_mpi_do_write (fd, method, buffer, buffer_offset);

            free (buffer);

            adios_clear_index_v1 (new_pg_root, new_vars_root);

            break;
        }

        case adios_mode_append:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = fd->base_offset + fd->offset;

            // build index
            adios_build_index_v1 (fd, &new_pg_root, &new_vars_root);
            // merge in old indicies
            adios_merge_index_v1 (&md->old_pg_root, &md->old_vars_root
                                 ,new_pg_root, new_vars_root);
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, md->old_pg_root, md->old_vars_root
                                 );
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_mpi_do_write (fd, method, buffer, buffer_offset);

            free (buffer);

            break;
        }

        case adios_mode_read:
        {
            // read the index to find the place to start reading
            adios_mpi_do_read (fd, method);
            struct adios_var_struct * v = fd->group->vars;
            while (v)
            {
                v->data = 0;
                v = v->next;
            }

            break;
        }

        default:
        {
            fprintf (stderr, "Unknown file mode: %d\n", fd->mode);
        }
    }

    if (md && md->fh)
        MPI_File_close (&md->fh);

    MPI_Comm_free (&md->group_comm);

    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->group_comm = MPI_COMM_NULL;

    adios_clear_index_v1 (md->old_pg_root, md->old_vars_root);
    md->old_pg_root = 0;
    md->old_vars_root = 0;
}

void adios_mpi_finalize (int mype, struct adios_method_struct * method)
{
// nothing to do here
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
