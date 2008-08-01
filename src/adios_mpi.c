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

static
void build_offsets (struct adios_bp_buffer_struct_v1 * b
                   ,MPI_Offset * offsets, int size, char * group_name
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

            offsets [pg_root->process_id * 2] = pg_root->offset_in_file;
            offsets [pg_root->process_id * 2 + 1] = size;
        }

        pg_root = pg_root->next;
    }
}

int adios_mpi_should_buffer (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            ,void * comm
                            )
{
    int i;
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                      method->method_data; 
    char * name;
    int err;
    int previous;
    int current;
    int next;

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
printf ("rank: %d mode_read\n", md->rank);
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

                    return 0;
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
printf ("*************\n");

                adios_init_buffer_read_index_offsets (&md->b);
                // already in the buffer
                adios_parse_index_offsets_v1 (&md->b);
printf ("*************\n");

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
printf ("*************\n");

                adios_init_buffer_read_vars_index (&md->b);
                MPI_File_seek (md->fh, md->b.vars_index_offset
                              ,MPI_SEEK_SET
                              );
                MPI_File_read (md->fh, md->b.buff, md->b.vars_size, MPI_BYTE
                              ,&md->status
                              );
                adios_parse_vars_index_v1 (&md->b, &md->old_vars_root);

printf ("*************\n");
                fd->base_offset = md->b.end_of_pgs;
            }
printf ("0 *************: %d\n", md->rank);

            if (md->group_comm != MPI_COMM_NULL && md->group_comm != MPI_COMM_SELF)
            {
printf ("1 *************: %d\n", md->rank);
                if (md->rank == 0)
                {
printf ("2 *************\n");
                    MPI_Offset * offsets = malloc (  sizeof (MPI_Offset)
                                                   * md->size * 2
                                                  );
                    memset (offsets, 0, sizeof (MPI_Offset) * md->size * 2);

                    // go through the pg index to build the offsets array
                    build_offsets (&md->b, offsets, md->size
                                  ,fd->group->name, md->old_pg_root
                                  );
                    MPI_Scatter (offsets, 2, MPI_LONG_LONG
                                ,offsets, 2, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
                    md->b.read_pg_offset = offsets [0];
                    md->b.read_pg_size = offsets [1];
                }
                else
                {
printf ("3 *************\n");
                    MPI_Offset offset [2];
                    offset [0] = offset [1] = 0;

printf ("3a *************\n");
                    MPI_Scatter (0, 0, MPI_LONG_LONG
                                ,offset, 2, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
printf ("3b *************\n");

                    md->b.read_pg_offset = offset [0];
                    md->b.read_pg_size = offset [1];
                }
            }
printf ("4 *************\n");

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
printf ("5 *************\n");
                // note rank 0 is already open
                // don't open it again here

                if (next != -1)
                {
                    MPI_Isend (&current, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
printf ("6 *************\n");
                MPI_Recv (&current, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&current, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_RDONLY
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
            }
printf ("7 *************\n");

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

                return 0;
            }
printf ("8 *************\n");

            break;
        }

        case adios_mode_write:
        {
            MPI_Offset * offsets = 0;
            MPI_Offset offset;

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
                    offsets [0] = 0;
                    for (i = 1; i < md->size - 1; i++)
                    {
                        offsets [i] = offsets [i - 1] + last_offset;
                        last_offset = offsets [i];
                    }
                    md->b.pg_index_offset =   offsets [md->size - 1]
                                            + last_offset;
for (i = 0; i < md->size; i++)
{
    printf ("offsets [%d]: %llu\n", i, offsets [i]);
}
printf ("index offset: %lld\n", md->b.pg_index_offset);
                    MPI_Scatter (offsets, 1, MPI_LONG_LONG
                                ,offsets, 1, MPI_LONG_LONG
                                ,0, md->group_comm
                                );
                    fd->base_offset = offsets [0];
                    free (offsets);
                }
                else
                {
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
                }
            }
            else
            {
printf ("here\n");
                md->b.pg_index_offset = fd->write_size_bytes;
            }
printf ("calculated index offset: %lld\n", md->b.pg_index_offset);

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                MPI_File_delete (name, MPI_INFO_NULL);  // make sure clean

                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY | MPI_MODE_CREATE
                                    ,MPI_INFO_NULL
                                    ,&md->fh
                                    );
                if (next != -1)
                {
                    MPI_Isend (&current, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&current, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&current, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
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

                return 0;
            }

            break;
        }

        case adios_mode_append:
        {
printf ("rank %d in append should buffer\n", md->rank);
            int old_file = 1;
            adios_buffer_struct_clear (&md->b);

            err = MPI_File_open (MPI_COMM_SELF, name, MPI_MODE_RDWR
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

                    return 0;
                }
            }

            if (old_file)
            {
printf ("rank %d in old_file\n", md->rank);
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
printf ("file size: %lld offset for version: %lld\n", md->b.file_size, md->b.length);
                    MPI_File_seek (md->fh, md->b.file_size - md->b.length
                                  ,MPI_SEEK_SET
                                  );
                    MPI_File_read (md->fh, md->b.buff, md->b.length, MPI_BYTE
                                  ,&md->status
                                  );
                    adios_parse_version (&md->b, &md->b.version);
printf ("version: %d\n", md->b.version);

                    adios_init_buffer_read_index_offsets (&md->b);
                    // already in the buffer
                    adios_parse_index_offsets_v1 (&md->b);
printf ("pg_index: %lld\n", md->b.pg_index_offset);
printf ("vars_index: %lld\n", md->b.vars_index_offset);

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

                    fd->base_offset = md->b.end_of_pgs;
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

printf ("ap write_size_bytes: %lld\n", fd->write_size_bytes);
                    uint64_t last_offset = offsets [0];
                    offsets [0] = fd->base_offset;
                    for (i = 1; i < md->size - 1; i++)
                    {
                        offsets [i] = offsets [i - 1] + last_offset;
                        last_offset = offsets [i];
                    }
                    md->b.pg_index_offset =   offsets [md->size - 1]
                                            + last_offset;
for (i = 0; i < md->size; i++)
{
    printf ("ap offsets [%d]: %llu\n", i, offsets [i]);
}
printf ("ap index offset: %lld\n", md->b.pg_index_offset);
                        MPI_Scatter (offsets, 1, MPI_LONG_LONG
                                    ,offsets, 1, MPI_LONG_LONG
                                    ,0, md->group_comm
                                    );
                        fd->base_offset = offsets [0];
                    }
                    else
                    {
                        MPI_Offset offset = fd->write_size_bytes;
printf ("rank: %d write size: %lld\n", md->rank, offset);

                        MPI_Gather (&offset, 1, MPI_LONG_LONG
                                   ,&offset, 1, MPI_LONG_LONG
                                   ,0, md->group_comm
                                   );

                        MPI_Scatter (&offset, 1, MPI_LONG_LONG
                                    ,&offset, 1, MPI_LONG_LONG
                                    ,0, md->group_comm
                                    );

                        fd->base_offset = offset;
                    }
                }
            }
            else
            {
                fd->base_offset = 0;
            }

            // cascade the opens to avoid trashing the metadata server
            if (previous == -1)
            {
                // already open on node 0 so don't need to open again

                if (next != -1)
                {
                    MPI_Isend (&current, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
            }
            else
            {
                MPI_Recv (&current, 1, MPI_INTEGER, previous, previous
                         ,md->group_comm, &md->status
                         );
                if (next != -1)
                {
                    MPI_Isend (&current, 1, MPI_INTEGER, next, current
                              ,md->group_comm, &md->req
                              );
                }
                err = MPI_File_open (MPI_COMM_SELF, name
                                    ,MPI_MODE_WRONLY
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

                return 0;
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

    switch (md->b.version)
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

            int i;

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

static void adios_mpi_do_write (struct adios_file_struct * fd
                               ,struct adios_method_struct * method
                               ,char * buffer
                               ,uint64_t buffer_size
                               )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                 method->method_data;

printf ("write offset: %lld\n", fd->base_offset);
    MPI_File_seek (md->fh, fd->base_offset, MPI_SEEK_SET);
printf ("writing data: %lld\n", fd->bytes_written);
    MPI_File_write (md->fh, fd->buffer, fd->bytes_written, MPI_BYTE
                   ,&md->status
                   );
printf ("write index offset: %lld\n", md->b.pg_index_offset);
    MPI_File_seek (md->fh, md->b.pg_index_offset, MPI_SEEK_SET);
printf ("writing index: %lld\n", buffer_size);
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
            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    int ranks_sent = 1; // assume that we have sent to ourselves
                    buffer_size = 100 * 1024; // try 100k to start
                    buffer = malloc (buffer_size);
                    buffer_offset = 0;
                    int count = 0;

                    while (ranks_sent < md->size)
                    {
                        MPI_Recv (buffer, buffer_size, MPI_BYTE, MPI_ANY_SOURCE
                                 ,0, md->group_comm, &md->status
                                 );
                        MPI_Get_count (&md->status, MPI_BYTE, &count);
                        if (buffer_size <= count)
                        {
                            fprintf (stderr, "100k buffer size too small.\n");
                        }
                        else
                        {
                            char * buffer_save = md->b.buff;
                            uint64_t buffer_size_save = md->b.length;
                            uint64_t offset_save = md->b.offset;
    
                            md->b.buff = buffer;
                            md->b.length = count;
                            md->b.offset = 0;
    
                            adios_parse_process_group_index_v1 (&md->b
                                                               ,&new_pg_root
                                                               );
                            adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                            adios_merge_index_v1 (&md->old_pg_root
                                                 ,&md->old_vars_root
                                                 ,new_pg_root, new_vars_root
                                                 );
                            new_pg_root = 0;
                            new_vars_root = 0;
                            md->b.buff = buffer_save;
                            md->b.length = buffer_size_save;
                            md->b.offset = offset_save;
                        }
                        ranks_sent++;
		    }
                }
                else
                {
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,index_start, new_pg_root
                                         ,new_vars_root
                                         );

                    MPI_Send (buffer, buffer_offset, MPI_BYTE, 0, 0
                             ,md->group_comm
                             );
                }
            }
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
printf ("start append write\n");
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;
            uint64_t index_start = md->b.pg_index_offset;

            // build index
            adios_build_index_v1 (fd, &new_pg_root, &new_vars_root);
            // merge in old indicies
            adios_merge_index_v1 (&md->old_pg_root, &md->old_vars_root
                                 ,new_pg_root, new_vars_root);
            new_pg_root = 0;
            new_vars_root = 0;
            if (md->group_comm != MPI_COMM_NULL)
            {
                if (md->rank == 0)
                {
                    int ranks_sent = 1;
                    buffer_size = 100 * 1024; // try 100k to start
                    buffer = malloc (buffer_size);
                    buffer_offset = 0;
                    int count = 0;

                    while (ranks_sent < md->size)
                    {
                        MPI_Recv (buffer, buffer_size, MPI_BYTE, MPI_ANY_SOURCE
                                 ,0, md->group_comm, &md->status
                                 );
                        MPI_Get_count (&md->status, MPI_BYTE, &count);
                        if (buffer_size <= count)
                        {
                            fprintf (stderr, "100k buffer size too small.\n");
                        }
                        else
                        {
                            char * buffer_save = md->b.buff;
                            uint64_t buffer_size_save = md->b.length;
                            uint64_t offset_save = md->b.offset;
    
                            md->b.buff = buffer;
                            md->b.length = count;
                            md->b.offset = 0;
    
                            adios_parse_process_group_index_v1 (&md->b
                                                               ,&new_pg_root
                                                               );
                            adios_parse_vars_index_v1 (&md->b, &new_vars_root);
                            adios_merge_index_v1 (&md->old_pg_root
                                                 ,&md->old_vars_root
                                                 ,new_pg_root, new_vars_root
                                                 );
                            new_pg_root = 0;
                            new_vars_root = 0;
                            md->b.buff = buffer_save;
                            md->b.length = buffer_size_save;
                            md->b.offset = offset_save;
                        }
                        ranks_sent++;
		    }
                }
                else
                {
                    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                         ,index_start, new_pg_root
                                         ,new_vars_root
                                         );

                    MPI_Send (buffer, buffer_offset, MPI_BYTE, 0, 0
                             ,md->group_comm
                             );
                    adios_clear_index_v1 (new_pg_root, new_vars_root);
                    new_pg_root = 0;
                    new_vars_root = 0;
                }
            }
            adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                 ,index_start, md->old_pg_root
                                 ,md->old_vars_root
                                 );
            adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
            adios_mpi_do_write (fd, method, buffer, buffer_offset);

            free (buffer);

            break;
printf ("end append write\n");
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

printf ("before close\n");
    if (md && md->fh)
        MPI_File_close (&md->fh);

printf ("before comm free\n");
    if (   md->group_comm != MPI_COMM_WORLD
        && md->group_comm != MPI_COMM_SELF
        && md->group_comm != MPI_COMM_NULL
       )
        MPI_Comm_free (&md->group_comm);

    md->fh = 0;
    md->req = 0;
    memset (&md->status, 0, sizeof (MPI_Status));
    md->group_comm = MPI_COMM_NULL;

printf ("before clear index\n");
    adios_clear_index_v1 (md->old_pg_root, md->old_vars_root);
    md->old_pg_root = 0;
    md->old_vars_root = 0;
printf ("end of close\n");
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
