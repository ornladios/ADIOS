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
#include <assert.h>

// xml parser
#include <mxml.h>

// see if we have MPI or other tools
#include "config.h"

#include "public/adios.h"
#include "public/adios_types.h"
#include "public/adios_mpi.h"
#include "public/adios_error.h"
#include "core/common_adios.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/adios_internals_mxml.h"
#include "core/buffer.h"
#include "core/util.h"
#include "core/adios_logger.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

extern int adios_errno; // in adios_error.c
extern struct adios_transport_struct * adios_transports; // in adios_internals.c

struct adios_staging_data_struct
{
    struct adios_file_struct * stagedf; // handle returned by selected method's adios_open() function
    int rank;   // dataspaces rank or MPI rank if MPI is available
    int size;
#if HAVE_MPI
    MPI_Comm group_comm; // for use in open..close
    MPI_Comm mpi_comm_init; // for use in init/finalize
#endif
    char * method_name; // actual staging method (dataspaces, etc.) used for the group
    char * params; // parameter string passed down to the actual staging method
    int method_initialized;
    uint64_t pg_global_size; // total size of data written by ALL processes
    struct adios_group_struct * group; // a copy of the original group to pass down
    struct adios_bp_buffer_struct_v1 b;
    struct adios_index_struct_v1 * index;

};

static struct adios_staging_data_struct * create_data_struct (struct adios_method_struct * method)
{
    struct adios_staging_data_struct * md = malloc (sizeof (struct adios_staging_data_struct));
    if (md != NULL)
    {
        md->rank = 0;
        md->size = 1;
#if HAVE_MPI
        md->mpi_comm_init = method->init_comm; // unused here, adios_open will set the current comm
        md->group_comm = method->init_comm; // temporarily set to a communicator, so that
                                            // adios_common_select_method_by_group_id() will
                                            // find a required_communicator
#endif
        md->method_name = NULL;
        md->params = NULL;
        md->method_initialized = 0;
        md->index = adios_alloc_index_v1(1); // with hashtables
        adios_buffer_struct_init (&md->b);
    }
    return md;
}

static void destroy_data_struct (struct adios_staging_data_struct ** md)
{
    adios_free_index_v1 ((*md)->index);
    adios_buffer_struct_clear (&(*md)->b);
    if ((*md)->method_name)
        free ((*md)->method_name);
    if ((*md)->params)
        free ((*md)->params);
    free (*md);
    *md = NULL;
}

/*static struct adios_group_struct * duplicate_group (const struct adios_group_struct * const gin, const char * const name)
{
    struct adios_group_struct * g = (struct adios_group_struct *)
        malloc (sizeof (struct adios_group_struct));
    if (g != NULL) {
        memcpy (g, gin, sizeof(struct adios_group_struct));
        g->name = strdup(name);
    }
    return g;
}*/

static void init_output_parameters(struct adios_method_struct * method, const PairStruct *params)
{
    const PairStruct *p = params;
    struct adios_staging_data_struct * md = (struct adios_staging_data_struct *)
                                                        method->method_data;
    while (p) {
        if (!strcasecmp (p->name, "method")) {
            errno = 0;
            md->method_name = strdup(p->value);
            // method->group is NULL at this point
            log_debug ("method set to %s for staging\n", md->method_name);
        } else if (!strcasecmp (p->name, "parameters")) {
            errno = 0;
            md->params = strdup(p->value);
            log_debug ("parameters set to %s for staging\n", md->params);
        } else {
            log_error ("Parameter name %s is not recognized by the staging method\n", p->name);
        }
        p = p->next;
    }
}

static int is_supported_method (const char * const method_name)
{
    if (!strcmp(method_name, "NULL") ||
        !strcmp(method_name, "DATASPACES") ||
        !strcmp(method_name, "DIMES")      ||
        !strcmp(method_name, "FLEXPATH")   ||
        !strcmp(method_name, "ICEE")
       )
    {
       return 1;
    }
    return 0;
}

void adios_staging_init (const PairStruct * parameters
                    ,struct adios_method_struct * method
                    )
{
    struct adios_staging_data_struct * md = create_data_struct(method);
    method->method_data = (void*) md;
    init_output_parameters(method, parameters);
    if (md->method_name && !is_supported_method(md->method_name)) {
        adios_error (err_unspecified, "Selected I/O method %s is not supported by the staging method\n", md->method_name);
        free(md->method_name);
        md->method_name = strdup("NULL");
    }
    if (!md->method_name) {
        /*int i;
        for (i = 0; i < ADIOS_METHOD_COUNT; i++) {
            if (is_supported_method (adios_transports[i].method_name))
            {
                md->method = strdup(adios_transports[i].method_name);
            }
        }*/
        md->method_name = strdup("MPI");
        if (md->method_name) {
            log_warn ("Actual staging method is not given for the staging method. "
                     "Will use %s\n", md->method_name);
        } else {
            adios_error (err_unspecified, "Actual staging method is not given for the staging method "
                    "and there is no supported method available in this ADIOS build\n");
        }
    }
    if (!md->params) {
        md->params = strdup ("");
    }
}

int adios_staging_open (struct adios_file_struct * fd,
                        struct adios_method_struct * method,
                        MPI_Comm comm
                       )
{
    struct adios_staging_data_struct * md = (struct adios_staging_data_struct *)
                                                    method->method_data;

    if (!md->method_name)
        return 0;

    adios_buffer_struct_clear (&md->b);

#if HAVE_MPI
    md->group_comm = comm;
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }
#endif

    fd->group->process_id = md->rank;

    if (!md->method_initialized)
    {
        char *gname = malloc(strlen(method->group->name)+2);
        gname[0] = '$';
        strcpy (gname+1, method->group->name);
        int64_t g;
        adios_common_declare_group (&g, gname, method->group->adios_host_language_fortran,
                                    "", "", "", method->group->stats_flag);
        adios_common_select_method_by_group_id (method->priority,
                                                md->method_name,
                                                md->params,
                                                g,
                                                method->base_path,
                                                method->iterations);
        md->group = (struct adios_group_struct *) g;
        /* The function above calls the selected method's adios_init() function */

        /* Flexpath needs to have the group's variables defined by the time of open,
         * so we need to define them here, instead postponing it to adios_close()
         */
        adios_common_define_var ((int64_t)md->group, "pg_ldim", "", adios_long, 0, 0, 0);
        adios_common_define_var ((int64_t)md->group, "pg_gdim", "", adios_long, 0, 0, 0);
        adios_common_define_var ((int64_t)md->group, "pg_offs", "", adios_long, 0, 0, 0);

        adios_common_define_var ((int64_t)md->group, "pg", "", adios_byte,
                                           "pg_ldim", "pg_gdim", "pg_offs");

        if (md->rank == 0)
        {
            adios_common_define_var ((int64_t)md->group, "index_size", "", adios_long, 0, 0, 0);

            adios_common_define_var ((int64_t)md->group, "index", "", adios_byte,
                                     "index_size", "index_size", "0");
        }


        md->method_initialized = 1;
        free (gname);
    }

    switch (fd->mode)
    {
        case adios_mode_write:
        case adios_mode_append:
        case adios_mode_update:
        {
            char file_mode[] = "w";
            if (fd->mode == adios_mode_append)
                file_mode[0] = 'a';
            else if (fd->mode == adios_mode_update)
                file_mode[0] = 'u';

            md->stagedf = common_adios_open (md->group->name, fd->name, file_mode, md->group_comm);
            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode, 
                         "Staging method: Unknown file mode requested: %d\n",
                         fd->mode);
            return 0;
        }
    }
    return 1;
}

static void build_file_offsets (struct adios_staging_data_struct *md,
                                struct adios_file_struct *fd)
{
    if (md->group_comm != MPI_COMM_NULL)
    {
        if (md->rank == 0)
        {
            // make one space for offset and one for size
            MPI_Offset * offsets = malloc(sizeof (MPI_Offset) * md->size);
            int i;

            offsets [0] = fd->bytes_written; // = the size of data in buffer on this processor
// mpixlc_r on Eugene doesn't support 64 bit mode. Therefore the following may have problem
// on Eugene for large data size since MPI_LONG_LONG is 32bit 
            MPI_Gather (&(fd->bytes_written), 1, MPI_LONG_LONG
                       ,offsets, 1, MPI_LONG_LONG
                       ,0, md->group_comm);


            uint64_t last_pgsize = offsets [0];
            offsets [0] = md->b.end_of_pgs; // 0 or where to append to existing data
            //printf ("offsets[%d] = {%llu", md->size, offsets[0]);
            for (i = 1; i < md->size; i++)
            {
                uint64_t this_offset = offsets [i];
                offsets [i] = offsets [i - 1] + last_pgsize;
                last_pgsize = this_offset;
                //printf ("  %llu", offsets[i]);
            }
            //printf (" }\n");
            md->b.pg_index_offset =   offsets [md->size - 1]
                                    + last_pgsize;
            //printf (" last_pgsize = %llu, pg_index_offset = %llu\n", last_pgsize, md->b.pg_index_offset);

            MPI_Scatter (offsets, 1, MPI_LONG_LONG
                        ,MPI_IN_PLACE, 1, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            md->pg_global_size = md->b.pg_index_offset;
            fd->current_pg->pg_start_in_file = offsets[0];
            free (offsets);
        }
        else
        {
            MPI_Offset offset[1];
            offset[0] = fd->bytes_written;

            MPI_Gather (offset, 1, MPI_LONG_LONG
                       ,0, 1, MPI_LONG_LONG
                       ,0, md->group_comm
                       );

            MPI_Scatter (0, 1, MPI_LONG_LONG
                        ,offset, 1, MPI_LONG_LONG
                        ,0, md->group_comm
                        );
            fd->current_pg->pg_start_in_file = offset[0];
        }
        MPI_Bcast(&md->pg_global_size, 1, MPI_LONG_LONG, 0, md->group_comm);
        //printf ("rank %d: pg_start_in_file = %llu\n", md->rank, fd->current_pg->pg_start_in_file);
    }
    else
    {
        md->b.pg_index_offset = fd->bytes_written;
        md->pg_global_size = md->b.pg_index_offset;
        fd->current_pg->pg_start_in_file = md->b.end_of_pgs; // 0 or where to append to existing data
    }
}

enum BUFFERING_STRATEGY adios_staging_should_buffer (struct adios_file_struct * fd
                                                ,struct adios_method_struct * method
                                                )
{
    return stop_on_overflow;
}

void adios_staging_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,const void * data
                     ,struct adios_method_struct * method
                     )
{
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
}

void adios_staging_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
    uint64_t mem_allowed;
    struct adios_staging_data_struct * md = (struct adios_staging_data_struct *)
                                                      method->method_data;

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
                         "MPI method, rank %d: cannot allocate %llu bytes for variable %s\n",
                         md->rank, *size, v->name);
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
                "MPI method, rank %d: OVERFLOW: Cannot get requested ADIOS buffer of %llu "
                "bytes for variable %s. Remaining buffer size was %llu\n",
                md->rank, *size, v->name, mem_allowed);
        *size = 0;
        *buffer = 0;
    }
}

void adios_staging_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )
{
    return;
}

void adios_staging_buffer_overflow (struct adios_file_struct * fd,
                                struct adios_method_struct * method)
{
    struct adios_staging_data_struct * md = (struct adios_staging_data_struct *)
                                                 method->method_data;
    adios_error (err_buffer_overflow, 
            "rank %d: MPI method only works with complete buffering of data between adios_open() "
            "and adios_close(). Portions of global arrays, that do not fit into the "
            "buffer on some processors will not be written by this method to %s\n", 
            md->rank, fd->name);
}


void adios_staging_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_staging_data_struct * md = (struct adios_staging_data_struct *)
                                                 method->method_data;
    if (!md->method_name)
        return;

    struct adios_index_process_group_struct_v1 * new_pg_root = 0;
    struct adios_index_var_struct_v1 * new_vars_root = 0;
    struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

    switch (fd->mode)
    {
        case adios_mode_write:
        case adios_mode_append:
        case adios_mode_update:
        {
            char * buffer = 0;
            uint64_t buffer_size = 0;
            uint64_t buffer_offset = 0;

            // figure out the offsets
            // before writing out the buffer and build the index based on target offsets
            build_file_offsets (md, fd);

            // build index appending to any existing index
            adios_build_index_v1 (fd, md->index);

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

                        adios_parse_process_group_index_v1 (&md->b, &new_pg_root, NULL);
                        adios_parse_vars_index_v1 (&md->b, &new_vars_root, NULL, NULL);
                        // do not merge attributes from other processes from 1.4
                        /*
                        adios_parse_attributes_index_v1 (&md->b
                                                        ,&new_attrs_root
                                                        );
                        */
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

                    uint32_t tmp_buffer_size = (uint32_t) buffer_offset;
                    MPI_Gather (&tmp_buffer_size, 1, MPI_INT, 0, 0, MPI_INT
                               ,0, md->group_comm
                               );
                    MPI_Gatherv (buffer, buffer_offset, MPI_BYTE
                                ,0, 0, 0, MPI_BYTE
                                ,0, md->group_comm
                                );
                    free (buffer);
                }
            }

            // everyone writes the BP buffer as their variable now
            /*
            char ldim[32], gdim[32], offset[32];
            sprintf (ldim,   "%"PRIu64, fd->bytes_written);
            sprintf (gdim,   "%"PRIu64, md->pg_global_size);
            sprintf (offset, "%"PRIu64, fd->current_pg->pg_start_in_file);

            int64_t var = adios_common_define_var ((int64_t)md->group, "pg", "", adios_unsigned_byte,
                                               ldim, gdim, offset);
            common_adios_write_byid (md->stagedf, (struct adios_var_struct *)var, fd->buffer);

             */
            struct adios_var_struct *v;
            v = adios_find_var_by_name (md->group, "pg_ldim");
            common_adios_write_byid (md->stagedf, v, &fd->bytes_written);
            v = adios_find_var_by_name (md->group, "pg_gdim");
            common_adios_write_byid (md->stagedf, v, &md->pg_global_size);
            v = adios_find_var_by_name (md->group, "pg_offs");
            common_adios_write_byid (md->stagedf, v, &fd->current_pg->pg_start_in_file) ;

            v = adios_find_var_by_name (md->group, "pg");
            common_adios_write_byid (md->stagedf, v, fd->buffer);


            /* Rank 0 writes the index */
            if (md->rank == 0)
            {
                adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
                                     ,md->b.pg_index_offset, md->index);
                adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

                /*
                sprintf (ldim,   "%"PRIu64, buffer_offset);
                sprintf (gdim,   "%"PRIu64, buffer_offset);
                int64_t idx = adios_common_define_var ((int64_t)md->group, "index", "",
                                                       adios_unsigned_byte, ldim, gdim, "0");
                common_adios_write_byid (md->stagedf, (struct adios_var_struct *)idx, buffer);
                */
                v = adios_find_var_by_name (md->group, "index_size");
                common_adios_write_byid (md->stagedf, v, &buffer_offset);
                v = adios_find_var_by_name (md->group, "index");
                common_adios_write_byid (md->stagedf, v, buffer);
                free (buffer);
            }
            break;
        }

        default:
        {
            adios_error (err_invalid_file_mode,
                    "Staging method: Unknown file mode: %d in adios_close()\n",
                    fd->mode);
        }
    }

    common_adios_close(md->stagedf);
    md->stagedf = NULL;
    md->group_comm = MPI_COMM_NULL;

    adios_clear_index_v1 (md->index);
}

void adios_staging_finalize (int mype, struct adios_method_struct * method)
{
    struct adios_staging_data_struct * md = (struct adios_staging_data_struct *)
                                                    method->method_data;
    if (!md->method_name)
        return;

    // common_adios_finalize() calls the actual method's finalize itself
    // BUT after this one because adios_select_method() at open added it to the
    // end of the list of existing methods
    adios_common_free_group ((int64_t)md->group);
    destroy_data_struct (&md);
    method->method_data = NULL;
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
