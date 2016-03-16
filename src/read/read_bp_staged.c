/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 *
 * Author: Qing Liu
 */


/******************************/
/* A read method for BP files */
/******************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <assert.h>
#include <errno.h>
#include <sys/time.h>
#include "public/adios_types.h"
#include "public/adios_read.h"
#include "public/adios_error.h"
#include "core/bp_utils.h"
#include "core/bp_types.h"
#include "core/adios_read_hooks.h"
#include "core/futils.h"
#include "core/common_read.h"
#include "core/adios_logger.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define TAG_DATA 0

static int chunk_buffer_size = -1;  // chunk size - must be set
static int poll_interval = 10; // 10 secs by default
static int show_hidden_attrs = 0; // don't show hidden attr by default
static int num_aggregators = -1; // number of aggregator - must be set

//private structure of BP_PROC
typedef struct bp_proc_private_struct
{
    int rank;
    int num_aggregators;
    int new_rank;
    int size;
    int groups;
    int group_size;
    int group;
    MPI_Comm group_comm;
    MPI_Comm new_comm;
    MPI_Comm new_comm2;
    int aggregator_rank;
    int aggregator_new_rank;
    read_request * split_read_request_list;
    int * aggregator_rank_array;
} bp_proc_pvt_struct;

//private structure of read_request
typedef struct read_request_private
{
    int rank;
    void * data;
    int file_idx;
    uint64_t offset;
    void * parent;
} rr_pvt_struct;

static void read_buffer (const ADIOS_FILE * fp,
                         uint64_t buffer_offset,
                         read_request * r,
                         read_request * s
                        );

static int isAggregator (BP_PROC * p)
{
    assert (p);

    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;

    return (pvt->rank == pvt->aggregator_rank);
}

static int rank_to_group_mapping (bp_proc_pvt_struct * pvt, int rank)
{
    int grp_size = pvt->size / pvt->groups;
    int remain = pvt->size - grp_size * pvt->groups;
    int g;

    if (remain == 0)
    {
        g = rank / grp_size;
    }
    else
    {
        if (rank < (grp_size + 1) * remain)
        {
            g = rank / (grp_size + 1);
        }
        else
        {
            g = remain + (rank - (grp_size + 1) * remain) / grp_size;
        }
    }

    return g;
}

static void buffer_write (void ** buffer, void * data, int size)
{
    memcpy (* buffer, data, size);
    * buffer = (char *)* buffer + size;
}

// *****************************************************************************
static void _buffer_write (char ** buffer, uint64_t * buffer_size
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
            fprintf (stderr, "Cannot allocate memory in buffer_write.  "
                             "Requested: %" PRIu64 "\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}

static void _buffer_write32 (char ** buffer, uint32_t * buffer_size
                            ,uint32_t * buffer_offset
                            ,const void * data, uint32_t size
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
            fprintf (stderr, "Cannot allocate memory in buffer_write.  "
                             "Requested: %u\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}

// *****************************************************************************
static void _buffer_read (char * buffer, uint64_t * buffer_offset
                         ,void * data, uint64_t size
                         )
{
    memcpy (data, buffer + *buffer_offset, size);
    *buffer_offset += size;
}

static int calc_data_size (BP_PROC * p)
{
    int size = 0;
    read_request * h = p->local_read_request_list;

    // count
    size += 4;

    while (h)
    {
        // varid + from_steps, nsteps, ndim + start + count
        size += 4 + 4 + 4 + 4 + h->sel->u.bb.ndim * (8 + 8) + 8 + 4 + 8;
        h = h->next;
    }

    return size;
}

/**************************************************
* Get the subfile index and associate data offset *
***************************************************/
static void get_data_addr (const ADIOS_FILE * fp, int varid,
                           const read_request * r,
                           int * file_idx,
                           uint64_t * offset,
                           uint64_t * payload_size
                          )
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    struct adios_index_var_struct_v1 * v;
    int j, idx, t;
    int start_idx, stop_idx;
    int ndim, time, file_is_fortran;
    uint64_t * dims = 0;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t * start, * count;
    int nsteps, is_global = 0, flag;

    file_is_fortran = is_fortran_file (fh);
    v = bp_find_var_byid (fh, varid);

    /* Get dimensions and flip if caller != writer language */
    bp_get_and_swap_dimensions (fp, v, file_is_fortran,
                                &ndim, &dims, &nsteps,
                                file_is_fortran);

    assert (ndim == r->sel->u.bb.ndim);

    start = r->sel->u.bb.start;
    count = r->sel->u.bb.count;

    /* Note fp->current_step is always 0 for file mode. */
    for (t = fp->current_step + r->from_steps;
         t < fp->current_step + r->from_steps + r->nsteps; t++
        )
    {

        if (!p->streaming)
        {
            time = get_time (v, t);
        }
        else
        {
            time = t + 1;
        }

        start_idx = get_var_start_index (v, time);
        stop_idx = get_var_stop_index (v, time);

        if (start_idx < 0 || stop_idx < 0)
        {
            adios_error (err_no_data_at_timestep,
                         "Variable (id=%d) has no data at %d time step in %s\n",
                         varid, t, __FUNCTION__);
            continue;
        }

        if (ndim == 0)
        {
            /* THIS IS A SCALAR VARIABLE */
            idx = 0;

            * file_idx = v->characteristics[start_idx + idx].file_index;
            * offset = v->characteristics[start_idx + idx].payload_offset;
            * payload_size = bp_get_type_size (v->type,
                                               v->characteristics[start_idx + idx].value);
            return;
        }

         /* READ AN ARRAY VARIABLE */
        int * idx_table = (int *) malloc (sizeof(int) * (stop_idx - start_idx + 1));

        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++)
        {
            idx_table[idx] = 1;
            /* Each pg can have a different sized array, so we need the actual dimensions from it */
            is_global = bp_get_dimension_characteristics_notime (&(v->characteristics[start_idx + idx]),
                                                                 ldims,
                                                                 gdims,
                                                                 offsets,
                                                                 file_is_fortran
                                                                );
            if (!is_global)
            {
                memcpy (gdims, ldims, ndim * 8);
            }

            for (j = 0; j < ndim; j++)
            {
                if ( (count[j] > gdims[j])
                  || (start[j] > gdims[j])
                  || (start[j] + count[j] > gdims[j]))
                {
                    adios_error (err_out_of_bound,
                                 "Error: Variable (id=%d) out of bound ("
                                 "the data in dimension %d to read is %" PRIu64 " elements from index %" PRIu64
                                 " but the actual data is [0,%" PRId64 "])",
                                 varid, j + 1, count[j], start[j], gdims[j] - 1);
                    return;
                }

                /* check if there is any data in this pg and this dimension to read in */
                flag = (offsets[j] >= start[j]
                        && offsets[j] < start[j] + count[j])
                    || (offsets[j] < start[j]
                        && offsets[j] + ldims[j] > start[j] + count[j])
                    || (offsets[j] + ldims[j] > start[j]
                        && offsets[j] + ldims[j] <= start[j] + count[j]);

                idx_table[idx] = idx_table[idx] && flag;
            }

            //FIXME
            if (idx_table[idx])
            {
                free (idx_table);
                if (dims)
                {
                    free (dims);
                }

                * file_idx = v->characteristics[start_idx + idx].file_index;
                * offset = v->characteristics[start_idx + idx].payload_offset;
                * payload_size = bp_get_type_size (v->type, v->characteristics[start_idx + idx].value);
                for (j = 0; j < ndim; j++)
                {
                    * payload_size *= ldims[j];
                }
                return;
            }
        }

        free (idx_table);
/*
        if (isTimeless (file_tdim))
            break;
*/
    } // end for (timestep ... loop over timesteps

    if (dims)
    {
        free (dims);
    }
}

static void sort_read_requests (BP_PROC * p)
{
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;
    rr_pvt_struct * rr_r, * rr_t; 
    int file_idx;
    uint64_t offset;
    read_request * r = pvt->split_read_request_list;
    read_request * n = 0, * t, * t_prev, * next;
    while (r)
    {
//printf ("[%d]: r->ra->offset = %" PRIu64 "\n", p->rank, r->ra->offset);
        rr_r = (rr_pvt_struct *) r->priv;
        t = n;
        t_prev = 0;
        next = r->next;

        file_idx = rr_r->file_idx;
        offset = rr_r->offset;

        //while (t && t->ra->file_idx <= file_idx)
        while (t)
        {
            rr_t = (rr_pvt_struct *) t->priv;
            if (rr_t->file_idx > file_idx ||
                (rr_t->file_idx == file_idx && rr_t->offset > offset)
               )
            {
                break;
            }

            t_prev = t;
            t = t->next;
        }

        if (!t_prev)
        {
            r->next = n;
            n = r;
        }
        else
        {
            t_prev->next = r;
            r->next = t;
        }
//if (p->rank == 1) list_print_readers (p, n);
        r = next;
    }

    pvt->split_read_request_list = n;
}

/*
static void send_read_data1 (BP_PROC * p)
{
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;
    read_request * r = p->local_read_request_list;
    rr_pvt_struct * rr_pvt;

    while (r)
    {
        rr_pvt = (rr_pvt_struct *) r->priv;
        assert (rr_pvt);

        if (pvt->rank != rr_pvt->rank)
        {
            MPI_Send (r->data, r->datasize, MPI_BYTE
                     ,rr_pvt->rank - pvt->aggregator_rank, TAG_DATA, pvt->new_comm
                     );
        }

        r = r->next;
    }
}
*/

static void send_read_data (BP_PROC * p)
{
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;
    int g, i;
    char * b = 0, * recv_buff = 0;
    uint32_t offset = 0, buffer_size = 0;
    int size, * sizes = 0, * offsets = 0;
    read_request * r = p->local_read_request_list;
    rr_pvt_struct * rr_pvt;

    MPI_Comm_size (pvt->new_comm, &size);

    sizes = (int *) malloc (size * 4);
    offsets = (int *) malloc (size * 4);
    assert (sizes);
    assert (offsets);

    for (i = 0; i < size; i++)
    {
        sizes[i] = 0;
        offsets[i] = -1;
    }

    while (r)    {
        rr_pvt = (rr_pvt_struct *) r->priv;
        assert (rr_pvt);

        g = rank_to_group_mapping (pvt, rr_pvt->rank);

        /*  g == p->group       -> requests are from processors that belong to th
e group
            p->rank != r->rank  -> requests are NOT from aggregator itself
        */

        if (g == pvt->group && pvt->rank != rr_pvt->rank)
        {
            assert (r->data);

            if (offsets[rr_pvt->rank - pvt->aggregator_rank] == -1)
            {
                offsets[rr_pvt->rank - pvt->aggregator_rank] = offset;
            }

            sizes[rr_pvt->rank - pvt->aggregator_rank] += r->datasize;

            _buffer_write32 (&b, &buffer_size, &offset, r->data, r->datasize);

            /* Free it immediately to avoid double buffering */
            free (r->data);
            r->data = 0;
        }

        r = r->next;
    }

    MPI_Scatterv (b, sizes, offsets, MPI_BYTE
                 ,recv_buff, 0, MPI_BYTE
                 ,0, pvt->new_comm
                 );

    free (b);
    free (sizes);
    free (offsets);
}

/*
static void get_read_data1 (BP_PROC * p)
{
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;
    MPI_Status status;
    read_request * r = p->local_read_request_list;

    if (pvt->rank != pvt->aggregator_rank)
    {
        while (r)
        {
            MPI_Recv (r->data, r->datasize, MPI_BYTE
                     ,MPI_ANY_SOURCE, MPI_ANY_TAG, pvt->new_comm
                     ,&status
                     );

            r = r->next;
        }
    }
}
*/

/* Receive read data from aggregator */
static void get_read_data (BP_PROC * p)
{
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;
    char * b = 0, * recv_buff = 0;
    int * sizes = 0, * offsets = 0;
    int size = 0;
    read_request * r = p->local_read_request_list;

    if (pvt->rank == pvt->aggregator_rank)
    {
        return;
    }

    r = p->local_read_request_list;
    while (r)
    {
        size += r->datasize;
        r = r->next;
    }

    recv_buff = malloc (size);
    if (recv_buff == 0)
    {
        printf ("Warning: the size of data is 0\n");
        return;
    }

    MPI_Scatterv (b, sizes, offsets, MPI_BYTE
                 ,recv_buff, size, MPI_BYTE
                 ,0, pvt->new_comm
                 );

    b = recv_buff;
    r = p->local_read_request_list;
    while (r)
    {
        memcpy (r->data, b, r->datasize);
        b += r->datasize;

        r = r->next;
    }

    free (recv_buff);
}

static void parse_buffer (BP_PROC * p, char * b, int src)
{
    int i, count;
    read_request * r;
    rr_pvt_struct * rr_pvt;

    // count
    count = * (uint32_t *) b;
    b += 4;

    for (i = 0; i < count; i++)
    {
        r = (read_request *) malloc (sizeof (read_request));
        assert (r);

        rr_pvt = (rr_pvt_struct *) malloc (sizeof (rr_pvt_struct));
        assert (rr_pvt);
        r->priv = rr_pvt;

        rr_pvt->rank = src;

        r->varid = * (uint32_t *) b;
        b += 4;

        r->from_steps = * (int32_t *) b;
        b += 4;

        r->nsteps = * (int32_t *) b;
        b += 4;

        r->sel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
        assert (r->sel);

        //FIXME: bb only for now
        r->sel->type = ADIOS_SELECTION_BOUNDINGBOX;
        r->sel->u.bb.ndim = * (uint32_t *) b;
        b += 4;

        r->sel->u.bb.start = (uint64_t *) malloc (r->sel->u.bb.ndim * 8);
        r->sel->u.bb.count = (uint64_t *) malloc (r->sel->u.bb.ndim * 8);
        assert (r->sel->u.bb.start);
        assert (r->sel->u.bb.count);

        memcpy (r->sel->u.bb.start, b, r->sel->u.bb.ndim * 8);
        b += r->sel->u.bb.ndim * 8;

        memcpy (r->sel->u.bb.count, b, r->sel->u.bb.ndim * 8);
        b += r->sel->u.bb.ndim * 8;

        r->datasize = * (uint64_t *) b;
        b += 8;
#if 0
        rr_pvt->file_idx = * (uint32_t *) b;
        b += 4;

        rr_pvt->offset = * (uint64_t *) b;
        b += 8;
#endif
        r->data = malloc (r->datasize);
        assert (r->data);

        // I am already a parent.
        rr_pvt->parent = 0;

        r->next = 0;

        list_insert_read_request_tail (&p->local_read_request_list, r);
    }
}

static read_request * split_read_requests (const ADIOS_FILE * fp, read_request * r)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    struct adios_index_var_struct_v1 * v;
    int i, j, idx, t, varid, nsteps, time, dummy;
    int start_idx, stop_idx;
    int ndim, file_is_fortran;
    uint64_t * dims = 0;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t * start, * count;
    int is_global = 0, flag;
    read_request * h = 0;
    rr_pvt_struct * rr_pvt = (rr_pvt_struct *) r->priv;
    assert (rr_pvt);

    varid = r->varid;
    start = r->sel->u.bb.start;
    count = r->sel->u.bb.count;
    file_is_fortran = is_fortran_file (fh);
    v = bp_find_var_byid (fh, varid);

    bp_get_and_swap_dimensions (fp, v, file_is_fortran,
                                &ndim, &dims, &nsteps,
                                file_is_fortran);

    if (futils_is_called_from_fortran ())
    {
        swap_order (ndim, start, &dummy);
        swap_order (ndim, count, &dummy);
    }
    for (t = fp->current_step + r->from_steps;
         t < fp->current_step + r->from_steps + r->nsteps; t++
        )
    {
        if (!p->streaming)
        {
            time = get_time (v, t);
        }
        else
        {
            time = t + 1;
        }

        start_idx = get_var_start_index (v, time);
        stop_idx = get_var_stop_index (v, time);

        if (start_idx < 0 || stop_idx < 0)
        {
            log_info ("Variable (id=%d) has no data at %d time step in %s\n",
                      varid, t, __FUNCTION__);
            continue;
        }

        if (ndim == 0)
        {
            /* THIS IS A SCALAR VARIABLE */
            idx = 0;

            read_request * n = (read_request *) malloc (sizeof (read_request));
            assert (n);

            rr_pvt_struct * n_pvt = (rr_pvt_struct *) malloc (sizeof (rr_pvt_struct));
            assert (n_pvt);
            n->priv = n_pvt;
           
            n->varid = r->varid;
            n->from_steps = r->from_steps;
            n->nsteps = r->nsteps;
            n_pvt->rank = rr_pvt->rank;
            n_pvt->file_idx = v->characteristics[start_idx + idx].file_index;
            n_pvt->offset = v->characteristics[start_idx + idx].payload_offset;
            n_pvt->parent = r;
            n->sel->u.bb.ndim = 0;
            n->sel->u.bb.start = 0;
            n->sel->u.bb.count = 0;
            n->next = 0;

            list_insert_read_request_tail (&h, n);

            continue;
        }

         /* READ AN ARRAY VARIABLE */
        int * idx_table = (int *) malloc (sizeof(int) * (stop_idx - start_idx + 1));
        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++)
        {
            idx_table[idx] = 1;
            is_global = bp_get_dimension_characteristics_notime (
                            &(v->characteristics[start_idx + idx]),
                            ldims,
                            gdims,
                            offsets,
                            file_is_fortran
                           );
            if (!is_global)
            {
                memcpy (gdims, ldims, ndim * 8);
            }
/*
            printf("ldims   = "); for (j = 0; j<ndim; j++) printf("%d ",ldims[j]); printf("\n");
            printf("gdims   = "); for (j = 0; j<ndim; j++) printf("%d ",gdims[j]); printf("\n");
            printf("offsets = "); for (j = 0; j<ndim; j++) printf("%d ",offsets[j]); printf("\n");
*/
            for (j = 0; j < ndim; j++)
            {
                if ( (count[j] > gdims[j])
                  || (start[j] > gdims[j])
                  || (start[j] + count[j] > gdims[j]))
                {
                    fprintf (stderr, "Error: Variable (id=%d, %s) out of bound ("
                        "the data in dimension %d to read is %" PRIu64 " elements from index %" PRIu64
                        " but the actual data is [0,%" PRId64 "])\n",
                        varid, v->var_name, j + 1, count[j], start[j], gdims[j] - 1);
                    return 0;
                }

                /* check if there is any data in this pg and this dimension to read in */
                flag = (offsets[j] >= start[j]
                        && offsets[j] < start[j] + count[j])
                    || (offsets[j] < start[j]
                        && offsets[j] + ldims[j] > start[j] + count[j])
                    || (offsets[j] + ldims[j] > start[j]
                        && offsets[j] + ldims[j] <= start[j] + count[j]);

                idx_table[idx] = idx_table[idx] && flag;
            }

            if (!idx_table[idx])
            {
                continue;
            }

            /* determined how many (fastest changing) dimensions can we read in in one read */
            int hole_break;
            for (i = ndim - 1; i > -1; i--)
            {
                if (offsets[i] == start[i] && ldims[i] == count[i])
                {
                }
                else
                {
                    break;
                }
            }

            hole_break = i;

            read_request * n = (read_request *) malloc (sizeof (read_request));
            assert (n);

            n->varid = r->varid;
            n->from_steps = r->from_steps;
            n->nsteps = r->nsteps;
            n->sel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
            assert (n->sel);
            //FIXME
            n->sel->type = ADIOS_SELECTION_BOUNDINGBOX;
            n->sel->u.bb.ndim = r->sel->u.bb.ndim;

            n->sel->u.bb.start = (uint64_t *) malloc (ndim * 8);
            assert (n->sel->u.bb.start);
            n->sel->u.bb.count = (uint64_t *) malloc (ndim * 8);
            assert (n->sel->u.bb.count);

            memcpy (n->sel->u.bb.start, start, ndim * 8);
            memcpy (n->sel->u.bb.count, count, ndim * 8);

            n->priv = (rr_pvt_struct *) malloc (sizeof (rr_pvt_struct));
            assert (n->priv);
            rr_pvt_struct * nrr_pvt = (rr_pvt_struct *) n->priv;
            
            nrr_pvt->rank = ((rr_pvt_struct *) r->priv)->rank;;
            nrr_pvt->file_idx = v->characteristics[start_idx + idx].file_index;
            nrr_pvt->offset = v->characteristics[start_idx + idx].payload_offset;
            nrr_pvt->parent = r;
            n->next = 0;

            if (hole_break == -1)
            {
            }
            else if (hole_break == 0)
            {
                int isize;

                isize = offsets[0] + ldims[0];
                if (start[0] >= offsets[0])
                {
                    // head is in
                    if (start[0] < isize)
                    {
                        if (start[0] + count[0] > isize)
                        {
                            n->sel->u.bb.count[0] = isize - start[0];
                        }
                        else
                        {
                            n->sel->u.bb.count[0] = count[0];
                        }
                        n->sel->u.bb.start[0] = start[0];
                    }
                }
                else
                {
                    // middle is in
                    if (isize < start[0] + count[0])
                    {
                        n->sel->u.bb.count[0] = ldims[0];
                    }
                    else
                    {
                    // tail is in
                        n->sel->u.bb.count[0] = count[0] + start[0] - offsets[0];
                    }
                    n->sel->u.bb.start[0] = offsets[0];
                }

            }
            else
            {
                int isize;
                uint64_t size_in_dset[10];
                uint64_t offset_in_dset[10];
                uint64_t offset_in_var[10];

                memset(size_in_dset, 0 , 10 * 8);
                memset(offset_in_dset, 0 , 10 * 8);
                memset(offset_in_var, 0 , 10 * 8);

                for (i = 0; i < ndim; i++)
                {
                    isize = offsets[i] + ldims[i];
                    if (start[i] >= offsets[i])
                    {
                        // head is in
                        if (start[i] < isize)
                        {
                            if (start[i] + count[i] > isize)
                            {
                                n->sel->u.bb.count[i] = isize - start[i];
                            }
                            else
                            {
                                n->sel->u.bb.count[i] = count[i];
                            }
                            n->sel->u.bb.start[i] = start[i];
                            offset_in_var[i] = 0;
                        }
                    }
                    else
                    {
                        // middle is in
                        if (isize < start[i] + count[i])
                        {
                            n->sel->u.bb.count[i] = ldims[i];
                        }
                        else
                        {
                            // tail is in
                            n->sel->u.bb.count[i] = count[i] + start[i] - offsets[i];
                        }
                        n->sel->u.bb.start[i] = offsets[i];
                        //offset_in_dset[i] = 0;
                        offset_in_var[i] = offsets[i] - start[i];
                    }
                }

            }

            n->datasize = bp_get_type_size (v->type, v->characteristics[start_idx + idx].value);
            for (i = 0; i < ndim; i++)
            {
                n->datasize *= n->sel->u.bb.count[i];
            }

            list_insert_read_request_tail (&h, n);
        }

        free (idx_table);
    } // end for (timestep ... loop over timesteps

    if (dims)
    {
        free (dims);
    }

    if (!h)
    {
        fprintf (stderr, "v = %s\n", v->var_name);
    }
    assert (h);

    return h;
}

static void process_read_requests (const ADIOS_FILE * fp)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;
    
    read_request * h = p->local_read_request_list, * n;

    while (h)
    {
        n = split_read_requests (fp, h);

        list_append_read_request_list (&pvt->split_read_request_list, n);

        h = h->next;
    }
}

static void read_chunk (const ADIOS_FILE * fp,
                        int file_idx, 
                        uint64_t chunk_offset, 
                        uint64_t size
                       )
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    MPI_File * sfh;
    MPI_Status status;
    int has_subfile;

    //struct timeval t0, t1;
    //gettimeofday (&t0, NULL);

    has_subfile = has_subfiles (fh);

    bp_realloc_aligned(fh->b, size);
    fh->b->offset = 0;

    if (has_subfile)
    {
        sfh = get_BP_subfile_handle (fh, file_idx);   

        if (!sfh)
        {
            int err;
            char * ch, * name_no_path, * name;
            struct BP_file_handle * new_h =
                  (struct BP_file_handle *) malloc (sizeof (struct BP_file_handle));

            new_h->file_index = file_idx;
            new_h->next = 0;
            if ( (ch = strrchr (fh->fname, '/')))
            {
                name_no_path = malloc (strlen (ch + 1) + 1);
                strcpy (name_no_path, ch + 1);
            }
            else
            {
                name_no_path = malloc (strlen (fh->fname) + 1);
                strcpy (name_no_path, fh->fname);
            }

            name = malloc (strlen (fh->fname) + 5 + strlen (name_no_path) + 1 + 10 + 1);
            sprintf (name, "%s.dir/%s.%d", fh->fname, name_no_path, new_h->file_index);

            err = MPI_File_open (MPI_COMM_SELF
                                ,name
                                ,MPI_MODE_RDONLY
                                ,MPI_INFO_NULL
                                ,&new_h->fh
                                );
            if (err != MPI_SUCCESS)
            {
                adios_error (err_file_open_error, "Can not open file %s\n", name);
                return;
            }

            add_BP_subfile_handle (fh, new_h);
            sfh = &new_h->fh;

            free (name_no_path);
            free (name);
        }

        MPI_File_seek (*sfh
                      ,(MPI_Offset)chunk_offset
                      ,MPI_SEEK_SET
                      );
        MPI_File_read (*sfh
                      ,fh->b->buff
                      ,size
                      ,MPI_BYTE
                      ,&status
                      );
    }
    else
    {
        MPI_File_seek (fh->mpi_fh
                      ,(MPI_Offset)chunk_offset
                      ,MPI_SEEK_SET
                      );
        MPI_File_read (fh->mpi_fh
                      ,fh->b->buff
                      ,size
                      ,MPI_BYTE
                      ,&status
                      );
    }

    fh->b->offset = 0;

    //gettimeofday (&t1, NULL);
    //printf ("read chunk time = %f \n", t1.tv_sec - t0.tv_sec + (double)(t1.tv_usec - t0.tv_usec)/1000000 );
}

/* 
The read request link list needs to be sorted beforehand
*/
static void do_read (const ADIOS_FILE * fp)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;
    int file_idx;
    uint64_t offset, payload_size;
    //struct timeval t0;
    //struct timeval t1;
    //double t2, t3, t4, t5;

    read_request * s = pvt->split_read_request_list, * f_start = s, * f_end = s;
    read_request * o_start = s, * o_end = s, * o_prev_end = 0, * parent = 0;

    //gettimeofday (&t0, NULL);

    //t2 = MPI_Wtime();

    while (f_start)
    {
        f_end = f_start;

        // Find a set of reqeusts that fall into the same file, i.e., [f_start, f_end)
        while (f_end
          && ((rr_pvt_struct *) f_end->priv)->file_idx == ((rr_pvt_struct *) f_start->priv)->file_idx)
        {
            f_end = f_end->next;
        }

        o_start = f_start;
        o_end = f_start;
        o_prev_end = 0;

        while (o_start != f_end)
        {
            // Find a set of requests that fall into chunk size, i.e., [o_start, o_end)
            while (o_end && o_end != f_end
                  && ((rr_pvt_struct *) o_end->priv)->offset - ((rr_pvt_struct *)o_start->priv)->offset <= chunk_buffer_size)
            {
                o_prev_end = o_end;
                o_end = o_end->next;
            }

            // Calculate the var payload size of the last request
            get_data_addr (fp,
                           o_prev_end->varid, 
                           o_prev_end, 
                           &file_idx, 
                           &offset, 
                           &payload_size
                          );
//printf ("o_start.offset = %" PRIu64 "\n", o_start->ra->offset);
//printf ("o_prev_end.offset = %" PRIu64 "\n", o_prev_end->ra->offset);

            //t4 = MPI_Wtime ();
            // read a chunk from file into internal buffer
            read_chunk (fp, 
                        ((rr_pvt_struct *) o_start->priv)->file_idx, 
                        ((rr_pvt_struct *) o_start->priv)->offset, 
                        ((rr_pvt_struct *) o_prev_end->priv)->offset
                         - ((rr_pvt_struct *) o_start->priv)->offset + payload_size
                       );

            //t5 = MPI_Wtime ();
//    printf ("read chunk = %f \n", t5 - t4);

            s = o_start;
            do
            {
                parent = ((rr_pvt_struct *) s->priv)->parent;
                // copy data from internal buffer to user buffer
                read_buffer (fp, ((rr_pvt_struct *) o_start->priv)->offset, parent, s);

                s = s->next;
            } while (s != o_end);

            o_start = o_end;
            o_prev_end = 0;
        }

        f_start = f_end;
    }

    //gettimeofday (&t1, NULL);
    //t3 = MPI_Wtime ();

    //printf ("while time = %f \n", t1.tv_sec - t0.tv_sec + (double)(t1.tv_usec - t0.tv_usec)/1000000 );
    //printf ("while time = %f \n", t3 - t2);
}

/* r - the original read rquest
 * s - the read request as the result of spliting
 */
static void read_buffer (const ADIOS_FILE * fp,
                         uint64_t buffer_offset,
                         read_request * r,
                         read_request * s
                        )
{
#define MAX_DIMS 32
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    struct adios_index_var_struct_v1 * v;
    uint64_t * s_start, * s_count;
    uint64_t * start, * count;
    int i, j, idx, t, time, ndim, dummy;
    int varid, start_idx, stop_idx, file_is_fortran;
    uint64_t ldims[MAX_DIMS], gdims[MAX_DIMS], offsets[MAX_DIMS];
    uint64_t datasize, dset_stride, var_stride, total_size = 0, items_read;
    int is_global = 0, size_unit, break_dim, idx_check1, idx_check2;
    uint64_t slice_offset, slice_size;
    void * data;

    varid = r->varid;
    v = bp_find_var_byid (fh, varid);
    data = r->data;
    start = r->sel->u.bb.start;
    count = r->sel->u.bb.count;
    s_start = s->sel->u.bb.start;
    s_count = s->sel->u.bb.count;
    file_is_fortran = is_fortran_file (fh);
    ndim = r->sel->u.bb.ndim;
    size_unit = bp_get_type_size (v->type,
                                  v->characteristics [0].value
                                 );

    /* We need to swap them here to read correctly in C order */
    if (futils_is_called_from_fortran ())
    {
        swap_order (ndim, start, &dummy);
        swap_order (ndim, count, &dummy);
    }

    for (i = 0, items_read = 1; i < ndim; i++)
    {
        items_read *= count[i];
    }

    /* For each timestep, do the reading separately.
     * Note fp->current_step is always 0 for file mode.
     */
    for (t = fp->current_step + r->from_steps;
         t < fp->current_step + r->from_steps + r->nsteps; t++
        )
    {
        if (!p->streaming)
        {
            time = get_time (v, t);
        }
        else
        {
            time = t + 1;
        }

        start_idx = get_var_start_index (v, time);
        stop_idx = get_var_stop_index (v, time);

        if (start_idx < 0 || stop_idx < 0)
        {
            log_info ("Variable (id=%d) has no data at %d time step in %s\n",
                      varid, t, __FUNCTION__
                     );
            continue;
        }

        if (ndim == 0)
        {
            /* READ A SCALAR VARIABLE */
            slice_size = 1 * size_unit;
            idx = 0;

            if (v->type == adios_string)
            {
                // strings are stored without \0 in file
                // size_of_type here includes \0 so decrease by one
                size_unit--;
            }

            slice_offset = v->characteristics[start_idx + idx].payload_offset;

            memcpy (data, fh->b->buff + slice_offset - buffer_offset, size_unit);
            if (fh->mfooter.change_endianness == adios_flag_yes)
            {
                change_endianness (data, size_unit, v->type);
            }

            if (v->type == adios_string)
            {
                // add \0 to the end of string
                // size_of_type here is the length of string
                // FIXME: how would this work for strings written over time?
                ((char *)data + total_size)[size_unit] = '\0';
            }

            total_size += size_unit;

            continue;
        }

         /* READ AN ARRAY VARIABLE */
        // loop over the list of pgs to read from one-by-one
        for (idx = 0; idx < stop_idx - start_idx + 1; idx++)
        {
            int flag1, flag2;
            uint64_t payload_size = size_unit;

            datasize = 1;
            var_stride = 1;
            dset_stride = 1;
            idx_check1 = 1;
            idx_check2 = 1;
            is_global = bp_get_dimension_characteristics_notime (
                                    &(v->characteristics[start_idx + idx]),
                                    ldims,
                                    gdims,
                                    offsets,
                                    file_is_fortran
                                   );

            if (!is_global)
            {
                memcpy (gdims, ldims, ndim * 8);
            }

            for (j = 0; j < ndim; j++)
            {
                payload_size *= ldims [j];

                if ( (count[j] > gdims[j])
                  || (start[j] > gdims[j])
                  || (start[j] + count[j] > gdims[j]))
                {
                    log_error ("Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %" PRIu64 " elements from index %" PRIu64
                        " but the actual data is [0,%" PRId64 "])",
                        varid, j+1, count[j], start[j], gdims[j] - 1);
                    return;
                }

                /* check if there is any data in this pg and this dimension to read in */
                flag1 = (offsets[j] >= start[j]
                        && offsets[j] < start[j] + count[j])
                    || (offsets[j] < start[j]
                        && offsets[j] + ldims[j] > start[j] + count[j])
                    || (offsets[j] + ldims[j] > start[j]
                        && offsets[j] + ldims[j] <= start[j] + count[j]);

                idx_check1 = idx_check1 && flag1;

                flag2 = (offsets[j] >= s_start[j]
                        && offsets[j] < s_start[j] + s_count[j])
                    || (offsets[j] < s_start[j]
                        && offsets[j] + ldims[j] > s_start[j] + s_count[j])
                    || (offsets[j] + ldims[j] > s_start[j]
                        && offsets[j] + ldims[j] <= s_start[j] + s_count[j]);

                idx_check2 = idx_check2 && flag2;
            }

            if (!idx_check1)
            {
                continue;
            }

            break_dim =  ndim - 1;
            while (break_dim > -1)
            {
                if (start[break_dim] == 0 &&
                    ldims[break_dim] == count[break_dim]
                   )
                {
                    datasize *= ldims[break_dim];
                }
                else
                    break;

                break_dim--;
            }

            slice_offset = 0;
            slice_size = 0;

            if (break_dim == -1)
            {
                slice_size = payload_size;

                slice_offset = v->characteristics[start_idx + idx].payload_offset;

                if (idx_check2)
                {
                    memcpy (data, fh->b->buff + slice_offset - buffer_offset, slice_size);
                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {
                        change_endianness (data, slice_size, v->type);
                    }
                }
            }
            else if (break_dim == 0)
            {
                int isize;
                uint64_t size_in_dset = 0, offset_in_dset = 0;
                uint64_t offset_in_var = 0, write_offset;

                isize = offsets[0] + ldims[0];
                if (start[0] >= offsets[0])
                {
                    // head is in
                    if (start[0] < isize)
                    {
                        if (start[0] + count[0] > isize)
                        {
                            size_in_dset = isize - start[0];
                        }
                        else
                        {
                            size_in_dset = count[0];
                        }
                        offset_in_var = 0;
                        offset_in_dset = start[0] - offsets[0];
                    }
                }
                else
                {
                    // middle is in
                    if (isize < start[0] + count[0])
                    {
                        size_in_dset = ldims[0];
                    }
                    else
                    {
                    // tail is in
                        size_in_dset = count[0] + start[0] - offsets[0];
                    }
                    offset_in_var = offsets[0] - start[0];
                    offset_in_dset = 0;
                }

                slice_size = size_in_dset * datasize * size_unit;
                slice_offset = v->characteristics[start_idx + idx].payload_offset
                                 + offset_in_dset * datasize * size_unit;

                write_offset = offset_in_var * size_unit;
                for (i = 1; i < ndim; i++)
                {
                    write_offset *= count[i];
                }

                if (idx_check2)
                {
                    memcpy ((char *) data + write_offset,
                            fh->b->buff + slice_offset - buffer_offset,
                            slice_size
                           );
                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {
                        change_endianness ((char *) data + write_offset, slice_size, v->type);
                    }
                }
            }
            else
            {
                int isize;
                uint64_t size_in_dset[MAX_DIMS];
                uint64_t offset_in_dset[MAX_DIMS];
                uint64_t offset_in_var[MAX_DIMS];

                memset (size_in_dset, 0, MAX_DIMS * 8);
                memset (offset_in_dset, 0, MAX_DIMS * 8);
                memset (offset_in_var, 0, MAX_DIMS * 8);

                for (i = 0; i < ndim; i++)
                {
                    isize = offsets[i] + ldims[i];
                    if (start[i] >= offsets[i])
                    {
                        // head is in
                        if (start[i] < isize)
                        {
                            if (start[i] + count[i] > isize)
                            {
                                size_in_dset[i] = isize - start[i];
                            }
                            else
                            {
                                size_in_dset[i] = count[i];
                            }
                            offset_in_dset[i] = start[i] - offsets[i];
                            offset_in_var[i] = 0;
                        }
                        else
                        {
                        }
                    }
                    else
                    {
                        // middle is in
                        if (isize < start[i] + count[i])
                        {
                            size_in_dset[i] = ldims[i];
                        }
                        else
                        {
                            // tail is in
                            size_in_dset[i] = count[i] + start[i] - offsets[i];
                        }
                        offset_in_dset[i] = 0;
                        offset_in_var[i] = offsets[i] - start[i];
                    }
                }

                datasize = 1;
                var_stride = 1;

                for (i = ndim - 1; i >= break_dim; i--)
                {
                    datasize *= size_in_dset[i];
                    dset_stride *= ldims[i];
                    var_stride *= count[i];
                }

                uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
                uint64_t var_offset = 0, dset_offset = 0;

                for (i = ndim - 1; i > -1; i--)
                {
                    start_in_payload += s * offset_in_dset[i] * size_unit;
                    end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_unit;
                    s *= ldims[i];
                }

                slice_size = end_in_payload - start_in_payload + 1 * size_unit;
                slice_offset =  v->characteristics[start_idx + idx].payload_offset
                                  + start_in_payload;

                for (i = 0; i < ndim; i++)
                {
                    offset_in_dset[i] = 0;
                }

                for (i = 0; i < ndim; i++)
                {
                    var_offset = offset_in_var[i] + var_offset * count[i];
                    dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
                }

                if (idx_check2)
                {
                    copy_data (data
                              ,fh->b->buff + fh->b->offset + slice_offset - buffer_offset
                              ,0
                              ,break_dim
                              ,size_in_dset
                              ,ldims
                              ,count
                              ,var_stride
                              ,dset_stride
                              ,var_offset
                              ,dset_offset
                              ,datasize
                              ,size_unit
                              ,fh->mfooter.change_endianness
                              ,v->type
                              );

                }
            }
        }  // for idx ... loop over pgs

        // shift target pointer for next read in
        data = (char *)data + (items_read * size_unit);
    } // for t
#undef MAX_DIMS
}

static void broadcast_fh_buffer (ADIOS_FILE * fp)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    struct bp_index_pg_struct_v1 * pgs_root = fh->pgs_root, * pg;
    struct adios_index_var_struct_v1 * vars_root = fh->vars_root, * v;
    char * buffer;
    uint64_t buffer_size, buffer_offset = 0;
    int i, j, nsteps;
    uint16_t len;
    uint8_t flag;
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;

    bp_realloc_aligned (fh->b, 0);
/*
    buffer = fh->b->buff;
*/
    buffer = 0;
    buffer_size = 0;
    buffer_offset = 0;
    if (isAggregator (p))
    {
        _buffer_write (&buffer, &buffer_size, &buffer_offset,
                       &pvt->num_aggregators, 4); // n_sf
        _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                       &fh->gvar_h->group_count, 2); //group_count 
        _buffer_write (&buffer, &buffer_size, &buffer_offset,
                       &fh->mfooter.pgs_count, 8); //vars_count 
        _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                       &fh->mfooter.vars_count, 4); //vars_count 
        _buffer_write (&buffer, &buffer_size, &buffer_offset,
                       &fh->mfooter.attrs_count, 4); //attrs_count 

        for (i = 0; i < fh->gvar_h->group_count; i++)
        {
            len = strlen (fh->gvar_h->namelist[i]);
            _buffer_write (&buffer, &buffer_size, &buffer_offset,
                           &len, 2); // namelist
            _buffer_write (&buffer, &buffer_size, &buffer_offset,
                           fh->gvar_h->namelist[i], len); // namelist
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &fh->gvar_h->var_counts_per_group[i], sizeof(uint32_t)); // var_counts_per_group
        }

        _buffer_write (&buffer, &buffer_size, &buffer_offset,
                       &fh->gattr_h->group_count, 2); //group_count 
        for (i = 0; i < fh->gattr_h->group_count; i++)
        {
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &fh->gattr_h->attr_counts_per_group[i], sizeof(uint32_t)); // attr_counts_per_group
        }

        for (i = 0; i < fh->mfooter.vars_count; i++)
        {
            len = strlen (fh->gvar_h->var_namelist[i]);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &len, 2); // namelist
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           fh->gvar_h->var_namelist[i], len); // namelist
        }

        for (i = 0; i < fh->mfooter.attrs_count; i++)
        {
            len = strlen (fh->gattr_h->attr_namelist[i]);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &len, 2); // namelist
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           fh->gattr_h->attr_namelist[i], len); // namelist
        }

        _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                       &fh->tidx_start, 4); // tidx_start
        _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                       &fh->tidx_stop, 4); // tidx_start

        pgs_root = fh->pgs_root;
        while (pgs_root)
        {
            len = strlen (pgs_root->group_name);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &len, 2); // group_name len
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           pgs_root->group_name, len); // group_name

            flag = (pgs_root->adios_host_language_fortran == adios_flag_yes) ? 'y' : 'n';
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &flag, 1);

            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &pgs_root->process_id, 4);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &pgs_root->time_index, 4);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &pgs_root->offset_in_file, 8);

            pgs_root = pgs_root->next;
        }

        vars_root = fh->vars_root;
        while (vars_root)
        {
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &vars_root->id, 4); // id

            len = strlen (vars_root->group_name);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &len, 2); // group_name len
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           vars_root->group_name, len); // group_name

            len = strlen (vars_root->var_name);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &len, 2); // var_name len
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           vars_root->var_name, len); // var_name

            len = strlen (vars_root->var_path);
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &len, 2); // var_path len
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           vars_root->var_path, len); // var_path
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &vars_root->type, 4); // type

            ADIOS_VARINFO * vi = bp_inq_var_byid (fp, vars_root->id - 1);
            assert (vi);

            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &vi->ndim, 4); // ndim
            if (vi->ndim)
            {
                _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                               vi->dims, vi->ndim * 8); // dims
            }

            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &vi->nsteps, 4); // ndim
            len = (vars_root->characteristics[0].value == 0) ?
                  0 : bp_get_type_size (vars_root->type, vars_root->characteristics[0].value);
            if (vars_root->type == adios_string)
            {
                len--;
            }
/*
int rank;
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
if (rank == 0)
fprintf (stderr, "bc %s bo 1 = %" PRIu64 ", bo 2 = %" PRIu64 ", len = %d\n", vars_root->var_name, bo, buffer_offset, len);
*/
            _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                           &len, 2);

            if (len)
            {
                _buffer_write (&buffer, &buffer_size, &buffer_offset, 
                               vars_root->characteristics[0].value, len); // ndim

            }

            common_read_free_varinfo (vi);

            vars_root = vars_root->next;
        }
    }

    MPI_Bcast (&buffer_offset, 8, MPI_BYTE, 0, pvt->new_comm);
    if (!isAggregator (p))
    {
/*
        bp_realloc_aligned (fh->b, buffer_offset);
        assert (fh->b->buff);

        buffer = fh->b->buff;
*/
        buffer = malloc (buffer_offset);
        assert (buffer);
    }

    MPI_Bcast (buffer, buffer_offset, MPI_BYTE, 0, pvt->new_comm);

    if (!isAggregator (p))
    {
        uint16_t len; 

        buffer_offset = 0;

        _buffer_read (buffer, &buffer_offset, &pvt->num_aggregators, 4);

        fh->gvar_h = (struct BP_GROUP_VAR *) malloc (sizeof (struct BP_GROUP_VAR));
        assert (fh->gvar_h);
        fh->gvar_h->time_index = 0;
        fh->gvar_h->var_offsets = 0;
        fh->gvar_h->pg_offsets = 0;

        _buffer_read (buffer, &buffer_offset, &fh->gvar_h->group_count, 2); //group_count 
        _buffer_read (buffer, &buffer_offset, &fh->mfooter.pgs_count, 8); //pgs_count 
        _buffer_read (buffer, &buffer_offset, &fh->mfooter.vars_count, 4); //vars_count 
        _buffer_read (buffer, &buffer_offset, &fh->mfooter.attrs_count, 4); //attrs_count 

        fh->gvar_h->namelist = (char **) malloc (fh->gvar_h->group_count * sizeof (char *));
        fh->gvar_h->var_counts_per_group = (uint32_t *) malloc (fh->gvar_h->group_count * sizeof(uint32_t));

        for (i = 0; i < fh->gvar_h->group_count; i++)
        {
            _buffer_read (buffer, &buffer_offset, &len, 2); // len
            fh->gvar_h->namelist[i] = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, fh->gvar_h->namelist[i], len); // namelist
            fh->gvar_h->namelist[i][len] = '\0';

            _buffer_read (buffer, &buffer_offset, &fh->gvar_h->var_counts_per_group[i], sizeof(uint32_t)); // var_counts_per_group
        }

        fh->gattr_h = (struct BP_GROUP_ATTR *) malloc (sizeof (struct BP_GROUP_ATTR));
        assert (fh->gattr_h);
        fh->gattr_h->attr_offsets = 0;

        _buffer_read (buffer, &buffer_offset, &fh->gattr_h->group_count, 2); //group_count 

        fh->gattr_h->attr_counts_per_group = (uint32_t *) malloc (fh->gattr_h->group_count * sizeof(uint32_t));
        fh->gattr_h->namelist = fh->gvar_h->namelist;

        for (i = 0; i < fh->gattr_h->group_count; i++)
        {
            _buffer_read (buffer, &buffer_offset, &fh->gattr_h->attr_counts_per_group[i], sizeof(uint32_t)); // attr_counts_per_group
        }

        fh->gvar_h->var_namelist = (char **) malloc (fh->mfooter.vars_count * sizeof (char *));
        for (i = 0; i < fh->mfooter.vars_count; i++)
        {
            _buffer_read (buffer, &buffer_offset, &len, 2); // len
            fh->gvar_h->var_namelist[i] = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, fh->gvar_h->var_namelist[i], len); // namelist
            fh->gvar_h->var_namelist[i][len] = '\0';
        }

        fh->gattr_h->attr_namelist = (char **) malloc (fh->mfooter.attrs_count * sizeof (char *));
        for (i = 0; i < fh->mfooter.attrs_count; i++)
        {
            _buffer_read (buffer, &buffer_offset, &len, 2); // len
            fh->gattr_h->attr_namelist[i] = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, fh->gattr_h->attr_namelist[i], len); // namelist
            fh->gattr_h->attr_namelist[i][len] = '\0';
        }

        _buffer_read (buffer, &buffer_offset, &fh->tidx_start, 4);
        _buffer_read (buffer, &buffer_offset, &fh->tidx_stop, 4);

        pgs_root = 0;
        int pgs_count = fh->mfooter.pgs_count;

        for (i = 0; i < pgs_count; i++)
        {
            pg = (struct bp_index_pg_struct_v1 *) malloc (sizeof (struct bp_index_pg_struct_v1));
            assert (pg);

            _buffer_read (buffer, &buffer_offset, &len, 2);

            pg->group_name = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, pg->group_name, len);
            pg->group_name[len] = '\0';

            _buffer_read (buffer, &buffer_offset, &flag, 1);
            pg->adios_host_language_fortran = (flag == 'y') ? adios_flag_yes : adios_flag_no;

            _buffer_read (buffer, &buffer_offset, &pg->process_id, 4);
            _buffer_read (buffer, &buffer_offset, &pg->time_index, 4);
            _buffer_read (buffer, &buffer_offset, &pg->offset_in_file, 8);

            pg->time_index_name = 0;
            pg->next = 0;

            if (!pgs_root)
            {
                pgs_root = pg;
                fh->pgs_root = pg;
            }
            else
            {
                pgs_root->next = pg;
                pgs_root = pg;
            }
        }

        fh->vars_table = (struct adios_index_var_struct_v1 **) malloc (fh->mfooter.vars_count * 8);

        vars_root = 0;
        for (i = 0; i < fh->mfooter.vars_count; i++)
        {
            v = (struct adios_index_var_struct_v1 *) malloc (sizeof (struct adios_index_var_struct_v1));
            assert (v);
            fh->vars_table[i] = v;

            _buffer_read (buffer, &buffer_offset, &v->id, 4);

            _buffer_read (buffer, &buffer_offset, &len, 2);
            v->group_name = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, v->group_name, len);
            v->group_name[len] = '\0';

            _buffer_read (buffer, &buffer_offset, &len, 2);
            v->var_name = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, v->var_name, len);
            v->var_name[len] = '\0';

            _buffer_read (buffer, &buffer_offset, &len, 2);
            v->var_path = (char *) malloc (len + 1);
            _buffer_read (buffer, &buffer_offset, v->var_path, len);
            v->var_path[len] = '\0';

            _buffer_read (buffer, &buffer_offset, &v->type, 4);

            v->characteristics_count = 1;
            v->characteristics = (struct adios_index_characteristic_struct_v1 *)
                                     malloc (sizeof (struct adios_index_characteristic_struct_v1));
            assert (v->characteristics);

            v->characteristics->stats = 0;
            _buffer_read (buffer, &buffer_offset, &(v->characteristics->dims.count), 4);

            v->characteristics->dims.dims = 0;
            if (v->characteristics->dims.count)
            {
                v->characteristics->dims.dims = (uint64_t *) malloc (v->characteristics->dims.count * 3 * 8);
                assert (v->characteristics->dims.dims);

                uint64_t * gdims = (uint64_t *) malloc (v->characteristics->dims.count * 8);
                assert (gdims);

                _buffer_read (buffer, &buffer_offset, gdims, v->characteristics->dims.count * 8);
                for (j = 0; j < v->characteristics->dims.count; j++)
                {
                    v->characteristics->dims.dims[j * 3] = gdims[j];
                    v->characteristics->dims.dims[j * 3 + 1] = gdims[j];
                    v->characteristics->dims.dims[j * 3 + 2] = 0;
                }
                free (gdims);
            }
            _buffer_read (buffer, &buffer_offset, &nsteps, 4);
            _buffer_read (buffer, &buffer_offset, &len, 2);

            if (len)
            {
                if (v->type == adios_string)
                {
                    v->characteristics->value = malloc (len + 1);
                }
                else
                {
                    v->characteristics->value = malloc (len);
                }

                _buffer_read (buffer, &buffer_offset, v->characteristics->value, len);

                if (v->type == adios_string)
                {
                    ((char * )(v->characteristics->value))[len] = '\0';
                }
            }
            else
            {
                v->characteristics->value = 0;
            }
            v->next = 0;

            if (!vars_root)
            {
                vars_root = v;
                fh->vars_root = v;
            }
            else
            {
                vars_root->next = v;
                vars_root = v;
            }
        }
    }

    if (buffer)
    {
        free (buffer);
    }
}

int adios_read_bp_staged_init_method (MPI_Comm comm, PairStruct * params)
{
    int global_rank;
    PairStruct * p = params;
    char * env_str;

    while (p)
    {
        if (!strcasecmp (p->name, "max_chunk_size"))
        {
            chunk_buffer_size = strtol(p->value, NULL, 10);
            if (chunk_buffer_size > 0)
            {
                log_debug ("max_chunk_size set to %dMB for the read method\n", chunk_buffer_size);
                chunk_buffer_size = chunk_buffer_size * 1024 * 1024;
            }
        }
        else if (!strcasecmp (p->name, "poll_interval"))
        {
            errno = 0;
            poll_interval = strtol(p->value, NULL, 10);
            if (poll_interval > 0 && !errno)
            {
                log_debug ("poll_interval set to %d secs for READ_BP read method\n",
                            poll_interval);
            }
            else
            {
                log_error ("Invalid 'poll_interval' parameter given to the READ_BP "
                            "read method: '%s'\n", p->value);
            }
        }
        else if (!strcasecmp (p->name, "show_hidden_attrs"))
        {
            show_hidden_attrs = 1;

            log_debug ("show_hidden_attrs is set\n");
        }
        else if (!strcasecmp (p->name, "num_aggregators"))
        {
            errno = 0;
            num_aggregators = strtol (p->value, NULL, 10);
            if (num_aggregators > 0 && !errno)
            {
                log_debug ("num_aggregators set to %d for STAGED_READ_BP read method",
                           num_aggregators);
            }
        }

        p = p->next;
    }

    /* if v1 read api is being used, the aggregation parameters can only
     * be set through environment variables and the follow
     * code deals with this.
     */
    MPI_Comm_rank (MPI_COMM_WORLD, &global_rank);
    if (num_aggregators <= 0)
    {
        env_str = getenv ("num_aggregators");
        if (!env_str)
        {
            adios_error (err_unspecified, 
                         "Environment variable \"num_aggregators\" hasn't been set.\n");
            exit(0);
        }

        num_aggregators = atoi (env_str);

        if (global_rank == 0)
        {
            printf ("%d aggregators are used.\n", num_aggregators);
        }
    }

    if (chunk_buffer_size <= 0)
    {
        env_str = getenv ("chunk_size");
        if (!env_str)
        {
            adios_error (err_unspecified,
                         "Environment variable \"chunk_size\" (in MB) hasn't been set.\n");
            exit(0);
        }

        chunk_buffer_size = 1024 * 1024 * atoi (env_str);
    }

    return 0;
}

int adios_read_bp_staged_finalize_method ()
{
    /* Set these back to default */
    chunk_buffer_size = -1;
    poll_interval = 10;
    show_hidden_attrs = 0;
    num_aggregators = -1;

    return 0;
}

ADIOS_FILE * adios_read_bp_staged_open (const char * fname, MPI_Comm comm, enum ADIOS_LOCKMODE lock_mode, float timeout_sec)
{
    log_error (" adios_read_open() is not supported in this method. Use the file-only adios_read_open_file().\n");
    return 0;
}

static void free_proc_struct (BP_PROC * p)
{
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;

    /* We need to free split_read_request_list cause
     * people might do two back-to-back sets of
     * schedule_read, schedule_read, perform_read
     */
    list_free_read_request (pvt->split_read_request_list);
    pvt->split_read_request_list = 0;

    list_free_read_request (p->local_read_request_list);
    p->local_read_request_list = 0;
}

static void init_read (BP_PROC * p)
{
    BP_FILE * fh = (BP_FILE *) p->fh;
    int i, remain, global_rank;
    int color1, color2;
    char * env_str;

    /* if v1 read api is being used, the aggregation parameters can only
     * be set through environment variables and the follow
     * code deals with this.
     */
    MPI_Comm_rank (MPI_COMM_WORLD, &global_rank);
    if (num_aggregators <= 0)
    {
        env_str = getenv ("num_aggregators");
        if (!env_str)
        {
            adios_error (err_unspecified,
                         "Environment variable \"num_aggregators\" hasn't been set.\n");
            exit(0);
        }

        num_aggregators = atoi (env_str);

        if (global_rank == 0)
        {
            printf ("%d aggregators are used.\n", num_aggregators);
        }
    }

    if (chunk_buffer_size <= 0)
    {
        env_str = getenv ("chunk_size");
        if (!env_str)
        {
            adios_error (err_unspecified,
                         "Environment variable \"chunk_size\" (in MB) hasn't been set.\n");
            exit(0);
        }

        chunk_buffer_size = 1024 * 1024 * atoi (env_str);
    }

    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) malloc (sizeof (bp_proc_pvt_struct));
    assert (pvt);

    p->priv = pvt;
    MPI_Comm_rank (fh->comm, &pvt->rank);
    MPI_Comm_size (fh->comm, &pvt->size);

    pvt->num_aggregators = num_aggregators;
    pvt->groups = (pvt->num_aggregators > pvt->size || pvt->num_aggregators <= 0) ? pvt->size : pvt->num_aggregators;
    pvt->group_size = pvt->size / pvt->groups;
    remain = pvt->size - pvt->group_size * pvt->groups;

    pvt->aggregator_rank_array = (int *) malloc (pvt->groups * sizeof (int));
    for (i = 0; i < pvt->groups; i++)
    {
        if (remain == 0)
        {
            pvt->aggregator_rank_array[i] = pvt->group_size * i;
        }
        else
        {
            if (i < remain)
            {
                pvt->aggregator_rank_array[i] = (pvt->group_size + 1) * i;
            }
            else
            {
                pvt->aggregator_rank_array[i] = remain * (pvt->group_size + 1) + (i - remain) * pvt->group_size;
            }
        }
    }

    if (remain == 0)
    {
        color1 = pvt->rank / pvt->group_size;
        color2 = pvt->rank % pvt->group_size;

        pvt->aggregator_rank = color1 * pvt->group_size;
    }
    else
    {
        if (pvt->rank < (pvt->group_size + 1) * remain)
        {
            color1 = pvt->rank / (pvt->group_size + 1);
            color2 = pvt->rank % (pvt->group_size + 1);

            pvt->aggregator_rank = color1 * (pvt->group_size + 1);
            pvt->group_size++;
        }
        else
        {
            color1 = remain + (pvt->rank - (pvt->group_size + 1) * remain) / pvt->group_size;
            color2 = (pvt->rank - (pvt->group_size + 1) * remain) % pvt->group_size;

            pvt->aggregator_rank = remain * (pvt->group_size + 1) + (color1 - remain) * pvt->group_size;
        }
    }

    pvt->group = color1;

    MPI_Comm_split (fh->comm, color1, pvt->rank, &pvt->new_comm);
    MPI_Comm_split (fh->comm, color2, pvt->rank, &pvt->new_comm2);
    MPI_Comm_rank (pvt->new_comm, &pvt->new_rank);

    pvt->aggregator_new_rank = 0;
    pvt->group_comm = fh->comm;
    pvt->split_read_request_list = 0;
    p->b = 0;

    return;
}

ADIOS_FILE * adios_read_bp_staged_open_file (const char * fname, MPI_Comm comm)
{
    int rank;
    BP_PROC * p;
    BP_FILE * fh;
    ADIOS_FILE * fp;
    bp_proc_pvt_struct * pvt;

    log_debug ("adios_read_bp_staged_open_file\n");

    MPI_Comm_rank (comm, &rank);

    fh = (BP_FILE *) malloc (sizeof (BP_FILE));
    assert (fh);
    fh->fname = (fname ? strdup (fname) : 0L);
    fh->mpi_fh = 0;
    fh->subfile_handles.n_handles = 0;
    fh->comm = comm;
    fh->gvar_h = 0;
    fh->pgs_root = 0;
    fh->vars_root = 0;
    fh->attrs_root = 0;
    fh->vars_table = 0;
    fh->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    assert (fh->b);
    adios_buffer_struct_init (fh->b);

    p = (BP_PROC *) malloc (sizeof (BP_PROC));
    assert (p);
    p->fh = fh;
    p->streaming = 0;
    p->varid_mapping = 0; // maps perceived id to real id
    p->local_read_request_list = 0;
    p->b = 0;
    p->priv = 0;
    init_read (p);

    fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    assert (fp);
    fp->fh = (uint64_t) p;

    pvt = (bp_proc_pvt_struct *) p->priv;
    if (isAggregator (p))
    {
        if (bp_open (fname, pvt->new_comm2, fh) < 0)
        {
            adios_error (err_file_open_error, "File open failed: %s\n", fname);
            return 0;
        }
    }

    broadcast_fh_buffer (fp);

    /* '-1' means that we want all steps. 
     * This will seek to the last step. So we need to set current_step back properly.
     * Usually bp_seek_to_step comes after release_step call, to first free up some
     * memory allocated by the previous step. This is the first seek call and, therefore,
     * no release_step.
     */
    bp_seek_to_step (fp, -1, show_hidden_attrs);
    /* It was agreed that, for file open the current step should be set to 0,
     * instead of the start time. The var_namelist and attr_namel
ist should
     * consist of all steps. For stream open, this is done differ
ently.
     * 07/2012 - Q.Liu
     */
    fp->current_step = 0;
    fp->last_step = fh->tidx_stop - fh->tidx_start;

    fp->path = strdup (fh->fname);
    fp->endianness = bp_get_endianness (fh->mfooter.change_endianness);
    fp->version = fh->mfooter.version & ADIOS_VERSION_NUM_MASK;
    fp->file_size = fh->mfooter.file_size;

    return fp;
}

int adios_read_bp_staged_close (ADIOS_FILE *fp)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;

    if (pvt->aggregator_rank_array)
    {
        free (pvt->aggregator_rank_array);
        pvt->aggregator_rank_array = 0;
    }

    free (pvt);
    p->priv = 0;

    if (p->fh)
    {
        bp_close (fh);
        p->fh = 0;
    }

    if (p->varid_mapping)
    {
        free (p->varid_mapping);
        p->varid_mapping = 0;
    }

    if (p->local_read_request_list)
    {
        list_free_read_request (p->local_read_request_list);
        p->local_read_request_list = 0;
    }

    free (p);

    if (fp->var_namelist)
    {
        a2s_free_namelist (fp->var_namelist, fp->nvars);
        fp->var_namelist = 0;
    }

    if (fp->attr_namelist)
    {
        a2s_free_namelist (fp->attr_namelist, fp->nattrs);
        fp->attr_namelist = 0;
    }

    if (fp->path)
    {
        free (fp->path);
        fp->path = 0;
    }
    // internal_data field is taken care of by common reader layer
    free (fp);

    return 0;
}

/* Staged read method doesn't support streaming read yet
 */
int adios_read_bp_staged_advance_step (ADIOS_FILE *fp, int last, float timeout_sec)
{
    log_error ("adios_advance_step() is not supported in this method.\n");
    return 0;
}

/* Staged read method doesn't support streaming read yet
 */
void adios_read_bp_staged_release_step (ADIOS_FILE *fp)
{
    log_error ("adios_release_step() is not supported in this method.\n");
}

ADIOS_VARINFO * adios_read_bp_staged_inq_var_byid (const ADIOS_FILE * fp, int varid)
{
    return adios_read_bp_inq_var_byid (fp, varid);
}

int adios_read_bp_staged_inq_var_stat (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo, int per_step_stat, int per_block_stat)
{
    return 0;
}

int adios_read_bp_staged_inq_var_blockinfo (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    return 0;
}

int adios_read_bp_staged_schedule_read_byid (const ADIOS_FILE * fp,
                                             const ADIOS_SELECTION * sel,
                                             int varid,
                                             int from_steps,
                                             int nsteps,
                                             void * data
                                            )
{
    // simply call the 'simple reader' scheudle read routine. 
    return adios_read_bp_schedule_read_byid (fp, sel, varid, from_steps, nsteps, data);
}

int adios_read_bp_staged_perform_reads (const ADIOS_FILE *fp, int blocking)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    bp_proc_pvt_struct * pvt = (bp_proc_pvt_struct *) p->priv;
    read_request * r, * h;
    int i, count, varid, ndims, total_size, size;
    void * buf;

    // First populate the read request private struct for
    // each local read request. Have to do it here since the
    // schedule_read_byid is simply calling what is in simple bp reader.
    h = p->local_read_request_list;
    while (h)
    {
        h->priv = malloc (sizeof (rr_pvt_struct));
        assert (h->priv);
        ((rr_pvt_struct *) h->priv)->rank = pvt->rank;
        h = h->next;
    }

    size = calc_data_size (p);
    p->b = malloc (size);
    assert (p->b);
    buf = p->b;

    // count
    count = list_get_length (p->local_read_request_list);
    buffer_write (&buf, &count, 4);

    h = p->local_read_request_list;
    while (h)
    {
        varid = h->varid;
        ndims = h->sel->u.bb.ndim;

        buffer_write (&buf, &varid, 4);
        buffer_write (&buf, &h->from_steps, 4);
        buffer_write (&buf, &h->nsteps, 4);
        buffer_write (&buf, &ndims, 4);
        //TODO: bb only for now
        buffer_write (&buf, h->sel->u.bb.start, ndims * 8);
        buffer_write (&buf, h->sel->u.bb.count, ndims * 8);
        buffer_write (&buf, &h->datasize, 8);

        h = h->next;
    }

    int * sizes = malloc (pvt->group_size * 4);
    int * offsets = malloc (pvt->group_size * 4);
    void * recv_buffer = NULL;

    MPI_Gather (&size, 1, MPI_INT
               ,sizes, 1, MPI_INT
               ,pvt->aggregator_new_rank, pvt->new_comm
               );

    if (isAggregator (p))
    {
        total_size = 0;
        offsets[0] = 0;

        for (i = 0; i < pvt->group_size; i++)
        {
            total_size += sizes[i];
            if (i > 0)
            {
                offsets[i] = offsets[i - 1] + sizes[i - 1];
            }
        }

        recv_buffer = malloc (total_size);
        assert (recv_buffer);
    }

    MPI_Gatherv (p->b, size, MPI_BYTE
                ,recv_buffer, sizes, offsets
                ,MPI_BYTE, pvt->aggregator_new_rank, pvt->new_comm
                );

    if (isAggregator (p))
    {
        for (i = 1; i < pvt->group_size; i++)
        {
            parse_buffer (p, (char *)recv_buffer + offsets[i], pvt->aggregator_rank + i);
        }
        free (recv_buffer);

        process_read_requests (fp);
    }

    free (p->b);
    p->b = 0;
    free (sizes);
    free (offsets);

    if (isAggregator (p))
    {
        sort_read_requests (p);

        do_read (fp);

        send_read_data (p);
    }
    else
    {
        get_read_data (p);
    }

    free_proc_struct (p);

    /* 1. prepare all reads */
    // check if all user memory is provided for blocking read
    if (blocking)
    {
        r = p->local_read_request_list;
        while (r)
        {
            rr_pvt_struct * rr_pvt = (rr_pvt_struct *) r->priv;
            if (rr_pvt->rank == pvt->rank && !r->data)
            {
                adios_error (err_operation_not_supported,
                    "[%d] Blocking mode at adios_perform_reads() requires that user "
                    "provides the memory for each read request. Request for "
                    "variable %d was scheduled without user-allocated me mory\n",
                    pvt->rank, r->varid);
                return err_operation_not_supported;
            }

            r = r->next;
        }
    }
    else
    {
        return 0;
    }

    return 0;
}

/* Note: staged read method doesn't support check_reads so far.
 * It only supports the sceanario that user allocates all memory and subsequently 
 * does perform_reads with blocking flagged.
 */
int adios_read_bp_staged_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk)
{
    log_error ("adios_check_reads() is not supported in this method.\n");

    return 0;
}

int adios_read_bp_staged_get_attr_byid (const ADIOS_FILE * fp, int attrid, enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    return 0;
}

int  adios_read_bp_staged_get_dimension_order (const ADIOS_FILE *fp)
{
    BP_FILE * fh = GET_BP_FILE (fp);
    return is_fortran_file (fh);
}

void adios_read_bp_staged_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{
    adios_read_bp_reset_dimension_order (fp, is_fortran);
}

void adios_read_bp_staged_get_groupinfo (const ADIOS_FILE *fp, int *ngroups, char ***group_namelist, uint32_t **nvars_per_group, uint32_t **nattrs_per_group)
{
    BP_PROC * p;
    BP_FILE * fh;
    int i, j, offset;

    p = (BP_PROC *) fp->fh;
    fh = (BP_FILE *) p->fh;

    * ngroups = fh->gvar_h->group_count;

    *group_namelist = (char **) malloc (sizeof (char *) * fh->gvar_h->group_count);
    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        (*group_namelist)[i] = malloc (strlen (fh->gvar_h->namelist[i]) + 1);
        assert ((*group_namelist)[i]);

        memcpy ((*group_namelist)[i], fh->gvar_h->namelist[i], strlen (fh->gvar_h->
namelist[i]) + 1);
    }

    * nvars_per_group = (uint32_t *) malloc (fh->gvar_h->group_count * sizeof (uint32_t));
    assert (* nvars_per_group);

    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        (* nvars_per_group)[i] = fh->gvar_h->var_counts_per_group[i];
    }

    * nattrs_per_group = (uint32_t *) malloc (fh->gattr_h->group_count * sizeof (uint32_t));
    assert (* nattrs_per_group);

    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        offset = 0;
        for (j = 0; j < i; j++)
        {
            offset += fh->gattr_h->attr_counts_per_group[j];
        }

        (* nattrs_per_group)[i] = 0;
        for (j = 0; j < fh->gattr_h->attr_counts_per_group[i]; j++)
        {
            if (!show_hidden_attrs && strstr (fh->gattr_h->attr_namelist[offset + j], "__adios__"))
            {
            }
            else
            {
                (* nattrs_per_group)[i] ++;
            }
        }
    }

    return;
}

/* This is simply a copy of simle bp reader */
int adios_read_bp_staged_is_var_timed (const ADIOS_FILE *fp, int varid)
{
    BP_PROC * p;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * v;
    struct adios_index_characteristic_struct_v1 ch;
    int retval = 0, ndim, k;
    uint64_t gdims[32];

    p = (BP_PROC *) fp->fh;
    fh = (BP_FILE *) p->fh;

    v = bp_find_var_byid (fh, varid);
    ch = v->characteristics[0];
    ndim = ch.dims.count; //ndim possibly has 'time' dimension

    log_debug ("adios_read_bp_is_var_timed: varid = %d, ndim = %d\n", varid, ndim);

    if (ndim == 0)
    {
        return 0;
    }

    for (k = 0; k < ndim; k++)
    {
        gdims[k] = ch.dims.dims[k * 3 + 1];
    }
/*
    if (is_fortran_file (fh))
    {
        swap_order (ndim, gdims, &dummy);
    }
*/

    if (gdims[ndim - 1] == 0) // with time
    {
        if (v->characteristics_count <= 1) {
            // a local array written once
            retval = 0;
        } else {
            retval = 1;
        }
        /* FIXME: This last test tests if the last l:g:o is only an 'l'.
           This is true for a variable over time but also 
           true for a 1D local array (which has no global dimension)
           The characteristics_count is 1 only if the local array is written
           from one process and only at one timestep.
           How do we identify local arrays written from many processes?
           And local arrays written several times?
        */
    }

    log_debug ("%s is_var_timed: = %d\n", v->var_name, retval);

    return retval;
}
