/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */


/****************************/
/* Read method for BP files */
/****************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include "public/adios_read.h"
#include "public/adios_error.h"
#include "public/adios_types.h"
#include "core/util.h"
#include "core/bp_utils.h"
#include "core/bp_types.h"
#include "core/adios_read_hooks.h"
#include "core/adios_read_transformed.h"
#include "core/futils.h"
#include "core/common_read.h"
#include "core/adios_logger.h"

// LAYERFIX
// NCSU ALACRITY-ADIOS - Include necessary headers
//#include "core/adios_transforms_common.h"
//#include "core/adios_transforms_read.h"
//#include "core/adios_transforms_util.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

static int chunk_buffer_size = 1024*1024*16;
static int poll_interval = 10; // 10 secs by default
static int show_hidden_attrs = 0; // don't show hidden attr by default

static ADIOS_VARCHUNK * read_var_bb (const ADIOS_FILE * fp, read_request * r);
static int bp_close (BP_FILE * fh);
static int adios_read_bp_get_endianness( uint32_t change_endianness);

// NCSU - For custom memory allocation
#define CALLOC(var, num, sz, comment)\
{\
    var = calloc (num, sz); \
    if (!var)    {\
        adios_error_at_line (err_no_memory, __FILE__, __LINE__, "Could not allocate memory for ", comment, " in common_read_get_characteristics"); \
        return 0; \
    }\
}

#define MALLOC(var,sz,comment)\
{\
    var = malloc (sz); \
    if (!var)    {\
        adios_error_at_line (err_no_memory, __FILE__, __LINE__, "Could not allocate memory for ", comment, " in common_read_get_characteristics"); \
        return 0; \
    }\
}\

// Search for the start var index.
static int64_t get_var_start_index (struct adios_index_var_struct_v1 * v, int t)
{
    int64_t i = 0;

    while (i < v->characteristics_count) {
        if (v->characteristics[i].time_index == t) {
            return i;
        }

        i++;
    }

    return -1;
}

// Search for the stop var index
static int64_t get_var_stop_index (struct adios_index_var_struct_v1 * v, int t)
{
    int64_t i = v->characteristics_count - 1;

    while (i > -1) {
        if (v->characteristics[i].time_index == t) {
            return i;
        }

        i--;
    }

    return -1;
}

MPI_File * get_BP_file_handle(struct BP_file_handle * l, uint32_t file_index)
{
    if (!l)
        return 0;

    while (l)
    {
        if (l->file_index == file_index)
            return &l->fh;

        l = l->next;
    }

    return 0;
}

void add_BP_file_handle (struct BP_file_handle ** l, struct BP_file_handle * n)
{
    if (!n)
        return;

    n->next = *l;
    *l = n;
}


void close_all_BP_files (struct BP_file_handle * l)
{
    struct BP_file_handle * n;

    while (l)
    {
        n = l->next;

        MPI_File_close (&l->fh);
        free (l);

        l = n;
    }
}

#define MPI_FILE_READ_OPS1                          \
        bp_realloc_aligned(fh->b, slice_size);      \
        fh->b->offset = 0;                          \
        MPI_File_seek (fh->mpi_fh                   \
                      ,(MPI_Offset)slice_offset     \
                      ,MPI_SEEK_SET                 \
                      );                            \
                                                    \
        MPI_File_read (fh->mpi_fh                   \
                      ,fh->b->buff                  \
                      ,slice_size                   \
                      ,MPI_BYTE                     \
                      ,&status                      \
                      );                            \
        fh->b->offset = 0;                          \

// To read subfiles
#define MPI_FILE_READ_OPS2                                                                  \
        bp_realloc_aligned(fh->b, slice_size);                                              \
        fh->b->offset = 0;                                                                  \
                                                                                            \
        MPI_File * sfh;                                                                     \
        sfh = get_BP_file_handle (fh->sfh                                                   \
                                 ,v->characteristics[start_idx + idx].file_index            \
                                 );                                                         \
        if (!sfh)                                                                           \
        {                                                                                   \
            int err;                                                                        \
            char * ch, * name_no_path, * name;                                              \
            MPI_Info info;                                                                  \
            struct BP_file_handle * new_h =                                                 \
                  (struct BP_file_handle *) malloc (sizeof (struct BP_file_handle));        \
            new_h->file_index = v->characteristics[start_idx + idx].file_index;             \
            new_h->next = 0;                                                                \
            if (ch = strrchr (fh->fname, '/'))                                              \
            {                                                                               \
                name_no_path = (char *) malloc (strlen (ch + 1) + 1);                       \
                strcpy (name_no_path, ch + 1);                                              \
            }                                                                               \
            else                                                                            \
            {                                                                               \
                name_no_path = (char *) malloc (strlen (fh->fname) + 1);                    \
                strcpy (name_no_path, fh->fname);                                           \
            }                                                                               \
                                                                                            \
            name = (char *) malloc (strlen (fh->fname) + 5 + strlen (name_no_path) + 1 + 10 + 1); \
            sprintf (name, "%s.dir/%s.%d", fh->fname, name_no_path, new_h->file_index);     \
            err = MPI_File_open (MPI_COMM_SELF                                                   \
                                ,name                                                       \
                                ,MPI_MODE_RDONLY                                            \
                                ,info                                                       \
                                ,&new_h->fh                                                 \
                                );                                                          \
                                                                                            \
           if (err)                                                                         \
           {                                                                                \
               fprintf (stderr, "can not open file %s\n", name);                            \
               return 0;                                                                    \
           }                                                                                \
           add_BP_file_handle (&fh->sfh                                                     \
                              ,new_h                                                        \
                              );                                                            \
           sfh = &new_h->fh;                                                                \
                                                                                            \
           free (name_no_path);                                                             \
           free (name);                                                                     \
        }                                                                                   \
                                                                                            \
        MPI_File_seek (*sfh                                                                 \
                      ,(MPI_Offset)slice_offset                                             \
                      ,MPI_SEEK_SET                                                         \
                      );                                                                    \
        MPI_File_read (*sfh                                                                 \
                      ,fh->b->buff                                                          \
                      ,slice_size                                                           \
                      ,MPI_BYTE                                                             \
                      ,&status                                                              \
                      );                                                                    \
        fh->b->offset = 0;                                                                  \

/* This routine release one step. It only frees the var/attr namelist. */
static void release_step (ADIOS_FILE *fp)
{
    BP_PROC * p = (BP_PROC *) fp->fh;

    if (p->varid_mapping)
    {
        free (p->varid_mapping);
        p->varid_mapping = 0;
    }

    if (fp->var_namelist)
    {
        free_namelist (fp->var_namelist, fp->nvars);
        fp->var_namelist = 0;
        fp->nvars = 0;
    }

    if (fp->attr_namelist)
    {
        free_namelist (fp->attr_namelist, fp->nattrs);
        fp->attr_namelist = 0;
        fp->nattrs = 0;
    }
}

/* This routin open a ADIOS-BP file with no timeout.
 * It first checks whether this is a valid BP file. This is done by
 * checking the validity on rank 0 and communicating to other ranks (avoiding
 * the situation where every one kicks the tire.
 * If the file is ok, then go ahead get the metadata.
 */
static BP_FILE * open_file (const char * fname, MPI_Comm comm)
{
    int i, rank, file_ok;
    BP_FILE * fh;

    MPI_Comm_rank (comm, &rank);

    if (rank == 0)
    {
        file_ok = check_bp_validity (fname);

        MPI_Bcast (&file_ok, 1, MPI_INT, 0, comm);
    }
    else
    {
        MPI_Bcast (&file_ok, 1, MPI_INT, 0, comm);
    }

    if (!file_ok)
    {
        return 0;
    }

    fh = (BP_FILE *) malloc (sizeof (BP_FILE));
    assert (fh);

    fh->fname = (fname ? strdup (fname) : 0L);
    fh->sfh = 0;
    fh->comm = comm;
    fh->gvar_h = 0;
    fh->pgs_root = 0;
    fh->vars_root = 0;
    fh->attrs_root = 0;
    fh->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    assert (fh->b);

    bp_open (fname, comm, fh);

    return fh;
}

/* This routine set ADIOS_FILE fields from BP_FILE struct */
void build_ADIOS_FILE_struct (ADIOS_FILE * fp, BP_FILE * fh)
{
    BP_PROC * p;
    int rank;

    log_debug ("build_ADIOS_FILE_struct is called\n");

    MPI_Comm_rank (fh->comm, &rank);

    p = (struct BP_PROC *) malloc (sizeof (struct BP_PROC));
    assert (p);
    p->rank = rank;
    p->fh = fh;
    p->streaming = 1;
    p->local_read_request_list = 0;
    p->priv = 0;

    fp->fh = (uint64_t) p;
    fp->file_size = fh->mfooter.file_size;
    fp->version = fh->mfooter.version;
    fp->endianness = adios_read_bp_get_endianness (fh->mfooter.change_endianness);
    fp->last_step = fh->tidx_stop - 1;

    /* Seek to the initial step. */
    release_step (fp);
    bp_seek_to_step (fp, 0, show_hidden_attrs);

    /* For file, the last step is tidx_stop */
    fp->last_step = fh->tidx_stop - 1;

    return;
}

static int get_new_step (ADIOS_FILE * fp, const char * fname, MPI_Comm comm, int last_step, float timeout_sec)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    int err, i;
    BP_FILE * fh;
    double t1 = adios_gettime();

    log_debug ("enter get_new_step\n");
    /* First check if the file has been updated with more steps. */
    /* While loop for handling timeout
       timeout > 0: wait up to this long to open the stream
       timeout = 0: return immediately
       timeout < 0: wait forever
    */
    int stay_in_poll_loop = 1;
    int found_stream = 0;

    while (stay_in_poll_loop)
    {
        /* Re-open the file */
        fh = open_file (fname, comm);
        if (!fh)
        {
            // file is bad so keeps polling.
            stay_in_poll_loop = 1;
        }
        else if (fh && (fh->tidx_stop - 1 == last_step))
        {
            // file is good but no new steps in it. Continue polling.
            bp_close (fh);
            stay_in_poll_loop = 1;
        }
        else
        {
            // the file looks good and there are new steps written.
            build_ADIOS_FILE_struct (fp, fh);
            stay_in_poll_loop = 0;
            found_stream = 1;
        }
        // check if we need to stay in loop
        if (stay_in_poll_loop)
        {
            if (timeout_sec == 0.0)
            {
                stay_in_poll_loop = 0;
            }
            else if (timeout_sec < 0.0)
            {
                stay_in_poll_loop = 1;
            }
            else if (timeout_sec > 0.0 && (adios_gettime () - t1 > timeout_sec))
            {
                log_debug ("Time is out in get_new_step()\n");
                stay_in_poll_loop = 0;
            }
            else
            {
                adios_nanosleep (poll_interval, 0);
            }
        }

    } // while (stay_in_poll_loop)

    log_debug ("exit get_new_step\n");

    return found_stream;
}

/* This routine processes a read request and returns data in ADIOS_VARCHUNK.
   If the selection type is not bounding box, convert it. The basic file reading
   functionality is implemented in read_var_bb() routine.
*/
static ADIOS_VARCHUNK * read_var (const ADIOS_FILE * fp, read_request * r)
{
    BP_PROC * p;
    BP_FILE * fh;
    int size_of_type;
    struct adios_index_var_struct_v1 * v;
    uint64_t i;
    read_request * nr;
    ADIOS_SELECTION * sel, * nsel;
    ADIOS_VARCHUNK * chunk;

    sel = r->sel;
    p = (BP_PROC *) fp->fh;
    fh = (BP_FILE *) p->fh;

    v = bp_find_var_byid (fh, r->varid);

    switch (sel->type)
    {
        case ADIOS_SELECTION_BOUNDINGBOX:
            chunk = read_var_bb (fp, r);
            break;
        case ADIOS_SELECTION_POINTS:
            /* The idea is we convert a point selection to bounding box section. */
            size_of_type = bp_get_type_size (v->type, v->characteristics [0].value);
            nr = (read_request *) malloc (sizeof (read_request));
            assert (nr);

            nr->rank = r->rank;
            nr->varid = r->varid;
            nr->from_steps = r->from_steps;
            nr->nsteps = r->nsteps;
            nr->data = r->data;
            nr->datasize  = r->datasize;
            nr->priv = r->priv;

            nsel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
            nr->sel = nsel;
            assert (nsel);

            nsel->type = ADIOS_SELECTION_BOUNDINGBOX;
            nsel->u.bb.ndim = sel->u.points.ndim;
            nsel->u.bb.start = (uint64_t *) malloc (nsel->u.bb.ndim * 8);
            nsel->u.bb.count = (uint64_t *) malloc (nsel->u.bb.ndim * 8);
            assert (nsel->u.bb.start && nsel->u.bb.count);

            for (i = 0; i < nsel->u.bb.ndim; i++)
            {
                nsel->u.bb.count[i] = 1;
            }

            for (i = 0; i < sel->u.points.npoints; i++)
            {
                memcpy (nsel->u.bb.start, sel->u.points.points + i * sel->u.points.ndim, sel->u.points.ndim * 8);

                chunk = read_var_bb (fp, nr);
                nr->data = (char *) nr->data + size_of_type;
            }

            free_selection (nsel);
            free (nr);
            break;
        case ADIOS_SELECTION_WRITEBLOCK:
            break;
        case ADIOS_SELECTION_AUTO:
            break;
    }

    return chunk;
}

/* Convert 'step' to time, which is used in ADIOS internals.
 * 'step' should start from 0.
 */
static int get_time (struct adios_index_var_struct_v1 * v, int step)
{
    int i = 0;
    int prev_ti = 0, counter = 0;

    while (i < v->characteristics_count)
    {
        if (v->characteristics[i].time_index != prev_ti)
        {
            counter ++;
            if (counter == (step + 1))
            {
                return v->characteristics[i].time_index;
            }
            prev_ti = v->characteristics[i].time_index;
        }

        i++;
    }

    return -1;

}

// LAYERFIX
/*
// NCSU ALACRITY-ADIOS - This state struct is given to the transform module when
//   reading a subvolume, to be passed to the file reading delegate function.
//   It is specific to the read_bp.c reader implementation.
typedef struct {
    BP_FILE *fh;
    int has_subfile;
} read_bp_read_state;

static void * read_bp_read_delegate(struct adios_index_var_struct_v1 *v, int time_index,
                                    uint64_t slice_offset, uint64_t slice_size, void *state_tmp) {
    read_bp_read_state *state = (read_bp_read_state *)state_tmp;
    BP_FILE *fh = state->fh;
    MPI_Status status; 			// For MPI_FILE_READ_OPS[12]
    int start_idx = time_index;	// For MPI_FILE_READ_OPS2
    int idx = 0;				// For MPI_FILE_READ_OPS2

    slice_offset += v->characteristics[start_idx + idx].payload_offset;
    if (!state->has_subfile)
    {
        MPI_FILE_READ_OPS1
    }
    else
    {
        MPI_FILE_READ_OPS2
    }

    return fh->b->buff + fh->b->offset;
}
*/

/* This routines read in data for bounding box selection.
   If the selection is not bounding box, it should be converted to it.
   The data returned is saved in ADIOS_VARCHUNK.
 */
static ADIOS_VARCHUNK * read_var_bb (const ADIOS_FILE *fp, read_request * r)
{
    BP_PROC * p;
    BP_FILE * fh;
    ADIOS_SELECTION * sel;
    struct adios_index_var_struct_v1 * v;
    int i, j, k, t, time, nsteps;
    int64_t start_idx, stop_idx, idx;
    int ndim, has_subfile, file_is_fortran;
    uint64_t size, * dims;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t datasize, dset_stride,var_stride, total_size=0, items_read;
    uint64_t * count, * start;
    void * data;
    int dummy = -1, is_global = 0, size_of_type;
    uint64_t slice_offset, slice_size;
    MPI_Status status;
    ADIOS_VARCHUNK * chunk;

    p = (BP_PROC *) fp->fh;
    fh = (BP_FILE *) p->fh;
    file_is_fortran = is_fortran_file (fh);
    has_subfile = has_subfiles (fh);

    sel = r->sel;
    start = sel->u.bb.start;
    count = sel->u.bb.count;
    data = r->data;

    v = bp_find_var_byid (fh, r->varid);

    /* Get dimensions and flip if caller != writer language */
    /* Note: ndim below does include time if there is any */
    // NCSU ALACRITY-ADIOS - Note: this function has been modified to return
    //   the "original" dimensions of the variable (i.e., not 1D byte array)
    bp_get_and_swap_dimensions (fh, v, file_is_fortran, &ndim, &dims, &nsteps, file_is_fortran);

    assert (ndim == sel->u.bb.ndim);
    ndim = sel->u.bb.ndim;

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if (futils_is_called_from_fortran ())
    {
        swap_order (ndim, start, &dummy);
        swap_order (ndim, count, &dummy);
    }

    /* items_read = how many data elements are we going to read in total (per timestep) */
    items_read = 1;
    for (i = 0; i < ndim; i++)
    {
        items_read *= count[i];
    }

    size_of_type = bp_get_type_size (v->type, v->characteristics [0].value);

    log_debug ("read_var_bb: from_steps = %d, nsteps = %d\n", r->from_steps, r->nsteps);

    /* Note fp->current_step is always 0 for files. */
    for (t = fp->current_step + r->from_steps; t < fp->current_step + r->from_steps + r->nsteps; t++)
    {

        if (!p->streaming)
        {
            time = get_time (v, t);
        }
        else
        {
            time = t + 1;
        }
//printf ("t = %d(%d,%d), time = %d\n", t, fp->current_step, r->from_steps, time);
//printf ("c = %d, f = %d, time = %d\n", fp->current_step, r->from_steps, time);
        start_idx = get_var_start_index (v, time);
        stop_idx = get_var_stop_index (v, time);

        if (start_idx < 0 || stop_idx < 0)
        {
            adios_error (err_no_data_at_timestep,"Variable %s has no data at %d time step\n",
                         v->var_name, t);
            continue;
        }

        if (ndim == 0)
        {
            /* READ A SCALAR VARIABLE */
            /* Prepare slice_offset, slice_size and idx for the later macro:
               MPI_FILE_READ_OPS1 and MPI_FILE_READ_OPS2
            */
            idx = 0;
            slice_offset = v->characteristics[start_idx + idx].payload_offset;
            slice_size = size_of_type;

            if (v->type == adios_string)
            {
                // Note: strings are stored without \0 in file
                // size_of_type above includes \0 so decrease by one
                size_of_type--;
            }

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS1
            }
            else
            {
                MPI_FILE_READ_OPS2
            }

            memcpy ((char *)data, fh->b->buff + fh->b->offset, size_of_type);

            if (fh->mfooter.change_endianness == adios_flag_yes)
            {
                change_endianness ((char *)data, size_of_type, v->type);
            }

            if (v->type == adios_string)
            {
                // add \0 to the end of string
                // size_of_type here is the length of string
                // FIXME: how would this work for strings written over time?
                ((char*)data)[size_of_type] = '\0';
            }

            data = (char *)data + (v->type == adios_string?  size_of_type + 1 : size_of_type);
        }
        else
        {
            /* READ AN ARRAY VARIABLE */
            int * idx_table = (int *) malloc (sizeof (int) * (stop_idx - start_idx + 1));
            uint64_t write_offset = 0;

            // loop over the list of pgs to read from one-by-one
            for (idx = 0; idx < stop_idx - start_idx + 1; idx++)
            {
                int flag;
                datasize = 1;
                var_stride = 1;
                dset_stride = 1;
                idx_table[idx] = 1;
                uint64_t payload_size = size_of_type;

                is_global = bp_get_dimension_characteristics_notime (&(v->characteristics[start_idx + idx]),
                                                                    ldims, gdims, offsets, file_is_fortran);
                if (!is_global)
                {
                    // we use gdims below, which is 0 for a local array; set to ldims here
                    for (j = 0; j < ndim; j++)
                    {
                        gdims[j] = ldims[j];
                    }
                    // we need to read only the first PG, not all, so let's prevent a second loop
                    stop_idx = start_idx;
                }
/*
                log_debug ("ldims   = "); for (j = 0; j<ndim; j++) log_debug_cont ("%d ",ldims[j]); log_debug_cont ("\n");
                log_debug ("gdims   = "); for (j = 0; j<ndim; j++) log_debug_cont ("%d ",gdims[j]); log_debug_cont ("\n");
                log_debug ("offsets = "); for (j = 0; j<ndim; j++) log_debug_cont ("%d ",offsets[j]); log_debug_cont ("\n");
                log_debug ("count   = "); for (j = 0; j<ndim; j++) log_debug_cont ("%d ",count[j]); log_debug_cont ("\n");
                log_debug ("start   = "); for (j = 0; j<ndim; j++) log_debug_cont ("%d ",start[j]); log_debug_cont ("\n");
*/
                for (j = 0; j < ndim; j++)
                {
                    payload_size *= ldims [j];

                    if ( (count[j] > gdims[j])
                      || (start[j] > gdims[j])
                      || (start[j] + count[j] > gdims[j]))
                    {
                        adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound 1("
                            "the data in dimension %d to read is %llu elements from index %llu"
                            " but the actual data is [0,%llu])\n",
                            r->varid, j + 1, count[j], start[j], gdims[j] - 1);
                        return 0;
                    }

                    /* check if there is any data in this pg and this dimension to read in */
                    flag = (offsets[j] >= start[j]
                            && offsets[j] < start[j] + count[j])
                        || (offsets[j] < start[j]
                            && offsets[j] + ldims[j] > start[j] + count[j])
                        || (offsets[j] + ldims[j] > start[j]
                            && offsets[j] + ldims[j] <= start[j] + count[j]);
                    idx_table [idx] = idx_table[idx] && flag;
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
                        datasize *= ldims[i];
                    }
                    else
                        break;
                }

                hole_break = i;
                slice_offset = 0;
                slice_size = 0;

/*
                log_debug ("hole_break = %d\n", hole_break);
*/

                // LAYERFIX
                /*
                // NCSU ALACRITY-ADIOS - If the variable is transformed, skip the normal read code
                //   and delegate to the transform method
                if (adios_transform_var_is_transformed(v))
                {
                    int it;
                    // Build the reader state
                    read_bp_read_state state;
                    state.fh = fh;
                    state.has_subfile = has_subfile;

                    assert(sel->type == ADIOS_SELECTION_BOUNDINGBOX);

//                    printf("Selection box: [");
//                    for (it = 0; it < sel->u.bb.ndim; it++)
//                        printf(" %llu", sel->u.bb.start[it]);
//                    printf(" ] [");
//                    for (it = 0; it < sel->u.bb.ndim; it++)
//                        printf(" %llu", sel->u.bb.count[it]);
//                    printf(" ]\n");

                    adios_subvolume_copy_spec copy_spec;
                    int intersects = adios_selection_to_copy_spec(&copy_spec, &sel->u.bb, ldims, offsets);
//                    printf(">>> copyspec dim1: subvolume %llu from src start:offset %llu:%llu to dst start:offset %llu:%llu\n",
//                            copy_spec.subv_dims[0], copy_spec.src_subv_offsets[0], copy_spec.src_dims[0], copy_spec.dst_subv_offsets[0], copy_spec.dst_dims[0]);
                    if (intersects) {
                        // Delegate reading to the transform code, which will in
                        // turn delegate reading to the supplied reader delegate
                        int err = adios_transform_retrieve_transformed_data_subvolume(v, data, start_idx + idx, &copy_spec, &state, read_bp_read_delegate, fh->mfooter.change_endianness);

                        // TODO: Doesn't work on non-contiguous subvolumes...should do this in copy_subvolume

                        // Reverse endianness, if required
                        //if (fh->mfooter.change_endianness == adios_flag_yes)
                        //{
                        //    change_endianness (data, slice_size, adios_transform_get_var_original_type(v));
                        //}

                    }
                    adios_subvolume_copy_spec_destroy(&copy_spec, 1);
                }
                else */
                if (hole_break == -1)
                {
                    /* The complete read happens to be exactly one pg, and the entire pg */
                    /* This means we enter this only once, and npg=1 at the end */
                    /* This is a rare case. FIXME: cannot eliminate this? */
                    slice_size = payload_size;

                    slice_offset = v->characteristics[start_idx + idx].payload_offset;
                    if (!has_subfile)
                    {
                        MPI_FILE_READ_OPS1
                    }
                    else
                    {
                        MPI_FILE_READ_OPS2
                    }

                    memcpy ((char *)data, fh->b->buff + fh->b->offset, slice_size);
                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {
                        change_endianness (data, slice_size, v->type);
                    }
                }
                else if (hole_break == 0)
                {
                    /* The slowest changing dimensions should not be read completely but
                       we still need to read only one block */
                    int isize;
                    uint64_t size_in_dset = 0;
                    uint64_t offset_in_dset = 0;
                    uint64_t offset_in_var = 0;

                    isize = offsets[0] + ldims[0];
                    if (start[0] >= offsets[0])
                    {
                        // head is in
                        if (start[0]<isize)
                        {
                            if (start[0] + count[0] > isize)
                                size_in_dset = isize - start[0];
                            else
                                size_in_dset = count[0];
                            offset_in_dset = start[0] - offsets[0];
                            offset_in_var = 0;
                        }
                    }
                    else
                    {
                        // middle is in
                        if (isize < start[0] + count[0])
                            size_in_dset = ldims[0];
                        else
                        // tail is in
                            size_in_dset = count[0] + start[0] - offsets[0];
                        offset_in_dset = 0;
                        offset_in_var = offsets[0] - start[0];
                    }

                    slice_size = size_in_dset * datasize * size_of_type;
                    write_offset = offset_in_var * datasize * size_of_type;

                    if (v->characteristics[start_idx + idx].payload_offset > 0)
                    {
                        slice_offset = v->characteristics[start_idx + idx].payload_offset
                                     + offset_in_dset * datasize * size_of_type;
                        if (!has_subfile)
                        {
                            MPI_FILE_READ_OPS1
                        }
                        else
                        {
                            MPI_FILE_READ_OPS2
                        }
                    }

                    memcpy ((char *)data + write_offset, fh->b->buff + fh->b->offset, slice_size);
                    if (fh->mfooter.change_endianness == adios_flag_yes)
                    {
                        change_endianness((char *)data + write_offset, slice_size, v->type);
                    }

                    //write_offset +=  slice_size;
                }
                else
                {
                    uint64_t stride_offset = 0;
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
                            if (start[i]<isize)
                            {
                                if (start[i] + count[i] > isize)
                                    size_in_dset[i] = isize - start[i];
                                else
                                    size_in_dset[i] = count[i];
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

                    for (i = ndim - 1; i >= hole_break; i--)
                    {
                        datasize *= size_in_dset[i];
                        dset_stride *= ldims[i];
                        var_stride *= count[i];
                    }

                    uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
                    for (i = ndim - 1; i > -1; i--)
                    {
                        start_in_payload += s * offset_in_dset[i] * size_of_type;
                        end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
                        s *= ldims[i];
                    }

                    slice_size = end_in_payload - start_in_payload + 1 * size_of_type;

                    slice_offset =  v->characteristics[start_idx + idx].payload_offset
                                  + start_in_payload;
                    if (!has_subfile)
                    {
                        MPI_FILE_READ_OPS1
                    }
                    else
                    {
                       MPI_FILE_READ_OPS2
                    }

                    for (i = 0; i < ndim; i++)
                    {
                        offset_in_dset[i] = 0;
                    }

                    uint64_t var_offset = 0;
                    uint64_t dset_offset = 0;

                    for (i = 0; i < ndim; i++)
                    {
                        var_offset = offset_in_var[i] + var_offset * count[i];
                        dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
                    }

                    copy_data (data
                              ,fh->b->buff + fh->b->offset
                              ,0
                              ,hole_break
                              ,size_in_dset
                              ,ldims
                              ,count
                              ,var_stride
                              ,dset_stride
                              ,var_offset
                              ,dset_offset
                              ,datasize
                              ,size_of_type
                              );
                }
            }  // end for (idx ... loop over pgs

            free (idx_table);

            total_size += items_read * size_of_type;
            // shift target pointer for next read in
            data = (char *)data + (items_read * size_of_type);

//FIXME
/*
        if (timedim == -1)
            break;
*/
        }
    } // end for t

    free (dims);

    chunk = (ADIOS_VARCHUNK *) malloc (sizeof (ADIOS_VARCHUNK));
    assert (chunk);

    chunk->varid = r->varid;
    chunk->type = v->type;
    chunk->sel = copy_selection (r->sel);
    chunk->data = data;

    return chunk;
}

/* Return 0: if file is little endian, 1 if file is big endian
 * We know if it is different from the current system, so here
 * we determine the current endianness and report accordingly.
 */
static int adios_read_bp_get_endianness( uint32_t change_endianness )
{
   int LE = 0;
   int BE = !LE;
   int i = 1;
   char *p = (char *) &i;
   int current_endianness;
   if (p[0] == 1) // Lowest address contains the least significant byte
       current_endianness = LE;
   else
       current_endianness = BE;
    if (change_endianness == adios_flag_yes)
        return !current_endianness;
    else
        return current_endianness;
}

int adios_read_bp_init_method (MPI_Comm comm, PairStruct * params)
{
    int  max_chunk_size;
    PairStruct * p = params;

    while (p)
    {
        if (!strcasecmp (p->name, "max_chunk_size"))
        {
            max_chunk_size = strtol(p->value, NULL, 10);
            if (max_chunk_size > 0)
            {
                log_debug ("max_chunk_size set to %dMB for the read method\n", max_chunk_size);
                chunk_buffer_size = max_chunk_size * 1024 * 1024;
            }
            else
            {
                log_error ("Invalid 'max_chunk_size' parameter given to the read method: '%s'\n", p->value);
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

        p = p->next;
    }

    return 0;
}

int adios_read_bp_finalize_method ()
{
    /* Set these back to default */
    chunk_buffer_size = 1024*1024*16;
    poll_interval = 10; // 10 secs by default
    show_hidden_attrs = 0; // don't show hidden attr by default

    return 0;
}

static int open_stream (ADIOS_FILE * fp, const char * fname,
                        MPI_Comm comm, float timeout_sec)
{
    int i, rank, ret;
    struct BP_PROC * p;
    BP_FILE * fh;
    int err, stay_in_poll_loop = 1;
    int file_ok = 0;
    double t1 = adios_gettime();

    MPI_Comm_rank (comm, &rank);
    // We need to first check if this is a valid ADIOS-BP file. This is done by
    // check whether there is 'ADIOS-BP' string written before the 28-bytes minifooter.
    // If it is valid, we will proceed with bp_open(). The potential issue is that before
    // calling bp_open, the next step could start writing and the footer will be corrupted.
    // This needs to be fixed later. Q. Liu, 06/2012
    /* While loop for handling timeout
       timeout > 0: wait up to this long to open the stream
       timeout = 0: return immediately
       timeout < 0: wait forever
    */

    // Only rank 0 does the poll
    if (rank == 0)
    {
        while (stay_in_poll_loop)
        {
            file_ok = check_bp_validity (fname);
            if (!file_ok)
            {
                // This stream does not exist yet
                log_debug ("Warning: file %s does not exsit!\n", fname);
            }

            if (stay_in_poll_loop)
            {
                // check if we need to stay in loop
                if (timeout_sec == 0.0)  //return immediately, which means check file once
                {
                    stay_in_poll_loop = 0;
                }
                else if (timeout_sec < 0.0) // check file until it arrives
                {
                    adios_nanosleep (poll_interval, 0);
                    stay_in_poll_loop = 1;
                }
                else if (timeout_sec > 0.0 && (adios_gettime () - t1 > timeout_sec))
                {
                    stay_in_poll_loop = 0;
                }
                else
                {
                    adios_nanosleep (poll_interval, 0);
                }
            }
        } // while (stay_in_poll_loop)

        if (!file_ok)
        {
            adios_error (err_file_open_error, "File not found: %s\n", fname);
        }

        MPI_Bcast (&file_ok, 1, MPI_INT, 0, comm);
    }
    else
    {
        MPI_Bcast (&file_ok, 1, MPI_INT, 0, comm);
    }

    if (!file_ok)
    {
        return err_file_not_found;
    }

    fh = (BP_FILE *) malloc (sizeof (BP_FILE));
    assert (fh);

    fh->fname = (fname ? strdup (fname) : 0L);
    fh->sfh = 0;
    fh->comm = comm;
    fh->gvar_h = 0;
    fh->pgs_root = 0;
    fh->vars_root = 0;
    fh->attrs_root = 0;
    fh->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    assert (fh->b);

    p = (struct BP_PROC *) malloc (sizeof (struct BP_PROC));
    assert (p);
    p->rank = rank;
    p->fh = fh;
    p->streaming = 1;
    p->varid_mapping = 0;
    p->local_read_request_list = 0;
    p->priv = 0;

    /* BP file open and gp/var/att parsing */
    bp_open (fname, comm, fh);

    fp->fh = (uint64_t) p;
    fp->file_size = fh->mfooter.file_size;
    fp->version = fh->mfooter.version;
    fp->path = strdup (fh->fname);
    fp->endianness = adios_read_bp_get_endianness (fh->mfooter.change_endianness);

    /* Seek to the initial step. */
    bp_seek_to_step (fp, 0, show_hidden_attrs);

    fp->current_step = 0;
    /* For file, the last step is tidx_stop */
    fp->last_step = fh->tidx_stop - 1;

    return 0;
}

/* As opposed to open_file, open_stream opens the first step in the file only.
 * The lock_mode for file reading is ignored for now.
 */
ADIOS_FILE * adios_read_bp_open_stream (const char * fname, MPI_Comm comm, enum ADIOS_LOCKMODE lock_mode, float timeout_sec)
{
    log_debug ("adios_read_bp_open_stream\n");

    ADIOS_FILE * fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    assert (fp);

    if (open_stream (fp, fname, comm, timeout_sec) < 0)
    {
        free (fp);
        fp = 0;
    }

    return fp;
}

ADIOS_FILE * adios_read_bp_open_file (const char * fname, MPI_Comm comm)
{
    int i, rank;
    struct BP_PROC * p;
    BP_FILE * fh;
    ADIOS_FILE * fp;

    log_debug ("adios_read_bp_open_file\n");

    MPI_Comm_rank (comm, &rank);

    fh = (BP_FILE *) malloc (sizeof (BP_FILE));
    assert (fh);

    fh->fname = (fname ? strdup (fname) : 0L);
    fh->sfh = 0;
    fh->comm = comm;
    fh->gvar_h = 0;
    fh->pgs_root = 0;
    fh->vars_root = 0;
    fh->attrs_root = 0;
    fh->b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    assert (fh->b);

    p = (struct BP_PROC *) malloc (sizeof (struct BP_PROC));
    assert (p);
    p->rank = rank;
    p->fh = fh;
    p->streaming = 0;
    p->varid_mapping = 0; // maps perceived id to real id
    p->local_read_request_list = 0;
    p->priv = 0;

    /* The ADIOS_FILE struct looks like the following */
#if 0
typedef struct {
        uint64_t fh;                /* File handler                                                   */
        int      nvars;             /* Number of variables in all groups (with full path)             */
        char     ** var_namelist;   /* Variable names in a char* array                                */
        int      nattrs;            /* Number of attributes in all groups                             */
        char     ** attr_namelist;  /* Attribute names in a char* array                               */

        /* Step */
        int      current_step;      /* The current step                                               */
        int      last_step;         /* The currently available latest step in the stream/file         */

        /* Information about file/stream */
        char     *path;             /* Full path file name (as passed at open)                        */
        int      endianness;        /* 0: little endian, 1: big endian                                */
                                    /*   the read API takes care of conversion automatically          */
        int      version;           /* Version of ADIOS-BP format                                     */
        uint64_t file_size;         /* Size of file in bytes not including subfiles                   */

        /* Internals */
        void     * internal_data;   /* Data for internal use                                          */
} ADIOS_FILE;
#endif
    fp = (ADIOS_FILE *) malloc (sizeof (ADIOS_FILE));
    assert (fp);

    if (bp_open (fname, comm, fh) < 0)
    {
        adios_error (err_file_open_error, "File open failed: %s", fname);
        return 0;
    }

    /* fill out ADIOS_FILE struct */
    fp->fh = (uint64_t) p;

    /* '-1' means that we want all steps.
     * This will seek to the last step. So we need to set current_step back properly.
     * Usually bp_seek_to_step comes after release_step call, to first free up some
     * memory allocated by the previous step. This is the first seek call and, therefore,
     * no release_step.
     */
    bp_seek_to_step (fp, -1, show_hidden_attrs);
    /* It was agreed that, for file open the current step should be set to 0,
     * instead of the start time. The var_namelist and attr_namelist should
     * consist of all steps. For stream open, this is done differently.
     * 07/2012 - Q.Liu
     */
    fp->current_step = 0;
    fp->last_step = fh->tidx_stop - fh->tidx_start;

    fp->path = strdup (fh->fname);
    fp->endianness = adios_read_bp_get_endianness (fh->mfooter.change_endianness);
    fp->version = fh->mfooter.version;
    fp->file_size = fh->mfooter.file_size;

    return fp;
}

int adios_read_bp_close (ADIOS_FILE * fp)
{
    struct BP_PROC * p;
    BP_FILE * fh;

    if (!fp)
    {
        return 0;
    }

    p = (struct BP_PROC *) fp->fh;
    assert (p);

    fh = p->fh;

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
        free_namelist (fp->var_namelist, fp->nvars);
        fp->var_namelist = 0;
    }

    if (fp->attr_namelist)
    {
        free_namelist (fp->attr_namelist, fp->nattrs);
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

int bp_close (BP_FILE * fh)
{
    struct BP_GROUP_VAR * gh = fh->gvar_h;
    struct BP_GROUP_ATTR * ah = fh->gattr_h;
    struct adios_index_var_struct_v1 * vars_root = fh->vars_root, *vr;
    struct adios_index_attribute_struct_v1 * attrs_root = fh->attrs_root, *ar;
    struct bp_index_pg_struct_v1 * pgs_root = fh->pgs_root, *pr;
    int i,j;
    MPI_File mpi_fh = fh->mpi_fh;

    adios_errno = 0;
    if (fh->mpi_fh)
        MPI_File_close (&mpi_fh);

    if (fh->sfh)
        close_all_BP_files (fh->sfh);

    if (fh->b) {
        adios_posix_close_internal (fh->b);
        free(fh->b);
    }

    /* Free variable structures */
    /* alloc in bp_utils.c: bp_parse_vars() */
    while (vars_root) {
        vr = vars_root;
        vars_root = vars_root->next;
        for (j = 0; j < vr->characteristics_count; j++) {
            // alloc in bp_utils.c:bp_parse_characteristics() <- bp_get_characteristics_data()
            if (vr->characteristics[j].dims.dims)
                free (vr->characteristics[j].dims.dims);
            if (vr->characteristics[j].value)
                free (vr->characteristics[j].value);
            // NCSU - Clearing up statistics
            if (vr->characteristics[j].stats)
            {
                uint8_t k = 0, idx = 0;
                uint8_t i = 0, count = adios_get_stat_set_count(vr->type);

                while (vr->characteristics[j].bitmap >> k)
                {
                    if ((vr->characteristics[j].bitmap >> k) & 1)
                    {
                        for (i = 0; i < count; i ++)
                        {
                            if (k == adios_statistic_hist)
                            {
                                struct adios_index_characteristics_hist_struct * hist = (struct adios_index_characteristics_hist_struct *) (vr->characteristics [j].stats[i][idx].data);
                                free (hist->breaks);
                                free (hist->frequencies);
                                free (hist);
                            }
                            else
                            free (vr->characteristics[j].stats [i][idx].data);
                        }
                        idx ++;
                    }
                    k ++;
                }

                for (i = 0; i < count; i ++)
                    free (vr->characteristics[j].stats [i]);

                free (vr->characteristics[j].stats);
                vr->characteristics[j].stats = 0;
            }

            // NCSU ALACRITY-ADIOS - Clear transform metadata
            adios_transform_clear_transform_characteristic(&vr->characteristics[j].transform);
        }
        if (vr->characteristics)
            free (vr->characteristics);
        if (vr->group_name)
            free (vr->group_name);
        if (vr->var_name)
            free (vr->var_name);
        if (vr->var_path)
            free (vr->var_path);
        free(vr);
    }

    fh->vars_root = 0;

    /* Free attributes structures */
    /* alloc in bp_utils.c bp_parse_attrs() */
    while (attrs_root) {
        ar = attrs_root;
        attrs_root = attrs_root->next;
        for (j = 0; j < ar->characteristics_count; j++) {
            if (ar->characteristics[j].value)
                free (ar->characteristics[j].value);
        }
        if (ar->characteristics)
            free (ar->characteristics);
        if (ar->group_name)
            free (ar->group_name);
        if (ar->attr_name)
            free (ar->attr_name);
        if (ar->attr_path)
            free (ar->attr_path);
        free(ar);
    }

    fh->attrs_root = 0;

    /* Free process group structures */
    /* alloc in bp_utils.c bp_parse_pgs() first loop */
    //printf ("pgs: %d\n", fh->mfooter.pgs_count);
    while (pgs_root) {
        pr = pgs_root;
        pgs_root = pgs_root->next;
        //printf("%d\tpg pid=%d addr=%x next=%x\n",i, pr->process_id, pr, pr->next);
        if (pr->group_name)
            free(pr->group_name);
        if (pr->time_index_name)
            free(pr->time_index_name);
        free(pr);
    }

    fh->pgs_root = 0;

    /* Free variable structures in BP_GROUP_VAR */
    if (gh) {
        for (j=0;j<2;j++) {
            for (i=0;i<gh->group_count;i++) {
                if (gh->time_index[j][i])
                    free(gh->time_index[j][i]);
            }
            if (gh->time_index[j])
                free(gh->time_index[j]);
        }
        free (gh->time_index);

        for (i=0;i<gh->group_count;i++) {
            if (gh->namelist[i])
                free(gh->namelist[i]);
        }
        if (gh->namelist)
            free (gh->namelist);

        for (i=0;i<fh->mfooter.vars_count;i++) {
            if (gh->var_namelist[i])
                free(gh->var_namelist[i]);
            if (gh->var_offsets[i])
                free(gh->var_offsets[i]);
        }
        if (gh->var_namelist)
            free (gh->var_namelist);

        if (gh->var_offsets)
            free(gh->var_offsets);

        if (gh->var_counts_per_group)
            free(gh->var_counts_per_group);

        if (gh->pg_offsets)
            free (gh->pg_offsets);

        free (gh);
    }

    fh->gvar_h = 0;

    /* Free attribute structures in BP_GROUP_ATTR */
    if (ah) {
        for (i = 0; i < fh->mfooter.attrs_count; i++) {
            if (ah->attr_offsets[i])
                free(ah->attr_offsets[i]);
            if (ah->attr_namelist[i])
                free(ah->attr_namelist[i]);
        }
        if (ah->attr_offsets)
            free(ah->attr_offsets);
        if (ah->attr_namelist)
            free(ah->attr_namelist);
        if (ah->attr_counts_per_group)
            free(ah->attr_counts_per_group);

        free(ah);
    }

    fh->gattr_h = 0;

    if (fh->fname)
    {
        free (fh->fname);
        fh->fname = 0;
    }

    if (fh)
        free (fh);

    return 0;
}

/* Since we have no way to know that the end of the stream has been reached, we cannot
 * block the call and expect new step will arrive. Therefore, if last == 0, we simply sleep
 * for a specified period of time and re-open the file to see if there is any new steps came in.
 * The trick is that when step is being advanced, it is likely that file has
 * already being appended with new steps. Therefore, we have to close and reopen
 * the file if the expected step is not found.
 * last - 0: next available step, !=0: newest available step
 *  RETURN: 0 OK, !=0 on error (also sets adios_errno)
 *
 *  Possible errors (adios_errno values):
 *       err_end_of_stream    Stream has ended, no more steps should be expected
 *       err_step_notready    The requested step is not yet available
 *       err_step_disappeared The requested step is not available anymore

 */
int adios_read_bp_advance_step (ADIOS_FILE * fp, int last, float timeout_sec)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    int last_step;
    MPI_Comm comm;
    char * fname;

    log_debug ("adios_read_bp_advance_step\n");

    //TODO: this part of code needs to cleaned up a bit. Some if-else branches can be merged. Q.Liu
    adios_errno = 0;
    if (last == 0) // read in the next step
    {
        if (fp->current_step < fp->last_step) // no need to re-open file. The next step is already in.
        {
            release_step (fp);
            bp_seek_to_step (fp, ++fp->current_step, show_hidden_attrs);
        }
        else // re-open to read in footer again. We should keep polling until there are new steps in OR
             // time out.
        {
            last_step = fp->last_step;
            fname = strdup (fh->fname);
            comm = fh->comm;

            if (p->fh)
            {
                bp_close (fh);
                p->fh = 0;
            }

            if (!get_new_step (fp, fname, comm, last_step, timeout_sec))
            {
                // With file reading, how can we tell it is the end of the streams?
                adios_errno = err_step_notready;
            }

            free (fname);

            log_debug ("Seek from step %d to step %d\n", last_step, last_step + 1);

            if (adios_errno == 0)
            {
                release_step (fp);
                bp_seek_to_step (fp, last_step + 1, show_hidden_attrs);
                fp->current_step = last_step + 1;
            }
        }
    }
    else // read in newest step. Re-open no matter whether current_step < last_step
    {
        last_step = fp->last_step;
        fname = strdup (fh->fname);
        comm = fh->comm;

        if (p->fh)
        {
            bp_close (fh);
            p->fh = 0;
        }

        // lockmode is currently not supported.
        if (!get_new_step (fp, fh->fname, fh->comm, last_step, timeout_sec))
        {
            adios_errno = err_step_notready;
        }

        free (fname);

        if (adios_errno == 0)
        {
            release_step (fp);
            bp_seek_to_step (fp, fp->last_step, show_hidden_attrs);
            fp->current_step = fp->last_step;
        }
    }

    return adios_errno;
}

/* Right now, this function does nothing, since locking hasn't been implemented yet */
void adios_read_bp_release_step (ADIOS_FILE *fp)
{
}

ADIOS_VARINFO * adios_read_bp_inq_var_byid (const ADIOS_FILE * fp, int varid)
{
    struct BP_PROC * p;
    BP_FILE * fh;
    ADIOS_VARINFO * varinfo;
    int file_is_fortran, i, k, timedim, size;
    struct adios_index_var_struct_v1 * v;

    adios_errno = 0;

    p = (struct BP_PROC *) fp->fh;
    fh = (BP_FILE *)p->fh;

    if (varid < 0 || varid >= fp->nvars)
    {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)\n", varid, fp->nvars);
        return 0;
    }

    v = bp_find_var_byid (fh, varid);

    varinfo = (ADIOS_VARINFO *) malloc (sizeof (ADIOS_VARINFO));
    assert (varinfo);

    /* Note: set varid as the real varid.
       Common read layer should convert it to the perceived id after the read method returns.
    */
    varinfo->varid = varid;
    // NCSU ALACRITY-ADIOS - Use transformed datatype if a transform is applied
    varinfo->type = v->type;
    //varinfo->type = adios_transform_get_var_original_type(v); // LAYERFIX

    file_is_fortran = is_fortran_file (fh);

    assert (v->characteristics_count);

    bp_get_and_swap_dimensions (fh, v, file_is_fortran,
                                &varinfo->ndim, &varinfo->dims,
                                &varinfo->nsteps,
                                file_is_fortran != futils_is_called_from_fortran()
                               );

    if (p->streaming)
    {
        varinfo->nsteps = 1;
    }

    // set value for scalar
    if (v->characteristics [0].value)
    {
        size = bp_get_type_size (v->type, v->characteristics [0].value);
        varinfo->value = (void *) malloc (size);
        assert (varinfo->value);

        memcpy (varinfo->value, v->characteristics [0].value, size);
    }
    else
    {
        varinfo->value = NULL;
    }

    varinfo->global = is_global_array (&(v->characteristics[0]));
    varinfo->nblocks = get_var_nblocks (v, varinfo->nsteps);
    assert (varinfo->nblocks);

    varinfo->sum_nblocks = v->characteristics_count;
    varinfo->statistics = 0;
    varinfo->blockinfo = 0;

    return varinfo;
}

/* Note: most of the code of this routine is copied from previous 1.3.1 implementation
 * of NCSU. One thing to fix is that any counter related to characteristic index should be of uint64_t,
 * instead of int.
 */
int adios_read_bp_inq_var_stat (const ADIOS_FILE *fp, ADIOS_VARINFO * varinfo, int per_step_stat, int per_block_stat)
{
#if 0
typedef struct {
        void     * min;            /* minimum value in an array variable, = value for a scalar       */
        void     * max;            /* maximum value of an array variable (over all steps)            */
        double   * avg;            /* average value of an array variable (over all steps)            */
        double   * std_dev;        /* standard deviation value of an array variable (over all steps) */

        struct ADIOS_STAT_STEP     /* per step statistics (if requested and recorded at writing) */
        {
            void     ** mins;      /* minimum per each step (array of 'nsteps' elements)             */
            void     ** maxs;      /* maximum per each step (array of 'nsteps' elements)             */
            double   ** avgs;      /* average per each step (array of 'nsteps' elements)             */
            double   ** std_devs;  /* standard deviation per each step (array of 'nsteps' elements)  */
        } *steps;

        struct ADIOS_STAT_BLOCK    /* per block statistics (if requested and recorded at writing) */
        {
            void     ** mins;      /* minimum per each block (array of 'nblocks' elements)         */
            void     ** maxs;      /* maximum per each block (array of 'nblocks' elements)         */
            double   ** avgs;      /* average per each block (array of 'nblocks' elements)         */
            double   ** std_devs;  /* std deviation per each block (array of 'nblocks' elements)   */
        } *blocks;

        struct ADIOS_HIST           /* Histogram if recorded at writing */
        {
            uint32_t    num_breaks;
            double      max;
            double      min;
            double *    breaks;
            uint32_t ** frequencies;
            uint32_t *  gfrequencies;
        } *histogram;

} ADIOS_VARSTAT;
#endif
    int i, j, c, count = 1, timestep;
    int size, sum_size, sum_type, nsteps, prev_timestep;
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    ADIOS_VARSTAT * vs;
    struct adios_index_var_struct_v1 * var_root;

    assert (varinfo);

    varinfo->statistics = vs = (ADIOS_VARSTAT *) malloc (sizeof (ADIOS_VARSTAT));
    assert (vs);

    vs->min = 0;
    vs->max = 0;
    vs->avg = 0;
    vs->std_dev = 0;

    vs->steps = (struct ADIOS_STAT_STEP *) malloc (sizeof (struct ADIOS_STAT_STEP));
    assert (vs->steps);
    vs->steps->mins = 0;
    vs->steps->maxs = 0;
    vs->steps->avgs = 0;
    vs->steps->std_devs = 0;

    //TODO
    vs->blocks = 0;
    vs->histogram = 0;

    uint64_t gcnt = 0, * cnts;

    double *gsum = NULL, *gsum_square = NULL;
    double **sums = NULL, **sum_squares = NULL;

    int16_t map[32];
    memset (map, -1, sizeof(map));

    var_root = bp_find_var_byid (fh, varinfo->varid);

    // Bitmap shows which statistical information has been calculated
    i = j = 0;
    while (var_root->characteristics[0].bitmap >> j)
    {
        if ((var_root->characteristics[0].bitmap >> j) & 1)
        {
            map [j] = i ++;
        }

        j ++;
    }

    nsteps = varinfo->nsteps;

    if (map[adios_statistic_min] != -1)
    {
        MALLOC(vs->steps->mins, nsteps * sizeof(void *), "minimum per timestep");
        for (i = 0; i < nsteps; i++)
        {
            vs->steps->mins[i] = 0;
        }
    }

    if (map[adios_statistic_max] != -1)
    {
        MALLOC(vs->steps->maxs, nsteps * sizeof(void *), "maximum per timestep");
        for (i = 0; i < nsteps; i++)
        {
            vs->steps->maxs[i] = 0;
        }
    }

    if (map[adios_statistic_sum] != -1)
    {
        MALLOC(sums, nsteps * sizeof(double *), "summation per timestep");
        MALLOC(vs->steps->avgs, nsteps * sizeof(double *), "average per timestep");

        for (i = 0; i < nsteps; i++)
        {
            sums[i] = vs->steps->avgs[i] = 0;
        }
        CALLOC(cnts, nsteps, sizeof(uint64_t), "count of elements per timestep");
    }

    if (map[adios_statistic_sum_square] != -1)
    {
        MALLOC(sum_squares, nsteps * sizeof(double *), "summation per timestep");
        MALLOC(vs->steps->std_devs, nsteps * sizeof(double *), "standard deviation per timestep");

        for (i = 0; i < nsteps; i++)
        {
            vs->steps->std_devs[i] = sum_squares[i] = 0;
        }
    }
/*
    if (map[adios_statistic_hist] != -1 && (var_root->characteristics[0].stats[0][map[adios_statistic_hist]].data))
    {
        struct adios_index_characteristics_stat_struct * stats = var_root->characteristics[0].stats[0];
        struct adios_index_characteristics_hist_struct * hist = stats[map[adios_statistic_hist]].data;
        int num_breaks = hist->num_breaks;

        MALLOC(vi->hist, sizeof(struct ADIOS_HIST), "histogram");
        MALLOC(vi->hist->breaks, num_breaks * sizeof(double), "break points of histogram");
        MALLOC(vi->hist->gfrequencies, (num_breaks + 1) * sizeof(uint32_t), "global frequencies of histogram");

        vi->hist->num_breaks = hist->num_breaks;
        vi->hist->min = hist->min;
        vi->hist->max = hist->max;

        memcpy(vi->hist->breaks, hist->breaks, num_breaks * sizeof(double));
        CALLOC(vi->hist->gfrequencies, (num_breaks + 1), bp_get_type_size(adios_unsigned_integer, ""), "global frequency");

        if (ntimes > 0)
        {
            MALLOC(vi->hist->frequenciess, (ntimes * sizeof(int32_t *)), "frequencies for timesteps");
            for(i = 0; i < ntimes; i++)
                CALLOC(vi->hist->frequenciess[i], (num_breaks + 1), bp_get_type_size(adios_unsigned_integer, ""), "frequency at timestep");
        }
    }
*/
    size = bp_get_type_size (var_root->type, "");
    sum_size = bp_get_type_size (adios_double, "");

    if (var_root->type == adios_complex || var_root->type == adios_double_complex)
    {
        int type;
        count = 3;
        timestep = -1;
        prev_timestep = 0;

        if (var_root->type == adios_complex)
        {
            type = adios_double;
        }
        else
        {
            type = adios_long_double;
        }

        // Only a double precision returned for all complex values
        size = bp_get_type_size (adios_double, "");

        for (i = 0; i < var_root->characteristics_count; i++)
        {
            // changes for 1.4.x. Q. Liu
            if (var_root->characteristics[i].time_index != prev_timestep)
            {
                timestep++;
                prev_timestep = var_root->characteristics[i].time_index;
            }

            assert (timestep < nsteps);

            if (!var_root->characteristics[i].stats)
                continue;

            struct adios_index_characteristics_stat_struct ** stats = var_root->characteristics[i].stats;

            if ((map[adios_statistic_finite] != -1) && (* ((uint8_t *) stats[0][map[adios_statistic_finite]].data) == 0))
                continue;

            if (map[adios_statistic_min] != -1 && stats[0][map[adios_statistic_min]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double((enum ADIOS_DATATYPES)type, stats[c][map[adios_statistic_min]].data);

                if(!vs->min)
                {
                    MALLOC (vs->min, count * size, "global minimum")
                    for (c = 0; c < count; c ++)
                           ((double * ) vs->min)[c] = data[c];

                }
                else
                {
                    for (c = 0; c < count; c ++)
                        if (data[c] < ((double *) vs->min)[c])
                               ((double * ) vs->min)[c] = data[c];
                }

                if(!vs->steps->mins[timestep])
                {
                    MALLOC (vs->steps->mins[timestep], count * size, "minimum per timestep")
                    for (c = 0; c < count; c ++)
                    {
                        ((double **) vs->steps->mins)[timestep][c] = data[c];
                    }
                }
                else
                {
                    for (c = 0; c < count; c ++)
                    {
                        if (data[c] < ((double **) vs->steps->mins)[timestep][c])
                        {
                            ((double **) vs->steps->mins)[timestep][c] = data[c];
                        }
                    }
                }
            }

            if (map[adios_statistic_max] != -1 && stats[0][map[adios_statistic_max]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double((enum ADIOS_DATATYPES)type, stats[c][map[adios_statistic_max]].data);

                if(!vs->max) {
                    MALLOC (vs->max, count * size, "global minimum")
                    for (c = 0; c < count; c ++)
                        ((double * ) vs->max)[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        if (data[c] > ((double *) vs->max)[c])
                            ((double * ) vs->max)[c] = data[c];
                }

                if(!vs->steps->maxs[timestep]) {
                    MALLOC (vs->steps->maxs[timestep], count * size, "minimum per timestep")
                    for (c = 0; c < count; c ++)
                        ((double **) vs->steps->maxs)[timestep][c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        if (data[c] > ((double **) vs->steps->maxs)[timestep][c])
                            ((double **) vs->steps->maxs)[timestep][c] = data[c];
                }
            }

            if (map[adios_statistic_sum] != -1 && stats[0][map[adios_statistic_sum]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double((enum ADIOS_DATATYPES)type, stats[c][map[adios_statistic_sum]].data);

                if(!gsum) {
                    MALLOC(gsum, count * sum_size, "global summation")
                    for (c = 0; c < count; c ++)
                           gsum[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        gsum[c] = gsum[c] + data[c];
                }

                if(!sums[timestep]) {
                    MALLOC(sums[timestep], count * sum_size, "summation per timestep")
                    for (c = 0; c < count; c ++)
                        sums[timestep][c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        sums[timestep][c] = sums[timestep][c] + data[c];
                }
            }

            if (map[adios_statistic_sum_square] != -1 && stats[0][map[adios_statistic_sum_square]].data)
            {
                double data[3];
                for (c = 0; c < count; c ++)
                    data[c] = bp_value_to_double((enum ADIOS_DATATYPES)type, stats[c][map[adios_statistic_sum_square]].data);

                if(!gsum_square) {
                    MALLOC(gsum_square, count * sum_size, "global summation of squares")
                    for (c = 0; c < count; c ++)
                        gsum_square[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                           gsum_square[c] = gsum_square[c] + data[c];
                }

                if(!sum_squares[timestep]) {
                    MALLOC(sum_squares[timestep], count * sum_size, "summation of square per timestep")
                    for (c = 0; c < count; c ++)
                        sum_squares[timestep][c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        sum_squares[timestep][c] = sum_squares[timestep][c] + data[c];
                }
            }

            if (map[adios_statistic_cnt] != -1 && stats[0][map[adios_statistic_cnt]].data)
            {
                cnts[timestep] += * ((uint32_t *) stats[0][map[adios_statistic_cnt]].data);
                gcnt += * (uint32_t *) stats[0][map[adios_statistic_cnt]].data;
            }
        }

        if(vs->min && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            // min, max, summation exists only for arrays
            // Calculate average / timestep

            for(timestep = 0; timestep < nsteps; timestep ++) {
                MALLOC(vs->steps->avgs[timestep], count * sum_size, "average per timestep")
                for (c = 0; c < count; c ++)
                    vs->steps->avgs[timestep][c] = sums[timestep][c] / cnts[timestep];

                MALLOC(vs->steps->std_devs[timestep], count * sum_size, "standard deviation per timestep")
                for (c = 0; c < count; c ++)
                    vs->steps->std_devs[timestep][c] = sqrt((sum_squares[timestep][c] / cnts[timestep]) - (vs->steps->avgs[timestep][c] * vs->steps->avgs[timestep][c]));

                free (sums[timestep]);
                free (sum_squares[timestep]);
            }
        }

        // Calculate global average
        if(vs->min && gsum && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
            MALLOC(vs->avg, count * sum_size, "global average")

            if(gcnt > 0)
                for (c = 0; c < count; c ++)
                    vs->avg[c] = gsum[c] / gcnt;
            else
                for (c = 0; c < count; c ++)
                    vs->avg[c] = gsum[c];

            MALLOC(vs->std_dev, count * sum_size, "global average")
            if(vs->avg && gcnt > 0)
                for (c = 0; c < count; c ++)
                    vs->std_dev[c] = sqrt(gsum_square[c] / gcnt - (vs->avg[c] * vs->avg[c]));
            else
                for (c = 0; c < count; c ++)
                    vs->std_dev[c] = 0;
        }
    }
    else
    {
        timestep = -1;
        prev_timestep = 0;

        for (i = 0; i < var_root->characteristics_count; i++)
        {
            // changes for 1.4.x. Q. Liu
            if (var_root->characteristics[i].time_index != prev_timestep)
            {
                timestep++;
                prev_timestep = var_root->characteristics[i].time_index;
            }

            assert (timestep < nsteps);

            if (!var_root->characteristics[i].stats)
            {
                continue;
            }

            struct adios_index_characteristics_stat_struct * stats = var_root->characteristics[i].stats[0];
            struct adios_index_characteristics_hist_struct * hist = stats[map[adios_statistic_hist]].data;

            if (map[adios_statistic_finite] != -1 && (* ((uint8_t *) stats[map[adios_statistic_finite]].data) == 0))
                continue;

            if (map[adios_statistic_min] != -1 && stats[map[adios_statistic_min]].data)
            {
                if(!vs->min)
                {
                    MALLOC (vs->min, size, "global minimum")
                    memcpy(vs->min, stats[map[adios_statistic_min]].data, size);

                }
                else if (adios_lt(var_root->type, stats[map[adios_statistic_min]].data, vs->min))
                {
                    memcpy(vs->min, stats[map[adios_statistic_min]].data, size);
                }

                if(!vs->steps->mins[timestep])
                {
                    MALLOC (vs->steps->mins[timestep], size, "minimum per timestep")
                    memcpy(vs->steps->mins[timestep], stats[map[adios_statistic_min]].data, size);
                }
                else if (adios_lt(var_root->type, stats[map[adios_statistic_min]].data, vs->steps->mins[timestep]))
                {
                    memcpy(vs->steps->mins[timestep], stats[map[adios_statistic_min]].data, size);
                }
            }

            if (map[adios_statistic_max] != -1 && stats[map[adios_statistic_max]].data)
            {
                if(!vs->max)
                {
                    MALLOC (vs->max, size, "global maximum")
                    memcpy(vs->max, stats[map[adios_statistic_max]].data, size);

                }
                else if (adios_lt(var_root->type, vs->max, stats[map[adios_statistic_max]].data))
                {
                    memcpy(vs->max, stats[map[adios_statistic_max]].data, size);
                }

                if(!vs->steps->maxs[timestep])
                {
                    MALLOC (vs->steps->maxs[timestep], size, "maximum per timestep")
                    memcpy(vs->steps->maxs[timestep], stats[map[adios_statistic_max]].data, size);
                }
                else if (adios_lt(var_root->type, vs->steps->maxs[timestep], stats[map[adios_statistic_max]].data))
                {
                    memcpy(vs->steps->maxs[timestep], stats[map[adios_statistic_max]].data, size);
                }
            }

            if (map[adios_statistic_sum] != -1 && stats[map[adios_statistic_sum]].data)
            {
                if(!gsum)
                {
                    MALLOC(gsum, sum_size, "global summation")
                    memcpy(gsum, stats[map[adios_statistic_sum]].data, sum_size);
                }
                else
                {
                    *gsum = *gsum + * ((double *) stats[map[adios_statistic_sum]].data);
                }

                if(!sums[timestep])
                {
                    MALLOC(sums[timestep], sum_size, "summation per timestep")
                    memcpy(sums[timestep], stats[map[adios_statistic_sum]].data, sum_size);
                }
                else
                {
                    *sums[timestep] = *sums[timestep] + * ((double *) stats[map[adios_statistic_sum]].data);
                }
            }

            if (map[adios_statistic_sum_square] != -1 && stats[map[adios_statistic_sum_square]].data)
            {
                if(!gsum_square)
                {
                    MALLOC(gsum_square, sum_size, "global summation of squares")
                    memcpy(gsum_square, stats[map[adios_statistic_sum_square]].data, sum_size);

                }
                else
                {
                    *gsum_square = *gsum_square + * ((double *) stats[map[adios_statistic_sum_square]].data);
                }

                if(!sum_squares[timestep])
                {
                    MALLOC(sum_squares[timestep], sum_size, "summation of square per timestep")
                    memcpy(sum_squares[timestep], stats[map[adios_statistic_sum_square]].data, sum_size);
                }
                else
                {
                    *sum_squares[timestep] = *sum_squares[timestep] + * ((double *) stats[map[adios_statistic_sum_square]].data);
                }
            }
//TODO
/*
            if(map[adios_statistic_hist] != -1 && stats[map[adios_statistic_hist]].data)
            {
                for(j = 0; j <= vi->hist->num_breaks; j++)
                {
                    uint32_t freq = hist->frequencies[j];
                    vi->hist->gfrequencies[j] += freq;
                    if (ntimes > 0)
                        vi->hist->frequenciess[timestep][j] += freq;
                }
            }
*/
            if (map[adios_statistic_cnt] != -1 && stats[map[adios_statistic_cnt]].data)
            {
                cnts[timestep] += * (uint32_t *) stats[map[adios_statistic_cnt]].data;
                gcnt += * (uint32_t *) stats[map[adios_statistic_cnt]].data;
            }
        }

//TODO
#if 0
        if(ntimes > 0 && vs->min && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1))
        {
            // min, max, summation exists only for arrays
            // Calculate average / timestep

            for(timestep = 0; timestep < ntimes; timestep ++)
            {
                MALLOC(vs->steps->avgs[timestep], sum_size, "average per timestep")
                *(vs->steps->avgs[timestep]) = *(sums[timestep]) / cnts[timestep];

                MALLOC(vs->steps->std_devs[timestep], sum_size, "standard deviation per timestep")
                *(vs->steps->std_devs[timestep]) = sqrt(*(sum_squares[timestep]) / cnts[timestep] - ((*(vs->steps->avgs[timestep]) * (*(vs->steps->avgs[timestep])))));

                free (sums[timestep]);
                free (sum_squares[timestep]);
            }
        }
#endif
        // Calculate global average
        if(vs->min && gsum && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1))
        {
            MALLOC(vs->avg, sum_size, "global average")
            if(gcnt > 0)
                *vs->avg = *gsum / gcnt;
            else
                vs->avg = gsum;

            MALLOC(vs->std_dev, sum_size, "global average")
            if(vs->avg && gcnt > 0)
                *vs->std_dev = sqrt(*gsum_square / gcnt - ((*(vs->avg)) * (*(vs->avg))));
            else
                *vs->std_dev = 0;
        }
    }

    if (!varinfo->value && vs->min)
    {
        varinfo->value = vs->min; // arrays have no value but we assign here the minimum
    }

    if(!vs->min)
    {
        vs->min = varinfo->value; // scalars have value but not min
    }

    if(!vs->max)
    {
        vs->max = varinfo->value; // scalars have value but not max
    }

    if (sums && gsum)
    {
        free (sums);
        free (gsum);
    }

    if (sum_squares && gsum_square)
    {
        free (sum_squares);
        free (gsum_square);
    }

    return 0;
}

// NCSU ALACRITY-ADIOS - Factored out VARBLOCK inquiry function to permit sourcing
static ADIOS_VARBLOCK * inq_var_blockinfo(const ADIOS_FILE * fp, const ADIOS_VARINFO * varinfo, int use_pretransform_dimensions) {
    struct BP_PROC * p = (struct BP_PROC *) fp->fh;
    int i, file_is_fortran, timedim;
    uint64_t * ldims, * gdims, * offsets;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * var_root;
    ADIOS_VARBLOCK *blockinfo;

    assert (varinfo);

    fh = (BP_FILE *) p->fh;
    file_is_fortran = is_fortran_file (fh);
    var_root = bp_find_var_byid (fh, varinfo->varid);

    if (use_pretransform_dimensions)
        assert(var_root->characteristics[0].transform.transform_type != adios_transform_none);

    blockinfo = (ADIOS_VARBLOCK *) malloc (varinfo->sum_nblocks * sizeof (ADIOS_VARBLOCK));
    assert(blockinfo);

    // NCSU ALACRITY-ADIOS - Use pre-transform dimensions if instructed to do so
    int dimcount;
    if (use_pretransform_dimensions) {
        dimcount = var_root->characteristics[0].transform.pre_transform_dimensions.count;
    } else {
        dimcount = var_root->characteristics[0].dims.count;
    }
    /* dim.count possibily include 'time' dim in it. */
    ldims = (uint64_t *) malloc (dimcount * sizeof(uint64_t));
    gdims = (uint64_t *) malloc (dimcount * sizeof(uint64_t));
    offsets = (uint64_t *) malloc (dimcount * sizeof(uint64_t));
    assert (ldims && gdims && offsets);

    for (i = 0; i < varinfo->sum_nblocks; i++)
    {
        blockinfo[i].start = (uint64_t *) malloc (dimcount * sizeof(uint64_t));
        blockinfo[i].count = (uint64_t *) malloc (dimcount * sizeof(uint64_t));
        assert (blockinfo[i].start && blockinfo[i].count);

        // NCSU ALACRITY-ADIOS - Fixed a bug with inq_var_blockinfo that was
        //   triggered in the presence of a global array with time dimension
        bp_get_dimension_generic_notime(use_pretransform_dimensions ?
                                    &var_root->characteristics[i].transform.pre_transform_dimensions :
                                    &var_root->characteristics[i].dims,
                                 ldims, gdims, offsets, file_is_fortran
                                 );

        /*
        if (dimcount != varinfo->ndim) // has time
        {
            if (file_is_fortran) // For Fortran file, time must be the last dim
            {
                memcpy (blockinfo[i].start, offsets, varinfo->ndim * sizeof(uint64_t));
                memcpy (blockinfo[i].count, ldims, varinfo->ndim * sizeof(uint64_t));
            }
            else // For C file, time must be the first dim
            {
                memcpy (blockinfo[i].start, offsets + 1, varinfo->ndim * sizeof(uint64_t));
                memcpy (blockinfo[i].count, ldims + 1, varinfo->ndim * sizeof(uint64_t));
            }
        }
        else // no time
        {
            memcpy (blockinfo[i].start, offsets, varinfo->ndim * sizeof(uint64_t));
            memcpy (blockinfo[i].count, ldims, varinfo->ndim * sizeof(uint64_t));
        }*/
        memcpy (blockinfo[i].start, offsets, varinfo->ndim * sizeof(uint64_t));
        memcpy (blockinfo[i].count, ldims, varinfo->ndim * sizeof(uint64_t));

        if (file_is_fortran != futils_is_called_from_fortran())
        {
            swap_order (varinfo->ndim, blockinfo[i].start, &timedim);
            swap_order (varinfo->ndim, blockinfo[i].count, &timedim);
        }
    }

    free (ldims);
    free (gdims);
    free (offsets);

    return blockinfo;
}

int adios_read_bp_inq_var_blockinfo (const ADIOS_FILE * fp, ADIOS_VARINFO * varinfo)
{
    // NCSU ALACRITY-ADIOS - Delegate to shared VARBLOCK loader
    varinfo->blockinfo = inq_var_blockinfo(fp, varinfo, 0); // 0 -> use true dimensions, not original dimensions
    assert(varinfo->blockinfo);
    return 0;
/*
    struct BP_PROC * p = (struct BP_PROC *) fp->fh;
    int i, file_is_fortran, timedim;
    uint64_t * ldims, * gdims, * offsets;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * var_root;

    assert (varinfo);

    fh = (BP_FILE *) p->fh;
    file_is_fortran = is_fortran_file (fh);
    var_root = bp_find_var_byid (fh, varinfo->varid);

    varinfo->blockinfo = (ADIOS_VARBLOCK *) malloc (varinfo->sum_nblocks * sizeof (ADIOS_VARBLOCK));
    assert (varinfo->blockinfo);

    // NCSU ALACRITY-ADIOS - Use the pre-transform dimensions when applicable
    int dimcount = var_root->characteristics[0].dims.count;//adios_transform_get_characteristic_original_num_dims(&var_root->characteristics[0]); // LAYERFIX
    // dim.count possibily include 'time' dim in it.
    ldims = (uint64_t *) malloc (dimcount);
    gdims = (uint64_t *) malloc (dimcount);
    offsets = (uint64_t *) malloc (dimcount);
    assert (ldims && gdims && offsets);

    for (i = 0; i < varinfo->sum_nblocks; i++)
    {
        varinfo->blockinfo[i].start = (uint64_t *) malloc (varinfo->ndim * 8);
        varinfo->blockinfo[i].count = (uint64_t *) malloc (varinfo->ndim * 8);
        assert (varinfo->blockinfo[i].start && varinfo->blockinfo[i].count);

        bp_get_dimension_characteristics (&(var_root->characteristics[i]),
                                          ldims, gdims, offsets
                                         );

        if (var_root->characteristics[0].dims.count != varinfo->ndim) // has time
        {
            if (file_is_fortran) // For Fortran file, time must be the last dim
            {
                memcpy (varinfo->blockinfo[i].start, offsets, varinfo->ndim * 8);
                memcpy (varinfo->blockinfo[i].count, ldims, varinfo->ndim * 8);
            }
            else // For C file, time must be the first dim
            {
                memcpy (varinfo->blockinfo[i].start, offsets + 1, varinfo->ndim * 8);
                memcpy (varinfo->blockinfo[i].count, ldims + 1, varinfo->ndim * 8);
            }
        }
        else // no time
        {
            memcpy (varinfo->blockinfo[i].start, offsets, varinfo->ndim * 8);
            memcpy (varinfo->blockinfo[i].count, ldims, varinfo->ndim * 8);
        }

        if (file_is_fortran != futils_is_called_from_fortran())
        {
            swap_order (varinfo->ndim, varinfo->blockinfo[i].start, &timedim);
            swap_order (varinfo->ndim, varinfo->blockinfo[i].count, &timedim);
        }
    }

    free (ldims);
    free (gdims);
    free (offsets);

    return 0;
*/
}

ADIOS_TRANSINFO * adios_read_bp_inq_var_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi) {
    struct BP_PROC * p = (struct BP_PROC *) fp->fh;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * var_root;
    int file_is_fortran;
    int dummy;
    ADIOS_TRANSINFO *transinfo;

    assert(vi);
    fh = (BP_FILE *) p->fh;
    file_is_fortran = is_fortran_file (fh);
    var_root = bp_find_var_byid(fh, vi->varid);
    assert(var_root);

    transinfo = malloc(sizeof(ADIOS_TRANSINFO));

    const struct adios_index_characteristic_transform_struct *transform = &var_root->characteristics[0].transform;

    transinfo->transform_type = transform->transform_type;
    if (transform->transform_type != adios_transform_none) {
        transinfo->orig_type = transform->pre_transform_type;

        // Load orig_ndims/orig_dims using the utility function
        bp_get_and_swap_dimensions_generic (fh, var_root, file_is_fortran,
                                            &transinfo->orig_ndim, &transinfo->orig_dims,
                                            &dummy,
                                            file_is_fortran != futils_is_called_from_fortran(),
                                            1); // 1 -> get based on pre-transform dimensions

        transinfo->orig_global = is_global_array_generic(&var_root->characteristics[0].transform.pre_transform_dimensions);
    } else {
        transinfo->orig_type = adios_unknown;
        transinfo->orig_ndim = 0;
        transinfo->orig_dims = 0;
        transinfo->orig_global = 0;
    }
    transinfo->orig_blockinfo = 0;

    return transinfo;
}
ADIOS_TRANSINFO * adios_read_bp_staged_inq_var_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi) {
    return NULL;
}
ADIOS_TRANSINFO * adios_read_bp_staged1_inq_var_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi) {
    return NULL;
}

int adios_read_bp_inq_var_trans_blockinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    struct BP_PROC * p = (struct BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    struct adios_index_var_struct_v1 * var_root;

    assert(vi);
    assert(ti);
    var_root = bp_find_var_byid(fh, vi->varid);
    assert(var_root);

    ti->orig_blockinfo = inq_var_blockinfo(fp, vi, 1); // 1 -> use original, pretransform dimensions
    return 0;
}

int adios_read_bp_staged_inq_var_trans_blockinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    return 1;
}
int adios_read_bp_staged1_inq_var_trans_blockinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    return 1;
}


/*
ADIOS_TRANSINFO * adios_read_bp_inq_var_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi) {
    struct BP_PROC * p = (struct BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    struct adios_index_var_struct_v1 * var_root;

    var_root = bp_find_var_byid(fh, vi->varid);

    ADIOS_TRANSINFO *retval = malloc(sizeof(ADIOS_TRANSINFO));

    struct adios_index_characteristic_transform_struct *transform = &var_root->characteristics[0]->transform;

    retval->transform_type = transform->transform_type;
    retval->

    //adios_characteristic_transform_info *transform_info = malloc(sizeof(adios_characteristic_transform_info));
    transform_info->transform_type = var_root->characteristics[0].transform_type;
    transform_info->pre_transform_type = var_root->characteristics[0].pre_transform_type;
    transform_info->pre_transform_dimensions.count = var_root->characteristics[0].pre_transform_dimensions.count;
}*/


/* Note: the varid isn't the perceived varid from the user */
/** Schedule reading a variable (slice) from the file.
 *  You need to allocate the memory for the data.
 *  You need to call adios_perform_reads() to do the reading of
 *  variables.
 *  IN:  fp         pointer to an (opened) ADIOS_FILE struct
 *       sel        selection created beforehand with adios_selection...().
 *                  sel=NULL means global selection (whole variable)
 *       varname    name of the variable
 *       from_step  Read the 'nsteps' consecutive steps from this
 *                  step of a file variable.
                    It is not used in case of a stream.
 *       nsteps     Read 'nsteps' consecutive steps from current step.
 *                  Must be 1 for a stream.
 *  OUT: data       pointer to the memory to hold data of the variable
 *                  In blocking read mode, the memory should be
 *                  pre-allocated. In non-blocking mode, memory can be
 *                  allocated or not, and that changes the behavior of
 *                  the chunked read. If memory is allocated,
 *                  adios_check_read() returns a variable if it is completed.
 *                  If memory is not allocated, the check returns any chunk
 *                  already available of a variable (in ADIOS own memory)
 *                  and the application has to rearrange the data. The user
 *                  has to process/copy the data before getting new chunks.
 *  RETURN: 0 OK, !=0 on error, sets adios_errno too
 */
int adios_read_bp_schedule_read_byid (const ADIOS_FILE * fp, const ADIOS_SELECTION * sel,
                                      int varid, int from_steps, int nsteps, void * data)
{
    BP_PROC * p;
    BP_FILE * fh;
    read_request * r;
    ADIOS_SELECTION * nullsel = 0;
    struct adios_index_var_struct_v1 * v;
    uint64_t datasize;
    int i, type_size, ndim, ns, file_is_fortran, mapped_varid;
    uint64_t * dims = 0;

    assert (fp);

    p = (BP_PROC *) fp->fh;
    fh = (BP_FILE *) p->fh;

    mapped_varid = p->varid_mapping[varid];
    v = bp_find_var_byid (fh, mapped_varid);
    file_is_fortran = is_fortran_file (fh);

    r = (read_request *) malloc (sizeof (read_request));
    assert (r);

    // If no selection was specified, generate a selection that includes the
    // entire variable (Drew)
    if (!sel)
    {
        bp_get_and_swap_dimensions (fh, v, file_is_fortran,
                                    &ndim, &dims,
                                    &ns,
                                    file_is_fortran != futils_is_called_from_fortran()
                                    );

        nullsel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
        assert (nullsel);

        nullsel->type == ADIOS_SELECTION_BOUNDINGBOX;
        nullsel->u.bb.ndim = ndim;
        nullsel->u.bb.start = (uint64_t *) malloc (nullsel->u.bb.ndim * 8);
        assert (nullsel->u.bb.start);
        nullsel->u.bb.count = (uint64_t *) malloc (nullsel->u.bb.ndim * 8);
        assert (nullsel->u.bb.count);

        for (i = 0; i < nullsel->u.bb.ndim; i++)
        {
            nullsel->u.bb.start[i] = 0;
            nullsel->u.bb.count[i] = dims[i];
        }

        free (dims);
    }

    r->rank = p->rank;
    /* copy selection since we don't want to operate on user memory.
     */
    r->sel = (!nullsel ? copy_selection (sel) : nullsel);
    r->varid = mapped_varid;
    if (!p->streaming)
    {
        r->from_steps = from_steps;
        r->nsteps = nsteps;
    }
    else
    {
        r->from_steps = 0;
        r->nsteps = 1;
    }
    r->data = data;
    //FIXME
    type_size = bp_get_type_size (v->type, 0);;

    if (sel->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        datasize = type_size;
        for (i = 0; i < sel->u.bb.ndim; i++)
        {
            datasize *=  sel->u.bb.count[i];
        }
    }
    else if (sel->type == ADIOS_SELECTION_POINTS)
    {
        datasize = type_size * sel->u.points.npoints;
    }

    r->priv = 0;
    r->next = 0;

    list_insert_read_request_next (&p->local_read_request_list, r);

    return 0;
}

int adios_read_bp_perform_reads (const ADIOS_FILE *fp, int blocking)
{
    struct BP_PROC * p;
    read_request * r;
    ADIOS_VARCHUNK * chunk;

    p = (struct BP_PROC *) fp->fh;

    // If blocking mode, ensure all read requests have a user memory buffer.
    // Otherwise, do nothing and return immediately (Drew)
    /* 1. prepare all reads */
    // check if all user memory is provided for blocking read
    if (blocking)
    {
        r = p->local_read_request_list;
        while (r)
        {
            if (!r->data)
            {
                adios_error (err_operation_not_supported,
                    "Blocking mode at adios_perform_reads() requires that user "
                    "provides the memory for each read request. Request for "
                    "variable %d was scheduled without user-allocated memory\n",
                    r->varid);
                return err_operation_not_supported;
            }

            r = r->next;
        }
    }
    else
    {
        return 0;
    }

    // Issue all read requests (Drew)
    while (p->local_read_request_list)
    {
        chunk = read_var (fp, p->local_read_request_list);

        // remove head from list
        r = p->local_read_request_list;
        p->local_read_request_list = p->local_read_request_list->next;
        free(r);

        common_read_free_chunk (chunk);
    }

    return 0;
}

int adios_read_bp_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk)
{
    BP_PROC * p;
    read_request * r;
    ADIOS_VARCHUNK * varchunk;
#if 0
 *  RETURN:         0: all chunks have been returned previously,
 *                     no need to call again (chunk is NULL, too)
 *                  1: some chunks are/will be available, call again
 *                  <0 on error, sets adios_errno too
#endif
    p = (struct BP_PROC *) fp->fh;

    if (p->local_read_request_list)
    {
        varchunk = read_var (fp, p->local_read_request_list);

        if (varchunk)
        {
            // remove head from list
            r = p->local_read_request_list;
            p->local_read_request_list = p->local_read_request_list->next;
            free(r);

            * chunk = varchunk;
            return 1;
        }
        else
        {
            return adios_errno;
        }
    }

    return 0;
}

int adios_read_bp_get_attr_byid (const ADIOS_FILE * fp, int attrid, enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    int i, offset, count;
    struct BP_GROUP * gh;
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    struct adios_index_attribute_struct_v1 * attr_root;
    struct adios_index_var_struct_v1 * var_root;
    int file_is_fortran, last_step = fp->last_step;
    uint64_t k, attr_c_index, var_c_index;

    adios_errno = 0;

    attr_root = fh->attrs_root; /* need to traverse the attribute list of the group */
    for (i = 0; i < attrid && attr_root; i++)
    {
        attr_root = attr_root->next;
    }

    if (i != attrid)
    {
        adios_error (err_corrupted_attribute, "Attribute id=%d is valid but was not found in internal data structures!\n",attrid);
        return adios_errno;
    }

    /* Look for the last step because some of the hidden attributes, such as last update time,
     * make sense for the most recent value. 07/2011 - Q.Liu
     */

    attr_c_index = -1;
    for (k = 0; k < attr_root->characteristics_count; k++)
    {
        if (attr_root->characteristics[k].time_index - 1 == last_step)
        {
            attr_c_index = k;
            break;
        }
    }

    if (attr_c_index == -1)
    {
        log_debug ("adios_read_bp_get_attr_byid: cannot find step : %d\n", last_step);
        attr_c_index = 0;
    }

    file_is_fortran = is_fortran_file (fh);

    // check the last version
    if (attr_root->characteristics[attr_c_index].value)
    {
        /* Attribute has its own value */
        *size = bp_get_type_size (attr_root->type, attr_root->characteristics[attr_c_index].value);
        *type = attr_root->type;
        *data = (void *) malloc (*size);
        assert (*data);

        memcpy(*data, attr_root->characteristics[attr_c_index].value, *size);
    }
    else if (attr_root->characteristics[attr_c_index].var_id)
    {
        /* Attribute is a reference to a variable */
        /* FIXME: var ids are not unique in BP. If a group of variables are written several
           times under different path using adios_set_path(), the id of a variable is always
           the same (should be different). As a temporary fix, we look first for a matching
           id plus path between an attribute and a variable. If not found, then we look for
           a match on the ids only.*/
        var_root = fh->vars_root;
        while (var_root)
        {
            if (var_root->id == attr_root->characteristics[attr_c_index].var_id
               && !strcmp(var_root->var_path, attr_root->attr_path)
               )
                break;
            var_root = var_root->next;
        }
        if (!var_root)
        {
            var_root = fh->vars_root;
            while (var_root)
            {
                if (var_root->id == attr_root->characteristics[attr_c_index].var_id)
                    break;
                var_root = var_root->next;
            }
        }

        if (!var_root)
        {
            adios_error (err_invalid_attribute_reference,
                   "Attribute %s/%s in group %s is a reference to variable ID %d, which is not found\n",
                   attr_root->attr_path, attr_root->attr_name, attr_root->group_name,
                   attr_root->characteristics[attr_c_index].var_id);
            return adios_errno;
        }

        /* default values in case of error */
        *data = NULL;
        *size = 0;
        *type = attr_root->type;

        var_c_index = -1;
        for (k = 0; k < var_root->characteristics_count; k++)
        {
            if (var_root->characteristics[k].time_index - 1 == last_step)
            {
                var_c_index = k;
                break;
            }
        }

        if (var_c_index == -1)
        {
            var_c_index = 0;
            log_debug ("adios_read_bp_get_attr_byid: cannot find step : %d\n", last_step);
        }
        /* FIXME: variable and attribute type may not match, then a conversion is needed. */
        /* Cases:
                1. attr has no type, var is byte array     ==> string
                2. attr has no type, var is not byte array ==> var type
                3. attr is string, var is byte array       ==> string
                4. attr type == var type                   ==> var type
                5. attr type != var type                   ==> attr type and conversion needed
        */
        /* Error check: attr cannot reference an array in general */
        if (var_root->characteristics[var_c_index].dims.count > 0)
        {
            if ( (var_root->type == adios_byte || var_root->type == adios_unsigned_byte) &&
                 (attr_root->type == adios_unknown || attr_root->type == adios_string) &&
                 (var_root->characteristics[var_c_index].dims.count == 1))
            {
                 ; // this conversions are allowed
            }
            else
            {
                adios_error (err_invalid_attribute_reference,
                    "Attribute %s/%s in group %s, typeid=%d is a reference to an %d-dimensional array variable "
                    "%s/%s of type %s, which is not supported in ADIOS\n",
                    attr_root->attr_path, attr_root->attr_name, attr_root->group_name, attr_root->type,
                    var_root->characteristics[var_c_index].dims.count,
                    var_root->var_path, var_root->var_name, common_read_type_to_string(var_root->type));
                return adios_errno;
            }
        }

        if ( (attr_root->type == adios_unknown || attr_root->type == adios_string) &&
             (var_root->type == adios_byte || var_root->type == adios_unsigned_byte) &&
             (var_root->characteristics[var_c_index].dims.count == 1) )
        {
            /* 1D byte arrays are converted to string */
            /* 1. read in variable */
            char varname[512];
            char *tmpdata;
            ADIOS_VARCHUNK *vc;
            read_request * r;
            uint64_t start, count;
            int status;

            start = 0;
            count = var_root->characteristics[var_c_index].dims.dims[0];
            snprintf(varname, 512, "%s/%s", var_root->var_path, var_root->var_name);
            tmpdata = (char *) malloc (count+1);
            assert (tmpdata);

            r = (read_request *) malloc (sizeof (read_request));
            assert (r);

            r->rank = p->rank;
            r->sel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
            r->sel->type = ADIOS_SELECTION_BOUNDINGBOX;
            r->sel->u.bb.ndim = 1;
            r->sel->u.bb.start = &start;
            r->sel->u.bb.count = &count;
            r->varid = var_root->id;
            r->from_steps = fp->last_step;
            r->nsteps = 1;
            r->data = tmpdata;
            r->datasize = count;
            r->priv = 0;
            r->next = 0;

            vc = read_var_bb (fp, r);

            free (r->sel);
            free (r);

            if (vc == 0)
            {
                char *msg = strdup(adios_get_last_errmsg());
                adios_error ((enum ADIOS_ERRCODES) adios_errno,
                      "Cannot read data of variable %s/%s for attribute %s/%s of group %s: %s\n",
                      var_root->var_path, var_root->var_name,
                      attr_root->attr_path, attr_root->attr_name, attr_root->group_name,
                      msg);
                free(tmpdata);
                free(msg);
                return adios_errno;
            }

            *type = adios_string;
            if (file_is_fortran)
            {
                /* Fortran byte array to C string */
                *data = futils_fstr_to_cstr( tmpdata, (int)count); /* FIXME: supports only 2GB strings... */
                *size = strlen( (char *)data );
                free(tmpdata);
            }
            else
            {
                /* C array to C string */
                tmpdata[count] = '\0';
                *size = count+1;
                *data = tmpdata;
            }

            free (vc->sel);
            free (vc);
        }
        else
        {
            /* other types are inherited */
            *type = var_root->type;
            *size = bp_get_type_size (var_root->type, var_root->characteristics[var_c_index].value);
            *data = (void *) malloc (*size);
            assert (*data);
            memcpy(*data, var_root->characteristics[var_c_index].value, *size);
        }
    }

    return 0;
}

void adios_read_bp_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{
    BP_PROC * p = (BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *)(p->fh);
    struct bp_index_pg_struct_v1 ** root = &(fh->pgs_root);
    struct bp_minifooter * mh = &(fh->mfooter);
    uint64_t i;

    for (i = 0; i < mh->pgs_count; i++) {
        is_fortran ? ((*root)->adios_host_language_fortran = adios_flag_yes)
               : ((*root)->adios_host_language_fortran = adios_flag_no);
        root = &(*root)->next;
    }
}

void adios_read_bp_get_groupinfo (const ADIOS_FILE *fp, int *ngroups, char ***group_namelist, int **nvars_per_group, int **nattrs_per_group)
{
    BP_PROC * p;
    BP_FILE * fh;
    int i, j, k, offset;

    p = (BP_PROC *) fp->fh;
    fh = (BP_FILE *) p->fh;

    * ngroups = fh->gvar_h->group_count;

    alloc_namelist (group_namelist, fh->gvar_h->group_count);
    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        (*group_namelist)[i] = malloc (strlen (fh->gvar_h->namelist[i]) + 1);
        assert ((*group_namelist)[i]);

        memcpy ((*group_namelist)[i], fh->gvar_h->namelist[i], strlen (fh->gvar_h->namelist[i]) + 1);
    }

    * nvars_per_group = (int *) malloc (fh->gvar_h->group_count * sizeof (int));
    assert (* nvars_per_group);

    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        (* nvars_per_group)[i] = fh->gvar_h->var_counts_per_group[i];
    }

    * nattrs_per_group = (int *) malloc (fh->gattr_h->group_count * sizeof (int));
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

/* Check if a variable is timed. This is solely done by checking whether
 * a variable is tagged with time in XML.
 */
int adios_read_bp_is_var_timed (const ADIOS_FILE *fp, int varid)
{
    BP_PROC * p;
    BP_FILE * fh;
    struct adios_index_var_struct_v1 * v;
    struct adios_index_characteristic_struct_v1 ch;
    int retval = 0, ndim, k, dummy, file_is_fortran;
    uint64_t gdims[32];

    p = (BP_PROC *) fp->fh;
    fh = (BP_FILE *) p->fh;

    v = bp_find_var_byid (fh, varid);
    // NCSU ALACRITY-ADIOS - Use the original dimensions of the variable
    //ch = v->characteristics[0];
    struct adios_index_characteristic_dims_struct_v1 *dims = &v->characteristics[0].dims;//adios_transform_get_characteristic_original_dims(&v->characteristics[0]);
    ndim = dims->count; //ndim possibly has 'time' dimension

    log_debug ("adios_read_bp_is_var_timed: varid = %d, ndim = %d\n", varid, ndim);

    if (ndim == 0)
    {
        return 0;
    }

    for (k = 0; k < ndim; k++)
    {
        gdims[k] = dims->dims[k * 3 + 1];
    }
/*
    if (is_fortran_file (fh))
    {
        swap_order (ndim, gdims, &dummy);
    }
*/
    if (gdims[ndim - 1] == 0) // with time
    {
        retval = 1;
    }

    log_debug ("%s is_var_timed: = %d\n", v->var_name, retval);

    return retval;
}

#if 0
int64_t adios_read_bp_read_local_var (ADIOS_GROUP * gp, const char * varname,
                                      int vidx, const uint64_t * start,
                                      const uint64_t * count, void * data)
{
    struct BP_GROUP      * gh;
    BP_FILE       * fh;
    struct adios_index_var_struct_v1 * var_root;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    int    i,j,k, t, varid, start_idx, idx;
    int    ndim, ndim_notime, has_subfile, file_is_fortran;
    uint64_t size, * dims;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t datasize, nloop, dset_stride,var_stride, total_size=0, items_read;
    uint64_t count_notime[32], start_notime[32];
    int timedim = -1, temp_timedim, is_global = 0, size_of_type;
    uint64_t slice_offset, slice_size, tmpcount = 0;
    uint64_t datatimeoffset = 0; // offset in data to write a given timestep
    MPI_Status status;

    adios_errno = 0;
    if (!gp)
    {
        adios_error (err_invalid_group_struct, "Null pointer passed as group to adios_read_var()");
        return -adios_errno;
    }

    gh = (struct BP_GROUP *) gp->gh;
    if (!gh)
    {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh group handle is NULL!");
        return -adios_errno;
    }

    fh = gh->fh;
    if (!fh)
    {
        adios_error (err_invalid_group_struct, "Invalid ADIOS_GROUP struct: .gh->fh file handle is NULL!");
        return -adios_errno;
    }

    varid = adios_read_bp_find_var(gp, varname);
    if (varid < 0 || varid >= gh->vars_count)
    {
        adios_error (err_invalid_varid, "Invalid variable id %d (allowed 0..%d)", varid, gh->vars_count);
        return -adios_errno;
    }

    /* Check if file is written out by Fortran or C */
    file_is_fortran = (fh->pgs_root->adios_host_language_fortran == adios_flag_yes);

    /* Check whether we need to handle subfiles */
    has_subfile = fh->mfooter.version & ADIOS_VERSION_HAVE_SUBFILE;

    var_root = gh->vars_root; /* first variable of this group. Need to traverse the list */
    for (i = 0; i< varid && var_root; i++)
    {
        var_root = var_root->next;
    }

    if (i != varid)
    {
        adios_error (err_corrupted_variable,
                     "Variable id=%d is valid but was not found in internal data structures!",
                     varid);
        return -adios_errno;
    }

    if (vidx < 0 || vidx >= var_root->characteristics_count)
    {
        adios_error (err_out_of_bound, "idx=%d is out of bound", vidx);
    }

    ndim = var_root->characteristics [vidx].dims.count;

    /* count_notime/start_notime are working copies of count/start */
    for (i = 0; i < ndim; i++)
    {
        count_notime[i] = count[i];
        start_notime[i] = start[i];
    }

    ndim_notime = ndim;

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if (futils_is_called_from_fortran())
    {
        timedim = -1;
        swap_order(ndim_notime, count_notime, &timedim);
        swap_order(ndim_notime, start_notime, &timedim);
    }

    /* items_read = how many data elements are we going to read in total */
    items_read = 1;
    for (i = 0; i < ndim_notime; i++)
        items_read *= count_notime[i];

    size_of_type = bp_get_type_size (var_root->type, var_root->characteristics [vidx].value);

    /* READ A SCALAR VARIABLE */
    if (ndim_notime == 0)
    {
        slice_size = size_of_type;
        start_idx = 0; // OPS macros below need it
        idx = vidx; // OPS macros below need it

        if (var_root->type == adios_string)
        {
            // strings are stored without \0 in file
            // size_of_type here includes \0 so decrease by one
            size_of_type--;
        }

        /* Old BP files don't have payload_offset characteristic */
        if (var_root->characteristics[vidx].payload_offset > 0)
        {
            slice_offset = var_root->characteristics[vidx].payload_offset;

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }
        }
        else
        {
            slice_offset = 0;
            MPI_FILE_READ_OPS1
        }

        memcpy((char *)data + total_size, fh->b->buff + fh->b->offset, size_of_type);

        if (fh->mfooter.change_endianness == adios_flag_yes)
            change_endianness((char *)data + total_size
                             ,size_of_type
                             ,var_root->type
                             );

        if (var_root->type == adios_string)
        {
            // add \0 to the end of string
            // size_of_type here is the length of string
            // FIXME: how would this work for strings written over time?
            ((char*)data + total_size)[size_of_type] = '\0';
        }

        total_size += size_of_type;

        return total_size;
    } /* READ A SCALAR VARIABLE END HERE */

    /* READ AN ARRAY VARIABLE */
    uint64_t write_offset = 0;
    int npg = 0;
    tmpcount = 0;
    int flag;
    datasize = 1;
    nloop = 1;
    var_stride = 1;
    dset_stride = 1;
    uint64_t payload_size = size_of_type;

    /* To get ldims for the index vidx */
    adios_read_bp_get_dimensioncharacteristics( &(var_root->characteristics[vidx]),
                                                ldims, gdims, offsets);

    /* Again, a Fortran written file has the dimensions in Fortran order we need to swap here */
    /* Only local dims are needed for reading local vars */
    if (file_is_fortran)
    {
        i=-1;
        swap_order(ndim, ldims, &(i));
    }

    /*
    printf("ldims   = "); for (j = 0; j < ndim; j++) printf("%d ",ldims[j]); printf("\n");
    printf("count_notime   = "); for (j = 0; j < ndim_notime; j++) printf("%d ",count_notime[j]); printf("\n");
    printf("start_notime   = "); for (j = 0; j < ndim_notime; j++) printf("%d ",start_notime[j]); printf("\n");
    */

    for (j = 0; j < ndim_notime; j++)
    {
        payload_size *= ldims [j];

        if ( (start_notime[j] > ldims[j])
            || (start_notime[j] + count_notime[j] > ldims[j]))
        {
                    adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound ("
                        "the data in dimension %d to read is %llu elements from index %llu"
                        " but the actual data is [0,%llu])",
                        varid, j+1, count_notime[j], start_notime[j], ldims[j] - 1);
                    return -adios_errno;
        }
    }

    /* determined how many (fastest changing) dimensions can we read in in one read */
    int break_dim =  ndim_notime - 1;
    while (break_dim > -1)
    {
        if (start_notime[break_dim] == 0 && ldims[break_dim] == count_notime[break_dim])
        {
            datasize *= ldims[break_dim];
        }
        else
            break;

        break_dim--;
    }

    slice_offset = 0;
    slice_size = 0;
    /* Note: MPI_FILE_READ_OPS  - for reading single BP file.
     *       MPI_FILE_READ_OPS2 - for reading those with subfiles.
     *       MPI_FILE_READ_OPS1 - for reading old version of BP files
     *                            which don't contain "payload_offset"
     * Whenever to use OPS macro, start_idx and idx variable needs to be
     * properly set.
     */

    start_idx = 0;
    idx = vidx;

    if (break_dim <= 0)
    {
        /* The slowest changing dimensions should not be read completely but
           we still need to read only one block */

        uint64_t size_in_dset = count_notime[0];
        uint64_t offset_in_dset = start_notime[0];

        slice_size = (break_dim == -1 ? datasize * size_of_type : size_in_dset * datasize * size_of_type);

        if (var_root->characteristics[start_idx + idx].payload_offset > 0)
        {
            slice_offset = var_root->characteristics[start_idx + idx].payload_offset
                         + offset_in_dset * datasize * size_of_type;

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }
        }
        else
        {
            slice_offset = 0;
            MPI_FILE_READ_OPS1
        }

        memcpy ((char *)data, fh->b->buff + fh->b->offset, slice_size);
        if (fh->mfooter.change_endianness == adios_flag_yes)
        {
            change_endianness((char *)data + write_offset, slice_size, var_root->type);
        }
    }
    else
    {
        uint64_t stride_offset = 0;
        uint64_t * size_in_dset, * offset_in_dset, * offset_in_var;
        uint64_t start_in_payload, end_in_payload, s;
        uint64_t var_offset;
        uint64_t dset_offset;

        size_in_dset = (uint64_t *) malloc (8 * ndim_notime);
        offset_in_dset = (uint64_t *) malloc (8 * ndim_notime);
        offset_in_var = (uint64_t *) malloc (8 * ndim_notime);

        if (size_in_dset == 0 || offset_in_dset == 0 || offset_in_var == 0)
        {
             adios_error (err_no_memory, "Malloc failed in %s at %d\n"
                         , __FILE__, __LINE__
                         );
             return -adios_errno;
        }

        for (i = 0; i < ndim_notime ; i++)
        {
            size_in_dset[i] = count_notime[i];
            offset_in_dset[i] = start_notime[i];
            offset_in_var[i] = 0;
        }

        datasize = 1;
        var_stride = 1;
        for (i = ndim_notime - 1; i >= break_dim; i--)
        {
            datasize *= size_in_dset[i];
            dset_stride *= ldims[i];
            var_stride *= count_notime[i];
        }

        /* Calculate the size of the chunk we are trying to read in */
        start_in_payload = 0;
        end_in_payload = 0;
        s = 1;
        for (i = ndim_notime - 1; i >= 0; i--)
        {
            start_in_payload += s * offset_in_dset[i] * size_of_type;
            end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
            s *= ldims[i];
        }
        slice_size = end_in_payload - start_in_payload + 1 * size_of_type;

        if (var_root->characteristics[start_idx + idx].payload_offset > 0)
        {
            slice_offset =  var_root->characteristics[start_idx + idx].payload_offset
                          + start_in_payload;
            if (!has_subfile)
            {
                MPI_FILE_READ_OPS
            }
            else
            {
                MPI_FILE_READ_OPS2
            }

            for ( i = 0; i < ndim_notime ; i++)
            {
                offset_in_dset[i] = 0;
            }
        }
        else
        {
            slice_offset =  start_in_payload;
            MPI_FILE_READ_OPS1
        }

        var_offset = 0;
        dset_offset = 0;
        for (i = 0; i < ndim_notime ; i++)
        {
            var_offset = offset_in_var[i] + var_offset * count_notime[i];
            dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
        }

        copy_data (data
                  ,fh->b->buff + fh->b->offset
                  ,0
                  ,break_dim
                  ,size_in_dset
                  ,ldims
                  ,count_notime
                  ,var_stride
                  ,dset_stride
                  ,var_offset
                  ,dset_offset
                  ,datasize
                  ,size_of_type
                  );

        free (size_in_dset);
        free (offset_in_dset);
        free (offset_in_var);
    }

    total_size += items_read * size_of_type;

    return total_size;
}

// NCSU - Timer series analysis, correlation
double adios_stat_cor (ADIOS_VARINFO * vix, ADIOS_VARINFO * viy, char * characteristic, uint32_t time_start, uint32_t time_end, uint32_t lag)
{
    int i,j;

    double avg_x = 0.0, avg_y = 0.0, avg_lag = 0.0;
    double var_x = 0.0, var_y = 0.0, var_lag = 0.0;
    double cov = 0;

    if (vix == NULL)
    {
        fprintf(stderr, "Variable not defined\n");
        return 0;
    }

    // If the vix and viy are not time series objects, return.
    if ((vix->timedim < 0) && (viy->timedim < 0))
    {
        fprintf(stderr, "Covariance must involve timeseries data\n");
        return 0;
    }

    uint32_t min = vix->dims[0] - 1;
    if (viy && (min > viy->dims[0] - 1))
        min = viy->dims[0] - 1;

    if(time_start == 0 && time_end == 0)
    { //global covariance
        if(viy == NULL) {
            fprintf(stderr, "Must have two variables for global covariance\n");
            return 0;
        }

        // Assign vix to viy, and calculate covariance
        viy = vix;
        time_start = 0;
        time_end = min;
    }
    // Check the bounds of time
    if (    (time_start >= 0) && (time_start <= min)
            &&      (time_end >= 0)   && (time_end <= min)
            &&  (time_start <= time_end))
    {
        if(viy == NULL) //user must want to run covariance against itself
        {
            if(! (time_end+lag) > min)
            {
                fprintf(stderr, "Must leave enough timesteps for lag\n");
                return 0;
            }

            if (strcmp(characteristic, "average") == 0 || strcmp(characteristic, "avg") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (adios_double, vix->avgs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (adios_double, vix->avgs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (adios_double, vix->avgs[i]);
                    double val_lag = bp_value_to_double (adios_double, vix->avgs[i + lag]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_lag += (val_lag - avg_lag) * (val_lag - avg_lag) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_lag - avg_lag) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "standard deviation") == 0 || strcmp(characteristic, "std_dev") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (adios_double, vix->std_devs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (adios_double, vix->std_devs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (adios_double, vix->std_devs[i]);
                    double val_lag = bp_value_to_double (adios_double, vix->std_devs[i + lag]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_lag += (val_lag - avg_lag) * (val_lag - avg_lag) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_lag - avg_lag) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "minimum") == 0 || strcmp(characteristic, "min") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (vix->type, vix->mins[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (vix->type, vix->mins[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (vix->type, vix->mins[i]);
                    double val_lag = bp_value_to_double (vix->type, vix->mins[i + lag]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_lag += (val_lag - avg_lag) * (val_lag - avg_lag) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_lag - avg_lag) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "maximum") == 0 || strcmp(characteristic, "max") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (vix->type, vix->maxs[i]);
                    double val_lag = bp_value_to_double (vix->type, vix->maxs[i + lag]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_lag += (val_lag - avg_lag) * (val_lag - avg_lag) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_lag - avg_lag) / (time_end - time_start + 1);
                }
            }
            else
            {
                fprintf (stderr, "Unknown characteristic\n");
                return 0;
            }
            return cov / (sqrt (var_x) * sqrt (var_lag));
        }
        else
        {
            if (strcmp(characteristic, "average") == 0 || strcmp(characteristic, "avg") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(adios_double, vix->avgs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(adios_double, viy->avgs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (adios_double, vix->avgs[i]);
                    double val_y = bp_value_to_double (adios_double, viy->avgs[i]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_y += (val_y - avg_y) * (val_y - avg_y) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_y - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "standard deviation") == 0 || strcmp(characteristic, "std_dev") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(adios_double, vix->std_devs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(adios_double, viy->std_devs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (adios_double, vix->std_devs[i]);
                    double val_y = bp_value_to_double (adios_double, viy->std_devs[i]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_y += (val_y - avg_y) * (val_y - avg_y) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_y - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "minimum") == 0 || strcmp(characteristic, "min") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(vix->type, vix->mins[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(viy->type, viy->mins[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (vix->type, vix->mins[i]);
                    double val_y = bp_value_to_double (viy->type, viy->mins[i]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_y += (val_y - avg_y) * (val_y - avg_y) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_y - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "maximum") == 0 || strcmp(characteristic, "max") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(vix->type, viy->maxs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    double val_x = bp_value_to_double (vix->type, vix->maxs[i]);
                    double val_y = bp_value_to_double (viy->type, viy->maxs[i]);
                    var_x += (val_x - avg_x) * (val_x - avg_x) / (time_end - time_start + 1);
                    var_y += (val_y - avg_y) * (val_y - avg_y) / (time_end - time_start + 1);
                    cov += (val_x - avg_x) * (val_y - avg_y) / (time_end - time_start + 1);
                }
            }
            else
            {
                fprintf (stderr, "Unknown characteristic\n");
                return 0;
            }
            return cov / (sqrt (var_x) * sqrt (var_y));
        }
    }
    else
    {
        fprintf (stderr, "Time values out of bounds\n");
        return 0;
    }
}

// NCSU - Time series analysis, covariance
//covariance(x,y) = sum(i=1,..N) [(x_1 - x_mean)(y_i - y_mean)]/N
double adios_stat_cov (ADIOS_VARINFO * vix, ADIOS_VARINFO * viy, char * characteristic, uint32_t time_start, uint32_t time_end, uint32_t lag)
{
    int i,j;

    double avg_x = 0.0, avg_y = 0.0, avg_lag = 0.0;
    double cov = 0;

    if (vix == NULL)
    {
        fprintf(stderr, "Variable not defined\n");
        return 0;
    }

    // If the vix and viy are not time series objects, return.
    if ((vix->timedim < 0) && (viy->timedim < 0))
    {
        fprintf(stderr, "Covariance must involve timeseries data\n");
        return 0;
    }

    uint32_t min = vix->dims[0] - 1;
    if (viy && (min > viy->dims[0] - 1))
        min = viy->dims[0] - 1;

    if(time_start == 0 && time_end == 0)
    { //global covariance
        if(viy == NULL) {
            fprintf(stderr, "Must have two variables for global covariance\n");
            return 0;
        }

        // Assign vix to viy, and calculate covariance
        viy = vix;
        time_start = 0;
        time_end = min;
    }
    // Check the bounds of time
    if (    (time_start >= 0) && (time_start <= min)
            &&      (time_end >= 0)   && (time_end <= min)
            &&  (time_start <= time_end))
    {
        if(viy == NULL) //user must want to run covariance against itself
        {
            if(! (time_end+lag) > min)
            {
                fprintf(stderr, "Must leave enough timesteps for lag\n");
                return 0;
            }

            if (strcmp(characteristic, "average") == 0 || strcmp(characteristic, "avg") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (adios_double, vix->avgs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (adios_double, vix->avgs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                    cov += (bp_value_to_double (adios_double, vix->avgs[i]) - avg_x) * (bp_value_to_double (adios_double, vix->avgs[i+lag]) - avg_lag) / (time_end - time_start + 1);
            }
            else if (strcmp(characteristic, "standard deviation") == 0 || strcmp(characteristic, "std_dev") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (adios_double, vix->std_devs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (adios_double, vix->std_devs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                    cov += (bp_value_to_double (adios_double, vix->std_devs[i]) - avg_x) * (bp_value_to_double (adios_double, vix->std_devs[i+lag]) - avg_lag) / (time_end - time_start + 1);
            }
            else if (strcmp(characteristic, "minimum") == 0 || strcmp(characteristic, "min") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (vix->type, vix->mins[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (vix->type, vix->mins[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                    cov += (bp_value_to_double (vix->type, vix->mins[i]) - avg_x) * (bp_value_to_double (vix->type, vix->mins[i+lag]) - avg_lag) / (time_end - time_start + 1);
            }
            else if (strcmp(characteristic, "maximum") == 0 || strcmp(characteristic, "max") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double (vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                    avg_lag += bp_value_to_double (vix->type, vix->maxs[i + lag]) / (time_end - time_start + 1);
                }

                for (i = time_start; i <= time_end; i ++)
                    cov += (bp_value_to_double (vix->type, vix->maxs[i]) - avg_x) * (bp_value_to_double (vix->type, vix->maxs[i+lag]) - avg_lag) / (time_end - time_start + 1);
            }
            else
            {
                fprintf (stderr, "Unknown characteristic\n");
                return 0;
            }
        }
        else
        {
            if (strcmp(characteristic, "average") == 0 || strcmp(characteristic, "avg") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(adios_double, vix->avgs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(adios_double, viy->avgs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    cov += (bp_value_to_double(adios_double, vix->avgs[i]) - avg_x) * (bp_value_to_double(adios_double, viy->avgs[i]) - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "standard deviation") == 0 || strcmp(characteristic, "std_dev") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(adios_double, vix->std_devs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(adios_double, viy->std_devs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    cov += (bp_value_to_double(adios_double, vix->std_devs[i]) - avg_x) * (bp_value_to_double(adios_double, viy->std_devs[i]) - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "minimum") == 0 || strcmp(characteristic, "min") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(vix->type, vix->mins[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(viy->type, viy->mins[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    cov += (bp_value_to_double(vix->type, vix->mins[i]) - avg_x) * (bp_value_to_double(viy->type, viy->mins[i]) - avg_y) / (time_end - time_start + 1);
                }
            }
            else if (strcmp(characteristic, "maximum") == 0 || strcmp(characteristic, "max") == 0)
            {
                for (i = time_start; i <= time_end; i ++)
                {
                    avg_x += bp_value_to_double(vix->type, vix->maxs[i]) / (time_end - time_start + 1);
                    avg_y += bp_value_to_double(vix->type, viy->maxs[i]) / (time_end - time_start + 1);
                }
                for (i = time_start; i <= time_end; i ++)
                {
                    cov += (bp_value_to_double(vix->type, vix->maxs[i]) - avg_x) * (bp_value_to_double(viy->type, viy->maxs[i]) - avg_y) / (time_end - time_start + 1);
                }
            }
            else
            {
                fprintf (stderr, "Unknown characteristic\n");
                return 0;
            }
        }
    }
    else
    {
        fprintf (stderr, "Time values out of bounds\n");
    }
    return cov;
}
#endif
