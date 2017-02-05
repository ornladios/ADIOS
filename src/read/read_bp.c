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
#include "core/adios_internals.h"
#include "core/util.h"
#include "core/bp_utils.h"
#include "core/bp_types.h"
#include "core/adios_read_hooks.h"
#include "core/futils.h"
#include "core/common_read.h"
#include "core/adios_logger.h"
#include "core/a2sel.h"
#include "core/adios_clock.h"
#include "core/adios_selection_util.h"

#include "core/transforms/adios_transforms_transinfo.h"
#include "core/transforms/adios_transforms_common.h" // NCSU ALACRITY-ADIOS

#ifdef DMALLOC
#include "dmalloc.h"
#endif

static int chunk_buffer_size = 1024*1024*16;
static int poll_interval_msec = 10000; // 10 secs by default
static int show_hidden_attrs = 0; // don't show hidden attr by default

static ADIOS_VARCHUNK * read_var_bb  (const ADIOS_FILE * fp, read_request * r);
static ADIOS_VARCHUNK * read_var_pts (const ADIOS_FILE * fp, read_request * r);
static ADIOS_VARCHUNK * read_var_wb  (const ADIOS_FILE * fp, read_request * r);

static int map_req_varid (const ADIOS_FILE * fp, int varid);
static int adios_wbidx_to_pgidx (const ADIOS_FILE * fp, read_request * r, int step_offset);

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

/* emulate MPI_File_read with support for data >2GiB */
#define MPI_FILE_READ64(fh, buff, size, type, status)                               \
    {                                                                               \
        int read_len = 0;                                                           \
        uint64_t total_read = 0;                                                    \
        uint64_t to_read = size;                                                    \
        int count;                                                                  \
        while (to_read > 0)                                                         \
        {                                                                           \
            read_len = (to_read > MAX_MPIWRITE_SIZE) ? MAX_MPIWRITE_SIZE : to_read; \
            MPI_File_read (fh                                                       \
                          ,(char*)buff + total_read                                 \
                          ,read_len                                                 \
                          ,type                                                     \
                          ,status                                                   \
                          );                                                        \
            MPI_Get_count (status, type, &count);                                   \
            int err;                                                                \
            if (count != read_len)                                                  \
            {                                                                       \
                log_error ("Need to do multi-read (tried: "                         \
                        "%d read: %d) errno %d\n",                                  \
                        read_len, count, errno);                                    \
                err = count;                                                        \
                break;                                                              \
            }                                                                       \
            total_read += count;                                                    \
            to_read -= count;                                                       \
        }                                                                           \
    }

#define MPI_FILE_READ_OPS1                          \
        bp_realloc_aligned(fh->b, slice_size);      \
        fh->b->offset = 0;                          \
                                                    \
        MPI_File_seek (fh->mpi_fh                   \
                      ,(MPI_Offset)slice_offset     \
                      ,MPI_SEEK_SET                 \
                      );                            \
        MPI_FILE_READ64 (fh->mpi_fh                 \
                        ,fh->b->buff                \
                        ,slice_size                 \
                        ,MPI_BYTE                   \
                        ,&status                    \
                        );                          \
        fh->b->offset = 0;                          \

// To read subfiles
#define MPI_FILE_READ_OPS2                                                                  \
        bp_realloc_aligned(fh->b, slice_size);                                              \
        fh->b->offset = 0;                                                                  \
                                                                                            \
        MPI_File * sfh;                                                                     \
        sfh = get_BP_subfile_handle (fh, v->characteristics[start_idx + idx].file_index);   \
        if (!sfh)                                                                           \
        {                                                                                   \
            int err;                                                                        \
            char * ch, * name_no_path, * name;                                              \
            MPI_Info info = MPI_INFO_NULL;                                                  \
            struct BP_file_handle * new_h =                                                 \
                  (struct BP_file_handle *) malloc (sizeof (struct BP_file_handle));        \
            new_h->file_index = v->characteristics[start_idx + idx].file_index;             \
            new_h->next = 0;                                                                \
            if ( (ch = strrchr (fh->fname, '/')) )                                          \
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
           add_BP_subfile_handle (fh, new_h);                                                  \
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
        MPI_FILE_READ64 (*sfh                                                               \
                        ,fh->b->buff                                                        \
                        ,slice_size                                                         \
                        ,MPI_BYTE                                                           \
                        ,&status                                                            \
                        );                                                                  \
        fh->b->offset = 0;                                                                  \

//We also need to be able to read old .bp which doesn't have 'payload_offset'
#define MPI_FILE_READ_OPS3                                                                  \
        MPI_File_seek (fh->mpi_fh                                                           \
                      ,(MPI_Offset) v->characteristics[start_idx + idx].offset       \
                      ,MPI_SEEK_SET);                                                       \
        MPI_File_read (fh->mpi_fh, fh->b->buff, 8, MPI_BYTE, &status);                      \
        tmpcount= *((uint64_t*)fh->b->buff);                                                \
                                                                                            \
        bp_realloc_aligned(fh->b, tmpcount + 8);                                            \
        fh->b->offset = 0;                                                                  \
                                                                                            \
        MPI_File_seek (fh->mpi_fh                                                           \
                      ,(MPI_Offset) (v->characteristics[start_idx + idx].offset)     \
                      ,MPI_SEEK_SET);                                                       \
        MPI_FILE_READ64 (fh->mpi_fh, fh->b->buff, tmpcount + 8, MPI_BYTE, &status);         \
        fh->b->offset = 0;                                                                  \
        adios_parse_var_data_header_v1 (fh->b, &var_header);                                \

// NCSU ALACRITY-ADIOS: After much pain and consideration, I've decided to implement a
//     2nd version of this function to avoid substantial wasted time in the writeblock method
#define MPI_FILE_READ_OPS1_BUF(buf)                 \
        MPI_File_seek (fh->mpi_fh                   \
                      ,(MPI_Offset)slice_offset     \
                      ,MPI_SEEK_SET                 \
                      );                            \
        MPI_FILE_READ64 (fh->mpi_fh                 \
                        ,(buf)                      \
                        ,slice_size                 \
                        ,MPI_BYTE                   \
                        ,&status                    \
                        );

// To read subfiles
#define MPI_FILE_READ_OPS2_BUF(buf)                                                         \
        MPI_File * sfh;                                                                     \
        sfh = get_BP_subfile_handle (fh, v->characteristics[start_idx + idx].file_index);   \
        if (!sfh)                                                                           \
        {                                                                                   \
            int err;                                                                        \
            char * ch, * name_no_path, * name;                                              \
            MPI_Info info = MPI_INFO_NULL;                                                  \
            struct BP_file_handle * new_h =                                                 \
                  (struct BP_file_handle *) malloc (sizeof (struct BP_file_handle));        \
            new_h->file_index = v->characteristics[start_idx + idx].file_index;             \
            new_h->next = 0;                                                                \
            if ((ch = strrchr (fh->fname, '/')))                                              \
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
           add_BP_subfile_handle (fh, new_h);                                                 \
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
        MPI_FILE_READ64 (*sfh                                                               \
                        ,(buf)                                                              \
                        ,slice_size                                                         \
                        ,MPI_BYTE                                                           \
                        ,&status                                                            \
                        );


/* This routine release one step. It only frees the var/attr namelist. */
static void release_step (ADIOS_FILE *fp)
{
    BP_PROC * p = GET_BP_PROC ((const ADIOS_FILE *)fp);

    if (p->varid_mapping)
    {
        free (p->varid_mapping);
        p->varid_mapping = 0;
    }

    if (fp->var_namelist)
    {
        a2s_free_namelist (fp->var_namelist, fp->nvars);
        fp->var_namelist = 0;
        fp->nvars = 0;
    }

    if (fp->attr_namelist)
    {
        a2s_free_namelist (fp->attr_namelist, fp->nattrs);
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
    int rank, file_ok;
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

    fh = BP_FILE_alloc (fname, comm);

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

    p = (BP_PROC *) malloc (sizeof (BP_PROC));
    assert (p);
    p->fh = fh;
    p->streaming = 1;
    p->varid_mapping = 0;
    p->local_read_request_list = 0;
    p->b = 0;
    p->priv = 0;

    fp->fh = (uint64_t) p;
    fp->file_size = fh->mfooter.file_size;
    fp->version = fh->mfooter.version & ADIOS_VERSION_NUM_MASK;
    fp->endianness = bp_get_endianness (fh->mfooter.change_endianness);
    fp->last_step = fh->tidx_stop - 1;

    /* Seek to the initial step. */
    release_step (fp);
    bp_seek_to_step (fp, 0, show_hidden_attrs);

    /* For file, the last step is tidx_stop */
    fp->last_step = fh->tidx_stop - 1;

    return;
}

static int get_new_step (ADIOS_FILE * fp, const char * fname, MPI_Comm comm, int last_tidx, float timeout_sec)
{
    BP_FILE * new_fh;
    double t1 = adios_gettime_double();

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
        new_fh = open_file (fname, comm);
        if (!new_fh)
        {
            // file is bad so keeps polling.
            stay_in_poll_loop = 1;
        }
        else if (new_fh && new_fh->tidx_stop == last_tidx)
        {
            // file is good but no new steps in it. Continue polling.
            bp_close (new_fh);
            stay_in_poll_loop = 1;
        }
        else
        {
            // the file looks good and there are new steps written.
            build_ADIOS_FILE_struct (fp, new_fh);
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
            else if (timeout_sec > 0.0 && (adios_gettime_double () - t1 > timeout_sec))
            {
                log_debug ("Time is out in get_new_step()\n");
                stay_in_poll_loop = 0;
            }
            else
            {
                adios_nanosleep (poll_interval_msec/1000,
                    (int)(((uint64_t)poll_interval_msec * 1000000L)%1000000000L));
            }
        }

    } // while (stay_in_poll_loop)

    log_debug ("exit get_new_step\n");

    return found_stream;
}


//
// returns # of elements in a bounding box' range (not size in bytes!)
//
static uint64_t adios_get_nelements_of_box (int ndim, uint64_t* start, uint64_t* count)
{
    int k=0;
    uint64_t bbsize = 1;
    for (k=0; k<ndim; k++) {
        bbsize *= count[k];
    }
    return bbsize;
}


//
// calculate the spanning N-dim bounding box for a list of N-dim points
//
static void mGetPointlistSpan(ADIOS_SELECTION_POINTS_STRUCT* pts, uint64_t* start, uint64_t* count)
{
    uint64_t i=0, idx=0;
    int d=0;
    uint64_t max[pts->ndim];
    for (i = 0; i < pts->npoints; i++)
    {
        idx = i * pts->ndim;
        //printf("%ldth = [%ld, %ld, %ld] \n", i, pts->points[idx], pts->points[idx+1], pts->points[idx+2]);
        for (d = 0; d < pts->ndim; d++) {
            if (i == 0) {
                start[d] = pts->points[d]; max[d] = pts->points[d];
                continue;
            }

            uint64_t curr = pts->points[idx+d];
            if ((start[d] > curr)) {
                start[d] = curr;
            }
            if (max[d] < curr) {
                max[d] = curr;
            }
        }
        //printf("start[0]=%ld  max[0]=%ld \n", start[0], max[0]);
    }
    for (d = 0; d < pts->ndim; d++) {
        count[d] = max[d] - start[d] + 1;
    }
    return;
}

/*
 *  Calculate the spanning N-dim bounding box for a list of 1-dim points in an N-dim box.
 *  Simply take the smallest and largest offset, then make a bounding box that fits all of them
 *  Note that it is not giving the smallest container box. It can only decrease the box in the
 *  slowest dimension without converting all points to N-dimension.
 */
static void mGetPointlistSpan1D(ADIOS_SELECTION_POINTS_STRUCT* pts, int ndim,
        uint64_t* boxstart, uint64_t* boxcount,
        uint64_t* spanstart, uint64_t* spancount)
{
    uint64_t i=0;
    int d=0;
    uint64_t span[2]; // min and max offsets
    span[0] = 0xFFFFFFFFFFFFFFFF; // 18446744073709551615, the min offset
    span[1] = 0; // the max offset

    // find min and max offsets
    uint64_t *pt = pts->points;
    for (i = 0; i < pts->npoints; i++)
    {
        if (*pt < span[0])
            span[0] = *pt;
        if (*pt > span[1])
            span[1] = *pt;
        pt++;
    }

    // convert them to N-dim
    uint64_t spanND[2*ndim];
    a2sel_points_1DtoND_box (2, span, ndim, boxstart, boxcount, 1, spanND);

    // correct sub-dimensions (some other points may be outside the naive span over two points
    spanstart[0] = spanND[0];
    spancount[0] = spanND[ndim]-spanND[0]+1;
    for (d = 1; d < ndim; d++) {
        spanstart[d] = boxstart[d];
        spancount[d] = boxcount[d];
    }
    return;
}

#if 0
//
// returns # of segments to divide for a bounding box' range
//
static int mGetRange(int ndim, uint64_t* start, uint64_t* count)
{
    int k=0;
    uint64_t bbsize = 1;
    for (k=0; k<ndim; k++) {
        bbsize *= count[k] - start[k]+1;
        log_debug ("... bb at %d dimention: [%" PRIu64 ", %" PRIu64 "]\n", k, count[k], start[k]);
    }

    const uint64_t BBSIZELIMIT = 1048576; /* 1M elements, 4-8MB data usually to read at once */
    uint64_t nBB = bbsize/BBSIZELIMIT ;

    if (nBB * BBSIZELIMIT != bbsize) {
        nBB += 1;
    }

    log_debug ("... nBB=%" PRIu64 "\n", nBB);
    //return bbsize;
    return (int) nBB;
}
#endif

/*
 * Pick the points of a point list from a contiguous array of data.
 * The point coordinates are relative to the original 'bndim' dimensional container
 * while data is in a sub-box of the original container:
 * contstart[] should be the starting offsets of  the original container (points are relative to this)
 * contcount[] should be the size of the original container
 * nelems is the sum of elements (= PROD(count[i], i=0,...,bndim-1))
 * substart[] should be the starting offsets of the sub-box in the original container (data is contained in this, with nelems values)
 * subcount[] should be the size of the sub-box
 * src is the data array containing data of the sub-box
 * dest is the output data array for the points
 * It returns the number of points that fall inside the sub-box, i.e. the number of copied elements
 */
static uint64_t pick_points_from_boundingbox (ADIOS_SELECTION *sel, int size_of_type,
                                              int bndim, uint64_t *contstart, uint64_t *contcount,
                                              uint64_t nelems, uint64_t *substart, uint64_t *subcount,
                                              char *src, char *dest)
{
    uint64_t npoints = 0;
    uint64_t j;
    int d;
    assert (sel->type == ADIOS_SELECTION_POINTS);
    int pndim = sel->u.points.ndim;
    assert (bndim==pndim || pndim==1);

    uint64_t subproduct[bndim+1]; // number of elements in subcount
    subproduct[bndim] = 1;
    for (d = bndim-1; d >= 0; d--) {
        subproduct[d] = subproduct[d+1] * subcount[d];
    }
    assert (nelems == subproduct[0]);

    uint64_t suboffs[bndim]; // N-D offsets of starting points of sub-box from the starting point of original container
    for (d = 0; d < bndim; d++) {
        assert (substart[d] >= contstart[d]);
        suboffs[d] = substart[d] - contstart[d];
    }

    uint64_t suboffset = 0; // 1D offset of starting point of sub-box from the starting point of original container
    for (d = bndim-1; d >= 0; d--) {
        suboffset += suboffs[d] * subproduct[d+1];
    }

    if (pndim == 1)
    {
        // 1D local points in N-D block
        uint64_t *pt = sel->u.points.points;

        for (j = 0; j < sel->u.points.npoints; j++)
        {
            if (suboffset <= *pt && *pt-suboffset < nelems)
            {
                memcpy (dest,
                        src + (*pt-suboffset)*size_of_type,
                        size_of_type);
                npoints++;
                dest += size_of_type;
            }
            pt++;
        }
    }
    else
    {
        uint64_t *pt = sel->u.points.points; // first N-dim point in list
        for (j = 0; j < sel->u.points.npoints; j++)
        {
            int64_t idxInBB = 0;
            for (d = 0; d < bndim; d++)
            {
                if (suboffs[d] <= pt[d] && pt[d]-suboffs[d] < subcount[d])
                {
                    idxInBB += (pt[d] - suboffs[d]) * subproduct[d+1];
                } else {
                    // this point is outside of the sub-box
                    idxInBB = -1;
                    break;
                }
            }

            if (idxInBB >= 0)
            {
                memcpy(dest, src+idxInBB*size_of_type, size_of_type);
                //printf(" checking: %.3f vs %.3f \n", ((double*)(nr->data))[idxInBB], ((double*)(r->data))[i]);
                //printf("checking: [%ld th bb]  [point %ld],  idxInBB=%ld value %.3f vs %.3f\n",j, i, idxInBB, ((double*)(nr->data))[idxInBB], ((double*)(r->data))[step]);
                npoints++;
                dest += size_of_type; // next slot for data in output data array
            }
            pt += bndim; // next N-dim point in list
        }
    }
    return npoints;
}

/* This routine processes a read request and returns data in ADIOS_VARCHUNK.
   If the selection type is not bounding box, convert it. The basic file reading
   functionality is implemented in read_var_bb() routine.
*/
static ADIOS_VARCHUNK * read_var (const ADIOS_FILE * fp, read_request * r)
{
    log_debug ("read_var()\n");

    ADIOS_VARCHUNK * chunk = NULL;
    ADIOS_SELECTION * sel = r->sel;

    switch (sel->type)
    {
        case ADIOS_SELECTION_BOUNDINGBOX:
            chunk = read_var_bb (fp, r);
            break;
        case ADIOS_SELECTION_POINTS:
            chunk = read_var_pts (fp, r);
            break;
        case ADIOS_SELECTION_WRITEBLOCK:
            chunk = read_var_wb (fp, r);
            break;
        case ADIOS_SELECTION_AUTO:
            break;
        default:
            log_debug ("ADIOS selection type is wrong\n");
            break;
    }

    return chunk;
}

/* This routine reads in data for bounding box selection.
   If the selection is not bounding box, it should be converted to it.
   The data returned is saved in ADIOS_VARCHUNK.
 */
static ADIOS_VARCHUNK * read_var_bb (const ADIOS_FILE *fp, read_request * r)
{
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);

    ADIOS_SELECTION * sel;
    struct adios_index_var_struct_v1 * v;
    int i, j, t, time, nsteps;
    int64_t start_idx, stop_idx, idx;
    int ndim, has_subfile, file_is_fortran;
    uint64_t * dims, tmpcount;
    uint64_t ldims[32], gdims[32], offsets[32];
    uint64_t datasize, dset_stride,var_stride, total_size=0, items_read;
    uint64_t * count, * start;
    void * data;
    int dummy = -1, is_global = 0, size_of_type;
    uint64_t slice_offset, slice_size;
    MPI_Status status;
    ADIOS_VARCHUNK * chunk;
    struct adios_var_header_struct_v1 var_header;
    //struct adios_var_payload_struct_v1 var_payload;

//    log_debug ("read_var_bb()\n");
    file_is_fortran = is_fortran_file (fh);
    has_subfile = has_subfiles (fh);

    sel = r->sel;
    start = sel->u.bb.start;
    count = sel->u.bb.count;
    data = r->data;

    v = bp_find_var_byid (fh, r->varid);

    /* Get dimensions and flip if caller != writer language */
    /* Note: ndim below doesn't include time if there is any */
    // NCSU ALACRITY-ADIOS - Note: this function has been modified to return
    //   the "raw" dimensions (i.e., 1D byte array)
    bp_get_and_swap_dimensions (fp, v, file_is_fortran, &ndim, &dims, &nsteps, file_is_fortran);

    assert (ndim == sel->u.bb.ndim);
    ndim = sel->u.bb.ndim;

    /* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
       We need to swap them here to read correctly in C order */
    if (futils_is_called_from_fortran ())
    {
        swap_order (ndim, start, &dummy);
        swap_order (ndim, count, &dummy);
    }
/*
    log_debug ("read_var_bb(): start(");
    for (i = 0; i < ndim; i++)
    {
        log_debug_cont ("%lu", start[i]);
        if (i != ndim - 1)
        {
            log_debug_cont (",");
        }
    }
    log_debug_cont (")\n");

    log_debug ("read_var_bb(): count(");
    for (i = 0; i < ndim; i++)
    {
        log_debug_cont ("%lu", count[i]);
        if (i != ndim - 1)
        {
            log_debug_cont (",");
        }
    }
    log_debug_cont (")\n");
*/
    /* items_read = how many data elements are we going to read in total (per timestep) */
    items_read = 1;
    for (i = 0; i < ndim; i++)
    {
        items_read *= count[i];
    }

    size_of_type = bp_get_type_size (v->type, v->characteristics [0].value);

//    log_debug ("read_var_bb: from_steps = %d, nsteps = %d\n", r->from_steps, r->nsteps);

    /* Note fp->current_step is always 0 for file mode. */
    for (t = fp->current_step + r->from_steps; t < fp->current_step + r->from_steps + r->nsteps; t++)
    {

        if (!p->streaming)
        {
            time = get_time (v, t);
        }
        else
        {
            // Fix: the assumption that for streaming mode, the time in file
            // always starts from 1 is not correct. So here we add fh->tidx_start to adjust
            // Q. Liu, 06/2013
            time = fh->tidx_start + t;
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

            if (v->characteristics[start_idx + idx].payload_offset > 0)
            {
                if (!has_subfile)
                {
                    MPI_FILE_READ_OPS1
                }
                else
                {
                    MPI_FILE_READ_OPS2
                }
            }
            else
            {
                slice_offset = 0;
                MPI_FILE_READ_OPS3
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
/*
                printf ("count   = "); for (j = 0; j<ndim; j++) printf ("%d ",count[j]); printf ("\n");
                printf ("start   = "); for (j = 0; j<ndim; j++) printf ("%d ",start[j]); printf ("\n");
*/
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
                printf ("ldims   = "); for (j = 0; j<ndim; j++) printf ("%d ",ldims[j]); printf ("\n");
                printf ("gdims   = "); for (j = 0; j<ndim; j++) printf ("%d ",gdims[j]); printf ("\n");
                printf ("offsets = "); for (j = 0; j<ndim; j++) printf ("%d ",offsets[j]); printf ("\n");
*/
                for (j = 0; j < ndim; j++)
                {
                    payload_size *= ldims [j];

                    if ( (count[j] > gdims[j])
                      || (start[j] > gdims[j])
                      || (start[j] + count[j] > gdims[j]))
                    {
                        adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound 1("
                            "the data in dimension %d to read is %" PRIu64 " elements from index %" PRIu64
                            " but the actual data is [0,%" PRId64 "])\n",
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

                if (hole_break == -1)
                {
                    /* The complete read happens to be exactly one pg, and the entire pg */
                    /* This means we enter this only once, and npg=1 at the end */
                    /* This is a rare case. FIXME: cannot eliminate this? */
                    slice_size = payload_size;

                    slice_offset = v->characteristics[start_idx + idx].payload_offset;
                    if (v->characteristics[start_idx + idx].payload_offset > 0)
                    {
                        if (!has_subfile)
                        {
                            MPI_FILE_READ_OPS1
                        }
                        else
                        {
                            MPI_FILE_READ_OPS2
                        }
                    }
                    else
                    {
                         slice_offset = 0;
                         MPI_FILE_READ_OPS3
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
                    uint64_t isize;
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
                    else
                    {
                        slice_offset = 0;
                        MPI_FILE_READ_OPS3
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
                    uint64_t isize;
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
                    if (v->characteristics[start_idx + idx].payload_offset > 0)
                    {
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
                    }
                    else
                    {
                        slice_offset =  start_in_payload;
                        MPI_FILE_READ_OPS3
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

                    adios_util_copy_data (data
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
                              ,fh->mfooter.change_endianness
                              ,v->type
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

//    log_debug ("read_var_bb(): build ADIOS_VARCHUNK\n");
    chunk = (ADIOS_VARCHUNK *) malloc (sizeof (ADIOS_VARCHUNK));
    assert (chunk);

    chunk->varid = r->varid;
    chunk->type = v->type;
    // NCSU ALACRITY-ADIOS - Added timestep information into varchunks
    chunk->from_steps = r->from_steps;
    chunk->nsteps = r->nsteps;
    chunk->sel = a2sel_copy (r->sel);
    chunk->data = r->data;
    return chunk;
}

/* This routine reads in data for point selection.
*/
static ADIOS_VARCHUNK * read_var_pts (const ADIOS_FILE *fp, read_request * r)
{
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);
    int size_of_type;

    uint64_t i;
    uint64_t nerr = 0; // number of out of bound errors (points outside of subselection block)
    ADIOS_SELECTION * sel = r->sel;
    assert (sel->type == ADIOS_SELECTION_POINTS);
    ADIOS_VARCHUNK * chunk = NULL;
    int pndim = sel->u.points.ndim; // dimensionality of points
    int bndim; // dimensionality of subselection block
    struct adios_index_var_struct_v1 * v;
    int step, d;
    uint64_t zeros[32] = {0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                          0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0};
    // for statistics
    uint64_t nelems_read = 0;
    uint64_t nelems_container = 0;
    uint64_t points_read = 0;
    int nreads_performed = 0;
    double t_span = 0.0, t_read = 0.0, t_pick = 0.0;
    double te, ttemp;
    double tb = MPI_Wtime();


    v = bp_find_var_byid (fh, r->varid);
    size_of_type = bp_get_type_size (v->type, v->characteristics [0].value);

    /* The idea is we convert a point selection to bounding box or writeblock selection. */
    read_request nrs;
    read_request *nr = &nrs;

    nr->varid = r->varid;
    nr->priv = r->priv;

    ADIOS_SELECTION * container = sel->u.points.container_selection;
    if (container &&
        container->type != ADIOS_SELECTION_WRITEBLOCK &&
        container->type != ADIOS_SELECTION_BOUNDINGBOX)
    {
        log_error ("The container selection of a point selection "
                "must be either NULL or a bounding box or a writeblock. "
                "We will use the default global bounding box\n");
        container = NULL;
    }

    int free_container = 0;
    if (!container)
    {
        // create full bounding box as container
        ADIOS_VARINFO * vinfo = adios_read_bp_inq_var_byid (fp, r->varid);
        container = a2sel_boundingbox (vinfo->ndim, zeros, vinfo->dims);
        free_container = 1;
        sel->u.points.container_selection = container; // save the container instead of NULL to be used in pick_points_from_boundingbox
        common_read_free_varinfo (vinfo);
    }

    char * dest = (char*)(r->data); // pointer to output buffer, will increase
        // continuously as we pick points

    if (container->type == ADIOS_SELECTION_WRITEBLOCK)
    {
        // Read the writeblock and then pick the points from memory
        ADIOS_SELECTION * wbsel = container;
        ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb = &wbsel->u.block;
        uint64_t ldims[32], gdims[32], offsets[32];

        for (step = 0; step < r->nsteps; step++)
        {
            // Get the size of the writeblock
            //  but first prepare fake request to call adios_wbidx_to_pgidx
            ttemp = MPI_Wtime();
            nr->from_steps = r->from_steps;
            nr->nsteps = r->nsteps;
            nr->sel = wbsel;
            int wbidx = wb->is_absolute_index && !p->streaming ?
                        wb->index :
                        adios_wbidx_to_pgidx (fp, nr, step);

            uint64_t nelems = 1;
            /* get ldims for the chunk and then calculate payload size */
            bndim = v->characteristics[wbidx].dims.count;
            bp_get_dimension_characteristics_notime(&(v->characteristics[wbidx]), ldims, gdims, offsets, is_fortran_file (fh));
            for (d = 0; d < bndim; d++)
            {
                nelems *= ldims [d];
            }
            t_span += MPI_Wtime() - ttemp;

            // Read the writeblock
            nr->data = malloc (nelems*size_of_type);
            if (nr->data)
            {
                nr->from_steps = r->from_steps+step;
                nr->nsteps = 1;
                nr->datasize = nelems*size_of_type;
                nr->sel = wbsel;
                ttemp = MPI_Wtime();
                ADIOS_VARCHUNK * wbchunk = read_var_wb (fp, nr);
                t_read += MPI_Wtime() - ttemp;

                // Get the points
                ttemp = MPI_Wtime();
                uint64_t np = pick_points_from_boundingbox (sel, size_of_type, bndim, zeros, ldims,
                                              nelems, zeros, ldims, nr->data, dest);
                if (np < sel->u.points.npoints) {
                    // for each step we need to move npoints, so fill up the output array
                    memset (dest, 0, (sel->u.points.npoints-np)*size_of_type);
                }
                t_pick += MPI_Wtime() - ttemp;
                points_read += np;
                dest += sel->u.points.npoints * size_of_type;

                free (nr->data);
                common_read_free_chunk (wbchunk);
                nelems_read += nelems;
                nelems_container += nelems;
                nreads_performed++;
            }
            else
            {
                adios_error (err_no_memory, "Could not allocate memory to read %" PRIu64
                        "bytes of writeblock %d\n", nelems_read*size_of_type, wbidx);
                memset (dest, 0, sel->u.points.npoints * size_of_type);
            }
        }
    }
    else if (container->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        bndim = container->u.bb.ndim;
        assert (pndim == bndim || pndim == 1);
        ADIOS_SELECTION * nsel;
        // make one or many bounding boxes around the points
        nsel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
        assert (nsel);

        nr->sel = nsel;
        nr->from_steps = r->from_steps;
        nr->nsteps = 1;
        nr->datasize  = size_of_type;

        nsel->type = ADIOS_SELECTION_BOUNDINGBOX;
        nsel->u.bb.ndim = bndim;
        nsel->u.bb.start = (uint64_t *) malloc (nsel->u.bb.ndim * sizeof(uint64_t));
        nsel->u.bb.count = (uint64_t *) malloc (nsel->u.bb.ndim * sizeof(uint64_t));
        assert (nsel->u.bb.start && nsel->u.bb.count);

        if ((sel->u.points.npoints < 5) || (r->nsteps > 1))
        {
            // few points to read, read them one by one
            for (d = 0; d < nsel->u.bb.ndim; d++)
            {
                nsel->u.bb.count[d] = 1;
            }

            nr->data = dest;

            for (step = 0; step < r->nsteps; step++)
            {
                for (i = 0; i < sel->u.points.npoints; i++)
                {
                    ttemp = MPI_Wtime();
                    if (pndim == bndim)
                    {
                        for (d = 0; d < bndim; d++)
                        {
                            nsel->u.bb.start[d] = sel->u.points.points[d + i*pndim] + container->u.bb.start[d];
                        }
                    } else
                    {
                        a2sel_points_1DtoND_box (1, &sel->u.points.points[i], bndim, container->u.bb.start,
                                                                container->u.bb.count, 1, nsel->u.bb.start);
                    }
                    t_pick += MPI_Wtime() - ttemp;
                    //memcpy (nsel->u.bb.start, sel->u.points.points + i * sel->u.points.ndim, sel->u.points.ndim * 8);
                    ttemp = MPI_Wtime();
                    chunk = read_var_bb (fp, nr);
                    t_read += MPI_Wtime() - ttemp;
                    nreads_performed++;
                    nr->data = (char *) nr->data + size_of_type;
                    common_read_free_chunk (chunk);
                }
                nr->from_steps++;
                nelems_container += adios_get_nelements_of_box (bndim, container->u.bb.start, container->u.bb.count);
                nelems_read += sel->u.points.npoints;
            }
            //dest += r->nsteps * sel->u.points.npoints * size_of_type; // == nr->data
        }
        else
        {
            /* Points in a Bounding Box container */
            // Create smaller bounding boxes and read them, then pick the points in memory
            uint64_t start[bndim], count[bndim];
            nelems_container = adios_get_nelements_of_box (bndim, container->u.bb.start, container->u.bb.count);

            ttemp = MPI_Wtime();
            if (bndim == pndim) {
                // get the smallest box the points fit in
                mGetPointlistSpan(&(sel->u.points), start, count);
                // adjust this box (relative to (0,0)) to the original box
                for (d = 0; d < bndim; d++) {
                    start[d] += container->u.bb.start[d];
                }
            } else {
                mGetPointlistSpan1D(&(sel->u.points), bndim, container->u.bb.start, container->u.bb.count, start, count);
                //memcpy (start, container->u.bb.start, bndim*sizeof(uint64_t));
                //memcpy (count, container->u.bb.count, bndim*sizeof(uint64_t));
            }
            t_span += MPI_Wtime() - ttemp;

            uint64_t nelems = adios_get_nelements_of_box (bndim, start, count);

            // allocate temporary array for reading block data
            int max_buffersize = 268435456; // 256MB max read
            if (nelems * size_of_type < max_buffersize)
                max_buffersize = nelems * size_of_type;
            void * readbuf = (void *) malloc (max_buffersize);
            while (!readbuf && max_buffersize > 1024)
            {
                max_buffersize = max_buffersize / 2;
                readbuf = (void *) malloc (nelems);
            }

            log_debug ("Allocated read buffer size = %d bytes (%d elements)\n", max_buffersize, max_buffersize/size_of_type);

            // we read maximum 'maxreadn' elements at once
            uint64_t maxreadn = max_buffersize/size_of_type;
            if (nelems < maxreadn)
                maxreadn = nelems;

            // determine strategy how to read in:
            //  - at once
            //  - loop over 1st dimension
            //  - loop over 1st & 2nd dimension
            //  - etc
            log_debug ("Read size strategy: \n");
            uint64_t readn[bndim];
            uint64_t sum = (uint64_t) 1;
            uint64_t actualreadn = (uint64_t) 1;
            for (d=bndim-1; d>=0; d--)
            {
                if (sum >= maxreadn) {
                    readn[d] = 1;
                } else {
                    readn[d] = maxreadn / sum; // sum is small for 4 bytes here
                    // this may be over the max count for this dimension
                    if (readn[d] > count[d])
                        readn[d] = count[d];
                }
                log_debug ("    dim %d: read %" PRId64 " elements\n", d, readn[d]);
                sum = sum * count[d];
                actualreadn = actualreadn * readn[d];
            }
            log_debug ("    read %" PRId64 " elements at once, %" PRId64 " in total (nelems=%" PRId64 ")\n", actualreadn, sum, nelems);

            uint64_t *s = nsel->u.bb.start;
            uint64_t *c = nsel->u.bb.count; // for block reading of smaller chunks
            // init s and c
            memcpy (s, start, bndim * 8);
            memcpy (c, readn, bndim * 8);

            nr->data = readbuf;
            if (nr->data)
            {
                // read until we have read all 'nelems' elements
                sum = 0;
                while (sum < nelems) {

                    // how many elements do we read in next?
                    actualreadn = 1;
                    for (d=0; d<bndim; d++)
                        actualreadn *= c[d];
                    /*
                    log_debug ("Read the next block ");
                    PRINT_DIMS64("  start", s, tdims, j);
                    PRINT_DIMS64("  count", c, tdims, j);
                    log_debug ("  read %d elems\n", actualreadn);
                     */

                    ttemp = MPI_Wtime();
                    // read a slice finally
                    chunk = read_var_bb (fp, nr);
                    t_read += MPI_Wtime() - ttemp;

                    if (chunk) {
                        // Get the points
                        ttemp = MPI_Wtime();
                        uint64_t npoints = pick_points_from_boundingbox (sel, size_of_type,
                                bndim, container->u.bb.start, container->u.bb.count,
                                actualreadn, s, c, nr->data, dest);

                        dest += npoints * size_of_type;
                        common_read_free_chunk (chunk);
                        nreads_performed++;
                        t_pick += MPI_Wtime() - ttemp;
                    }

                    // prepare for next read
                    sum += actualreadn;
                    int incdim=1; // largest dim should be increased
                    for (d=bndim-1; d>=0; d--) {
                        if (incdim) {
                            if (s[d]+c[d] == start[d]+count[d]) {
                                // reached the end of this dimension
                                s[d] = start[d];
                                c[d] = readn[d];
                                incdim = 1; // next smaller dim can increase too
                            } else {
                                // move up in this dimension up to total count
                                s[d] += readn[d];
                                if (s[d]+c[d] >             start[d]+count[d]) {
                                    // do not reach over the limit
                                    c[d] = start[d]+count[d]-s[d];
                                }
                                incdim = 0;
                            }
                        }
                    }
                } // end while sum < nelems
                nelems_read += nelems;
                free(readbuf);
            }
            else
            {
                adios_error (err_no_memory, "Could not allocate memory to read %" PRIu64
                        "bytes of a bounding box\n", maxreadn*size_of_type);
                memset (dest, 0, sel->u.points.npoints * size_of_type);
            }
        }

        a2sel_free (nsel);
    }

    if (nerr > 0)
    {
        adios_error(err_out_of_bound, "%" PRIu64 " points were out of bound of the subselection\n", nerr);
    }

    chunk = (ADIOS_VARCHUNK *) malloc (sizeof (ADIOS_VARCHUNK));
    assert (chunk);

    chunk->varid = r->varid;
    chunk->type = v->type;
    chunk->from_steps = r->from_steps;
    chunk->nsteps = r->nsteps;
    chunk->sel = a2sel_copy (r->sel);
    chunk->data = r->data;

    if (free_container)
    {
        a2sel_free(container);
        sel->u.points.container_selection = NULL;
    }

    te = MPI_Wtime();
    log_info ("Point selection reading: number of points = %" PRIu64
            ", total container elements = %" PRIu64
            ", number of reads = %d"
            ", elements read from file = %" PRIu64
            "\n       Time total = %6.2fs  span = %6.2fs  read = %6.2fs  picking = %6.2fs\n",
            sel->u.points.npoints*r->nsteps, nelems_container, nreads_performed, nelems_read,
            te-tb, t_span, t_read, t_pick);
    return chunk;
}

int adios_read_bp_init_method (MPI_Comm comm, PairStruct * params)
{
    int  max_chunk_size, pollinterval;
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
            pollinterval = strtol(p->value, NULL, 10);
            if (pollinterval > 0 && !errno)
            {
                log_debug ("poll_interval set to %d secs for READ_BP read method\n",
                            pollinterval);
                poll_interval_msec = pollinterval;
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
    poll_interval_msec = 10000; // 10 secs by default
    show_hidden_attrs = 0; // don't show hidden attr by default

    return 0;
}

static int open_stream (ADIOS_FILE * fp, const char * fname,
                        MPI_Comm comm, float timeout_sec)
{
    int rank;
    BP_PROC * p;
    BP_FILE * fh;
    int stay_in_poll_loop = 1;
    int file_ok = 0;
    double t1 = adios_gettime_double();

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
            adios_errno = err_no_error; // clear previous intermittent error
            file_ok = check_bp_validity (fname);
            if (!file_ok)
            {
                // This stream does not exist yet
                log_debug ("file %s is not a valid file for streaming read."
                           "One possible reason is it's a VERY old BP file,"
                           "which doesn't allow reader to check its validity.\n", fname);

                if (stay_in_poll_loop)
                {
                    // check if we need to stay in loop
                    if (timeout_sec == 0.0)  //return immediately, which means check file once
                    {
                        stay_in_poll_loop = 0;
                    }
                    else if (timeout_sec < 0.0) // check file until it arrives
                    {
                        adios_nanosleep (poll_interval_msec/1000,
                            (int)(((uint64_t)poll_interval_msec * 1000000L)%1000000000L));
                        stay_in_poll_loop = 1;
                    }
                    else if (timeout_sec > 0.0 && (adios_gettime_double () - t1 > timeout_sec))
                    {
                        stay_in_poll_loop = 0;
                    }
                    else
                    {
                        adios_nanosleep (poll_interval_msec/1000,
                            (int)(((uint64_t)poll_interval_msec * 1000000L)%1000000000L));
                    }
                }
            }
            else
            {
                stay_in_poll_loop = 0;
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

    fh = BP_FILE_alloc (fname, comm);

    p = (BP_PROC *) malloc (sizeof (BP_PROC));
    assert (p);
    p->fh = fh;
    p->streaming = 1;
    p->varid_mapping = 0;
    p->local_read_request_list = 0;
    p->b = 0;
    p->priv = 0;

    /* BP file open and gp/var/att parsing */
    bp_open (fname, comm, fh);

    fp->fh = (uint64_t) p;
    fp->file_size = fh->mfooter.file_size;
    fp->version = fh->mfooter.version & ADIOS_VERSION_NUM_MASK;
    fp->path = strdup (fh->fname);
    fp->endianness = bp_get_endianness (fh->mfooter.change_endianness);

    /* Seek to the initial step. */
    bp_seek_to_step (fp, 0, show_hidden_attrs);

    fp->current_step = 0;
    /* For file, the last step is tidx_stop */
    fp->last_step = fh->tidx_stop - fh->tidx_start;

    return 0;
}

/* As opposed to open_file, open_stream opens the first step in the file only.
 * The lock_mode for file reading is ignored for now.
 */
ADIOS_FILE * adios_read_bp_open (const char * fname, MPI_Comm comm, enum ADIOS_LOCKMODE lock_mode, float timeout_sec)
{
    log_debug ("adios_read_bp_open\n");

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
    int rank;
    BP_PROC * p;
    BP_FILE * fh;
    ADIOS_FILE * fp;

    log_debug ("adios_read_bp_open_file\n");

    MPI_Comm_rank (comm, &rank);

    fh = BP_FILE_alloc (fname, comm);

    p = (BP_PROC *) malloc (sizeof (BP_PROC));
    assert (p);
    p->fh = fh;
    p->streaming = 0;
    p->varid_mapping = 0; // maps perceived id to real id
    p->local_read_request_list = 0;
    p->b = 0;
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
        adios_error (err_file_open_error, "File open failed: %s\n", fname);
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
    fp->endianness = bp_get_endianness (fh->mfooter.change_endianness);
    fp->version = fh->mfooter.version & ADIOS_VERSION_NUM_MASK;
    fp->file_size = fh->mfooter.file_size;

    return fp;
}

int adios_read_bp_close (ADIOS_FILE * fp)
{
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);

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
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);
    int last_tidx;
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
            last_tidx = fh->tidx_stop;
            fname = strdup (fh->fname);
            comm = fh->comm;

            if (p->fh)
            {
                bp_close (fh);
                p->fh = 0;
            }

            if (!get_new_step (fp, fname, comm, last_tidx, timeout_sec))
            {
                // With file reading, how can we tell it is the end of the streams?
                adios_errno = err_step_notready;
            }

            free (fname);

            if (adios_errno == 0)
            {
                release_step (fp);
                bp_seek_to_step (fp, fp->last_step + 1, show_hidden_attrs);
                fp->current_step = fp->last_step + 1;
            }
        }
    }
    else // read in newest step. Re-open no matter whether current_step < last_step
    {
        last_tidx = fh->tidx_stop;
        fname = strdup (fh->fname);
        comm = fh->comm;

        if (p->fh)
        {
            bp_close (fh);
            p->fh = 0;
        }

        // lockmode is currently not supported.
        if (!get_new_step (fp, fh->fname, fh->comm, last_tidx, timeout_sec))
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
    ADIOS_VARINFO * varinfo;

    int mapped_id = map_req_varid (fp, varid);;

    adios_errno = 0;

    /* this call sets varinfo->varid as the real mapped id.
     * Therefore, we need to set it back to perceived id.
     */
    varinfo = bp_inq_var_byid (fp, mapped_id);
    varinfo->varid = varid;

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
    int size, sum_size, nsteps, prev_timestep;
    int nb; // total number of blocks (varinfo->sum_nblocks)
    BP_FILE * fh = GET_BP_FILE (fp);

    ADIOS_VARSTAT * vs;
    struct adios_index_var_struct_v1 * var_root;

    assert (varinfo);

    varinfo->statistics = vs = (ADIOS_VARSTAT *) malloc (sizeof (ADIOS_VARSTAT));
    assert (vs);

    vs->min = NULL;
    vs->max = NULL;
    vs->avg = NULL;
    vs->std_dev = NULL;

    if (per_step_stat) {
        vs->steps = (struct ADIOS_STAT_STEP *) malloc (sizeof (struct ADIOS_STAT_STEP));
        assert (vs->steps);
        vs->steps->mins = NULL;
        vs->steps->maxs = NULL;
        vs->steps->avgs = NULL;
        vs->steps->std_devs = NULL;
    } else {
        vs->steps = NULL;
    }

    if (per_block_stat) {
        vs->blocks = (struct ADIOS_STAT_BLOCK *) malloc (sizeof (struct ADIOS_STAT_BLOCK));
        assert (vs->blocks);
        vs->blocks->mins = NULL;
        vs->blocks->maxs = NULL;
        vs->blocks->avgs = NULL;
        vs->blocks->std_devs = NULL;
    } else {
        vs->blocks = NULL;
    }

    //TODO
    vs->histogram = NULL;

    uint64_t gcnt = 0, *cnts=NULL, *bcnts = NULL;

    double *gsum = NULL, *gsum_square = NULL;
    double **sums = NULL,  **sum_squares = NULL;
    double **bsums = NULL, **bsum_squares = NULL;

    int16_t map[32];
    memset (map, -1, sizeof(map));

    nsteps = varinfo->nsteps;
    nb = varinfo->sum_nblocks;
    int time = adios_step_to_time (fp, varinfo->varid, 0);

    var_root = bp_find_var_byid (fh, varinfo->varid);
    // will loop from var_root->characteristics[from_ch..to_ch-1]
    int from_ch = 0; 
    int to_ch = var_root->characteristics_count;

    if (fp->is_streaming) {
        // find the current timestep in characteristics array, because when streaming from a file,
        // var_root.characteristics contains many timesteps
        from_ch = get_var_start_index (var_root, time); 
        to_ch = from_ch + nb;
        assert(from_ch < var_root->characteristics_count);
        assert(to_ch  <= var_root->characteristics_count);
    }

    // Bitmap shows which statistical information has been calculated
    i = j = 0;
    while (var_root->characteristics[from_ch].bitmap >> j)
    {
        if ((var_root->characteristics[from_ch].bitmap >> j) & 1)
        {
            map [j] = i ++;
        }

        j ++;
    }

    if (map[adios_statistic_min] != -1)
    {
        if (per_step_stat) {
            MALLOC(vs->steps->mins, nsteps * sizeof(void *), "minimum per timestep");
            for (i = 0; i < nsteps; i++)
            {
                vs->steps->mins[i] = NULL;
            }
        }
        if (per_block_stat) {
            MALLOC(vs->blocks->mins, nb * sizeof(void *), "minimum per writeblock");
            for (i = 0; i < nb; i++)
            {
                vs->blocks->mins[i] = NULL;
            }
        }
    }

    if (map[adios_statistic_max] != -1)
    {
        if (per_step_stat) {
            MALLOC(vs->steps->maxs, nsteps * sizeof(void *), "maximum per timestep");
            for (i = 0; i < nsteps; i++)
            {
                vs->steps->maxs[i] = NULL;
            }
        }
        if (per_block_stat) {
            MALLOC(vs->blocks->maxs, nb * sizeof(void *), "maximum per writeblock");
            for (i = 0; i < nb; i++)
            {
                vs->blocks->maxs[i] = NULL;
            }
        }
    }

    if (map[adios_statistic_sum] != -1)
    {
        if (per_step_stat) {
            MALLOC(sums, nsteps * sizeof(double *), "summation per timestep");
            MALLOC(vs->steps->avgs, nsteps * sizeof(double *), "average per timestep");

            for (i = 0; i < nsteps; i++)
            {
                sums[i] = vs->steps->avgs[i] = NULL;
            }
            CALLOC(cnts, nsteps, sizeof(uint64_t), "count of elements per timestep");
        }
        if (per_block_stat) {
            MALLOC(bsums, nb * sizeof(double *), "summation per writeblock");
            MALLOC(vs->blocks->avgs, nb * sizeof(double *), "average per writeblock");

            for (i = 0; i < nb; i++)
            {
                bsums[i] = vs->blocks->avgs[i] = NULL;
            }
            CALLOC(bcnts, nb, sizeof(uint64_t), "count of elements per writeblock");
        }
    }

    if (map[adios_statistic_sum_square] != -1)
    {
        if (per_step_stat) {
            MALLOC(sum_squares, nsteps * sizeof(double *), "summation per timestep");
            MALLOC(vs->steps->std_devs, nsteps * sizeof(double *), "standard deviation per timestep");

            for (i = 0; i < nsteps; i++)
            {
                vs->steps->std_devs[i] = sum_squares[i] = NULL;
            }
        }
        if (per_block_stat) {
            MALLOC(bsum_squares, nb * sizeof(double *), "summation per writeblock");
            MALLOC(vs->blocks->std_devs, nb * sizeof(double *), "standard deviation per writeblock");

            for (i = 0; i < nb; i++)
            {
                vs->blocks->std_devs[i] = bsum_squares[i] = NULL;
            }
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
    enum ADIOS_DATATYPES original_var_type = var_root->type;

    if (var_root->characteristics[from_ch].transform.transform_type != adios_transform_none) {
        original_var_type = var_root->characteristics[from_ch].transform.pre_transform_type;
    }

    size = bp_get_type_size (original_var_type, "");
    sum_size = bp_get_type_size (adios_double, "");

    if (original_var_type == adios_complex || original_var_type == adios_double_complex)
    {
        int type;
        int idx; // block array index, = i-from_ch in loops
        int tidx; // step array index, = timestep-fp->current_step in loops
        count = 3;
        timestep = fp->current_step; // 0..
        prev_timestep = time; // 1..

        if (original_var_type == adios_complex)
        {
            type = adios_double;
        }
        else
        {
            type = adios_long_double;
        }

        // Only a double precision returned for all complex values
        size = bp_get_type_size (adios_double, "");

        for (i = from_ch; i < to_ch; i++)
        {
            // changes for 1.4.x. Q. Liu
            if (var_root->characteristics[i].time_index != prev_timestep)
            {
                timestep++;
                prev_timestep = var_root->characteristics[i].time_index;
            }

            idx = i - from_ch;
            tidx = timestep - fp->current_step;

            assert (tidx < nsteps);
            assert (idx < nb);

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

                if(!vs->min) {
                    MALLOC (vs->min, count * size, "global minimum")
                    for (c = 0; c < count; c ++)
                           ((double * ) vs->min)[c] = data[c];

                } else {
                    for (c = 0; c < count; c ++)
                        if (data[c] < ((double *) vs->min)[c])
                               ((double * ) vs->min)[c] = data[c];
                }

                if (per_step_stat) {
                    if(!vs->steps->mins[tidx]) {
                        MALLOC (vs->steps->mins[tidx], count * size, "minimum per timestep")
                        for (c = 0; c < count; c ++) {
                            ((double **) vs->steps->mins)[tidx][c] = data[c];
                        }
                    } else {
                        for (c = 0; c < count; c ++) {
                            if (data[c] < ((double **) vs->steps->mins)[tidx][c]) {
                                ((double **) vs->steps->mins)[tidx][c] = data[c];
                            }
                        }
                    }
                }

                if (per_block_stat) {
                    if(!vs->blocks->mins[idx]) {
                        MALLOC (vs->blocks->mins[idx], count * size, "minimum per writeblock")
                        for (c = 0; c < count; c ++) {
                            ((double **) vs->blocks->mins)[idx][c] = data[c];
                        }
                    } else {
                        for (c = 0; c < count; c ++) {
                            if (data[c] < ((double **) vs->blocks->mins)[idx][c]) {
                                ((double **) vs->blocks->mins)[idx][c] = data[c];
                            }
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

                if (per_step_stat) {
                    if(!vs->steps->maxs[tidx]) {
                        MALLOC (vs->steps->maxs[tidx], count * size, "maximum per timestep")
                        for (c = 0; c < count; c ++)
                            ((double **) vs->steps->maxs)[tidx][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            if (data[c] > ((double **) vs->steps->maxs)[tidx][c])
                                ((double **) vs->steps->maxs)[tidx][c] = data[c];
                    }
                }

                if (per_block_stat) {
                    if(!vs->blocks->maxs[idx]) {
                        MALLOC (vs->blocks->maxs[idx], count * size, "maximum per writeblock")
                        for (c = 0; c < count; c ++)
                            ((double **) vs->blocks->maxs)[idx][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            if (data[c] > ((double **) vs->blocks->maxs)[idx][c])
                                ((double **) vs->blocks->maxs)[idx][c] = data[c];
                    }
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

                if (per_step_stat) {
                    if(!sums[tidx]) {
                        MALLOC(sums[tidx], count * sum_size, "summation per timestep")
                        for (c = 0; c < count; c ++)
                            sums[tidx][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            sums[tidx][c] = sums[tidx][c] + data[c];
                    }
                }

                if (per_block_stat) {
                    if(!bsums[idx]) {
                        MALLOC(bsums[idx], count * sum_size, "summation per writeblock")
                        for (c = 0; c < count; c ++)
                            bsums[idx][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            bsums[idx][c] = bsums[idx][c] + data[c];
                    }
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

                if (per_step_stat) {
                    if(!sum_squares[tidx]) {
                        MALLOC(sum_squares[tidx], count * sum_size, "summation of square per timestep")
                        for (c = 0; c < count; c ++)
                            sum_squares[tidx][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            sum_squares[tidx][c] = sum_squares[tidx][c] + data[c];
                    }
                }

                if (per_block_stat) {
                    if(!bsum_squares[idx]) {
                        MALLOC(bsum_squares[idx], count * sum_size, "summation of square per writeblock")
                        for (c = 0; c < count; c ++)
                            bsum_squares[idx][c] = data[c];

                    } else {
                        for (c = 0; c < count; c ++)
                            bsum_squares[idx][c] = bsum_squares[idx][c] + data[c];
                    }
                }
            }

            if (map[adios_statistic_cnt] != -1 && stats[0][map[adios_statistic_cnt]].data)
            {
                if (per_step_stat) {
                    cnts[tidx] += * ((uint32_t *) stats[0][map[adios_statistic_cnt]].data);
                }
                if (per_block_stat) {
                    bcnts[idx] += * ((uint32_t *) stats[0][map[adios_statistic_cnt]].data);
                }
                gcnt += * (uint32_t *) stats[0][map[adios_statistic_cnt]].data;
            }
        }

        if (per_step_stat) {
            if(vs->min && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) {
                // min, max, summation exists only for arrays
                // Calculate average / timestep

                for(timestep = 0; timestep < nsteps; timestep ++) {
                    MALLOC(vs->steps->avgs[timestep], count * sum_size, "average per timestep")
                    for (c = 0; c < count; c ++)
                        vs->steps->avgs[timestep][c] = sums[timestep][c] / cnts[timestep];

                    MALLOC(vs->steps->std_devs[timestep], count * sum_size, "standard deviation per timestep")
                    for (c = 0; c < count; c ++)
                        vs->steps->std_devs[timestep][c] = 
                            sqrt((sum_squares[timestep][c] / cnts[timestep]) - 
                            (vs->steps->avgs[timestep][c] * vs->steps->avgs[timestep][c]));

                    free (sums[timestep]);
                    free (sum_squares[timestep]);
                }
            }
        }

        // Calculate per-block average 
        if (per_block_stat) {
            if(vs->min && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1)) 
            {
                for (i = 0; i < to_ch-from_ch; i++)
                {
                    MALLOC(vs->blocks->avgs[i], count * sum_size, "average per writeblock")
                    MALLOC(vs->blocks->std_devs[i], count * sum_size, "standard deviation per writeblock")

                    if(bcnts[i]) {
                        for (c = 0; c < count; c ++)
                            vs->blocks->avgs[i][c] = bsums[i][c] / bcnts[i];

                        for (c = 0; c < count; c ++)
                            vs->blocks->std_devs[i][c] = 
                                sqrt((bsum_squares[i][c] / bcnts[i]) - 
                                        (vs->blocks->avgs[i][c] * vs->blocks->avgs[i][c]));
                    } else {
                        // this block is an empty block (0 size) in file
                        for (c = 0; c < count; c ++) {
                            vs->blocks->avgs[i][c] = 0.0;
                            vs->blocks->std_devs[i][c] = 0.0;
                        }
                    }
                    free (bsums[i]);
                    free (bsum_squares[i]);
                }
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
        int idx; // array index, = i-from_ch in loops
        int tidx; // step array index, = timestep-fp->current_step in loops
        timestep = fp->current_step; // 0..
        prev_timestep = time; // 1..

        for (i = from_ch; i < to_ch; i++)
        {
            //printf ("i = %3d, time_index = %d, prev = %d, count = %lld\n", 
            //   i,  var_root->characteristics[i].time_index, prev_timestep, var_root->characteristics_count);

            // changes for 1.4.x. Q. Liu
            if (var_root->characteristics[i].time_index != prev_timestep)
            {
                timestep++;
                prev_timestep = var_root->characteristics[i].time_index;
            }

            idx = i - from_ch;
            tidx = timestep - fp->current_step;
            //printf ("        tidx = %d, idx = %d, nsteps = %d, nb=%d\n", tidx, idx, nsteps, nb); 

            assert (tidx < nsteps);
            assert (idx < nb);

            if (!var_root->characteristics[i].stats)
            {
                continue;
            }

            struct adios_index_characteristics_stat_struct * stats = var_root->characteristics[i].stats[0];

            if (map[adios_statistic_finite] != -1 && (* ((uint8_t *) stats[map[adios_statistic_finite]].data) == 0))
                continue;

            if (map[adios_statistic_min] != -1 && stats[map[adios_statistic_min]].data)
            {
                if(!vs->min)
                {
                    MALLOC (vs->min, size, "global minimum")
                    memcpy(vs->min, stats[map[adios_statistic_min]].data, size);

                }
                else if (adios_lt(original_var_type, stats[map[adios_statistic_min]].data, vs->min))
                {
                    memcpy(vs->min, stats[map[adios_statistic_min]].data, size);
                }

                if (per_step_stat) {
                    if(!vs->steps->mins[tidx])
                    {
                        MALLOC (vs->steps->mins[tidx], size, "minimum per timestep")
                        memcpy(vs->steps->mins[tidx], stats[map[adios_statistic_min]].data, size);
                    }
                    else if (adios_lt(original_var_type, stats[map[adios_statistic_min]].data, vs->steps->mins[tidx]))
                    {
                        memcpy(vs->steps->mins[tidx], stats[map[adios_statistic_min]].data, size);
                    }
                }

                if (per_block_stat) {
                    if(!vs->blocks->mins[idx])
                    {
                        MALLOC (vs->blocks->mins[idx], size, "minimum per writeblock")
                        memcpy(vs->blocks->mins[idx], stats[map[adios_statistic_min]].data, size);
                    }
                    else if (adios_lt(original_var_type, stats[map[adios_statistic_min]].data, vs->blocks->mins[idx]))
                    {
                        memcpy(vs->blocks->mins[idx], stats[map[adios_statistic_min]].data, size);
                    }
                }
            }

            if (map[adios_statistic_max] != -1 && stats[map[adios_statistic_max]].data)
            {
                if(!vs->max)
                {
                    MALLOC (vs->max, size, "global maximum")
                    memcpy(vs->max, stats[map[adios_statistic_max]].data, size);

                }
                else if (adios_lt(original_var_type, vs->max, stats[map[adios_statistic_max]].data))
                {
                    memcpy(vs->max, stats[map[adios_statistic_max]].data, size);
                }

                if (per_step_stat) {
                    if(!vs->steps->maxs[tidx])
                    {
                        MALLOC (vs->steps->maxs[tidx], size, "maximum per timestep")
                        memcpy(vs->steps->maxs[tidx], stats[map[adios_statistic_max]].data, size);
                    }
                    else if (adios_lt(original_var_type, vs->steps->maxs[tidx], stats[map[adios_statistic_max]].data))
                    {
                        memcpy(vs->steps->maxs[tidx], stats[map[adios_statistic_max]].data, size);
                    }
                }

                if (per_block_stat) {
                    if(!vs->blocks->maxs[idx])
                    {
                        MALLOC (vs->blocks->maxs[idx], size, "maximum per writeblock")
                        memcpy(vs->blocks->maxs[idx], stats[map[adios_statistic_max]].data, size);
                    }
                    else if (adios_lt(original_var_type, stats[map[adios_statistic_max]].data, vs->blocks->maxs[idx]))
                    {
                        memcpy(vs->blocks->maxs[idx], stats[map[adios_statistic_max]].data, size);
                    }
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

                if (per_step_stat) {
                    if(!sums[tidx])
                    {
                        MALLOC(sums[tidx], sum_size, "summation per timestep")
                        memcpy(sums[tidx], stats[map[adios_statistic_sum]].data, sum_size);
                    }
                    else
                    {
                        *sums[tidx] = *sums[tidx] + * ((double *) stats[map[adios_statistic_sum]].data);
                    }
                }

                if (per_block_stat) {
                    if(!bsums[idx])
                    {
                        MALLOC(bsums[idx], sum_size, "summation per writeblock")
                        memcpy(bsums[idx], stats[map[adios_statistic_sum]].data, sum_size);
                    }
                    else
                    {
                        *bsums[idx] = *bsums[idx] + * ((double *) stats[map[adios_statistic_sum]].data);
                    }
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

                if (per_step_stat) {
                    if(!sum_squares[tidx])
                    {
                        MALLOC(sum_squares[tidx], sum_size, "summation of square per timestep")
                        memcpy(sum_squares[tidx], stats[map[adios_statistic_sum_square]].data, sum_size);
                    }
                    else
                    {
                        *sum_squares[tidx] = *sum_squares[tidx] + * ((double *) stats[map[adios_statistic_sum_square]].data);
                    }
                }

                if (per_block_stat) {
                    if(!bsum_squares[idx])
                    {
                        MALLOC(bsum_squares[idx], sum_size, "summation of square per writeblock")
                        memcpy(bsum_squares[idx], stats[map[adios_statistic_sum_square]].data, sum_size);
                    }
                    else
                    {
                        *bsum_squares[idx] = *bsum_squares[idx] + * ((double *) stats[map[adios_statistic_sum_square]].data);
                    }
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
                        vi->hist->frequenciess[tidx][j] += freq;
                }
            }
*/
            if (map[adios_statistic_cnt] != -1 && stats[map[adios_statistic_cnt]].data)
            {
                if (per_step_stat) {
                    cnts[tidx] += * (uint32_t *) stats[map[adios_statistic_cnt]].data;
                }
                if (per_block_stat) {
                    bcnts[idx] = * (uint32_t *) stats[map[adios_statistic_cnt]].data;
                }
                gcnt += * (uint32_t *) stats[map[adios_statistic_cnt]].data;
            }
        }

        if (per_step_stat) {
            if(nsteps > 0 && vs->min
                    && (map[adios_statistic_sum] != -1)
                    && (map[adios_statistic_sum_square] != -1)
              )
            {
                // min, max, summation exists only for arrays
                // Calculate average / timestep
                for(timestep = 0; timestep < nsteps; timestep ++)
                {
                    MALLOC(vs->steps->avgs[timestep], sum_size, "average per timestep")
                    if(cnts[timestep]) {
                        *(vs->steps->avgs[timestep]) = *(sums[timestep]) / cnts[timestep];
                    } else {
                        // no summation for this timestep (e.g. constant NAN array)
                        *(vs->steps->avgs[timestep]) = 0.0;
                    }

                    MALLOC(vs->steps->std_devs[timestep], sum_size, "standard deviation per timestep")
                    if(cnts[timestep]) {
                        *(vs->steps->std_devs[timestep]) = sqrt(*(sum_squares[timestep]) / cnts[timestep]
                                - ((*(vs->steps->avgs[timestep]) * (*(vs->steps->avgs[timestep])))));
                    } else {
                        *(vs->steps->std_devs[timestep]) = 0.0;
                    }

                    free (sums[timestep]);
                    free (sum_squares[timestep]);
                }
            }
        }

        // Calculate per-block  average
        if (per_block_stat) {
            if(vs->min && (map[adios_statistic_sum] != -1) && (map[adios_statistic_sum_square] != -1))
            {
                for (i = 0; i < to_ch-from_ch; i++)
                {
                    MALLOC(vs->blocks->avgs[i], sum_size, "average per writeblock")
                    MALLOC(vs->blocks->std_devs[i], sum_size, "standard deviation per writeblock")
                    if(bcnts[i]) {
                        *(vs->blocks->avgs[i]) = *(bsums[i]) / bcnts[i];

                        *(vs->blocks->std_devs[i]) = sqrt(*(bsum_squares[i]) / bcnts[i]
                                    - ((*(vs->blocks->avgs[i]) * (*(vs->blocks->avgs[i])))));
                    } else {
                        // this block is an empty block (0 size) in file
                        *(vs->blocks->avgs[i]) = 0.0;
                        *(vs->blocks->std_devs[i]) = 0.0;
                    }

                    free (bsums[i]);
                    free (bsum_squares[i]);
                }

            }
        }

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

    if (sums) free(sums);
    if (bsums) free(bsums);
    if (gsum) free(gsum);

    if (sum_squares) free (sum_squares);
    if (bsum_squares) free (bsum_squares);
    if (gsum_square) free (gsum_square);

    if (cnts) free (cnts);
    if (bcnts) free (bcnts);

    return 0;
}

// NCSU ALACRITY-ADIOS - Factored out VARBLOCK inquiry function to permit sourcing
static ADIOS_VARBLOCK * inq_var_blockinfo(const ADIOS_FILE * fp, const ADIOS_VARINFO * varinfo, int use_pretransform_dimensions) {
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);
    int i, j, file_is_fortran, nblks, time;
    uint64_t * ldims, * gdims, * offsets;
    int dummy = -1;
    struct adios_index_var_struct_v1 * var_root;
    struct bp_index_pg_struct_v1 * pgs = fh->pgs_root;
    ADIOS_VARBLOCK *blockinfo;
    // variables to heuristicly calculate the process_id of a PG in a subfile
    uint32_t current_process_id = pgs->process_id;
    uint32_t deduced_file_index = 0;
    int64_t current_offset = -1;

    assert (varinfo);

    file_is_fortran = is_fortran_file (fh);
    // For file mode: return all blocks info;
    // For streaming mode: return all blocks within the current step
    // 08/14/2014 Q. Liu
    nblks = (p->streaming ? varinfo->nblocks[0] : varinfo->sum_nblocks);

    // Perform variable ID mapping, since the input to this function is user-perceived
    int mapped_id = map_req_varid (fp, varinfo->varid);
    var_root = bp_find_var_byid (fh, mapped_id);

    blockinfo = (ADIOS_VARBLOCK *) malloc (nblks * sizeof (ADIOS_VARBLOCK));
    assert (blockinfo);

    const struct adios_index_characteristic_struct_v1 *root_characteristic = &var_root->characteristics[0];

    // NCSU ALACRITY-ADIOS - Use pre-transform dimensions if instructed to do so
    int dimcount;
    if (use_pretransform_dimensions && root_characteristic->transform.transform_type != adios_transform_none) {
        dimcount = var_root->characteristics[0].transform.pre_transform_dimensions.count;
    } else {
        dimcount = var_root->characteristics[0].dims.count;
    }
    /* dim.count possibily include 'time' dim in it. */
    ldims = (uint64_t *) malloc (dimcount * 8);
    gdims = (uint64_t *) malloc (dimcount * 8);
    offsets = (uint64_t *) malloc (dimcount * 8);
    assert (ldims && gdims && offsets);

    time = adios_step_to_time (fp, varinfo->varid, 0);

    j = 0; 
    for (i = 0; i < nblks; i++)
    {
        int k; /* to save i or j for the process_id determination step below */
        int has_oldschool_time_index = 0; // old BP file with time encoded as dimension
        blockinfo[i].start = (uint64_t *) malloc (dimcount * 8);
        blockinfo[i].count = (uint64_t *) malloc (dimcount * 8);
        assert (blockinfo[i].start && blockinfo[i].count);

        if (!p->streaming)
        {
            // NCSU ALACRITY-ADIOS
            const struct adios_index_characteristic_struct_v1 *blk_characteristic = &var_root->characteristics[i];
        	// Only use pre-transform dimensions if A) pre-transform dimensions were
        	// requested, and B) this varblock is actually transformed. Use normal
        	// dimensions otherwise
            const struct adios_index_characteristic_dims_struct_v1 *blk_dims =
            		use_pretransform_dimensions && blk_characteristic->transform.transform_type != adios_transform_none ?
            				&blk_characteristic->transform.pre_transform_dimensions :
            				&blk_characteristic->dims;

            bp_get_dimension_generic_notime(blk_dims, ldims, gdims, offsets, file_is_fortran, &has_oldschool_time_index);
            k = i;
        }
        else
        {
            while (j < var_root->characteristics_count && var_root->characteristics[j].time_index != time)
            {
                j++;
            }

            if (j < var_root->characteristics_count)
            {
                // NCSU ALACRITY-ADIOS
                const struct adios_index_characteristic_struct_v1 *blk_characteristic = &var_root->characteristics[j];
                // Only use pre-transform dimensions if A) pre-transform dimensions were
            	// requested, and B) this varblock is actually transformed. Use normal
            	// dimensions otherwise
                const struct adios_index_characteristic_dims_struct_v1 *blk_dims =
                		use_pretransform_dimensions && blk_characteristic->transform.transform_type != adios_transform_none ?
                				&blk_characteristic->transform.pre_transform_dimensions :
                				&blk_characteristic->dims;

                bp_get_dimension_generic_notime(blk_dims, ldims, gdims, offsets, file_is_fortran, &has_oldschool_time_index);
                k = j;
                j++;
            }
            else
            {
                // shoudn't be here.
            }
        }

        // NCSU ALACRITY-ADIOS - If a time dimension was removed above, update
        // dimcount so that dimension copy/swapping works below
        //if (dimcount > 0 && ldims[dimcount-1] == 0 && gdims[dimcount-1] != 0)
        if (has_oldschool_time_index && dimcount > 0)
            dimcount--;

        /*Fix: the function above swaps the dimensions to C order in any case. 
         * For Fortran callers, we have to swap it back here */
        if (futils_is_called_from_fortran ())
        {
            swap_order (dimcount, ldims, &dummy);
            swap_order (dimcount, offsets, &dummy);
        }

        memcpy (blockinfo[i].start, offsets, dimcount * 8);
        memcpy (blockinfo[i].count, ldims, dimcount * 8);

        // NCSU ALACRITY-ADIOS - This code was left over in the Transforms branch after the merge. Not sure if it's needed; preserved in case it represented a valid bugfix
//        if (file_is_fortran != futils_is_called_from_fortran())
//        {
//            swap_order (dimcount, blockinfo[i].start, &timedim);
//            swap_order (dimcount, blockinfo[i].count, &timedim);
//        }
        
        /* Find the process ID */
        //blockinfo[i].process_id = (uint32_t)-1;
        /*
        // old routine fine for single bp file (no subfiles)
        while (pgs != NULL && 
               pgs->offset_in_file <= var_root->characteristics[k].offset) 
        {
            current_process_id = pgs->process_id;
            pgs = pgs->next;
        }
        blockinfo[i].process_id = current_process_id;
        blockinfo[i].time_index = var_root->characteristics[k].time_index;
        */

        /* sub-files' PGs start from 0 offset again and again
           unfortunately the pgs don't have info on subfile index, which is 
           only stored in the variable characteristics. 
           Assumption: process_ids go from 0..n, and all the pgs are ordered
           incrementally in subfiles according to process_ids.
           This is true so far by all writing methods.
        */
        if (pgs)
            current_process_id = pgs->process_id;
            // if pgs==NULL, keep the current process id from the last PG

        while (pgs != NULL) {
            if ((int64_t)pgs->offset_in_file <= current_offset) {
                deduced_file_index++;
            }
            if ((int32_t)deduced_file_index > (int32_t)var_root->characteristics[k].file_index) {
                deduced_file_index--; 
                /* pgs and current_offset does not change anymore and we enter the while 
                   loop again for the next block and will increase this counter again */
                break;
            }
            if (deduced_file_index == var_root->characteristics[k].file_index &&
                    pgs->offset_in_file > var_root->characteristics[k].offset) {
                break;
            }
            current_offset = pgs->offset_in_file;
            current_process_id = pgs->process_id;
            pgs = pgs->next;
        }
        blockinfo[i].process_id = current_process_id;
        blockinfo[i].time_index = var_root->characteristics[k].time_index;
    }

    free (ldims);
    free (gdims);
    free (offsets);
    return blockinfo;
}

// NCSU ALACRITY-ADIOS - Delegate to shared VARBLOCK loader
int adios_read_bp_inq_var_blockinfo (const ADIOS_FILE * fp, ADIOS_VARINFO * varinfo)
{
    varinfo->blockinfo = inq_var_blockinfo(fp, varinfo, 0); // 0 -> use true dimensions, not original dimensions
    assert(varinfo->blockinfo);
    return 0;

}

// NCSU ALACRITY-ADIOS - Adding an inq function to get the new transform metadata from storage
ADIOS_TRANSINFO * adios_read_bp_inq_var_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi) {
    BP_FILE * fh = GET_BP_FILE (fp);
    struct adios_index_var_struct_v1 * var_root;
    int file_is_fortran;
    int dummy;
    ADIOS_TRANSINFO *transinfo;
    assert(vi);
    file_is_fortran = is_fortran_file (fh);

    // Perform variable ID mapping, since the input to this function is user-perceived
    int mapped_id = map_req_varid (fp, vi->varid);
    var_root = bp_find_var_byid(fh, mapped_id);
    assert(var_root);

    transinfo = malloc(sizeof(ADIOS_TRANSINFO));

    const struct adios_index_characteristic_transform_struct *transform = &var_root->characteristics[0].transform;

    transinfo->transform_type = transform->transform_type;
    if (transform->transform_type != adios_transform_none) {
        transinfo->orig_type = transform->pre_transform_type;

        // Load orig_ndims/orig_dims using the utility function
        bp_get_and_swap_dimensions_generic (fp, var_root, file_is_fortran,
                                            &transinfo->orig_ndim, &transinfo->orig_dims,
                                            &dummy,
                                            file_is_fortran != futils_is_called_from_fortran(),
                                            1); // 1 -> get based on pre-transform dimensions

        transinfo->orig_global = is_global_array_generic(&var_root->characteristics[0].transform.pre_transform_dimensions);

        transinfo->transform_metadata_len = transform->transform_metadata_len;
        transinfo->transform_metadata = transform->transform_metadata;
        transinfo->should_free_transform_metadata = 0;
    } else {
        transinfo->orig_type = adios_unknown;
        transinfo->orig_ndim = 0;
        transinfo->orig_dims = 0;
        transinfo->orig_global = 0;
        transinfo->transform_metadata_len = 0;
        transinfo->transform_metadata = 0;
        transinfo->should_free_transform_metadata = 0;
    }
    transinfo->orig_blockinfo = 0;
    transinfo->transform_metadatas = 0;

    return transinfo;
}

// NCSU ALACRITY-ADIOS - Adding an inq function to get original (pre-transform) blockinfo for variables from storage
int adios_read_bp_inq_var_trans_blockinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
	assert(fp);
	assert(vi);
	assert(ti);

	struct BP_PROC * p = (struct BP_PROC *) fp->fh;
    BP_FILE * fh = (BP_FILE *) p->fh;
    struct adios_index_var_struct_v1 * var_root;
    int i;

    // Perform variable ID mapping, since the input to this function is user-perceived
    int mapped_id = map_req_varid (fp, vi->varid);
    var_root = bp_find_var_byid (fh, mapped_id);

    ti->orig_blockinfo = inq_var_blockinfo(fp, vi, 1); // 1 -> use original, pretransform dimensions
    assert(ti->orig_blockinfo);

    // In streaming mode, we need to offset the transform metadata and length
    // arrays to start at the current timestep. For file mode, no such translation
    // is needed.
    int streaming_block_offset;
    if (p->streaming) {
    	int time = adios_step_to_time_v1(fp, var_root, 0);
    	streaming_block_offset = get_var_start_index(var_root, time);
    } else {
    	streaming_block_offset = 0;
    }

    assert(streaming_block_offset < var_root->characteristics_count);
    assert(streaming_block_offset + vi->sum_nblocks <= var_root->characteristics_count);

    // Allocate and fill the transform_metadatas array
    ti->transform_metadatas = (ADIOS_TRANSFORM_METADATA*)malloc(vi->sum_nblocks * sizeof(ADIOS_TRANSFORM_METADATA));
    assert(ti->transform_metadatas);
    for (i = 0; i < vi->sum_nblocks; i++) {
    	const struct adios_index_characteristic_transform_struct *transform_char = &var_root->characteristics[streaming_block_offset + i].transform;

    	ti->transform_metadatas[i] = (ADIOS_TRANSFORM_METADATA){
    		.length = transform_char->transform_metadata_len,
    		.content = transform_char->transform_metadata,
    	};
    }

    return 0;
}

// NCSU ALACRITY-ADIOS - Dummy function stubs for the staged and staged1 read transports; move to those files at some point
ADIOS_TRANSINFO * adios_read_bp_staged_inq_var_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi) {
    return NULL;
}
ADIOS_TRANSINFO * adios_read_bp_staged1_inq_var_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi) {
    return NULL;
}
int adios_read_bp_staged_inq_var_trans_blockinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    return 1;
}
int adios_read_bp_staged1_inq_var_trans_blockinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO *ti) {
    return 1;
}



uint64_t get_req_datasize (const ADIOS_FILE * fp, read_request * r, struct adios_index_var_struct_v1 * v)
{
    ADIOS_SELECTION * sel = r->sel;
    uint64_t datasize = bp_get_type_size (v->type, "");
    int i, pgidx, ndims;
    const struct BP_PROC * p = (struct BP_PROC *) fp->fh;

    if (sel->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        for (i = 0; i < sel->u.bb.ndim; i++)
        {
            datasize *=  sel->u.bb.count[i];
        }
    }
    else if (sel->type == ADIOS_SELECTION_POINTS)
    {
        datasize *= sel->u.points.npoints;
    }
    else if (sel->type == ADIOS_SELECTION_WRITEBLOCK)
    {
        //pgidx = adios_wbidx_to_pgidx (fp, r);
        // NCSU ALACRITY-ADIOS: Adding absolute PG indexing, but *only* in non-streaming
    	// mode (absolute writeblocks are interpreted as timestep-relative when in
    	// streaming mode)
        pgidx = sel->u.block.is_absolute_index && !p->streaming ?
                    sel->u.block.index :
                    adios_wbidx_to_pgidx (fp, r, 0);
        // NCSU ALACRITY-ADIOS: Adding sub-PG writeblock read support
        if (sel->u.block.is_sub_pg_selection) {
            datasize = sel->u.block.nelements;
        } else {
            // NCSU ALACRITY-ADIOS: This used to not be in this else block
            ndims = v->characteristics[pgidx].dims.count;
            for (i = 0; i < ndims; i++)
            {
                datasize *= v->characteristics[pgidx].dims.dims[i * 3];
            }
        }
    }

    return datasize;
}

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
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);

    read_request * r;
    ADIOS_SELECTION * nullsel = 0;
    struct adios_index_var_struct_v1 * v;
    int i, ndim, ns, file_is_fortran, mapped_varid;
    uint64_t * dims = 0;

    mapped_varid = p->varid_mapping[varid];
    v = bp_find_var_byid (fh, mapped_varid);
    file_is_fortran = is_fortran_file (fh);

    r = (read_request *) malloc (sizeof (read_request));
    assert (r);

    if (!sel)
    {
        bp_get_and_swap_dimensions (fp, v, file_is_fortran,
                                    &ndim, &dims,
                                    &ns,
                                    file_is_fortran != futils_is_called_from_fortran()
                                    );

        nullsel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
        assert (nullsel);

        nullsel->type = ADIOS_SELECTION_BOUNDINGBOX;
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

    /* copy selection since we don't want to operate on user memory.
     */
    r->sel = (!nullsel ? a2sel_copy (sel) : nullsel);
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
    r->datasize = get_req_datasize (fp, r, v);
    r->priv = 0;
    r->next = 0;

    list_insert_read_request_next (&p->local_read_request_list, r);

    return 0;
}

int adios_read_bp_perform_reads (const ADIOS_FILE *fp, int blocking)
{
    BP_PROC * p = GET_BP_PROC (fp);
    read_request * r;
    ADIOS_VARCHUNK * chunk;

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

    while (p->local_read_request_list)
    {
        chunk = read_var (fp, p->local_read_request_list);

        // remove head from list
        r = p->local_read_request_list;
        p->local_read_request_list = p->local_read_request_list->next;
        a2sel_free (r->sel);
        r->sel = NULL;
        free(r);

        common_read_free_chunk (chunk);
    }

    return 0;
}

/* This routine split a 'big' request into smaller ones which can be fit into
 * buffer_size bytes of memory.
 */
static read_request * split_req (const ADIOS_FILE * fp, const read_request * r, int buffer_size)
{
    BP_FILE * fh = GET_BP_FILE (fp);

    read_request * h = 0;
    ADIOS_SELECTION * sel = r->sel;
    struct adios_index_var_struct_v1 * v;
    int type_size, n_elements, ndim;
    int i, j, varid, remain, done;
    uint64_t pos[32], subbb[32], start[32], count[32];

    log_debug ("split_req()\n");
    varid = r->varid; //map_req_varid (fp, r->varid); // NCSU ALACRITY-ADIOS: Bugfix: r->varid has already been mapped
    v = bp_find_var_byid (fh, varid);
    type_size = bp_get_type_size (v->type, "");
    assert (type_size);

    n_elements = buffer_size / type_size;

    printf ("n_elements = %d\n", n_elements);
    //TODO: handle string
    if (sel->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        ndim = sel->u.bb.ndim;
        // convert chunk size to position within the bounding box
        for (i = ndim - 1; i > -1; i--)
        {
            pos[i] = n_elements % sel->u.bb.count[i];
            assert (sel->u.bb.count[i]);
            n_elements /= sel->u.bb.count[i];
        }

        log_debug ("pos = ");
        for (i = 0; i < ndim; i++)
        {
            log_debug_cont ("%" PRIu64 " ", pos[i]);
        }
        log_debug_cont ("\n");

        // calculate sub-bounding-box
        for (i = ndim - 1; i > -1; i--)
        {
            if (pos[i] != sel->u.bb.count[i] - 1)
            {
                j = 0;
                while (j <= i && pos[j] == 0)
                {
                    subbb[j] = 0;
                    j++;
                }

                if (j <= i)
                {
                    subbb[j] = pos[j];
                    j++;
                    while (j <= i)
                    {
                        subbb[j] = sel->u.bb.count[j];
                        j++;
                    }
                }
                break;
            }
            else
            {
                subbb[i] = sel->u.bb.count[i];
            }
        }

        log_debug ("subbb = ");
        for (i = 0; i < ndim; i++)
        {
            log_debug_cont ("%" PRIu64 " ", subbb[i]);
        }
        log_debug_cont ("\n");

        memcpy (start, sel->u.bb.start, ndim * 8);
        memcpy (count, sel->u.bb.count, ndim * 8);

        while (1)
        {
            read_request * newreq = (read_request *) malloc (sizeof (read_request));
            assert (newreq);

            newreq->sel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
            assert (newreq->sel);
            newreq->sel->type = ADIOS_SELECTION_BOUNDINGBOX;
            newreq->sel->u.bb.ndim = ndim;
            newreq->sel->u.bb.start = malloc (ndim * 8);
            newreq->sel->u.bb.count = malloc (ndim * 8);
            assert (newreq->sel->u.bb.start);
            assert (newreq->sel->u.bb.count);

            memcpy (newreq->sel->u.bb.start,
                    start,
                    ndim * 8
                   );

            // check whether the start + count will be out of bound
            for (i = 0; i < ndim; i++)
            {
                if (start[i] + subbb[i] > sel->u.bb.start[i] + sel->u.bb.count[i])
                {
                    count[i] = sel->u.bb.start[i] + sel->u.bb.count[i] - start[i];
                }
                else
                {
                    count[i] = subbb[i];
                }
            }

            memcpy (newreq->sel->u.bb.count,
                    count,
                    ndim * 8
                   );

            log_debug ("bb: (");
            for (i = 0; i < ndim; i++)
            {
                log_debug_cont ("%" PRIu64 "", newreq->sel->u.bb.start[i]);
                if (i != ndim - 1)
                {
                    log_debug_cont (",");
                }
            }
            log_debug_cont (") (");
            for (i = 0; i < ndim; i++)
            {
                log_debug_cont ("%" PRId64 "", newreq->sel->u.bb.start[i] + newreq->sel->u.bb.count[i] - 1);
                if (i != ndim - 1)
                {
                    log_debug_cont (",");
                }
            }
            log_debug_cont (")\n");

            done = 0;
            for (i = ndim - 1; i > -1; i--)
            {
                // This dimension is finished.
                if (start[i] + count[i] == sel->u.bb.start[i] + sel->u.bb.count[i])
                {
                    start[i] = sel->u.bb.start[i];
                }
                else
                {
                    start[i] += count[i];
                    break;
                }
            }

            if (i == -1)
            {
                done = 1;
            }

            newreq->varid = r->varid;
            newreq->from_steps = r->from_steps;
            newreq->nsteps = r->nsteps;
            newreq->data = r->data;
            newreq->datasize = type_size;
            for (i = 0; i < ndim; i++)
            {
                newreq->datasize *= count[i];
            }

            newreq->priv = r->priv;
            newreq->next = 0;

            list_insert_read_request_next (&h, newreq);

            // all dimensions are finished and we are done
            if (done)
            {
                break;
            }
        }
    }
    else if (sel->type == ADIOS_SELECTION_POINTS)
    {
        remain = sel->u.points.npoints;
        while (remain)
        {
            read_request * newreq = (read_request *) malloc (sizeof (read_request));
            assert (newreq);

            newreq->sel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
            assert (newreq->sel);
            newreq->sel->type = ADIOS_SELECTION_POINTS;
            newreq->sel->u.points.ndim = sel->u.points.ndim;
            newreq->sel->u.points.npoints = (remain > n_elements ? n_elements : remain);
            newreq->sel->u.points.points = malloc (newreq->sel->u.points.npoints * newreq->sel->u.points.ndim * 8);
            assert (newreq->sel->u.points.points);
            memcpy (newreq->sel->u.points.points,
                    sel->u.points.points + (sel->u.points.npoints - remain) * sel->u.points.ndim,
                    newreq->sel->u.points.npoints * sel->u.points.ndim * 8
                   );

            newreq->varid = r->varid;
            newreq->from_steps = r->from_steps;
            newreq->nsteps = r->nsteps;
            newreq->data = r->data;
            newreq->datasize = type_size * newreq->sel->u.points.npoints;
            newreq->priv = r->priv;
            newreq->next = 0;

            list_insert_read_request_next (&h, newreq);

            remain -= n_elements;
        }
    }
    else if (sel->type == ADIOS_SELECTION_WRITEBLOCK)
    {

    }

    return h;
}

int adios_read_bp_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk)
{
    BP_PROC * p = GET_BP_PROC (fp);

    read_request * r;
    ADIOS_VARCHUNK * varchunk;
/*
 *  RETURN:         0: all chunks have been returned previously,
 *                     no need to call again (chunk is NULL, too)
 *                  1: some chunks are/will be available, call again
 *                  <0 on error, sets adios_errno too
 */
    log_debug ("adios_read_bp_check_reads()\n");

    if (!p->local_read_request_list)
    {
        return 0;
    }

    // if memory is pre-allocated
    if (p->local_read_request_list->data)
    {
        log_debug ("adios_read_bp_check_reads(): memory is pre-allocated\n");
        varchunk = read_var (fp, p->local_read_request_list);

        if (varchunk)
        {
            // remove head from list
            r = p->local_read_request_list;
            p->local_read_request_list = p->local_read_request_list->next;
            a2sel_free (r->sel);
            r->sel = NULL;
            free(r);

            * chunk = varchunk;
            return 1;
        }
        else
        {
            return adios_errno;
        }
    }
    else // if memory is not pre-allocated
    {
        log_debug ("adios_read_bp_check_reads(): memory is not pre-allocated\n");
        // memory is large enough to contain the data
        if (chunk_buffer_size >= p->local_read_request_list->datasize)
        {
            log_debug ("adios_read_bp_check_reads(): memory is large enough to contain the data (%" PRIu64 ")\n",
                       p->local_read_request_list->datasize);
            assert (p->local_read_request_list->datasize);
            p->b = realloc (p->b, p->local_read_request_list->datasize);
            p->local_read_request_list->data = p->b;

            varchunk = read_var (fp, p->local_read_request_list);

            if (varchunk)
            {
                // remove head from list
                r = p->local_read_request_list;
                p->local_read_request_list = p->local_read_request_list->next;
                a2sel_free (r->sel);
                r->sel = NULL;
                free(r);

                * chunk = varchunk;
                return 1;
            }
            else
            {
                return adios_errno;
            }
        }
        else // memory is smaller than what it takes to read the entire thing in.
        {
            log_debug ("adios_read_bp_check_reads(): memory is not large enough to contain the data (%" PRIu64 ")\n",
                       p->local_read_request_list->datasize);
            read_request * subreqs = split_req (fp, p->local_read_request_list, chunk_buffer_size);
            assert (subreqs);

            // remove head from list
            r = p->local_read_request_list;
            p->local_read_request_list = p->local_read_request_list->next;
            a2sel_free (r->sel);
            r->sel = NULL;
            free(r);

            r = subreqs;
            while (r->next)
            {
                r = r->next;
            }

            r->next = p->local_read_request_list;
            p->local_read_request_list = subreqs;

            p->b = realloc (p->b, p->local_read_request_list->datasize);
            p->local_read_request_list->data = p->b;

            varchunk = read_var (fp, p->local_read_request_list);

            if (varchunk)
            {
                // remove head from list
                r = p->local_read_request_list;
                p->local_read_request_list = p->local_read_request_list->next;
                a2sel_free (r->sel);
                r->sel = NULL;
                free(r);

                * chunk = varchunk;
                return 1;
            }
            else
            {
                return adios_errno;
            }
        }
    }

    return 0;
}

int adios_read_bp_get_attr_byid (const ADIOS_FILE * fp, int attrid, enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    int i;
    BP_FILE * fh = GET_BP_FILE (fp);
    struct adios_index_attribute_struct_v1 * attr_root;
    struct adios_index_var_struct_v1 * var_root, * v1;
    int file_is_fortran, last_step = fp->last_step, show_hidden_attrs;
    uint64_t k, attr_c_index, var_c_index;

    adios_errno = 0;

    show_hidden_attrs = 0;
    for (i = 0; i < fp->nattrs; i++)
    {
        if (strstr (fp->attr_namelist[i], "__adios__"))
        {
            show_hidden_attrs = 1;
            break;
        }
    }

    attr_root = fh->attrs_root; /* need to traverse the attribute list of the group */
    i = 0;

    if (show_hidden_attrs)
    {
        while (i < attrid && attr_root)
        {
            i++;
            attr_root = attr_root->next;
        }
    }
    else
    {
        while (i < attrid && attr_root)
        {
            if (strstr (attr_root->attr_path, "__adios__"))
            {
            }
            else
            {
                i++;
            }

            attr_root = attr_root->next;
        }

        while (attr_root && strstr (attr_root->attr_path, "__adios__"))
        {
            attr_root = attr_root->next;
        }
    }

    assert (attr_root);

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
        *type = attr_root->type;
        int type_size = bp_get_type_size (attr_root->type, attr_root->characteristics[attr_c_index].value);
        if (*type == adios_string) {
            *size = type_size;
        } else {
            *size = attr_root->nelems * type_size;
        }
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
               && !strcmp(var_root->group_name, attr_root->group_name)
               )
                break;
            var_root = var_root->next;
        }

        if (!var_root)
        {
            var_root = fh->vars_root;
            while (var_root)
            {
                if (var_root->id == attr_root->characteristics[attr_c_index].var_id
                   && !strcmp(var_root->group_name, attr_root->group_name))
                    break;
                var_root = var_root->next;
            }
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
            int varid = 0;
            v1 = fh->vars_root;
            while (v1 && v1 != var_root)
            {
                v1 = v1->next;
                varid++;
            }

            start = 0;
            count = var_root->characteristics[var_c_index].dims.dims[0];
            snprintf(varname, 512, "%s/%s", var_root->var_path, var_root->var_name);
            tmpdata = (char *) malloc (count+1);
            assert (tmpdata);

            r = (read_request *) malloc (sizeof (read_request));
            assert (r);

            r->sel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
            r->sel->type = ADIOS_SELECTION_BOUNDINGBOX;
            r->sel->u.bb.ndim = 1;
            r->sel->u.bb.start = &start;
            r->sel->u.bb.count = &count;
            r->varid = varid;
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

int  adios_read_bp_get_dimension_order (const ADIOS_FILE *fp)
{
    BP_FILE * fh = GET_BP_FILE (fp);
    return is_fortran_file (fh);
}

void adios_read_bp_reset_dimension_order (const ADIOS_FILE *fp, int is_fortran)
{
    BP_FILE * fh = GET_BP_FILE (fp);
    struct bp_index_pg_struct_v1 ** root = &(fh->pgs_root);
    struct bp_minifooter * mh = &(fh->mfooter);
    uint64_t i;

    for (i = 0; i < mh->pgs_count; i++) {
        is_fortran ? ((*root)->adios_host_language_fortran = adios_flag_yes)
               : ((*root)->adios_host_language_fortran = adios_flag_no);
        root = &(*root)->next;
    }
}

void adios_read_bp_get_groupinfo (const ADIOS_FILE *fp, int *ngroups, char ***group_namelist, uint32_t **nvars_per_group, uint32_t **nattrs_per_group)
{
    BP_FILE * fh = GET_BP_FILE (fp);
    int i, j, offset;

    * ngroups = fh->gvar_h->group_count;

    *group_namelist = (char **) malloc (sizeof (char *) * fh->gvar_h->group_count);
    for (i = 0; i < fh->gvar_h->group_count; i++)
    {
        (*group_namelist)[i] = malloc (strlen (fh->gvar_h->namelist[i]) + 1);
        assert ((*group_namelist)[i]);

        memcpy ((*group_namelist)[i], fh->gvar_h->namelist[i], strlen (fh->gvar_h->namelist[i]) + 1);
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

/* Check if a variable is timed. This is solely done by checking whether
 * a variable is tagged with time in XML.
 */
int adios_read_bp_is_var_timed (const ADIOS_FILE *fp, int varid)
{
    BP_FILE * fh = GET_BP_FILE (fp);

    struct adios_index_var_struct_v1 * v;
    //struct adios_index_characteristic_struct_v1 ch;
    int retval = 0, ndim, k;
    uint64_t gdims[32];

    v = bp_find_var_byid (fh, varid);
    //ch = v->characteristics[0];
    //ndim = ch.dims.count; //ndim possibly has 'time' dimension
    // NCSU ALACRITY-ADIOS - An optimization. Not sure why it was originally added, but it works
    struct adios_index_characteristic_dims_struct_v1 *dims = &v->characteristics[0].dims;
    ndim = dims->count; //ndim possibly has 'time' dimension

    log_debug ("adios_read_bp_is_var_timed: varid = %d, ndim = %d\n", varid, ndim);

    if (ndim == 0)
    {
        return 0;
    }

    for (k = 0; k < ndim; k++)
    {
        // NCSU ALACRITY-ADIOS - An optimization
        gdims[k] = dims->dims[k * 3 + 1]; //ch.dims.dims[k * 3 + 1];
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

static int map_req_varid (const ADIOS_FILE * fp, int varid)
{
    BP_PROC * p = GET_BP_PROC (fp);

    return p->varid_mapping[varid];
}

/* This routine converts the write block index, which is of a particular step,
 * to the adios internal PG index.
 */
static int adios_wbidx_to_pgidx (const ADIOS_FILE * fp, read_request * r, int step_offset)
{
    BP_FILE * fh = GET_BP_FILE (fp);

    int time, start_idx, stop_idx, c, idx;
    int mapped_varid, ridx;
    struct adios_index_var_struct_v1* v;

    if (r->sel->type != ADIOS_SELECTION_WRITEBLOCK)
    {
        return -1;
    }

    time = adios_step_to_time (fp, r->varid, r->from_steps + step_offset);
    mapped_varid = r->varid; //map_req_varid (fp, r->varid); // NCSU ALACRITY-ADIOS: Bugfix: r->varid has already been mapped
    v = bp_find_var_byid (fh, mapped_varid);

    start_idx = get_var_start_index (v, time);
    stop_idx = get_var_stop_index (v, time);
    if (start_idx < 0 || stop_idx < 0)
    {
        adios_error (err_no_data_at_timestep,
                     "No data at step %d\n",
                     r->from_steps);
    }

    ridx =  r->sel->u.block.index;
    c = -1;
    idx = start_idx;
    while (idx <= stop_idx)
    {
        if (v->characteristics[idx].time_index == time)
        {
            c++;
        }

        if (c < ridx)
        {
            idx++;
        }
        else
        {
            break;
        }
    }

    if (c != ridx)
    {
        log_debug ("Error in adios_wbidx_to_pgidx().\n");
    }

    return idx;
}

/* This routine reads a write block. The 'index' value in the selection is
 * the block index within the context of the current step. Therefore, we
 * need to translate it to an absolute index.
 */
static ADIOS_VARCHUNK * read_var_wb (const ADIOS_FILE * fp, read_request * r)
{
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);

    struct adios_index_var_struct_v1 * v;
    int i, j, varid, start_idx, idx;
    int ndim, has_subfile;
    uint64_t ldims[32], gdims[32], offsets[32];
    int size_of_type;
    uint64_t slice_offset, slice_size;
    void * data;
    ADIOS_VARCHUNK * chunk;
    MPI_Status status;
    const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb;// NCSU ALACRITY-ADIOS

    adios_errno = 0;

    has_subfile = has_subfiles (fh);
    data = r->data;
    varid = r->varid; //varid = map_req_varid (fp, r->varid); // NCSU ALACRITY-ADIOS: Bugfix: r->varid has already been mapped
    v = bp_find_var_byid (fh, varid);

    // NCSU ALACRITY-ADIOS: Add support for absolute PG index for efficiency
    //time = adios_step_to_time (fp, r->varid, r->from_steps);
    //idx = adios_wbidx_to_pgidx (fp, r);
    assert(r->sel->type == ADIOS_SELECTION_WRITEBLOCK);
    wb = &r->sel->u.block;

    for (i = 0; i < r->nsteps; i++)
    {
        // NCSU ALACRITY-ADIOS: Adding absolute PG indexing, but *only* in non-streaming
    	// mode (absolute writeblocks are interpreted as timestep-relative when in
    	// streaming mode)
        idx = wb->is_absolute_index && !p->streaming ?
                  wb->index :
                  adios_wbidx_to_pgidx (fp, r, i);
        //if (!wb->is_absolute_index) printf("Timestep-relative writeblock index used!\n");
        assert (idx >= 0);

        ndim = v->characteristics [idx].dims.count;
        size_of_type = bp_get_type_size (v->type, v->characteristics [idx].value);

        if (ndim == 0)
        {
            r->datasize = size_of_type;
            slice_size = size_of_type;
            start_idx = 0; // OPS macros below need it

            if (v->type == adios_string)
            {
                size_of_type--;
            }

            slice_offset = v->characteristics[idx].payload_offset;

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS1
            }
            else
            {
                MPI_FILE_READ_OPS2
            }

            memcpy((char *)data, fh->b->buff + fh->b->offset, size_of_type);

            if (fh->mfooter.change_endianness == adios_flag_yes)
            {
                 change_endianness ((char *)data,
                                    size_of_type,
                                    v->type
                                   );
            }

            if (v->type == adios_string)
            {
                ((char*)data)[size_of_type] = '\0';
            }

            data = (char *) data + size_of_type;
        }
        else
        {
            // NCSU ALACRITY-ADIOS: Added sub-PG writeblock selection support
            // If this is a sub-PG selection, use nelements to compute slice_size
            // instead
            if (wb->is_sub_pg_selection) {
                // The start and end of the sub-PG selection must fall within the PG
                slice_size = wb->nelements * size_of_type;
            } else {
                // NCSU ALACRITY-ADIOS: This used to not be inside an else block
                // Else, do the old method of computing PG size from bounds
                slice_size = size_of_type;

                /* To get ldims for the chunk and then calculate payload size */
                bp_get_dimension_characteristics(&(v->characteristics[idx]),
                                             ldims, gdims, offsets);

                for (j = 0; j < ndim; j++)
                {
                    slice_size *= ldims [j];
                }
            }

            r->datasize = slice_size;
            /* Note: MPI_FILE_READ_OPS1 - for reading single BP file.
             *       MPI_FILE_READ_OPS2 - for reading those with subfiles.
             * Whenever to use OPS macro, start_idx and idx variable needs to be
             * properly set.
             */
            start_idx = 0;
            slice_offset = v->characteristics[idx].payload_offset;

            // NCSU ALACRITY-ADIOS: Added sub-PG writeblock selection support
            // If this is a sub-PG read, add the element_offset within the PG to the base offset in the file
            if (wb->is_sub_pg_selection) {
                slice_offset += wb->element_offset * size_of_type;
            }

            if (!has_subfile)
            {
                MPI_FILE_READ_OPS1_BUF(data) // NCSU ALACRITY-ADIOS: Read data directly to user buffer
            }
            else
            {
                MPI_FILE_READ_OPS2_BUF(data) // NCSU ALACRITY-ADIOS: Read data directly to user buffer
            }

            // NCSU ALACRITY-ADIOS: Reading directly to user buffer eliminates the need for this memcpy (profiling revealed it was hurting performance for transformed data)
            //memcpy ((char *)data, fh->b->buff + fh->b->offset, slice_size);
            if (fh->mfooter.change_endianness == adios_flag_yes)
            {
                change_endianness ((char *)data, slice_size, v->type);
            }
        }
    }

    chunk = (ADIOS_VARCHUNK *) malloc (sizeof (ADIOS_VARCHUNK));
    assert (chunk);

    chunk->varid = r->varid;
    chunk->type = v->type;
    // NCSU ALACRITY-ADIOS - Added timestep information into varchunks
    chunk->from_steps = r->from_steps;
    chunk->nsteps = r->nsteps;
    chunk->sel = a2sel_copy (r->sel);
    chunk->data = data;

    return chunk;
}

#if 0
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
