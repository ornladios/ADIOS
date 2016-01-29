#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>

#include "common_query.h"
#include "adios_query_hooks.h"
#include "public/adios_error.h"
#include "core/common_read.h"
#include "core/adios_logger.h"
#include "query_utils.h"

// get the number of elements in a block (not size in bytes)
// size in bytes = get_nelements(v,i) * adios_type_size (v->type, NULL)
static uint64_t get_nelements (int ndim, uint64_t *dims)
{
    int i;
    //int size = adios_type_size (v->type, NULL);
    uint64_t size = 1;
    for (i = 0; i < ndim; i++) {
        size *= dims[i];
    }
    return size;
}

// get the size of a block in bytes
/*
static uint64_t get_blocksize_bytes (int ndim, uint64_t *dims, enum ADIOS_DATATYPES type)
{
    return adios_type_size (type, NULL) * get_nelements (ndim, dims);
}
*/

static ADIOS_VARINFO * adios_query_find_varinfo (ADIOS_FILE *f, ADIOS_QUERY *q, const char *varname)
{
    if (!q->left && !q->right) {
        if (!strcmp(q->varName, varname) && f == q->file) {
            return q->varinfo;
        }
    }
    if (q->left) {
       ADIOS_VARINFO *v = adios_query_find_varinfo (f, q->left, varname);
       if (v)
           return v;
    }
    if (q->right) {
       ADIOS_VARINFO *v = adios_query_find_varinfo (f, q->right, varname);
       if (v)
           return v;
    }
    return NULL;
}

/* Calculate the 1D contiguous index position from the N dimensional coordinates.
 * lcoords are assumed at this point, that they fall inside the
 * {0..dims[0]} x ... x {0...dims[ndim-1]}  N-dimensional cube
 */
static uint64_t adios_query_calc_position (int ndim, uint64_t * dims, uint64_t * lcoords)
{
    int n;
    uint64_t pos = lcoords[ndim-1];
    uint64_t slice_size = dims[ndim-1];
    for (n = ndim-2; n >= 0; n--) {
        pos += lcoords[n] * slice_size;
        slice_size *= dims[ndim-1];
    }
    return pos;
}

/* Copy data from pointvalues[i] to data[X], i=0..npoints-1, where
 * X is calculated from the local position of the point in the bounding box.
 * Copy only if the point falls inside the box.
 */
static void adios_query_copy_points_to_bb (
                ADIOS_SELECTION_POINTS_STRUCT * pointsel,
                char * pointvalues,
                int elemsize,
                ADIOS_SELECTION * bb,
                char * data
            )
{
    assert (bb->type == ADIOS_SELECTION_BOUNDINGBOX);
    assert (pointsel->ndim == bb->u.bb.ndim);
    assert (pointsel->ndim <= 32);
    uint64_t npoints = pointsel->npoints;
    int ndim = pointsel->ndim;
    uint64_t n, i, coord_idx, coord;
    int falls_outside;
    uint64_t lcoords[32];

    for (n=0; n < npoints; n++) {
        // point coordinates = pointsel->u.points.points[n*ndim..(n+1)*ndim-1] into data[]
        falls_outside = 0;
        coord_idx = n*ndim;
        for (i=0; i < ndim; i++) {
            coord = pointsel->points[coord_idx];
            // check if point is in the box in the first place
            if (coord < bb->u.bb.start[n] || bb->u.bb.start[n]+bb->u.bb.count[n] <= coord) {
                falls_outside = 1;
                break;
            }
            // calculate the local point coordinate in bb
            lcoords[i] = coord - bb->u.bb.start[i];
            coord_idx++;
        }
        if (!falls_outside) {
            // copy elemsize bytes from &pointvalues[n*elemsize] into data[X]
            // where X is the local coordinate of point n in the boundingbox bb
            coord_idx = adios_query_calc_position(ndim, bb->u.bb.count, lcoords);
            memcpy (data+coord_idx*elemsize, pointvalues+n*elemsize, elemsize);
        }
    }
}

/* Copy data from block_data[] to bb_data[].
 * 'block_data' has 'nelements' of data points of 'elemsize' bytes.
 * Position of data in block is calculated from varinfo->blockinfo[block index] and block_start_offset.
 * 'bb_data' is user allocated, assumed to cover the size of bounding box 'bb'.
 * Copy only that portion that falls inside the bb box.
 * block_start_offset is the starting position in the original writeblock from which 'block' contains 'nelements'
 * data elements. It is used in calculating position in 'bb' but not in accessing data in 'block'.
 */
static void adios_query_copy_block_to_bb (ADIOS_VARINFO * vi, int block_index, char * block_data,
        uint64_t nelements, int elemsize,
        uint64_t block_start_offset,
        ADIOS_SELECTION * bb, char * bb_data)
{
    assert (bb->type == ADIOS_SELECTION_BOUNDINGBOX);
    assert (bb->u.bb.ndim > 0);
    assert (bb->u.bb.ndim <= 32);
    int ndim = bb->u.bb.ndim;

    // FIXME: block_start_offset is not supported yet. The calculations below do not use it. New design of the copy
    // algorithm is needed
    if (block_start_offset > 0) {
        log_error ("ADIOS QUERY ERROR: The %s function does not support partial writeblock "
                "selecions yet. Demand implementation from the ADIOS developers.\n", __func__);
        return;
    }

    // dimension and offset of current writeblock in global array
    uint64_t *wboffs = vi->blockinfo[block_index].start;
    uint64_t *wbdims = vi->blockinfo[block_index].count;

#if 0
    if (ndim == 1)
    {
        // 1D array, simple contiguous copy
        // starting coordinate of writeblock in the bounding box is wboffs[]-offs[]
        // will copy nelements data points from block[block_start..] -> data[bb_start]
        uint64_t block_start = 0;
        int64_t nelems = nelements;
        int64_t bb_start = wboffs[0] + block_start_offset - bb->u.bb.start[0];
        if (bb_start < 0) {
            // outside of left side of bounding box
            nelems += bb_start;    // decrease nelements
            block_start -= bb_start;  // positive integer
            bb_start = 0;
        }
        uint64_t bb_end = bb_start + nelems;
        if (bb_start + nelems > bb->u.bb.count[0]) {
            // reaching beyond the right side of the bounding box
            nelems = bb->u.bb.count[0] - bb_start;
        }
        if (nelems > 0) {
            memcpy (data+bb_start*elemsize, block+block_start*elemsize, nelems*elemsize);
        }
    }
    else if (ndim == 2) // FIXME remove this incomplete branch
    {
        // calculate 'relative to original bounding box' coordinate from 'starting point in writeblock'
        // starting coordinate of writeblock in the bounding box is wboffs[]-offs[]
        int64_t i, j, j_min, j_max, pos, n;
        // FIXME: incorrect starting point calculation with block_start_offset
        i = wboffs[0] - bb->u.bb.start[0] + block_start_offset / wbdims[1];
        j_min = wboffs[1] - bb->u.bb.start[1] + block_start_offset % wbdims[1];
        j = j_min;

        // calculate 1D position in output array
        pos = i*bb->u.bb.count[1]+j_min;

        // copy the elements from contiguous 'block' into 'data' in a non-contiguous manner
        j_max = wboffs[1] - bb->u.bb.start[1] + wbdims[1]; // j runs [j_min,j_max)
        n = 0;
        while (n < nelements) {
            //if (j==j_min) {
            //    print ("rank %d: Copy %dth element to v1[%d,%d] (=v1[%lld])\n", rank, n, i, j, pos);
            //}
            memcpy (data+pos*elemsize, block+n*elemsize, (j_max-j)*elemsize);
            n += j_max - j;
            j = j_min;
            i++;
            pos = i*bb->u.bb.count[1]+j_min;
        }
    }
    else {}
#endif

    /* check if there is any intersection */
    int flag;
    int i;
    for (i = 0; i < ndim; i++)
    {
        flag =   (wboffs[i] >= bb->u.bb.start[i]
                  && wboffs[i] < bb->u.bb.start[i] + bb->u.bb.count[i])
              || (wboffs[i] <  bb->u.bb.start[i]
                  && wboffs[i] + wbdims[i] > bb->u.bb.start[i] + bb->u.bb.count[i])
              || (wboffs[i] + wbdims[i] > bb->u.bb.start[i]
                  && wboffs[i] + wbdims[i] <= bb->u.bb.start[i] + bb->u.bb.count[i]);

        if (!flag) {
            return;
        }
    }

    /* determine how many (fastest changing) dimensions of the block can we copy in one swoop */
    uint64_t n_cont_elems = 1; // number of elements that can be contiguously copied
    for (i = ndim - 1; i > -1; i--)
    {
        if (wboffs[i] == bb->u.bb.start[i] && wbdims[i] == bb->u.bb.count[i])
        {
            n_cont_elems *= wbdims[i];
        }
        else
            break;
    }
    int hole_break = i;
    log_debug ("%s: hole_break = %d\n", __func__, hole_break);

    /* Handle different cases */
    if (hole_break == -1)
    {
        /* The complete bb happens to be exactly the writeblock, and the entire writeblock.
         * Just copy the block into the data (shifted by 'block_start_offset' elements)
         * This is a rare case. FIXME: can we eliminate this?
         */
        assert (n_cont_elems == nelements + block_start_offset);
        memcpy (bb_data+block_start_offset*elemsize, block_data, nelements*elemsize);
    }
    else if (hole_break == 0)
    {
        /* Block should not be copied entirely in the slowest changing dimension but
         * we still need to copy only one contiguous block.
         * Let's call the n-1 dimensional sub-block "row" here.
         */
        uint64_t block_nrows = 0; // number of "rows" in slowest dimension in block thats inside bb
        uint64_t block_startrow = 0; // start copying from this "row" in block
        uint64_t bb_startrow = 0; // starting "row" in bb to copy to

        uint64_t x = bb->u.bb.start[0] + bb->u.bb.count[0];
        if (wboffs[0] >= bb->u.bb.start[0])
        {
            // head of block is at or after the head of bb in dimension,
            //  block may fit into bb fully, partially or not at all
            block_startrow = 0;
            bb_startrow = wboffs[0] - bb->u.bb.start[0];
            if (wboffs[0] < x) {
                // head of the writeblock is inside the bb
                if (wboffs[0] + wbdims[0] > x)
                {
                    // the tail is outside of the target
                    block_nrows = x - wboffs[0];
                }
                else
                {
                    // the whole block fits inside the bb
                    block_nrows = wbdims[0];
                }
            }
            else
            {
                // else the whole writeblock is outside of (after) bb so nothing to do
                block_nrows = 0;
            }
        }
        else
        {
            // head of block is outside (before) bb, may cover bb, partially cover or not at all
            block_startrow = bb->u.bb.start[0] - wboffs[0];
            bb_startrow = 0;
            if (wboffs[0] + wbdims[0] > bb->u.bb.start[0])
            {
                // tail of block is inside bb, so there is coverage
                if (wboffs[0] + wbdims[0] < x)
                {
                    // tail of block is inside bb, partial cover
                    block_nrows = wboffs[0] + wbdims[0] - bb->u.bb.start[0];
                }
                else
                {
                    // block covers bb entirely
                    block_nrows = bb->u.bb.count[0];
                }
            }
            else
            {
                // else the whole writeblock is outside of (before) bb so nothing to do
                block_nrows = 0;
            }
        }

        log_debug ("%s: number of rows to copy = %" PRIu64 "\n", __func__, block_nrows);
        if (block_nrows)
        {
            uint64_t slice_size = block_nrows * n_cont_elems * elemsize;
            //uint64_t read_offset = block_startrow * n_cont_elems * elemsize;
            uint64_t write_offset = bb_startrow * n_cont_elems * elemsize;
            // write_offset += block_start_offset * elemsize; THIS IS WRONG
            log_debug ("%s: Copy %" PRIu64 " bytes from block to bb at offset = %" PRIu64 "\n", __func__, slice_size, write_offset);
            if (slice_size > 0) {
                memcpy (bb_data + write_offset, block_data, slice_size);
            }
        }
    }
    else
    {
        /* Block should not be copied entirely in the more than one dimension (hole+1 dimensions) so
         * we need to determine the largest contiguous block of the fastest N-hole-1 dimensions, and then
         * the number of such blocks to copy in a loop.
         */
        uint64_t block_nrows[32]; // number of "rows" in each non-contiguous dimension in block thats inside bb
        uint64_t block_startrow[32]; // start copying from this "row" in block
        uint64_t bb_startrow[32]; // starting "row" in bb to copy to

        memset(block_nrows, 0 , 32 * 8);
        memset(block_startrow, 0 , 32 * 8);
        memset(bb_startrow, 0 , 32 * 8);

        uint64_t x;

        for (i = 0; i < ndim; i++)
        {
            // for all dimension > hole, wboffs=bb.start and wbdims=bb.count, so
            // the same code will fill block_nrow=wbdims and block_startrow=0
            x = bb->u.bb.start[i] + bb->u.bb.count[i];
            if (wboffs[i] >= bb->u.bb.start[i])
            {
                // head of block is at or after the head of bb in dimension,
                //  block may fit into bb fully, partially or not at all
                block_startrow[i] = 0;
                bb_startrow[i] = wboffs[i] - bb->u.bb.start[i];
                if (wboffs[i] < x) {
                    // head of the writeblock is inside the bb
                    if (wboffs[i] + wbdims[i] > x)
                    {
                        // the tail is outside of the target
                        block_nrows[i] = x - wboffs[i];
                    }
                    else
                    {
                        // the whole block fits inside the bb
                        block_nrows[i] = wbdims[i];
                    }
                }
                else
                {
                    // else the whole writeblock is outside of (after) bb so nothing to do
                    block_nrows[i] = 0;
                }
            }
            else
            {
                // head of block is outside (before) bb, may cover bb, partially cover or not at all
                block_startrow[i] = bb->u.bb.start[i] - wboffs[i];
                bb_startrow[i] = 0;
                if (wboffs[i] + wbdims[i] > bb->u.bb.start[i])
                {
                    // tail of block is inside bb, so there is coverage
                    if (wboffs[0] + wbdims[i] < x)
                    {
                        // tail of block is inside bb, partial cover
                        block_nrows[i] = wboffs[i] + wbdims[i] - bb->u.bb.start[i];
                    }
                    else
                    {
                        // block covers bb entirely
                        block_nrows[i] = bb->u.bb.count[i];
                    }
                }
                else
                {
                    // else the whole writeblock is outside of (before) bb so nothing to do
                    block_nrows[i] = 0;
                }
            }
        }

        n_cont_elems = 1; // the contiguous piece in the fastest dimension(s) that we copy at once
        uint64_t block_stride = 1;
        uint64_t bb_stride = 1;

        for (i = ndim - 1; i >= hole_break; i--)
        {
            n_cont_elems *= block_nrows[i];
            bb_stride *= bb->u.bb.count[i];
            block_stride *= wbdims[i];
        }

        log_debug ("%s: Block calculation:\n", __func__);
        for (i = 0; i < ndim; i++)
        {
            log_debug ("   block_nrows[%d]=%" PRIu64 "\tblock_startrow[%d]=%" PRIu64 "\tbb_startrow=[%d]=%" PRIu64 "\n",
                    i, block_nrows[i], i, block_startrow[i], i, bb_startrow[i]);
            block_startrow[i] = 0;
        }

        uint64_t block_offset = 0; // not the same thing as block_start_offset!
        uint64_t bb_offset = 0;

        for (i = 0; i < ndim; i++)
        {
            block_offset = block_startrow[i] + block_offset * wbdims[i];
            bb_offset = bb_startrow[i] + bb_offset * bb->u.bb.count[i];
        }

        log_debug ("%s: Copy from block to bb, cont_elems = %" PRIu64 ", block_stride = %" PRIu64 ", bb_stride = %" PRIu64 "\n",
                __func__, n_cont_elems, block_stride, bb_stride);

        adios_util_copy_data (bb_data
                  ,block_data
                  ,0
                  ,hole_break
                  ,block_nrows
                  ,wbdims
                  ,bb->u.bb.count
                  ,bb_stride
                  ,block_stride
                  ,bb_offset
                  ,block_offset
                  ,n_cont_elems
                  ,elemsize
                  ,adios_flag_no
                  ,vi->type
                  );
    }
}


int common_query_read_boundingbox (
        ADIOS_FILE *f,
        ADIOS_QUERY *q,
        const char *varname,
        int timestep,
        unsigned int nselections,
        ADIOS_SELECTION *selections,
        ADIOS_SELECTION *bb,
        void *data
   )
{
    assert (q);
    assert (f);
    assert (varname);
    assert (bb->type == ADIOS_SELECTION_BOUNDINGBOX);
    assert (data);
    if (nselections == 0 || selections == NULL)
        return 0;
    assert (selections[0].type == ADIOS_SELECTION_POINTS ||
            selections[0].type == ADIOS_SELECTION_WRITEBLOCK);

    // Get var's varinfo from query if it has already been created during query evaluation
    int free_varinfo = 0;
    ADIOS_VARINFO *vinfo = adios_query_find_varinfo (f, q, varname);
    if (!vinfo) {
        vinfo = adios_inq_var (f, varname);
        free_varinfo = 1;
    }
    if (!vinfo) {
        adios_error (err_corrupted_variable,
                "Corrupted variable in file. Could not get information about variable %s.\n",
                varname);
        return adios_errno;
    }
    if (!vinfo->blockinfo) {
        adios_inq_var_blockinfo(f, vinfo);
    }

    int elemsize = adios_type_size(vinfo->type, NULL);

    int n;
    for (n = 0; n < nselections; n++)
    {
        uint64_t nelements, element_offset;
        if (selections[n].type == ADIOS_SELECTION_WRITEBLOCK) {
            if (selections[n].u.block.is_sub_pg_selection) {
                nelements = selections[n].u.block.nelements;
                element_offset = selections[n].u.block.element_offset;
            } else {
                nelements = get_nelements(vinfo->ndim, vinfo->blockinfo[selections[n].u.block.index].count);
                element_offset = 0;
            }
        } else {  // ADIOS_SELECTION_POINTS
            nelements = selections[n].u.points.npoints;
            element_offset = 0; // not used in this case at all
        }

        char *d = (char *) malloc (nelements * elemsize);
        adios_schedule_read (f, &selections[n], varname, timestep, 1, d);
        adios_perform_reads (f, 1);

        // place data into user's allocated data described by the bounding box
        if (selections[n].type == ADIOS_SELECTION_WRITEBLOCK)
        {
            adios_query_copy_block_to_bb(vinfo, selections[n].u.block.index, d,
                    nelements, elemsize, element_offset, bb, data);
        }
        else // ADIOS_SELECTION_POINTS
        {
            adios_query_copy_points_to_bb(&(selections[n].u.points), d, elemsize, bb, data);
        }

        free (d);
    }

    if (free_varinfo)
        adios_free_varinfo(vinfo);
    return 0;
}

