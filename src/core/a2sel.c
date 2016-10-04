/*
 * a2sel.c
 *
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>

#include "public/adios_error.h"
#include "core/a2sel.h"

extern int adios_errno;

ADIOS_SELECTION * a2sel_boundingbox (int ndim, const uint64_t *start, const uint64_t *count)
{
    adios_errno = err_no_error;
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_BOUNDINGBOX;
        sel->u.bb.ndim = ndim;
        sel->u.bb.start = (uint64_t *) malloc (ndim * sizeof(uint64_t));
        sel->u.bb.count = (uint64_t *) malloc (ndim * sizeof(uint64_t));
        memcpy (sel->u.bb.start, start, ndim * sizeof(uint64_t));
        memcpy (sel->u.bb.count, count, ndim * sizeof(uint64_t));
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for bounding box selection\n");
    }
    return sel;
}

ADIOS_SELECTION * a2sel_points (int ndim, uint64_t npoints, const uint64_t *points,
                                ADIOS_SELECTION * container, int free_points_on_delete)
{
    adios_errno = err_no_error;
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_POINTS;
        sel->u.points.ndim = ndim;
        sel->u.points.npoints = npoints;
        sel->u.points.points = (uint64_t *) points;
        sel->u.points.container_selection = container;
        sel->u.points._free_points_on_delete = free_points_on_delete;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for points selection\n");
    }
    return sel;
}

ADIOS_SELECTION * a2sel_writeblock (int index)
{
    adios_errno = err_no_error;
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_WRITEBLOCK;
        sel->u.block.index = index;
        // NCSU ALACRITY-ADIOS: Set the writeblock selection to be a full-PG selection by default
        sel->u.block.is_absolute_index = 0;
        sel->u.block.is_sub_pg_selection = 0;
        sel->u.block.element_offset = 0;
        sel->u.block.nelements = 0;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for writeblock selection\n");
    }
    return sel;
}

ADIOS_SELECTION * a2sel_auto (char *hints)
{
    adios_errno = err_no_error;
    ADIOS_SELECTION * sel = (ADIOS_SELECTION *) malloc (sizeof(ADIOS_SELECTION));
    if (sel) {
        sel->type = ADIOS_SELECTION_AUTO;
        sel->u.autosel.hints = hints;
    } else {
        adios_error(err_no_memory, "Cannot allocate memory for auto selection\n");
    }
    return sel;
}

void a2sel_free (ADIOS_SELECTION *sel)
{
    if (!sel)
        return;
    if (sel->type == ADIOS_SELECTION_POINTS)
    {
        if (sel->u.points.container_selection != NULL)
        {
            a2sel_free (sel->u.points.container_selection);
            sel->u.points.container_selection = NULL;
        }
        if (sel->u.points._free_points_on_delete)
        {
            free (sel->u.points.points);
        }
    }
    else if (sel->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        if (sel->u.bb.start) {
            free (sel->u.bb.start);
            sel->u.bb.start = NULL;
        }
        if (sel->u.bb.count) {
            free (sel->u.bb.count);
            sel->u.bb.count = NULL;
        }
    }
    free(sel);
}

ADIOS_SELECTION * a2sel_copy (const ADIOS_SELECTION * sel)
{
    ADIOS_SELECTION * nsel;

    nsel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
    assert (nsel);

    nsel->type = sel->type;

    if (sel->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        nsel->u.bb.ndim = sel->u.bb.ndim;
        nsel->u.bb.start = (uint64_t *) malloc (sel->u.bb.ndim * 8);
        nsel->u.bb.count = (uint64_t *) malloc (sel->u.bb.ndim * 8);
        assert (nsel->u.bb.start && nsel->u.bb.count);

        memcpy (nsel->u.bb.start, sel->u.bb.start, sel->u.bb.ndim * 8);
        memcpy (nsel->u.bb.count, sel->u.bb.count, sel->u.bb.ndim * 8);
    }
    else if (sel->type == ADIOS_SELECTION_POINTS)
    {
        nsel->u.points.ndim = sel->u.points.ndim;
        nsel->u.points.npoints = sel->u.points.npoints;
        if (sel->u.points.container_selection) {
            nsel->u.points.container_selection = a2sel_copy (sel->u.points.container_selection);
        } else {
            nsel->u.points.container_selection = NULL;
        }
        nsel->u.points.points = (uint64_t *) malloc (nsel->u.points.npoints * nsel->u.points.ndim * 8);
	nsel->u.points._free_points_on_delete = 1; // junmin
        assert (nsel->u.points.points);

        memcpy (nsel->u.points.points, sel->u.points.points, sel->u.points.npoints * sel->u.points.ndim * 8);
    }
    else if (sel->type == ADIOS_SELECTION_WRITEBLOCK)
    {
        nsel->u.block.index = sel->u.block.index;
        // NCSU ALACRITY-ADIOS: Copy the new fields
        nsel->u.block.is_absolute_index = sel->u.block.is_absolute_index;
        nsel->u.block.is_sub_pg_selection = sel->u.block.is_sub_pg_selection;
        nsel->u.block.element_offset = sel->u.block.element_offset;
        nsel->u.block.nelements = sel->u.block.nelements;
    }
    else if (sel->type == ADIOS_SELECTION_AUTO)
    {
        //TODO
    }
    else
    {
        //adios_error (err_invalid_argument, "Wrong ADIOS selection type.\n");
    }

    return nsel;
}


void a2sel_points_1DtoND_box (uint64_t npoints, uint64_t *pts1d,
                                             int ndim, uint64_t *start, uint64_t *count, int global,
                                             uint64_t *ptsNd)
{
    int n, d;
    assert (ndim > 0);

    uint64_t product[ndim];
    product[ndim-1] = count[ndim-1];
    for (d = ndim-2; d >= 0; d--) {
        product[d] = product[d+1] * count[d];
    }
    // Note, product[0] is never used

    // if global conversion, add start[] to each coordinate
    uint64_t extraoffs[ndim];
    for (d = 0; d < ndim; d++)
    {
        extraoffs[d] = (global ? start[d] : 0);
    }

    uint64_t *pN = ptsNd;
    uint64_t *p1 = pts1d;
    uint64_t rem;
    for (n = 0; n < npoints; n++)
    {
        rem = *p1;
        for (d = 0; d < ndim-1; d++)
        {
            *pN = rem / product[d+1] + extraoffs[d];
            rem = rem % product[d+1];
            pN++;
        }
        *pN = rem + extraoffs[ndim-1]; // last dimension is just the remainder
        pN++;
        p1++;
    }
 }

ADIOS_SELECTION * a2sel_points_1DtoND (ADIOS_SELECTION * pointsinbox1D, int global)
{
    if (!pointsinbox1D)
    {
        adios_error (err_invalid_selection, "in adios_selection_points_1DtoND(): NULL selection provided\n");
        return NULL;
    }

    if (pointsinbox1D->type != ADIOS_SELECTION_POINTS ||
        !pointsinbox1D->u.points.container_selection)
    {
        adios_error (err_invalid_selection, "in adios_selection_points_1DtoND(): "
                "Only point selections with a container selection can be converted\n");
        return NULL;
    }

    if (pointsinbox1D->u.points.container_selection->type != ADIOS_SELECTION_BOUNDINGBOX)
    {
        adios_error (err_invalid_selection, "in adios_selection_points_1DtoND(): "
                "Point selection's container can only be a bounding box\n");
                return NULL;
    }

    if (pointsinbox1D->u.points.ndim != 1)
    {
        adios_error (err_invalid_selection, "in adios_selection_points_1DtoND(): "
                "Only 1D points can be converted\n");
                return NULL;
    }

    uint64_t *ptsNd = (uint64_t *) malloc (pointsinbox1D->u.points.container_selection->u.bb.ndim * pointsinbox1D->u.points.npoints * sizeof(uint64_t));
    if (!ptsNd)
    {
        adios_error (err_no_memory, "in adios_selection_points_1DtoND(): "
                "Not enough memory to allocate %d-dimensional point selection for %" PRIu64 "points\n",
                pointsinbox1D->u.points.container_selection->u.bb.ndim, pointsinbox1D->u.points.npoints);
                return NULL;
    }

    ADIOS_SELECTION * container = a2sel_copy (pointsinbox1D->u.points.container_selection);

    a2sel_points_1DtoND_box (pointsinbox1D->u.points.npoints, pointsinbox1D->u.points.points,
                                            container->u.bb.ndim, container->u.bb.start, container->u.bb.count, global,
                                            ptsNd);
    int ndim = container->u.bb.ndim;
    if (global)
    {
        a2sel_free(container);
        container = NULL;
    }
    ADIOS_SELECTION * result = a2sel_points(ndim, pointsinbox1D->u.points.npoints, ptsNd,
                                            container, 1);
    return result;
}
