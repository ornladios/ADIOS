/*
 * a2sel.h  ADIOS_SELECTION functions
 *
 */

#ifndef A2SEL_H_
#define A2SEL_H_

#include <public/adios_selection.h>

ADIOS_SELECTION * a2sel_boundingbox (int ndim, const uint64_t *start, const uint64_t *count);
ADIOS_SELECTION * a2sel_points (int ndim, uint64_t npoints, const uint64_t *points,
                                ADIOS_SELECTION * container, int free_points_on_delete);
ADIOS_SELECTION * a2sel_writeblock (int index);
ADIOS_SELECTION * a2sel_auto (char *hints);
void a2sel_free (ADIOS_SELECTION *sel);

ADIOS_SELECTION * a2sel_copy (const ADIOS_SELECTION * sel);

/** Create a list of N-dimensional global points from a list of 1D offsets in a bounding box defined by 'start' and 'count'.
 *  Box must be ndim dimensional
 *  pts1d is an array of npoints elements
 *  ptsNd is an array of ndim * npoints elements, pre-allocated by caller
 *  If global is 0 the result points are local in the bounding box,
 *  otherwise they are global (i.e. the 'start' is added to all coordinates.
 */
void a2sel_points_1DtoND_box (uint64_t npoints, uint64_t *pts1d,
                                             int ndim, uint64_t *start, uint64_t *count, int global,
                                             uint64_t *ptsNd);

ADIOS_SELECTION * a2sel_points_1DtoND (ADIOS_SELECTION * pointsinbox1D, int global);

#endif /* A2SEL_H_ */
