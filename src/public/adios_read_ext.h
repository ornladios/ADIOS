/*
 * adios_read_ext.h
 *
 *  Created on: Apr 25, 2014
 *      Author: David A. Boyuka II
 */
#ifndef ADIOS_READ_EXT_H_
#define ADIOS_READ_EXT_H_

#ifdef ADIOS_USE_READ_API_1
#  error "Extended ADIOS Read API is not available in Read API v1"
#endif

#include <stdint.h>

#include "public/adios_read_v2.h"
#include "public/adios_selection.h"

// An opaque type defining a particular view of the data.
// Currently, there are only two possible values: LOGICAL_DATA_VIEW and PHYSICAL_DATA_VIEW
typedef void* data_view_t;

// LOGICAL_DATA_VIEW: the default, ADIOS presents the same view of the data as it was written to the
//   file (e.g., if processes wrote to a 3D global array of doubles, the user API will present a 3D global
//   array of doubles to the user when in this view mode).
extern const data_view_t LOGICAL_DATA_VIEW;

// PHYSICAL_DATA_VIEW: ADIOS will present the raw transport layer view of the data. If a variable is not
//   transformed, its presentation is equivalent to that under LOGICAL_DATA_VIEW. If a variable is
//   transformed, it will be presented as a 1D byte array, and reads will be answered directly from the
//   transformed data with out any de-transformation applied.
extern const data_view_t PHYSICAL_DATA_VIEW;

// An identifier for a particular transform type (e.g., identity, zlib, etc.)
// constant NO_TRANSFORM indicates the absence of any data transform being applied.
typedef int adios_transform_type_t;

extern const adios_transform_type_t NO_TRANSFORM;

typedef struct {
	void *content;
	uint64_t length;
} ADIOS_TRANSFORM_METADATA;

// A transform information structure describing how a particular variable has been transformed
typedef struct {
	const adios_transform_type_t transform_type; /* The data transform applied to this variable */
	const int should_free_transform_metadata; // Used internally for free

	/* An array of transform plugin-specific metadata buffers, one for each
       varblock in this file (number of varblocks == ADIOS_VARINFO.sum_nblocks).
	   Only needed by advanced applications requiring direct manipulation
	   of transformed data. */
	const ADIOS_TRANSFORM_METADATA *transform_metadatas;
} ADIOS_VARTRANSFORM;

// Sets the "data view" for this ADIOS file, which determines how ADIOS presents variables through
// adios_inq_var*, and how reads are evaluated in adios_schedule_reads/adios_check_reads calls.
// Currently, the choice is between a logical and physical view of the data, which only differ for
// transformed variables; a logical view of a transformed variable presents the data as it was
// originally written (this is the default), whereas a physical view presents the transformed data
// as it actually exists on disk.
void adios_read_set_data_view(ADIOS_FILE *fp, data_view_t vt);

// Populates data transform information about a given variable into an ADIOS_VARTRANSFORM struct
// Return NULL if failed
ADIOS_VARTRANSFORM *  adios_inq_var_transform(const ADIOS_FILE *fp, const ADIOS_VARINFO *varinfo);

// Frees an ADIOS_VARTRANSFORM previous received from a call to adios_inq_var_transform
void adios_free_var_transform(const ADIOS_VARINFO *vi, ADIOS_VARTRANSFORM *vartransform);

// Creates a writeblock selection that only retrieves elements [start_elem, start_elem + num_elems)
// within a variable. An element is a single value of whatever the varaible's datatype is (i.e.,
// 1 element = 1 double if the variable type is double, 1 byte if the variable type is byte, etc.)
ADIOS_SELECTION * adios_selection_writeblock_bounded(int index, uint64_t start_elem, uint64_t num_elems);

#endif /* ADIOS_READ_EXT_H_ */
