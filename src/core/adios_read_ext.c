/*
 * adios_read_ext.c
 *
 *  Created on: Apr 28, 2014
 *      Author: xczou
 */

#include <stdlib.h>

#include "public/adios_read_ext.h"
#include "core/transforms/adios_transforms_common.h"

// Sets whether the transform layer is enabled for the given ADIOS file descriptor. The default is
// enabled (true) for newly-opened file descriptors, in which case all Read API functions operate
// on the user's view of the data.
//
// If the transform layer is disabled, the internal, 1D-byte-array nature of transformed variables,
// as seen by the at the transport layer, will be exposed via Read API functions. adios_inq_var
// and adios_inq_var_blockinfo will return this internal view of transformed variables, and
// adios_schedule_read will operate on the 1D byte local array coordinate space, meaning only
// ADIOS_WRITEBLOCK selections are useful, and raw, transformed data will be returned into the user
// buffer by these read functions.
//
// Note: Read API functions operating on non-transformed variables are completely unaffected by this flag
void adios_read_set_transforms_enabled(ADIOS_FILE *fp, int enabled);

// Creates a writeblock selection that only retrieves elements [start_elem, start_elem + num_elems)
// within a variable. An element is a single value of whatever the varaible's datatype is (i.e.,
// 1 element = 1 double if the variable type is double, 1 byte if the variable type is byte, etc.)
ADIOS_SELECTION * adios_selection_writeblock_bounded(int index, uint64_t start_elem, uint64_t num_elems);


// Returns the user's view of a given variable, regardless of whether transforms are enabled for the given ADIOS file descriptor
ADIOS_TRANSINFO * adios_inq_transinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi);

// Stores the user's view of a given variable's PG bounds into ti->orig_blockinfo, regardless of
// whether transforms are enabled for the given ADIOS file descriptor
// Returns >0 on success
int adios_inq_trans_blockinfo(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO * ti);

// Returns the transform-specific metadata associated with all PGs of a variable into
// ti->transform_metadatas and ti->transform_metadatas
// Returns >0 on success
int adios_inq_trans_metadata(const ADIOS_FILE *fp, const ADIOS_VARINFO *vi, ADIOS_TRANSINFO * ti);
