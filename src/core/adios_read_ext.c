/*
 * adios_read_ext.c
 *
 *  Created on: Apr 28, 2014
 *      Author: xczou
 */

#include <stdlib.h>

#include "public/adios_read_ext.h"
#include "core/transforms/adios_transforms_common.h"

// Sets the "data view" for this ADIOS file, which determines how ADIOS presents variables through
// adios_inq_var*, and how reads are evaluated in adios_schedule_reads/adios_check_reads calls.
// Currently, the choice is between a logical and physical view of the data, which only differ for
// transformed variables; a logical view of a transformed variable presents the data as it was
// originally written (this is the default), whereas a physical view presents the transformed data
// as it actually exists on disk.
void adios_read_set_data_view(ADIOS_FILE *fp, data_view_t vt) {
	fp->data_view = vt;
}

// Populates data transform information about a given variable into an ADIOS_VARTRANSFORM struct
// Return NULL if failed
ADIOS_VARTRANSFORM *  adios_inq_var_transform(const ADIOS_FILE *fp, const ADIOS_VARINFO *varinfo){
	ADIOS_TRANSINFO* info = adios_read_bp_inq_var_transinfo(fp, varinfo);
	if (info == NULL)
		return NULL;
	ADIOS_VARTRANSFORM *transform = (ADIOS_VARTRANSFORM*) malloc(sizeof(ADIOS_VARTRANSFORM));

	transform->transform_type = info->transform_type;

//	int sumBlocks = varinfo->sum_nblocks;

	transform->transform_metadatas = info->transform_metadatas;
	return transform;

}

void adios_free_var_transform(ADIOS_VARTRANSFORM *vartransform) {

}

// Creates a writeblock selection that only retrieves elements [start_elem, start_elem + num_elems)
// within a variable. An element is a single value of whatever the varaible's datatype is (i.e.,
// 1 element = 1 double if the variable type is double, 1 byte if the variable type is byte, etc.)
ADIOS_SELECTION * adios_selection_writeblock_bounded(int index, uint64_t start_elem, uint64_t num_elems);
