/*
 * adios_read_ext.c
 *
 *  Created on: Apr 28, 2014
 *      Author: xczou
 */

#include <stdlib.h>

#include "public/adios_read_ext.h"
#include "core/transforms/adios_transforms_common.h"

// Ensure unique pointer-based values for each one
const data_view_t LOGICAL_DATA_VIEW = &LOGICAL_DATA_VIEW;
const data_view_t PHYSICAL_DATA_VIEW = &PHYSICAL_DATA_VIEW;

extern const adios_transform_type_t NO_TRANSFORM = adios_transform_none;

// Sets the "data view" for this ADIOS file, which determines how ADIOS presents variables through
// adios_inq_var*, and how reads are evaluated in adios_schedule_reads/adios_check_reads calls.
// Currently, the choice is between a logical and physical view of the data, which only differ for
// transformed variables; a logical view of a transformed variable presents the data as it was
// originally written (this is the default), whereas a physical view presents the transformed data
// as it actually exists on disk.
void adios_read_set_data_view(ADIOS_FILE *fp, data_view_t data_view) {
	common_read_set_data_view(fp, data_view);
}

// Populates data transform information about a given variable into an ADIOS_VARTRANSFORM struct
// Return NULL if failed
ADIOS_VARTRANSFORM *  adios_inq_var_transform(const ADIOS_FILE *fp, const ADIOS_VARINFO *varinfo){
	ADIOS_TRANSINFO* info = common_read_inq_transinfo(fp, varinfo);
	if (info == NULL)
		return NULL;
	ADIOS_VARTRANSFORM *vartransform = (ADIOS_VARTRANSFORM*) malloc(sizeof(ADIOS_VARTRANSFORM));

	*vartransform = (ADIOS_VARTRANSFORM){
		.varid = varinfo->varid,
		.sum_nblocks = varinfo->sum_nblocks,
		.transform_type = info->transform_type,
		.should_free_transform_metadata = info->should_free_transform_metadata,
		.transform_metadatas = info->transform_metadatas
	};

	//	vartransform->transform_type = info->transform_type;

	// Transfer ownership of the transform_metadatas array to the new struct
	//	vartransform->should_free_transform_metadata = info->should_free_transform_metadata;
	//	vartransform->transform_metadatas = info->transform_metadatas; // TODO: Make this work without casting
	info->transform_metadatas = NULL;
	common_read_free_transinfo(varinfo, info);
	return vartransform;
}

#define MYFREE(p) {if (p){ free((void*)(p)); (p) = NULL; }}
void adios_free_var_transform(ADIOS_VARTRANSFORM *vartransform) {
	if (vartransform->transform_metadatas) {
		if (vartransform->should_free_transform_metadata) {
			int i;
			for (i = 0; i < vartransform->sum_nblocks; i++)
				MYFREE(vartransform->transform_metadatas[i].content);
		}
		MYFREE(vartransform->transform_metadatas);
	}
	MYFREE(vartransform);
}
#undef MYFREE

// Creates a writeblock selection that only retrieves elements [start_elem, start_elem + num_elems)
// within a variable. An element is a single value of whatever the varaible's datatype is (i.e.,
// 1 element = 1 double if the variable type is double, 1 byte if the variable type is byte, etc.)
ADIOS_SELECTION * adios_selection_writeblock_bounded(int index, uint64_t start_elem, uint64_t num_elems, int is_timestep_relative) {
	ADIOS_SELECTION *sel = common_read_selection_writeblock(index);
	sel->u.block.is_absolute_index = !is_timestep_relative;
	sel->u.block.is_sub_pg_selection = 1;
	sel->u.block.element_offset = start_elem;
	sel->u.block.nelements = num_elems;
	return sel;
}

