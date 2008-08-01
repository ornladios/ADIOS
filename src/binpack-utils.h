#ifndef _BINPACK_UTILS_H
#define _BINPACK_UTILS_H

#include <stdint.h>

/*
 binpack-utils.h
 2006
 */

#include "binpack-general.h"

#define NUM_GP 24
#define MAX_PRIMITIVE_SIZE 16

enum tagbp_t { binpack_group, mesh_group, mesh_info, coord_data, connect_data,
             timestep_group, timestep_info, node_group, node_info,
             nodedata_group, nodedata_info, var_data, basic_group,
             scalar_attr, string_attr};

enum cell_t { point_cell, line_cell, tri_cell, quad_cell, tet_cell, hex_cell,
              pyr_cell, prism_cell };

int bp_getsize (enum vartype_t, void * val);
int bp_getCellInfo (enum cell_t, char *);
int bp_tagPeek (FILE *);
char ** bp_dirparser (char *, int *);

int adios_should_use_data (int element, int rank
                          ,struct adios_bp_dimension_struct * dims
                          ,int * position
                          );
const char * adios_tag_to_string (int tag);
void adios_var_element_count (int rank
                             ,struct adios_bp_dimension_struct * dims
                             ,uint64_t * use_count
                             ,uint64_t * total_count
                             );
#endif
