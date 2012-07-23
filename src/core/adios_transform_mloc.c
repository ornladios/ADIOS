#include <stdint.h>
#include <assert.h>

#include "adios_logger.h"
#include "adios_transforms_hooks.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_write.h"
#include "public/adios_error.h"

#ifdef MLOC
#include "mloc.h"

#define MAX_POSSIBLE_BINS 65536
static const int MAX_PART_METADATA_SIZE =
    sizeof(uint64_t) +
    sizeof(unsigned short int) +
    MAX_POSSIBLE_BINS * ( sizeof(unsigned short int) +  2 * sizeof(uint64_t) + sizeof(unsigned char) );


uint64_t adios_transform_mloc_vars_size(uint64_t orig_size, int num_vars) {
    return num_vars * MAX_PART_METADATA_SIZE +	// For the metadata
           orig_size * 5/4;						// For the index + data
}

// *apply function is in adios_transform_mloc_write.c, since it requires
// access to functions from adios_internals.c, and therefore cannot link in
// read-only libraries.

typedef enum {
    DISJOINT,			// Bounding box and variable PG are completely disjoint
    PARTIALLY_COVERED,	// Bounding box only partially covers variable PG
    COVERED,			// Bounding box completely covers variable PG
    UNKNOWN				// Some unknown selection method
} bb_intersetion_type;

enum ADIOS_ERRCODES adios_transform_mloc_retrieve_subvolume(
                    struct adios_index_var_struct_v1 *var, void *out,
                    uint64_t *out_size, void *read_state,
                    adios_transform_var_read_delegate read_delegate,
                    ADIOS_SELECTION *sel, int time_index) {

    assert(var->characteristics[time_index].transform_type == adios_transform_mloc);
    struct adios_index_characteristic_dims_struct_v1 *orig_dims =
            &var->characteristics[time_index].pre_transform_dimensions;

 
    return 0;
}






#else
DEFINE_FNS_UNIMPL(mloc);

#endif
