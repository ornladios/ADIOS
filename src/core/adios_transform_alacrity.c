#include <stdint.h>
#include <assert.h>

#include "adios_logger.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_hooks_read.h"
#include "adios_subvolume.h"
#include "public/adios_error.h"

#ifdef ALACRITY
#include "alacrity.h"

// *apply function is in adios_transform_alacrity_write.c, since it requires
// access to functions from adios_internals.c, and therefore cannot link in
// read-only libraries.

typedef enum {
    DISJOINT,			// Bounding box and variable PG are completely disjoint
    PARTIALLY_COVERED,	// Bounding box only partially covers variable PG
    COVERED,			// Bounding box completely covers variable PG
    UNKNOWN				// Some unknown selection method
} bb_intersetion_type;

// TODO: Support endianness change!!!
enum ADIOS_ERRCODES adios_transform_alacrity_retrieve_subvolume(
                    struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
                    const adios_subvolume_copy_spec *copy_spec,
                    void *read_state, adios_transform_var_read_delegate read_delegate,
                    enum ADIOS_FLAG swap_endianness) {

    assert(var->characteristics[time_index].transform.transform_type == adios_transform_alacrity);

    // Read the variable in its entirety
    void *buf = read_delegate(var, time_index, 0, adios_transform_var_get_transformed_size(var, time_index), read_state);

    // Deserialize the ALACRITY data
    ALPartitionData part;
    memstream_t ms = memstreamInitReturn(buf);
    ALDeserializePartitionData(&part, &ms);
    memstreamDestroy(&ms, false);

    // Decode the ALACRITY data
    uint64_t out_size;
    void *decode_buf = malloc(ALGetDecodeLength(&part));
    ALDecode(&part, decode_buf, &out_size);
    ALPartitionDataDestroy(&part);

    copy_subvolume_with_spec(global_out, decode_buf, copy_spec, var->characteristics[0].transform.pre_transform_type, swap_endianness);

    return err_no_error;
}

#else
DECLARE_TRANSFORM_READ_METHOD_UNIMPL(alacrity)
#endif
