
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "adios_bp_v1.h"
#include "common_adios.h"
#include "adios_logger.h"
#include "adios_internals.h"
#include "public/adios_types.h"
#include "util.h"

#include "adios_transforms_hooks_write.h"
#include "adios_transforms_common.h"
#include "adios_transforms_write.h"
#include "adios_transforms_specparse.h"

////////////////////////////////////////
// adios_group_size support
////////////////////////////////////////
uint64_t adios_transform_worst_case_transformed_group_size(uint64_t group_size, struct adios_file_struct *fd)
{
    // New group size is always at least the original group size
    uint64_t max_transformed_group_size = group_size;
    int non_scalar_transformed_var_count = 0;

    struct adios_var_struct *cur_var;
    uint64_t transformed_group_size;
    int transform_type;

    // Table of what transform types have been seen so far.
    // Allocated on stack; no dynamic memory to clean up.
    int transform_type_seen[num_adios_transform_types];
    memset(transform_type_seen, 0, num_adios_transform_types * sizeof(int));

    // Identify all transform methods used, and count the number of non-scalar
    // variables
    for (cur_var = fd->group->vars; cur_var; cur_var = cur_var->next)
    {
        // Skip unknown/none transform types and scalar variables
        if (cur_var->transform_type == adios_transform_none ||
            !cur_var->dimensions)
        {
            continue;
        }

        transform_type_seen[cur_var->transform_type] = 1;
        non_scalar_transformed_var_count++;
    }

    // For each transform type, get a worst-case group size estimate, and
    // record the worst of the worst cases
    for (transform_type = adios_transform_none + 1; transform_type < num_adios_transform_types; transform_type++) {
        if (!transform_type_seen[transform_type])
            continue;

        transformed_group_size = adios_transform_calc_vars_transformed_size(transform_type, group_size, non_scalar_transformed_var_count);

        if (transformed_group_size > max_transformed_group_size) {
            max_transformed_group_size = transformed_group_size;
        }
    }

    // Return the maximum worst case for the group size. Note that this is
    // always at least group_size, since it is initialized to that value,
    // and never decreases.
    return max_transformed_group_size;
}

////////////////////////////////////////
// Variable conversion to byte array (preparation for transform)
////////////////////////////////////////
static struct adios_dimension_struct * new_dimension() {
    struct adios_dimension_struct *dim = (struct adios_dimension_struct *)malloc(sizeof(struct adios_dimension_struct));
    dim->dimension.rank = 0;
    dim->dimension.var = NULL;
    dim->dimension.time_index = adios_flag_no;
    dim->global_dimension.rank = 0;
    dim->global_dimension.var = NULL;
    dim->global_dimension.time_index = adios_flag_unknown;
    dim->local_offset.rank = 0;
    dim->local_offset.var = NULL;
    dim->local_offset.time_index = adios_flag_unknown;
    dim->next = 0;
    return dim;
}

static int find_time_dimension(struct adios_dimension_struct *dim, struct adios_dimension_struct **time_dim, enum ADIOS_FLAG fortran_order_flag) {
    struct adios_dimension_struct *cur_dim;
    int i;
    for (i = 0, cur_dim = dim; cur_dim; cur_dim = cur_dim->next, i++ ) {
        if (cur_dim->dimension.time_index == adios_flag_yes) {
            if (time_dim) *time_dim = cur_dim;
            return i;
        }

        // If we're at the last dimension, and we haven't found a time
        // dimension yet, check whether the global array is zero; if so, we
        // have a time dimension, and must infer the location based on whether
        // this is FORTRAN or not
        if (!cur_dim->next) {
            if (cur_dim->global_dimension.var == NULL && cur_dim->global_dimension.rank == 0) {
                if (fortran_order_flag == adios_flag_yes) {
                    if (time_dim) *time_dim = cur_dim;
                    return i;
                } else if (fortran_order_flag == adios_flag_no) {
                    if (time_dim) *time_dim = dim;
                    return 0;
                }
            }
        }
    }

    if (time_dim) *time_dim = 0;
    return -1;
}

// If there is a time dimension the final dimensions will look like this:
//   local:  t  l1 l2 l3 (or l1 l2 l3 t if fortran order)
//   global: g1 g2 g3 0
//   offset: o1 o2 o3 0
static void adios_transform_attach_byte_array_dimensions(struct adios_group_struct *grp, struct adios_var_struct *var) {
    int i, new_ndim, new_time_dim_pos;
    uint64_t ldims[3];
    uint64_t gdims[3];
    uint64_t odims[3];

    int orig_ndim = count_dimensions(var->pre_transform_dimensions);
    int orig_time_dim_pos = find_time_dimension(var->pre_transform_dimensions, NULL, grp->adios_host_language_fortran);

    assert(orig_time_dim_pos == -1 || orig_time_dim_pos == 0 || orig_time_dim_pos == orig_ndim - 1); // Time dimension is either first, last, or non-existant

    ldims[0] = 1; // 1 PG ID in length
    ldims[1] = 0; // unknown bytes in length
    gdims[0] = UINT64_MAX >> 1; // Infinite max PGs
    gdims[1] = UINT64_MAX >> 1; // Infinite max bytes
    odims[0] = 0; // unknown PG ID
    odims[1] = 0; // 0 byte offset
    ldims[2] = gdims[2] = odims[2] = 0; // No 3rd dimension, yet

    // If we are writing from a FORTRAN file, we need to reverse the raw dimensions, so that when
    // it is read back it is properly swapped to C order
    if (grp->adios_host_language_fortran == adios_flag_yes) {
        int dummy = -1;
        swap_order(2, ldims, &dummy);
        swap_order(2, gdims, &dummy);
        swap_order(2, odims, &dummy);
    }

    // Add the time dimension
    if (orig_time_dim_pos == 0) {
        ldims[2] = ldims[1];
        ldims[1] = ldims[0];
        ldims[0] = 1;
        new_ndim = 3;
        new_time_dim_pos = 0;
    } else if (orig_time_dim_pos == orig_ndim - 1) {
        ldims[2] = 1;
        new_ndim = 3;
        new_time_dim_pos = 2;
    } else {
        new_ndim = 2;
        new_time_dim_pos = -1;
    }

    // Construct the dimension linked list
    for (i = 0; i < new_ndim; i++) {
        struct adios_dimension_struct *new_dim = new_dimension();
        new_dim->dimension.time_index = (i == new_time_dim_pos) ? adios_flag_yes : adios_flag_no;
        new_dim->dimension.rank = ldims[i];
        new_dim->global_dimension.rank = gdims[i];
        new_dim->local_offset.rank = odims[i];
        adios_append_dimension(&var->dimensions, new_dim);
    }
}

static void adios_transform_convert_var_to_byte_array(struct adios_group_struct *grp, struct adios_var_struct *var) {
    // Save old metadata
    var->pre_transform_type = var->type;
    var->pre_transform_dimensions = var->dimensions;

    // Convert the type to byte array and clear the dimensions (since they were
    // moved and shouldn't be double-referenced)
    var->type = adios_byte;
    var->dimensions = 0;

    // Attach the new dimension to the variable
    adios_transform_attach_byte_array_dimensions(grp, var);
}

////////////////////////////////////////
// Definition phase - set up transform parameters
////////////////////////////////////////

struct adios_var_struct * adios_transform_define_var(struct adios_group_struct *orig_var_grp,
                                                     struct adios_var_struct *orig_var,
                                                     struct adios_transform_spec *transform_spec) {
    log_debug("Transforming variable %s with type %d\n", orig_var->name, transform_spec->transform_type);

    // Set transform type and spec
    orig_var->transform_type = transform_spec->transform_type;
    orig_var->transform_spec = transform_spec;

    // If there is no transform, nothing else to do
    if (transform_spec->transform_type == adios_transform_none)
        return orig_var;

    // If we get here, transform_type is an actual transformation, so prepare
    // the variable. This entails 1) adding a new dimension variable for
    // the variable (it will become a 1D byte array), and 2) making the
    // variable into a 1D byte array.

    // Convert variable to 1D byte array
    adios_transform_convert_var_to_byte_array(orig_var_grp, orig_var);
    log_debug("Converted variable %s into byte array\n", orig_var->name);

    // Allocate the transform-specific metadata buffer
    orig_var->transform_metadata_len = adios_transform_get_metadata_size(transform_spec);
    if (orig_var->transform_metadata_len)
        orig_var->transform_metadata = malloc(orig_var->transform_metadata_len);

    orig_var->bitmap = 0; // Disable statistics

    // Return the modified variable
    return orig_var;
}

////////////////////////////////////////
// Write phase - transformed byte array dimension management
////////////////////////////////////////

// NCSU ALACRITY-ADIOS - Compute the pre-transform size of a variable, in bytes
// Precondition: var is a "dimensioned" variable; that is, not a scalar, and not a string
uint64_t adios_transform_get_pre_transform_var_size(struct adios_group_struct *group, struct adios_var_struct *var) {
    assert(var->dimensions);
    assert(var->type != adios_string);
    return adios_get_type_size(var->pre_transform_type, NULL) *
           adios_get_dimension_space_size(var,
                                          var->pre_transform_dimensions,
                                          group);
}

static inline uint64_t generate_unique_block_id(const struct adios_file_struct * fd, const struct adios_var_struct *var) {
    return ((uint64_t)fd->group->process_id << 32) + (uint64_t)var->write_count;
}

static int adios_transform_store_transformed_length(struct adios_file_struct * fd, struct adios_var_struct *var, uint64_t transformed_len) {
    struct adios_dimension_struct *dim1, *dim2, *dim3;
    struct adios_dimension_item_struct *pg_id_offset, *byte_length_ldim;

    const uint64_t pg_id = generate_unique_block_id(fd, var);//fd->pg_start_in_file; // Use the current file offset as a unique ID for this PG

    // Get the first two dimensions (which always exist)
    dim1 = var->dimensions;
    assert(dim1);
    dim2 = dim1->next;
    assert(dim2);

    // Find appropriate dimension items
    if (fd->group->adios_host_language_fortran == adios_flag_yes)
        pg_id_offset = &dim2->local_offset;
    else
        pg_id_offset = &dim1->local_offset;

    if (dim1->dimension.time_index == adios_flag_yes) {
        // If the first dimension is a time dimension, then dimension is
        // upshifted, but only for ->dimension
        dim3 = dim2->next;
        assert(dim3);
        if (fd->group->adios_host_language_fortran == adios_flag_yes)
            byte_length_ldim = &dim2->dimension;
        else
            byte_length_ldim = &dim3->dimension;
    } else {
        if (fd->group->adios_host_language_fortran == adios_flag_yes)
            byte_length_ldim = &dim1->dimension;
        else
            byte_length_ldim = &dim2->dimension;
    }

    // Finally, insert the values into the dimension items
    pg_id_offset->rank = pg_id;
    byte_length_ldim->rank = transformed_len;

    //printf(">>> Statistics bitmap at store-time: %08lx\n", var->bitmap);

    return 1;
}

int adios_transform_variable_data(struct adios_file_struct * fd,
                                  struct adios_var_struct *var,
                                  int use_shared_buffer,
                                  int *wrote_to_shared_buffer) {
    //printf("[TRANSFORM] Would be doing transform ID %d on variable %s here, if it were implemented\n", var->transform_type, var->name);
    //printf("[TRANSFORM] Apparent type of variable is %s, but was originally %s\n", adios_type_to_string_int(var->type), adios_type_to_string_int(var->pre_transform_type));
    //printf("[TRANSFORM] Original size of variable data is %llu\n", adios_transform_get_pre_transform_var_size(fd->group, var));

    assert(fd);
    assert(var);

    if (var->transform_type == adios_transform_none) {
        // If no transform is required, do nothing, and delegate payload
        // writing to the caller by not using the shared buffer
        *wrote_to_shared_buffer = 0;
        return 1;
    }

    assert(var->type == adios_byte); // Assume byte array
    assert(var->transform_type != adios_transform_none);

#if defined(WITH_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
    timer_start ("adios_transform_apply");
#endif
    // Transform the data, get the new length
    uint64_t transformed_len;
    int success = adios_transform_apply(fd, var, &transformed_len, use_shared_buffer, wrote_to_shared_buffer);
#if defined(WITH_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
    timer_stop ("adios_transform_apply");
#endif

    if (!success)
        return 0;

    // Store the new length in the metadata
    adios_transform_store_transformed_length(fd, var, transformed_len);

    return 1;
}



//////////////////////////////////////////////////////////////////////////
// Characteristic support, used in adios_internals.c so cannot be
// in the adios_transforms_common.c, as that isn't always linked
// with adios_internals.c.
//////////////////////////////////////////////////////////////////////////

// TODO: There are 3 exact copies of this function:
//       1) here, 2) core/adios_internals.c, 3) write/adios_adaptive.c
//       Someone should extract to a single copy and make it public via
//       a header file.
static void buffer_write (char ** buffer, uint64_t * buffer_size
                         ,uint64_t * buffer_offset
                         ,const void * data, uint64_t size
                         )
{
    if (*buffer_offset + size > *buffer_size || *buffer == 0)
    {
        char * b = realloc (*buffer, *buffer_offset + size + 1000);
        if (b)
        {
            *buffer = b;
            *buffer_size = (*buffer_offset + size + 1000);
        }
        else
        {
            fprintf (stderr, "Cannot allocate memory in buffer_write.  "
                             "Requested: %llu\n", *buffer_offset + size + 1000);

            return;
        }
    }

    memcpy (*buffer + *buffer_offset, data, size);
    *buffer_offset += size;
}

// Init
int adios_transform_init_transform_var(struct adios_var_struct *var) {
    var->transform_type = adios_transform_none;
    var->transform_spec = 0;
    var->pre_transform_dimensions = 0;
    var->pre_transform_type = adios_unknown;
    //var->transform_type_param_len = 0;
    //var->transform_type_param = 0;
    var->transform_metadata_len = 0;
    var->transform_metadata = 0;
    return 1;
}

// Serialize
static void adios_transform_dereference_dimensions_characteristic(struct adios_index_characteristic_dims_struct_v1 *dst_char_dims, const struct adios_dimension_struct *src_var_dims) {
    uint8_t i;
    uint8_t c = count_dimensions(src_var_dims);

    dst_char_dims->count = c;
    dst_char_dims->dims = malloc(3 * 8 * c); // (local, global, local offset) * count
    assert(dst_char_dims->dims);
    uint64_t *ptr = dst_char_dims->dims;
    for (i = 0; i < c; i++)
    {
        //  Casts to eliminate const-ness problems
        ptr[0] = adios_get_dim_value((struct adios_dimension_item_struct *)&src_var_dims->dimension);
        ptr[1] = adios_get_dim_value((struct adios_dimension_item_struct *)&src_var_dims->global_dimension);
        ptr[2] = adios_get_dim_value((struct adios_dimension_item_struct *)&src_var_dims->local_offset);
        src_var_dims = src_var_dims->next;
        ptr += 3; // Go to the next set of 3
    }
}

static void adios_transform_dereference_dimensions_var(struct adios_dimension_struct **dst_var_dims, const struct adios_dimension_struct *src_var_dims) {
    uint8_t i;
    uint8_t c = count_dimensions(src_var_dims);

    for (i = 0; i < c; i++) {
        struct adios_dimension_struct * d_new =
            (struct adios_dimension_struct *)malloc(sizeof (struct adios_dimension_struct));

        // de-reference dimension id
        d_new->dimension.var = NULL;
        d_new->dimension.rank = adios_get_dim_value(&src_var_dims->dimension);
        d_new->dimension.time_index = src_var_dims->dimension.time_index;
        d_new->global_dimension.var = NULL;
        d_new->global_dimension.rank = adios_get_dim_value(&src_var_dims->global_dimension);
        d_new->global_dimension.time_index = src_var_dims->global_dimension.time_index;
        d_new->local_offset.var = NULL;
        d_new->local_offset.rank = adios_get_dim_value(&src_var_dims->local_offset);
        d_new->local_offset.time_index = src_var_dims->local_offset.time_index;
        d_new->next = 0;

        adios_append_dimension(dst_var_dims, d_new);

        src_var_dims = src_var_dims->next;
    }
}

static void adios_transform_clean_dimensions(struct adios_index_characteristic_dims_struct_v1 *dst_char_dims) {
    dst_char_dims->count = 0;
    if (dst_char_dims->dims)
        free(dst_char_dims->dims);
    dst_char_dims->dims = 0;
}

static uint8_t adios_transform_serialize_transform(enum ADIOS_TRANSFORM_TYPE transform_type,
                                                   enum ADIOS_DATATYPES pre_transform_type,
                                                   const struct adios_index_characteristic_dims_struct_v1 *pre_transform_dimensions,
                                                   uint16_t transform_metadata_len,
                                                   void *transform_metadata,
                                                   uint64_t *write_length, char **buffer, uint64_t *buffer_size, uint64_t *buffer_offset) {

    // Either there is no metadata, or the metadata buffer is non-null
    assert(!transform_metadata_len || transform_metadata);

    *write_length = 0;

    // No transform case
    if (transform_type == adios_transform_none)
        return 0;

    // From this point on, this is an actual transform type

    // Write characteristic flag
    uint8_t flag = (uint8_t)adios_characteristic_transform_type;
    buffer_write(buffer, buffer_size, buffer_offset, &flag, 1);
    *write_length += 1;

    // Write transform type
    buffer_write(buffer, buffer_size, buffer_offset, &transform_type, 1);
    *write_length += 1;

    // Write the pre-transform datatype
    buffer_write(buffer, buffer_size, buffer_offset, &pre_transform_type, 1);
    *write_length += 1;

    // Write the number of pre-transform dimensions
    buffer_write (buffer, buffer_size, buffer_offset, &pre_transform_dimensions->count, 1);
    *write_length += 1;

    // Write the length of pre-transform dimension data about to be written
    uint16_t len = 3 * 8 * pre_transform_dimensions->count;
    buffer_write (buffer, buffer_size, buffer_offset, &len, 2);
    *write_length += 2;

    // Write the pre-transform dimensions
    buffer_write (buffer, buffer_size, buffer_offset, pre_transform_dimensions->dims, len);
    *write_length += len;

    buffer_write (buffer, buffer_size, buffer_offset, &transform_metadata_len, 2);
    *write_length += 2;

    if (transform_metadata_len) {
        buffer_write (buffer, buffer_size, buffer_offset, transform_metadata, transform_metadata_len);
        *write_length += transform_metadata_len;
    }

    return 1; // Return that we wrote 1 characteristic flag
}

uint8_t adios_transform_serialize_transform_characteristic(const struct adios_index_characteristic_transform_struct *transform, uint64_t *write_length,
                                                           char **buffer, uint64_t *buffer_size, uint64_t *buffer_offset) {
    return adios_transform_serialize_transform(
                transform->transform_type, transform->pre_transform_type, &transform->pre_transform_dimensions,
                transform->transform_metadata_len, transform->transform_metadata,
                write_length, buffer, buffer_size, buffer_offset
           );
}

uint8_t adios_transform_serialize_transform_var(const struct adios_var_struct *var, uint64_t *write_length,
                                                char **buffer, uint64_t *buffer_size, uint64_t *buffer_offset) {

    // In this case, we are going to actually serialize the dimensions as a
    // adios_index_characteristic_dims_struct_v1, but it is currently in the
    // form of an adios_dimension_struct. We must convert here before passing
    // to the common serialization routine.

    struct adios_index_characteristic_dims_struct_v1 tmp_dims;
    adios_transform_dereference_dimensions_characteristic(&tmp_dims, var->pre_transform_dimensions);

    // Perform the serialization using the common function with the temp dimension structure
    uint8_t char_write_count =
            adios_transform_serialize_transform(
                var->transform_type, var->pre_transform_type, &tmp_dims,
                var->transform_metadata_len, var->transform_metadata,
                write_length, buffer, buffer_size, buffer_offset
            );

    // Free the temp dimension structure
    adios_transform_clean_dimensions(&tmp_dims);

    return char_write_count;
}

// Clear
int adios_transform_clear_transform_var(struct adios_var_struct *var) {
    var->transform_type = adios_transform_none;
    if (var->transform_spec)
        adios_transform_free_spec(&var->transform_spec); // Also clears to 0

    var->pre_transform_type = 0;

    // Frees and zeros-out pre_transform_dimensions (since last ->next is 0)
    while (var->pre_transform_dimensions)
    {
        struct adios_dimension_struct * dimensions = var->pre_transform_dimensions->next;
        free(var->pre_transform_dimensions);
        var->pre_transform_dimensions = dimensions;
    }

    // Free/clear transform-specific metadata
    var->transform_metadata_len = 0;
    if (var->transform_metadata)
        free(var->transform_metadata);
    var->transform_metadata = 0;

    return 1; // Return success
}

// Copy
int adios_transform_copy_transform_characteristic(struct adios_index_characteristic_transform_struct *dst_transform, const struct adios_var_struct *src_var) {
    adios_transform_init_transform_characteristic(dst_transform);

    dst_transform->transform_type = src_var->transform_type;
    dst_transform->pre_transform_type = src_var->pre_transform_type;

    adios_transform_dereference_dimensions_characteristic(&dst_transform->pre_transform_dimensions, src_var->pre_transform_dimensions);

    dst_transform->transform_metadata_len = src_var->transform_metadata_len;
    if (src_var->transform_metadata_len) {
        dst_transform->transform_metadata = malloc(src_var->transform_metadata_len);
        memcpy(dst_transform->transform_metadata, src_var->transform_metadata, src_var->transform_metadata_len);
    } else {
        dst_transform->transform_metadata = 0;
    }

    return 1;
}

int adios_transform_copy_var_transform(struct adios_var_struct *dst_var, const struct adios_var_struct *src_var) {
    adios_transform_init_transform_var(dst_var);

    dst_var->transform_type = src_var->transform_type;
    dst_var->pre_transform_type = src_var->pre_transform_type;

    adios_transform_dereference_dimensions_var(&dst_var->pre_transform_dimensions, src_var->pre_transform_dimensions);

    // for parameter
    dst_var->transform_spec = adios_transform_spec_copy(src_var->transform_spec);

    dst_var->transform_metadata_len = src_var->transform_metadata_len;
    if (src_var->transform_metadata_len) {
        dst_var->transform_metadata = malloc(src_var->transform_metadata_len);
        memcpy(dst_var->transform_metadata, src_var->transform_metadata, src_var->transform_metadata_len);
    } else {
        dst_var->transform_metadata = 0;
    }

    return 1;
}

// Calculate overhead
uint64_t adios_transform_calc_transform_characteristic_overhead(struct adios_var_struct *var) {
    if (var->transform_type == adios_transform_none) {
        return 0; // No overhead needed, since characteristic won't be written
    } else {
        return 1 +    // For characterstic flag
               1 +    // For transform_type field
               1 +  // For pre_transform_type field
               adios_calc_var_characteristics_dims_overhead(var->pre_transform_dimensions) + // For pre-transform dimensions field
               2 +    // For transform_metadata_len
               var->transform_metadata_len;    // For transform_metadata
    }
}
