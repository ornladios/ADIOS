/*
 * adios_transforms_hooks_write.h
 *
 *  Created on: Jul 25, 2012
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_TRANSFORMS_HOOKS_WRITE_H_
#define ADIOS_TRANSFORMS_HOOKS_WRITE_H_

#include <stdint.h>
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "adios_transforms_common.h"
#include "public/adios_error.h"

// Initialize the transform system for adios read/write libraries
void adios_transform_init();

// Delegation functions
uint16_t adios_transform_get_metadata_size(struct adios_transform_spec *transform_spec);
uint64_t adios_transform_calc_vars_transformed_size(enum ADIOS_TRANSFORM_TYPE transform_type, uint64_t orig_size, int num_vars);
int adios_transform_apply(
        struct adios_file_struct *fd, struct adios_var_struct *var,
        uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer);

/////////////////////////////////////////////////
// Transform write method registry/registration
/////////////////////////////////////////////////

// Transform write method registry entry
typedef struct {
    uint16_t (*transform_get_metadata_size)(struct adios_transform_spec *transform_spec);

    uint64_t (*transform_calc_vars_transformed_size)(
                enum ADIOS_TRANSFORM_TYPE type,
                uint64_t orig_size, int num_vars);

    int (*transform_apply)(
            struct adios_file_struct *fd, struct adios_var_struct *var,
            uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer);
} adios_transform_write_method;

//
// Every transform plugin has a set of functions that must go through three stages:
// * Declaration: as with a C header
// * Definition: the functions must be defined with bodies (or defined as
//   unimplemented using the DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL utility macro)
// * Registration: loading pointers to the functions into a callback table
//

// Transform method function declarations
#define DECLARE_TRANSFORM_WRITE_METHOD(tmethod) \
    uint16_t adios_transform_##tmethod##_get_metadata_size(struct adios_transform_spec *transform_spec); \
    uint64_t adios_transform_##tmethod##_calc_vars_transformed_size(enum ADIOS_TRANSFORM_TYPE type,      \
                                                                    uint64_t orig_size, int num_vars);   \
    int adios_transform_##tmethod##_apply(struct adios_file_struct *fd, struct adios_var_struct *var,    \
                                          uint64_t *transformed_len,                                     \
                                          int use_shared_buffer, int *wrote_to_shared_buffer);

// Transform method function registration
#define TRANSFORM_WRITE_METHOD_HOOK_LIST(tmethod) \
    adios_transform_##tmethod##_get_metadata_size, \
    adios_transform_##tmethod##_calc_vars_transformed_size, \
    adios_transform_##tmethod##_apply

#define REGISTER_TRANSFORM_WRITE_METHOD_HOOKS(ttable, tmethod, method_type) \
    ttable[method_type] = (adios_transform_write_method){ TRANSFORM_WRITE_METHOD_HOOK_LIST(tmethod) };

// Transform method function helper definitions for unimplemented methods
#define UNIMPL_TRANSFORM_WRITE_FN(tmethod, func) \
    adios_error(err_operation_not_supported,                                \
                "Transform method %s is not supported for write in this "   \
                "configuration of ADIOS (function missing: %s)\n",          \
                #tmethod, func);

// Note: this is actually a "definition" in the language-semantic sense, but this detail is
//  irrelevant to users, so we name it similarly to DECLARE_TRANSFORM_WRITE_METHOD
#define DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(tmethod)                                       \
        uint16_t adios_transform_##tmethod##_get_metadata_size(struct adios_transform_spec *transform_spec) { \
            UNIMPL_TRANSFORM_WRITE_FN(tmethod, __FUNCTION__);                                \
            return 0;                                                                        \
        }                                                                                    \
        uint64_t adios_transform_##tmethod##_calc_vars_transformed_size(enum ADIOS_TRANSFORM_TYPE type,     \
                                                                        uint64_t orig_size, int num_vars) { \
            UNIMPL_TRANSFORM_WRITE_FN(tmethod, __FUNCTION__);                                \
            return 0;                                                                        \
        }                                                                                    \
        int adios_transform_##tmethod##_apply(struct adios_file_struct *fd,                  \
                                              struct adios_var_struct *var,                  \
                                              uint64_t *transformed_len,                     \
                                              int use_shared_buffer, int *wrote_to_shared_buffer) {  \
            UNIMPL_TRANSFORM_WRITE_FN(tmethod, __FUNCTION__);                                \
            return 0;                                                                        \
        }

#endif /* ADIOS_TRANSFORMS_HOOKS_WRITE_H_ */
