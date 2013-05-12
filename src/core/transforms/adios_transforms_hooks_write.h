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
uint64_t adios_transform_calc_vars_transformed_size(struct adios_transform_spec *transform_spec, uint64_t orig_size, int num_vars);
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
                struct adios_transform_spec *transform_spec,
                uint64_t orig_size, int num_vars);

    int (*transform_apply)(
            struct adios_file_struct *fd, struct adios_var_struct *var,
            uint64_t *transformed_len, int use_shared_buffer, int *wrote_to_shared_buffer);
} adios_transform_write_method;

// Transform write method registry
extern adios_transform_write_method TRANSFORM_WRITE_METHODS[num_adios_transform_types];

#define DECLARE_TRANSFORM_WRITE_METHOD(tmethod) \
    uint16_t adios_transform_##tmethod##_get_metadata_size(struct adios_transform_spec *transform_spec); \
    uint64_t adios_transform_##tmethod##_calc_vars_transformed_size(struct adios_transform_spec *transform_spec, \
                                                                    uint64_t orig_size, int num_vars);	 \
    int adios_transform_##tmethod##_apply(struct adios_file_struct *fd, struct adios_var_struct *var,	 \
                                          uint64_t *transformed_len,									 \
                                          int use_shared_buffer, int *wrote_to_shared_buffer);

#define UNIMPL_TRANSFORM_WRITE_FN(tmethod, func) \
    adios_error(err_operation_not_supported,								\
                "Transport method %s is not supported for write in this "	\
                "configuration of ADIOS (function missing: %s)\n",			\
                #tmethod, func);

// Note: this is actually a "definition" in the language-semantic sense, but this detail is
//  irrelevant to users, so we name it similarly to DECLARE_TRANSFORM_WRITE_METHOD
#define DECLARE_TRANSFORM_WRITE_METHOD_UNIMPL(tmethod) 										\
        uint16_t adios_transform_##tmethod##_get_metadata_size(struct adios_transform_spec *transform_spec) { \
            UNIMPL_TRANSFORM_WRITE_FN(tmethod, __FUNCTION__);								\
            return 0;																		\
        }																					\
        uint64_t adios_transform_##tmethod##_calc_vars_transformed_size(struct adios_transform_spec *transform_spec, \
                                                                        uint64_t orig_size, int num_vars) {	\
            UNIMPL_TRANSFORM_WRITE_FN(tmethod, __FUNCTION__);								\
            return 0;																		\
        }																					\
        int adios_transform_##tmethod##_apply(struct adios_file_struct *fd,							\
                                              struct adios_var_struct *var,							\
                                              uint64_t *transformed_len,							\
                                              int use_shared_buffer, int *wrote_to_shared_buffer) {	\
            UNIMPL_TRANSFORM_WRITE_FN(tmethod, __FUNCTION__);										\
            return 0;																				\
        }

#define REGISTER_TRANSFORM_WRITE_METHOD(tmethod, method_type) \
    TRANSFORM_WRITE_METHODS[method_type].transform_get_metadata_size = adios_transform_##tmethod##_get_metadata_size;		\
    TRANSFORM_WRITE_METHODS[method_type].transform_calc_vars_transformed_size = adios_transform_##tmethod##_calc_vars_transformed_size;	\
    TRANSFORM_WRITE_METHODS[method_type].transform_apply = adios_transform_##tmethod##_apply;

#endif /* ADIOS_TRANSFORMS_HOOKS_WRITE_H_ */
