/*
 * adios_transforms_hooks.h
 *
 *  Created on: Jul 12, 2012
 *      Author: Drew
 */

#ifndef ADIOS_TRANSFORMS_HOOKS_H_
#define ADIOS_TRANSFORMS_HOOKS_H_

#include <stdint.h>
#include "adios_bp_v1.h"
#include "adios_logger.h"
#include "adios_internals.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_write.h"
#include "adios_transforms_util.h"
#include "public/adios_selection.h"
#include "public/adios_error.h"

// Transform registry entry
typedef struct {
    uint16_t (*transform_get_metadata_size)();

    uint64_t (*transform_calc_vars_transformed_size)(uint64_t orig_size, int num_vars);

    int (*transform_apply)(
            struct adios_file_struct *fd, struct adios_var_struct *var,
            uint64_t *transformed_len, int *use_shared_buffer, int *wrote_to_shared_buffer);

    enum ADIOS_ERRCODES (*transform_retrieve_subvolume)(
            struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
            const adios_subvolume_copy_spec *subv_spec,
            void *read_state, adios_transform_var_read_delegate read_delegate,
            enum ADIOS_FLAG swap_endianness);
} adios_transform_method;

/*
 * The transform registry, containing all necessary function pointers for
 * each transform method.
 */
extern adios_transform_method TRANSFORM_METHODS[num_adios_transform_types];

// Initialize the transform system for adios read-only libraries
void adios_transform_read_init();
// Initialize the transform system for adios read/write libraries
void adios_transform_init();

// Delegation functions
uint16_t adios_transform_get_metadata_size(enum ADIOS_TRANSFORM_TYPE transform_type);
uint64_t adios_transform_calc_vars_transformed_size(enum ADIOS_TRANSFORM_TYPE transform_type, uint64_t orig_size, int num_vars);
int adios_transform_apply(
        enum ADIOS_TRANSFORM_TYPE transform_type,
        struct adios_file_struct *fd, struct adios_var_struct *var,
        uint64_t *transformed_len, int *use_shared_buffer, int *wrote_to_shared_buffer);

enum ADIOS_ERRCODES adios_transform_retrieve_subvolume(
        enum ADIOS_TRANSFORM_TYPE transform_type,
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *subv_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness);

// Transport method function declaration/assignment macros
#define FORWARD_DECLARE(tmethod) \
    uint16_t adios_transform_##tmethod##_get_metadata_size();											\
    uint64_t adios_transform_##tmethod##_calc_vars_transformed_size(uint64_t orig_size, int num_vars);	\
    int adios_transform_##tmethod##_apply(struct adios_file_struct *fd, struct adios_var_struct *var,	\
                                          uint64_t *transformed_len,									\
                                          int *use_shared_buffer, int *wrote_to_shared_buffer);			\
    enum ADIOS_ERRCODES adios_transform_##tmethod##_retrieve_subvolume(									\
            struct adios_index_var_struct_v1 *var, void *global_out, int time_index,					\
            const adios_subvolume_copy_spec *copy_spec,													\
            void *read_state, adios_transform_var_read_delegate read_delegate,							\
            enum ADIOS_FLAG swap_endianness);

#define ASSIGN_FNS(tmethod, method_type) \
    TRANSFORM_METHODS[method_type].transform_get_metadata_size = adios_transform_##tmethod##_get_metadata_size;		\
    TRANSFORM_METHODS[method_type].transform_calc_vars_transformed_size = adios_transform_##tmethod##_calc_vars_transformed_size;	\
    TRANSFORM_METHODS[method_type].transform_apply = adios_transform_##tmethod##_apply;								\
    TRANSFORM_METHODS[method_type].transform_retrieve_subvolume = adios_transform_##tmethod##_retrieve_subvolume;

#define ASSIGN_READ_FNS(tmethod, method_type) \
    TRANSFORM_METHODS[method_type].transform_get_metadata_size = 0;													\
    TRANSFORM_METHODS[method_type].transform_calc_vars_transformed_size = adios_transform_##tmethod##_calc_vars_transformed_size;	\
    TRANSFORM_METHODS[method_type].transform_apply = 0;																\
    TRANSFORM_METHODS[method_type].transform_retrieve_subvolume = adios_transform_##tmethod##_retrieve_subvolume;

void unimplemented_transform_function(const char *tmethod, const char *func);

#define DEFINE_FNS_UNIMPL(tmethod) 															\
        uint16_t adios_transform_##tmethod##_get_metadata_size() {							\
            unimplemented_transform_function(#tmethod, __FUNCTION__);						\
            return 0;																		\
        }																					\
        uint64_t adios_transform_##tmethod##_calc_vars_transformed_size(uint64_t orig_size, int num_vars) {	\
            unimplemented_transform_function(#tmethod, __FUNCTION__);						\
            return 0;																		\
        }																					\
        int adios_transform_##tmethod##_apply(struct adios_file_struct *fd,							\
                                              struct adios_var_struct *var,							\
                                              uint64_t *transformed_len,							\
                                              int *use_shared_buffer, int *wrote_to_shared_buffer) {	\
            unimplemented_transform_function(#tmethod, __FUNCTION__);								\
            return 0;																				\
        }																							\
        enum ADIOS_ERRCODES adios_transform_##tmethod##_retrieve_subvolume(							\
                struct adios_index_var_struct_v1 *var, void *global_out, int time_index,			\
                const adios_subvolume_copy_spec *subv_spec,											\
                void *read_state, adios_transform_var_read_delegate read_delegate,			\
                enum ADIOS_FLAG swap_endianness) {													\
            unimplemented_transform_function(#tmethod, __FUNCTION__);								\
            return err_operation_not_supported;														\
        }

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//// SETUP YOUR NEW TRANSFORM METHODS                                      ////
//// 1. Add an entry to enum adios_transform_type in adios_bp_v1.h (be     ////
////    sure to update num_adios_transform_types accordingly               ////
//// 2. Add a FOWARD_DECLARE line below (assuming standard naming)         ////
//// 3. Add an ASSIGN_FNS line to adios_transforms_hooks.c (assuming       ////
////    standard naming)                                                   ////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

FORWARD_DECLARE(none);
FORWARD_DECLARE(identity);
FORWARD_DECLARE(alacrity);

#undef FORWARD_DECLARE

#endif /* ADIOS_TRANSFORMS_HOOKS_H_ */
