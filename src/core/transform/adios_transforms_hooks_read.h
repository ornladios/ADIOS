/*
 * adios_transforms_hooks_read.h
 *
 *  Created on: Jul 25, 2012
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_TRANSFORMS_HOOKS_READ_H_
#define ADIOS_TRANSFORMS_HOOKS_READ_H_

#include <stdint.h>
#include "adios_bp_v1.h"
#include "adios_subvolume.h"
#include "adios_transforms_common.h"
#include "adios_transforms_read.h"
#include "adios_transforms_reqgroup.h"
#include "adios_read_transformed.h"
#include "public/adios_read.h"
#include "public/adios_selection.h"
#include "public/adios_error.h"

enum ADIOS_READ_RESULT_MODE {
    adios_read_noreturn,
    adios_read_return_complete,
    adios_read_return_partial
};

// Initialize the transform system for adios read-only libraries
void adios_transform_read_init();

// Delegation functions
adios_transform_read_reqgroup * adios_transform_generate_read_reqgroup(const ADIOS_VARINFO *vi, const ADIOS_TRANSINFO* ti, const ADIOS_FILE *fp,
                                                                       const ADIOS_SELECTION *sel, int from_steps, int nsteps, void *data);

ADIOS_VARCHUNK * adios_transform_subrequest_completed(adios_transform_read_reqgroup *reqgroup,
                                                      adios_transform_pg_reqgroup *pg_reqgroup,
                                                      adios_transform_read_subrequest *completed_subreq,
                                                      enum ADIOS_READ_RESULT_MODE mode);

ADIOS_VARCHUNK * adios_transform_pg_reqgroup_completed(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *completed_pg_reqgroup,
                                                       enum ADIOS_READ_RESULT_MODE mode);

ADIOS_VARCHUNK * adios_transform_read_reqgroup_completed(adios_transform_read_reqgroup *completed_reqgroup,
                                                         enum ADIOS_READ_RESULT_MODE mode);

enum ADIOS_ERRCODES adios_transform_retrieve_subvolume(
        enum ADIOS_TRANSFORM_TYPE transform_type,
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *subv_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness);



////////////////////////////////////////////////
// Transform read method registry/registration
////////////////////////////////////////////////

// Transform read method registry entry
typedef struct {
    enum ADIOS_ERRCODES (*transform_retrieve_subvolume)(
            struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
            const adios_subvolume_copy_spec *subv_spec,
            void *read_state, adios_transform_var_read_delegate read_delegate,
            enum ADIOS_FLAG swap_endianness);

    int (*transform_generate_read_subrequests)(
            adios_transform_read_reqgroup *reqgroup,
            adios_transform_pg_reqgroup *pg_reqgroup);

    ADIOS_VARCHUNK * (*transform_subrequest_completed)(
            adios_transform_read_reqgroup *reqgroup,
            adios_transform_pg_reqgroup *pg_reqgroup,
            adios_transform_read_subrequest *completed_subreq,
            enum ADIOS_READ_RESULT_MODE mode);

    ADIOS_VARCHUNK * (*transform_pg_reqgroup_completed)(
            adios_transform_read_reqgroup *reqgroup,
            adios_transform_pg_reqgroup *completed_pg_reqgroup,
            enum ADIOS_READ_RESULT_MODE mode);

    ADIOS_VARCHUNK * (*transform_reqgroup_completed)(
            adios_transform_read_reqgroup *completed_reqgroup,
            enum ADIOS_READ_RESULT_MODE mode);
} adios_transform_read_method;

// Transform read method registry
extern adios_transform_read_method TRANSFORM_READ_METHODS[num_adios_transform_types];

// Transport method function declaration/assignment macros
#define DECLARE_TRANSFORM_READ_METHOD(tmethod) \
    enum ADIOS_ERRCODES adios_transform_##tmethod##_retrieve_subvolume(					\
            struct adios_index_var_struct_v1 *var, void *global_out, int time_index,	\
            const adios_subvolume_copy_spec *copy_spec,									\
            void *read_state, adios_transform_var_read_delegate read_delegate,			\
            enum ADIOS_FLAG swap_endianness);											\
            int adios_transform_##tmethod##_generate_read_subrequests(	\
                    adios_transform_read_reqgroup *reqgroup,			\
                    adios_transform_pg_reqgroup *pg_reqgroup);			\
            ADIOS_VARCHUNK * adios_transform_##tmethod##_subrequest_completed(	\
                    adios_transform_read_reqgroup *reqgroup,					\
                    adios_transform_pg_reqgroup *pg_reqgroup,					\
                    adios_transform_read_subrequest *completed_subreq,			\
                    enum ADIOS_READ_RESULT_MODE mode);							\
            ADIOS_VARCHUNK * adios_transform_##tmethod##_pg_reqgroup_completed(	\
                    adios_transform_read_reqgroup *reqgroup,					\
                    adios_transform_pg_reqgroup *completed_pg_reqgroup,			\
                    enum ADIOS_READ_RESULT_MODE mode);							\
            ADIOS_VARCHUNK * adios_transform_##tmethod##_reqgroup_completed(	\
                    adios_transform_read_reqgroup *completed_reqgroup,			\
                    enum ADIOS_READ_RESULT_MODE mode);

#define UNIMPL_TRANSFORM_READ_FN(tmethod, func) \
    adios_error(err_operation_not_supported,								\
                "Transport method %s is not supported for read in this "	\
                "configuration of ADIOS (function missing: %s)\n",			\
                #tmethod, func);

#define DECLARE_TRANSFORM_READ_METHOD_UNIMPL(tmethod) 									\
    enum ADIOS_ERRCODES adios_transform_##tmethod##_retrieve_subvolume(					\
            struct adios_index_var_struct_v1 *var, void *global_out, int time_index,	\
            const adios_subvolume_copy_spec *subv_spec,									\
            void *read_state, adios_transform_var_read_delegate read_delegate,			\
            enum ADIOS_FLAG swap_endianness) {											\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);								\
        return err_operation_not_supported;												\
    }																					\
    int adios_transform_##tmethod##_generate_read_subrequests(	\
            adios_transform_read_reqgroup *reqgroup,			\
            adios_transform_pg_reqgroup *pg_reqgroup) {			\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);		\
        return adios_errno;										\
    }															\
    ADIOS_VARCHUNK * adios_transform_##tmethod##_subrequest_completed(	\
            adios_transform_read_reqgroup *reqgroup,					\
            adios_transform_pg_reqgroup *pg_reqgroup,					\
            adios_transform_read_subrequest *completed_subreq,			\
            enum ADIOS_READ_RESULT_MODE mode) {							\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);				\
        return NULL;													\
    }																	\
    ADIOS_VARCHUNK * adios_transform_##tmethod##_pg_reqgroup_completed(	\
            adios_transform_read_reqgroup *reqgroup,					\
            adios_transform_pg_reqgroup *completed_pg_reqgroup,			\
            enum ADIOS_READ_RESULT_MODE mode) {							\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);				\
        return NULL;													\
    }																	\
    ADIOS_VARCHUNK * adios_transform_##tmethod##_reqgroup_completed(	\
            adios_transform_read_reqgroup *completed_reqgroup,			\
            enum ADIOS_READ_RESULT_MODE mode) {							\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);				\
        return NULL;													\
    }

#define REGISTER_TRANSFORM_READ_METHOD(tmethod, method_type) \
        TRANSFORM_READ_METHODS[method_type].transform_retrieve_subvolume = adios_transform_##tmethod##_retrieve_subvolume;					\
        TRANSFORM_READ_METHODS[method_type].transform_generate_read_subrequests = adios_transform_##tmethod##_generate_read_subrequests;	\
        TRANSFORM_READ_METHODS[method_type].transform_subrequest_completed = adios_transform_##tmethod##_subrequest_completed;				\
        TRANSFORM_READ_METHODS[method_type].transform_pg_reqgroup_completed = adios_transform_##tmethod##_pg_reqgroup_completed;			\
        TRANSFORM_READ_METHODS[method_type].transform_reqgroup_completed = adios_transform_##tmethod##_reqgroup_completed;

#endif /* ADIOS_TRANSFORMS_HOOKS_READ_H_ */
