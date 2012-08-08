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


typedef struct {
    enum ADIOS_DATATYPES elem_type;
    const ADIOS_SELECTION *bounds;
    uint64_t ragged_offset;
    void *data;
} adios_datablock;

// Initialize the transform system for adios read-only libraries
void adios_transform_read_init();

// Datablock management
adios_datablock * adios_datablock_new(
        enum ADIOS_DATATYPES elem_type,
        const ADIOS_SELECTION *bounds,
        void *data);

adios_datablock * adios_datablock_new_ragged(
        enum ADIOS_DATATYPES elem_type,
        const ADIOS_SELECTION *bounds,
        const uint64_t *ragged_offsets, void *data);

adios_datablock * adios_datablock_new_ragged_offset(
        enum ADIOS_DATATYPES elem_type,
        const ADIOS_SELECTION *bounds,
        uint64_t ragged_offset, void *data);

void adios_datablock_free(adios_datablock *datablock, int free_data);

// Delegation functions
adios_transform_read_reqgroup * adios_transform_generate_read_reqgroup(const ADIOS_VARINFO *vi, const ADIOS_TRANSINFO* ti, const ADIOS_FILE *fp,
                                                                       const ADIOS_SELECTION *sel, int from_steps, int nsteps, void *data);

adios_datablock * adios_transform_subrequest_completed(adios_transform_read_reqgroup *reqgroup,
                                                       adios_transform_pg_reqgroup *pg_reqgroup,
                                                       adios_transform_read_subrequest *completed_subreq);

adios_datablock * adios_transform_pg_reqgroup_completed(adios_transform_read_reqgroup *reqgroup,
                                                        adios_transform_pg_reqgroup *completed_pg_reqgroup);

adios_datablock * adios_transform_read_reqgroup_completed(adios_transform_read_reqgroup *completed_reqgroup);


////////////////////////////////////////////////
// Transform read method registry/registration
////////////////////////////////////////////////

// Transform read method registry entry
typedef struct {
    int (*transform_generate_read_subrequests)(
            adios_transform_read_reqgroup *reqgroup,
            adios_transform_pg_reqgroup *pg_reqgroup);

    adios_datablock * (*transform_subrequest_completed)(
            adios_transform_read_reqgroup *reqgroup,
            adios_transform_pg_reqgroup *pg_reqgroup,
            adios_transform_read_subrequest *completed_subreq);

    adios_datablock * (*transform_pg_reqgroup_completed)(
            adios_transform_read_reqgroup *reqgroup,
            adios_transform_pg_reqgroup *completed_pg_reqgroup);

    adios_datablock * (*transform_reqgroup_completed)(
            adios_transform_read_reqgroup *completed_reqgroup);
} adios_transform_read_method;

// Transform read method registry
extern adios_transform_read_method TRANSFORM_READ_METHODS[num_adios_transform_types];

// Transport method function declaration/assignment macros
#define DECLARE_TRANSFORM_READ_METHOD(tmethod)							\
    int adios_transform_##tmethod##_generate_read_subrequests(			\
            adios_transform_read_reqgroup *reqgroup,					\
            adios_transform_pg_reqgroup *pg_reqgroup);					\
   adios_datablock * adios_transform_##tmethod##_subrequest_completed(	\
            adios_transform_read_reqgroup *reqgroup,					\
            adios_transform_pg_reqgroup *pg_reqgroup,					\
            adios_transform_read_subrequest *completed_subreq);			\
   adios_datablock * adios_transform_##tmethod##_pg_reqgroup_completed(	\
            adios_transform_read_reqgroup *reqgroup,					\
            adios_transform_pg_reqgroup *completed_pg_reqgroup);		\
   adios_datablock * adios_transform_##tmethod##_reqgroup_completed(		\
            adios_transform_read_reqgroup *completed_reqgroup);

#define UNIMPL_TRANSFORM_READ_FN(tmethod, func) \
    adios_error(err_operation_not_supported,								\
                "Transport method %s is not supported for read in this "	\
                "configuration of ADIOS (function missing: %s)\n",			\
                #tmethod, func);

#define DECLARE_TRANSFORM_READ_METHOD_UNIMPL(tmethod) 					\
    int adios_transform_##tmethod##_generate_read_subrequests(			\
            adios_transform_read_reqgroup *reqgroup,					\
            adios_transform_pg_reqgroup *pg_reqgroup) {					\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);				\
        return adios_errno;												\
    }																	\
    adios_datablock * adios_transform_##tmethod##_subrequest_completed(	\
            adios_transform_read_reqgroup *reqgroup,					\
            adios_transform_pg_reqgroup *pg_reqgroup,					\
            adios_transform_read_subrequest *completed_subreq) {		\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);				\
        return NULL;													\
    }																	\
    adios_datablock * adios_transform_##tmethod##_pg_reqgroup_completed(	\
            adios_transform_read_reqgroup *reqgroup,					\
            adios_transform_pg_reqgroup *completed_pg_reqgroup) {		\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);				\
        return NULL;													\
    }																	\
    adios_datablock * adios_transform_##tmethod##_reqgroup_completed(	\
            adios_transform_read_reqgroup *completed_reqgroup) {		\
        UNIMPL_TRANSFORM_READ_FN(tmethod, __FUNCTION__);				\
        return NULL;													\
    }

#define REGISTER_TRANSFORM_READ_METHOD(tmethod, method_type) \
        TRANSFORM_READ_METHODS[method_type].transform_generate_read_subrequests = adios_transform_##tmethod##_generate_read_subrequests;	\
        TRANSFORM_READ_METHODS[method_type].transform_subrequest_completed = adios_transform_##tmethod##_subrequest_completed;				\
        TRANSFORM_READ_METHODS[method_type].transform_pg_reqgroup_completed = adios_transform_##tmethod##_pg_reqgroup_completed;			\
        TRANSFORM_READ_METHODS[method_type].transform_reqgroup_completed = adios_transform_##tmethod##_reqgroup_completed;

#endif /* ADIOS_TRANSFORMS_HOOKS_READ_H_ */
