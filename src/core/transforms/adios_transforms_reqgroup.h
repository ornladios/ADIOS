/*
 * adios_transformed_reqgroup.h
 *
 *  Created on: Jul 30, 2012
 *      Author: David A. Boyuka II
 */

#ifndef ADIOS_TRANSFORMS_REQGROUP_H_
#define ADIOS_TRANSFORMS_REQGROUP_H_

#include <stdint.h>
#include "adios_transforms_common.h"
#include "adios_read_transformed.h"
#include "public/adios_read.h"
#include "public/adios_types.h"

typedef struct _adios_transform_read_subrequest {
    int				completed;

    ADIOS_SELECTION *sel;
    void 			*data;

    void 			*transform_internal;

    struct _adios_transform_read_subrequest *next;
} adios_transform_read_subrequest;

typedef struct _adios_transform_pg_reqgroup {
    int completed;

    int blockidx;
    uint64_t raw_var_length;			// Transformed variable data length, in bytes
    const ADIOS_VARBLOCK *raw_varblock;	// Points into adios_transform_read_reqgroup->varinfo->blockinfo; do not free here
    const ADIOS_VARBLOCK *orig_varblock;// Points into adios_transform_read_reqgroup->transinfo->orig_blockinfo; do not free here

    ADIOS_SELECTION *pg_selection;
    adios_subvolume_copy_spec *pg_intersection_to_global_copyspec;

    int num_subreqs;
    int num_completed_subreqs;
    adios_transform_read_subrequest *subreqs;

    void *transform_internal;

    struct _adios_transform_pg_reqgroup *next;

} adios_transform_pg_reqgroup;

typedef struct _adios_transform_read_reqgroup {
    int completed;
    ADIOS_VARCHUNK *lent_varchunk;

    const ADIOS_FILE      	*fp;

    const ADIOS_VARINFO 	*raw_varinfo;
    const ADIOS_TRANSINFO	*transinfo;
    enum ADIOS_FLAG 		swap_endianness;

    int						from_steps;
    int                     nsteps;
    const ADIOS_SELECTION 	*orig_sel;
    void                  	*orig_data;

    int num_pg_reqgroups;
    int num_completed_pg_reqgroups;
    adios_transform_pg_reqgroup *pg_reqgroups;

    void *transform_internal;

    struct _adios_transform_read_reqgroup *next;
} adios_transform_read_reqgroup;

// adios_transform_read_subrequest manipulation
adios_transform_read_subrequest * adios_transform_new_subreq(ADIOS_SELECTION *sel, void *data);
adios_transform_read_subrequest * adios_transform_new_subreq_byte_segment(const ADIOS_VARBLOCK *raw_varblock, uint64_t start, uint64_t count, void *data);
adios_transform_read_subrequest * adios_transform_new_subreq_whole_pg(const ADIOS_VARBLOCK *raw_varblock, void *data);
void adios_transform_subreq_mark_complete(adios_transform_read_reqgroup *regroup_parent, adios_transform_pg_reqgroup *parent_pg_reqgroup,
                                          adios_transform_read_subrequest *subreq);
void adios_transform_free_subreq(adios_transform_read_subrequest **subreq_ptr);

// adios_transform_pg_reqgroup manipulation
adios_transform_pg_reqgroup * adios_transform_new_pg_reqgroup(int blockidx, const ADIOS_VARBLOCK *orig_varblock,
                                                              const ADIOS_VARBLOCK *raw_varblock, ADIOS_SELECTION *pg_sel,
                                                              adios_subvolume_copy_spec *pg_intersection_to_global_cs);
void adios_transform_pg_reqgroup_append_subreq(adios_transform_pg_reqgroup *pg_reqgroup, adios_transform_read_subrequest *subreq);
int adios_transform_pg_reqgroup_find_subreq(const adios_transform_pg_reqgroup *pg_reqgroup, const ADIOS_VARCHUNK *chunk,
                                            int skip_completed, adios_transform_read_subrequest **matching_subreq);
int adios_transform_pg_reqgroup_remove_subreq(adios_transform_pg_reqgroup *pg_reqgroup, adios_transform_read_subrequest *subreq);
adios_transform_read_subrequest * adios_transform_pg_reqgroup_pop_subreq(adios_transform_pg_reqgroup *pg_reqgroup);
void adios_transform_free_pg_reqgroup(adios_transform_pg_reqgroup **pg_reqgroup_ptr);

// adios_transform_read_reqgroup manipulation
adios_transform_read_reqgroup * adios_transform_new_read_reqgroup(const ADIOS_FILE *fp, const ADIOS_VARINFO *varinfo,
                                                                  const ADIOS_TRANSINFO *transinfo,
                                                                  const ADIOS_SELECTION *sel, int from_steps, int nsteps,
                                                                  void *data, enum ADIOS_FLAG swap_endianness);
void adios_transform_read_reqgroup_append_pg_reqgroup(adios_transform_read_reqgroup *reqgroup, adios_transform_pg_reqgroup *pg_reqgroup);
int adios_transform_read_reqgroup_remove_pg_reqgroup(adios_transform_read_reqgroup *reqgroup, adios_transform_pg_reqgroup *pg_reqgroup);
adios_transform_pg_reqgroup * adios_transform_read_reqgroup_pop_pg_reqgroup(adios_transform_read_reqgroup *reqgroup);
int adios_transform_read_reqgroup_find_subreq(const adios_transform_read_reqgroup *reqgroup, const ADIOS_VARCHUNK *chunk, int skip_completed,
                                         adios_transform_pg_reqgroup **matching_pg_reqgroup, adios_transform_read_subrequest **matching_subreq);
void adios_transform_read_reqgroup_append(adios_transform_read_reqgroup **head, adios_transform_read_reqgroup *new_reqgroup);
adios_transform_read_reqgroup * adios_transform_read_reqgroup_remove(adios_transform_read_reqgroup **head, adios_transform_read_reqgroup *reqgroup);
adios_transform_read_reqgroup * adios_transform_read_reqgroup_pop(adios_transform_read_reqgroup **head);
void adios_transform_free_read_reqgroup(adios_transform_read_reqgroup **reqgroup_ptr);

#endif /* ADIOS_TRANSFORMS_REQGROUP_H_ */
