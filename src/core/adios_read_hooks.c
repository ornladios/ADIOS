/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include "adios_read_hooks.h"

#define MATCH_STRING_TO_METHOD(b,d,r) \
if (!strcasecmp (buf,b)) \
{*method=d;*requires_group_comm=r;return 1;}

#define ASSIGN_FNS(a,b) \
(*t) [b].adios_init_method_fn = adios_read_##a##_init_method; \
(*t) [b].adios_finalize_method_fn = adios_read_##a##_finalize_method; \
(*t) [b].adios_open_stream_fn = adios_read_##a##_open_stream; \
(*t) [b].adios_open_file_fn = adios_read_##a##_open_file; \
(*t) [b].adios_close_fn = adios_read_##a##_close; \
(*t) [b].adios_advance_step_fn = adios_read_##a##_advance_step; \
(*t) [b].adios_release_step_fn = adios_read_##a##_release_step; \
(*t) [b].adios_inq_var_byid_fn = adios_read_##a##_inq_var_byid; \
(*t) [b].adios_inq_var_stat_fn = adios_read_##a##_inq_var_stat; \
(*t) [b].adios_inq_var_blockinfo_fn = adios_read_##a##_inq_var_blockinfo; \
(*t) [b].adios_schedule_read_byid_fn = adios_read_##a##_schedule_read_byid; \
(*t) [b].adios_perform_reads_fn = adios_read_##a##_perform_reads; \
(*t) [b].adios_check_reads_fn = adios_read_##a##_check_reads; \
(*t) [b].adios_get_attr_byid_fn = adios_read_##a##_get_attr_byid; \
(*t) [b].adios_reset_dimension_order_fn = adios_read_##a##_reset_dimension_order; \
(*t) [b].adios_get_groupinfo_fn = adios_read_##a##_get_groupinfo; \

void adios_read_hooks_init (struct adios_read_hooks_struct ** t)
{
    static int did_init = 0;
    // we need to init only once in the lifetime of an application
    // called from common_read.c/common_read_init_method() and 
    // from common_read.c/common_read_open_*() 
    if (!did_init) {
        *t = (struct adios_read_hooks_struct *)
               calloc (ADIOS_READ_METHOD_COUNT, sizeof (struct adios_read_hooks_struct));

        ASSIGN_FNS(bp,ADIOS_READ_METHOD_BP)
#ifndef __MPI_DUMMY_H__
        ASSIGN_FNS(bp_staged,ADIOS_READ_METHOD_BP_STAGED)
        ASSIGN_FNS(bp_staged1,ADIOS_READ_METHOD_BP_STAGED1)
#endif
#if HAVE_DART
        ASSIGN_FNS(dart,ADIOS_READ_METHOD_DART)
#endif

#if HAVE_DIMES
        ASSIGN_FNS(dimes,ADIOS_READ_METHOD_DIMES)
#endif

#if HAVE_PHDF5
        //ASSIGN_FNS(hdf5,ADIOS_READ_METHOD_HDF5)
#endif

#if HAVE_NSSI
        //ASSIGN_FNS(nssi,ADIOS_READ_METHOD_NSSI)
#endif

#if HAVE_DATATAP
        ASSIGN_FNS(datatap,ADIOS_READ_METHOD_DATATAP)
#endif

        //printf("%s: adios_read_hooks = %x\n",__func__,*t);
        did_init = 1;
    }

}

