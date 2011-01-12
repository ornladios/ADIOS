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
(*t) [b].adios_init_fn = adios_read_##a##_init; \
(*t) [b].adios_finalize_fn = adios_read_##a##_finalize; \
(*t) [b].adios_fopen_fn = adios_read_##a##_fopen; \
(*t) [b].adios_fclose_fn = adios_read_##a##_fclose; \
(*t) [b].adios_gopen_fn = adios_read_##a##_gopen; \
(*t) [b].adios_gopen_byid_fn = adios_read_##a##_gopen_byid; \
(*t) [b].adios_gclose_fn = adios_read_##a##_gclose; \
(*t) [b].adios_inq_var_fn = adios_read_##a##_inq_var; \
(*t) [b].adios_inq_var_byid_fn = adios_read_##a##_inq_var_byid; \
(*t) [b].adios_read_var_fn = adios_read_##a##_read_var; \
(*t) [b].adios_read_local_var_fn = adios_read_##a##_read_local_var; \
(*t) [b].adios_read_var_byid_fn = adios_read_##a##_read_var_byid; \
(*t) [b].adios_get_attr_fn = adios_read_##a##_get_attr; \
(*t) [b].adios_get_attr_byid_fn = adios_read_##a##_get_attr_byid;

void adios_read_hooks_init (struct adios_read_hooks_struct ** t)
{
    static int did_init = 0;
    // we need to init only once in the lifetime of an application
    // called from common_read.c/common_read_fopen()
    if (!did_init) {
        *t = (struct adios_read_hooks_struct *)
               calloc (ADIOS_READ_METHOD_COUNT, sizeof (struct adios_read_hooks_struct));

        ASSIGN_FNS(bp,ADIOS_READ_METHOD_BP)
        ASSIGN_FNS(bp_subfile,ADIOS_READ_METHOD_BP_SUBFILE)

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
        ASSIGN_FNS(nssi,ADIOS_READ_METHOD_NSSI)
#endif

#if HAVE_DATATAP
        ASSIGN_FNS(datatap,ADIOS_READ_METHOD_DATATAP)
#endif

        //printf("%s: adios_read_hooks = %x\n",__func__,*t);
        did_init = 1;
    }

}

