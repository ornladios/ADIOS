/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_INTERNALS_MXML_H
#define ADIOS_INTERNALS_MXML_H

int adios_parse_config (const char * config);
int adios_local_config ();
int adios_common_select_method (int priority, const char * method
                               ,const char * parameters, const char * group 
                               ,const char * base_path, int iters
                               );
int adios_common_select_method_by_group_id (int priority, const char * method
                                           ,const char * parameters, int64_t group_id
                                           ,const char * base_path, int iters
                                           );
void adios_cleanup ();

int adios_set_buffer_size (void);

uint64_t adios_method_buffer_alloc (uint64_t size);
int adios_method_buffer_free (uint64_t size);

#endif
