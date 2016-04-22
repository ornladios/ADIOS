/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef ADIOS_BUFFER_H
#define ADIOS_BUFFER_H

#include "public/adios_types.h"
#include "core/adios_internals.h"

/* Set the maximum buffer size usable by one adios_open()...adios_close() operation */
void  adios_databuffer_set_max_size (uint64_t v);

/* Return a size with which the buffer can be extended up to the maximum.
   It returns a default size unless the existing buffer size plus the default size is
   greater than the max size. In that case it returns max size - current size.
   It does not mean that realloc will succeed.
   It does not return current size + something, just the something.
*/
uint64_t adios_databuffer_get_extension_size (struct adios_file_struct *fd);

/* Resize (or create) the buffer to 'size' if size is less than maximum. 
   It does NOT resize the buffer up to the maximum if size is greater than the maximum
*/
int adios_databuffer_resize (struct adios_file_struct *fd, uint64_t size);
void adios_databuffer_free (struct adios_file_struct *fd);



/*
   OBSOLETE functions, cannot remove until disfunctional method functions are revised 
*/
void      adios_buffer_size_requested_set (uint64_t v);
uint64_t  adios_buffer_size_requested_get (void);
void      adios_buffer_size_remaining_set (uint64_t v);
void      adios_buffer_alloc_percentage_set (int v);
void      adios_buffer_alloc_when_set (enum ADIOS_BUFFER_ALLOC_WHEN v);

enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when_get (void);

int adios_set_buffer_size (void);
uint64_t adios_method_buffer_alloc (uint64_t size);
int adios_method_buffer_free (uint64_t size);

#endif
