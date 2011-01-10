/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
 *  Metadata (index) for DART and DIMES methods
 * 
 *  Data structures to store the index (aka footer for BP format) that can be
 *  put into DataSpaces to describe the variables/attributes there
 */
#ifndef __SPACES_INDEX_H__
#define __SPACES_INDEX_H__

#include <stdint.h>
#include "adios_types.h"
#include "adios_read.h"


/* functions currently defined in adios_dart.c and read_dart.c */

/* Tell the DataSpaces order of dimensions for a 1-3 dim array written from Fortran or C.
   unpack=1: the reverse of packing (to retrieve the original order).
   didx should be an int [3] array in any case.
 */
void ds_dimension_ordering(int ndims, int is_app_fortran, int unpack, int *didx);

void ds_pack_file_info (int time, int nvars, int nattrs, int group_index_len, char * groupname,
                        /*OUT*/char **buf, /*OUT*/int *buf_len);

//ADIOS_FILE * ds_unpack_file_info (char * buf, int buf_len,
//                                  /* OUT */ struct adios_read_dart_data_struct * ds)

void ds_pack_group_info (struct adios_file_struct *fd
        ,struct adios_method_struct * method
        ,struct adios_index_var_struct_v1 *vars_root
        ,struct adios_index_attribute_struct_v1 * attrs_root
        ,char ** buffer, int *buffer_size, int *nvars, int *nattrs
        );

//ADIOS_GROUP * ds_unpack_group_info (char * buf, struct adios_read_dart_group_struct * group);


/***********************
  FILE info buffer: 

SI = sizeof(int) = 4

CURRENT structure
        128 fix bytes, 1 group
        bytes   content
         SI      = length of this buffer (=128 bytes fixed right now)
         SI      = time index
         SI      = number of variables in file (actually in the group)
         SI      = number of attributes in file (actually in the group)
         SI      = length of group index (in GROUP@fn/gn variable)
         SI      = N length of name  
         N      = group name (no \0)
GOAL structure (needs read from space before writing)
        bytes   content
         SI      = length of this buffer (=128 bytes fixed right now)
         SI      = time index
         SI      = number of variables in file 
         SI      = number of attributes in file
         SI      = G number of groups
    G*  ( 
         SI      = length of group index (in GROUP@fn/gn variable)
         SI      = N length of name 
         N      = group name (no \0)
        )
*/


/*********************
  GROUP info buffer:

   SI = sizeof(int) = 4

   Size = N bytes, N is found in FILE buffer
        bytes   content
         SI      = length of this buffer  (=N above)
         SI      = NV number of variables in file (actually in the group)
         SI      = NA number of attributes in file (actually in the group)
    NV*(
         SI      = L length of name  
         L      = variable name (includes path) (no \0)
         SI      = type of variable
         SI      = has time dimension (0 or 1) (time dim is not stored)
         SI      = D number of dimensions
      D*(
         8      = ith dimension
        )
         K      = value of scalar variable (if D = 0), K is according to type
                  string represented as (SI bytes for len, then string w/o \0) 
       )
    NA*(
         SI      = L length of name  
         L      = attribute name (includes path) (no \0)
         SI      = type of attribute
         K      = value of attribute, K is according to type
                  string represented as (SI bytes for len, then string w/o \0) 
       )

*/


#endif
