/* 
 * Staged write of ADIOS files using a staging method
 *
 * Copyright (c) 2008 - 2012.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef __UTILS_H_
#define __UTILS_H_

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "adios_types.h"

#define MAX3(a,b,c) (a > b ? (a > c ? a : c) : (b > c ? b : c))
#define print(...) fprintf (stderr, __VA_ARGS__); 
#define print0(...) if (!rank) fprintf (stderr, __VA_ARGS__); 

#define bool int
#define false 0
#define true 1

void ints_to_str (int n, int *values, char *s);
void int64s_to_str (int n, uint64_t *values, char *s);
const char * value_to_string (enum ADIOS_DATATYPES type, void * data, int idx);

/* Get basename and dirname from path.
   Allocates memory for both strings, they should be freed after use.
   May return NULL for dirname. "" for dirname means path is like /file.

   dir1/dir2/file   dir=dir1/dir2   base=file
   dir1/file        dir=dir1        base=file
   file             dir=null        base=file
   /file            dir=""          base=file
   dir1/            dir=null        base=dir1
   /                dir=""          base=""
   ""               dir=null        base=""
*/
void getbasename (char *path, char **dirname, char **basename);

bool file_exists (char *path);  // true if stat(path) succeeds, 
                                // i.e. it is an accessible item on file system
bool is_dir (char *path);       // true if path is a directory (and stat() succeeds on it)

/** mkdir -r
  * return: 0 on success, otherwise mkdir() syscall's return value
  */
int createdir_recursive( char* path);

#endif

