#ifndef __DIRUTUL_H_
#define __DIRUTUL_H_

#include "globals.h"

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
