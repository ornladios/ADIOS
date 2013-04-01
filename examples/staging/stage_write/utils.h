/* 
 * Staged write of ADIOS files using a staging method
 *
 * Copyright (c) 2008 - 2012.  UT-BATTELLE, LLC. All rights reserved.
 */

#ifndef __UTILS_H_
#define __UTILS_H_

#include <stdio.h>
#include <stdint.h>
#include "adios_types.h"

#define MAX3(a,b,c) (a > b ? (a > c ? a : c) : (b > c ? b : c))
#define print(...) fprintf (stderr, __VA_ARGS__); 
#define print0(...) if (!rank) fprintf (stderr, __VA_ARGS__); 

#define bool int
#define false 0
#define true 1

static inline void ints_to_str (int n, int *values, char *s)
{
    int i;
    char v[32];
    if (!n) {
        s[0] = '\0';
        return;
    }
    sprintf(s,"%d", values[0]);
    for (i=1; i<n; i++)
    {
        sprintf (v,",%d", values[i]);
        strcat (s,v);
    }
}
static inline void int64s_to_str (int n, uint64_t *values, char *s)
{
    int i;
    char v[32];
    if (!n) {
        s[0] = '\0';
        return;
    }
    sprintf(s,"%llu", values[0]);
    for (i=1; i<n; i++)
    {
        sprintf (v,",%llu", values[i]);
        strcat (s,v);
    }
}

static inline const char * value_to_string (enum ADIOS_DATATYPES type, void * data, int idx)
    static char s [100];
    s [0] = 0;


    switch (type)
    {
        case adios_unsigned_byte:
            sprintf (s, "%u", ((uint8_t *) data)[idx]);
            break;

        case adios_byte:
            sprintf (s, "%d", ((int8_t *) data)[idx]);
            break;

        case adios_short:
            sprintf (s, "%hd", ((int16_t *) data)[idx]);
            break;

        case adios_unsigned_short:
            sprintf (s, "%hu", ((uint16_t *) data)[idx]);
            break;

        case adios_integer:
            sprintf (s, "%d", ((int32_t *) data)[idx]);
            break;

        case adios_unsigned_integer:
            sprintf (s, "%u", ((uint32_t *) data)[idx]);
            break;

        case adios_long:
            sprintf (s, "%lld", ((int64_t *) data)[idx]);
            break;

        case adios_unsigned_long:
            sprintf (s, "%llu", ((uint64_t *) data)[idx]);
            break;

        case adios_real:
            sprintf (s, "%g", ((float *) data)[idx]);
            break;

        case adios_double:
            sprintf (s, "%lg", ((double *) data)[idx]);
            break;

        case adios_long_double:
            sprintf (s, "%Lg", ((long double *) data)[idx]);
            break;

        case adios_string:
            return (char*) ((char *)data+idx);
            break;

        case adios_complex:
            sprintf (s, "(%g, %g)", 
                    ((float *) data)[2*idx], ((float *) data)[2*idx+1]);
            break;

        case adios_double_complex:
            sprintf (s, "(%lg, %lg)", 
                    ((double *) data)[2*idx], ((double *) data)[2*idx+1]);
            break;
        
        default:
            sprintf (s, "unknown");
            break;
    }

    return s;
}

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
void getbasename (char *path, char **dirname, char **basename)
{
    char *work, *ptr;

    work = strdup(path);
    if (work[strlen(work)-1] == '/' && strlen(work)>1)
        work[strlen(work)-1] = '\0';

    ptr = rindex(work, '/');
    if (ptr && ptr != work) {
        // found a / and but not the beginning / of a full path
        ptr[0] = 0;
        *dirname = strdup(work);
        ptr[0] = '/';
        *basename = strdup(ptr+1);
    } else if (ptr == work) {
        // found / as the first character 
        *dirname = strdup("");
        *basename = strdup(ptr+1);
    } else {
        *dirname = NULL; //strdup(".");
        *basename = strdup(work);
    }
    free(work);
}

bool file_exists (char *path)  
{
// true if stat(path) succeeds, 
// i.e. it is an accessible item on file system
    struct stat sb;
    int i = stat ( path, &sb );
    if ( i == 0 )
	/* File found */
	return true;
    return false;
}
bool is_dir (char *path)
{
       // true if path is a directory (and stat() succeeds on it)
     struct stat sb;
     if ( !stat( path, &sb) ) {
         if ( sb.st_mode & S_IFDIR ) {
             //if (verbose) printf("isDir( %s ) = true\n", path);
             return  true;
         }
     } /* else {
         // for whatever reason, stat cannot be retrieved
         // so just say here, it is not a directory 
        retval = false;
     }*/
     return false;
}

/** mkdir -r
  * return: 0 on success, otherwise mkdir() syscall's return value
  */
int createdir_recursive( char* path)
{
    int res;
    char *dirs, *base;
    /*printf(" called mkdir %s\n", path);*/
    getbasename(path, &dirs, &base);
    if ( dirs != NULL)  {
        if ( !is_dir(dirs) ) {
	    /*printf(" dirs=%s  base=%s\n", dirs, base);*/
            res = createdir_recursive(dirs);
	} else res = 0;
        free(dirs);
        if (res) return res;
    }
    if (file_exists(path))
        res = 0;
    else
        res = mkdir(path, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
    if (res) {
        if (errno == EEXIST) {
            print("Concurrency error: createdir %s failed: %s\n", path, strerror(errno));
            print("Some other process created this dir already during this function call\n");
            res = 0;
        } else {
            print("createdir %s failed: %s\n", path, strerror(errno));
        }
    }
    //else
    //   print("createdir %s succeeded\n", path);

    free(base);
    return (res);
}


#endif

