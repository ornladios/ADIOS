/* 
 * Staged write of ADIOS files using a staging method
 *
 * Copyright (c) 2008 - 2012.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

#include "utils.h"

void ints_to_str (int n, int *values, char *s)
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

void int64s_to_str (int n, uint64_t *values, char *s)
{
    int i;
    char v[32];
    if (!n) {
        s[0] = '\0';
        return;
    }
    sprintf(s,"%" PRIu64, values[0]);
    for (i=1; i<n; i++)
    {
        sprintf (v,",%" PRIu64, values[i]);
        strcat (s,v);
    }
}


const char * value_to_string (enum ADIOS_DATATYPES type, void * data, int idx)
{
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
            sprintf (s, "%" PRId64, ((int64_t *) data)[idx]);
            break;

        case adios_unsigned_long:
            sprintf (s, "%" PRIu64, ((uint64_t *) data)[idx]);
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


bool file_exists (char * path)
{
    struct stat sb;
    int i = stat ( path, &sb );
    if ( i == 0 )
	/* File found */
	return true;
    return false;
}

bool is_dir(char *path) {
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

