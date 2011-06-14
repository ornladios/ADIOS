#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

#include "globals.h"  // bool, false, true
#include "dirutil.h"


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
            fprintf(stderr, "Concurrency error: createdir %s failed: %s\n", path, strerror(errno));
            fprintf(stderr, "Some other process created this dir already during this function call\n");
            res = 0;
        } else {
            fprintf(stderr, "createdir %s failed: %s\n", path, strerror(errno));
        }
    }
    //else
    //   printf("createdir %s succeeded\n", path);

    free(base);
    return (res);
}

