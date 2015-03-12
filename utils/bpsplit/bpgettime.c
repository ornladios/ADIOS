/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
   BPGETTIME: Determine number of timesteps in a timestepped BP file

   Author: Norbert Podhorszki, pnorbert@ornl.gov

**/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef __USE_LARGEFILE64
#define __USE_LARGEFILE64
#endif

#include <sys/types.h> 
#include <unistd.h>  
#include <fcntl.h>  // open64
#include <sys/stat.h> // S_IRUSR and alike
#include <errno.h>

#ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#endif

#include <string.h>
#include <getopt.h>

#include "adios_types.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

#ifndef strndup   //HAVE_STRNDUP
#  define strndup(str,len) strdup(str)
#endif

#define DIVIDER "========================================================\n"

#ifndef bool
typedef int bool;
#endif
#define true 1
#define false 0

/** Prototypes */
int bpgettime(char *filein);

/** Global variables */
int verbose   = 0;            // 1: print log to stdout, 2: debug 
bool simple_output  = false;  // true: simple print for scripts to use: min max format


struct option options[] = {
    {"help",        no_argument,       NULL, 'h'},
    {"verbose",     no_argument,       NULL, 'v'},
    {"simple",      no_argument,       NULL, 's'},
    {NULL,          0,                 NULL, 0}
};

static const char *optstring = "+hvs";

char *prgname; /* argv[0] */


void display_help() {
   printf(
"Usage: %s [-h | --help] [-v | --verbose] [-s | --simple] bpfile \n"
"\n"
"Print the number of timesteps available in the BP file\n"
"\n"
"  --help    | -h  Print this help.\n"
"  --verbose | -v  Log activity about what this program is doing.\n"
"                  Extra -v increases the number of messages.\n"
"  --simple  | -s  Simple output for scripts:  min max\n"
"\n"
,prgname
);

}


/** Main */
int main( int argc, char *argv[] ) {
    int excode;
    char *filein  = NULL;

    prgname = strdup(argv[0]);

    /* other variables */
    int c;
    //int last_c='_';
    int idx = 1;
    /* Process the arguments */
    while ((c = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
        switch (c) {
        case 'h':
            display_help();
            return 0;
        case 'v':
            verbose++;
            break;
        case 's':
            simple_output=true;
            break;
        case 1:
            /* This means a field is unknown, could be multiple arg (we do not have such)
               or bad arg*/
            /*
            if (last_c=='z') {
                // process additional argument; 
            }*/
            /* else { Unknown param } */
            fprintf(stderr, "Unrecognized argument: %s\n", optarg);
            return 1;
        default:
            printf ("\nError: Unknown command line option: %s\n\n", argv[optind]);
            display_help();
            exit (1);
            break;
        } /* end switch */
        //last_c = c;
        idx++;
    } /* end while */

    if (argc >= optind+1) {
        // Read the required arguments
        filein  = strndup(argv[optind], 256);
    } else {
        display_help();
        return 1;
    }

    /* Do the work */
    excode = bpgettime(filein);

    return excode; 
}


/** Read the indexes from the input file and parse them into the in_* structs 
 */
int bpgettime(char *filename) {
    struct adios_bp_buffer_struct_v1           * in_bp         = 0;
    struct adios_index_process_group_struct_v1 * in_pg_root    = 0;
    struct adios_index_process_group_struct_v1 * pg = 0;
    uint32_t ts_min=0, ts_max=0, ngroups=0;
    int excode = 0;
    uint32_t version = 0;
    int rc = 0;

    // init bp struct for reading
    in_bp = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    adios_buffer_struct_init (in_bp);

    // open infile 
    rc = adios_posix_open_read_internal (filename, "", in_bp);
    if (!rc) {
        fprintf (stderr, "%s: file not found: %s\n", prgname, filename);
        return 2;
    }

    // read version
    adios_posix_read_version (in_bp);
    adios_parse_version (in_bp, &version);

    // read and parse process group index
    adios_posix_read_index_offsets (in_bp);
    adios_parse_index_offsets_v1 (in_bp);
    adios_posix_read_process_group_index (in_bp);
    adios_parse_process_group_index_v1 (in_bp, &in_pg_root, NULL);

    // close file
    adios_posix_close_internal (in_bp);
    
    // get the smallest timestep (= timestep of first process group)
    // get the largest timestep (= timestep of last process group)
    pg = in_pg_root;
    ts_min = pg->time_index;
    ngroups = 0;
    while (pg) {
        ts_max = pg->time_index;
        ngroups++;
        pg = pg->next; 
    }
    
    if (!simple_output) { 
        printf("Smallest timestep: %u\n", ts_min);
        printf("Largest timestep: %u\n", ts_max);
        if (verbose) {
            printf("Number of groups: %u\n", ngroups);
        }
    } else {
        printf("%u %u\n", ts_min, ts_max); 
    }

    return excode;
}

