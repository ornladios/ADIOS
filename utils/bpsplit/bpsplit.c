/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
   BPSPLIT: Carve out some timesteps from a timestepped BP file

   Author: Norbert Podhorszki, pnorbert@ornl.gov

**/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

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

#if defined(__APPLE__) || defined(__WIN32__)  || defined(__CYGWIN__) 
#    define lseek64 lseek
#    define open64  open
#endif

#define DIVIDER "========================================================\n"

#ifndef bool
typedef int bool;
#endif
#define true 1
#define false 0

/** Prototypes */
int bpsplit(char *filein, char *fileout, char *recordfile, int from_in, int to_in, bool skiplast);
void cleanup(void);

/** Global variables */
int verbose   = 0;            // 1: print log to stdout, 2: debug 


struct option options[] = {
    {"from",        required_argument, NULL, 'n'},
    {"to",          required_argument, NULL, 'm'},
    {"recordfile",  required_argument, NULL, 'r'},
    {"skiplast",    no_argument,       NULL, 's'},
    {"help",        no_argument,       NULL, 'h'},
    {"verbose",     no_argument,       NULL, 'v'},
    {NULL,          0,                 NULL, 0}
};

static const char *optstring = "+hvn:m:r:s";

char *prgname; /* argv[0] */


void display_help() {
   printf(
"Usage: %s [-h | --help] [-v | --verbose] [--from N | -n N] [--to M | -m M]\n"
"          [--recordfile path | -r path ] [ --skiplast | -s ] \n"
"          inputfile outputfile\n"
"\n"
"Copy some timesteps from a time-stepped BP file into another file.\n"
"\n"
"  --help | -h     Print this help.\n"
"  --verbose | -v  Log activity about what this program is doing.\n"
"                  Extra -v increases the number of messages.\n"
"  --from | -n N   Start from Nth timestep. \n"
"                  Note: BP file timesteps starts from 1.\n"
"  --to | -m M     Finish at Mth timestep.\n"
"                  -k means last+k+1 (-1=last)\n"
"  --recordfile | \n"
"   -r path        write out last timestep into file <path>\n"
"  --skiplast      Skip the last timestep (like -m -2)\n"
"\n"
"  Default behavior: -n -1 -m -1 \n"
"\n"
"  Do not use -r and -n together. If -r is used, bpsplit will examine the\n" 
"  content of file 'last' and splits all timesteps from inputfile that\n"
"  are newer (greater) then that timestep\n"
"\n"
"  --skiplast would do the same up to the one before last timestep.\n"
"  This is useful if you expect a process may be still writing in.bp\n"
,prgname
);

}


/** Main */
int main( int argc, char *argv[] ) {
    int excode;
    int from      = -1;           // split from  'from'. -1=last
    int to        = -1;           // split until 'to'. -1=last
    bool skiplast = false;
    char *filein  = NULL;
    char *fileout = NULL;       
    char *recordfile = NULL;      // contains last timestep splitted before
    bool n_was_set = false;       // test if both n and recordfile is specified 

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
        case 'n':
            errno = 0; 
            from = strtol(optarg, (char **)NULL, 0);
            if (errno) {
                fprintf(stderr, "Error: could not convert --from/-n option's value: %s\n", optarg);
                return 1;
            }
            n_was_set=true;
            break;
        case 'm':
            errno = 0; 
            to = strtol(optarg, (char **)NULL, 0);
            if (errno) {
                fprintf(stderr, "Error: could not convert --to/-m option's value: %s\n", optarg);
                return 1;
            }
            break;
        case 'r':
            recordfile = strndup(optarg,256);
            break;
        case 's':
            skiplast = true;
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
            printf ("\nError: Unknown command line option after: %s\n\n", argv[optind]);
            display_help();
            exit (1);
            break;
        } /* end switch */
        //last_c = c;
        idx++;
    } /* end while */

    if (argc >= optind+2) {
        // Read the required arguments
        filein  = strndup(argv[optind], 256);
        fileout = strndup(argv[optind+1], 256);
    } else {
        display_help();
        return 1;
    }

    if (n_was_set && recordfile) {
        fprintf(stderr, "Error: -n/--from index and -r/--recordfile should not be used together!\n");
        return 2;
    }

    /* Do the work */
    excode = bpsplit(filein, fileout, recordfile, from, to, skiplast);
    if (excode)
        cleanup();

    return excode; 
}

// global vars for the functions below
//   filein's data structs
struct adios_bp_buffer_struct_v1           * in_bp         = 0;
struct adios_index_struct_v1               * idx           = 0;
struct adios_index_process_group_struct_v1 * in_pg_root    = 0;
//   fileout's data structs 
struct adios_bp_buffer_struct_v1           * out_bp         = 0;
//struct adios_index_process_group_struct_v1 * out_pg_root    = 0; // use idx->pg_root
uint64_t out_offset_start = 0;  // the beginning offset of group data in in_bp to write out to out_bp
uint64_t out_offset_end = 0;    // the end offset of group data in in_bp to write out to out_bp
//   aux variables to contain tail of index chains (for cleanup only)
struct adios_index_process_group_struct_v1 * tail_pg_root    = 0;
//   recordfile
uint32_t lasttime; // timestep written previously into the recordfile

/** Read the indexes from the input file and parse them into the in_* structs 
 */
int read_indexes(char *filename) {
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
    
    // read and parse variable index
    adios_posix_read_vars_index (in_bp);
    adios_parse_vars_index_v1 (in_bp, &idx->vars_root, NULL, NULL);

    // read and parse attribute index
    adios_posix_read_attributes_index (in_bp);
    adios_parse_attributes_index_v1 (in_bp, &idx->attrs_root);

    if (verbose>1) {
        printf (DIVIDER);
        printf ("Process Groups Index:\n");
        printf ("End of process groups       = %" PRIu64 "\n", in_bp->end_of_pgs);
        printf ("Process Groups Index Offset = %" PRIu64 "\n", in_bp->pg_index_offset);
        printf ("Process Groups Index Size   = %" PRIu64 "\n", in_bp->pg_size);
        printf ("Variable Index Offset       = %" PRIu64 "\n", in_bp->vars_index_offset);
        printf ("Variable Index Size         = %" PRIu64 "\n", in_bp->vars_size);
        printf ("Attribute Index Offset      = %" PRIu64 "\n", in_bp->attrs_index_offset);
        printf ("Attribute Index Size        = %" PRIu64 "\n", in_bp->attrs_size);
    }

    return excode;
}

/** get the largest timestep in the input file (= timestep of last process group) 
 */
uint32_t timesteps_available() {
    struct adios_index_process_group_struct_v1 * pg = in_pg_root;
    uint32_t ts;
    while (pg) {
        ts = pg->time_index;
        pg = pg->next; 
    }
    if (verbose) printf("Largest timestep in input file: %u\n", ts);
    return ts;
}

/** get last timestep written into the recordfile (if exist) 
 *  lasttime is set to the value in the file or 0
 *  return value is error code, 
 *  0 on success, 
 *  -1 if file does not exists (not an error!) 
 *  1 for failure to read existing file
 */
int get_last_time(char *recordfile) {
    int excode = 0;
    int n;
    FILE * rf;
    lasttime = 0;  // there is at least 1 group with Time=1 in a bp file
    if (recordfile == NULL) 
        return 0;
    if (verbose>1) printf("Try to read recordfile: %s\n", recordfile);
    rf = fopen( recordfile, "r");
    if (rf == NULL) 
        return -1; // no problem, just no such file (yet)
    if (verbose) printf("Read recordfile: %s\n", recordfile);
    n = fscanf(rf, "%u\n", &lasttime);
    if (n != 1) {
        fprintf(stderr, "Recordfile %s does not contain an unsigned integer value in the first line: %s\n",
           recordfile, strerror(errno));
        excode = 1;
        lasttime = 0;
    }
    if (verbose && n==1) printf("Last time from recordfile: %u\n", lasttime);
    fclose(rf);

    return excode;
}

/** Write out last timestep into the recordfile.
 */
int record_last_time(char *recordfile, int laststep) {
    FILE * rf;
    if (recordfile == NULL) 
        return 0;
    if (verbose>1) printf("Try to write last step into recordfile: %s\n", recordfile);
    rf = fopen( recordfile, "w");
    if (rf == NULL) {
        fprintf(stderr, "Could not create recordfile %s: %s\n", recordfile, strerror(errno));
        return 1;
    }
    if (verbose) printf("Write step %d into recordfile: %s\n", laststep, recordfile);
    fprintf( rf, "%u\n", laststep);
    fclose(rf);
    return 0;
}

/** Handle -maxt+1 ... -1 and 1 ... maxt indexes.
    Return 0 on any other (invalid) indexes
 */
uint32_t handle_negative_index( int t, uint32_t maxt, char *var) {
    int j = t;
    if (j < 0) {
        j = maxt + j + 1; // -1 = last, -2 = last-1...
    }
    if (j <= 0 || j > maxt) {
        fprintf(stderr, "Error: %s index %u should be in 1..%u or -%u..-1\n",
            var, t, maxt, maxt-1);
        j = 0;
    }
    return j;
}

/** 1. Split in_pg_root index chain to 3 parts (1..from-1, from..to, to+1..end), so that 
 *    in_pg_root   index chain contains indexes from the beginning
 *    out_pg_root  index chain contains indexes from 'from' to 'to' timesteps
 *    tail_pg_root index chain contains anything greater than 'to' timestep
 *  We assume here that 
 *     1 <= from <= to 
 *     the chain is ordered by time
 *
 *  2. Also recalculate the offsets of groups to be written out to their future
 *  position in the output.
 *
 *  3. Set the start and end offsets in the input file that contains all groups
 *  to be copied in out_offset_start and out_offset_end
 *      out_offset_start   points to the first byte to be copied.
 *      out_offset_end     points to a byte which is not copied!
 */
void split_pg_index( uint32_t from, uint32_t to) {
    struct adios_index_process_group_struct_v1 * pg = in_pg_root;
    struct adios_index_process_group_struct_v1 * pg_prev = NULL;
    int    section=0; // 0=in, 1=out, 2=tail
    out_offset_start = 0; // default to return if there is nothing to write out too
    out_offset_end = 0;
    if (verbose) printf("Split process group index chain\n");
    while (pg) {
        if (section == 0 && pg->time_index >= from) {
            // reached from..to section
            // start out_pg_root index chain
            idx->pg_root = pg;
            // this is the starting offset from which data should be copied
            out_offset_start = pg->offset_in_file;
            // unlink previous->next pointer to this item
            if (pg_prev) pg_prev->next = NULL;
            else in_pg_root = NULL;
            section = 1;
            if (verbose>1) printf("  section out starts at time %d\n", pg->time_index);
        }
        else if (section < 2 && pg->time_index > to) {
            // reached to+1.. section
            // start tail_pg_root index chain
            tail_pg_root = pg;
            // unlink previous->next pointer to this item
            if (pg_prev) pg_prev->next = NULL;
            else in_pg_root = NULL;
            if (verbose>1) printf("  section tail starts at time %d\n", pg->time_index);
            break; // nothing more to look for
        }

        if (section == 1)  {
            // recalculate the offset of the outgoing group to the position in the output
            if (verbose>1) printf("    group time %d offset %" PRId64 " -> %" PRId64 "\n",
                            pg->time_index, pg->offset_in_file, pg->offset_in_file - out_offset_start);
            pg->offset_in_file -= out_offset_start;
        }

        pg_prev = pg;
        pg = pg->next; 
    }

    // determine the out_offset
    if (tail_pg_root) 
        out_offset_end = tail_pg_root->offset_in_file; // end points to a byte which is not copied!
    else
        out_offset_end = in_bp->pg_index_offset; // pg index is right after all the pg data in the bp file
    if (verbose>1) printf("  offset start = %" PRIu64 "  end = %" PRIu64 "\n", out_offset_start, out_offset_end);
}

/** Determine the start and beginning offsets in input file (of groups) that should be
 *  copied into output file.
 *  out_offset_start points to the first byte to be copied.
 *  out_offset_end points to a byte which is not copied!
 */
/*
void determine_pg_offsets() {
    // determine offsets
    if (verbose) printf("Determine process group offsets\n");
    if (idx->pg_root) {
        out_offset_start = idx->pg_root->offset_in_file;
        if (tail_pg_root) 
            out_offset_end = tail_pg_root->offset_in_file; // end points to a byte which is not copied!
        else
            out_offset_end = in_bp->pg_index_offset; // pg index is right after all the pg data in the bp file
    } else {
        out_offset_start = 0;
        out_offset_end = 0;
    }
    if (verbose>1) printf("  offset start = %" PRIu64 "  end = %" PRIu64 "\n", out_offset_start, out_offset_end);
}
*/

/** Remove all characteristics from variable and attr indexes that fall outside 
 *  the offsets to be written out. Also recalculate the offsets of characteristics
 *  to their future position in the output. 
 *  Call after split_pg_index()
 */
void weed_out_indexes(void) {
    struct adios_index_var_struct_v1       * vg = idx->vars_root;
    struct adios_index_var_struct_v1       * vg_prev = NULL;
    struct adios_index_attribute_struct_v1 * ag = idx->attrs_root;
    struct adios_index_attribute_struct_v1 * ag_prev = NULL;
    int i, start, count;
    
    // process variables
    if (verbose) printf("Weed out characteristics from variables\n");
    while (vg) {
        // look for those characteristics that are within our offsets
        start = 0;
        count = 0; 
        for (i=0; i < vg->characteristics_count; i++) {
            if (vg->characteristics[i].offset < out_offset_start)
                start++;
            else if (vg->characteristics[i].offset >= out_offset_start && 
                     vg->characteristics[i].offset <= out_offset_end) {
                // recalc offset for output offset
                vg->characteristics[i].offset -= out_offset_start; 
                count++;
            } else
                break;
        }
        if (verbose>1)
            printf("    Var %s/%s: characteristics start=%d, count=%d to be kept\n", 
                    vg->var_path, vg->var_name, start, count);
        if (start > 0) {
            for (i=0; i<count; i++)   // shift valid characteristics to the beginning of the array
                vg->characteristics[i] = vg->characteristics[start+i];
        }
        vg->characteristics_count = count;
        if (count == 0) {
            // no characteristics <=> this variable is not contained in the output slice
            // take it out from the chain
            if (vg_prev != NULL) vg_prev->next = vg->next;
            else idx->vars_root = vg->next;
        } else {
            vg_prev = vg; // advance prev only if this variable is kept in chain
        }
        vg = vg->next;
    }

    // process attributes
    if (verbose) printf("Weed out characteristics from attributes\n");
    while (ag) {
        // look for those characteristics that are within our offsets
        start = 0;
        count = 0; 
        for (i=0; i < ag->characteristics_count; i++) {
            if (ag->characteristics[i].offset < out_offset_start)
                start++;
            else if (ag->characteristics[i].offset >= out_offset_start && 
                     ag->characteristics[i].offset <= out_offset_end) {
                // recalc offset for output offset
                ag->characteristics[i].offset -= out_offset_start; 
                count++;
            } else
                break;
        }
        if (verbose>1)
            printf("    Attr %s/%s: characteristics start=%d, count=%d to be kept\n", 
                ag->attr_path, ag->attr_name, start, count);
        if (start > 0) {
            for (i=0; i<count; i++) 
                ag->characteristics[i] = ag->characteristics[start+i];
        }
        ag->characteristics_count = count;
        if (count == 0) {
            // no characteristics <=> this attribute is not contained in the output slice
            // take it out from the chain
            if (ag_prev != NULL) ag_prev->next = ag->next;
            else idx->attrs_root = ag->next;
        } else {
            ag_prev = ag; // advance prev only if this attribute is kept in chain
        }
        ag = ag->next;
    }
}

/** Write the output file.
 *   copy data from input file from min offset to max offset
 *   write out indexes
 *   write out version
 *
 *  FIXME: how should we copy if endianness should be changed?
 */
#define COPYBUFSIZE 1024
int write_out( const char *fileout, const char *filein) {
    int f;
    // open file
    if (verbose) printf("Write data to output file %s\n", fileout);
    errno = 0;
    f = open64( fileout, O_WRONLY | O_CREAT | O_TRUNC,
                S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH
              );
    if (f == -1) {
        fprintf(stderr, "Cannot open output file %s: %s\n", fileout, strerror(errno));
        return 1;
    }
   
    // copy data
    char buf[COPYBUFSIZE];
    ssize_t bytes_read, bytes_written;
    uint64_t bytes_to_copy = out_offset_end - out_offset_start; // end byte should not be copied
    uint64_t bytes_copied  = 0;
    if (verbose>1) printf("  seek in input file to %" PRIu64 "\n", out_offset_start);
    lseek64 (in_bp->f, out_offset_start, SEEK_SET);
    while (bytes_copied < bytes_to_copy) {
        bytes_read = read( in_bp->f, (void *)buf, COPYBUFSIZE);
        if (bytes_read < 0) {
            fprintf(stderr, "Error at reading input file %s from offset %" PRIu64 ": %s\n",
                    filein, out_offset_start+bytes_copied, strerror(errno));
            close(f);
            return 2;
        } else if (bytes_read == 0) {
            fprintf(stderr, "Error unexpected end of input file %s at offset %" PRIu64 ": %s\n",
                    filein, out_offset_start+bytes_copied, strerror(errno));
            close(f);
            return 3;
        }
        // check to not to copy too much
        if (bytes_copied + bytes_read > bytes_to_copy)
            bytes_read = bytes_to_copy - bytes_copied;

        bytes_written = write( f, (void *) buf, bytes_read);
        if (bytes_written != bytes_read) {
            fprintf(stderr, "Error: could not write %lld bytes to output file %s at offset %" PRIu64 ": %s\n",
                    (long long)bytes_read, fileout, out_offset_start+bytes_copied, strerror(errno));
            close(f);
            return 4;
        }

        bytes_copied += bytes_written;
    }
    if (verbose>1) printf("  written %" PRIu64 " %" PRIx64 " bytes of data into %s\n", bytes_copied, bytes_copied, fileout);

    // write indexes and version into a buffer
    char * buffer = NULL;
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;
    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset, bytes_copied, idx);
    if (verbose>1) printf("  index size %" PRIu64 " 0x%" PRIx64 "\n", buffer_offset, buffer_offset);
    adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

    // write buffer out
    if (verbose>1) printf("  write %" PRIu64 " 0x%" PRIx64 " bytes of indexes into %s\n", buffer_offset, buffer_offset, fileout);
    bytes_written = write (f, buffer, buffer_offset);

    if (verbose>1) printf("  written %zu 0x%zx bytes of indexes into %s\n", bytes_written, bytes_written, fileout);

    // clean up
    free(buffer);
    close(f);

    return 0;
}

int bpsplit(char *filein, char *fileout, char *recordfile, int from_in, int to_in, bool skiplast) {
    int excode = 0;
    uint32_t from, to;
    uint32_t maxtime = 1;   // at least there is time=1 (single group) in a bp file
    idx = adios_alloc_index_v1(0);

    // open input file, read and parse indexes 
    excode = read_indexes( filein );
    if (excode) 
        return excode; // finish here on error
    maxtime = timesteps_available();

    // get last timestep from recordfile (if exists)
    excode = get_last_time(recordfile);
    if (excode > 0) // -1 is not an error
        return excode; // finish here on error
    if (excode == -1)
        from = 1;  // nothing in recordfile, so we start from beginning
    if (lasttime >= maxtime) {
        fprintf(stderr, "Error: last timestep recorded in recordfile %s (%u) >= maxtime %u. No output file is created\n",
            recordfile, lasttime, maxtime);
        return 5; // finish here on error
    }
    // if recordfile was given, it's content is stronger than from, to!
    if (lasttime > 0) {
        from = lasttime+1;
    }


    // set from and to correctly now (handle negative values)
    //  type is converted from int to uint32_t at this point
    if (!recordfile) 
        from = handle_negative_index(from_in, maxtime, "from");
    to  = handle_negative_index(to_in, maxtime, "to");
    // handle skiplast option
    if (skiplast && to == maxtime) {
        if (verbose) printf("Skip-last option applied.\n");
        to--;
    }
    if (verbose) printf("Indexes from %u to %u\n", from, to);
    if (from == 0 || to == 0)
        return 4; // finish here on error
    if (from > to) {
        fprintf(stderr, "Error: From index evaluated %d > to index evaluated %d.\n", from, to);
        return 4; // finish here on error
    }

    if (lasttime >= to) {
        fprintf(stderr, "Error: To index %u <=  what was recorded in recordfile %s (%u).\n",
            to, recordfile, lasttime);
        return 6; // exit with error, because no output file will be created.
    }

    // split in_pg_root index chain to in*, out* and tail*  (1..from-1, from..to, to+1..end)
    // also get min and max offset we have to copy
    split_pg_index( from, to);
    if (out_offset_start == out_offset_end) {
        fprintf(stderr, "Error: Nothing to write out. File is not created.\n");
        return 7;
    }
    weed_out_indexes();     // remove things from vars and attrs referring outside the offsets
    
    // write output file
    excode = write_out(fileout, filein);  // filein just for error messages
    if (excode > 0) 
        return excode; // finish here on error

    // write recordfile (if exists)
    excode = record_last_time( recordfile, to);

    return excode;
}

void cleanup(void) {
    if (!in_bp)
         adios_posix_close_internal (in_bp);

}
