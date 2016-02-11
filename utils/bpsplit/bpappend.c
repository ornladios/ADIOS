/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/*
   BPAPPEND: append a BP file to another 

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
int bpappend(char *filein, char *fileout);
void cleanup(void);

/** Global variables */
int verbose   = 0;            // 1: print log to stdout, 2: debug 


struct option options[] = {
    {"help",        no_argument,       NULL, 'h'},
    {"verbose",     no_argument,       NULL, 'v'},
    {NULL,          0,                 NULL, 0}
};

static const char *optstring = "+hv";

char *prgname; /* argv[0] */


void display_help() {
   printf(
"Usage: %s [-h | --help] [-v | --verbose] splitfile appendfile\n"
"\n"
"Append timesteps into a time-stepped BP file from another file BP file splitted by bpsplit.\n"
"splitfile will be appended to appendfile.\n"
"\n"
"  --help | -h     Print this help.\n"
"  --verbose | -v  Log activity about what this program is doing.\n"
"                  Extra -v increases the number of messages.\n"
"\n"
,prgname
);

}


/** Main */
int main( int argc, char *argv[] ) {
    int excode;
    char *filein  = NULL;
    char *fileout = NULL;       

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

    if (argc >= optind+2) {
        // Read the required arguments
        filein  = strndup(argv[optind], 256);
        fileout = strndup(argv[optind+1], 256);
    } else {
        display_help();
        return 1;
    }

    /* Do the work */
    excode = bpappend(filein, fileout);
    if (excode)
        cleanup();

    return excode; 
}

// global vars for the functions below
//   filein's data structs
struct adios_bp_buffer_struct_v1           * in_bp         = 0;
struct adios_index_process_group_struct_v1 * in_pg_root    = 0;
struct adios_index_var_struct_v1           * in_vars_root  = 0; 
struct adios_index_attribute_struct_v1     * in_attrs_root = 0; 
//   fileout's data structs 
struct adios_bp_buffer_struct_v1           * out_bp         = 0;
struct adios_index_process_group_struct_v1 * out_pg_root    = 0;
struct adios_index_var_struct_v1           * out_vars_root  = 0; 
struct adios_index_attribute_struct_v1     * out_attrs_root = 0; 
int    outf;  // file handle for output file (we write raw since we copy groups, not vars and attrs)

uint64_t append_offset_start = 0;  // the beginning offset in out_bp to write group data of in_bp

//   aux variables to contain tail of index chains (for cleanup only)
struct adios_index_process_group_struct_v1 * tail_pg_root    = 0;

/** Read the indexes from the input file and parse them into the in_* structs 
 */
int read_indexes(char *filename, bool input) {
    struct adios_bp_buffer_struct_v1           * bp;
    struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_var_struct_v1           * vars_root = 0;
    struct adios_index_attribute_struct_v1     * attrs_root = 0;
    int excode = 0;
    uint32_t version = 0;
    int rc = 0;

    if (verbose) printf("Read indexes from %s\n", filename);

    // init bp struct for reading
    bp = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    adios_buffer_struct_init (bp);

    // open infile 
    rc = adios_posix_open_read_internal (filename, "", bp);
    if (!rc) {
        if (input) fprintf (stderr, "%s: file not found: %s\n", prgname, filename);
        return 128;
    }

    // read version
    adios_posix_read_version (bp);
    adios_parse_version (bp, &version);

    // read and parse process group index
    adios_posix_read_index_offsets (bp);
    adios_parse_index_offsets_v1 (bp);
    adios_posix_read_process_group_index (bp);
    adios_parse_process_group_index_v1 (bp, &pg_root, NULL);
    
    // read and parse variable index
    adios_posix_read_vars_index (bp);
    adios_parse_vars_index_v1 (bp, &vars_root, NULL, NULL);

    // read and parse attribute index
    adios_posix_read_attributes_index (bp);
    adios_parse_attributes_index_v1 (bp, &attrs_root);

    if (verbose>1) {
        printf (DIVIDER);
        printf ("%s file %s:\n", (input ? "Input" : "Output"), filename);
        printf ("Process Groups Index:\n");
        printf ("Process Groups Index Offset = %" PRIu64 "\n", bp->pg_index_offset);
        printf ("Process Groups Index Size   = %" PRIu64 "\n", bp->pg_size);
        printf ("Variable Index Offset       = %" PRIu64 "\n", bp->vars_index_offset);
        printf ("Variable Index Size         = %" PRIu64 "\n", bp->vars_size);
        printf ("Attribute Index Offset      = %" PRIu64 "\n", bp->attrs_index_offset);
        printf ("Attribute Index Size        = %" PRIu64 "\n", bp->attrs_size);
    }

    // set global pointers 
    if (input) {
        in_bp         = bp;
        in_pg_root    = pg_root;
        in_vars_root  = vars_root;
        in_attrs_root = attrs_root;
    } else {
        out_bp         = bp;
        out_pg_root    = pg_root;
        out_vars_root  = vars_root;
        out_attrs_root = attrs_root;
        // we do not need the readonly file handle anymore
        close(bp->f);
        bp->f = -1;
    }

    return excode;
}

/** Copy the input file to the output file.
 *  FIXME: how should we copy if endianness should be changed?
 */
#define COPYBUFSIZE 1024
int copy_file( const char *filein, const char *fileout) {
    int inf, outf;
    if (verbose) printf("Copy input %s to output %s\n", filein, fileout);

    // open files
    errno = 0;
    inf = open64( filein, O_RDONLY);  
    if (inf == -1) {
        fprintf(stderr, "Cannot open input file %s: %s\n", filein, strerror(errno));
        return 1;
    }
    errno = 0;
    outf = open64( fileout, O_WRONLY | O_CREAT | O_TRUNC,
                S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH
              );
    if (outf == -1) {
        fprintf(stderr, "Cannot open output file %s: %s\n", fileout, strerror(errno));
        close(inf);
        return 1;
    }
   
    // copy data
    char buf[COPYBUFSIZE];
    ssize_t bytes_read, bytes_written;
    uint64_t bytes_copied  = 0;
    bytes_read = read( inf, (void *)buf, COPYBUFSIZE);
    while (bytes_read > 0) {
        bytes_written = write( outf, (void *) buf, bytes_read);
        if (bytes_written != bytes_read) {
            fprintf(stderr, "Error: could not write %lld bytes to output file %s at offset 0: %s\n",
                    (long long)bytes_read, fileout, strerror(errno));
            close(inf);
            close(outf);
            return 4;
        }
        bytes_copied += bytes_written;
        bytes_read = read( inf, (void *)buf, COPYBUFSIZE);
    }
    if (bytes_read < 0) {
        fprintf(stderr, "Error at reading input file %s from offset %" PRIu64 ": %s\n",
                filein, bytes_copied, strerror(errno));
        close(inf);
        close(outf);
        return 2;
    }
    if (verbose>1) printf("  copied %" PRIu64 " (0x%" PRIx64 ") bytes of data into %s\n", bytes_copied, bytes_copied, fileout);

    close(inf);
    close(outf);

    return 0;
}


/** Return error (-1) if timesteps in input is NOT >= than timesteps in output.
 */
int check_timesteps(void) {
    struct adios_index_process_group_struct_v1 * pg;
    uint32_t ts_start_in, ts_end_in;
    uint32_t ts_start_out, ts_end_out;
    if (verbose>1) printf("Check timesteps\n");

    if (out_pg_root == NULL) return 0;

    pg = out_pg_root;
    ts_start_out = pg->time_index;
    while (pg->next) pg = pg->next; 
    ts_end_out = pg->time_index;

    pg = in_pg_root;
    ts_start_in = pg->time_index;
    while (pg->next) pg = pg->next; 
    ts_end_in = pg->time_index;

    if (verbose) { 
        printf("Timesteps in output file: %u - %u\n", ts_start_out, ts_end_out);
        printf("Timesteps in input  file: %u - %u\n", ts_start_in, ts_end_in);
    }

    if (ts_start_in < ts_end_out) {
        fprintf(stderr, "Error: Timesteps in input file should be >= than last timestep in output file\n");
        return 1;
    }

    return 0;
}

/** Recalculate offsets of input to the position in output. 
  * That is, add out_bp->pg_index_offset to all offsets.
  */
void recalc_offsets(void) {
    struct adios_index_process_group_struct_v1 * pg = in_pg_root;
    struct adios_index_var_struct_v1       * vg = in_vars_root;
    struct adios_index_attribute_struct_v1 * ag = in_attrs_root;
    int i;

    // process group offsets
    if (verbose) printf("Recalc offsets for new process groups\n");
    while (pg) {
        pg->offset_in_file += out_bp->pg_index_offset;
        pg = pg->next;
    }

    // variable characteristics offsets
    if (verbose) printf("Recalc offsets for new variable characteristics\n");
    while (vg) {
        for (i=0; i < vg->characteristics_count; i++) {
            vg->characteristics[i].offset += out_bp->pg_index_offset;
        }
        vg = vg->next;
    }

    // attribute characteristics offsets
    if (verbose) printf("Recalc offsets for new attribute characteristics\n");
    while (ag) {
        for (i=0; i < ag->characteristics_count; i++) {
            ag->characteristics[i].offset += out_bp->pg_index_offset;
        }
        ag = ag->next;
    }
}

/** Recover output file if something went bad during the append.
  */
int recover(int f) {
    // Append input indexes into the output indexes
    char * buffer = 0;
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;
    uint64_t index_start = out_bp->pg_index_offset;
    off_t pos;
    ssize_t bytes_written;

    fprintf(stderr, "Recover original output after error in append...\n");

    // seek back to the point where the original index was in the output file
    pos = lseek64 (f, out_bp->pg_index_offset, SEEK_SET); 
    if ( pos == (off_t)-1 ) {
        fprintf(stderr, "  Error at seeking to the original index position in the output file: %s\n", strerror(errno));
        return 1;
    }

    struct adios_index_struct_v1 * idx = adios_alloc_index_v1(0);
    idx->pg_root = out_pg_root;
    idx->vars_root = out_vars_root;
    idx->attrs_root = out_attrs_root;

    // write old index into a buffer
    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset, index_start, 
                          idx);
    if (verbose>1) fprintf(stderr, "  original index size %" PRIu64 " 0x%" PRIx64 "\n", buffer_offset, buffer_offset);
    adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

    adios_clear_index_v1 (idx);
    adios_free_index_v1 (idx);

    // write index buffer out
    if (verbose>1) fprintf(stderr, "  write %" PRIu64 " 0x%" PRIx64 " bytes of indexes...\n", buffer_offset, buffer_offset);
    bytes_written = write (f, buffer, buffer_offset);
    free(buffer);
    if (bytes_written == -1) {
        fprintf(stderr, "  Error at writing original index into output file: %s\n", strerror(errno));
        return 2;
    }
    if (bytes_written != buffer_offset) {
        fprintf(stderr, "  Error: could not write complete original index into output file: %zu bytes from %" PRIu64 "\n"
                "Index in output file will be damaged.\n",
                bytes_written, buffer_offset);
        return 2;
    }

    // truncate file back to original size
    if ( (ftruncate( f, out_bp->pg_index_offset + bytes_written)) == -1) {
        fprintf(stderr, "  Error: Could not truncate output file back to original size: %s\n", strerror(errno));
        return 3;
    }

    return 0;
}

/** Write the output file.
 *   copy data from input file into the output
 *   write out merged indexes
 *   write out version
 *
 *  FIXME: how should we copy if endianness should be changed?
 */
#define COPYBUFSIZE 1024
int append_in_to_out( const char *fileout, const char *filein) {
    int f;
    // open file
    if (verbose) printf("Write data to output file %s\n", fileout);
    errno = 0;
    f = open64( fileout, O_WRONLY);
    if (f == -1) {
        fprintf(stderr, "Cannot open output file %s: %s\n", fileout, strerror(errno));
        return 1;
    }
   
    // copy data
    char buf[COPYBUFSIZE];
    ssize_t bytes_read, bytes_written;
    uint64_t bytes_to_copy = in_bp->pg_index_offset; // all groups in input but no indexes
    uint64_t bytes_copied  = 0;
    off_t    pos;
    lseek64 (in_bp->f, 0, SEEK_SET);
    if (verbose>1) 
        printf("  seek in output to end of group data, %" PRIu64 " (0x%" PRIx64 ") \n",
                out_bp->pg_index_offset, out_bp->pg_index_offset);
    pos = lseek64 (f, out_bp->pg_index_offset, SEEK_SET); // seek behind groups, overwrite old index
    if (pos != out_bp->pg_index_offset) {
            fprintf(stderr, "Error: cannot seek to offset %" PRIu64 " (0x%" PRIx64 ") in file %s: %s\n",
                    out_bp->pg_index_offset, out_bp->pg_index_offset, filein, strerror(errno));
            close(f);
            return 2;
    }
    if (verbose>1) printf("  copy data from input file, %" PRIu64 " bytes\n", in_bp->pg_index_offset);
    while (bytes_copied < bytes_to_copy) {
        bytes_read = read( in_bp->f, (void *)buf, COPYBUFSIZE);
        if (bytes_read < 0) {
            fprintf(stderr, "Error at reading input file %s from offset %" PRIu64 ": %s\n",
                    filein, bytes_copied, strerror(errno));
            recover(f);
            close(f);
            return 3;
        } else if (bytes_read == 0) {
            fprintf(stderr, "Error unexpected end of input file %s at offset %" PRIu64 ": %s\n",
                    filein, bytes_copied, strerror(errno));
            recover(f);
            close(f);
            return 4;
        }
        // check to not to copy too much
        if (bytes_copied + bytes_read > bytes_to_copy)
            bytes_read = bytes_to_copy - bytes_copied;

        bytes_written = write( f, (void *) buf, bytes_read);
        if (bytes_written != bytes_read) {
            fprintf(stderr, "Error: could not write %lld bytes to output file %s at offset %" PRIu64 ": %s\n",
                    (long long)bytes_read, fileout, out_bp->pg_index_offset+bytes_copied, strerror(errno));
            recover(f);
            close(f);
            return 4;
        }

        bytes_copied += bytes_written;
    }
    if (verbose>1) printf("  written %" PRIu64 " (0x%" PRIx64 ") bytes of data into %s\n", bytes_copied, bytes_copied, fileout);

    // Append input indexes into the output indexes
    char * buffer = 0;
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;
    uint64_t index_start =  in_bp->pg_index_offset + out_bp->pg_index_offset;

    struct adios_index_struct_v1 * idx = adios_alloc_index_v1(0);
    idx->pg_root = out_pg_root;
    idx->vars_root = out_vars_root;
    idx->attrs_root = out_attrs_root;

    if (verbose>1) printf("  index starts at %" PRIu64 " (0x%" PRIx64 ")\n", index_start, index_start);

    // merge in old indicies
    adios_merge_index_v1 (idx,
                          in_pg_root, in_vars_root, in_attrs_root, 0);
    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset, index_start, 
                          idx);
    if (verbose>1) printf("  index size %" PRIu64 " 0x%" PRIx64 "\n", buffer_offset, buffer_offset);
    adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);

    // write index buffer out
    if (verbose>1) printf("  write %" PRIu64 " 0x%" PRIx64 " bytes of indexes into %s\n", buffer_offset, buffer_offset, fileout);
    bytes_written = write (f, buffer, buffer_offset);

    if (verbose>1) printf("  written %zu 0x%zx bytes of indexes into %s\n", bytes_written, bytes_written, fileout);

    // clean up
    free(buffer);
    close(f);
    adios_clear_index_v1 (idx);
    adios_free_index_v1 (idx);

    return 0;
}

int bpappend(char *filein, char *fileout) {
    int excode = 0;
    // open output file if exists, read and parse indexes 
    excode = read_indexes( fileout, false );
    if (excode == 128) {
        // output file does not exist, so just copy input to output
        return copy_file(filein, fileout);
    } else if (excode) {
        return excode; // finish here on error
    }
    // we have to append input to output

    // open input file, read and parse indexes 
    excode = read_indexes( filein, true );
    if (excode) 
        return excode; // finish here on error

    // check if timesteps of input >= timesteps already in output 
    excode = check_timesteps();
    if (excode)
        return excode; // finish here on error

    // recalculate offsets for the new data
    recalc_offsets();

    // append input to output file
    excode = append_in_to_out(fileout, filein);  // filein just for error messages

    return excode;
}

void cleanup(void) {
    if (!in_bp)
         adios_posix_close_internal (in_bp);
    if (!outf)
         close (outf);

}
