/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include "config.h"

#ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <glob.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include "adios_types.h"
#include "adios_internals.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "bp_utils.h"
#include "adios_transforms_common.h" // NCSU ALACRITY-ADIOS
#include "adios_transforms_read.h" // NCSU ALACRITY-ADIOS

#if HAVE_PTHREAD
#   include "pthread.h"
#endif

#define DIVIDER "========================================================\n"

// User arguments
int verbose=0;   // 1: print summary, 2: print indexes 3: print working info
int nthreads=1;  // Number of threads to use (main counts as 1 thread)
char * filename; // process 'filename'.dir/'filename'.NNN subfiles and 
                 //   generate metadata file 'filename'
int nsubfiles=0; // number of subfiles to process

struct option options[] = {
    {"help",                 no_argument,          NULL,    'h'},
    {"verbose",              no_argument,          NULL,    'v'},
    {"nsubfiles",            required_argument,    NULL,    'n'},
#if HAVE_PTHREAD
    {"nthreads",             required_argument,    NULL,    't'},
#endif
    {NULL,                   0,                    NULL,    0}
};

#if HAVE_PTHREAD
static const char *optstring = "hvn:t:";
#else
static const char *optstring = "hvn:";
#endif

// help function
void display_help() {
    printf ("usage: bpmeta [OPTIONS] -n <N> <filename>\n"
            "\nbpmeta processes <filename>.dir/<filename>.<nnn> subfiles and\n"
            "generates a metadata file <filename>.\n"
            "\nIt is used to generate the missing metadata file after using\n"
            "the MPI_AGGREGATE output method with 'have_metada_file=0' option.\n" 
            "\n"
            "  --nsubfiles | -n <N>   The number of subfiles to process in\n"
            "                           <filename>.dir\n"
#if HAVE_PTHREAD
            "  --nthreads  | -t <T>   Parallel reading with <T> threads.\n"
            "                           The main thread is counted in.\n"
#endif
            "\n"
            "Help options\n"
            "  --help      | -h       Print this help.\n"
            "  --verbose   | -v       Print log about what this program is doing.\n"
            "                           Use multiple -v to increase logging level.\n"
            "Typical use: bpmeta -t 16 -n 1024 mydata.bp\n"
           );
}


/* Global variables among threads */
struct adios_bp_buffer_struct_v1 ** b = 0;
  /* sub-index structure variables */
struct adios_index_struct_v1 ** subindex;

int process_subfiles (int tid, int startidx, int endidx);
int write_index (struct adios_index_struct_v1 * index, char * fname);
int get_nsubfiles (char *filename);
void print_pg_index ( int tid, struct adios_index_process_group_struct_v1 * pg_root);
void print_variable_index (int tid, struct adios_index_var_struct_v1 * vars_root);
void print_attribute_index (int tid,  struct adios_index_attribute_struct_v1 * attrs_root);

#if HAVE_PTHREAD
struct thread_args 
{
    int tid;
    int startidx;
    int endidx;
};

void * thread_main (void *arg)
{
    struct thread_args *targ = (struct thread_args *) arg;
    process_subfiles (targ->tid, targ->startidx, targ->endidx);
    pthread_exit(NULL);
    return NULL; // just to avoid compiler warning
}
#endif


int main (int argc, char ** argv)
{
    long int tmp;
    int c;
    while ((c = getopt_long(argc, argv, optstring, options, NULL)) != -1) {
        switch (c) {
            case 'n':
            case 't':
                errno = 0;
                tmp = strtol(optarg, (char **)NULL, 0);
                if (errno) {
                    fprintf(stderr, "Error: could not convert -%c value: %s\n", 
                            c, optarg);
                    return 1;
                }
                if (c == 'n')
                    nsubfiles=tmp;
                else 
                    nthreads=tmp;
                break;

            case 'h':
                display_help();
                return 0;
                break;

            case 'v':
                verbose++;
                break;

            default:
                printf("Unrecognized argument: %s\n", optarg);
                break;

        }
    }

    /* Check if we have a file defined */
    if (optind >= argc) {
        printf ("Missing file name\n");
        display_help();
        return 1;
    }

    filename = strdup(argv[optind++]);

    if (nsubfiles < 1)
        nsubfiles = get_nsubfiles (filename);
    if (nsubfiles < 1) {
        printf ("Cannot determine the number of subfiles. To avoid this problem, "
                "provide the number of subfiles manually with the -n <N> option.\n");
        return -1;
    }

    if (nthreads < 1) 
        nthreads = 1;
    if (nthreads > nsubfiles) {
        printf ("Warning: asked for processing %d subfiles using %d threads. "
                "We will utilize only %d threads.\n", 
                nsubfiles, nthreads, nsubfiles);
        nthreads = nsubfiles;
    }

    if (verbose>1)
        printf ("Create metadata file %s from %d subfiles using %d threads\n", 
                filename, nsubfiles, nthreads);

    /* Initialize global variables */
    b = malloc (nsubfiles * sizeof (struct adios_bp_buffer_struct_v1*));
    subindex = malloc (nthreads * sizeof (struct adios_index_struct_v1*));

    /* Split the processing work among T threads */
    int tid;

#if HAVE_PTHREAD

    pthread_t *thread = (pthread_t *) malloc (nthreads * sizeof(pthread_t));
    struct thread_args *targs = (struct thread_args*) 
                                  malloc (nthreads * sizeof(struct thread_args));
    int K = nsubfiles/nthreads; // base number of files to be processed by one thread
    int L = nsubfiles%nthreads; // this many threads processes one more subfile
    int startidx, endidx;
    int rc;
    //printf ("K=%d L=%d\n", K, L);
    endidx = -1;
    for (tid=0; tid<nthreads; tid++)
    {
        startidx = endidx + 1;
        endidx = startidx + K - 1;
        if (tid < L) {
            endidx++;
        }
        targs[tid].tid = tid;
        targs[tid].startidx = startidx;
        targs[tid].endidx = endidx;
        if (verbose)
            printf ("Process subfiles from %d to %d with thread %d\n", 
                    targs[tid].startidx, targs[tid].endidx, targs[tid].tid);

        if (tid < nthreads-1) {
            /* Start worker thread. */
            pthread_attr_t attr;
            pthread_attr_init (&attr);
            pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
            rc = pthread_create (&thread[tid], &attr, thread_main, &targs[tid]);
            if (rc) {
                printf ("ERROR: Thread %d: Cannot create thread, err code = %d\n", 
                        tid, rc);
            }
            pthread_attr_destroy(&attr);
        } 
        else
        {
            // last "thread" is the main thread
            process_subfiles (tid, startidx, endidx);
        }
    }
    // wait here for everyone to finish
    for (tid=0; tid<nthreads-1; tid++)
    {
        void *status;
        rc = pthread_join (thread[tid], &status);
        if (rc) {
            printf ("ERROR: Thread %d: Cannot join thread, err code = %d\n", tid, rc);
        } else {
            if (verbose>1)
                printf ("Thread %d: Joined thread.\n", tid);
        }
    }
    free (targs);
    free (thread);

#else /* non-threaded version */

    nthreads = 1;
    process_subfiles (0, 0, nsubfiles-1);

#endif

    /* Merge the T indexes into the global output index */
    struct adios_index_struct_v1 * globalindex;
    globalindex = adios_alloc_index_v1(1);
    for (tid=0; tid<nthreads; tid++)
    {
        adios_merge_index_v1 (globalindex, 
                              subindex[tid]->pg_root, 
                              subindex[tid]->vars_root, 
                              subindex[tid]->attrs_root, 1); 
    }
    write_index (globalindex, filename);

    /* Clean-up */
    adios_clear_index_v1 (globalindex);
    /*... already cleaned-up by globalindex clearing
    for (tid=0; tid<nthreads; tid++) {
        adios_clear_index_v1 (subindex[tid]);
        free (subindex[tid]);
    }
    */
    free (subindex);
    free (globalindex);
    free (b);
    return 0;
}

int write_index (struct adios_index_struct_v1 * index, char * fname)
{
    /* Write out the global index */
    char * buffer = 0;
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;
    uint64_t index_start = 0; 
    int f;
    uint16_t flag = 0;
    ssize_t bytes_written = 0;

    flag |= ADIOS_VERSION_HAVE_SUBFILE;
    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
            ,index_start, index);
    if (verbose>2) {
        printf ("buffer=%p size=%" PRId64 " offset=%" PRId64 "\n", buffer, buffer_size, buffer_offset);
    }

    adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);
    if (verbose>2) {
        printf ("buffer=%p size=%" PRId64 " offset=%" PRId64 "\n", buffer, buffer_size, buffer_offset);
    }

    f = open (fname, O_CREAT | O_RDWR, 0644);
    if (f == -1)
    {
        fprintf (stderr, "Failed to create metadata file %s: %s\n", 
                 fname, strerror(errno));
        return -1;
    }

    bytes_written = write (f, buffer, (size_t)buffer_offset);
    if (bytes_written == -1) 
    {
        fprintf (stderr, "Failed to write metadata to file %s: %s\n", 
                 fname, strerror(errno));
    } 
    else if (bytes_written != (ssize_t) buffer_offset) 
    {
        fprintf (stderr, "Failed to write metadata of %" PRId64 " bytes to file %s. "
                "Only wrote %lld bytes\n", buffer_offset, fname, (long long)bytes_written);
    }
    close(f);
    return 0;
}

int process_subfiles (int tid, int startidx, int endidx)
{
    char fn[256];
    uint32_t version = 0;
    int idx;
    int rc = 0;

    subindex[tid] = adios_alloc_index_v1(1);

    for (idx=startidx; idx<=endidx; idx++) 
    {
        b[idx] = malloc (sizeof (struct adios_bp_buffer_struct_v1));
        adios_buffer_struct_init (b[idx]);

        snprintf (fn, 256, "%s.dir/%s.%d", filename, filename, idx);
        rc = adios_posix_open_read_internal (fn, "", b[idx]);
        if (!rc)
        {
            fprintf (stderr, "bpmeta: file not found: %s\n", fn);
            return -1;
        }

        adios_posix_read_version (b[idx]);
        adios_parse_version (b[idx], &version);
        version = version & ADIOS_VERSION_NUM_MASK;
        if (verbose) {
            //printf (DIVIDER);
            printf ("Thread %d: Metadata of %s:\n", tid, fn);
            printf ("Thread %d: BP format version: %d\n", tid, version);
        }
        if (version < 2)
        {
            fprintf (stderr, "bpmeta: This version of bpmeta can only work with BP format version 2 and up. "
                    "Use an older bpmeta from adios 1.6 to work with this file.\n");
            adios_posix_close_internal (b[idx]);
            return -1;
        }

        struct adios_index_process_group_struct_v1 * new_pg_root = 0;
        struct adios_index_var_struct_v1 * new_vars_root = 0;
        struct adios_index_attribute_struct_v1 * new_attrs_root = 0;

        adios_posix_read_index_offsets (b[idx]);
        adios_parse_index_offsets_v1 (b[idx]);

        /*
           printf ("End of process groups       = %" PRIu64 "\n", b->end_of_pgs);
           printf ("Process Groups Index Offset = %" PRIu64 "\n", b->pg_index_offset);
           printf ("Process Groups Index Size   = %" PRIu64 "\n", b->pg_size);
           printf ("Variable Index Offset       = %" PRIu64 "\n", b->vars_index_offset);
           printf ("Variable Index Size         = %" PRIu64 "\n", b->vars_size);
           printf ("Attribute Index Offset      = %" PRIu64 "\n", b->attrs_index_offset);
           printf ("Attribute Index Size        = %" PRIu64 "\n", b->attrs_size);
         */

        adios_posix_read_process_group_index (b[idx]);
        adios_parse_process_group_index_v1 (b[idx], &new_pg_root, NULL);
        print_pg_index (tid, new_pg_root);

        adios_posix_read_vars_index (b[idx]);
        adios_parse_vars_index_v1 (b[idx], &new_vars_root, NULL, NULL);
        print_variable_index (tid, new_vars_root);

        adios_posix_read_attributes_index (b[idx]);
        adios_parse_attributes_index_v1 (b[idx], &new_attrs_root);
        print_attribute_index (tid, new_attrs_root);

        adios_merge_index_v1 (subindex[tid], new_pg_root, new_vars_root, new_attrs_root, 1); 

        adios_posix_close_internal (b[idx]);
        adios_shared_buffer_free (b[idx]);

    }

    if (verbose>1) {
        //printf (DIVIDER);
        printf ("Thread %d: End of reading all subfiles\n", tid);
    }


    return 0;
}


int get_nsubfiles (char *filename)
{
    char pattern[256];
    glob_t g;
    int err,ret;

    snprintf (pattern, 256, "%s.dir/%s.*", filename, filename);
    err = glob (pattern, GLOB_ERR | GLOB_NOSORT, NULL, &g);
    if (!err) {
        ret = g.gl_pathc;
    } else {
        switch (err) {
            case GLOB_NOMATCH:
                printf ("ERROR: No matching file found for %s\n", pattern);
                break;

            case GLOB_NOSPACE:
                printf ("ERROR: Not enough memory for running glob for pattern %s\n", pattern);
                break;

            case GLOB_ABORTED:
                printf ("ERROR: glob was aborted by a reading error for pattern %s\n", pattern);
                printf ("errno = %d: %s\n", errno, strerror(errno));
                break;
        }
        ret = 0;
    }
    return ret;
}

void print_pg_index (int tid, struct adios_index_process_group_struct_v1 * pg_root)
{
    unsigned int npg=0;
    if (verbose>1) {
        //printf (DIVIDER);
        printf ("Thread %d:   Process Groups Index:\n", tid);
    }
    while (pg_root)
    {
        if (verbose>1) {
                printf ("Thread %d:   Group: %s\n", tid, pg_root->group_name);
                printf ("Thread %d:   \tProcess ID: %d\n", tid, pg_root->process_id);
                printf ("Thread %d:   \tTime Name: %s\n", tid, pg_root->time_index_name);
                printf ("Thread %d:   \tTime: %d\n", tid, pg_root->time_index);
                printf ("Thread %d:   \tOffset in File: %" PRIu64 "\n", tid, pg_root->offset_in_file);
        }
        pg_root = pg_root->next;
        npg++;
    }
    if (verbose==1) {
        printf ("Thread %d: Number of process groups: %u\n", tid, npg);
    }
}

void print_variable_index (int tid, struct adios_index_var_struct_v1 * vars_root)
{
    unsigned int nvars=0;
    if (verbose>1) {
        //printf (DIVIDER);
        printf ("Thread %d: Variable Index:\n", tid);
    }
    while (vars_root)
    {
        if (verbose>1) 
        {
            if (!strcmp (vars_root->var_path, "/"))
            {
                printf ("Thread %d:   Var (Group) [ID]: /%s (%s) [%d]\n", 
                        tid, vars_root->var_name,vars_root->group_name, vars_root->id
                       );
            }
            else
            {
                printf ("Thread %d:   Var (Group) [ID]: %s/%s (%s) [%d]\n", 
                        tid, vars_root->var_path, 
                        vars_root->var_name, vars_root->group_name, vars_root->id
                       );
            }
            const char * typestr = adios_type_to_string_int (vars_root->type);
            printf ("Thread %d: \tDatatype: %s\n", tid, typestr);
            printf ("Thread %d: \tVars Characteristics: %" PRIu64 "\n",
                    tid, vars_root->characteristics_count
                   );
            uint64_t i;
            for (i = 0; i < vars_root->characteristics_count; i++)
            {
                printf ("Thread %d: \tOffset(%" PRIu64 ")", tid, vars_root->characteristics [i].offset);
                printf ("\tPayload Offset(%" PRIu64 ")", vars_root->characteristics [i].payload_offset);
                printf ("\tFile Index(%d)", vars_root->characteristics [i].file_index);
                printf ("\tTime Index(%d)", vars_root->characteristics [i].time_index);

                /* NCSU - Print min, max */
                if (vars_root->type == adios_complex || vars_root->type == adios_double_complex)
                {
                    uint8_t type;

                    if (vars_root->type == adios_complex)
                        type = adios_double;
                    else
                        type = adios_long_double;


                    if (vars_root->characteristics [i].stats && vars_root->characteristics [i].stats[0][adios_statistic_min].data)
                    {
                        printf ("\tMin(%s)", bp_value_to_string (type
                                    ,vars_root->characteristics [i].stats[0][adios_statistic_min].data
                                    )
                               );
                    }
                    if (vars_root->characteristics [i].stats && vars_root->characteristics [i].stats[0][adios_statistic_max].data)
                    {
                        printf ("\tMax(%s)", bp_value_to_string (type
                                    ,vars_root->characteristics [i].stats[0][adios_statistic_max].data
                                    )
                               );
                    }
                }
                else
                {
                    if (vars_root->characteristics [i].stats && vars_root->characteristics [i].stats[0][adios_statistic_min].data)
                    {
                        printf ("\tMin(%s)", bp_value_to_string (vars_root->type
                                    ,vars_root->characteristics [i].stats[0][adios_statistic_min].data
                                    )
                               );
                    }
                    if (vars_root->characteristics [i].stats && vars_root->characteristics [i].stats[0][adios_statistic_max].data)
                    {
                        printf ("\tMax(%s)", bp_value_to_string (vars_root->type
                                    ,vars_root->characteristics [i].stats[0][adios_statistic_max].data
                                    )
                               );
                    }
                }
                //*/

                if (vars_root->characteristics [i].value)
                {
                    if (vars_root->type != adios_string)
                        printf ("\tValue(%s)", bp_value_to_string (vars_root->type,
                                    vars_root->characteristics [i].value));
                }
                if (vars_root->characteristics [i].dims.count != 0)
                {
                    int j;

                    printf ("\tDims (l:g:o): (");
                    for (j = 0; j < vars_root->characteristics [i].dims.count; j++)
                    {
                        if (j != 0)
                            printf (",");
                        if (  vars_root->characteristics [i].dims.dims [j * 3 + 1]
                                != 0
                           )
                        {
                            printf ("%" PRIu64 ":%" PRIu64 ":%" PRIu64
                                    ,vars_root->characteristics [i].dims.dims [j * 3 + 0]
                                    ,vars_root->characteristics [i].dims.dims [j * 3 + 1]
                                    ,vars_root->characteristics [i].dims.dims [j * 3 + 2]
                                   );
                        }
                        else
                        {
                            printf ("%" PRIu64
                                    ,vars_root->characteristics [i].dims.dims [j * 3 + 0]
                                   );
                        }
                    }
                    printf (")");
                }

                // NCSU ALACRITY-ADIOS - Print transform info
                if (vars_root->characteristics[i].transform.transform_type != adios_transform_none) {
                    struct adios_index_characteristic_transform_struct *transform = &vars_root->characteristics[i].transform;
                    struct adios_index_characteristic_dims_struct_v1 *dims = &transform->pre_transform_dimensions;
                    int j;

                    printf ("\tTransform type: %s (ID = %hhu)", adios_transform_plugin_desc(transform->transform_type), transform->transform_type);
                    printf ("\tPre-transform datatype: %s", adios_type_to_string_int(transform->pre_transform_type));
                    printf ("\tPre-transform dims (l:g:o): (");
                    for (j = 0; j < dims->count; j++)
                    {
                        if (j != 0)
                            printf (",");
                        if (  dims->dims [j * 3 + 1]
                                != 0
                           )
                        {
                            printf ("%" PRIu64 ":%" PRIu64 ":%" PRIu64
                                    ,dims->dims [j * 3 + 0]
                                    ,dims->dims [j * 3 + 1]
                                    ,dims->dims [j * 3 + 2]
                                   );
                        }
                        else
                        {
                            printf ("%" PRIu64
                                    ,dims->dims [j * 3 + 0]
                                   );
                        }
                    }
                    printf (")");
                }

                printf ("\n");
            }

        }
        vars_root = vars_root->next;
        nvars++;
    }
    if (verbose==1) {
        printf ("Thread %d: Number of variables: %u\n", tid, nvars);
    }
}

void print_attribute_index (int tid, struct adios_index_attribute_struct_v1 * attrs_root)
{
    unsigned int nattrs=0;
    if (verbose>1) {
        //printf (DIVIDER);
        printf ("Thread %d: Attribute Index:\n", tid);
    }
    while (attrs_root)
    {
        if (verbose>1)
        {
            if (!strcmp (attrs_root->attr_path, "/"))
            {
                printf ("Thread %d:   Attribute (Group) [ID]: /%s (%s) [%d]\n",
                        tid, attrs_root->attr_name, attrs_root->group_name,
                        attrs_root->id
                       );
            }
            else
            {
                printf ("Thread %d:   Attribute (Group) [ID]: %s/%s (%s) [%d]\n",
                        tid, attrs_root->attr_path, attrs_root->attr_name, 
                        attrs_root->group_name, attrs_root->id
                       );
            }
            printf ("Thread %d: \tDatatype: %s\n", tid, 
                        adios_type_to_string_int (attrs_root->type));
            printf ("Thread %d: \tAttribute Characteristics: %" PRIu64 "\n", tid,
                        attrs_root->characteristics_count
                   );
            uint64_t i;
            for (i = 0; i < attrs_root->characteristics_count; i++)
            {
                printf ("Thread %d: \t\tOffset(%" PRIu64 ")", tid,
                            attrs_root->characteristics [i].offset);
                printf ("\t\tPayload Offset(%" PRIu64 ")", attrs_root->characteristics [i].payload_offset);
                printf ("\t\tFile Index(%d)", attrs_root->characteristics [i].file_index);
                printf ("\t\tTime Index(%d)", attrs_root->characteristics [i].time_index);

                /* NCSU - Print min, max  */
                if (attrs_root->type == adios_complex || attrs_root->type == adios_double_complex)
                {
                    uint8_t type;
                    if (attrs_root->type == adios_complex)
                        type = adios_double;
                    else
                        type = adios_long_double;

                    if (attrs_root->characteristics [i].stats && attrs_root->characteristics [i].stats[0][adios_statistic_min].data)
                    {
                        printf ("\tMin(%s)", bp_value_to_string (type
                                    ,attrs_root->characteristics [i].stats[0][adios_statistic_min].data
                                    )
                               );
                    }
                    if (attrs_root->characteristics [i].stats && attrs_root->characteristics [i].stats[0][adios_statistic_max].data)
                    {
                        printf ("\tMax(%s)", bp_value_to_string (type
                                    ,attrs_root->characteristics [i].stats[0][adios_statistic_max].data
                                    )
                               );
                    }
                }
                else
                {
                    if (attrs_root->characteristics [i].stats && attrs_root->characteristics [i].stats[0][adios_statistic_min].data)
                    {
                        printf ("\tMin(%s)", bp_value_to_string (attrs_root->type
                                    ,attrs_root->characteristics [i].stats[0][adios_statistic_min].data
                                    )
                               );
                    }
                    if (attrs_root->characteristics [i].stats && attrs_root->characteristics [i].stats[0][adios_statistic_max].data)
                    {
                        printf ("\tMax(%s)", bp_value_to_string (attrs_root->type
                                    ,attrs_root->characteristics [i].stats[0][adios_statistic_max].data
                                    )
                               );
                    }
                }

                if (attrs_root->characteristics [i].value)
                {
                    printf ("\t\tValue(%s)", bp_value_to_string (attrs_root->type
                                ,attrs_root->characteristics [i].value
                                )
                           );
                }
                if (attrs_root->characteristics [i].var_id)
                {
                    printf ("\t\tVar(%u)", attrs_root->characteristics [i].var_id);
                }
                if (attrs_root->characteristics [i].dims.count != 0)
                {
                    int j;

                    printf ("\t\tDims (l:g:o): (");
                    for (j = 0; j < attrs_root->characteristics [i].dims.count; j++)
                    {
                        if (j != 0)
                            printf (",");
                        if (  attrs_root->characteristics [i].dims.dims [j * 3 + 1]
                                != 0
                           )
                        {
                            printf ("%" PRIu64 ":%" PRIu64 ":%" PRIu64
                                    ,attrs_root->characteristics [i].dims.dims [j * 3 + 0]
                                    ,attrs_root->characteristics [i].dims.dims [j * 3 + 1]
                                    ,attrs_root->characteristics [i].dims.dims [j * 3 + 2]
                                   );
                        }
                        else
                        {
                            printf ("%" PRIu64
                                    ,attrs_root->characteristics [i].dims.dims [j * 3 + 0]
                                   );
                        }
                    }
                    printf (")");
                }
                printf ("\n");
            }
        }

        attrs_root = attrs_root->next;
        nattrs++;
    }
    if (verbose==1) {
        printf ("Thread %d: Number of attributes: %u\n", tid, nattrs);
    }
}
