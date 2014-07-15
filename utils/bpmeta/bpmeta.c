/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
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

#define DIVIDER "========================================================\n"


void print_process_group_index (
                         struct adios_index_process_group_struct_v1 * pg_root
                         );
void print_vars_index (struct adios_index_var_struct_v1 * vars_root);
void print_attributes_index
                         (struct adios_index_attribute_struct_v1 * attrs_root);

int verbose=1; // 1: print summary, 2: print indexes 3: print working info

int main (int argc, char ** argv)
{
    char * filename;
    int rc = 0;
    int nsubfiles=1;

    if (argc < 3)
    {
        fprintf (stderr, "usage: %s <filename> <number of subfiles>\n" ,argv [0]);
        return -1;
    }

    filename = argv [1];
    errno = 0;
    nsubfiles = strtol (argv[2], NULL, 10);
    if (errno != 0 || nsubfiles < 1) {
        fprintf (stderr, "Invalid number of subfiles given\n" ,argv [2]);
        return -1;
    }
    printf ("Create metadata file %s from %d subfiles\n", filename, nsubfiles);

    struct adios_bp_buffer_struct_v1 ** b = 0;
    uint32_t version = 0;
    int idx;
    char fn[256];

    /* Output index structure variables */
    struct adios_index_struct_v1 * index;
    index = adios_alloc_index_v1(1);

    b = malloc (nsubfiles * sizeof (struct adios_bp_buffer_struct_v1*));
    for (idx=0; idx<nsubfiles; idx++) 
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
            printf (DIVIDER);
            printf ("Metadata of %s:\n", fn);
            printf ("BP format version: %d\n", version);
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
           printf ("End of process groups       = %llu\n", b->end_of_pgs);
           printf ("Process Groups Index Offset = %llu\n", b->pg_index_offset);
           printf ("Process Groups Index Size   = %llu\n", b->pg_size);
           printf ("Variable Index Offset       = %llu\n", b->vars_index_offset);
           printf ("Variable Index Size         = %llu\n", b->vars_size);
           printf ("Attribute Index Offset      = %llu\n", b->attrs_index_offset);
           printf ("Attribute Index Size        = %llu\n", b->attrs_size);
         */

        adios_posix_read_process_group_index (b[idx]);
        adios_parse_process_group_index_v1 (b[idx], &new_pg_root);
        print_process_group_index (new_pg_root);

        adios_posix_read_vars_index (b[idx]);
        adios_parse_vars_index_v1 (b[idx], &new_vars_root, NULL, NULL);
        print_vars_index (new_vars_root);

        adios_posix_read_attributes_index (b[idx]);
        adios_parse_attributes_index_v1 (b[idx], &new_attrs_root);
        print_attributes_index (new_attrs_root);

        adios_merge_index_v1 (index, new_pg_root, new_vars_root, new_attrs_root); 

        adios_posix_close_internal (b[idx]);
        adios_shared_buffer_free (b[idx]);

    }

    if (verbose>1) {
        printf (DIVIDER);
        printf ("End of reading all subfiles\n");
    }

    /* Write out the global index */
    char * buffer = 0;
    uint64_t buffer_size = 0;
    uint64_t buffer_offset = 0;
    uint64_t index_start = 0; 
    int err;
    int f;
    uint16_t flag = 0;
    ssize_t bytes_written = 0;

    flag |= ADIOS_VERSION_HAVE_SUBFILE;
    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset
            ,index_start, index);
    if (verbose>2) {
        printf ("buffer=%p size=%lld offset=%lld\n", buffer, buffer_size, buffer_offset);
    }

    adios_write_version_flag_v1 (&buffer, &buffer_size, &buffer_offset, flag);
    if (verbose>2) {
        printf ("buffer=%p size=%lld offset=%lld\n", buffer, buffer_size, buffer_offset);
    }

    f = open (filename, O_CREAT | O_RDWR, 0644);
    if (f == -1)
    {
        fprintf (stderr, "Failed to create metadata file %s: %s\n", 
                 filename, strerror(errno));
        return -1;
    }

    bytes_written = write (f, buffer, (size_t)buffer_offset);
    if (bytes_written == -1) 
    {
        fprintf (stderr, "Failed to write metadata to file %s: %s\n", 
                 filename, strerror(errno));
    } 
    else if (bytes_written != (ssize_t) buffer_offset) 
    {
        fprintf (stderr, "Failed to write metadata of %lld bytes to file %s. "
                "Only wrote %lld bytes\n", buffer_offset, filename, (long long)bytes_written);
    }
    close(f);

    free (b);

    return 0;
}



void print_process_group_index (
                         struct adios_index_process_group_struct_v1 * pg_root
                         )
{
    unsigned int npg=0;
    if (verbose>1) {
        printf (DIVIDER);
        printf ("Process Groups Index:\n");
    }
    while (pg_root)
    {
        if (verbose>1) {
            printf ("Group: %s\n", pg_root->group_name);
            printf ("\tProcess ID: %d\n", pg_root->process_id);
            printf ("\tTime Name: %s\n", pg_root->time_index_name);
            printf ("\tTime: %d\n", pg_root->time_index);
            printf ("\tOffset in File: %llu\n", pg_root->offset_in_file);
        }
        pg_root = pg_root->next;
        npg++;
    }
    if (verbose==1) {
        printf ("Number of process groups: %u\n", npg);
    }
}

void print_vars_index (struct adios_index_var_struct_v1 * vars_root)
{
    unsigned int nvars=0;
    if (verbose>1) {
        printf (DIVIDER);
        printf ("Variable Index:\n");
    }
    while (vars_root)
    {
        if (verbose>1) 
        {
            if (!strcmp (vars_root->var_path, "/"))
            {
                printf ("Var (Group) [ID]: /%s (%s) [%d]\n", vars_root->var_name
                        ,vars_root->group_name, vars_root->id
                       );
            }
            else
            {
                printf ("Var (Group) [ID]: %s/%s (%s) [%d]\n", vars_root->var_path
                        ,vars_root->var_name, vars_root->group_name, vars_root->id
                       );
            }
            const char * typestr = adios_type_to_string_int (vars_root->type);
            printf ("\tDatatype: %s\n", typestr);
            //printf ("\tDatatype: %s\n", adios_type_to_string_int (vars_root->type));
            printf ("\tVars Characteristics: %llu\n"
                    ,vars_root->characteristics_count
                   );
            uint64_t i;
            for (i = 0; i < vars_root->characteristics_count; i++)
            {
                printf ("\tOffset(%llu)", vars_root->characteristics [i].offset);
                printf ("\tPayload Offset(%llu)", vars_root->characteristics [i].payload_offset);
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
                    else
                        printf ("\tValue(\"%s\")", bp_value_to_string (vars_root->type,
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
                            printf ("%llu:%llu:%llu"
                                    ,vars_root->characteristics [i].dims.dims [j * 3 + 0]
                                    ,vars_root->characteristics [i].dims.dims [j * 3 + 1]
                                    ,vars_root->characteristics [i].dims.dims [j * 3 + 2]
                                   );
                        }
                        else
                        {
                            printf ("%llu"
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
                            printf ("%llu:%llu:%llu"
                                    ,dims->dims [j * 3 + 0]
                                    ,dims->dims [j * 3 + 1]
                                    ,dims->dims [j * 3 + 2]
                                   );
                        }
                        else
                        {
                            printf ("%llu"
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
        printf ("Number of variables: %u\n", nvars);
    }
}

void print_attributes_index
                          (struct adios_index_attribute_struct_v1 * attrs_root)
{
    unsigned int nattrs=0;
    if (verbose>1) {
        printf (DIVIDER);
        printf ("Attribute Index:\n");
    }
    while (attrs_root)
    {
        if (verbose>1)
        {
            if (!strcmp (attrs_root->attr_path, "/"))
            {
                printf ("Attribute (Group) [ID]: /%s (%s) [%d]\n"
                        ,attrs_root->attr_name, attrs_root->group_name
                        ,attrs_root->id
                       );
            }
            else
            {
                printf ("Attribute (Group) [ID]: %s/%s (%s) [%d]\n"
                        ,attrs_root->attr_path
                        ,attrs_root->attr_name, attrs_root->group_name
                        ,attrs_root->id
                       );
            }
            printf ("\tDatatype: %s\n", adios_type_to_string_int (attrs_root->type));
            printf ("\tAttribute Characteristics: %llu\n"
                    ,attrs_root->characteristics_count
                   );
            uint64_t i;
            for (i = 0; i < attrs_root->characteristics_count; i++)
            {
                printf ("\t\tOffset(%llu)", attrs_root->characteristics [i].offset);
                printf ("\t\tPayload Offset(%llu)", attrs_root->characteristics [i].payload_offset);
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
                            printf ("%llu:%llu:%llu"
                                    ,attrs_root->characteristics [i].dims.dims [j * 3 + 0]
                                    ,attrs_root->characteristics [i].dims.dims [j * 3 + 1]
                                    ,attrs_root->characteristics [i].dims.dims [j * 3 + 2]
                                   );
                        }
                        else
                        {
                            printf ("%llu"
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
        printf ("Number of attributes: %u\n", nattrs);
    }
}
