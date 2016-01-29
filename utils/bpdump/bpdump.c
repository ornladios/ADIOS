/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "adios_types.h"
#include "adios_internals.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "bp_utils.h"
#include "adios_transforms_common.h" // NCSU ALACRITY-ADIOS
#include "adios_transforms_read.h" // NCSU ALACRITY-ADIOS
//#include "adios_internals.h"

#define DIVIDER "========================================================\n"

struct dump_struct
{
    int do_dump;
    char * dump_var;
    enum ADIOS_FLAG host_language_fortran;
};

struct var_dim
{
    uint16_t id;
    uint64_t rank;
};

void print_process_group_header (uint64_t num
                      ,struct adios_process_group_header_struct_v1 * pg_header
                      );
void print_vars_header (struct adios_vars_header_struct_v1 * vars_header);
void print_var_header (struct adios_var_header_struct_v1 * var_header);
void print_var_payload (struct adios_var_header_struct_v1 * var_header
                       ,struct adios_var_payload_struct_v1 * var_payload
                       ,struct dump_struct * dump
                       ,int var_dims_count
                       ,struct var_dim * var_dims
                       );
void print_attrs_header (
                      struct adios_attributes_header_struct_v1 * attrs_header
                      );
void print_attribute (struct adios_attribute_struct_v1 * attribute);
void print_process_group_index (
                         struct adios_index_process_group_struct_v1 * pg_root
                         );
void print_vars_index (struct adios_index_var_struct_v1 * vars_root);
void print_attributes_index
                         (struct adios_index_attribute_struct_v1 * attrs_root);

/*int print_dataset (int type, int ranks, struct adios_bp_dimension_struct * dims
                  ,void * data
                  );*/

//const char * value_to_string (enum ADIOS_DATATYPES type, void * data);
const char * value_to_string_ptr (enum ADIOS_DATATYPES type, void * data, uint64_t element);

int have_subfiles = 0; // global flag indicating if we variables are stored in subfiles (bpdump does not process subfiles)

int main (int argc, char ** argv)
{
    char * filename;
    int rc = 0;
    struct dump_struct dump;

    if (argc < 2 || argc > 4)
    {
        fprintf (stderr, "usage: %s [-d [var]|--dump [var]] <filename>\n"
                ,argv [0]
                );

        return -1;
    }

    if (argv [1][0] && argv [1][0] == '-')
    {
        if (   !strcmp (argv [1], "-d")
            || !strcmp (argv [1], "--dump")
           )
        {
            dump.do_dump = 1;
            if (argc > 3)
            {
                dump.dump_var = argv [2];
                filename = argv [3];
                printf("%s %s\n",dump.dump_var,filename);
            }
            else
            {
                dump.dump_var = 0;
                filename = argv [2];
                printf("ALLVARS %s\n",filename);
            }
        }
        else
        {
            fprintf (stderr, "usage: %s [-d [var]|--dump [var]] <filename>\n"
                    ,argv [0]
                    );

            return -1;
        }
    }
    else
    {
        filename = argv [1];
        dump.do_dump = 0;
        dump.dump_var = 0;
    }

    have_subfiles = 0;
    struct adios_bp_buffer_struct_v1 * b = 0;
    uint32_t version = 0;

    b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    adios_buffer_struct_init (b);

    rc = adios_posix_open_read_internal (filename, "", b);
    if (!rc)
    {
        fprintf (stderr, "bpdump: file not found: %s\n", filename);

        return -1;
    }

    adios_posix_read_version (b);
    adios_parse_version (b, &version);
    version = version & ADIOS_VERSION_NUM_MASK;
    printf ("BP format version: %d\n", version);
    if (version < 2)
    {
        fprintf (stderr, "bpdump: This version of bpdump can only dump BP format version 2. "
                 "Use an older bpdump from adios 1.6 to dump this file.\n");
        adios_posix_close_internal (b);
        return -1;
    }

    struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_process_group_struct_v1 * pg = 0;
    struct adios_index_var_struct_v1 * vars_root = 0;
    struct adios_index_attribute_struct_v1 * attrs_root = 0;

    printf (DIVIDER);
    printf ("Process Groups Index:\n");
    adios_posix_read_index_offsets (b);
    adios_parse_index_offsets_v1 (b);

    /*
    printf ("End of process groups       = %" PRIu64 "\n", b->end_of_pgs);
    printf ("Process Groups Index Offset = %" PRIu64 "\n", b->pg_index_offset);
    printf ("Process Groups Index Size   = %" PRIu64 "\n", b->pg_size);
    printf ("Variable Index Offset       = %" PRIu64 "\n", b->vars_index_offset);
    printf ("Variable Index Size         = %" PRIu64 "\n", b->vars_size);
    printf ("Attribute Index Offset      = %" PRIu64 "\n", b->attrs_index_offset);
    printf ("Attribute Index Size        = %" PRIu64 "\n", b->attrs_size);
    */

    adios_posix_read_process_group_index (b);
    adios_parse_process_group_index_v1 (b, &pg_root, NULL);
    print_process_group_index (pg_root);

    printf (DIVIDER);
    printf ("Vars Index:\n");
    adios_posix_read_vars_index (b);
    adios_parse_vars_index_v1 (b, &vars_root, NULL, NULL);
    print_vars_index (vars_root);

    printf (DIVIDER);
    printf ("Attributes Index:\n");
    adios_posix_read_attributes_index (b);
    adios_parse_attributes_index_v1 (b, &attrs_root);
    print_attributes_index (attrs_root);

    if (version & ADIOS_VERSION_HAVE_SUBFILE)
    {
        printf (DIVIDER);
        return 0;
    }


    uint64_t pg_num = 0;
    pg = pg_root;
    while (pg)
    {

        pg_num++;
        /* Avoid processing PG's whose offset is beyond the start of index (they are actually in subfiles) */
        /* Note: only variables have subfile index, PG's don't so we need to be this indirect in the test */
        if (pg->offset_in_file >= b->pg_index_offset) 
        {
            printf ("Process Group %" PRIu64 " offset is beyond the footer.", pg_num);
            if (have_subfiles) {
                printf (" It is probably in a subfile but bpdump does not process subfiles.\n");
            } else {
                printf (" Since there are no subfiles, this probably means the BP file is corrupt.\n"
                        "offset=%" PRIu64 "  footer starts at %" PRIu64 "\n", pg->offset_in_file, b->pg_index_offset);
            }
            pg = pg->next;
            continue;
        }

        int var_dims_count = 0;
        struct var_dim * var_dims = 0;

        printf (DIVIDER);

        struct adios_process_group_header_struct_v1 pg_header;
        struct adios_vars_header_struct_v1 vars_header;
        struct adios_attributes_header_struct_v1 attrs_header;

        struct adios_var_header_struct_v1 var_header;
        struct adios_var_payload_struct_v1 var_payload;
        struct adios_attribute_struct_v1 attribute;

        // setup where to read the process group from (and size)
        b->read_pg_offset = pg->offset_in_file;
        if (pg->next)
        {
            b->read_pg_size =   pg->next->offset_in_file
                              - pg->offset_in_file;
        }
        else
        {
            b->read_pg_size =   b->pg_index_offset
                              - pg->offset_in_file;
        }

        adios_posix_read_process_group (b);
        adios_parse_process_group_header_v1 (b, &pg_header);
        print_process_group_header (pg_num, &pg_header);
        //printf ("\tSize of group in fil: %" PRIu64 " bytes\n",  b->read_pg_size);

        adios_parse_vars_header_v1 (b, &vars_header);
        print_vars_header (&vars_header);

        dump.host_language_fortran = pg_header.host_language_fortran;
        int i;
        for (i = 0; i < vars_header.count; i++)
        {
            var_payload.payload = 0;

            adios_parse_var_data_header_v1 (b, &var_header);
            print_var_header (&var_header);

            if ( var_header.dims == 0 || 
                 ( dump.do_dump &&
                   dump.dump_var &&
                   !strcasecmp (dump.dump_var, var_header.name)
                 )
               )
            {
                // add one for string null terminators
                var_payload.payload = malloc (var_header.payload_size + 1);
                adios_parse_var_data_payload_v1 (b, &var_header, &var_payload
                                                ,var_header.payload_size
                                                );
            }
            else
            {
                adios_parse_var_data_payload_v1 (b, &var_header, NULL, 0);
            }

            if (var_header.is_dim == adios_flag_yes)
            {
                var_dims = realloc (var_dims,   (var_dims_count + 1)
                                              * sizeof (struct var_dim)
                                   );

                var_dims [var_dims_count].id = var_header.id;
                var_dims [var_dims_count].rank = *(unsigned int *)
                                                        var_payload.payload;
                var_dims_count++;
            }

            if (dump.do_dump)
            {
                // make sure the buffer is big enough or send in null
                print_var_payload (&var_header, &var_payload, &dump
                                  ,var_dims_count, var_dims
                                  );
            }
            if (var_payload.payload)
            {
                free (var_payload.payload);
                var_payload.payload = 0;
            }
            printf ("\n");
        }

        adios_parse_attributes_header_v1 (b, &attrs_header);
        print_attrs_header (&attrs_header);

        for (i = 0; i < attrs_header.count; i++)
        {
            adios_parse_attribute_v1 (b, &attribute);
            print_attribute (&attribute);
            printf ("\n");
        }

        var_dims_count = 0;
        if (var_dims)
            free (var_dims);
        pg = pg->next;
    }
    printf (DIVIDER);
    printf ("End of %s\n", filename);

    adios_posix_close_internal (b);

    return 0;
}

/* In src/bp_utils.c
const char * value_to_string (enum ADIOS_DATATYPES type, void * data)
{
    static char s [100];
    s [0] = 0;


    switch (type)
    {
        case adios_unsigned_byte:
            sprintf (s, "%u", *(((uint8_t *) data)));
            break;

        case adios_byte:
            sprintf (s, "%d", *(((int8_t *) data)));
            break;

        case adios_short:
            sprintf (s, "%hd", *(((int16_t *) data)));
            break;

        case adios_unsigned_short:
            sprintf (s, "%uh", *(((uint16_t *) data)));
            break;

        case adios_integer:
            sprintf (s, "%d", *(((int32_t *) data)));
            break;

        case adios_unsigned_integer:
            sprintf (s, "%u", *(((uint32_t *) data)));
            break;

        case adios_long:
            sprintf (s, "%" PRId64, *(((int64_t *) data)));
            break;

        case adios_unsigned_long:
            sprintf (s, "%" PRIu64, *(((uint64_t *) data)));
            break;

        case adios_real:
            sprintf (s, "%f", *(((float *) data)));
            break;

        case adios_double:
            sprintf (s, "%le", *(((double *) data)));
            break;

        case adios_long_double:
            sprintf (s, "%Le", *(((long double *) data)));
            break;

        case adios_string:
            sprintf (s, "%s", ((char *) data));
            break;

        case adios_complex:
            sprintf (s, "(%f %f)", *(((float *) data) + 0)
                                 , *(((float *) data) + 1)
                    );
            break;

        case adios_double_complex:
            sprintf (s, "(%lf %lf)", *(((double *) data) + 0)
                                   , *(((double *) data) + 1)
                    );
            break;
    }

    return s;
}
*/

const char * value_to_string_ptr (enum ADIOS_DATATYPES type, void * data, uint64_t element)
{
    static char s [100];
    s [0] = 0;


    switch (type)
    {
        case adios_unsigned_byte:
        {
            uint8_t * p = (uint8_t *) data;
            sprintf (s, "%u", p [element]);
            break;
        }

        case adios_byte:
        {
            int8_t * p = (int8_t *) data;
            sprintf (s, "%d", p [element]);
            break;
        }

        case adios_short:
        {
            int16_t * p = (int16_t *) data;
            sprintf (s, "%hd", p [element]);
            break;
        }

        case adios_unsigned_short:
        {
            uint16_t * p = (uint16_t *) data;
            sprintf (s, "%uh", p [element]);
            break;
        }

        case adios_integer:
        {
            int32_t * p = (int32_t *) data;
            sprintf (s, "%d", p [element]);
            break;
        }

        case adios_unsigned_integer:
        {
            uint32_t * p = (uint32_t *) data;
            sprintf (s, "%u", p [element]);
            break;
        }

        case adios_long:
        {
            int64_t * p = (int64_t *) data;
            sprintf (s, "%" PRId64, p [element]);
            break;
        }

        case adios_unsigned_long:
        {
            uint64_t * p = (uint64_t *) data;
            sprintf (s, "%" PRIu64, p [element]);
            break;
        }

        case adios_real:
        {
            float * p = (float *) data;
            sprintf (s, "%f", p [element]);
            break;
        }

        case adios_double:
        {
            double * p = (double *) data;
            sprintf (s, "%le", p [element]);
            break;
        }

        case adios_long_double:
        {
            long double * p = (long double *) data;
            sprintf (s, "%Le", p [element]);
            break;
        }

        case adios_string:
        case adios_string_array:
        {
            //char * p = (char *) data;
            //sprintf (s, "%s", p [element]);
            fprintf (stderr, "arrays of strings not fully supported\n");
            break;
        }

        case adios_complex:
        {
            float * p = (float *) data;
            sprintf (s, "(%g %g)", p [element * 2 + 0]
                                 , p [element * 2 + 1]
                    );
            break;
        }

        case adios_double_complex:
        {
            double * p = (double *) data;
            sprintf (s, "(%lf %lf)", p [element * 2 + 0]
                                   , p [element * 2 + 1]
                    );
            break;
        }

        case adios_unknown:
            break;
    }

    return s;
}
        
void print_process_group_header (uint64_t num
                      ,struct adios_process_group_header_struct_v1 * pg_header
                      )
{
    struct adios_method_info_struct_v1 * m;

    printf ("Process Group: %" PRIu64 "\n", num);
    printf ("\tGroup Name: %s\n", pg_header->name);
    printf ("\tHost Language Fortran?: %c\n"
           ,(pg_header->host_language_fortran == adios_flag_yes ? 'Y' : 'N')
           );
    printf ("\tCoordination Var Member ID: %d\n", pg_header->coord_var_id);
    printf ("\tTime Name: %s\n", pg_header->time_index_name);
    printf ("\tTime: %d\n", pg_header->time_index);
    printf ("\tMethods used in output: %d\n", pg_header->methods_count);
    m = pg_header->methods;
    while (m)
    {
        printf ("\t\tMethod ID: %d\n", m->id);
        printf ("\t\tMethod Parameters: %s\n", m->parameters);
 
        m = m->next;
    }
}

void print_vars_header (struct adios_vars_header_struct_v1 * vars_header)
{
    printf ("\tVars Count: %u\n", vars_header->count);
}

void print_characteristic_dims(struct adios_index_characteristic_dims_struct_v1 *dims) {
    int j;

    for (j = 0; j < dims->count; j++)
    {
        if (j != 0)
            printf (",");
        if (dims->dims [j * 3 + 1] != 0)
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
}

void print_var_header (struct adios_var_header_struct_v1 * var_header)
{
    int i = 0;

    printf ("\t\tVar Name (ID): %s (%d)\n", var_header->name, var_header->id);
    printf ("\t\tVar Path: %s\n", var_header->path);
    printf ("\t\tDatatype: %s\n", adios_type_to_string_int (var_header->type));
    printf ("\t\tIs Dimension: %c\n"
           ,(var_header->is_dim == adios_flag_yes ? 'Y' : 'N')
           );
    if (var_header->dims)
    {
        struct adios_dimension_struct_v1 * d = var_header->dims;
        printf ("\t\tDimensions:\n");
        while (d)
        {
            printf ("\t\t\tDim %d l:g:o: ", i++);
            if (d->dimension.var_id == 0)
            {
                if (d->dimension.is_time_index == adios_flag_yes)
                {
                    printf ("TIME");
                }
                else
                {
                    printf ("R(%" PRIu64 ")", d->dimension.rank);
                }
            }
            else
            {
                printf ("V(%u)", d->dimension.var_id);
            }

            if (   d->global_dimension.var_id != 0
                || d->global_dimension.rank != 0
               )
            {
                if (d->global_dimension.var_id == 0)
                {
                    if (d->global_dimension.is_time_index == adios_flag_yes)
                    {
                        printf (":TIME");
                    }
                    else
                    {
                        printf (":R(%" PRIu64 ")", d->global_dimension.rank);
                    }
                }
                else
                {
                    printf (":V(%u)", d->global_dimension.var_id);
                }
                if (d->local_offset.var_id == 0)
                {
                    if (d->local_offset.is_time_index == adios_flag_yes)
                    {
                        printf (":TIME");
                    }
                    else
                    {
                        printf (":R(%" PRIu64 ")\n", d->local_offset.rank);
		    }
                }
                else
                {
                    printf (":V(%u)\n", d->local_offset.var_id);
                }
            }
            printf ("\n");

            d = d->next;
	}
    }
    printf ("\t\tCharacteristics:\n");
    printf ("\t\t\tOffset(%" PRIu64 ")", var_header->characteristics.offset);

    /* NCSU - Print min, max */

	if (var_header->type == adios_complex || var_header->type == adios_double_complex)
	{
		uint8_t type;

		if (var_header->type == adios_complex)
			type = adios_double;
		else
			type = adios_long_double;

        if (var_header->characteristics.stats && var_header->characteristics.stats[0][adios_statistic_min].data)
        {
            printf ("\tMin(%s)", bp_value_to_string (type
                                       ,var_header->characteristics.stats[0][adios_statistic_min].data
                                       )
                   );
        }
        if (var_header->characteristics.stats && var_header->characteristics.stats[0][adios_statistic_max].data)
        {
            printf ("\tMax(%s)", bp_value_to_string (type
                                       ,var_header->characteristics.stats[0][adios_statistic_max].data
                                       )
                   );
        }

	}
	else
	{
    	if (var_header->characteristics.stats && var_header->characteristics.stats[0][adios_statistic_min].data)
    	{
    	    printf ("\tMin(%s)", bp_value_to_string (var_header->type
    	                               ,var_header->characteristics.stats[0][adios_statistic_min].data
    	                               )
    	           );
    	}
    	if (var_header->characteristics.stats && var_header->characteristics.stats[0][adios_statistic_max].data)
    	{
    	    printf ("\tMax(%s)", bp_value_to_string (var_header->type
    	                               ,var_header->characteristics.stats[0][adios_statistic_max].data
    	                               )
    	           );
    	}
	}
	
    

    if (var_header->characteristics.value)
    {
        printf ("\t\t\tValue(%s)", bp_value_to_string (var_header->type
                                          ,var_header->characteristics.value
                                          )
               );
    }
    if (var_header->characteristics.dims.count != 0)
    {
        printf ("\t\t\tDims (l:g:o): (");
        print_characteristic_dims(&var_header->characteristics.dims);
        printf(")");
        }

    // NCSU ALACRITY-ADIOS - Adding printing of transform type
    printf ("\t\t\tTransform-type(%hhu = %s)",
            var_header->characteristics.transform.transform_type,
            adios_transform_plugin_primary_xml_alias(var_header->characteristics.transform.transform_type));
    if (var_header->characteristics.transform.transform_type != adios_transform_none) {
        printf ("\t\t\tPre-transform-datatype(%s)", adios_type_to_string_int(var_header->characteristics.transform.pre_transform_type));
        printf ("\t\t\tPre-transform-dims(l:g:o = ");
        print_characteristic_dims(&var_header->characteristics.transform.pre_transform_dimensions);
        printf (")");
        printf ("\t\t\tTransform-metadata-length(%hu)", var_header->characteristics.transform.transform_metadata_len);
    }

    printf ("\n");
}

#if 0
void adios_var_element_count (int rank

                             ,struct adios_bp_dimension_struct * dims
                             ,uint64_t * use_count
                             ,uint64_t * total_count
                             )
{
    int i;

    *use_count = 1;
    *total_count = 1;

    for (i = 0; i < rank; i++)
    {
        *use_count *= dims [i].local_bound;
        *total_count *= dims [i].local_bound;
        int use_size = dims [i].use_upper_bound - dims [i].use_lower_bound + 1;
        int total_size = dims [i].upper_bound - dims [i].lower_bound + 1;

        // adjust for the stride
        if (dims [i].stride <= use_size / 2)
        {
            if (use_size % dims [i].stride == 1) // correct fencepost error
                use_size = use_size / dims [i].stride + 1;
            else
                use_size = use_size / dims [i].stride;
        }
        else
        {
            if (dims [i].stride >= use_size)
                use_size = 1;
            else
                use_size = use_size / dims [i].stride + 1;  // maybe always 2?
        }

        // need to correct for empty/bad arrays
        if (   dims [i].use_upper_bound < dims [i].use_lower_bound
            || dims [i].upper_bound < dims [i].lower_bound
           )
        {
            use_size = 0;
            total_size = 0;
        }

        // number of items in the array
        *use_count *= use_size;
        *total_count *= total_size;
    }
}
#endif

    // for writing out bits in a global way, we would need this piece
static
int increment_dimension (enum ADIOS_FLAG host_language_fortran
                        ,uint64_t element
                        ,int ranks
                        ,uint64_t * dims
                        ,uint64_t * position
                        )
{
    int i;
    int done = 0;
    
    if (element == 0)
    {
        for (i = 0; i < ranks; i++)
        {
            position [i] = 0;
        }   
        done = 1;
    }   
    else  // increment our position
    {
//        if (host_language_fortran == adios_flag_yes)
        {
            i = 0;
            while (!done && i < ranks)
            {
                // if less than max, just increment this dim
                if (position [i] < dims [i] - 1)
                {
                    position [i]++;
                    done = 1;
                }
                else  // reset dim and move to next to increment
                {
                    position [i] = 0;
                    i++;
                }
            }
        }
#if 0
        else
        {
            i = ranks - 1;
            while (!done && i >= 0)
            {
                // if less than max, just increment this dim
                if (position [i] < dims [i])
                {
                    position [i]++;
                    done = 1;
                }
                else  // reset dim and move to next to increment
                {
                    position [i] = 0;
                    i--;
                }
            }
        }
#endif
    }

    return done;

#if 0
    // for writing out bits in a global way, we would need this piece
    // check against bounds
    for (i = 0; i < rank; i++)
    {
        if (   position [i] < dims [i].use_lower_bound
            || position [i] > dims [i].use_upper_bound
           )
        {
            return 0;
        }
        else
        {
            // (pos - use lower) mod stride == 0 == use this element
            if (((position [i] - dims [i].use_lower_bound) % dims [i].stride) != 0)
                return 0;
        }
    }

    return 1;  // we only get here if we are within all bounds
#endif
}

static
int dims_not_max (uint64_t * position, uint64_t * dims, int ranks)
{
    int i;

    for (i = 0; i < ranks; i++)
    {
        if (position [i] != dims [i] - 1)
            return 1;
    }

    return 0;
}

void print_var_payload (struct adios_var_header_struct_v1 * var_header
                       ,struct adios_var_payload_struct_v1 * var_payload
                       ,struct dump_struct * dump
                       ,int var_dims_count
                       ,struct var_dim * var_dims
                       )
{
    if (   dump->do_dump
    	&& dump->dump_var
        && var_header->dims
        && !strcasecmp (dump->dump_var, var_header->name)
       )
    {
        if (var_header->dims)
        {
            uint64_t element = 0;
            int ranks = 0;
            struct adios_dimension_struct_v1 * d = var_header->dims;
            int c = 0;
            uint64_t * position = 0;
            uint64_t * dims = 0;
            int i = 0;

            while (d)
            {
                ranks++;
                d = d->next;
            }

            position = (uint64_t *) malloc (8 * ranks);
            memset (position, 0, 8 * ranks);
            dims = (uint64_t *) malloc (8 * ranks);
            memset (dims, 0, 8 * ranks);

            d = var_header->dims;
            uint64_t * dims_t = dims;

            while (d)
            {
                if (d->dimension.var_id != 0)
                {
                    for (i = 0; i < var_dims_count; i++)
                    {
                        if (var_dims [i].id == d->dimension.var_id)
                        {
                            *dims_t = var_dims [i].rank;
                        }
                    }
                }
                else
                {
                    if (d->dimension.is_time_index == adios_flag_yes)
                    {
                        *dims_t = 1;
                    }
                    else
                    {
                        *dims_t = d->dimension.rank;
                    }
                }

                d = d->next;
                dims_t++;
            }

            while (dims_not_max (position, dims, ranks))
            {
                increment_dimension (dump->host_language_fortran
                                    ,element
                                    ,ranks
                                    ,dims
                                    ,position
                                    );
                if (c > 65)
                {
                    printf ("\n");
                    c = 0;
                }

                c += printf ("[");
                for (i = 0; i < ranks; i++)
                {
                    if (i > 0)
                        c += printf (",%" PRIu64, position [i]);
                    else
                        c += printf ("%" PRIu64, position [i]);
                }
                c += printf ("] ");
                c += printf ("%s ", value_to_string_ptr (var_header->type
                                                 ,var_payload->payload, element
                                                 )
                       );

                element++;
            }
            printf ("\n");

            if (position)
                free (position);
            if (dims)
                free (dims);
            position = 0;
            dims = 0;
        }
    }
    if (!var_header->dims)
    {
        if (var_header->type != adios_string) 
            printf ("\t\t\tValue: %s\n", bp_value_to_string (var_header->type ,var_payload->payload));
        else
            printf ("\t\t\tValue: \"%s\"\n", bp_value_to_string (var_header->type ,var_payload->payload));
    }
}

void print_attrs_header (
                      struct adios_attributes_header_struct_v1 * attrs_header
                      )
{
    printf ("\tAttributes Count: %u\n", attrs_header->count);
}

void print_attribute (struct adios_attribute_struct_v1 * attribute)
{
    int i;
    printf ("\t\tAttribute Name (ID): %s (%d)\n"
           ,attribute->name, attribute->id
           );
    printf ("\t\tAttribute Path: %s\n", attribute->path);
    if (attribute->is_var == adios_flag_yes)
    {
        printf ("\t\tAssociated Var ID: %d\n", attribute->var_id);
    }
    else
    {
        printf ("\t\tDatatype: %s\n", adios_type_to_string_int (attribute->type));
        printf ("\t\t# of elements:   %d\n", attribute->nelems);
        printf ("\t\tLength in bytes: %d\n", attribute->length);
        char * p = (char *) attribute->value;
        int elemsize = (int) adios_get_type_size (attribute->type, attribute->value);
        if (attribute->nelems == 1) {
            printf ("\t\tValue: %s\n", bp_value_to_string (attribute->type ,attribute->value));
        } else {
            printf ("\t\tValues: %s", bp_value_to_string (attribute->type, p));
            for (i=1; i<attribute->nelems; i++) {
                p += elemsize;
                printf (", %s", bp_value_to_string (attribute->type, p));
            }
            printf ("\n");
        }
    }
}

void print_process_group_index (
                         struct adios_index_process_group_struct_v1 * pg_root
                         )
{
    while (pg_root)
    {
        printf ("Group: %s\n", pg_root->group_name);
        printf ("\tProcess ID: %d\n", pg_root->process_id);
        printf ("\tTime Name: %s\n", pg_root->time_index_name);
        printf ("\tTime: %d\n", pg_root->time_index);
        printf ("\tOffset in File: %" PRIu64 "\n", pg_root->offset_in_file);

        pg_root = pg_root->next;
    }
}

void print_vars_index (struct adios_index_var_struct_v1 * vars_root)
{
    while (vars_root)
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
        printf ("\tVars Characteristics: %" PRIu64 "\n"
               ,vars_root->characteristics_count
               );
        uint64_t i;
        for (i = 0; i < vars_root->characteristics_count; i++)
        {
            printf ("\tOffset(%" PRIu64 ")", vars_root->characteristics [i].offset);
            printf ("\tPayload Offset(%" PRIu64 ")", vars_root->characteristics [i].payload_offset);
            printf ("\tFile Index(%d)", vars_root->characteristics [i].file_index);
            printf ("\tTime Index(%d)", vars_root->characteristics [i].time_index);

            if (vars_root->characteristics [i].file_index != (uint32_t)-1) { 
                have_subfiles = 1;
            }

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

        vars_root = vars_root->next;
    }
}

void print_attributes_index
                          (struct adios_index_attribute_struct_v1 * attrs_root)
{
    while (attrs_root)
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
        printf ("\tAttribute Characteristics: %" PRIu64 "\n"
               ,attrs_root->characteristics_count
               );
        uint64_t i;
        for (i = 0; i < attrs_root->characteristics_count; i++)
        {
            printf ("\t\tOffset(%" PRIu64 ")", attrs_root->characteristics [i].offset);
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
                if (attrs_root->nelems == 1) {
                    printf ("\t\tValue(%s)", bp_value_to_string (attrs_root->type, attrs_root->characteristics [i].value));
                } else {
                    char * p = (char *) attrs_root->characteristics [i].value;
                    int elemsize = (int) adios_get_type_size (attrs_root->type, p);
                    int k;
                    printf ("\t\tValues(%s", bp_value_to_string (attrs_root->type, p));
                    for (k=1; k<attrs_root->nelems; k++) {
                        p += elemsize;
                        printf (", %s", bp_value_to_string (attrs_root->type, p));
                    }
                    printf (")\n");
                }
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

        attrs_root = attrs_root->next;
    }
}
