//#include "binpack-general.h"
//#include "br-utils.h"
#include <sys/types.h>
#include "adios_types.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

#define DIVIDER "========================================================\n"

struct dump_struct
{
    int do_dump;
    char * dump_var;
};

void print_process_group_header (uint64_t num
                      ,struct adios_process_group_header_struct_v1 * pg_header
                      );
void print_vars_header (struct adios_vars_header_struct_v1 * vars_header);
void print_var_header (struct adios_var_header_struct_v1 * var_header);
void print_var_payload (struct adios_var_header_struct_v1 * var_header
                       ,struct adios_var_payload_struct_v1 * var_payload
                       ,struct dump_struct * dump
                       );
void print_attrs_header (
                      struct adios_attributes_header_struct_v1 * attrs_header
                      );
void print_attribute (struct adios_attribute_struct_v1 * attribute);
void print_process_group_index (
                         struct adios_index_process_group_struct_v1 * pg_root
                         );
void print_vars_index (struct adios_index_var_struct_v1 * vars_root);

int print_dataset (int type, int ranks, struct adios_bp_dimension_struct * dims
                  ,void * data
                  );

const char * value_to_string (enum ADIOS_DATATYPES type, void * data, uint64_t element);

int main (int argc, char ** argv)
{
    char * filename;
    char * var;
    int i = 0;
    int rc = 0;
    uint64_t element_size = 0;
    struct adios_bp_element_struct * element = NULL;
    struct dump_struct dump;

    if (argc < 2)
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
            if (argc > 2)
            {
                dump.dump_var = argv [2];
                filename = argv [3];
                printf("%s %s\n",dump.dump_var,filename);
            }
            else
            {
                dump.dump_var = 0;
                filename = argv [2];
                printf("%s %s\n",dump.dump_var,filename);
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

    struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_process_group_struct_v1 * pg = 0;
    struct adios_index_var_struct_v1 * vars_root = 0;

    printf (DIVIDER);
    printf ("Process Groups Index:\n");
    adios_posix_read_index_offsets (b);
    adios_parse_index_offsets_v1 (b);

    adios_posix_read_process_group_index (b);
    adios_parse_process_group_index_v1 (b, &pg_root);
    print_process_group_index (pg_root);

    printf (DIVIDER);
    printf ("Vars Index:\n");
    adios_posix_read_vars_index (b);
    adios_parse_vars_index_v1 (b, &vars_root);
    print_vars_index (vars_root);

    uint64_t element_num = 1;
    pg = pg_root;
    while (pg)
    {
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
        print_process_group_header (element_num++, &pg_header);

        adios_parse_vars_header_v1 (b, &vars_header);
        print_vars_header (&vars_header);

        int i;
        for (i = 0; i < vars_header.count; i++)
        {
            adios_parse_var_data_header_v1 (b, &var_header);
            print_var_header (&var_header);

            if (dump.do_dump)
            {
                // make sure the buffer is big enough or send in null
                adios_parse_var_data_payload_v1 (b, &var_header, &var_payload);
                print_var_payload (&var_header, &var_payload, &dump);
            }
            else
            {
                adios_parse_var_data_payload_v1 (b, &var_header, NULL);
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

        pg = pg->next;
    }
    printf (DIVIDER);
    printf ("End of %s\n", filename);

    adios_posix_close_internal (b);

    return 0;
}

#if 0
int print_dataset (int type, int ranks, struct adios_bp_dimension_struct * dims
                  ,void * data
                  )
{
    printf("My Value: \n");
    int use_element_count = 1;
    int total_element_count = 1;
    int e = 0; // which element we are printing
    int i,j;
    int position[ranks];
    for (i = 0; i < ranks; i++)
    {
        use_element_count *= (  dims [i].local_bound
                             );
        position[i]=dims[i].local_bound; 
        total_element_count *= dims [i].local_bound;
    } 
    printf ("DIMS: dims[%d][%d]\n",position[0],position[1]);
    printf ("ranks=%d\n",ranks);
    for (j = 0; j < total_element_count; j++)
    {
        printf ("%s ", value_to_string (type, data, e));
            switch (type)
            {
                case bp_uchar:
                    printf ("%c ", (((unsigned char *) data) [e]));
                    break;

                case bp_char:
                    printf ("%c ", (((char *) data) [e]));
                    break;

                case bp_int:
                    printf ("%d ", (((int *) data) [e]));
                    break;

                case bp_float:
                    printf ("(%d: %e) ", j,(((float *) data) [e]));
                    break;

                case bp_double:
                    printf ("(%d: %le) ",j, (((double *) data) [e]));
                    break;

                case bp_string:
                    break;

                case bp_longlong: // adios_long
                    printf ("%lld ", (((int64_t *) data) [e]));
                    break;

                case bp_complex: // adios_complex
                    printf ("(%lf %lf) ", ((double *) data) [e * 2 + 0]
                                        , ((double *) data) [e * 2 + 1]
                           );
                    break;
            }
            e++;
    }
}
#endif

const char * value_to_string (enum ADIOS_DATATYPES type, void * data, uint64_t element)
{
    static char s [100];
    s [0] = 0;

    switch (type)
    {
        case adios_unsigned_byte:
            //sprintf (s, "%u", *(((uint8_t *) data) + element));
            sprintf (s, "%u", *(((uint8_t *) &data) + element));
            break;

        case adios_byte:
            //sprintf (s, "%d", *(((int8_t *) data) + element));
            sprintf (s, "%d", *(((int8_t *) &data) + element));
            break;

        case adios_short:
            //sprintf (s, "%hd", *(((int8_t *) data) + element));
            sprintf (s, "%hd", *(((int8_t *) &data) + element));
            break;

        case adios_unsigned_short:
            //sprintf (s, "%uh", *(((int8_t *) data) + element));
            sprintf (s, "%uh", *(((int8_t *) &data) + element));
            break;

        case adios_integer:
            //sprintf (s, "%d", *(((int32_t *) data) + element));
            sprintf (s, "%d", *(((int32_t *) &data) + element));
            break;

        case adios_unsigned_integer:
            //sprintf (s, "%u", *(((uint32_t *) data) + element));
            sprintf (s, "%u", *(((uint32_t *) &data) + element));
            break;

        case adios_long:
            //sprintf (s, "%lld", *(((int64_t *) data) + element));
            sprintf (s, "%lld", *(((int64_t *) &data) + element));
            break;

        case adios_unsigned_long:
            //sprintf (s, "%llu", *(((uint64_t *) data) + element));
            sprintf (s, "%llu", *(((uint64_t *) &data) + element));
            break;

        case adios_real:
            //sprintf (s, "%e", *(((float *) data) + element));
            sprintf (s, "%e", *(((float *) &data) + element));
            break;

        case adios_double:
            //sprintf (s, "%le", *(((double *) data) + element));
            sprintf (s, "%le", *(((double *) &data) + element));
            break;

        case adios_long_double:
            //sprintf (s, "%Le", *(((long double *) data) + element));
            sprintf (s, "%Le", *(((long double *) &data) + element));
            break;

        case adios_string:
            //sprintf (s, "%s", ((char *) data) + element);
            sprintf (s, "%s", ((char *) data) + element);
            break;

        case adios_complex:
            sprintf (s, "(%f %f)", *((float *) data) + (element * 2 + 0)
                                 , *((float *) data) + (element * 2 + 1)
                    );
            break;

        case adios_double_complex:
            sprintf (s, "(%lf %lf)", *((double *) data) + (element * 2 + 0)
                                   , *((double *) data) + (element * 2 + 1)
                    );
            break;
    }

    return s;
}

        
void print_process_group_header (uint64_t num
                      ,struct adios_process_group_header_struct_v1 * pg_header
                      )
{
    int i;
    struct adios_method_info_struct_v1 * m;

    printf ("Process Group: %llu\n", num);
    printf ("\tGroup Name: %s\n", pg_header->name);
    printf ("\tHost Language Fortran?: %c\n"
           ,(pg_header->host_language_fortran == adios_flag_yes ? 'Y' : 'N')
           );
    printf ("\tCoordination Comm Member ID: %d\n", pg_header->coord_comm_id);
    printf ("\tCoordination Var Member ID: %d\n", pg_header->coord_var_id);
    printf ("\tTimestep Member ID: %d\n", pg_header->timestep_id);
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

void print_var_header (struct adios_var_header_struct_v1 * var_header)
{
    int i = 0;

    printf ("\t\tVar Name (ID): %s (%d)\n", var_header->name, var_header->id);
    printf ("\t\tVar Path: %s\n", var_header->path);
    printf ("\t\tDatatype: %s\n", adios_type_to_string (var_header->type));
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
                printf ("R(%llu):", d->dimension.rank);
            }
            else
            {
                printf ("V(%hu):", d->dimension.var_id);
            }
            if (d->global_dimension.var_id == 0)
            {
                printf ("R(%llu):", d->global_dimension.rank);
            }
            else
            {
                printf ("V(%hu):", d->global_dimension.var_id);
            }
            if (d->local_offset.var_id == 0)
            {
                printf ("R(%llu)\n", d->local_offset.rank);
            }
            else
            {
                printf ("V(%hu)\n", d->local_offset.var_id);
            }

            d = d->next;
	}
    }
}

void print_var_payload (struct adios_var_header_struct_v1 * var_header
                       ,struct adios_var_payload_struct_v1 * var_payload
                       ,struct dump_struct * dump
                       )
{
    printf ("\t\tMin: %s\n", value_to_string (var_header->type
                                             ,var_payload->min, 0
                                             )
           );
    printf ("\t\tMax: %s\n", value_to_string (var_header->type
                                             ,var_payload->max, 0
                                             )
           );
    if (dump->do_dump)
    {
        printf ("dump the data\n");
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
        printf ("\t\tDatatype: %s\n", adios_type_to_string (attribute->type));
        printf ("\t\tValue: %s\n", value_to_string (attribute->type
                                                   ,attribute->value, 0
                                                   )
               );
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
        printf ("\tTimestep: %d\n", pg_root->timestep);
        printf ("\tOffset in File: %llu\n", pg_root->offset_in_file);

        pg_root = pg_root->next;
    }
}

void print_vars_index (struct adios_index_var_struct_v1 * vars_root)
{
    while (vars_root)
    {
        if (!strcmp (vars_root->var_path, "/"))
        {
            printf ("Var (Group): /%s (%s)\n", vars_root->var_name
                   ,vars_root->group_name
                   );
        }
        else
	{
            printf ("Var (Group): %s/%s (%s)\n", vars_root->var_path
                   ,vars_root->var_name, vars_root->group_name
                   );
	}
        printf ("\tDatatype: %s\n", adios_type_to_string (vars_root->type));
        printf ("\tVars Entries: %llu\n", vars_root->entries_count);
        uint64_t i;
        printf ("\t\tOffset\t\tMin\t\tMax\n");
        for (i = 0; i < vars_root->entries_count; i++)
        {
            printf ("\t\t%s\t\t", value_to_string (adios_long
                                  ,(void *) (vars_root->entries [i].offset), 0)
                                  );
            printf ("%s\t\t", value_to_string (vars_root->type
                                           ,vars_root->entries [i].min, 0)
                                           );
            printf ("%s\n", value_to_string (vars_root->type
                                           ,vars_root->entries [i].max, 0)
                                           );
        }

        vars_root = vars_root->next;
    }
}
