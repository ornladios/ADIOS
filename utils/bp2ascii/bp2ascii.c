/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "adios_internals.h"
#include "adios_bp_v1.h"
#include "adios_read.h"

struct dump_struct
{
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
void print_process_group_index (
                         struct adios_index_process_group_struct_v1 * pg_root
                         );
void print_var_header (struct adios_var_header_struct_v1 * var_header);

int main (int argc, char ** argv)
{
    char * filename;
    char newfilename [256];
    int i = 0;
    int rc = 0;
    struct dump_struct dump;

    if (argv [1][0] && argv [1][0] == '-' && argc > 3)
    {
        if (   !strcmp (argv [1], "-v")
            || !strcmp (argv [1], "--var")
           )
        {
            dump.dump_var = argv [2];
            filename = argv [3];
            printf("%s %s\n",dump.dump_var,filename);
	    if (argc > 4)
	        strcpy (newfilename, argv[4]);
	    else
            {
		strcpy (newfilename, dump.dump_var);
		strcat (newfilename, ".dat");
            }
        }
        else
        {
                fprintf (stderr, "usage: %s [-v [var]|--var [var]] <filename> [newfilename]\n"
                        ,argv [0]
                        );
                return -1;
        }
    }    
    else
    {
        fprintf (stderr, "usage: %s [-v [var]|--var [var]] <filename> [newfilename]\n"
                ,argv [0]
                );
        return -1;
    }

    struct adios_bp_buffer_struct_v1 * b = 0;
    uint32_t version = 0;
    FILE * outf;
    outf = fopen (newfilename, "w");
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
    struct adios_index_var_struct_v1 * vars_root = 0;
    //struct adios_index_attribute_struct_v1 * attrs_root = 0;

    adios_posix_read_index_offsets (b);
    adios_parse_index_offsets_v1 (b);

    adios_posix_read_process_group_index (b);
    adios_parse_process_group_index_v1 (b, &pg_root, NULL);
    print_process_group_index (pg_root);

    adios_posix_read_vars_index (b);
    adios_parse_vars_index_v1 (b, &vars_root, NULL, NULL);

    while (vars_root)
    {
        //printf("%s\n",vars_root->var_name);
        if (!strcasecmp (dump.dump_var, vars_root->var_name) )
            break;
        else
            vars_root = vars_root->next;
    }

    uint64_t element_counts = vars_root->characteristics_count;
    //struct adios_process_group_header_struct_v1 pg_header;
    struct adios_var_header_struct_v1 var_header;
    struct adios_var_payload_struct_v1 var_payload;
    uint64_t offset, var_len;
    printf("characteristics count: %" PRIu64 "\n", element_counts);

    for (i = 0; i < element_counts; i++)
    {
        offset = vars_root->characteristics[i].offset;
        b->read_pg_offset = pg_root->offset_in_file;

        printf("offset: %" PRIu64 " read_pg_offset=%" PRIu64 "\n", offset, b->read_pg_offset);
        if (pg_root->next)
        {
            b->read_pg_size =   pg_root->next->offset_in_file
                              - pg_root->offset_in_file;
        }
        else
        {
            b->read_pg_size =   b->pg_index_offset
                              - pg_root->offset_in_file;
        }
        
        adios_init_buffer_read_process_group (b);
        lseek(b->f, b->read_pg_offset+offset,SEEK_SET);

        read(b->f,b->buff,8);
        var_len = *(uint64_t *) b->buff;

        read (b->f,b->buff+8, var_len);

        printf("var length: %" PRIu64 " offset: %" PRIu64 "\n",var_len, offset);
        printf ("payload_size %" PRIu64 "x\n", var_header.payload_size);

        adios_parse_var_data_header_v1 (b, &var_header);
        print_var_header (&var_header);
        
        var_payload.payload = malloc (var_header.payload_size);   
	adios_parse_var_data_payload_v1 (b, &var_header, &var_payload
                                        ,var_header.payload_size
                                        );
        int j;
        switch (var_header.type) {
            case adios_long_double:
                for (j=0; j<var_header.payload_size/16;j++) 
                    fprintf(outf, "%Lf ", *((long double*)var_payload.payload+j));
                    break;
            case adios_double:
                for (j=0; j<var_header.payload_size/8;j++) 
                    fprintf(outf, "%f ", *((double*)var_payload.payload+j));
                    break;
            case adios_real:
                for (j=0; j<(var_header.payload_size)/4;j++) 
                    fprintf(outf, "%f ", *((float *)var_payload.payload+j));
                    break;
            case adios_long:
                for (j=0; j<(var_header.payload_size)/8;j++) 
                    fprintf(outf, "%ld ", *((long *)var_payload.payload+j));
                    break;
            case adios_unsigned_long:
                for (j=0; j<(var_header.payload_size)/8;j++) 
                    fprintf(outf, "%lu ", *((unsigned long*)var_payload.payload+j));
                    break;
            case adios_integer:
                for (j=0; j<(var_header.payload_size)/4;j++) 
                    fprintf(outf, "%d ", *((int *)var_payload.payload+j));
                    break;
            case adios_unsigned_integer:
                for (j=0; j<(var_header.payload_size)/4;j++) 
                    fprintf(outf, "%u ", *((unsigned int *)var_payload.payload+j));
                    break;
            case adios_unsigned_short:
                for (j=0; j<(var_header.payload_size)/2;j++) 
                    fprintf(outf, "%d ", *((unsigned short*)var_payload.payload+j));
                    break;
            case adios_short:
                for (j=0; j<(var_header.payload_size)/2;j++) 
                    fprintf(outf, "%u ", *((short *)var_payload.payload+j));
                    break;
            case adios_unsigned_byte:
                for (j=0; j<var_header.payload_size;j++) 
                    fprintf(outf, "%c ", *((unsigned char *)var_payload.payload+j));
                    break;
            case adios_byte:
                for (j=0; j<var_header.payload_size;j++) 
                    fprintf(outf, "%c ", *((char *)var_payload.payload+j));
                    break;
            case adios_string:
                    fprintf(outf, "%s ", (char *)var_payload.payload);
                    break;
            case adios_string_array:
                    fprintf(outf, "%s ", *(char **)var_payload.payload);
                    break;
            case adios_complex:
                for (j=0; j<(var_header.payload_size)/8;j++) 
                    fprintf(outf, "%f + %fi", *((float  *)var_payload.payload+2*j),*((float  *)var_payload.payload+j*2+1));
                    break;
            case adios_double_complex:
                for (j=0; j<(var_header.payload_size)/16;j++) 
                    fprintf(outf, "%f + %fi", *((double *)var_payload.payload+2*j),*((double *)var_payload.payload+j*2+1));
                    break;
            case adios_unknown:
                    break;
         }        
         fprintf(outf,"\n");
    }
    fclose(outf);
    adios_posix_close_internal (b);
    return 0;
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
void print_var_header (struct adios_var_header_struct_v1 * var_header)
{
    printf ("\t\tVar Name (ID): %s (%d)\n", var_header->name, var_header->id);
    printf ("\t\tVar Path: %s\n", var_header->path);
    printf ("\t\tDatatype: %s\n", adios_type_to_string_int (var_header->type));
    printf ("\t\tIs Dimension: %c\n"
           ,(var_header->is_dim == adios_flag_yes ? 'Y' : 'N')
           );
}
