#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/stat.h>
#include "adios_types.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"

void adios_buffer_struct_init (struct adios_bp_buffer_struct_v1 * b)
{
    b->f = -1;
    b->buff = 0;
    b->length = 0;
    b->change_endianness = adios_flag_unknown;
    b->offset = 0;
    b->end_of_pgs = 0;
    b->pg_index_offset = 0;
    b->pg_size = 0;
    b->vars_index_offset = 0;
    b->vars_size = 0;
    b->file_size = 0;
    b->read_pg_offset = 0;
    b->read_pg_size = 0;
}

void adios_buffer_struct_clear (struct adios_bp_buffer_struct_v1 * b)
{
    if (b->buff)
        free (b->buff);
    adios_buffer_struct_init (b);
}

//*****************************************************************************
// buff must be 4 bytes
int adios_parse_version (struct adios_bp_buffer_struct_v1 * b
                        ,uint32_t * version
                        )
{
    // if high bit set, big endian
    uint64_t test = 1;

    if (b->length < 4)
    {
        fprintf (stderr, "adios_parse_version requires a buffer of at least "
                         "4 bytes.  Only %lld were provided\n"
                ,b->length
                );

        return 1;
    }

    *version = ntohl (*(uint32_t *) (b->buff + b->offset));

    if (   (!*(char *) version && !*(char *) &test)
        || (*(char *) version && *(char *) &test)
       )
    {
        b->change_endianness = adios_flag_no;
    }
    else
    {
        b->change_endianness = adios_flag_yes;
    }

    *version = *version & 0x7fffffff;

    return 0;
}

// buff must be 16 bytes
int adios_parse_index_offsets_v1 (struct adios_bp_buffer_struct_v1 * b)
{
    if (b->length - b->offset < 16)
    {
        fprintf (stderr, "adios_parse_index_offsets_v1 requires a buffer of "
                         "at least 16 bytes.  Only %lld were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    uint64_t vars_end = b->file_size - 20;

    b->pg_index_offset = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;
    b->vars_index_offset = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    b->end_of_pgs = b->pg_index_offset;
    b->pg_size = b->vars_index_offset - b->pg_index_offset;
    b->vars_size = vars_end - b->vars_index_offset;

    return 0;
}

int adios_parse_process_group_index_v1 (struct adios_bp_buffer_struct_v1 * b
                         ,struct adios_index_process_group_struct_v1 ** pg_root
                         )
{
    struct adios_index_process_group_struct_v1 ** root;
    if (b->length - b->offset < 16)
    {
        fprintf (stderr, "adios_parse_process_group_index_v1 requires a buffer "
                         "of at least 16 bytes.  Only %lld were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    root = pg_root;

    uint64_t process_groups_count;
    uint64_t process_groups_length;

    process_groups_count = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;
    process_groups_length = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    // validate remaining length

    uint64_t i;
    for (i = 0; i < process_groups_count; i++)
    {
        uint16_t length_of_group;
        // validate remaining length
        length_of_group = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;

        if (!*root)
        {
            *root = (struct adios_index_process_group_struct_v1 *)
                   malloc (sizeof(struct adios_index_process_group_struct_v1));
            (*root)->next = 0;
        }
        uint16_t length_of_name;

        length_of_name = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->group_name = (char *) malloc (length_of_name + 1);
        (*root)->group_name [length_of_name] = '\0';
        memcpy ((*root)->group_name, b->buff + b->offset, length_of_name);
        b->offset += length_of_name;

        (*root)->process_id = *(uint32_t *) (b->buff + b->offset);
        b->offset += 4;

        (*root)->timestep = *(uint32_t *) (b->buff + b->offset);
        b->offset += 4;

        (*root)->offset_in_file = *(uint64_t *) (b->buff + b->offset);
        b->offset += 8;

        root = &(*root)->next;
    }

    return 0;
}

int adios_parse_vars_index_v1 (struct adios_bp_buffer_struct_v1 * b
                              ,struct adios_index_var_struct_v1 ** vars_root
                              )
{
    struct adios_index_var_struct_v1 ** root;

    if (b->length - b->offset < 10)
    {
        fprintf (stderr, "adios_parse_vars_index_v1 requires a buffer "
                         "of at least 10 bytes.  Only %lld were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    root = vars_root;

    uint16_t vars_count;
    uint64_t vars_length;

    vars_count = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    vars_length = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    // validate remaining length

    int i;
    for (i = 0; i < vars_count; i++)
    {
        if (!*root)
        {
            *root = (struct adios_index_var_struct_v1 *)
                          malloc (sizeof (struct adios_index_var_struct_v1));
            (*root)->next = 0;
        }
        uint32_t var_entry_length;
        uint16_t len;
        uint64_t offsets_count;

        var_entry_length = *(uint32_t *) (b->buff + b->offset);
        b->offset += 4;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->group_name = (char *) malloc (len + 1);
        (*root)->group_name [len] = '\0';
        strncpy ((*root)->group_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->var_name = (char *) malloc (len + 1);
        (*root)->var_name [len] = '\0';
        strncpy ((*root)->var_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->var_path = (char *) malloc (len + 1);
        (*root)->var_path [len] = '\0';
        strncpy ((*root)->var_path, b->buff + b->offset, len);
        b->offset += len;

        (*root)->type = (enum ADIOS_DATATYPES) *(b->buff + b->offset);
        b->offset += 1;

        offsets_count = *(uint64_t *) (b->buff + b->offset);
        (*root)->entries_count = offsets_count;
        (*root)->entries_allocated = offsets_count;
        b->offset += 8;

        // validate remaining length: offsets_count * (8 + 2 * (size of type))
        uint64_t j;
        (*root)->entries = malloc (offsets_count
                            * sizeof (struct adios_index_var_entry_struct_v1));
        for (j = 0; j < offsets_count; j++)
        {
            (*root)->entries [j].offset = *(uint64_t *) (b->buff + b->offset);
            b->offset += 8;
            char * str = "";
            void * x = adios_dupe_data_scalar (adios_string, str);
//                             ((*root)->type, (void *) (b->buff + b->offset));
            (*root)->entries [j].min = adios_dupe_data_scalar
                                         ((*root)->type, b->buff + b->offset);
            b->offset += adios_get_type_size ((*root)->type, "");
            (*root)->entries [j].max = adios_dupe_data_scalar
                                         ((*root)->type, b->buff + b->offset);
            b->offset += adios_get_type_size ((*root)->type, "");
        }

        root = &(*root)->next;
    }

    return 0;
}

int adios_parse_process_group_header_v1 (struct adios_bp_buffer_struct_v1 * b
                       ,struct adios_process_group_header_struct_v1 * pg_header
                            )
{
    if (b->length - b->offset < 16)
    {
        fprintf (stderr, "adios_parse_process_group_header_v1 requires a "
                         "buffer of at least 16 bytes.  "
                         "Only %lld were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    uint64_t size;
    size = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    pg_header->host_language_fortran =
                             (*(b->buff + b->offset) == 'y' ? adios_flag_yes
                                                            : adios_flag_no
                             );
    b->offset += 1;

    uint16_t len;
    uint8_t count;
    len = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    pg_header->name = (char *) malloc (len + 1);
    pg_header->name [len] = '\0';
    memcpy (pg_header->name, b->buff + b->offset, len);
    b->offset += len;

    pg_header->coord_comm_id = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    pg_header->coord_var_id = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    pg_header->timestep_id = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;

    pg_header->methods_count = *(b->buff + b->offset);
    b->offset += 1;

    // length of methods section
    len = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;

    int i;
    struct adios_method_info_struct_v1 ** root;
    
    pg_header->methods = 0;
    root = &pg_header->methods;
    for (i = 0; i < pg_header->methods_count; i++)
    {
        if (!*root)
        {
            *root = (struct adios_method_info_struct_v1 *)
                        malloc (sizeof (struct adios_method_info_struct_v1));
            (*root)->next = 0;
        }

        (*root)->id = (enum ADIOS_IO_METHOD) *(b->buff + b->offset);
        b->offset += 1;

        len = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        (*root)->parameters = (char *) malloc (len + 1);
        (*root)->parameters [len] = '\0';
        strncpy ((*root)->parameters, b->buff + b->offset, len);
        b->offset += len;

        root = &(*root)->next;
    }

    return 0;
}

int adios_parse_vars_header_v1 (struct adios_bp_buffer_struct_v1 * b
                               ,struct adios_vars_header_struct_v1 * vars_header
                               )
{
    if (b->length - b->offset < 12)
    {
        fprintf (stderr, "adios_parse_var_header_v1 requires a "
                         "buffer of at least 12 bytes.  "
                         "Only %lld were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    vars_header->count = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    vars_header->length = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    return 0;
}

int adios_parse_var_data_header_v1 (struct adios_bp_buffer_struct_v1 * b
                               ,struct adios_var_header_struct_v1 * var_header
                               )
{
    if (b->length - b->offset < 21)
    {
        fprintf (stderr, "adios_parse_var_data_header_v1 requires a "
                         "buffer of at least 21 bytes.  "
                         "Only %lld were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    uint64_t initial_offset = b->offset;  // save to calc payload size
    uint64_t length_of_var;
    uint16_t len;

    length_of_var = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    //validate remaining length

    var_header->id = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    
    len = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    var_header->name = (char *) malloc (len + 1);
    var_header->name [len] = '\0';
    memcpy (var_header->name, b->buff + b->offset, len);
    b->offset += len;

    len = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    var_header->path = (char *) malloc (len + 1);
    var_header->path [len] = '\0';
    memcpy (var_header->path, b->buff + b->offset, len);
    b->offset += len;

    var_header->type = (enum ADIOS_DATATYPES) *(b->buff + b->offset);
    b->offset += 1;

    var_header->is_dim = (*(b->buff + b->offset) == 'y' ? adios_flag_yes
                                                        : adios_flag_no
                         );
    b->offset += 1;

    int i;
    uint8_t dims_count;
    uint16_t dims_length;

    // validate remaining length

    dims_count = *(b->buff + b->offset);
    b->offset += 1;
    dims_length = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;

    var_header->dims = 0;
    struct adios_dimension_struct_v1 ** root = &var_header->dims;

    for (i = 0; i < dims_count; i++)
    {
        uint8_t flag;

        if (!*root)
        {
            *root = (struct adios_dimension_struct_v1 *)
                           malloc (sizeof (struct adios_dimension_struct_v1));
            (*root)->next = 0;
        }

        flag = *(b->buff + b->offset);
        b->offset += 1;
        if (flag == 'y')
        {
            (*root)->dimension.rank = 0;
            (*root)->dimension.var_id = *(uint16_t *) (b->buff + b->offset);
            b->offset += 2;
        }
        else
        {
            (*root)->dimension.rank = *(uint64_t *) (b->buff + b->offset);
            (*root)->dimension.var_id = 0;
            b->offset += 8;
        }

        flag = *(b->buff + b->offset);
        b->offset += 1;
        if (flag == 'y')
        {
            (*root)->global_dimension.rank = 0;
            (*root)->global_dimension.var_id = *(uint16_t *)
                                                         (b->buff + b->offset);
            b->offset += 2;
        }
        else
        {
            (*root)->global_dimension.rank = *(uint64_t *)
                                                         (b->buff + b->offset);
            (*root)->global_dimension.var_id = 0;
            b->offset += 8;
        }

        flag = *(b->buff + b->offset);
        b->offset += 1;
        if (flag == 'y')
        {
            (*root)->local_offset.rank = 0;
            (*root)->local_offset.var_id = *(uint16_t *) (b->buff + b->offset);
            b->offset += 2;
        }
        else
        {
            (*root)->local_offset.rank = *(uint64_t *) (b->buff + b->offset);
            (*root)->local_offset.var_id = 0;
            b->offset += 8;
        }
    }

    var_header->payload_size = length_of_var - (b->offset - initial_offset);

    return 0;
}

int adios_parse_var_data_payload_v1 (struct adios_bp_buffer_struct_v1 * b
                             ,struct adios_var_header_struct_v1 * var_header
                             ,struct adios_var_payload_struct_v1 * var_payload
                             )
{
    if (b->length - b->offset < var_header->payload_size)
    {
        fprintf (stderr, "adios_parse_var_data_payload_v1 requires a "
                         "buffer of at least %lld bytes.  "
                         "Only %lld were provided\n"
                ,var_header->payload_size, b->length - b->offset
                );

        return 1;
    }

    uint64_t size = adios_get_type_size (var_header->type, "");

    if (var_payload)
    {
        memcpy (&var_payload->min, (b->buff + b->offset), size);
        b->offset += size;

        memcpy (&var_payload->max, (b->buff + b->offset), size);
        b->offset += size;

        if (var_payload->payload)
        {
            memcpy (var_payload->payload, (b->buff + b->offset)
                   ,var_header->payload_size - (2 * size)
                   );
            b->offset += var_header->payload_size - (2 * size);
        }
        else
        {
            b->offset += var_header->payload_size - (2 * size);
        }
    }
    else
    {
        b->offset += var_header->payload_size;
    }

    return 0;
}

int adios_parse_attributes_header_v1 (struct adios_bp_buffer_struct_v1 * b
                      ,struct adios_attributes_header_struct_v1 * attrs_header
                      )
{
    if (b->length - b->offset < 10)
    {
        fprintf (stderr, "adios_parse_attribute_header_v1 requires a "
                         "buffer of at least 10 bytes.  "
                         "Only %lld were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    attrs_header->count = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    attrs_header->length = *(uint64_t *) (b->buff + b->offset);
    b->offset += 8;

    return 0;
}

int adios_parse_attribute_v1 (struct adios_bp_buffer_struct_v1 * b
                             ,struct adios_attribute_struct_v1 * attribute
                             )
{
    if (b->length - b->offset < 15)
    {
        fprintf (stderr, "adios_parse_attribute_data_payload_v1 requires a "
                         "buffer of at least 15 bytes.  "
                         "Only %lld were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    uint32_t attribute_length;
    uint16_t len;

    attribute_length = *(uint32_t *) (b->buff + b->offset);
    b->offset += 4;
    attribute->id = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;

    len = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    attribute->name = (char *) malloc (len + 1);
    attribute->name [len] = '\0';
    strncpy (attribute->name, b->buff + b->offset, len);
    b->offset += len;

    len = *(uint16_t *) (b->buff + b->offset);
    b->offset += 2;
    attribute->path = (char *) malloc (len + 1);
    attribute->path [len] = '\0';
    strncpy (attribute->path, b->buff + b->offset, len);
    b->offset += len;

    attribute->is_var = (*(b->buff + b->offset) == 'y' ? adios_flag_yes
                                                       : adios_flag_no
                        );
    b->offset += 1;
    if (attribute->is_var == adios_flag_yes)
    {
        attribute->var_id = *(uint16_t *) (b->buff + b->offset);
        b->offset += 2;
        attribute->type = adios_unknown;
        attribute->length = 0;
        attribute->value = 0;
    }
    else
    {
        attribute->var_id = 0;
        attribute->type = (enum ADIOS_DATATYPES) *(b->buff + b->offset);
        b->offset += 1;
        attribute->length = *(uint32_t *) (b->buff + b->offset);
        b->offset += 4;

        attribute->value = malloc (attribute->length + 1);
        ((char *) attribute->value) [attribute->length] = '\0';
        memcpy (attribute->value, (b->buff + b->offset), attribute->length);
        b->offset += attribute->length;
    }

    return 0;
}

void * adios_dupe_data_scalar (enum ADIOS_DATATYPES type, void * in)
{
    void * out;
    int element_size = adios_get_type_size (type, in);

    void * d;

    switch (type)
    {
        case adios_byte:
        case adios_short:
        case adios_integer:
        case adios_long:
        case adios_unsigned_byte:
        case adios_unsigned_short:
        case adios_unsigned_integer:
        case adios_unsigned_long:
        case adios_real:
        case adios_double:
        case adios_long_double:
        case adios_complex:
        case adios_double_complex:
            d = malloc (element_size);
            if (!d)
            {
                fprintf (stderr, "cannot allocate %d bytes to copy scalar\n"
                        ,element_size
                        );

                return 0;
            }

            memcpy ((char *) d, (char *) in, element_size);
            break;

        case adios_string:
            d = malloc (element_size + 1);
            if (!d)
            {
                fprintf (stderr, "cannot allocate %d bytes to copy scalar\n"
                        ,element_size + 1
                        );

                return 0;
            }

            memcpy ((char *) d, (char *) in, element_size + 1);
            break;

        default:
            d = 0;
            break;
    }

    return d;
}

//*****************************************************************************
void adios_init_buffer_read_version (struct adios_bp_buffer_struct_v1 * b)
{
    if (!b->buff)
    {
        b->buff = malloc (20); // all we need for version and indexes
        if (!b->buff)
            fprintf(stderr, "could not allocate 20 bytes\n");
        b->length = 20;
        b->offset = 16;
    }
}

// last 4 bytes of file
void adios_posix_read_version (struct adios_bp_buffer_struct_v1 * b)
{
    uint64_t buffer_size;
    uint64_t start;
    uint64_t r;

    adios_init_buffer_read_version (b);

    lseek64 (b->f, b->file_size - 20, SEEK_SET);

    r = read (b->f, b->buff, 20);
    if (r != 20)
        fprintf (stderr, "could not read 20 bytes. read only: %lld\n", r);
}

void adios_posix_read_index_offsets (struct adios_bp_buffer_struct_v1 * b)
{
    b->offset = 0; // just move to the start of the buffer
}

void adios_init_buffer_read_process_group_index (
                                          struct adios_bp_buffer_struct_v1 * b
                                          )
{
    b->buff = realloc (b->buff, b->pg_size);
    b->offset = 0;
}

void adios_posix_read_process_group_index (struct adios_bp_buffer_struct_v1 * b)
{
    adios_init_buffer_read_process_group_index (b);

    lseek (b->f, b->pg_index_offset, SEEK_SET);
    read (b->f, b->buff, b->pg_size);
}

void adios_init_buffer_read_vars_index (struct adios_bp_buffer_struct_v1 * b)
{
    b->buff = realloc (b->buff, b->vars_size);
    b->offset = 0;
}

void adios_posix_read_vars_index (struct adios_bp_buffer_struct_v1 * b)
{
    adios_init_buffer_read_vars_index (b);

    lseek (b->f, b->vars_index_offset, SEEK_SET);
    uint64_t r = read (b->f, b->buff, b->vars_size);
    if (r != b->vars_size)
        fprintf (stderr, "reading vars_index: wanted %llu, read: %llu\n"
                ,b->vars_size, r
                );
}

void adios_init_buffer_read_procss_group (struct adios_bp_buffer_struct_v1 * b
                          ,uint64_t pg_number
                          ,struct adios_index_process_group_struct_v1 * pg_root
                          )
{
    while (pg_root && --pg_number > 0)
    {
        pg_root = pg_root->next;
    }
    if (pg_root)
    {
        uint64_t size = 0;

        if (pg_root->next)
        {
            size = pg_root->next->offset_in_file - pg_root->offset_in_file;
        }
        else
        {
            size = b->pg_index_offset - pg_root->offset_in_file;
        }

        b->buff = realloc (b->buff, size);
        b->length = size;
        b->offset = 0;
        b->read_pg_offset = pg_root->offset_in_file;
        b->read_pg_size = size;
    }
    else
    {
        fprintf (stderr, "invalid process group number %llu\n", pg_number);
    }
}

uint64_t adios_posix_read_process_group (struct adios_bp_buffer_struct_v1 * b
                          ,uint64_t pg_number
                          ,struct adios_index_process_group_struct_v1 * pg_root
                          )
{
    uint64_t pg_size = 0;

    if (!pg_root)
    {
        fprintf (stderr, "just start at the front and go from there (save "
                         "position for faster reading)\n"
                );
    }
    else
    {
        adios_init_buffer_read_process_group (b, pg_number, pg_root);
        lseek (b->f, b->read_pg_offset, SEEK_SET);
        pg_size = read (b->f, b->buff, b->read_pg_size);
        if (pg_size != b->read_pg_size)
        {
            fprintf (stderr, "adios_read_process_group: "
                             "Tried to read: %llu, but only got: %llu\n"
                    ,b->read_pg_size, pg_size
                    );

            pg_size = 0;
        }
    }

    return pg_size;
}

int adios_posix_open_read_internal (const char * filename
                                   ,const char * base_path
                                   ,struct adios_bp_buffer_struct_v1 * b
                                   )
{
    char name [STR_LEN];
    struct stat s;

    sprintf (name, "%s%s", base_path, filename);

    if (stat (name, &s) == 0)
        b->file_size = s.st_size;

    b->f = open64 (name, O_RDONLY);
    if (b->f == -1)
    {
        fprintf (stderr, "ADIOS POSIX: file not found: %s\n", name);

        return 0;
    }

    return 1;
}

void adios_posix_close_internal (struct adios_bp_buffer_struct_v1 * b)
{
    if (b->f != -1)
    {
        close (b->f);
    }

    b->f = -1;
    adios_buffer_struct_clear (b);
}
