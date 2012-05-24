/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "adios_types.h"
#include "adios_internals.h"
#include "adios_bp_v1.h"
#include "adios_endianness.h"

#if defined(__APPLE__)
#    define O_LARGEFILE 0
#endif

#define BYTE_ALIGN 8

#ifdef DMALLOC
#include "dmalloc.h"
#endif

static void alloc_aligned (struct adios_bp_buffer_struct_v1 * b, uint64_t size)
{
    b->allocated_buff_ptr = malloc (size + BYTE_ALIGN - 1);
    if (!b->allocated_buff_ptr)
    {
        fprintf (stderr, "Cannot allocate: %llu\n", size);

        b->buff = 0;
        b->length = 0;

        return;
    }
    uint64_t p = (uint64_t) b->allocated_buff_ptr;
    b->buff = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
    b->length = size;
}

static void realloc_aligned (struct adios_bp_buffer_struct_v1 * b
                            ,uint64_t size
                            )
{
    b->allocated_buff_ptr = realloc (b->allocated_buff_ptr
                                    ,size + BYTE_ALIGN - 1
                                    );
    if (!b->allocated_buff_ptr)
    {
        fprintf (stderr, "Cannot allocate: %llu\n", size);

        b->buff = 0;
        b->length = 0;

        return;
    }
    uint64_t p = (uint64_t) b->allocated_buff_ptr;
    b->buff = (char *) ((p + BYTE_ALIGN - 1) & ~(BYTE_ALIGN - 1));
    b->length = size;
}

void adios_shared_buffer_free (struct adios_bp_buffer_struct_v1 * b)
{
    if (b->allocated_buff_ptr)
        free (b->allocated_buff_ptr);
    b->allocated_buff_ptr = 0;
    b->buff = 0;
    b->offset = 0;
    b->length = 0;
}

void adios_buffer_struct_init (struct adios_bp_buffer_struct_v1 * b)
{
    b->f = -1;
    b->allocated_buff_ptr = 0;
    b->buff = 0;
    b->length = 0;
    b->change_endianness = adios_flag_unknown;
    b->version = 0;
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
    if (b->allocated_buff_ptr)
        free (b->allocated_buff_ptr);
    adios_buffer_struct_init (b);
}

// *****************************************************************************
// buff must be 4 bytes
int adios_parse_version (struct adios_bp_buffer_struct_v1 * b
                        ,uint32_t * version
                        )
{
    // if high bit set, big endian
    uint32_t test = 1;

    if (b->length < 4)
    {
        fprintf (stderr, "adios_parse_version requires a buffer of at least "
                         "4 bytes.  Only %llu were provided\n"
                ,b->length
                );

        return 1;
    }

    *version = ntohl (*(uint32_t *) (b->buff + b->offset));
    char *v = (char *) version;
    if (   (*v && !*(char *) &test)       // both writer and this machine are big endian
        || (!*(v+3) && *(char *) &test)   // both are little endian
       )
    {
        b->change_endianness = adios_flag_no;
        //fprintf(stderr, "no need to change endianness\n");
    }
    else
    {
        //fprintf(stderr, "change endianness\n");
        b->change_endianness = adios_flag_yes;
    }

    *version = *version & 0x7fffffff;

    return 0;
}

// buff must be 16 bytes
int adios_parse_index_offsets_v1 (struct adios_bp_buffer_struct_v1 * b)
{
    if (b->length - b->offset < 24)
    {
        fprintf (stderr, "adios_parse_index_offsets_v1 requires a buffer of "
                         "at least 24 bytes.  Only %llu were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    uint64_t attrs_end = b->file_size - 28;
    int i;

    char * t;

    t = b->buff + b->offset;
    b->pg_index_offset = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(b->pg_index_offset);
    }

    b->offset += 8;

    b->vars_index_offset = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(b->vars_index_offset);
    }

    t = b->buff + b->offset;
    b->offset += 8;

    b->attrs_index_offset = *(uint64_t *) (b->buff + b->offset);
    t = b->buff + b->offset;
    if(b->change_endianness == adios_flag_yes) {
        swap_64(b->attrs_index_offset);
    }

    b->offset += 8;

    b->end_of_pgs = b->pg_index_offset;
    b->pg_size = b->vars_index_offset - b->pg_index_offset;
    b->vars_size = b->attrs_index_offset - b->vars_index_offset;
    b->attrs_size = attrs_end - b->attrs_index_offset;

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
                         "of at least 16 bytes.  Only %llu were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    root = pg_root;

    uint64_t process_groups_count;
    uint64_t process_groups_length;

    process_groups_count = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(process_groups_count);
    }
    b->offset += 8;

    process_groups_length = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(process_groups_length);
    }
    b->offset += 8;

    // validate remaining length

    uint64_t i;
    for (i = 0; i < process_groups_count; i++)
    {
        uint16_t length_of_group;
        // validate remaining length
        length_of_group = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(length_of_group);
        }
        b->offset += 2;

        if (!*root)
        {
            *root = (struct adios_index_process_group_struct_v1 *)
                   malloc (sizeof(struct adios_index_process_group_struct_v1));
            (*root)->next = 0;
        }
        uint16_t length_of_name;

        length_of_name = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(length_of_name);
        }
        b->offset += 2;
        (*root)->group_name = (char *) malloc (length_of_name + 1);
        (*root)->group_name [length_of_name] = '\0';
        memcpy ((*root)->group_name, b->buff + b->offset, length_of_name);
        b->offset += length_of_name;

        (*root)->adios_host_language_fortran =
                             (*(b->buff + b->offset) == 'y' ? adios_flag_yes
                                                            : adios_flag_no
                             );
        b->offset += 1;

        (*root)->process_id = *(uint32_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_32((*root)->process_id);
        }
        b->offset += 4;

        length_of_name = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(length_of_name);
        }
        b->offset += 2;
        (*root)->time_index_name = (char *) malloc (length_of_name + 1);
        (*root)->time_index_name [length_of_name] = '\0';
        memcpy ((*root)->time_index_name, b->buff + b->offset, length_of_name);
        b->offset += length_of_name;

        (*root)->time_index = *(uint32_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_32((*root)->time_index);
        }
        b->offset += 4;

        (*root)->offset_in_file = *(uint64_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_64((*root)->offset_in_file);
        }
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
                         "of at least 10 bytes.  Only %llu were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    root = vars_root;

    uint16_t vars_count;
    uint64_t vars_length;

    vars_count = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(vars_count);
    }
    b->offset += 2;

    vars_length = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(vars_length);
    }
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
        uint8_t flag;
        uint32_t var_entry_length;
        uint16_t len;
        uint64_t characteristics_sets_count;
        uint64_t type_size;

        var_entry_length = *(uint32_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_32(var_entry_length);
        }
        b->offset += 4;

        (*root)->id = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16((*root)->id);
        }
        b->offset += 2;

        len = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(len);
        }
        b->offset += 2;
        (*root)->group_name = (char *) malloc (len + 1);
        (*root)->group_name [len] = '\0';
        strncpy ((*root)->group_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(len);
        }
        b->offset += 2;
        (*root)->var_name = (char *) malloc (len + 1);
        (*root)->var_name [len] = '\0';
        strncpy ((*root)->var_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(len);
        }
        b->offset += 2;
        (*root)->var_path = (char *) malloc (len + 1);
        (*root)->var_path [len] = '\0';
        strncpy ((*root)->var_path, b->buff + b->offset, len);
        b->offset += len;

        flag = *(b->buff + b->offset);
        (*root)->type = (enum ADIOS_DATATYPES) flag;
        b->offset += 1;
        type_size = adios_get_type_size ((*root)->type, "");

        characteristics_sets_count = *(uint64_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_64(characteristics_sets_count);
        }
        (*root)->characteristics_count = characteristics_sets_count;
        (*root)->characteristics_allocated = characteristics_sets_count;
        b->offset += 8;

        // validate remaining length: offsets_count * (8 + 2 * (size of type))
        uint64_t j;
        (*root)->characteristics = malloc (characteristics_sets_count
                         * sizeof (struct adios_index_characteristic_struct_v1)
                        );
        memset ((*root)->characteristics, 0
               ,  characteristics_sets_count
                * sizeof (struct adios_index_characteristic_struct_v1)
               );
        for (j = 0; j < characteristics_sets_count; j++)
        {
            uint8_t characteristic_set_count;
            uint32_t characteristic_set_length;
            uint8_t item = 0;

            (*root)->characteristics [j].stats = 0;
            characteristic_set_count = (uint8_t) *(b->buff + b->offset);
            b->offset += 1;

            characteristic_set_length = *(uint32_t *) (b->buff + b->offset);
            if(b->change_endianness == adios_flag_yes) {
                swap_32(characteristic_set_length);
            }
            b->offset += 4;

            while (item < characteristic_set_count)
            {
                uint8_t flag;
                enum ADIOS_CHARACTERISTICS c;
                flag = *(b->buff + b->offset);
                c = (enum ADIOS_CHARACTERISTICS) flag;
                b->offset += 1;

                switch (c)
                {
                    case adios_characteristic_min:
                    case adios_characteristic_max:
                    case adios_characteristic_value:
                    {
                        uint16_t data_size;
                        void * data = 0;

                        if ((*root)->type == adios_string)
                        {
                            data_size = *(uint16_t *) (b->buff + b->offset);
                            if(b->change_endianness == adios_flag_yes) {
                                swap_16(data_size);
                            }
                            b->offset += 2;
                        }
                        else
                        {
                            data_size = adios_get_type_size ((*root)->type, "");
                        }

                        switch ((*root)->type)
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
                                data = malloc (data_size);

                                if (!data)
                                {
                                    fprintf (stderr, "cannot allocate %d bytes "
                                                     "to copy scalar %s\n"
                                            ,data_size
                                            ,(*root)->var_name
                                            );

                                    return 1;
                                }

                                memcpy (data, (b->buff + b->offset), data_size);
                                if(b->change_endianness == adios_flag_yes) {
                                    if((*root)->type == adios_complex) {
                                        // TODO
                                    }
                                    else if((*root)->type == adios_double_complex) {
                                        // TODO
                                    }
                                    else {
                                        switch(data_size)
                                        {
                                            case 2:
                                                swap_16_ptr(data);
                                                break;
                                            case 4:
                                                swap_32_ptr(data);
                                                break;
                                            case 8:
                                                swap_64_ptr(data);
                                                break;
                                            case 16:
                                                swap_128_ptr(data);
                                                break;
                                       }
                                    }
                                }
                                b->offset += data_size;
                                break;

                            case adios_string:
                                data = malloc (data_size + 1);

                                if (!data)
                                {
                                    fprintf (stderr, "cannot allocate %d bytes "
                                                     "to copy scalar %s\n"
                                            ,data_size
                                            ,(*root)->var_name
                                            );

                                    return 1;
                                }

                                ((char *) data) [data_size] = '\0';
                                memcpy (data, (b->buff + b->offset), data_size);
                                b->offset += data_size;
                                break;

                            default:
                                data = 0;
                                break;
                        }

                        switch (c)
                        {
                            case adios_characteristic_value:
                                (*root)->characteristics [j].value = data;
                                break;

                            // NCSU - reading older bp files
                            // adios_characteristic_min, max are not used anymore. If this is encountered it is an older bp file format
                            // Code below reads min and min, and sets the bitmap for those 2 alone
                            case adios_characteristic_min:
                                if (!(*root)->characteristics [j].stats)
                                {
                                    (*root)->characteristics [j].stats = malloc (sizeof(struct adios_index_characteristics_stat_struct *));
                                    (*root)->characteristics [j].stats[0] = malloc (2 * sizeof(struct adios_index_characteristics_stat_struct));
                                    (*root)->characteristics [j].bitmap = 0;
                                }
                                (*root)->characteristics [j].stats[0][adios_statistic_min].data = data;
                                (*root)->characteristics [j].bitmap |= (1 << adios_statistic_min);
                                break;

                            case adios_characteristic_max:
                                if (!(*root)->characteristics [j].stats)
                                {
                                    (*root)->characteristics [j].stats = malloc (sizeof(struct adios_index_characteristics_stat_struct *));
                                    (*root)->characteristics [j].stats[0] = malloc (2 * sizeof(struct adios_index_characteristics_stat_struct));
                                    (*root)->characteristics [j].bitmap = 0;
                                }
                                (*root)->characteristics [j].stats[0][adios_statistic_max].data = data;
                                (*root)->characteristics [j].bitmap |= (1 << adios_statistic_max);
                                break;
                        }
                        break;
                    }

                    // NCSU - Statistics - Parsing stat related info from bp file based on the bitmap
                    case adios_characteristic_stat:
                    {
                        uint8_t k, c, idx;
                        uint64_t count = adios_get_stat_set_count((*root)->type);
                        uint16_t characteristic_size;

                        (*root)->characteristics [j].stats = malloc (count * sizeof(struct adios_index_characteristics_stat_struct *));

                        for (c = 0; c < count; c ++)
                        {
                            (*root)->characteristics [j].stats[c] = calloc(ADIOS_STAT_LENGTH, sizeof(struct adios_index_characteristics_stat_struct));

                            k = idx = 0;
                            while ((*root)->characteristics[j].bitmap >> k)
                            {
                                (*root)->characteristics [j].stats[c][k].data = 0;

                                if (((*root)->characteristics[j].bitmap >> k) & 1)
                                {
                                    if (k == adios_statistic_hist)
                                    {
                                        struct adios_index_characteristics_hist_struct * hist = malloc(sizeof(struct adios_index_characteristics_hist_struct));
                                        uint32_t bi, num_breaks;

                                        (*root)->characteristics [j].stats[c][idx].data = hist;

                                        // Getting the number of breaks of histogram
                                        hist->num_breaks = * (uint32_t *) (b->buff + b->offset);
                                        if(b->change_endianness == adios_flag_yes) {
                                            swap_32(hist->num_breaks);
                                        }
                                        b->offset += 4;

                                        num_breaks = hist->num_breaks;

                                        // Getting the min of histogram
                                        hist->max = *(double *) (b->buff + b->offset);
                                        if(b->change_endianness == adios_flag_yes) {
                                            swap_64(hist->min);
                                        }
                                        b->offset += 8;

                                        // Getting the max of histogram
                                        hist->max = *(double *) (b->buff + b->offset);
                                        if(b->change_endianness == adios_flag_yes) {
                                            swap_64(hist->max);
                                        }
                                        b->offset += 8;

                                        // Getting the frequencies of the histogram
                                        hist->frequencies = malloc ((num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
                                        memcpy(hist->frequencies, (b->buff + b->offset), (num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));

                                        if(b->change_endianness == adios_flag_yes) {
                                            for(bi = 0; bi <= num_breaks; bi ++) {
                                                swap_32(hist->frequencies[bi]);
                                            }
                                        }
                                        b->offset += 4 * (num_breaks + 1);

                                        // Getting the breaks of the histogram
                                        hist->breaks = malloc (num_breaks * adios_get_type_size(adios_double, ""));
                                        memcpy(hist->breaks, (b->buff + b->offset), num_breaks * adios_get_type_size(adios_double, ""));
                                        if(b->change_endianness == adios_flag_yes) {
                                            for(bi = 0; bi < num_breaks; bi ++)
                                                swap_64(hist->breaks[bi]);
                                        }
                                        b->offset += 8 * num_breaks;
                                    }
                                    else
                                    {
                                        // NCSU - Generic for non-histogram data
                                        characteristic_size = adios_get_stat_size((*root)->characteristics [j].stats[c][idx].data, (*root)->type, k);
                                        (*root)->characteristics [j].stats[c][idx].data = malloc (characteristic_size);

                                        void * data = (*root)->characteristics [j].stats[c][idx].data;
                                        memcpy (data, (b->buff + b->offset), characteristic_size);
                                        b->offset += characteristic_size;

                                        if(b->change_endianness == adios_flag_yes)
                                            swap_ptr(data, characteristic_size * 8);
                                    }
                                    idx ++;
                                }
                                k ++;
                            }
                        }
                        break;
                    }

                    // NCSU - Reading bitmap value
                    case adios_characteristic_bitmap:
                    {
                        (*root)->characteristics [j].bitmap =
                                            *(uint32_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_32((*root)->characteristics [j].bitmap);
                        }
                        b->offset += 4;
                        break;
                    }

                    case adios_characteristic_offset:
                    {
                        uint64_t size = adios_get_type_size ((*root)->type, "");
                        (*root)->characteristics [j].offset =
                                            *(uint64_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_64((*root)->characteristics [j].offset);
                        }
                        b->offset += 8;

                        break;
                    }

                    case adios_characteristic_payload_offset:
                    {
                        uint64_t size = adios_get_type_size ((*root)->type, "");
                        (*root)->characteristics [j].payload_offset =
                                            *(uint64_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_64((*root)->characteristics [j].payload_offset);
                        }
                        b->offset += 8;

                        break;
                    }

                    case adios_characteristic_file_index:
                    {
                        (*root)->characteristics [j].file_index =
                                            *(uint32_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_32((*root)->characteristics [j].file_index);
                        }
                        b->offset += 4;

                        break;
                    }

                    case adios_characteristic_time_index:
                    {
                        (*root)->characteristics [j].time_index =
                                            *(uint32_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_32((*root)->characteristics [j].time_index);
                        }
                        b->offset += 4;

                        break;
                    }

                    case adios_characteristic_dimensions:
                    {
                        uint16_t dims_length;

                        (*root)->characteristics [j].dims.count =
                                           *(uint8_t *) (b->buff + b->offset);
                        b->offset += 1;

                        dims_length = *(uint16_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_16(dims_length);
                        }
                        b->offset += 2;

                       (*root)->characteristics [j].dims.dims = (uint64_t *)
                                                         malloc (dims_length);
                       memcpy ((*root)->characteristics [j].dims.dims
                              ,(b->buff + b->offset)
                              ,dims_length
                              );
                        if(b->change_endianness == adios_flag_yes) {
                            uint16_t di = 0;
                            uint16_t dims_num = dims_length / 8;
                            for (di = 0; di < dims_num; di ++) {
                                swap_64(((*root)->characteristics [j].dims.dims)[di]);
                            }
                        }
                        b->offset += dims_length;
                    }
                }
                item++;
            }
        }

        root = &(*root)->next;
    }

    return 0;
}

int adios_parse_attributes_index_v1 (struct adios_bp_buffer_struct_v1 * b
                          ,struct adios_index_attribute_struct_v1 ** attrs_root
                          )
{
    struct adios_index_attribute_struct_v1 ** root;

    if (b->length - b->offset < 10)
    {
        fprintf (stderr, "adios_parse_attributes_index_v1 requires a buffer "
                         "of at least 10 bytes.  Only %llu were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    root = attrs_root;

    uint16_t attrs_count;
    uint64_t attrs_length;

    attrs_count = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(attrs_count);
    }
    b->offset += 2;

    attrs_length = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(attrs_length);
    }
    b->offset += 8;

    // validate remaining length

    int i;
    for (i = 0; i < attrs_count; i++)
    {
        if (!*root)
        {
            *root = (struct adios_index_attribute_struct_v1 *)
                      malloc (sizeof (struct adios_index_attribute_struct_v1));
            (*root)->next = 0;
        }
        uint8_t flag;
        uint32_t attr_entry_length;
        uint16_t len;
        uint64_t characteristics_sets_count;
        uint64_t type_size;

        attr_entry_length = *(uint32_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_32(attr_entry_length);
        }
        b->offset += 4;

        (*root)->id = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16((*root)->id);
        }
        b->offset += 2;

        len = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(len);
        }
        b->offset += 2;
        (*root)->group_name = (char *) malloc (len + 1);
        (*root)->group_name [len] = '\0';
        strncpy ((*root)->group_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(len);
        }
        b->offset += 2;
        (*root)->attr_name = (char *) malloc (len + 1);
        (*root)->attr_name [len] = '\0';
        strncpy ((*root)->attr_name, b->buff + b->offset, len);
        b->offset += len;

        len = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(len);
        }
        b->offset += 2;
        (*root)->attr_path = (char *) malloc (len + 1);
        (*root)->attr_path [len] = '\0';
        strncpy ((*root)->attr_path, b->buff + b->offset, len);
        b->offset += len;

        flag = *(b->buff + b->offset);
        (*root)->type = (enum ADIOS_DATATYPES) flag;
        b->offset += 1;
        type_size = adios_get_type_size ((*root)->type, "");

        characteristics_sets_count = *(uint64_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_64(characteristics_sets_count);
        }
        (*root)->characteristics_count = characteristics_sets_count;
        (*root)->characteristics_allocated = characteristics_sets_count;
        b->offset += 8;

        // validate remaining length: offsets_count * (8 + 2 * (size of type))
        uint64_t j;
        (*root)->characteristics = malloc (characteristics_sets_count
                       * sizeof (struct adios_index_characteristic_struct_v1)
                      );
        memset ((*root)->characteristics, 0
               ,  characteristics_sets_count
                * sizeof (struct adios_index_characteristic_struct_v1)
               );
        for (j = 0; j < characteristics_sets_count; j++)
        {
            uint8_t characteristic_set_count;
            uint32_t characteristic_set_length;
            uint8_t item = 0;

            characteristic_set_count = (uint8_t) *(b->buff + b->offset);
            b->offset += 1;

            characteristic_set_length = *(uint32_t *) (b->buff + b->offset);
            if(b->change_endianness == adios_flag_yes) {
                swap_32(characteristic_set_length);
            }
            b->offset += 4;

            while (item < characteristic_set_count)
            {
                uint8_t flag;
                enum ADIOS_CHARACTERISTICS c;
                flag = *(b->buff + b->offset);
                c = (enum ADIOS_CHARACTERISTICS) flag;
                b->offset += 1;

                switch (c)
                {
                    case adios_characteristic_value:
                    {
                        uint16_t data_size;
                        void * data = 0;

                        if ((*root)->type == adios_string)
                        {
                            data_size = *(uint16_t *) (b->buff + b->offset);
                            if(b->change_endianness == adios_flag_yes) {
                                swap_16(data_size);
                            }
                            b->offset += 2;
                        }
                        else
                        {
                            data_size = adios_get_type_size ((*root)->type, "");
                        }

                        data = malloc (data_size + 1);
                        ((char *) data) [data_size] = '\0';

                        if (!data)
                        {
                            fprintf (stderr, "cannot allocate %d bytes to "
                                             "copy scalar %s\n"
                                    ,data_size
                                    ,(*root)->attr_name
                                    );

                            return 1;
                        }

                        switch ((*root)->type)
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
                                memcpy (data, (b->buff + b->offset), data_size);

                                if(b->change_endianness == adios_flag_yes) {
                                    if((*root)->type == adios_complex) {
                                        // TODO
                                    }
                                    else if((*root)->type == adios_double_complex) {
                                        // TODO
                                    }
                                    else {
                                        switch(data_size)
                                        {
                                            case 2:
                                                swap_16_ptr(data);
                                                break;
                                            case 4:
                                                swap_32_ptr(data);
                                                break;
                                            case 8:
                                                swap_64_ptr(data);
                                                break;
                                            case 16:
                                                swap_128_ptr(data);
                                                break;
                                        }
                                    }
                                }

                                b->offset += data_size;
                                break;
                            case adios_string:
                                memcpy (data, (b->buff + b->offset), data_size);
                                b->offset += data_size;
                                break;

                            default:
                                free (data);
                                data = 0;
                                break;
                        }

                        (*root)->characteristics [j].value = data;

                        break;
                    }

                    case adios_characteristic_offset:
                    {
                        uint64_t size = adios_get_type_size ((*root)->type, "");
                        (*root)->characteristics [j].offset =
                                            *(uint64_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_64((*root)->characteristics [j].offset);
                        }
                        b->offset += 8;

                        break;
                    }

                    case adios_characteristic_bitmap:
                    {
                        uint32_t size = adios_get_type_size ((*root)->type, "");
                        (*root)->characteristics [j].offset =
                                            *(uint32_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_32((*root)->characteristics [j].offset);
                        }
                        b->offset += 4;

                        break;
                    }

                    case adios_characteristic_payload_offset:
                    {
                        uint64_t size = adios_get_type_size ((*root)->type, "");
                        (*root)->characteristics [j].payload_offset =
                                            *(uint64_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_64((*root)->characteristics [j].payload_offset);
                        }
                        b->offset += 8;

                        break;
                    }

                    case adios_characteristic_file_index:
                    {
                        (*root)->characteristics [j].file_index =
                                            *(uint32_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_32((*root)->characteristics [j].file_index);
                        }
                        b->offset += 4;

                        break;
                    }

                    case adios_characteristic_time_index:
                    {
                        (*root)->characteristics [j].time_index =
                                            *(uint32_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_32((*root)->characteristics [j].time_index);
                        }
                        b->offset += 4;

                        break;
                    }

                    case adios_characteristic_var_id:
                    {
                        (*root)->characteristics [j].var_id =
                                            *(uint16_t *) (b->buff + b->offset);
                        if(b->change_endianness == adios_flag_yes) {
                            swap_16((*root)->characteristics [j].var_id);
                        }
                        b->offset += 2;

                        break;
                    }
                }
                item++;
            }
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
                         "Only %llu were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    uint64_t size;
    size = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(size);
    }
    b->offset += 8;

    pg_header->host_language_fortran =
                             (*(b->buff + b->offset) == 'y' ? adios_flag_yes
                                                            : adios_flag_no
                             );
    b->offset += 1;

    uint16_t len;
    uint8_t count;
    len = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(len);
    }
    b->offset += 2;
    pg_header->name = (char *) malloc (len + 1);
    pg_header->name [len] = '\0';
    memcpy (pg_header->name, b->buff + b->offset, len);
    b->offset += len;

    pg_header->coord_var_id = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(pg_header->coord_var_id);
    }
    b->offset += 2;
    len = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(len);
    }
    b->offset += 2;
    pg_header->time_index_name = (char *) malloc (len + 1);
    pg_header->time_index_name [len] = '\0';
    memcpy (pg_header->time_index_name, b->buff + b->offset, len);
    b->offset += len;

    pg_header->time_index = *(uint32_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_32(pg_header->time_index);
    }
    b->offset += 4;

    pg_header->methods_count = *(b->buff + b->offset);
    b->offset += 1;

    // length of methods section
    len = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(len);
    }
    b->offset += 2;

    int i;
    struct adios_method_info_struct_v1 ** root;

    pg_header->methods = 0;
    root = &pg_header->methods;
    for (i = 0; i < pg_header->methods_count; i++)
    {
        uint8_t flag;
        if (!*root)
        {
            *root = (struct adios_method_info_struct_v1 *)
                        malloc (sizeof (struct adios_method_info_struct_v1));
            (*root)->next = 0;
        }

        flag = *(b->buff + b->offset);
        (*root)->id = (enum ADIOS_IO_METHOD) flag;
        b->offset += 1;

        len = *(uint16_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_16(len);
        }
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

    if (b->length - b->offset < 10)
    {
        fprintf (stderr, "adios_parse_vars_header_v1 requires a "
                         "buffer of at least 10 bytes.  "
                         "Only %llu were provided\n"
                ,b->length - b->offset
                );

        vars_header->count = 0;
        vars_header->length = 0;

        return 1;
    }

    vars_header->count = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(vars_header->count);
    }
    b->offset += 2;
    vars_header->length = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(vars_header->length);
    }
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
                         "Only %llu were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    uint64_t initial_offset = b->offset;  // save to calc payload size
    uint64_t length_of_var;
    uint16_t len;
    uint8_t flag;

    length_of_var = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(length_of_var);
    }
    b->offset += 8;

    //validate remaining length

    var_header->id = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(var_header->id);
    }
    b->offset += 2;

    len = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(len);
    }
    b->offset += 2;
    var_header->name = (char *) malloc (len + 1);
    var_header->name [len] = '\0';
    memcpy (var_header->name, b->buff + b->offset, len);
    b->offset += len;

    len = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(len);
    }
    b->offset += 2;
    var_header->path = (char *) malloc (len + 1);
    var_header->path [len] = '\0';
    memcpy (var_header->path, b->buff + b->offset, len);
    b->offset += len;

    flag = *(b->buff + b->offset);
    var_header->type = (enum ADIOS_DATATYPES) flag;
    b->offset += 1;

    var_header->is_dim = (*(b->buff + b->offset) == 'y' ? adios_flag_yes
                                                        : adios_flag_no
                         );
    b->offset += 1;

    int i;
    uint8_t dims_count;
    uint16_t dims_length;
    uint8_t characteristics_count;
    uint32_t characteristics_length;

    // validate remaining length

    dims_count = *(b->buff + b->offset);
    b->offset += 1;
    dims_length = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(dims_length);
    }
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

        (*root)->dimension.rank = 0;
        (*root)->dimension.var_id = 0;
        (*root)->dimension.time_index = adios_flag_no;

        (*root)->global_dimension.rank = 0;
        (*root)->global_dimension.var_id = 0;
        (*root)->global_dimension.time_index = adios_flag_no;

        (*root)->local_offset.rank = 0;
        (*root)->local_offset.var_id = 0;
        (*root)->local_offset.time_index = adios_flag_no;

        flag = *(b->buff + b->offset);
        b->offset += 1;
        if (flag == 'y')
        {
            (*root)->dimension.rank = 0;
            (*root)->dimension.var_id = *(uint16_t *) (b->buff + b->offset);
            if(b->change_endianness == adios_flag_yes) {
                swap_16((*root)->dimension.var_id);
            }
            if ((*root)->dimension.var_id == 0)
                (*root)->dimension.time_index = adios_flag_yes;
            b->offset += 2;
        }
        else
        {
            (*root)->dimension.rank = *(uint64_t *) (b->buff + b->offset);
            if(b->change_endianness == adios_flag_yes) {
                swap_64((*root)->dimension.rank);
            }
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
            if(b->change_endianness == adios_flag_yes) {
                swap_16((*root)->global_dimension.var_id);
            }
            if ((*root)->global_dimension.var_id == 0)
                (*root)->global_dimension.time_index = adios_flag_yes;
            b->offset += 2;
        }
        else
        {
            (*root)->global_dimension.rank = *(uint64_t *)
                                                         (b->buff + b->offset);
            if(b->change_endianness == adios_flag_yes) {
                swap_64((*root)->global_dimension.rank);
            }
            (*root)->global_dimension.var_id = 0;
            b->offset += 8;
        }

        flag = *(b->buff + b->offset);
        b->offset += 1;
        if (flag == 'y')
        {
            (*root)->local_offset.rank = 0;
            (*root)->local_offset.var_id = *(uint16_t *) (b->buff + b->offset);
            if(b->change_endianness == adios_flag_yes) {
                swap_16((*root)->local_offset.var_id);
            }
            if ((*root)->local_offset.var_id == 0)
                (*root)->local_offset.time_index = adios_flag_yes;
            b->offset += 2;
        }
        else
        {
            (*root)->local_offset.rank = *(uint64_t *) (b->buff + b->offset);
            if(b->change_endianness == adios_flag_yes) {
                swap_64((*root)->local_offset.rank);
            }
            (*root)->local_offset.var_id = 0;
            b->offset += 8;
        }

        root = &(*root)->next;
    }

    characteristics_count = *(b->buff + b->offset);
    b->offset += 1;
    characteristics_length = *(uint32_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_32(characteristics_length);
    }
    b->offset += 4;

    uint64_t size = adios_get_type_size (var_header->type, "");

    var_header->characteristics.offset = 0;
    var_header->characteristics.payload_offset = 0;
    var_header->characteristics.stats = 0;
    var_header->characteristics.value = 0;
    var_header->characteristics.dims.count = 0;
    var_header->characteristics.dims.dims = 0;
    for (i = 0; i < characteristics_count; i++)
    {
        uint8_t flag;
        enum ADIOS_CHARACTERISTICS c;

        flag = *(b->buff + b->offset);
        b->offset += 1;
        c = (enum ADIOS_CHARACTERISTICS) flag;

        switch (c)
        {
            case adios_characteristic_offset:
                var_header->characteristics.offset = *(uint64_t *)
                                                        (b->buff + b->offset);
                if(b->change_endianness == adios_flag_yes) {
                    swap_64(var_header->characteristics.offset);
                }
                b->offset += 8;
                break;

            case adios_characteristic_payload_offset:
                var_header->characteristics.payload_offset = *(uint64_t *)
                                                        (b->buff + b->offset);
                if(b->change_endianness == adios_flag_yes) {
                    swap_64(var_header->characteristics.payload_offset);
                }
                b->offset += 8;
                break;

            //NCSU - Read in bitmap
            case adios_characteristic_bitmap:
                var_header->characteristics.bitmap = *(uint32_t *)
                                                        (b->buff + b->offset);
                if(b->change_endianness == adios_flag_yes) {
                    swap_32(var_header->characteristics.bitmap);
                }
                b->offset += 4;
                break;

            //NCSU - Read in the statistics
            case adios_characteristic_stat:
            {
                uint8_t j = 0, idx = 0;
                uint8_t c = 0, count = adios_get_stat_set_count (var_header->type);
                uint64_t characteristic_size;

                var_header->characteristics.stats = malloc(count * sizeof(struct adios_index_characteristics_stat_struct *));

                for (c = 0; c < count; c ++)
                {
                    idx = j = 0;
                    var_header->characteristics.stats[c] = calloc(ADIOS_STAT_LENGTH, sizeof(struct adios_index_characteristics_stat_struct));

                    while (var_header->characteristics.bitmap >> j)
                    {
                        if ((var_header->characteristics.bitmap >> j) & 1)
                        {
                            if (j == adios_statistic_hist)
                            {
                                struct adios_index_characteristics_hist_struct * hist = malloc (sizeof(struct adios_index_characteristics_hist_struct));
                                uint32_t bi, num_breaks;

                                var_header->characteristics.stats[c][idx].data = hist;

                                // Getting the number of bins
                                hist->num_breaks = * (uint32_t *) (b->buff + b->offset);

                                if(b->change_endianness == adios_flag_yes) {
                                    swap_32(hist->num_breaks);
                                }
                                b->offset += 4;

                                num_breaks = hist->num_breaks;

                                // Getting the min hist value
                                hist->min = * (double *) (b->buff + b->offset);
                                if(b->change_endianness == adios_flag_yes) {
                                    swap_64(hist->min);
                                }
                                b->offset += 8;

                                // Getting the max hist value
                                hist->max = * (double *) (b->buff + b->offset);

                                if(b->change_endianness == adios_flag_yes) {
                                    swap_64(hist->max);
                                }
                                b->offset += 8;

                                // Getting the frequencies
                                hist->frequencies = malloc((num_breaks + 1)  * adios_get_type_size(adios_unsigned_integer, ""));
                                memcpy (hist->frequencies
                                       ,(b->buff + b->offset), (num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, "")
                                       );

                                if(b->change_endianness == adios_flag_yes) {
                                    for(bi = 0; bi <= num_breaks; bi++)
                                        swap_32(hist->frequencies[bi]);
                                }
                                b->offset += 4 * (num_breaks + 1);

                                // Getting the break points
                                hist->breaks = malloc(num_breaks * adios_get_type_size(adios_double, ""));
                                memcpy (hist->breaks, (b->buff + b->offset), num_breaks * 8);

                                if(b->change_endianness == adios_flag_yes) {
                                    for(bi = 0; bi < num_breaks; bi++)
                                        swap_64(hist->breaks[bi]);
                                }
                                b->offset += 8 * num_breaks;
                            }
                            else
                            {
                                characteristic_size = adios_get_stat_size(var_header->characteristics.stats[c][idx].data, var_header->type, j);

                                var_header->characteristics.stats[c][idx].data = malloc(characteristic_size);
                                memcpy (var_header->characteristics.stats[c][idx].data, (b->buff + b->offset)
                                   ,characteristic_size
                                   );

                                if(b->change_endianness == adios_flag_yes) {
                                    swap_adios_type(var_header->characteristics.stats[c][idx].data, var_header->type);
                                }

                                b->offset += characteristic_size;
                            }

                            idx ++;
                        }
                        j ++;
                    }
                }
                break;
            }

            case adios_characteristic_dimensions:
            {
                uint8_t dim_count;
                uint16_t dim_length;

                dim_count = *(b->buff + b->offset);
                b->offset += 1;
                dim_length = *(uint16_t *) (b->buff + b->offset);
                if(b->change_endianness == adios_flag_yes) {
                    swap_16(dim_length);
                }
                b->offset += 2;

                var_header->characteristics.dims.dims = malloc (dim_length);
                memcpy (var_header->characteristics.dims.dims
                       ,(b->buff + b->offset), dim_length
                       );
                if(b->change_endianness == adios_flag_yes) {
                    uint16_t di = 0;
                    uint16_t dim_num = dim_length / 8;
                    for (di = 0; di < dim_num; di ++) {
                         swap_64((var_header->characteristics.dims.dims)[di]);
                    }
                }
                b->offset += dim_length;
                break;
            }

            // NCSU - Adding backward compatibility
            // Reading min and max here implies older bp file
            case adios_characteristic_min:
                if (!var_header->characteristics.stats)
                {
                    var_header->characteristics.stats = malloc(sizeof(struct adios_index_characteristics_stat_struct *));
                    var_header->characteristics.stats[0] = malloc(2 * sizeof(struct adios_index_characteristics_stat_struct));
                    var_header->characteristics.bitmap = 0;
                }
                var_header->characteristics.stats[0][adios_statistic_min].data = malloc (size);
                var_header->characteristics.bitmap |= (1 << adios_statistic_min);
                memcpy (var_header->characteristics.stats[0][adios_statistic_min].data, (b->buff + b->offset)
                       ,size
                       );

                if(b->change_endianness == adios_flag_yes) {
                    swap_adios_type(var_header->characteristics.stats[0][adios_statistic_min].data, var_header->type);
                }
                b->offset += size;
                break;

            case adios_characteristic_max:
                if (!var_header->characteristics.stats)
                {
                    var_header->characteristics.stats = malloc(sizeof(struct adios_index_characteristics_stat_struct *));
                    var_header->characteristics.stats[0] = malloc(2 * sizeof(struct adios_index_characteristics_stat_struct));
                    var_header->characteristics.bitmap = 0;
                }
                var_header->characteristics.stats[0][adios_statistic_max].data = malloc (size);
                var_header->characteristics.bitmap |= (1 << adios_statistic_max);
                memcpy (var_header->characteristics.stats[0][adios_statistic_max].data, (b->buff + b->offset)
                       ,size
                       );

                if(b->change_endianness == adios_flag_yes) {
                    swap_adios_type(var_header->characteristics.stats[0][adios_statistic_max].data, var_header->type);
                }
                b->offset += size;
                break;

            case adios_characteristic_value:
            {
                uint16_t val_size = *(uint16_t *) (b->buff + b->offset);
                if(b->change_endianness == adios_flag_yes) {
                    swap_16(val_size);
                }
                b->offset += 2;

                var_header->characteristics.value = malloc (val_size + 1);
                ((char *) var_header->characteristics.value) [val_size] = '\0';
                memcpy (var_header->characteristics.value, (b->buff + b->offset)
                       ,val_size
                       );

                // TODO: do we need byte-swap here?
                b->offset += val_size;
                break;
            }
        }
    }

    var_header->payload_size = length_of_var - (b->offset - initial_offset);

    return 0;
}

int adios_clear_process_group_header_v1 (
                       struct adios_process_group_header_struct_v1 * pg_header)
{
    pg_header->host_language_fortran = adios_flag_unknown;
    if (pg_header->name)
    {
        free (pg_header->name);
        pg_header->name = 0;
    }
    pg_header->coord_var_id = 0;
    if (pg_header->time_index_name)
    {
        free (pg_header->time_index_name);
        pg_header->time_index_name = 0;
    }
    pg_header->time_index = 0;
    while (pg_header->methods)
    {
        struct adios_method_info_struct_v1 * t = pg_header->methods->next;
        pg_header->methods->id = (enum ADIOS_IO_METHOD) 0;
        if (pg_header->methods->parameters)
        {
            free (pg_header->methods->parameters);
            pg_header->methods->parameters = 0;
        }
        free (pg_header->methods);
        pg_header->methods = t;
    }
    pg_header->methods_count = 0;

    return 0;
}

int adios_clear_var_header_v1 (struct adios_var_header_struct_v1 * var_header)
{
    if (var_header->name)
    {
        free (var_header->name);
        var_header->name = 0;
    }
    if (var_header->path)
    {
        free (var_header->path);
        var_header->path = 0;
    }
    while (var_header->dims)
    {
        struct adios_dimension_struct_v1 * d = var_header->dims->next;
        free (var_header->dims);
        var_header->dims = d;
    }
    struct adios_index_characteristic_struct_v1 * c
                                                = &var_header->characteristics;

    c->offset = 0;

    // NCSU - Clear stat
    if (c->stats)
    {
        uint8_t j = 0, idx = 0;
        uint8_t i = 0, count = adios_get_stat_set_count(var_header->type);

        while (c->bitmap >> j)
        {
            if ((c->bitmap >> j) & 1)
            {
                for (i = 0; i < count; i ++)
				{
					if (j == adios_statistic_hist)
					{
						struct adios_index_characteristics_hist_struct * hist = (struct adios_index_characteristics_hist_struct *) (c->stats[i][idx].data);
						free (hist->breaks);
						free (hist->frequencies);
						free (hist);
					}
					else
                    	free (c->stats[i][idx].data);
				}
                idx ++;
            }
            j ++;
        }

        for (i = 0; i < count; i ++)
            free (c->stats [i]);

        free (c->stats);
        c->stats = 0;
    }
    if (c->dims.dims)
    {
        free (c->dims.dims);
        c->dims.count = 0;
        c->dims.dims = 0;
    }
    if (c->value)
    {
        free (c->value);
        c->value = 0;
    }
    c->var_id = 0;

    return 0;
}

int adios_parse_var_data_payload_v1 (struct adios_bp_buffer_struct_v1 * b
                             ,struct adios_var_header_struct_v1 * var_header
                             ,struct adios_var_payload_struct_v1 * var_payload
                             ,uint64_t payload_buffer_size
                             )
{
    if (b->length - b->offset < var_header->payload_size)
    {
        fprintf (stderr, "adios_parse_var_data_payload_v1 for name %s "
                         "path %s requires a "
                         "buffer of at least %llu bytes.  "
                         "Only %llu were provided\n"
                ,var_header->name, var_header->path
                ,var_header->payload_size, b->length - b->offset
                );

        b->offset += var_header->payload_size;

        return 1;
    }
    if (   payload_buffer_size < var_header->payload_size
        && var_payload && var_payload->payload
       )
    {
        fprintf (stderr, "reading var name %s path %s requires a "
                         "buffer of at least %llu bytes.  "
                         "Only %llu were provided\n"
                ,var_header->name, var_header->path
                ,var_header->payload_size, payload_buffer_size
                );

        b->offset += var_header->payload_size;

        return 1;
    }

    uint64_t size = adios_get_type_size (var_header->type, "");

    if (var_payload)
    {
        if (var_payload->payload)
        {
            memcpy (var_payload->payload, (b->buff + b->offset)
                   ,var_header->payload_size
                   );
            if(b->change_endianness == adios_flag_yes) {
                swap_adios_type_array(var_payload->payload, var_header->type, var_header->payload_size);
            }
            if (var_header->type == adios_string) {
                ((char*)var_payload->payload)[var_header->payload_size] = '\0';
            }
            b->offset += var_header->payload_size;
        }
        else
        {
            b->offset += var_header->payload_size;
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
        fprintf (stderr, "adios_parse_attributes_header_v1 requires a "
                         "buffer of at least 10 bytes.  "
                         "Only %llu were provided\n"
                ,b->length - b->offset
                );

        attrs_header->count = 0;
        attrs_header->length = 0;

        return 1;
    }

    attrs_header->count = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(attrs_header->count);
    }
    b->offset += 2;

    attrs_header->length = *(uint64_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_64(attrs_header->length);
    }
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
                         "Only %llu were provided\n"
                ,b->length - b->offset
                );

        return 1;
    }

    uint32_t attribute_length;
    uint16_t len;

    attribute_length = *(uint32_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_32(attribute_length);
    }
    b->offset += 4;
    attribute->id = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(attribute->id);
    }
    b->offset += 2;

    len = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(len);
    }
    b->offset += 2;
    attribute->name = (char *) malloc (len + 1);
    attribute->name [len] = '\0';
    strncpy (attribute->name, b->buff + b->offset, len);
    b->offset += len;

    len = *(uint16_t *) (b->buff + b->offset);
    if(b->change_endianness == adios_flag_yes) {
        swap_16(len);
    }
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
        if(b->change_endianness == adios_flag_yes) {
            swap_16(attribute->var_id);
        }
        b->offset += 2;
        attribute->type = adios_unknown;
        attribute->length = 0;
        attribute->value = 0;
    }
    else
    {
        uint8_t flag;
        attribute->var_id = 0;
        flag = *(b->buff + b->offset);
        attribute->type = (enum ADIOS_DATATYPES) flag;
        b->offset += 1;
        attribute->length = *(uint32_t *) (b->buff + b->offset);
        if(b->change_endianness == adios_flag_yes) {
            swap_32(attribute->length);
        }
        b->offset += 4;

        attribute->value = malloc (attribute->length + 1);
        ((char *) attribute->value) [attribute->length] = '\0';
        memcpy (attribute->value, (b->buff + b->offset), attribute->length);

        // TODO: do we need byte-swap here?
        if(b->change_endianness == adios_flag_yes) {
            swap_adios_type(attribute->value, attribute->type);
        }
        b->offset += attribute->length;
    }

    return 0;
}

int adios_clear_attribute_v1 (struct adios_attribute_struct_v1 * attribute)
{
    attribute->id = 0;
    if (attribute->name)
    {
        free (attribute->name);
        attribute->name = 0;
    }
    if (attribute->path)
    {
        free (attribute->path);
        attribute->path = 0;
    }
    attribute->is_var = adios_flag_unknown;
    attribute->var_id = 0;
    attribute->type = adios_unknown;
    attribute->length = 0;
    if (attribute->value)
    {
        free (attribute->value);
        attribute->value = 0;
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

// *****************************************************************************
void adios_init_buffer_read_version (struct adios_bp_buffer_struct_v1 * b)
{
    if (!b->buff)
    {
        alloc_aligned (b, 28);
        memset (b->buff, 0, 28);
        if (!b->buff)
            fprintf(stderr, "could not allocate 28 bytes\n");
        b->offset = 24;
    }
}

// last 4 bytes of file
void adios_posix_read_version (struct adios_bp_buffer_struct_v1 * b)
{
    uint64_t buffer_size;
    uint64_t start;
    uint64_t r;

    adios_init_buffer_read_version (b);

    lseek (b->f, b->file_size - 28, SEEK_SET);

    r = read (b->f, b->buff, 28);
    if (r != 28)
        fprintf (stderr, "could not read 28 bytes. read only: %llu\n", r);
}

void adios_init_buffer_read_index_offsets (struct adios_bp_buffer_struct_v1 * b)
{
    b->offset = 0; // just move to the start of the buffer
}

void adios_posix_read_index_offsets (struct adios_bp_buffer_struct_v1 * b)
{
    adios_init_buffer_read_index_offsets (b);
}

void adios_init_buffer_read_process_group_index (
                                          struct adios_bp_buffer_struct_v1 * b
                                          )
{
    realloc_aligned (b, b->pg_size);
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
    realloc_aligned (b, b->vars_size);
    b->offset = 0;
}

void adios_posix_read_vars_index (struct adios_bp_buffer_struct_v1 * b)
{
    adios_init_buffer_read_vars_index (b);

    uint64_t r;

    lseek (b->f, b->vars_index_offset, SEEK_SET);
    r = read (b->f, b->buff, b->vars_size);

    if (r != b->vars_size)
        fprintf (stderr, "reading vars_index: wanted %llu, read: %llu\n"
                ,b->vars_size, r
                );
}

void adios_init_buffer_read_attributes_index
                                        (struct adios_bp_buffer_struct_v1 * b)
{
    realloc_aligned (b, b->attrs_size);
    b->offset = 0;
}

void adios_posix_read_attributes_index (struct adios_bp_buffer_struct_v1 * b)
{
    adios_init_buffer_read_attributes_index (b);

    uint64_t r;

    lseek (b->f, b->attrs_index_offset, SEEK_SET);
    r = read (b->f, b->buff, b->attrs_size);

    if (r != b->attrs_size)
        fprintf (stderr, "reading attributess_index: wanted %llu, read: %llu\n"
                ,b->attrs_size, r
                );
}

void adios_init_buffer_read_process_group (struct adios_bp_buffer_struct_v1 * b)
{
    realloc_aligned (b, b->read_pg_size);
    b->offset = 0;
}

uint64_t adios_posix_read_process_group (struct adios_bp_buffer_struct_v1 * b)
{
    uint64_t pg_size = 0;

    adios_init_buffer_read_process_group (b);

    do
    {
        lseek (b->f, b->read_pg_offset + pg_size, SEEK_SET);
        pg_size += read (b->f, b->buff + pg_size, b->read_pg_size - pg_size);
    } while (errno && pg_size != b->read_pg_size);

    if (pg_size != b->read_pg_size)
    {
        fprintf (stderr, "adios_read_process_group: "
                         "Tried to read: %llu, but only got: %llu error: %s\n"
                ,b->read_pg_size, pg_size, strerror (errno)
                );

        pg_size = 0;
    }

    return pg_size;
}

int adios_posix_open_read_internal (const char * filename
                                   ,const char * base_path
                                   ,struct adios_bp_buffer_struct_v1 * b
                                   )
{
    char * name;
    struct stat s;

    name = malloc (strlen (base_path) + strlen (filename) + 1);
    sprintf (name, "%s%s", base_path, filename);

    if (stat (name, &s) == 0)
        b->file_size = s.st_size;

    b->f = open (name, O_RDONLY | O_LARGEFILE);
    if (b->f == -1)
    {
        fprintf (stderr, "ADIOS POSIX: file not found: %s\n", name);

        free (name);

        return 0;
    }

    free (name);

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

uint64_t adios_get_type_size (enum ADIOS_DATATYPES type, void * var)
{
    switch (type)
    {
        case adios_byte:
        case adios_unsigned_byte:
            return 1;

        case adios_string:
            if (!var)
                return 0;
            else
                return strlen ((char *) var);

        case adios_short:
        case adios_unsigned_short:
            return 2;

        case adios_integer:
        case adios_unsigned_integer:
            return 4;

        case adios_long:
        case adios_unsigned_long:
            return 8;

        case adios_real:
            return 4;

        case adios_double:
            return 8;

        case adios_long_double:
            return 16;

        case adios_complex:
            return 2 * 4;

        case adios_double_complex:
            return 2 * 8;

        default:
            return -1;
    }
}

//NCSU - Get the size of the statistic id
uint64_t adios_get_stat_size (void * data, enum ADIOS_DATATYPES type, enum ADIOS_STAT stat_id)
{
    uint64_t size = 0;

    switch (type)
    {
        case adios_complex:
            switch (stat_id)
            {
                case adios_statistic_min:
                case adios_statistic_max:
                case adios_statistic_sum:
                case adios_statistic_sum_square:
                    return adios_get_type_size(adios_double, "");

                case adios_statistic_finite:
                    return adios_get_type_size(adios_byte, "");

                case adios_statistic_cnt:
                    return adios_get_type_size(adios_unsigned_integer, "");
            }
        case adios_double_complex:
            switch (stat_id)
            {
                case adios_statistic_min:
                case adios_statistic_max:
                case adios_statistic_sum:
                case adios_statistic_sum_square:
                    return adios_get_type_size(adios_long_double, "");

                case adios_statistic_finite:
                    return adios_get_type_size(adios_byte, "");

                case adios_statistic_cnt:
                    return adios_get_type_size(adios_unsigned_integer, "");
            }
        default:
        {
            switch (stat_id)
            {
                case adios_statistic_min:
                case adios_statistic_max:
                    return adios_get_type_size(type, "");

                case adios_statistic_sum:
                case adios_statistic_sum_square:
                    return adios_get_type_size(adios_double, "");

                case adios_statistic_finite:
                    return adios_get_type_size(adios_byte, "");

                case adios_statistic_cnt:
                    return adios_get_type_size(adios_unsigned_integer, "");

                case adios_statistic_hist:
                {
                    struct adios_hist_struct * hist = (struct adios_hist_struct *) data;
                    size += adios_get_type_size(adios_unsigned_integer, ""); // Number of breaks
                    size += adios_get_type_size(adios_double, ""); // Min
                    size += adios_get_type_size(adios_double, ""); // Max
                    size += ((hist->num_breaks + 1) * adios_get_type_size(adios_unsigned_integer, ""));
                    size += (hist->num_breaks * adios_get_type_size(adios_double, ""));
                    return size;
                }
                default:
                    return -1;
            }
        }
    }
}

//NCSU - Returns count of the set of characteristcs for a variable
uint8_t adios_get_stat_set_count (enum ADIOS_DATATYPES type)
{
    // NCSU - Complex numbers have the statistic related values repeated.
    // NCSU - Holds values for real, imaginary, complex parts
    if (type == adios_complex || type == adios_double_complex)
        return 3;
    return 1;
}

