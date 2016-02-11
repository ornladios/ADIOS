/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "adios_types.h"
#include "adios_internals.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "bp_utils.h"
#include "adios_transforms_common.h" // NCSU ALACRITY-ADIOS
#include "adios_transforms_read.h" // NCSU ALACRITY-ADIOS
//#include "adios_internals.h"

#define DIVIDER "========================================================\n"

int do_write_index = 0; // write recovered index at the end, default is no

void print_process_group_header (uint64_t num
                      ,struct adios_process_group_header_struct_v1 * pg_header
                      );
void print_vars_header (struct adios_vars_header_struct_v1 * vars_header);
void print_var_header (struct adios_var_header_struct_v1 * var_header);
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

/* Temporarily store scalar's values to be used to determine
   an array's dimensions. Indexed by variable id (which is the reference 
   in the array dimension to the scalar)
*/
#define MAX_DIMENSION_INDEX 1024
uint64_t dim_value[MAX_DIMENSION_INDEX];
const uint64_t INVALID_DIM = (uint64_t) -1L;

void init_dimensions ()
{
    int i;
    for (i=0; i < MAX_DIMENSION_INDEX; i++) {
        dim_value[i] = INVALID_DIM;
    }
}

/* Store a scalar variable's values temporarily while we process the 
   dimensions of the arrays in the same PG */
void store_scalar_dimensions (
        struct adios_var_header_struct_v1 * var_header,
        struct adios_var_payload_struct_v1 * var_payload)
{
    if (var_header->is_dim == adios_flag_yes)
    {
        if (var_header->id < MAX_DIMENSION_INDEX) {
            uint64_t d = 0L;
            switch (var_header->type)
            {
                case adios_byte:
                    d = *(signed char *) var_payload->payload;
                    break;
                case adios_unsigned_byte:
                    d = *(unsigned char *) var_payload->payload;
                    break;
                case adios_short:
                    d = *(signed short *) var_payload->payload;
                    break;
                case adios_unsigned_short:
                    d = *(unsigned short *) var_payload->payload;
                    break;
                case adios_integer:
                    d = *(signed int *) var_payload->payload;
                    break;
                case adios_unsigned_integer:
                    d = *(unsigned int *) var_payload->payload;
                    break;
                case adios_long:
                    d = *(signed long long *) var_payload->payload;
                    break;
                case adios_unsigned_long:
                    d = *(unsigned long long *) var_payload->payload;
                    break;

                default:
                    printf ("  !!! ERROR !!! The variable %s/%s appears to be used "
                            "as dimension but it has non-integer type.\n",
                           var_header->path, var_header->name);
                    break;
            }
            dim_value [var_header->id] = d;


        } else {
            printf ("  !!! ERROR !!! This tool can only handle scalar variables"
                    " used as dimensions if they appear as the first %d variables. "
                    "Here we encountered a scalar named %s/%s with id=%d\n",
                    MAX_DIMENSION_INDEX, var_header->path, var_header->name, var_header->id);
        }
    }
}

/* Return the actual dimension/offset value.
   If it points to a scalar variable, get its value.
   If it's time return the value provided by the caller 
     (should be 1 for a dimension, 0 for an offset)
 */
uint64_t get_dimension (struct adios_dimension_item_struct_v1 * d, int return_for_time)
{
    int id = d->var_id; 
    uint64_t dim = 0L;
    if (id == 0)
    {
        if (d->is_time_index == adios_flag_yes)
        {
            dim = return_for_time;
        }
        else
        {
            dim = d->rank;
        }
    } 
    else 
    {
        if (id < 0 || id > MAX_DIMENSION_INDEX) 
        {
            printf ("  !!! ERROR !!! Reference to dimension of variable id %d cannot be handled."
                    " We only store scalar variables up to id %d\n",
                    id, MAX_DIMENSION_INDEX);
            dim = INVALID_DIM;
        }
        else if (dim_value[id] == INVALID_DIM) 
        {
            printf ("  !!! ERROR !!! Reference to dimension of variable id %d but "
                    " we haven't stored this scalar variable. Result is invalid dimension\n",
                    id);
            dim = INVALID_DIM;
        }
        else 
        {
            dim = dim_value[id];
        }
    }
    return dim;
}


int MAX_GROUP_NAME_LENGTH = 64; // could be 65535 but who does that??? 

/* print char with printable char or its \xxx code
   return 1 if its ASCII 32..126, usable for names in ADIOS
*/
int print_namechar (char c)
{
    if (32 <= c && c <= 126) {
        // SPACE, symbols, alphanumeric
        printf ("%c",c);
        return 1;
    } else {
        printf (" \\%3.3hhu", c);
        return 0;
    }
}

int looks_like_PG (const char * buf, int blen)
{
    printf ("  Check if it looks like a PG...\n");
    int offset = 8; // skip now the PG size

    // host_language_fortran flag, single char
    char fort = buf[offset]; 
    offset += 1;
    printf ("       Fortran flag char should be 'y' or 'n': '");
    print_namechar (fort);
    printf ("'\n");
    if (fort != 'y' && fort != 'n')
        return 0; 

    // group name length, 16 bit
    uint16_t namelen =  *(uint16_t *) (buf + offset);
    offset += 2;
    printf ("       Group name length expected to be less than %d characters: %hu\n", 
            MAX_GROUP_NAME_LENGTH, namelen);
    if (namelen > MAX_GROUP_NAME_LENGTH)
        return 0;

    char gname[MAX_GROUP_NAME_LENGTH];
    memcpy (gname, buf+offset, namelen);
    offset += namelen;

    int i;
    int validname = 1;
    printf ("       Group name: \"");
    for (i = 0; i < namelen; i++)
    {
        validname &= print_namechar(gname[i]);
    }
    printf ("\"\n");

    if (!validname) {
        printf ("       Group name contains characters that are invalid for a name\n");
        return 0;
    }

    return 1;
}

/* Look for a PG at a given offset. If it looks like a PG, 
   return 1 and also return the reported PG size.
   pgsize_reported is changed only when returning 1.
*/
int find_pg (int fd, uint64_t offset, uint64_t file_size, uint64_t * pgsize_reported) 
{

    const int N_READ_AHEAD = 1024;
    char buf[N_READ_AHEAD]; // temporary buffer to read data in and parse for info

    printf ("Look for a Process Group (PG) at offset %" PRIu64 "\n", offset);
    lseek (fd, offset, SEEK_SET);
    // read a few bytes in
    errno = 0;
    int blen = read (fd, buf, N_READ_AHEAD);

    if (blen < 8) {
        printf ("  === Could not read even 8 bytes. Finish PG reading\n");
        if (errno) {
            printf ("  Error when reading: %s\n", strerror(errno));
        }
        return 0;
    }

    uint64_t pgsize;
    pgsize = *(uint64_t*) buf;  // first 8 bytes is size of PG
    printf ("  PG reported size: %" PRIu64 "\n", pgsize);

    if (pgsize < 28) {
        printf ("   === PG reported size is too small. This is not a (good) PG.\n");
        return 0;
    }

    if (pgsize + offset > file_size + 1024*1024*1024 /* a GB index??? */ ) {
        printf ("   === Offset + PG reported size >> file size. This is not a (good) PG.\n");
        return 0;
    }

    if (!looks_like_PG (buf, blen)) {
        return 0;
    }

    // All tests passed, it seems to be a PG
    *pgsize_reported = pgsize;
    return 1;
}


void add_pg_to_index (
        struct adios_index_struct_v1 *  index, 
        struct adios_process_group_header_struct_v1 * pg_header, 
        uint64_t curr_offset)
{
    /* similar to code from adios_internals.c:adios_build_index_v1() */
    struct adios_index_process_group_struct_v1 * g_item;
    g_item = (struct adios_index_process_group_struct_v1 *)
        malloc (sizeof (struct adios_index_process_group_struct_v1));
    g_item->group_name = strdup (pg_header->name);
    g_item->adios_host_language_fortran = pg_header->host_language_fortran;

    /* FIXME: process id (rank) is only recorded in the index, but not in the PG header.
       This should be the original MPI rank of the writer. */
    g_item->process_id = 0; 

    g_item->time_index_name = strdup (pg_header->time_index_name);
    g_item->time_index = pg_header->time_index;
    g_item->offset_in_file = curr_offset;
    g_item->next = 0;

    // add to the groups
    index_append_process_group_v1 (index, g_item);
}


/* Get the dimensions of the variable and add it to the index.
   This is complicated because of references to other variables */
void process_dimensions (
        struct adios_var_header_struct_v1 * var_header,
        struct adios_index_var_struct_v1 * v_index
        )
{
    // Get number of dimensions first
    struct adios_dimension_struct_v1 * d = var_header->dims;
    v_index->characteristics [0].dims.count = 0;
    while (d)
    {
        v_index->characteristics [0].dims.count++;
        d = d->next;
    }

    if (v_index->characteristics [0].dims.count)
    {
        // (local, global, local offset)
        v_index->characteristics [0].dims.dims = malloc
            (3 * 8 * v_index->characteristics [0].dims.count);

        int j;
        d = var_header->dims;
        for (j = 0; j < v_index->characteristics [0].dims.count; j++)
        {
            v_index->characteristics [0].dims.dims [j * 3 + 0] =
                get_dimension (&d->dimension, 1);
            v_index->characteristics [0].dims.dims [j * 3 + 1] =
                get_dimension (&d->global_dimension, 0);
            v_index->characteristics [0].dims.dims [j * 3 + 2] =
                get_dimension (&d->local_offset, 0);

            d = d->next;
        }
    } 
    else
    {
        v_index->characteristics [0].dims.dims = NULL;
    }
}

void add_var_to_index (
        struct adios_index_struct_v1 *  index, 
        struct adios_process_group_header_struct_v1 * pg_header, 
        struct adios_var_header_struct_v1 * var_header,
        struct adios_var_payload_struct_v1 * var_payload)
{

    /* similar to code from adios_internals.c:adios_build_index_v1() */
    struct adios_index_var_struct_v1 * v_index;
    v_index = malloc (sizeof (struct adios_index_var_struct_v1));
    v_index->characteristics = malloc (
            sizeof (struct adios_index_characteristic_struct_v1)
            );

    v_index->id = var_header->id;
    v_index->group_name = strdup (pg_header->name);
    v_index->var_name = (var_header->name ? strdup (var_header->name) : 0L);
    v_index->var_path = (var_header->path ? strdup (var_header->path) : 0L);
    v_index->type = var_header->type;

    v_index->characteristics_count = 1;
    v_index->characteristics_allocated = 1;
    v_index->characteristics [0].offset = var_header->characteristics.offset;
    v_index->characteristics [0].payload_offset = var_header->characteristics.payload_offset;

    v_index->characteristics [0].file_index = -1; // It works only with single BP files
    v_index->characteristics [0].time_index = pg_header->time_index;

    v_index->characteristics [0].value = 0;

    /* Determine the dimensions from actual values or references of scalars */
    process_dimensions (var_header, v_index);


    // NCSU - Initializing stat related info in index
    v_index->characteristics [0].bitmap = 0;
    v_index->characteristics [0].stats = 0;

    memcpy (&v_index->characteristics[0].transform,
            &var_header->characteristics.transform,
            sizeof (struct adios_index_characteristic_transform_struct));

    // Save scalar's value
    if (var_payload->payload)
    {
        uint64_t size = adios_get_type_size (var_header->type, var_payload->payload);
        if (var_header->type == adios_string) {
            v_index->characteristics [0].value = malloc (size + 1);
            memcpy (v_index->characteristics [0].value, var_payload->payload, size);
            ((char *) (v_index->characteristics [0].value)) [size] = 0;
        } else {
            v_index->characteristics [0].value = malloc (size);
            memcpy (v_index->characteristics [0].value, var_payload->payload, size);
        }
    }

    /* FIXME: Missing statistics and transformation statistics */

    v_index->next = 0;

    // this fn will either take ownership for free
    //printf ("  add index var %s/%s\n", v_index->var_path, v_index->var_name);
    index_append_var_v1 (index, v_index, 1);
}

void write_index (int fd, uint64_t file_offset, struct adios_index_struct_v1 *  index)
{
    char * buffer = NULL;
    uint64_t buffer_size = 0L;
    uint64_t buffer_offset = 0L;
    adios_write_index_v1 (&buffer, &buffer_size, &buffer_offset, file_offset, index);
    adios_write_version_v1 (&buffer, &buffer_size, &buffer_offset);
    printf ("Index metadata size is %" PRIu64 " bytes\n", buffer_offset);
    if  (do_write_index) 
    {
        printf ("Ready to write index to file offset %" PRIu64 ". Size is %" PRIu64 " bytes\n",
                file_offset, buffer_offset);
        lseek (fd, file_offset, SEEK_SET);
        ssize_t s = write (fd, buffer, buffer_offset);
        if (s != buffer_offset) {
            printf ("ERROR: Only wrote %llu bytes of index data instead of the attempted %" PRIu64 " bytes\n",
                    (unsigned long long)s, buffer_offset);
        } else {
            printf ("Wrote index to file offset %" PRIu64 ". Size is %" PRIu64 " bytes\n", file_offset, buffer_offset);
        }
        printf ("Truncate file to size %" PRIu64 "\n", file_offset+s);
        ftruncate (fd, file_offset+s);
    }
    else 
    {
        printf ("To actually recover the file by rewriting the index, use -f option.\n"
                "Final file will be truncted to %" PRIu64 " bytes\n", file_offset+buffer_offset);
    }
}

void print_usage (int argc, char ** argv)
{
    printf ("Usage: %s [-f | --force] <filename>\n"
            "  -f:  do write the recovered index to the end of file\n",
            argv [0]); 
    printf (
"This recovery tool parses the data blocks in the file, reconstructs the "
"index metadata and writes it to the end. It is useful only if the original "
"index is damaged somehow. It only works with a single BP file. Subfiles are "
"not handled. Also, the tool is very limited. It does not work if the file "
"has variables with transformations (compression). It does not recover attributes "
"nor statistics. It does not recover the file beyond the first data corruption in "
"the middle. So use it with care, copy the corrupted file before "
"using this tool. Without -f option, you can test if the processing goes well "
"without changing the file.\n"
    );
}

int main (int argc, char ** argv)
{
    char * filename;

    if (argc < 2 || argc > 3)
    {
        print_usage (argc, argv);
        return -1;
    }

    if (argv [1][0] && argv [1][0] == '-')
    {
        if (   !strcmp (argv [1], "-f")
            || !strcmp (argv [1], "--force")
           )
        {
            do_write_index = 1;
            filename = argv [2];
        }
        else
        {
            print_usage (argc, argv);
            return -1;
        }
    }
    else
    {
        filename = argv [1];
        do_write_index = 0;
    }

    have_subfiles = 0;
    struct adios_bp_buffer_struct_v1 * b = 0;

    b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    adios_buffer_struct_init (b);

    int flags = O_RDONLY;
    if (do_write_index)
        flags = O_RDWR;

    int fd = open (filename, flags);
    if (fd < 0) {
        fprintf (stderr, "recover: cannot open file %s\n", filename);
        if (errno)
            fprintf (stderr, "%s\n", strerror(errno));
        return -1;
    }

    struct stat statbuf; 
    fstat (fd, &statbuf);
    uint64_t file_size = statbuf.st_size;

    b->f = fd; 

    printf ("File size in bytes: %" PRIu64 "\n", file_size);

    /* Variables to build new index */
    struct adios_index_struct_v1 * index = adios_alloc_index_v1(1);
    //struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_process_group_struct_v1 * pg = 0;
    //struct adios_index_var_struct_v1 * vars_root = 0;
    //struct adios_index_attribute_struct_v1 * attrs_root = 0;


    uint64_t pg_num = 0L;
    uint64_t curr_offset = 0L;
    uint64_t new_offset = 0L;
    uint64_t pgsize_reported; // size of current pg (as indicated in PG header (wrongly))
    uint64_t pgsize_actual; // size of current pg based on processing (accurate)
    int found_pg = 0;

    printf (DIVIDER);
    found_pg = find_pg (fd, new_offset, file_size, &pgsize_reported);

    // pass over the PGs from beginning of file
    while (found_pg)
    {
        curr_offset = new_offset;
        pg_num++;
        printf ("PG %" PRIu64 " found at offset %" PRIu64 "\n", pg_num, curr_offset);


        /* Let's process the group */
        /* Allocate PG index struct here to allow to be used below */
        pg = (struct adios_index_process_group_struct_v1 *)
                malloc (sizeof(struct adios_index_process_group_struct_v1));

        // setup where to read the process group from (and size)
        pg->offset_in_file = curr_offset;
        //b->read_pg_offset = pg->offset_in_file;
        b->read_pg_offset = curr_offset;
        b->read_pg_size = pgsize_reported;

        /* Temporary variables for parsing one PG */
        struct adios_process_group_header_struct_v1 pg_header;
        struct adios_vars_header_struct_v1 vars_header;
        struct adios_attributes_header_struct_v1 attrs_header;

        struct adios_var_header_struct_v1 var_header;
        struct adios_var_payload_struct_v1 var_payload;
        struct adios_attribute_struct_v1 attribute;

        init_dimensions ();  // store scalar values from this PG temporarily

        /* Read the whole PG into a buffer and start parsing */
        adios_posix_read_process_group (b);
        adios_parse_process_group_header_v1 (b, &pg_header);
        print_process_group_header (pg_num, &pg_header);

        add_pg_to_index (index, &pg_header, curr_offset);

        adios_parse_vars_header_v1 (b, &vars_header);
        print_vars_header (&vars_header);

        int i;
        for (i = 0; i < vars_header.count; i++)
        {
            var_payload.payload = 0;

            adios_parse_var_data_header_v1 (b, &var_header);
            print_var_header (&var_header);

            if ( var_header.dims == 0)
            {
                // Load scalars to save them for handling as dimension values
                var_payload.payload = malloc (var_header.payload_size + 1);
                adios_parse_var_data_payload_v1 (b, &var_header, &var_payload
                                                ,var_header.payload_size
                                                );
            }
            else
            {
                // Just parse to move the offset in buffer, don't read data
                adios_parse_var_data_payload_v1 (b, &var_header, NULL, 0);
            }

            store_scalar_dimensions (&var_header, &var_payload);
            add_var_to_index (index, &pg_header, &var_header, &var_payload);

            if (var_payload.payload)
            {
                free (var_payload.payload);
                var_payload.payload = 0;
            }
            //printf ("\n");

        }

        adios_parse_attributes_header_v1 (b, &attrs_header);
        print_attrs_header (&attrs_header);

        for (i = 0; i < attrs_header.count; i++)
        {
            adios_parse_attribute_v1 (b, &attribute);
            //print_attribute (&attribute);
            //printf ("\n");
        }

        pgsize_actual = b->offset;
        printf ("Actual size of group by processing: %" PRIu64 " bytes\n", pgsize_actual);

        pg = pg->next;
        

        printf (DIVIDER);
        found_pg = 0;
        if (curr_offset + pgsize_actual < file_size) 
        {
            new_offset =  curr_offset + pgsize_actual;
            found_pg = find_pg (fd, curr_offset+pgsize_actual, file_size, &pgsize_reported);
        }
        if (!found_pg && 
            pgsize_actual != pgsize_reported &&
            curr_offset + pgsize_reported < file_size) 
        {
            new_offset =  curr_offset + pgsize_reported;
            found_pg = find_pg (fd, curr_offset+pgsize_reported, file_size, &pgsize_reported);
        }
    }

    // The end of the last successfully processed PG
    // This will be the start of the index data
    curr_offset += pgsize_actual;

    printf (DIVIDER);
    printf ("Found %" PRIu64 " PGs to be processable\n", pg_num);

    write_index (fd, curr_offset, index);
    adios_posix_close_internal (b); // will close fd

    return 0;
}


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

    //printf ("Process Group: %" PRIu64 "\n", num);
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
    return;
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

/*
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
*/

/*
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
*/

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
    return;
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
