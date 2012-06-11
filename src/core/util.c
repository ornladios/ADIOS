#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <errno.h>

#include "core/util.h"


static char * remove_whitespace (char *start, char *end) 
{
    char *s = start;
    char *e = end;
    int orig_len = (int) (e-s);
    int final_len;
    char *res;
    // remove front whitespace (but do not go far beyond the end)
    while (s <= e && 
           (*s==' ' || *s=='\t' || *s=='\n')
          ) s++;
    if (s <= e) { // there is some text 
        // remove tail whitespace
        while (s <= e && 
               (*e==' ' || *e=='\t' || *e=='\n')
              ) e--;
        // create result 
        final_len = e - s + 1; //  length of result 
        if (final_len > 0) {
            res = (char *) malloc (final_len + 1); // allocate space s..e and \0
            memcpy(res, s, final_len);
            res[final_len] = 0;
        } else {
            // "   = something" patterns end here
            res = NULL;
        }
    } else {
        // no non-whitespace character found
        res = NULL;
    }
    return res;
}


/* Split a line at = sign into name and value pair
   Remove " ", TAB and Newline from around names and values
   Return NULL for name and value if there is no = sign in line
   Return newly allocated strings otherwise
   Used by: esimmon_internal_text_to_name_value_pairs
 */
static void splitnamevalue (const char * line, int linelen,  char **name, char **value)
{
    char *equal; // position of first = sign in line

    equal = strchr (line, '=');
    if (equal && equal != line) {
        /* 1. name */
        // from first char to before =
        *name = remove_whitespace ((char*)line, equal-1);
        //printf ("      --name=[%s]\n", *name);
        /* 2. value */
        // from after = to the last character of line
        *value = remove_whitespace (equal+1, (char*)line+linelen-1);
        //printf ("      --value=[%s]\n", *value);

    } else if (equal != line) {
        /* check if it as name without = value statement */
        *name = remove_whitespace ((char*)line, (char*)line+linelen-1);
        //printf ("      --name only=[%s]\n", *name);
        *value = NULL;
    } else { 
        // funny text starting with =. E.g. "=value" 
        *name = NULL;
        *value = NULL;
    }
}

PairStruct * text_to_name_value_pairs (const char * text)
{
    /* Process a multi-line and/or ;-separated text and create a list
       of name=value pairs from each line which has a 
           name = value 
       pattern. Whitespaces are removed. 
         "X = 1
          Y = 2"  
       is not valid because of missing ';', but
          "X=1; Y=5;
          Z=apple"  
       is valid
    */
    char *name, *value; 
    char *item, *delim;
    int len;
    char line[256];
    PairStruct *res = NULL, *last = NULL, *pair;

    if (!text) return res;

    item  = text; 
    while (item) {
        delim = strchr (item, ';');
        if (delim) 
            len = (int) (delim-item); 
        else 
            len = strlen (item);

        strncpy (line, item, len);
        line[len] = '\0';

        //printf ("    --Line=[%s]\n", line);
        splitnamevalue(line, len, &name, &value);
        if (name) {
            pair = (PairStruct *) malloc (sizeof(PairStruct));
            pair->name = name;
            pair->value = value;
            pair->next = NULL;
            if (last) {
                last->next = pair;
                last = pair;
            } else {
                res = pair; 
                last = pair;
            }
        }
        if (delim && delim+1 != 0)
            item = delim+1;
        else
            item = NULL;
    }
    return res;
}


void free_name_value_pairs (PairStruct * pairs)
{
    PairStruct *p;
    while (pairs) {
        free(pairs->name);
        free(pairs->value);
        p = pairs;
        pairs=pairs->next;
        free(p);
    }
}

/* Reverse the order in an array in place.
   use swapping from Fortran/column-major order to ADIOS-read-api/C/row-major order and back
*/
void swap_order(int n, uint64_t *array, int *timedim)
{
    int i;
    uint64_t tmp;
    for (i=0; i<n/2; i++) {
        tmp = array[i];
        array[i] = array[n-1-i];
        array[n-1-i] = tmp;
    }
    if (*timedim > -1)
        *timedim = (n-1) - *timedim; // swap the time dimension too
}

/* Change endianness of each element in an array */
/* input: array, size in bytes(!), size of one element */
void change_endianness( void *data, uint64_t slice_size, enum ADIOS_DATATYPES type)
{
    int size_of_type = bp_get_type_size(type, "");
    uint64_t n = slice_size / size_of_type;
    uint64_t i;
    char *ptr = (char *) data;

    if (slice_size % size_of_type != 0) {
        fprintf(stderr, "Adios error in bp_utils.c:change_endianness(): "
                    "An array's endianness is to be converted but the size of array "
                    "is not dividable by the size of the elements: "
                    "size = %lld, element size = %d\n", slice_size, size_of_type);
    }

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
            switch (size_of_type) {
                /* case 1: nothing to do */
                case 2:
                    for (i=0; i < n; i++) {
                        swap_16_ptr(ptr);
                        ptr += size_of_type;
                    }
                    break;
                case 4:
                    for (i=0; i < n; i++) {
                        swap_32_ptr(ptr);
                        ptr += size_of_type;
                    }
                    break;
                case 8:
                    for (i=0; i < n; i++) {
                        swap_64_ptr(ptr);
                        ptr += size_of_type;
                    }
                    break;
                case 16:
                    for (i=0; i < n; i++) {
                        swap_128_ptr(ptr);
                        ptr += size_of_type;
                    }
                    break;
            }
            break;

        case adios_complex:
            for (i=0; i < n; i++) {
                swap_32_ptr(ptr);   // swap REAL part 4 bytes 
                swap_32_ptr(ptr+4); // swap IMG part 4 bytes
                ptr += size_of_type;
            }
            break;

        case adios_double_complex:
            for (i=0; i < n; i++) {
                swap_64_ptr(ptr);   // swap REAL part 8 bytes 
                swap_64_ptr(ptr+8); // swap IMG part 8 bytes
                ptr += size_of_type;
            }
            break;

        case adios_string:
        default:
            /* nothing to do */
            break;
    }
}

void copy_data (void *dst, void *src,
        int idim,
        int ndim,
        uint64_t* size_in_dset,
        uint64_t* ldims,
        const uint64_t * readsize,
        uint64_t dst_stride,
        uint64_t src_stride,
        uint64_t dst_offset,
        uint64_t src_offset,
        uint64_t ele_num,
        int      size_of_type
        )
{
    unsigned int i, j;
    uint64_t dst_offset_new=0;
    uint64_t src_offset_new=0;
    uint64_t src_step, dst_step;
    if (ndim-1==idim) {
        for (i=0;i<size_in_dset[idim];i++) {
            memcpy ((char *)dst + (i*dst_stride+dst_offset)*size_of_type,
                    (char *)src + (i*src_stride+src_offset)*size_of_type,
                    ele_num*size_of_type);
        }
        return;
    }

    for (i = 0; i<size_in_dset[idim];i++) {
        // get the different step granularity 
        // for each different reading pattern broke
        src_step = 1;
        dst_step = 1;
        for (j = idim+1; j <= ndim-1;j++) {
            src_step *= ldims[j];
            dst_step *= readsize[j];
        }
        src_offset_new =src_offset + i * src_stride * src_step;
        dst_offset_new = dst_offset + i * dst_stride * dst_step;
        copy_data ( dst, src, idim+1, ndim, size_in_dset,
                ldims,readsize,
                dst_stride, src_stride,
                dst_offset_new, src_offset_new,
                ele_num, size_of_type);
    }
}

void alloc_namelist (char ***namelist, int length)
{
        int j;

        *namelist = (char **) malloc(length*sizeof(char*));
        for (j=0;j<length;j++)
                (*namelist)[j] = (char *) malloc(255);

        return;
}

void free_namelist (char **namelist, int length)
{
        int i;
        if (namelist) {
                for (i=0;i<length;i++) {
                        if(namelist[i])
                                free(namelist[i]);
                }
                free(namelist);
        }
        return;
}

void list_insert_reader (read_request ** h, read_request * q)
{
    read_request * head;
    if (!h || !q)
    {
        printf ("Error: list_insert_reader ()\n");
        return;
    }

    head = * h;
    if (!head)
    {
        * h = q;
        q->next = NULL;

        return;
    }

    while (head->next)
    {
        head = head->next;
    }

    head->next = q;
    q->next = NULL;

    return;
}

ADIOS_SELECTION * copy_selection (ADIOS_SELECTION * sel)
{
    ADIOS_SELECTION * nsel;

    nsel = (ADIOS_SELECTION *) malloc (sizeof (ADIOS_SELECTION));
    assert (nsel);

    nsel->type = sel->type;

    if (sel->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        nsel->u.bb.ndim = sel->u.bb.ndim;
        nsel->u.bb.start = (uint64_t *) malloc (sel->u.bb.ndim * 8);
        nsel->u.bb.count = (uint64_t *) malloc (sel->u.bb.ndim * 8);
        assert (nsel->u.bb.start && nsel->u.bb.count);

        memcpy (nsel->u.bb.start, sel->u.bb.start, sel->u.bb.ndim * 8);
        memcpy (nsel->u.bb.count, sel->u.bb.count, sel->u.bb.ndim * 8);
    }
    else if (sel->type == ADIOS_SELECTION_POINTS)
    {
        nsel->u.points.ndim = sel->u.points.ndim;
        nsel->u.points.npoints = sel->u.points.npoints;
        nsel->u.points.points = (uint64_t *) malloc (sel->u.points.npoints * 8);
        assert (nsel->u.points.points);

        memcpy (nsel->u.points.points, sel->u.points.points, sel->u.points.npoints * 8);
    }
    else if (sel->type == ADIOS_SELECTION_WRITEBLOCK)
    {
        nsel->u.block.index = sel->u.block.index;
    }
    else if (sel->type == ADIOS_SELECTION_AUTO)
    {
        //FIXME
    }
    else
    {
        //FIXME
    }

    return nsel;
}

