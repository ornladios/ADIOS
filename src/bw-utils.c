#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include "bw-utils.h"

void bw_fopen_ (char * filename, long long * file_hd)
{
    int fbp = -1;
    int pos = 10;

    fbp = open (filename, O_RDWR);
    if (fbp == -1)
    {
        //fprintf (stderr, "File does not exisit!\n");
        fbp = open (filename, O_CREAT | O_RDWR
                   ,  S_IRUSR | S_IWUSR
                    | S_IRGRP | S_IWGRP
                    | S_IROTH | S_IWOTH
                   );

    }
    lseek (fbp, 0, SEEK_END);
    off_t size = lseek (fbp, 0, SEEK_CUR);
    if (size > 0)
    {
        lseek (fbp, size - 4, SEEK_SET);
        read (fbp, &pos, sizeof (int) * 1);
        if (pos == 0)
            lseek (fbp, size - 4, SEEK_SET);
        else
            return;
    }
    *file_hd = (long long) fbp;
}

void bw_fclose_ (long long * fptr)
{
    if (*fptr)
    {
        close (*fptr);
        *fptr = (long long) -1;
    }
}

void bw_fwrite_ (long long * fptr, void * buffer, size_t * number)
{
    write (*(int *) fptr, buffer, sizeof (char) * *number);
}

// write to a buffer
void bwrite (void * val, int size, int number, void  * ptrbuf, int * ptridx)
{
    memcpy (ptrbuf + *ptridx, val, size * number);
    *ptridx = *ptridx + size * number;
}

// write to a file directly
void bwrite_f (void * val, int size, int number, void * fptr, int * ptridx)
{
    size_t count = size * number;
    bw_fwrite_ (fptr, val, &count);
}

// we use this to decide how to write (buffer or file);
typedef void (* WRITE_FN) (void *, int, int, void *, int *);

static WRITE_FN BWRITE = bwrite;

void bw_set_write (int a)
{
    if (a == 1)
        BWRITE = bwrite;
    else
        BWRITE = bwrite_f;
}

// rank = number of dimension entries
// dims = dimension entries (5 ints wide):
//            lower bound, upper bound, use lower bound, use upper bound, stride
// g_dims = global dimension entries (5 ints wide): 
//            lower bound, upper bound, use lower bound, use upper bound, stride
void bw_dset_ (void * buf, int * startidx, int * endidx, char * str, char * name
              ,void * val, int * type, int * rank, int * dims
              )
{
    bw_dset (buf, *startidx, endidx, str, name, val
            ,(enum vartype_t) *type, *rank
            ,(struct adios_bp_dimension_struct *) dims
            );    
}

void bw_scalar_ (void * buf, int * startidx, int * endidx, char * str, char * name
                ,void * val, int * type
                )
{
    bw_scalar (buf, *startidx, endidx, str, name, val, (enum vartype_t) *type);
}

void bw_an_d_ (void * buf, int * startidx, int * endidx, char * str, char * aname
              ,void * aval, int * type
              )
{
    bw_attr_num_ds (buf, *startidx, endidx, str, aname, aval
                   ,(enum vartype_t) *type
                   );
}

void bw_an_g_ (void * buf, int * startidx, int * endidx, char * str, char * aname
              ,void * aval, int * type
              )
{
    bw_attr_num_gp (buf, *startidx, endidx, str, aname, aval
                   ,(enum vartype_t) *type
                   );
}

void bw_as_d_ (void * buf, int * startidx, int * endidx, char * str, char * aname
              ,char * aval
              )
{
    bw_attr_str_ds (buf, *startidx, endidx, str, aname, aval);    
}

void bw_as_g_ (void * buf, int * startidx, int * endidx, char * str, char * aname
              ,char * aval
              )
{
    bw_attr_str_gp (buf, *startidx, endidx, str, aname, aval);
}

// rank = number of dimension entries
void bw_dset (void * fbp, int startidx, int * endidx, char * path, char * name
             ,void * val, enum vartype_t type, int rank
             ,struct adios_bp_dimension_struct * dims
             )
{
    int var_step;
    int localidx = startidx;
    enum TAG_t tag;

    if ((*(int *) fbp) == -1)
        fprintf (stderr, "Failed to open a file!\n");

    var_step = bcalsize_dset (path, name, type, rank, dims);
    *endidx = startidx + var_step;
    BWRITE (&var_step, sizeof (int), 1, fbp, &localidx);

    tag = DST_TAG;
    bw_stringtag (fbp, &localidx, tag, name);

    tag = DIR_TAG;
    bw_stringtag (fbp, &localidx, tag, path);

    tag = DSTVAL_TAG;
    bw_dsettag (fbp, &localidx, tag, val, type, rank, dims);
}

void bw_scalar (void * fbp, int startidx, int * endidx, char * str, char * name
               ,void * val, enum vartype_t type
               )
{
    int var_step;
    int localidx = startidx;
    enum TAG_t tag;

    if ((*(int *) fbp) == -1)
        fprintf (stderr, "Failed to open a file!\n");

    var_step = bcalsize_scalar (str, name, type, val);
    *endidx = startidx + var_step;
    BWRITE (&var_step, sizeof (int), 1, fbp, &localidx);

    tag = SCR_TAG;
    bw_stringtag (fbp, &localidx, tag, name);

    tag = DIR_TAG;
    bw_stringtag (fbp, &localidx, tag, str);
    
    tag = VAL_TAG;
    bw_scalartag (fbp, &localidx, tag, val, type);
}

void bw_attr_num_ds (void * fbp, int startidx, int * endidx, char * str
                    ,char * aname, void * val, enum vartype_t type
                    )
{
    int var_step;
    enum TAG_t tag;
    int localidx = startidx;    

    var_step = bcalsize_attr_num (str, aname, type, val);
    *endidx = startidx + var_step;
    BWRITE (&var_step, sizeof (int), 1, fbp, &localidx);

    tag = DSTATRN_TAG;
    bw_stringtag (fbp, &localidx, tag, aname);

    tag = DIR_TAG;
    bw_stringtag (fbp, &localidx, tag, str);
    
    tag = VAL_TAG;
    bw_scalartag (fbp, &localidx, tag, val, type);
}

void bw_attr_num_gp (void * fbp, int startidx, int * endidx, char * str
                    ,char * aname, void * val, enum vartype_t type
                    )
{
    int var_step;
    enum TAG_t tag;
    int localidx = startidx;

    var_step = bcalsize_attr_num (str, aname, type, val);
    *endidx = startidx + var_step;
    BWRITE (&var_step, sizeof (int), 1, fbp, &localidx);

    tag = GRPATRN_TAG;
    bw_stringtag (fbp, &localidx, tag, aname);

    tag = DIR_TAG;
    bw_stringtag (fbp, &localidx, tag, str);
    
    tag = VAL_TAG;
    bw_scalartag (fbp, &localidx, tag, val, type);
}

void bw_attr_str_gp (void * fbp, int startidx, int * endidx, char * str
                    ,char * aname, char * aval
                    )
{
    int var_step;
    enum TAG_t tag;
    enum vartype_t type = bp_string;
    int localidx = startidx;    

    var_step = bcalsize_attr_str (str, aname, aval);
    *endidx = startidx + var_step;
    BWRITE (&var_step, sizeof (int), 1, fbp, &localidx);

    tag = GRPATRS_TAG;
    bw_stringtag (fbp, &localidx, tag, aname);

    tag = DIR_TAG;
    bw_stringtag (fbp, &localidx, tag, str);
    
    tag = VAL_TAG;
    bw_scalartag (fbp, &localidx, tag, aval, type);
}

void bw_attr_str_ds (void * fbp, int startidx, int * endidx, char * str
                    ,char * aname, char * aval
                    )
{
    int var_step;
    enum TAG_t tag;
    int localidx = startidx;

    var_step = bcalsize_attr_str (str, aname, aval);
    *endidx = startidx + var_step;
    BWRITE (&var_step, sizeof (int), 1, fbp, &localidx);

    tag = DSTATRS_TAG;
    bw_stringtag (fbp, &localidx, tag, aname);

    tag = DIR_TAG;
    bw_stringtag (fbp, &localidx, tag, str);
    
    tag = VAL_TAG;
    bw_scalartag (fbp, &localidx, tag, aval, bp_string);
}

int bp_calsize_stringtag (char * name)
{
    return  sizeof (int)
          + sizeof (int)
          + strlen (name);
}

void bw_stringtag (void * fbp, int * ptridx, enum TAG_t tag, char * name)
{ 
    int size = strlen (name);

    BWRITE (&tag, sizeof (int), 1, fbp, ptridx);
    BWRITE (&size, sizeof (int), 1, fbp, ptridx);
    BWRITE (name, sizeof (char), size, fbp, ptridx);
}

int bp_calsize_scalartag (enum vartype_t type, void * val)
{
    int unit_size = 0;

    unit_size = bp_getsize (type, val);
    if (unit_size < 0)
    {
        fprintf (stderr, "Invalid type %s.  Defaulting to 4 bytes\n"
                ,adios_type_to_string (type)
                );

        unit_size = 4;
    }

    return  sizeof (int)
          + sizeof (int)
          + sizeof (int)
          + unit_size * sizeof (char);
}

void bw_scalartag (void * fbp, int * ptridx, enum TAG_t tag, void * val
                  ,enum vartype_t type
                  )
{
    int size;

    size = bp_getsize (type, val);
    if (size < 0)
    {
        fprintf (stderr, "Invalid type %s.  Defaulting to 4 bytes\n"
                ,adios_type_to_string (type)
                );

        size = 4;
    }

    BWRITE (&tag, sizeof (int), 1, fbp, ptridx);
    BWRITE (&size, sizeof (int), 1, fbp, ptridx);
    BWRITE (&type, sizeof (int), 1, fbp, ptridx);
    BWRITE (val, sizeof (char), size, fbp, ptridx);
}

// rank = number of dimension entries
int bp_calsize_dsettag (enum vartype_t type, int rank
                       ,struct adios_bp_dimension_struct * dims
                       )
{
    char * tmp = "";  // dummy since we don't expect strings in arrays
    int subsize = bp_getsize (type, tmp);
    unsigned long long use_count = 0;
    unsigned long long total_count = 0;

    if (subsize < 0)
    {
        fprintf (stderr, "Invalid type %s.  Defaulting to 4 bytes\n"
                ,adios_type_to_string (type)
                );

        subsize = 4;
    }

    adios_var_element_count (rank, dims, &use_count, &total_count);
    subsize *= use_count;

    return  sizeof (int)
          + sizeof (int)
          + sizeof (int)
          + sizeof (struct adios_bp_dimension_struct) * rank
          + sizeof (int)
          + sizeof (char) * subsize;
}

// rank = number of dimension entries
void bw_dsettag (void * fbp, int * ptridx, enum TAG_t tag, void * val
                ,enum vartype_t type, int rank
                ,struct adios_bp_dimension_struct * dims
                )
{
    int total_size;
    int element_size;
    unsigned long long use_count = 0;
    unsigned long long total_count = 0;
    int full_write = 1;

    element_size = bp_getsize (type, val);
    if (element_size < 0)
    {
        fprintf (stderr, "Invalid type %s.  Defaulting to 4 bytes\n"
                ,adios_type_to_string (type)
                );

        element_size = 4;
    }

    adios_var_element_count (rank, dims, &use_count, &total_count);

    if (use_count == total_count)
        full_write = 1;
    else
        full_write = 0;

    total_size = use_count * element_size;

    BWRITE (&tag, sizeof (int), 1, fbp, ptridx);
    BWRITE (&total_size, sizeof (int), 1, fbp, ptridx);
    BWRITE (&rank, sizeof (int), 1, fbp, ptridx);
    BWRITE (dims, sizeof (struct adios_bp_dimension_struct), rank, fbp, ptridx);
    BWRITE (&type, sizeof (int), 1, fbp, ptridx);
    BWRITE (val, sizeof (char), total_size, fbp, ptridx);
}

// rank = number of dimension entries
int bcalsize_dset (char * dir_name, char * name, enum vartype_t type, int rank
                  ,struct adios_bp_dimension_struct * dims
                  )
{
    int size = 4;                           // overall element size

    size += bp_calsize_stringtag (name);
    size += bp_calsize_stringtag (dir_name);
    size += bp_calsize_dsettag (type, rank, dims);

    return size;
}

int bw_calsize_dset_ (int * subsize, char * dir_name, char * name
                     ,enum vartype_t * type, int * rank, int * dims
                     )
{
    *subsize = bcalsize_dset (dir_name, name, *type, *rank
                             ,(struct adios_bp_dimension_struct *) dims
                             );

    return *subsize;
}

int bcalsize_scalar (char *dir_name, char* name,enum vartype_t type,void * val)
{
    int size = 4;                          // overall element size

    size += bp_calsize_stringtag (name);
    size += bp_calsize_stringtag (dir_name);
    size += bp_calsize_scalartag (type, val);

    return size;
}

void bw_calsize_scalar_ (int *subsize,char *dir_name, char* name
                        ,enum vartype_t *type, void * val
                        )
{
    *subsize = bcalsize_scalar (dir_name, name, *type, val);
}

int bcalsize_attr_str (char * dir_name, char * aname, char * aval)
{
    int size = 4;                       // overall element size

    size += bp_calsize_stringtag (aname);
    size += bp_calsize_stringtag (dir_name);
    size += bp_calsize_scalartag (bp_string, aval);

    return size;
}

int bcalsize_attr_num (char * dir_name, char * aname, enum vartype_t type
                      ,void * val
                      )
{
    int size = 4;                           // overall element size

    size += bp_calsize_stringtag (aname);
    size += bp_calsize_stringtag (dir_name);
    size += bp_calsize_scalartag (type, val);

    return size;
}
