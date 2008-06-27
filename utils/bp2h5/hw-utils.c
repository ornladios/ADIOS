#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include "hw-utils.h"
#include "br-utils.h"
#include "binpack-utils.h"
#define  MAX_RANK 10 * 5
#define STRLEN 1000 

/*
 * Global config variables to set the convenstion of conversion
 */

/*
 * array to dataset convention (transpose array dimension order if using Fortran)
 */
enum lang_convention array_dim_order_fortran = USE_FORTRAN;

/*
 * string-typed dataset attribute write convention
 */
enum lang_convention attr_str_ds_fortran = USE_FORTRAN;

/*
 * string-typed group attribute write convention
 */
enum lang_convention attr_str_gp_fortran = USE_FORTRAN;

/* 
 * scalar var/attribute write convention
 */
enum scalar_convention scalar_dataspace = USE_SCALAR;

/*
 * initialization flag
 */
int bp2h5_initialized = 0;

/*
 * verbose level flag
 */
enum verbose_level verbose = NO_INFO;

/*
 * buffer used to read value of scalars or datasets from bp file
 */
struct buffer_struct
{
    unsigned long long DATALEN;
    void * val;
} read_buffer;


/*
 * Initialization
 * initialize_bp2h5() sets config variables and allocate
 * internal read buffer.
 * It returns 0 if no error and -1 otherwise
 */
int initialize_bp2h5(enum lang_convention array_dim_order,
                     enum lang_convention attr_str_ds,
                     enum lang_convention attr_str_gp,
                     enum scalar_convention scalar_ds,
                     unsigned long long buf_size,
                     enum verbose_level verb
                     )
{
    if(!bp2h5_initialized) {
        array_dim_order_fortran = array_dim_order;
        attr_str_ds_fortran = attr_str_ds;
        attr_str_gp_fortran = attr_str_gp;
        scalar_dataspace = scalar_ds;
        verbose = verb;
        
        // allocate initial buffer size
        if(buf_size <= 0) {
            fprintf (stderr, "Error in initialization: invalid buffer size\n");
            return -1;            
        }
        read_buffer.DATALEN = buf_size;
        read_buffer.val = (char *) malloc (read_buffer.DATALEN);
        if (!read_buffer.val)
        {
            fprintf (stderr, "Error in initialization: cannot allocate memory for read buffer\n");
            return -1;
        }

        bp2h5_initialized = 1;
    }

    return 0;
}

/*
 * Helper function for buffer reallocation
 * code borrowed from GNU Libc Manual.
 */
void *xrealloc (void *ptr, size_t size)
{
    register void *value = realloc (ptr, size);
    if (value == 0) {
        fprintf(stderr, "Virtual memory exhausted\n");
        exit(-1);
    }
    return value;
}

/*
 * Allocate buffer if necessary to accommodate the data value or dataset being read 
 */
static void pre_element_fetch (struct adios_bp_element_struct * element
                              ,void ** buffer, long long * buffer_size
                              ,void * private_data
                              )
{
    struct buffer_struct * data = (struct buffer_struct *) private_data;

    // check if the buffer size is big enough
    if(data->DATALEN < element->size) {
        data->val = xrealloc(data->val, element->size);
        data->DATALEN = element->size;
    }

    element->data = data->val;
    *buffer = data->val;
    *buffer_size = data->DATALEN;
}

void hw_free2D (char ** grp_name, int level)
{
    int i;

    for (i = 0; i < level; i++)
    {
        free (grp_name [i]);
        grp_name [i] = NULL;
    }

    free (grp_name);
}

/*
 * Convert bp file fnamein to a h5 file fnameout
 */
int hw_makeh5 (char * fnamein, char * fnameout)
{
    long long handle = 0;
    char * tmpstr;
    int size;
    struct adios_bp_element_struct * element = NULL;
    long long element_size = 0;

    handle = br_fopen (fnamein);
    if (!handle)
    {
        fprintf (stderr, "Error in open bp file: cannot open %s file!\n", fnamein);

        return -1;
    }

    // initialization
    if(!bp2h5_initialized && 
        (initialize_bp2h5(USE_FORTRAN, USE_FORTRAN, USE_FORTRAN, USE_SCALAR, 100 * 1024 * 1024, NO_INFO))) {
        return -1;
    }

    // xxx.bp --> xxx.h5
    tmpstr = strdup (fnamein);
    size = strlen (fnamein);
    tmpstr [size - 2] = 'h';
    tmpstr [size - 1] = '5';
    if(!fnameout) {
        fnameout=tmpstr;
    }

    hid_t h5file_id;
    hid_t root_id;
    h5file_id = H5Fcreate (fnameout, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    root_id = H5Gopen (h5file_id, "/");

    // Use br_get_next_element_specific instead of br_get_next_element_general
    // in order to handle arbitrary-sized data. 
    while (element_size = br_get_next_element_specific (handle
                                                       ,pre_element_fetch
                                                       ,0
                                                       ,&read_buffer
                                                       ,&element
                                                       )
          )
    {
        if(verbose >= DEBUG_INFO) {
            fprintf(stderr, "Write element %-20s : Tag: %-10s Size: %lu Bytes\n"
                   ,element->name 
                   ,adios_tag_to_string(element->tag)
                   ,element->size
                   );
        }
        switch (element->tag)
        {
            case DSTATRS_TAG:
                hw_attr_str_ds (root_id, element->path, element->name, read_buffer.val);
                break;

            case DSTATRN_TAG:
                hw_attr_num_ds (root_id, element->path, element->name, read_buffer.val
                               ,element->type
                               );
                break;

            case GRPATRS_TAG:
                hw_attr_str_gp (root_id, element->path, element->name, read_buffer.val);
                break;

            case GRPATRN_TAG:
                hw_attr_num_gp (root_id, element->path, element->name, read_buffer.val
                               ,element->type
                               );
                break;

            case SCR_TAG:
                hw_scalar (root_id, element->path, element->name, read_buffer.val
                          ,element->type, 0
                          );
                break;

            case DST_TAG:
                if(verbose >= DEBUG_INFO) { 
                    fprintf (stderr, "  Ranks: %u\n", element->ranks);
                    if (element->dims) {
                        for (int i = 0; i < element->ranks; i++) {                                      
                            fprintf(stderr, "  Dim(%d): %d(%d):%d\n"
                                   ,i
                                   ,element->dims [i].local_bound
                                   ,element->dims [i].global_bound
                                   ,element->dims [i].global_offset
                                   );
                        }
                    }
                }

                hw_dset (root_id, element->path, element->name, read_buffer.val
                        ,element->type, element->ranks, element->dims
                        );
                break;

            default:
                break;
        }

        br_free_element (element);
    }

    H5Gclose (root_id);
    H5Fclose (h5file_id);
    br_fclose (handle);

    free (read_buffer.val);

    return 0;
}

/*
 * Write dataset to h5 file
 */
void hw_dset (hid_t root_id, char * dirstr, char * name, void * data
             ,enum vartype_t type, int rank
             ,struct adios_bp_dimension_struct * dims
             )
{
    H5Eset_auto (NULL,NULL);
    char ** grp_name;
    int level;
    int i;
    hid_t grp_id [NUM_GP + 1];
    
    grp_id [0] = root_id;
    
    grp_name = bp_dirparser (dirstr, &level);
    
    for (i = 0; i < level; i++)
    {
        grp_id [i + 1] = H5Gopen (grp_id [i],grp_name [i]);
        if (grp_id [i + 1] < 0)
        {
            grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
        }
    }

    if (dims [0].global_bound)
    {
        hsize_t * global_h5dims;
        hsize_t * local_h5dims;
        hsize_t * start;
        hsize_t * stride;
        hsize_t * count;

        hid_t dataspace;
        hid_t memspace;
        hid_t dataset;
        hid_t type_id;

        herr_t h5_status;

        h5_status = bp_getH5TypeId (type, &type_id, data);

        global_h5dims = (hsize_t *) malloc (rank * sizeof (hsize_t));
        local_h5dims = (hsize_t *) malloc (rank * sizeof (hsize_t));
        start = (hsize_t *) malloc (rank * sizeof (hsize_t));
        stride = (hsize_t *) malloc (rank * sizeof (hsize_t));
        count = (hsize_t *) malloc (rank * sizeof (hsize_t));

        if(array_dim_order_fortran == USE_FORTRAN) { 
            // transpose dimension order for Fortran arrays
            for (int i = 0; i < rank; i++)
            {
                int reverse_i = rank - 1 - i;

                global_h5dims [reverse_i] = dims [i].global_bound;

                local_h5dims [reverse_i] = dims [i].local_bound;

                start [reverse_i] = dims [i].global_offset;

                stride [reverse_i] = 1;

                count [reverse_i] = dims [i].local_bound;
            }
        }
        else if(array_dim_order_fortran == USE_C) {
            for (int i = 0; i < rank; i++)
            {
                global_h5dims [i] = dims [i].global_bound;

                local_h5dims [i] = dims [i].local_bound;

                start [i] = dims [i].global_offset;

                stride [i] = 1;

                count [i] = dims [i].local_bound;
            }
        }

        if(verbose >= DEBUG_INFO) {
            for (int i = 0; i < rank; i++)
                fprintf(stderr, "  Hyperslab index = %d local bound = %d global bound = %d start = %d\n"
                       ,i
                       ,local_h5dims[i]
                       ,global_h5dims[i]
                       ,start[i]
                       );
        }

        dataspace = H5Screate_simple (rank, global_h5dims, NULL);

        h5_status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET
                                  ,start, stride, count, NULL
                                  );

        memspace = H5Screate_simple (rank, local_h5dims, NULL);

        dataset = H5Dopen (grp_id [level], name);
        if (dataset < 0)
        {
            dataset = H5Dcreate (grp_id [level], name, type_id
                                ,dataspace, H5P_DEFAULT
                                );
        }

        h5_status = H5Dwrite (dataset, type_id, memspace, dataspace
                             ,H5P_DEFAULT, data
                             );

        H5Sclose (memspace);
        H5Sclose (dataspace);
        H5Dclose (dataset);

        free (global_h5dims);
        free (local_h5dims);
        free (start);
        free (stride);
        free (count);
    }
    else
    {
        hsize_t * h5dims;
        h5dims = (hsize_t *) malloc (rank * sizeof (hsize_t));    
        if(array_dim_order_fortran == USE_FORTRAN) { 
            // transpose dimension order for Fortran arrays
            for (int i = 0; i < rank; i++)
            {
                h5dims [rank-1-i] = dims [i].local_bound;
            }
        }
        else if(array_dim_order_fortran == USE_C) {
            for (int i = 0; i < rank; i++)
            {
                h5dims [i] = dims [i].local_bound;
            }
        }

        if(verbose >= DEBUG_INFO) {
            for (int i = 0; i < rank; i++) {
                fprintf(stderr, "  Dataspace index = %d dimension = %d \n"
                       ,i
                       ,h5dims[i]
                       );
            }
        }
        
        hw_dataset (grp_id [level], name, data, type, rank, h5dims);
        for (i = 1; i < level + 1; i++)
            H5Gclose (grp_id [i]);

        free (h5dims);
    }
    hw_free2D (grp_name, level);
}

/*
 * Write a scalar var to h5 file
 */
void hw_scalar (hid_t root_id, char * dirstr, char * name, void * data
               ,enum vartype_t type, int append
               )
{
    H5Eset_auto (NULL, NULL);
    char ** grp_name;
    int level;
    int i;
    int ndims = 1;
    hid_t grp_id [NUM_GP + 1];
    hsize_t dims [] = {1};
    
    grp_id [0] = root_id;
    grp_name = bp_dirparser (dirstr, &level);
    
    for (i = 0; i < level; i++)
    {
        grp_id [i + 1] = H5Gopen (grp_id [i], grp_name [i]);
        if (grp_id [i + 1] < 0)
            grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
    }

    if(scalar_dataspace == USE_SCALAR) {
        // write in a scalar
        hw_scalar_as_scalar (grp_id [level], name, data, type, append);
    }
    else if(scalar_dataspace == USE_SINGLE_ELE_ARRAY) {
        // write in a single-element array
        hw_scalar_as_array (grp_id [level], name, data, type, ndims, dims, append);
    }

    for (i = 1; i < level + 1; i++)
        H5Gclose (grp_id [i]);

    hw_free2D (grp_name, level);
}

/*
 * Write a scalar var to h5 file as a scalar
 */
void hw_scalar_as_scalar (hid_t parent_id, char * name, void * data
                 ,enum vartype_t type, int append
                 )
{
    hid_t h5_dataset_id;
    hid_t h5_dataspace_id;
    hid_t h5_type_id;
    int status;
    herr_t h5_status;

    status = bp_getH5TypeId (type, &h5_type_id, data);
    if(status == 0 && h5_type_id > 0)
    {
        h5_dataspace_id = H5Screate(H5S_SCALAR);
        if (h5_dataspace_id > 0) 
        {
            h5_dataset_id = H5Dcreate (parent_id, name, h5_type_id
                                      ,h5_dataspace_id, H5P_DEFAULT
                                      );
            if (h5_dataset_id > 0)
            {
                
                status = H5Dwrite (h5_dataset_id, h5_type_id, H5S_ALL
                                  ,H5S_ALL, H5P_DEFAULT, data
                                  );
                status = H5Dclose (h5_dataset_id);
            }
            else
            {
                if (append)
                {
                    h5_dataset_id = H5Dopen(parent_id, name);
                    if (h5_dataset_id > 0)
                    {
                        h5_status = H5Dwrite (h5_dataset_id, h5_type_id, H5S_ALL
                                             ,H5S_ALL, H5P_DEFAULT, data
                                             );
                        status = H5Sclose (h5_dataset_id);
                    }
                }
            }
            status = H5Sclose (h5_dataspace_id);
        }
        H5Tclose (h5_type_id);
    }

    return;
}

/*
 * Write a scalar var to h5 file as a single-element array
 */
void hw_scalar_as_array (hid_t parent_id, char * name, void * data
                 ,enum vartype_t type, int ndims, hsize_t * dims
                 ,int append
                 )
{
    hid_t h5_dataset_id;
    hid_t h5_dataspace_id;
    hid_t h5_type_id;
    int status;
    herr_t h5_status;

    status = bp_getH5TypeId (type, &h5_type_id, data);
    if(status == 0 && h5_type_id > 0)
    {
        h5_dataspace_id = H5Screate_simple (ndims, dims, NULL);
        if (h5_dataspace_id > 0) 
        {
            h5_dataset_id = H5Dcreate (parent_id, name, h5_type_id
                                      ,h5_dataspace_id, H5P_DEFAULT
                                      );
            if (h5_dataset_id > 0)
            {
                
                status = H5Dwrite (h5_dataset_id, h5_type_id, H5S_ALL
                                  ,H5S_ALL, H5P_DEFAULT, data
                                  );
                status = H5Dclose (h5_dataset_id);
            }
            else
            {
                if (append)
                {
                    h5_dataset_id = H5Dopen(parent_id, name);
                    if (h5_dataset_id > 0)
                    {
                        h5_status = H5Dwrite (h5_dataset_id, h5_type_id, H5S_ALL
                                             ,H5S_ALL, H5P_DEFAULT, data
                                             );
                        status = H5Sclose (h5_dataset_id);
                    }
                }
            }
            status = H5Sclose (h5_dataspace_id);
        }
        H5Tclose (h5_type_id);
    }

    return;
}

/*
 * hw_attr_str_gp() writes string-typed group attribute 
 * hw_attr_str_gp() finds the group to which the attribute is attached and
 * writes the attribute. 
 */
void hw_attr_str_gp (hid_t root_id, char * dirstr, char * aname, char * aval)
{
    H5Eset_auto (NULL, NULL);
    char ** grp_name;
    int level;
    int i;
    hid_t grp_id [NUM_GP + 1];

    grp_id [0] = root_id;
    grp_name = bp_dirparser (dirstr, &level);
    
    for (i = 0; i < level; i++)
    {
        grp_id [i + 1] = H5Gopen (grp_id [i], grp_name [i]);
        if (grp_id [i + 1] < 0)
        {
            grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
        }
    }

    hw_string_attr_gp_internal (grp_id [level], aname, aval);

    for (i = 1; i < level; i++)
        H5Gclose (grp_id [i]);

    hw_free2D (grp_name, level);
}

/*
 * hw_attr_str_ds() writes string-typed dataset attribute 
 * hw_attr_str_ds() finds the dataset to which the attribute is attached and
 * writes the attribute. 
 */
void hw_attr_str_ds (hid_t root_id, char * dirstr, char * aname, char * aval)
{
    H5Eset_auto (NULL, NULL);
    char ** grp_name;
    int level;
    int i;
    hid_t grp_id [NUM_GP + 1];
    hid_t dataset_id;
    hid_t type_id;
    hid_t space_id;
    hsize_t dims [] = {1};

    grp_id [0] = root_id;
    
    grp_name = bp_dirparser (dirstr, &level);
    
    for (i = 0; i < level - 1; i++)
    {
        grp_id [i + 1] = H5Gopen (grp_id [i], grp_name [i]);
        if (grp_id [i + 1] < 0)
        {
            grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
        }
    }
    
    dataset_id = H5Dopen (grp_id [level - 1], grp_name [level - 1]);
    
    type_id = H5Tcopy (H5T_NATIVE_INT);
    if (dataset_id == -1)
    {
        // BUG:
        // create a dataset if the dataset specified by path doesn't exist
        // actually this code never get executed because ADIOS always writes
        // all vars before writing attributes
        // in case the dataset does not exist, then there is no way to ensure
        // the dataspace is single element array
        // FIX:
        // delay attribute write until corresponding dataset is written
        fprintf(stderr, "Warning in writing h5 file: you hit a bug (%s: %d)\n", __FILE__, __LINE__);

        space_id = H5Screate_simple (1, dims, NULL);
        dataset_id = H5Dcreate (grp_id [level - 1], grp_name [level - 1]
                               ,type_id, space_id, H5P_DEFAULT
                               );
    }

    hw_string_attr_ds_internal (dataset_id, aname, aval);

    H5Dclose (dataset_id);
    for (i = 1; i < level; i++)
        H5Gclose (grp_id [i]);

    hw_free2D (grp_name, level);
}

/*
 * String-typed group attribute is converted to hdf5 dataset attribute 
 * in C convention or Fortran convention. 
 */
void hw_string_attr_gp_internal (hid_t parent_id, const char *name,const char *value)
{
    if(attr_str_gp_fortran == USE_C) {
        hw_string_attr_c(parent_id, name, value);   
    }
    else if(attr_str_gp_fortran = USE_FORTRAN) {
        hw_string_attr_f(parent_id, name, value);    
    }
}

/*
 * String-typed dataset attribute is converted to hdf5 dataset attribute 
 * in C convention or Fortran convention. 
 */
void hw_string_attr_ds_internal (hid_t parent_id, const char *name,const char *value)
{
    if(attr_str_ds_fortran == USE_C) {
        hw_string_attr_c(parent_id, name, value);   
    }
    else if(attr_str_ds_fortran = USE_FORTRAN) {
        hw_string_attr_f(parent_id, name, value);    
    }
}

/*
 * Write string-typed attribute in C convention
 */
void hw_string_attr_c (hid_t parent_id, const char *name,const char *value)
{
    hid_t dspace_id,dtype_id,attr_id;
    dspace_id = H5Screate(H5S_SCALAR);
    if(dspace_id>0)
    {
        dtype_id = H5Tcopy(H5T_C_S1_g);
        if(dtype_id>0)
        {
            H5Tset_size(dtype_id,strlen(value)+1);
            attr_id = H5Acreate(parent_id,name,dtype_id,dspace_id,H5P_DEFAULT);
            if(attr_id>0)
            {
                H5Awrite(attr_id,dtype_id,value);
                H5Aclose(attr_id);
            }
            H5Tclose(dtype_id);
        }
        H5Sclose(dspace_id);
    }

    return;
}

/*
 * Write string-typed attribute in Fortran convention
 */
void hw_string_attr_f (hid_t parent_id, const char *name,const char *value)
{
    hid_t dspace_id, dtype_id, attr_id;
    hsize_t adims[1] = {1};

    dspace_id = H5Screate_simple(1, adims, NULL);
    if(dspace_id > 0)
    {
        dtype_id = H5Tcopy(H5T_FORTRAN_S1); // Fortran string
        if(dtype_id > 0)
        {
            H5Tset_size(dtype_id, strlen(value));
            attr_id = H5Acreate(parent_id, name, dtype_id, dspace_id, H5P_DEFAULT);
            if(attr_id > 0)
            {
                H5Awrite(attr_id, dtype_id, value);
                H5Aclose(attr_id);
            }
            H5Tclose(dtype_id);
        }
        H5Sclose(dspace_id);
    }

    return;
}

/*
 * hw_attr_num_gp() writes numeric-typed group attribute 
 * hw_attr_num_gp() finds the group to which the attribute is attached and
 * writes the attribute. 
 */
void hw_attr_num_gp(hid_t root_id, char *dirstr, char *aname, void *avalue, enum vartype_t type)
{
    H5Eset_auto(NULL,NULL);
    char **grp_name;
    int level,i;
    hid_t grp_id[NUM_GP+1];

    grp_id[0] = root_id;
    grp_name = bp_dirparser(dirstr,&level);
    
    for(i=0;i<level;i++)
    {
        grp_id[i+1] = H5Gopen(grp_id[i],grp_name[i]);
        if(grp_id[i+1]<0)
        {
            grp_id[i+1] = H5Gcreate(grp_id[i],grp_name[i], 0);
        }
    }

    hw_scalar_attr(grp_id[level], aname, avalue,type);

    for(i=1;i<level;i++)
        H5Gclose(grp_id[i]);

    hw_free2D(grp_name,level);
}

/*
 * hw_attr_num_ds() writes numeric-typed dataset attribute 
 * hw_attr_num_ds() finds the dataset to which the attribute is attached and
 * writes the attribute. 
 */
void hw_attr_num_ds(hid_t root_id, char *dirstr, char *aname, void *avalue, enum vartype_t type)
{
    H5Eset_auto(NULL,NULL);
    char **grp_name;
    int level,i;
    hid_t grp_id[NUM_GP+1],dataset_id,type_id, space_id;
    hsize_t dims[]={1};

    grp_id[0] = root_id;
    
    grp_name = bp_dirparser(dirstr,&level);
    
    for(i=0;i<level-1;i++)
    {
        grp_id[i+1] = H5Gopen(grp_id[i],grp_name[i]);
        if(grp_id[i+1]<0)
        {
            grp_id[i+1] = H5Gcreate(grp_id[i],grp_name[i], 0);
        }
    }

    dataset_id = H5Dopen(grp_id[level-1],grp_name[level-1]);
    
    if(dataset_id==-1)
    {
        // BUG:
        // create a dataset if the dataset specified by path doesn't exist
        // actually this code never get executed because ADIOS always writes
        // all vars before writing attributes
        // in case the dataset does not exist, then there is no way to ensure
        // the dataspace is single element array
        // FIX:
        // delay attribute write until corresponding dataset is written
        fprintf(stderr, "Warning in writing h5 file: you hit a bug (%s: %d)\n", __FILE__, __LINE__);

        type_id = H5Tcopy(H5T_NATIVE_INT);
        space_id = H5Screate_simple(1,dims,NULL);
        dataset_id = H5Dcreate(grp_id[level-1],grp_name[level-1],type_id,space_id,H5P_DEFAULT);
    }

    hw_scalar_attr(dataset_id, aname, avalue,type);

    H5Dclose(dataset_id);
    for(i=1;i<level;i++)
        H5Gclose(grp_id[i]);

    hw_free2D(grp_name,level);
}

/*
 * Write an integer attribute as a scalar or a single-element-array
 */
void hw_scalar_attr( hid_t parent_id, const char *name, void *value,enum vartype_t type)
{
    if(scalar_dataspace == USE_SCALAR) {
        hw_scalar_attr_scalar (parent_id, name, value, type);
    }
    else if(scalar_dataspace == USE_SINGLE_ELE_ARRAY){
        hw_scalar_attr_array (parent_id, name, value, type);
    }
}

/*
 * Write an integer attribute to a h5 file as a scalar
 *
 * Edited from the AVS example file src/hdf5/examp/write_struct.c
 */
void hw_scalar_attr_scalar ( hid_t parent_id, const char *name, void *value, enum vartype_t type)
{
    hid_t h5_dspace_id, h5_attr_id, h5_type_id;
    int status;
    herr_t h5_status;
  
    h5_dspace_id = H5Screate( H5S_SCALAR );
    if( h5_dspace_id > 0 ) 
    {
        status = bp_getH5TypeId( type, &h5_type_id, value );
        if ( status == 0 && h5_type_id > 0 ) 
        {
            h5_attr_id = H5Acreate( parent_id, name, h5_type_id,h5_dspace_id, H5P_DEFAULT );
            if( h5_attr_id > 0 ) 
            {
                h5_status = H5Awrite( h5_attr_id, h5_type_id, value );
                if (h5_status < 0)
                    fprintf(stderr, "Failed to write an attribute (%s) to the h5 file!\n",name);
                H5Aclose( h5_attr_id );     /* close attribute */
            }
            H5Sclose( h5_dspace_id );       /* close dataspace */
        }
        H5Tclose( h5_type_id );
    }
    return;
}

/*
 * Write an integer attribute to a h5 file as a single-element-array
 */
void hw_scalar_attr_array ( hid_t parent_id, const char *name, void *value, enum vartype_t type)
{
    hid_t h5_dspace_id, h5_attr_id, h5_type_id;
    int status;
    herr_t h5_status;
    hsize_t dims[1] = {1};
  
    h5_dspace_id = H5Screate_simple(1, dims, NULL);
    if( h5_dspace_id > 0 ) 
    {
        status = bp_getH5TypeId( type, &h5_type_id, value );
        if ( status == 0 && h5_type_id > 0 ) 
        {
            h5_attr_id = H5Acreate( parent_id, name, h5_type_id, h5_dspace_id, H5P_DEFAULT );
            if( h5_attr_id > 0 ) 
            {
                h5_status = H5Awrite( h5_attr_id, h5_type_id, value );
                if (h5_status < 0)
                    fprintf(stderr, "Failed to write an attribute (%s) to the h5 file!\n", name);
                H5Aclose( h5_attr_id );     /* close attribute */
            }
            H5Sclose( h5_dspace_id );       /* close dataspace */
        }
        H5Tclose( h5_type_id );
    }
    return;
}

/*
 * Write an array as a dataset to a h5 file 
 */
void hw_dataset(hid_t parent_id, char* name, void* data,enum vartype_t type, int rank, hsize_t* dims)
{
    hid_t dataset_id, dataspace_id, cparms, type_id,filespace;    
    int i,rank_old;
    herr_t h5_status;
    hsize_t *offset;
    hsize_t *maxdims;
    maxdims = (hsize_t*)malloc(sizeof(hsize_t)*rank);
    offset = (hsize_t*)malloc(sizeof(hsize_t)*rank);
    for(i=0;i<rank;i++)
    {
        maxdims[i] = H5S_UNLIMITED;    
        offset[i] = 0;
    }
    cparms = H5Pcreate(H5P_DATASET_CREATE);
    h5_status = H5Pset_chunk(cparms,rank,dims);
    h5_status = bp_getH5TypeId(type, &type_id, data);
    if(h5_status == 0 && type_id>0)
    {
        dataset_id = H5Dopen(parent_id,name);
        if(dataset_id<0) 
        {
            dataspace_id = H5Screate_simple(rank, dims, NULL);
            if(dataspace_id>0 && h5_status==0) 
            {
                dataset_id = H5Dcreate(parent_id, name, type_id, dataspace_id,cparms);
                if(dataset_id<0)
                {
                    H5Sclose(dataspace_id);
                    return;
                }
                h5_status = H5Dextend (dataset_id, dims);
                filespace = H5Dget_space(dataset_id);
                h5_status = H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,dims,NULL);
                h5_status = H5Dwrite(dataset_id,type_id,dataspace_id,filespace,H5P_DEFAULT,data);
            }
        }
        else
        {
            filespace = H5Dget_space(dataset_id);
            rank_old = H5Sget_simple_extent_ndims(filespace);
            if(rank_old!=rank && filespace>0)
            {
                H5Sclose(filespace);
                return;
            }
            h5_status = H5Sget_simple_extent_dims(filespace,maxdims,NULL);
            
            offset[0] = maxdims[0];
            maxdims [0] += dims [0];
            h5_status = H5Dextend (dataset_id, maxdims);
            dataspace_id = H5Screate_simple(rank,dims,NULL);
            h5_status = H5Sget_simple_extent_dims(filespace,maxdims,NULL);
            filespace = H5Dget_space(dataset_id);

            if(verbose >= DEBUG_INFO) {
//                printf("parent_id=%d,dataset_id=%d, name=%s,filespace=%d\n",\
//                        parent_id,dataset_id, name,filespace);
            }

            h5_status = H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,dims,NULL);
            dataspace_id = H5Screate_simple(rank, dims, NULL);
            h5_status = H5Dwrite(dataset_id,type_id,dataspace_id,filespace,H5P_DEFAULT,data);
        }
        if(dataset_id>0)
            h5_status = H5Dclose(dataset_id);
        if(dataspace_id>0)
            h5_status = H5Sclose(dataspace_id);
        if(filespace>0)
            h5_status = H5Sclose(filespace);
        
    }
    H5Tclose(type_id);
    H5Pclose(cparms);
    free(maxdims);
    free(offset);
    return;
}

/*
 * Maps bp datatypes to h5 datatypes 
 */
int bp_getH5TypeId(enum vartype_t type, hid_t* h5_type_id, void * val)
{
    int size, status=0;
    switch (type)
    {
    case bp_float:
        *h5_type_id = H5Tcopy(H5T_NATIVE_FLOAT);
        break;
    case bp_double:
        *h5_type_id = H5Tcopy(H5T_NATIVE_DOUBLE);
        break;
    case bp_longdouble:
        *h5_type_id = H5Tcopy(H5T_NATIVE_LDOUBLE);
        break;
    case bp_uchar:
    case bp_ushort:
    case bp_uint:
    case bp_ulong:
    case bp_ulonglong:
#ifdef BP_USEUNSIGNED
        size = bp_getsize(type, val);
        if (size == 1)
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT8);
        else if (size == 2)
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT16);
        else if (size == 4)
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT32);
         else if (size == 8)
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT64);
        else
            status = -1;
        break;
#endif
    case bp_char:
    case bp_string:
    case bp_short:
    case bp_int:
    case bp_long:
    case bp_longlong:
        size = bp_getsize(type, val);
        if (size == 1)
            *h5_type_id = H5Tcopy(H5T_NATIVE_INT8);
        else if (size == 2)
            *h5_type_id = H5Tcopy(H5T_NATIVE_INT16);
        else if (size == 4)
            *h5_type_id = H5Tcopy(H5T_NATIVE_INT32);
        else if (size == 8)
            *h5_type_id = H5Tcopy(H5T_NATIVE_INT64);
        else
            status = -1;
            break;
    default:
        status = -1;
    }
    return status;
}

