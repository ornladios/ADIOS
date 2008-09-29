#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stdint.h"
#include <signal.h>
#include <sys/types.h>
#include "adios_types.h"
#include "adios_transport_hooks.h"
#include "adios_bp_v1.h"
#include "adios_internals.h"
#include "hw-utils.h"

#define NUM_GP 24
#define  MAX_RANK 10 * 5
#define STRLEN 1000 

struct var_dim
{
    uint16_t id;
    uint64_t rank;
};

/*
 * Global config variables to set the convenstion of conversion
 */

/*
 * array to dataset convention (transpose array dimension order if using Fortran)
 */
enum lang_convention array_dim_order_fortran = USE_FORTRAN;

/*
 * string-typed scalar/dataset write convention
 */
enum lang_convention var_str_fortran = USE_FORTRAN;

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
const char * value_to_string (enum ADIOS_DATATYPES type, void * data);
enum verbose_level verbose = NO_INFO;
void set_lang_convention(enum ADIOS_FLAG host_language_fortran)
{
    if(host_language_fortran == adios_flag_yes) {
        array_dim_order_fortran = USE_FORTRAN;
        var_str_fortran = USE_FORTRAN;
        attr_str_ds_fortran = USE_FORTRAN;
        attr_str_gp_fortran = USE_FORTRAN;
    }
    else {
        // otherwise assume host language is C
        array_dim_order_fortran = USE_C;
        var_str_fortran = USE_C;
        attr_str_ds_fortran = USE_C;
        attr_str_gp_fortran = USE_C;
    }
}

/*
 * Parse path of var/attribute
 *
 * (copied from binpack-utils.c; this function
 * is no longer included in libadios.a or 
 * libadios_int.a)
 */ 
char** bp_dirparser(char *str, int *nLevel)
{
  char **grp_name;
  char *pch;
  int idx = 0, len=0;
  char *tmpstr;
  tmpstr= (char *)malloc(1*(strlen(str)+1));
  strcpy(tmpstr,str);
  pch = strtok(tmpstr,"/");
  grp_name = (char **)malloc(NUM_GP*1);
  while(pch!=NULL && *pch!=' ')
  {

     len = strlen(pch);
     grp_name[idx]  = (char *)malloc((len+1)*1);
     grp_name[idx][0]='\0';
     strcat(grp_name[idx],pch);
     pch=strtok(NULL,"/");
     idx=idx+1;
  }
  *nLevel = idx;
  free(tmpstr);
  return grp_name;
}

/*
 * Initialization
 * initialize_bp2h5() sets config variables and allocate
 * internal read buffer.
 * It returns 0 if no error and -1 otherwise
 */
int initialize_bp2h5(enum lang_convention array_dim_order,
                     enum lang_convention var_str,
                     enum lang_convention attr_str_ds,
                     enum lang_convention attr_str_gp,
                     enum scalar_convention scalar_ds,
                     enum verbose_level verb
                     )
{
    if(!bp2h5_initialized) {
        array_dim_order_fortran = array_dim_order;
        var_str_fortran = var_str;
        attr_str_ds_fortran = attr_str_ds;
        attr_str_gp_fortran = attr_str_gp;
        scalar_dataspace = scalar_ds;
        verbose = verb;
        
        bp2h5_initialized = 1;
    }

    return 0;
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
    char * tmpstr;
    int size;
    int i;
    int rc;
    uint64_t element_size = 0;
    struct adios_bp_element_struct * element = NULL;

    // initialization
    if(!bp2h5_initialized && 
        (initialize_bp2h5(USE_FORTRAN, USE_FORTRAN, USE_FORTRAN, USE_FORTRAN, USE_SCALAR, NO_INFO))) {
        return -1;
    }

    // open bp file for read
    struct adios_bp_buffer_struct_v1 * b = 0;
    uint32_t version = 0;

    b = malloc (sizeof (struct adios_bp_buffer_struct_v1));
    adios_buffer_struct_init (b);

    rc = adios_posix_open_read_internal (fnamein, "", b);
    if (!rc)
    {
        fprintf (stderr, "Error in open bp file: cannot open %s file!\n", fnamein);

        return -1;
    }

    tmpstr = strdup (fnamein);
    size = strlen (fnamein);
    tmpstr [size - 2] = 'h';
    tmpstr [size - 1] = '5';
    if(!fnameout) {
        fnameout=tmpstr;
    }

    // create h5 file
    hid_t h5file_id;
    hid_t root_id;
    h5file_id = H5Fcreate (fnameout, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    root_id = H5Gopen (h5file_id, "/");

    // read and parse footer
    adios_posix_read_version (b);
    adios_parse_version (b, &version);

    struct adios_index_process_group_struct_v1 * pg_root = 0;
    struct adios_index_process_group_struct_v1 * pg = 0;
    struct adios_index_var_struct_v1 * vars_root = 0;
    struct adios_index_attribute_struct_v1 * attrs_root = 0;

    // read and parse index
    adios_posix_read_index_offsets (b);
    adios_parse_index_offsets_v1 (b);

    adios_posix_read_process_group_index (b);
    adios_parse_process_group_index_v1 (b, &pg_root);

    adios_posix_read_vars_index (b);
    adios_parse_vars_index_v1 (b, &vars_root);
    uint64_t dimval;
    int var_dims_count = 0;
    struct var_dim * var_dims = 0;
    while (vars_root)
    {
       if (vars_root->characteristics_count==1 && vars_root->characteristics[0].dims.count==0) {
          var_dims = realloc (var_dims,   (var_dims_count + 1)
                             * sizeof (struct var_dim)
                             );
          var_dims [var_dims_count].id = vars_root->id;
          switch (vars_root->type ) {
              case adios_integer:
                   dimval =(uint64_t) *(int *) vars_root->characteristics[0].value;
                   printf("%s %d %llu\n", vars_root->var_name,vars_root->id,dimval); 
                   break;
              case adios_long: 
                   dimval =(uint64_t) *(long *) vars_root->characteristics[0].value;
                   printf("%s %d %llu\n", vars_root->var_name,vars_root->id,dimval); 
                   break;
              default:
                   printf("ERROR: unsupported var adios type %s %d\n",vars_root->var_name,vars_root->type); 
                   break;
           }
           var_dims [var_dims_count].rank = dimval;
           var_dims_count++;
       }
       vars_root = vars_root->next;
    } 
    adios_posix_read_attributes_index (b);
    adios_parse_attributes_index_v1 (b, &attrs_root);
    while (attrs_root)
    {
       if  (  attrs_root->characteristics_count==1 
           && attrs_root->characteristics[0].dims.count==0
           && attrs_root->type != adios_string) {
          var_dims = realloc (var_dims,   (var_dims_count + 1)
                             * sizeof (struct var_dim)
                             );
          var_dims [var_dims_count].id = attrs_root->id;
          if (attrs_root->characteristics[0].var_id !=0) {
              printf("attribute:%s %d\n",attrs_root->attr_name, attrs_root->characteristics[0].var_id); 
              for (i = 0; i<var_dims_count;i++) {
                   if (var_dims[i].id == attrs_root->characteristics[0].var_id) {
                       dimval = var_dims[i].rank;
                       printf("attribute:%s %d\n",attrs_root->attr_name, dimval); 
                       hw_attr_num_ds (root_id, attrs_root->attr_path,attrs_root->attr_name, &dimval
                                       ,adios_long
                                      );
                   }
              }
          }
          else {
              switch (attrs_root->type ) {
                  case adios_integer:
                      dimval =(uint64_t) *(int *) attrs_root->characteristics[0].value;
                      //printf("%s %d %llu\n", attrs_root->attr_name,attrs_root->id,dimval); 
                      break;
                  case adios_long: 
		      dimval =(uint64_t) *(long *) attrs_root->characteristics[0].value;
		      //printf("%s %d %llu\n", attrs_root->attr_name,attrs_root->id,dimval); 
		      break;
		  default:
		      printf("ERROR: unsupported adios type: %s \n",attrs_root->attr_name);//,vars_root->type); 
		      //value_to_string (attrs_root->type,attrs_root->characteristics[0].value)); 
		      break ;
              }
          }
           var_dims [var_dims_count].rank = dimval;
           var_dims_count++;
       }
       attrs_root = attrs_root->next;
    }
/*
    for (i = 0; i< var_dims_count; i++)
       printf("total: %d %d %d\n",var_dims_count, var_dims[i].id, var_dims[i].rank); 
*/
    // xxx.bp --> xxx.h5
    // parse element from bp file and write to hdf5 file
    uint64_t element_num = 1;
    pg = pg_root;
    while (pg)
    {

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
        //print_process_group_header (element_num++, &pg_header);

        adios_parse_vars_header_v1 (b, &vars_header);

        set_lang_convention(pg_header.host_language_fortran);

        // process each var in current process group
        int i;
        for (i = 0; i < vars_header.count; i++)
        {
            var_payload.payload = 0;

            adios_parse_var_data_header_v1 (b, &var_header);
            //print_var_header (&var_header);

            //if (var_header.is_dim == adios_flag_yes)
            //{
                var_payload.payload = malloc (var_header.payload_size);
                adios_parse_var_data_payload_v1 (b, &var_header, &var_payload
                                                ,var_header.payload_size
                                                );
            //}
            //else
            //{
            //    adios_parse_var_data_payload_v1 (b, &var_header, NULL, 0);
            //}


            // write to h5 file
            // make sure the buffer is big enough or send in null
            if (var_header.dims) 
            { 
                // dataset
                    
                uint64_t element = 0;
                int ranks = 0;
                struct adios_dimension_struct_v1 * d = var_header.dims;
                int c = 0;
                uint64_t * dims;
                uint64_t * global_dims;
                uint64_t * offsets;
                int i = 0;

                while (d)
                {
                    ranks++;
                    d = d->next;
                }

                // get local dimensions
                dims = (uint64_t *) malloc (8 * ranks);
                memset (dims, 0, 8 * ranks);

                d = var_header.dims;
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
                        *dims_t = d->dimension.rank;
                    }

                    d = d->next;
                    dims_t++;
                }
                    
                // get global dimensions
                global_dims = (uint64_t *) malloc (8 * ranks);
                memset (global_dims, 0, 8 * ranks);

                d = var_header.dims;
                dims_t = global_dims;

                while (d)
                {
                    if (d->global_dimension.var_id != 0)
                    {
                        for (i = 0; i < var_dims_count; i++)
                        {
                            if (var_dims [i].id == d->global_dimension.var_id)
                            {
                                *dims_t = var_dims [i].rank;
                            }
                        }
                    }
                    else
                    {
                        *dims_t = d->global_dimension.rank;
                    }
                    d = d->next;
                    dims_t++;
                }

                // get offsets
                offsets = (uint64_t *) malloc (8 * ranks);
                memset (offsets, 0, 8 * ranks);

                d = var_header.dims;
                dims_t = offsets;
                while (d) {
                    if (d->local_offset.var_id != 0)
                    {
                        for (i = 0; i < var_dims_count; i++)
                        {
                            if (var_dims [i].id == d->local_offset.var_id)
                            {
                                *dims_t = var_dims [i].rank;
                            }
                        }
                    }
                    else
                    {
                        *dims_t = d->local_offset.rank;
                    }

                    d = d->next;
                    dims_t++;
                }
                // now ready to write dataset to h5 file
                hw_dset (root_id, var_header.path, var_header.name, var_payload.payload
                        ,var_header.type, ranks, dims, global_dims, offsets
                        );        
                free(dims);
                free(global_dims);
                free(offsets);
            }
            else {
                // scalar var
                hw_scalar (root_id, var_header.path, var_header.name
                          ,var_payload.payload
                          ,var_header.type
                          ,0
                          );     
            }
  
            if (var_payload.payload)
            {
                free (var_payload.payload);
            }
        }

        adios_parse_attributes_header_v1 (b, &attrs_header);
        // process each attribute in current process group
        for (i = 0; i < attrs_header.count; i++)
        {
            adios_parse_attribute_v1 (b, &attribute);
            //printf("%s_%s %s\n",attribute.path,attribute.name, attribute.value);
            // write to h5 file
            if(attribute.is_var == adios_flag_no) {
                switch(attribute.type) 
                {
                    case adios_string:
                        hw_attr_str_ds (root_id, attribute.path, attribute.name, attribute.value);
                        break;
                    case adios_byte:
                    case adios_real:
                    case adios_double:
                    case adios_long_double:
                    case adios_unsigned_byte:
                    case adios_unsigned_short:
                    case adios_unsigned_integer:
                    case adios_unsigned_long:
                    case adios_short:
                    case adios_integer:
                    case adios_long:
                         hw_attr_num_ds (root_id, attribute.path, attribute.name, attribute.value
                                       ,attribute.type
                                       );
                        break;
                    case adios_complex:
                    case adios_double_complex:
                        fprintf(stderr, "Error in mapping ADIOS Data Types to HDF5: complex not supported yet.\n");
                        break;
                    case adios_unknown:
                    default:
                        fprintf(stderr, "Error in mapping ADIOS Data Types to HDF5: unknown data type.\n");
                }
            }
            else {
                // var 
                //hw_attr_num_ds (root_id, attribute.path, attribute.name, attribute.value
                //                ,attribute.type
                //                );
                 
            }
        }

        var_dims_count = 0;
        
        if (var_dims) {
            free (var_dims);
        }

        pg = pg->next;
    }

    adios_posix_close_internal (b);

    H5Gclose (root_id);
    H5Fclose (h5file_id);

    return 0;
}

/*
 * Write dataset to h5 file
 */
void hw_dset (hid_t root_id, char * dirstr, char * name, void * data
             ,enum ADIOS_DATATYPES type, int rank
             ,uint64_t *dims
             ,uint64_t *global_dims
             ,uint64_t *offsets
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

    if (global_dims[0])
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

                global_h5dims [reverse_i] = global_dims[i];

                local_h5dims [reverse_i] = dims[i];

                start [reverse_i] = offsets[i];

                stride [reverse_i] = 1;

                count [reverse_i] = dims[i];
            }
        }
        else if(array_dim_order_fortran == USE_C) {
            for (int i = 0; i < rank; i++)
            {
                global_h5dims [i] = global_dims[i];

                local_h5dims [i] = dims[i];

                start [i] = offsets[i];

                stride [i] = 1;

                count [i] = dims[i];
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
                h5dims [rank-1-i] = dims[i];
            }
        }
        else if(array_dim_order_fortran == USE_C) {
            for (int i = 0; i < rank; i++)
            {
                h5dims [i] = dims[i];
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
               ,enum ADIOS_DATATYPES type, int append
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
        if (grp_id [i + 1] < 0) {
            grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
        }
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
                 ,enum ADIOS_DATATYPES type, int append
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
                 ,enum ADIOS_DATATYPES type, int ndims, hsize_t * dims
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
    
    for (i = 0; i < level; i++)
    {
        grp_id [i + 1] = H5Gopen (grp_id [i], grp_name [i]);
        if (grp_id [i + 1] < 0)
        {
            grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
        }
    }
    if (grp_id[level] > 0) {
        hw_string_attr_gp_internal (grp_id [level], aname, aval);

        for (i = 1; i < level; i++)
            H5Gclose (grp_id [i]);

        hw_free2D (grp_name, level);
        return;
    }
    for (i = 0; i < level - 1; i++) {
        grp_id [i + 1] = H5Gopen (grp_id [i], grp_name [i]);
	if (grp_id [i + 1] < 0) {
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
        //fprintf(stderr, "Warning in writing h5 file: you hit a bug (%s: %d)\n", __FILE__, __LINE__);

        space_id = H5Screate_simple (1, dims, NULL);
        dataset_id = H5Dcreate (grp_id [level - 1], grp_name [level - 1]
                               ,type_id, space_id, H5P_DEFAULT
                               );
        if(dataset_id==-1)
           hw_attr_str_gp (root_id, dirstr, aname, aval);
    }
   
    if(dataset_id>0) hw_string_attr_ds_internal (dataset_id, aname, aval);

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
void hw_attr_num_gp(hid_t root_id, char *dirstr, char *aname, void *avalue, enum ADIOS_DATATYPES type)
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
void hw_attr_num_ds(hid_t root_id, char *dirstr, char *aname, void *avalue, enum ADIOS_DATATYPES type)
{
    H5Eset_auto(NULL,NULL);
    char **grp_name;
    int level,i;
    hid_t grp_id[NUM_GP+1],dataset_id,type_id, space_id;
    hsize_t dims[]={1};

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

    if (grp_id[level] > 0) {    
        hw_scalar_attr(grp_id[level], aname, avalue,type);
        for(i=1;i<level;i++)
            H5Gclose(grp_id[i]);
        hw_free2D(grp_name,level);
        return; 
    }
    
    for (i = 0; i < level - 1; i++) {
        grp_id [i + 1] = H5Gopen (grp_id [i], grp_name [i]);
	if (grp_id [i + 1] < 0) {
            grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
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
        if(dataset_id==-1)
        {
           printf("hit group att!%d\n",type);
           if(type==adios_string)
              hw_attr_str_gp (root_id, dirstr, aname, avalue);
           else
              hw_attr_num_gp (root_id, dirstr, aname, avalue,type);
        }
    }

    if(dataset_id>0)
    {
       hw_scalar_attr(dataset_id, aname, avalue,type);
      //printf("%s: dataset_id:%d, %s\n",grp_name[0],dataset_id,aname);
       H5Dclose(dataset_id);
    }


    for(i=1;i<level;i++)
        H5Gclose(grp_id[i]);

    hw_free2D(grp_name,level);
}

/*
 * Write an integer attribute as a scalar or a single-element-array
 */
void hw_scalar_attr( hid_t parent_id, const char *name, void *value,enum ADIOS_DATATYPES type)
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
void hw_scalar_attr_scalar ( hid_t parent_id, const char *name, void *value, enum ADIOS_DATATYPES type)
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
void hw_scalar_attr_array ( hid_t parent_id, const char *name, void *value, enum ADIOS_DATATYPES type)
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
void hw_dataset(hid_t parent_id, char* name, void* data,enum ADIOS_DATATYPES type, int rank, hsize_t* dims)
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
 *
 * The Mapping is according to HDF5 Reference Manual
 * (http://hdf.ncsa.uiuc.edu/HDF5/doc1.6/Datatypes.html)
 */
int bp_getH5TypeId(enum ADIOS_DATATYPES type, hid_t* h5_type_id, void * val)
{
    int size, status=0;

    switch (type)
    {
        case adios_byte:
            *h5_type_id = H5Tcopy(H5T_NATIVE_CHAR);
            break;
        case adios_string:
            if(var_str_fortran = USE_FORTRAN) {
                *h5_type_id = H5Tcopy(H5T_FORTRAN_S1);
            }
            else { // otherwise assume C
                *h5_type_id = H5Tcopy(H5T_C_S1_g);
            }
            break;
        case adios_real:
            *h5_type_id = H5Tcopy(H5T_NATIVE_FLOAT);
            break;
        case adios_double:
            *h5_type_id = H5Tcopy(H5T_NATIVE_DOUBLE);
            break;
        case adios_long_double:
            *h5_type_id = H5Tcopy(H5T_NATIVE_LDOUBLE);
            break;
        case adios_unsigned_byte:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT8);
            break;
        case adios_unsigned_short:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT16);
            break;
        case adios_unsigned_integer:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT32);
            break;
        case adios_unsigned_long:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT64);
            break;
        case adios_short:
            *h5_type_id = H5Tcopy(H5T_NATIVE_INT16);
            break;
        case adios_integer:
            *h5_type_id = H5Tcopy(H5T_NATIVE_INT32);
            break;
        case adios_long:
            *h5_type_id = H5Tcopy(H5T_NATIVE_INT64);
            break;
        case adios_complex:
        case adios_double_complex:
            fprintf(stderr, "Error in mapping ADIOS Data Types to HDF5: complex not supported yet.\n");
            status = -1;
            break;
        case adios_unknown:
        default:
            fprintf(stderr, "Error in mapping ADIOS Data Types to HDF5: unknown data type.\n");
            status = -1;
    }
    return status;
}

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
            sprintf (s, "%lld", *(((int64_t *) data)));
            break;

        case adios_unsigned_long:
            sprintf (s, "%llu", *(((uint64_t *) data)));
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
