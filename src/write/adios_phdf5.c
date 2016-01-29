/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mpi.h"
#include "public/adios_types.h"
#include "core/adios_bp_v1.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_internals.h"
#include "core/util.h"
#ifdef PHDF5 
#include "hdf5.h"
#endif

static const int debug=0;

#define NUM_GP 24
void adios_phdf5_end_iteration (struct adios_method_struct * method)
{
}

void adios_phdf5_start_calculation (struct adios_method_struct * method)
{
}

void adios_phdf5_stop_calculation (struct adios_method_struct * method)
{
}
void adios_phdf5_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
}


#ifndef PHDF5
void adios_phdf5_init(const PairStruct * parameters
                     ,struct adios_method_struct * method
                     ){}
void adios_phdf5_finalize (int mype, struct adios_method_struct * method){}
enum BUFFERING_STRATEGY adios_phdf5_should_buffer (
                     struct adios_file_struct * fd
                    ,struct adios_method_struct * method
                    ){ return no_buffering; }
int adios_phdf5_open(struct adios_file_struct *fd
                    ,struct adios_method_struct * method
                    ,MPI_Comm comm
                    ){ return -1; }
void adios_phdf5_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     ){}
void adios_phdf5_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,const void * data
                       ,struct adios_method_struct * method
                       ){}
void adios_phdf5_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffersize
                    ,struct adios_method_struct * method
                    ){}
#else

///////////////////////////
// Function Declarations
///////////////////////////
int adios_getsize(enum ADIOS_DATATYPES type, void * val);

// adios_flag determine whether it is dataset or group

void hw_gopen (hid_t root_id, char * path, hid_t * grp_id, int * level,enum ADIOS_FLAG *grpflag);
void hw_gclose ( hid_t * grp_ids, int level, enum ADIOS_FLAG grpflag);

int getH5TypeId ( enum ADIOS_DATATYPES type, hid_t* h5_type_id, enum ADIOS_FLAG);
hsize_t parse_dimension(struct adios_var_struct *pvar_root
                       ,struct adios_attribute_struct *patt_root 
                       ,struct adios_dimension_item_struct * dim); 
int hw_var (hid_t root_id
           ,struct adios_var_struct *pvar_root
           ,struct adios_attribute_struct *patt
           ,struct adios_var_struct *pvar
           ,enum ADIOS_FLAG fortran_flag 
           ,int myrank, int nproc);
int hr_var (hid_t root_id
           ,struct adios_var_struct *pvar_root
           ,struct adios_attribute_struct *patt_root
           ,struct adios_var_struct *pvar
           ,enum ADIOS_FLAG fortran_flag
           ,int myrank, int nproc);
           
int hw_attribute ( hid_t root_id
                  ,struct adios_var_struct *pvar_root
                  ,struct adios_attribute_struct *pvar
                  ,enum ADIOS_FLAG fortran_flag 
                  ,int myrank, int nproc);

struct adios_phdf5_data_struct
{
  hid_t fh;
  hid_t root_id;
  hid_t dxpl_id;
  MPI_Comm group_comm;
  int rank;
  int size;
};
int adios_phdf5_initialized = 0;
void adios_phdf5_init(const PairStruct * parameters
                     ,struct adios_method_struct * method
                     )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct *)
                                                    method->method_data;
    if (!adios_phdf5_initialized)
    {
        adios_phdf5_initialized = 1;
    }
    method->method_data = malloc (sizeof (struct adios_phdf5_data_struct));
    md = (struct adios_phdf5_data_struct *) method->method_data;
    md->fh = 0;
    md->root_id = 0;
    md->rank = -1;
    md->size = 0;
    md->group_comm = MPI_COMM_NULL;
}
enum BUFFERING_STRATEGY adios_phdf5_should_buffer (struct adios_file_struct * fd
                                                  ,struct adios_method_struct * method
                                                  )
{
    return no_buffering;
}
 
int adios_phdf5_open(struct adios_file_struct *fd
                    ,struct adios_method_struct * method, MPI_Comm comm
                    )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct *)
                                                    method->method_data;

    char * name;
    MPI_Info info = MPI_INFO_NULL;
    hid_t fapl_id;
    fapl_id = H5P_DEFAULT;
    // no shared buffer 
    //fd->shared_buffer = adios_flag_no;

    // create a new file. If file exists its contents will be overwritten. //
    md->group_comm = comm;
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }
    else 
       md->group_comm=MPI_COMM_SELF;

    fd->group->process_id = md->rank;
    name = malloc (strlen (method->base_path) + strlen (fd->name) + 1);
    sprintf(name, "%s%s", method->base_path, fd->name);

    H5Eset_auto ( NULL, NULL);

    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    //fprintf (stderr, "************\n%s: comm=%x\n**************\n", __func__, md->group_comm);
    H5Pset_fapl_mpio(fapl_id,md->group_comm,info);

    switch (fd->mode) {
        case adios_mode_read:
        {
            md->fh = H5Fopen (name, H5F_ACC_RDONLY, fapl_id);
            if (md->fh <= 0)
            {
                fprintf (stderr, "ADIOS PHDF5: file not found: %s\n", fd->name);
                free (name);
                return adios_flag_no;
            } 
            break;
        }
        case adios_mode_write:
        case adios_mode_append:
        case adios_mode_update:
            md->fh = H5Fcreate (name, H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
            if (md->fh < 0)
            {
                md->fh = H5Fopen (name, H5F_ACC_RDWR, fapl_id);
                //md->fh = H5Fcreate (name, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
            }
            if (md->fh < 0)
            {
                fprintf (stderr, "ADIOS PHDF5: file not create/append: %s\n", fd->name);
                free (name);
                return adios_flag_no;
            } 
            break;
    }

    md->root_id = H5Gopen(md->fh,"/");
    if(md->root_id < 0)
        md->root_id = H5Gcreate(md->fh,"/",0);
    H5Pclose(fapl_id);
    free (name); 
    return 1;
}

void adios_phdf5_write (struct adios_file_struct * fd
                       ,struct adios_var_struct * v
                       ,const void * data
                       ,struct adios_method_struct * method
                       )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct *)
                                                       method->method_data;
    if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        if (debug && md->rank==0) {
            fprintf(stderr, "-------------------------\n");
            fprintf(stderr, "write var: %s start!\n", v->name);
        }
        hw_var (md->root_id, fd->group->vars, fd->group->attributes
                 ,v, fd->group->adios_host_language_fortran,md->rank,md->size);
        MPI_Barrier(md->group_comm);
    }
    else
    {
       //fprintf(stderr, "entering unknown phdf5 mode %d!\n", fd->mode);
    }
    if (debug && md->rank==0) {
        fprintf(stderr, "write var: %s end!\n", v->name);
        //fprintf(stderr, "-------------------------\n");
    }
}


void adios_phdf5_read (struct adios_file_struct * fd
                      ,struct adios_var_struct * v
                      ,void * buffer
                      ,uint64_t buffersize
                      ,struct adios_method_struct * method
                      )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct *)
                                                       method->method_data;
    if(fd->mode == adios_mode_read)
    {
       v->adata = buffer;
       v->data_size = buffersize;
       
       if (md->rank==0) {
           fprintf(stderr, "-------------------------\n");
           fprintf(stderr, "read var: %s! start\n", v->name);
       }
       hr_var (md->root_id
              ,fd->group->vars
              ,fd->group->attributes
              ,v
              ,fd->group->adios_host_language_fortran
              ,md->rank
              ,md->size);
       v->adata = 0;
       if (md->rank==0) {
           fprintf(stderr, "read var: %s! end\n", v->name);
           //fprintf(stderr, "-------------------------\n");
       }
       
    }
}

/*static void adios_phdf5_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
// This function is not useful for phdf5 since adios_read/write do real read/write 
}*/


void adios_phdf5_buffer_overflow (struct adios_file_struct * fd, 
                                  struct adios_method_struct * method)
{
    // this never happens without shared buffering
}

void adios_phdf5_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct*)
                                                  method->method_data;
      
    struct adios_attribute_struct * a = fd->group->attributes;

    if (fd->mode == adios_mode_read) {
        if (debug && md->rank==0) {
           fprintf(stderr, "-------------------------\n");
           fprintf(stderr, "reading done, phdf5 file is closed;\n");
           fprintf(stderr, "-------------------------\n");
        }
    }
    else if (fd->mode == adios_mode_write || fd->mode == adios_mode_append)
    {
        //fprintf(stderr, "entering phdf5 write attribute mode!\n");
        while(a)
        {
            //fprintf(stderr, "Write attribute [%s]/[%s]\n", a->path, a->name);
            if (strcmp(a->path, "/__adios__")) {

                hw_attribute ( md->root_id, fd->group->vars, a
                        ,fd->group->adios_host_language_fortran
                        ,md->rank
                        ,md->size);
            }
            a = a->next;
        }
        if (debug && md->rank==0) {
           fprintf(stderr, "-------------------------\n");
           fprintf(stderr, "writing done, phdf5 file is closed;\n");
           fprintf(stderr, "-------------------------\n");
        }
    }

    if (md && md->fh && md->root_id) {
        H5Gclose (md->root_id);
    }
    H5Fclose (md->fh);
    md->group_comm = MPI_COMM_NULL;
    md->fh = 0;
    md->rank = -1;
    md->size = 0;
}

void adios_phdf5_finalize (int mype, struct adios_method_struct * method)
{
    // nothing to do here
    if (adios_phdf5_initialized)
        adios_phdf5_initialized = 0;
}

int hw_attribute ( hid_t root_id
                ,struct adios_var_struct *pvar_root
                ,struct adios_attribute_struct *patt
                ,enum ADIOS_FLAG fortran_flag
                ,int myrank
                ,int nproc) {

    H5Eset_auto ( NULL, NULL);
    int i, rank = 0, level;
    hid_t   h5_plist_id, h5_type_id, h5_dataspace_id;
    hid_t h5_attribute_id, grp_ids[NUM_GP];
    hsize_t * h5_localdims;
    struct adios_dimension_struct * dims;
    struct adios_var_struct * var_linked;
    enum ADIOS_FLAG flag = adios_flag_unknown;
     
    h5_plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_plist_id=H5P_DEFAULT;  
 
    H5Pset_dxpl_mpio(h5_plist_id, H5FD_MPIO_COLLECTIVE); 
    char *path = patt->path;
    /*
    if (!strcmp(patt->path, "/__adios__")) 
        path = strdup ("ADIOSINFO");
    fprintf (stderr, "%s: write attribute %s under path %s\n", __func__, patt->name, path);
    */
    hw_gopen (root_id, path, grp_ids, &level, &flag);
    int err_code = 0;
    //printf("patt->type=%d patt->name : %s\n", patt->type, patt->name);
    if (patt->type == -1) {
        var_linked = patt->var;
        if (!var_linked || (var_linked && !var_linked->data)) {
                fprintf (stderr, "PHDF5 ERROR: invalid data in var_linked"
                        " (in attribute write): %s(%d)\n"
                        ,var_linked->name, var_linked->id);
                err_code = -2;  
                H5Pclose(h5_plist_id);
                hw_gclose (grp_ids,level, flag);
                return err_code;
        }
        else
            dims = var_linked->dimensions;
        getH5TypeId (var_linked->type, &h5_type_id, fortran_flag);
        // Scalar variable as attribute
        if (!dims) {
            h5_dataspace_id = H5Screate ( H5S_SCALAR);
            if (h5_dataspace_id > 0) {
                h5_attribute_id = H5Aopen_name ( grp_ids[level], patt->name);
                if (h5_attribute_id < 0) {
                    h5_attribute_id = H5Acreate ( grp_ids[level], patt->name
                                          ,h5_type_id,h5_dataspace_id,0);
                }
                if (h5_attribute_id > 0) {
                        if (myrank == 0)        
                            H5Awrite ( h5_attribute_id, h5_type_id, var_linked->data);
                }
                H5Aclose (h5_attribute_id);
                H5Sclose (h5_dataspace_id);
            }
            else {
                fprintf (stderr, "PHDF5 ERROR in h5_dataspace_id in hw_attribute\n");  
                err_code = -2;
            }
         }
         else {
             while (dims) {
                ++rank;
                dims = dims->next;
             }
            
             h5_localdims = (hsize_t *) malloc (rank * sizeof(hsize_t));
             dims = var_linked->dimensions;
             for ( i = 0; i < rank; i++) {
                 if ( dims->dimension.var) { 
                     h5_localdims [i] = *(int *)dims->dimension.var->data;
                 }
                 else if ( dims->dimension.attr) { 
                     if ( dims->dimension.attr->var)
                         h5_localdims [i] = *(int *)dims->dimension.attr->var->data;
                     else 
                         h5_localdims [i] = *(int *)dims->dimension.attr->value;
                 } else {
                     h5_localdims [i] = dims->dimension.rank;
                 }
             }
             h5_dataspace_id = H5Screate_simple(rank,h5_localdims, NULL);
             h5_attribute_id = H5Aopen_name ( grp_ids[level], patt->name);
             if (h5_attribute_id < 0) {
                 h5_attribute_id = H5Acreate ( grp_ids[level], patt->name
                                          ,h5_type_id,h5_dataspace_id,0);
                 if (h5_attribute_id < 0) {
                     fprintf (stderr, "PHDF5 ERROR: getting negative attribute_id "
                                      "in hw_attribute: %s\n", patt->name);
                     err_code = -2; 
                  }
             }
             if (h5_attribute_id > 0) {
                 if (myrank == 0 && var_linked->data)
                     H5Awrite ( h5_attribute_id, h5_type_id, var_linked->data);
                 H5Aclose ( h5_attribute_id);
             }
             H5Sclose ( h5_dataspace_id);
             free (h5_localdims);
        }
    }
    if (patt->type > 0)
    {
        getH5TypeId (patt->type, &h5_type_id, fortran_flag);
        if (h5_type_id > 0) {
            switch (patt->type) {
            case adios_string:
                h5_dataspace_id = H5Screate ( H5S_SCALAR);
                H5Tset_size ( h5_type_id, (strlen((char *)patt->value)+1));
                h5_attribute_id = H5Aopen_name ( grp_ids[level], patt->name);
                if (h5_attribute_id < 0) { 
                    h5_attribute_id = H5Acreate ( grp_ids[level], patt->name 
                                        ,h5_type_id, h5_dataspace_id, 0);
                    if (h5_attribute_id > 0) { 
                       if (myrank == 0)        
                          H5Awrite(h5_attribute_id, h5_type_id, patt->value);
                     }
                }
                H5Aclose(h5_attribute_id);
                H5Sclose (h5_dataspace_id); 
                break; 
            case adios_integer: 
            case adios_long: 
                break; 
            default: 
                break;
           }
        }
    }
    H5Pclose (h5_plist_id); 
    hw_gclose (grp_ids,level, flag);
    return err_code;
} 
int hr_var (hid_t root_id
           ,struct adios_var_struct *pvar_root
           ,struct adios_attribute_struct *patt_root
           ,struct adios_var_struct *pvar
           ,enum ADIOS_FLAG fortran_flag
           ,int myrank
           ,int nproc) {

    H5Eset_auto ( NULL, NULL);
    int i, rank = 0, level, err_code = -2;
    hid_t   h5_plist_id, h5_type_id, h5_dataspace_id, h5_dataset_id, h5_memspace_id, grp_ids[NUM_GP];
    struct adios_dimension_struct * dims = pvar->dimensions;
    enum ADIOS_FLAG flag_yes = adios_flag_yes;

    h5_plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_plist_id = H5P_DEFAULT;
    H5Pset_dxpl_mpio(h5_plist_id, H5FD_MPIO_INDEPENDENT);
    H5Pclose (h5_plist_id);
    getH5TypeId (pvar->type, &h5_type_id, fortran_flag);
    if (h5_type_id <=0 )
    {
        fprintf (stderr, "ERROR in getH5TypeId in hr_var!\n");
        return -2;
    }

    //name = (char *) malloc (strlen(pvar->path)+strlen(pvar->name)-1);
    //sprintf (name, "%s", pvar->path);
    
    //if (pvar->path[strlen(pvar->path)-1]=='/')
    //    sprintf (name, "%s%s", name, pvar->name);
    //else
    //    sprintf (name, "%s/%s", name, pvar->name);

    if (pvar->path)
        hw_gopen (root_id, pvar->path, grp_ids, &level, &flag_yes);
    //printf("root_id=%d, grd_ids[%d]=%d\n", root_id, level, grp_ids[level]); 
    // variable is scalar only need to read by every processor 

    if (!dims) {
       
        h5_dataspace_id = H5Screate(H5S_SCALAR);
        h5_dataset_id = H5Dopen (grp_ids[level], pvar->name);
        if ( h5_dataset_id > 0) {
            H5Dread (h5_dataset_id, h5_type_id, H5S_ALL
                          ,H5S_ALL, h5_plist_id, pvar->adata
                          );
            H5Dclose (h5_dataset_id);
            err_code = 0;
        }
        else
            fprintf (stderr, "PHDF5 ERROR: can not open dataset: %s in hr_var\n", pvar->name);
        H5Sclose (h5_dataspace_id);
        H5Tclose (h5_type_id);
        hw_gclose (grp_ids,level, flag_yes);
        return err_code; 
    }
// variable is dataset
    while (dims) {
        ++rank;
        dims = dims->next;
    }
    dims = pvar->dimensions;

    if ( dims->global_dimension.rank || 
         dims->global_dimension.var  || 
         dims->global_dimension.attr) 
    {

        hsize_t * h5_globaldims, * h5_localdims, * h5_offsets, * h5_strides; 
        hsize_t h5_gbstrides[2],h5_gbglobaldims[2], h5_gblocaldims[2], h5_gboffsets[2], *h5_gbdims;
        char name[256];
        if(debug && myrank==0)printf("\tenter global writing!\n");
        h5_gbdims = (hsize_t *)malloc(rank * 3 * sizeof(hsize_t));
        h5_strides = (hsize_t *) malloc (rank * sizeof(hsize_t));

        h5_gbstrides [0] = 1;
        h5_gbstrides [1] = 1;
        h5_gbglobaldims [0] = nproc;
        h5_gbglobaldims [1] = rank * 3;
        h5_gboffsets [0] = myrank;
        h5_gboffsets [1] = 0;
        h5_gblocaldims [0] = 1; 
        h5_gblocaldims [1] = rank * 3;

        // get the global/local/offset arrays
        for (i = 0; i < rank; i++) {
            h5_strides [i] = 1;
            //dims = dims->next;
        } // end of for

        // save global-bounds information into the hdf5 files
        h5_dataspace_id = H5Screate_simple (2, h5_gbglobaldims, NULL);
        h5_memspace_id = H5Screate_simple (2, h5_gblocaldims, NULL);
        H5Sselect_hyperslab (h5_dataspace_id, H5S_SELECT_SET
                            ,h5_gboffsets, h5_gbstrides, h5_gblocaldims,0);
        sprintf(name, "_%s_gbdims", pvar->name);
         
        h5_dataset_id  = H5Dopen ( grp_ids[level], name);
        if (h5_dataset_id > 0) {
            H5Dread (h5_dataset_id, H5T_STD_I64LE, h5_memspace_id
                    ,h5_dataspace_id, h5_plist_id, h5_gbdims);
            H5Dclose (h5_dataset_id);
        }
        h5_globaldims= h5_gbdims;
        h5_localdims = h5_gbdims+rank;
        h5_offsets = h5_gbdims+rank*2;

        for (i=0;i<rank;i++) {
            if(myrank==0) 
                printf("\tDIMS var:%s dim[%d]:  %llu %llu %llu\n",pvar->name
                         ,i, h5_globaldims[i], h5_localdims[i], h5_offsets[i]);
        }
        H5Sclose(h5_dataspace_id); 
        H5Sclose(h5_memspace_id);
          
        //free(h5_gbdims);
        //free (h5_strides); 
        //return 0;
        
        h5_dataspace_id = H5Screate_simple (rank, h5_globaldims, NULL);
        err_code = -2;
        if (h5_dataspace_id > 0) {
            H5Sselect_hyperslab (h5_dataspace_id, H5S_SELECT_SET
                                ,h5_offsets, h5_strides, h5_localdims, 0
                                );
            h5_memspace_id = H5Screate_simple (rank, h5_localdims, NULL);
            if (h5_memspace_id > 0) {
                h5_dataset_id  = H5Dopen ( grp_ids[level], pvar->name);
                if (h5_dataset_id > 0) {
                    H5Dread (h5_dataset_id, h5_type_id, h5_memspace_id
                            ,h5_dataspace_id, h5_plist_id, pvar->adata
                            );
                    H5Dclose (h5_dataset_id);
                    err_code = 0;
                }
                else
                    fprintf (stderr, "PHDF5 ERROR: dataset %s does not existed!\n"\
                        ,pvar->name);
                
                H5Sclose (h5_memspace_id); 
            }
            else
                fprintf (stderr, "PHDF5 ERROR: out of memory, cannot create local space in hr_var: %s\n"\
                        ,pvar->name);
            H5Sclose (h5_dataspace_id);
        }
        else
            fprintf (stderr, "PHDF5 ERROR: out of memory, cannot create global space in hr_var: %s\n"\
                    ,pvar->name);
         
        free(h5_gbdims);
        free (h5_strides); 
    }  // end of writing dataset with global-bounds

    else {
        hsize_t *h5_localdims;
        h5_localdims = (hsize_t *) malloc (rank * sizeof(hsize_t));
        for (i = 0; i < rank; i++) {
            h5_localdims [i] = parse_dimension(pvar_root, patt_root, &dims->dimension);
            //printf("\t var name: %s, dims[%d] = %d\n",pvar->name, i, h5_localdims [i]);
            dims=dims->next;
        }
        h5_dataspace_id = H5Screate_simple (rank, h5_localdims, NULL);
        if ( h5_dataspace_id > 0) {
            h5_dataset_id  = H5Dopen ( grp_ids[level], pvar->name);
            if ( h5_dataset_id > 0) {
                H5Dread (h5_dataset_id, h5_type_id, H5S_ALL
                        ,H5S_ALL, h5_plist_id, pvar->adata
                        );
                H5Dclose (h5_dataset_id);
                err_code = 0;
            }
            else
                fprintf ( stderr, "PHDF5 ERROR:  cannot create dataset id for var: %s\n", pvar->name);
            H5Sclose (h5_dataspace_id);
        }
        else
            fprintf (stderr, "PHDF5 ERROR: cannot create dataset space %s for var!\n", pvar->name);
        free (h5_localdims);
    }
    hw_gclose(grp_ids, level, adios_flag_yes);
    H5Tclose (h5_type_id);
    H5Pclose (h5_plist_id);
    return err_code;
}

int hw_var (hid_t root_id
           ,struct adios_var_struct *pvar_root
           ,struct adios_attribute_struct *patt_root
           ,struct adios_var_struct *pvar
           ,enum ADIOS_FLAG fortran_flag
           ,int myrank
           ,int nproc) {

    H5Eset_auto ( NULL, NULL);
    int i, rank = 0, level;
    hid_t   h5_plist_id, h5_type_id, h5_dataspace_id, h5_dataset_id, h5_memspace_id, grp_ids[NUM_GP];
    struct adios_dimension_struct * dims = pvar->dimensions;
    enum ADIOS_FLAG flag_yes = adios_flag_yes;
     
    h5_plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_plist_id=H5P_DEFAULT;   
    //fprintf (stderr, "H5Pset_dxpl_mpio not found on ewok so commented out\n");
    H5Pset_dxpl_mpio(h5_plist_id, H5FD_MPIO_INDEPENDENT); 

    getH5TypeId (pvar->type, &h5_type_id, fortran_flag);
    if ( h5_type_id <= 0) {
       fprintf (stderr, "PHDF5 ERROR in getH5TypeId in hw_var\n");  
       return -2;
    }
/*//////////////////////////////////////////////////////////////////////////
// This is for the HDF5/1.8.1
    name = ( char *) malloc ( strlen (pvar->path) + strlen (pvar->name) -1);
    if ( pvar->path[strlen(pvar->path)-1]=='/')
        sprintf ( name, "%s%s", pvar->path, pvar->name);
    else
        sprintf ( name, "%s/%s", pvar->path, pvar->name);
*/////////////////////////////////////////////////////////////////////////////
    //for(i=0;i<nproc;i++)
    {
    if (pvar->path)
        hw_gopen (root_id, pvar->path, grp_ids, &level, &flag_yes);
    }
//    printf("root_id=%d, grp_id=%d\n", root_id, grp_ids[level]);

    if (!dims) {
        h5_dataspace_id = H5Screate(H5S_SCALAR);
        h5_dataset_id = H5Dopen (grp_ids[level], pvar->name);
        if (h5_dataset_id <= 0) {
            h5_dataset_id = H5Dcreate (grp_ids[level], pvar->name, h5_type_id
                                      ,h5_dataspace_id, H5P_DEFAULT
                                      );
            if (h5_dataset_id <= 0) 
                fprintf (stderr, "PHDF5 ERROR: can not create scalar %s in hw_var!\n", pvar->name); 
        }
        if (h5_dataset_id>0 ) {
        	if ( myrank==1)
        		H5Dwrite (h5_dataset_id, h5_type_id, H5S_ALL
        				  ,H5S_ALL, h5_plist_id, pvar->data
        		         );
            //printf("groupid=%d level=%d datasetid=%d\n",grp_ids[level],level,h5_dataset_id);
            //printf("write dataset: name=%s/%s myrank=%d\n"
            //         , pvar->path,pvar->name,myrank);
            H5Dclose (h5_dataset_id); 
        }
        H5Sclose (h5_dataspace_id); 
        H5Tclose (h5_type_id);
        H5Pclose (h5_plist_id);
        hw_gclose (grp_ids,level, flag_yes);
        return 0;
    }// end of scalar read

    while (dims) {
        ++rank;
        dims = dims->next;
    }
    dims = pvar->dimensions;
    hsize_t * h5_globaldims, * h5_localdims, * h5_offsets, * h5_strides, * h5_gbdims; 
    if ( dims->dimension.rank || 
         dims->dimension.var  ||
         dims->dimension.attr ) 
    {

        hsize_t h5_gbstrides[2], h5_gbglobaldims[2], h5_gblocaldims[2], h5_gboffsets[2];
        char name[256];
        if(debug&&myrank==0)printf("\tenter global reading!\n");

        h5_gbstrides [0] = 1;
        h5_gbstrides [1] = 1;
        h5_gbglobaldims [0] = nproc;
        h5_gbglobaldims [1] = rank * 3;
        h5_gboffsets [0] = myrank;
        h5_gboffsets [1] = 0;
        h5_gblocaldims [0] = 1; 
        h5_gblocaldims [1] = rank * 3;

        h5_gbdims = (hsize_t *) malloc (rank * 3 * sizeof(hsize_t));
        h5_strides = (hsize_t *) malloc (rank * sizeof(hsize_t));

        h5_globaldims = h5_gbdims;
        h5_localdims = h5_gbdims + rank;
        h5_offsets = h5_gbdims + 2 * rank;
 
        // get the global/local/offset arrays       
        for ( i = 0; i < rank; i++) {
             h5_strides [i] = 1;
             h5_globaldims[i] = parse_dimension (pvar_root, patt_root, &dims->global_dimension);
             h5_localdims[i] = parse_dimension (pvar_root, patt_root, &dims->dimension);
             h5_offsets[i] = parse_dimension (pvar_root, patt_root, &dims->local_offset);
             if (dims)
                 dims = dims -> next;
             if (debug && myrank==0) {
                 printf("\t%s[%d]: g(%llu):l(%llu):o(%llu)\n"
                       ,pvar->name,i
                       ,h5_globaldims[i]
                       ,h5_localdims[i]
                       ,h5_offsets[i]);
             }
        }
        sprintf(name, "_%s_gbdims", pvar->name);
        // save the global bounds information into hdf5 file

        h5_dataspace_id = H5Screate_simple (2, h5_gbglobaldims, NULL);
        h5_memspace_id = H5Screate_simple (2, h5_gblocaldims, NULL);

        H5Sselect_hyperslab (h5_dataspace_id, H5S_SELECT_SET
                            ,h5_gboffsets, h5_gbstrides, h5_gblocaldims,0);
    
        h5_dataset_id  = H5Dopen ( grp_ids[level], name);
        if (h5_dataset_id < 0) {
             h5_dataset_id = H5Dcreate (grp_ids[level]
                                       ,name
                                       ,H5T_STD_I64LE
                                       ,h5_dataspace_id
                                       ,H5P_DEFAULT);
        }
        if (h5_dataset_id > 0) {
            H5Dwrite (h5_dataset_id, H5T_STD_I64LE, h5_memspace_id
                     ,h5_dataspace_id, h5_plist_id, h5_gbdims);
            H5Dclose (h5_dataset_id);
        }
        H5Sclose(h5_dataspace_id); 
        H5Sclose(h5_memspace_id);
          
        h5_dataspace_id = H5Screate_simple (rank, h5_globaldims, NULL);
        if ( h5_dataspace_id < 0 && rank ==2) { 
            fprintf (stderr, "PHDF5 ERROR: cannot create dataspace for var: %s %d %llu %llu\n"\
                    ,pvar->name, rank, h5_globaldims[0], h5_globaldims[1]);
            return -1; 
         } 

         H5Sselect_hyperslab (h5_dataspace_id, H5S_SELECT_SET
                                  ,h5_offsets, h5_strides, h5_localdims, 0 
                                  );

         h5_dataset_id  = H5Dopen ( grp_ids[level], pvar->name);
         if ( h5_dataset_id < 0) { 
             h5_dataset_id = H5Dcreate (grp_ids[level]
                                       ,pvar->name
                                       ,h5_type_id
                                       ,h5_dataspace_id
                                       ,H5P_DEFAULT);
             if ( h5_dataset_id < 0) {
                fprintf (stderr, "PHDF5 ERROR: can not create dataset: %s!\n"
                       , pvar->name); 
                return -2; 
            } 
        } 
        h5_memspace_id = H5Screate_simple (rank, h5_localdims, NULL);
        if ( h5_memspace_id < 0) { 
            fprintf ( stderr, "PHDF5 ERROR: can not create h5_dataset for"
                      " var: %s\n"
                     ,pvar->name);
            return -1; 
        }
        H5Dwrite (h5_dataset_id, h5_type_id, h5_memspace_id
                           ,h5_dataspace_id, h5_plist_id, pvar->data
                           );
        H5Dclose (h5_dataset_id);
        H5Sclose (h5_dataspace_id);
        H5Sclose (h5_memspace_id);
        free (h5_gbdims);  
        free (h5_strides);  
    }
    else {
        h5_localdims = (hsize_t *) malloc (rank * sizeof(hsize_t));
        hid_t h5p_dset_id;
        enum ADIOS_FLAG is_timeindex = adios_flag_no;
        int  dimindex = 0;
        for ( i = 0; i < rank; i++) {
            h5_localdims [i] = parse_dimension(pvar_root, patt_root, &dims->dimension);
            if ( dims->dimension.is_time_index == adios_flag_yes) {
                  is_timeindex = adios_flag_yes;
                  dimindex = i;
            }
            dims = dims -> next;
        }
        h5_dataset_id  = H5Dopen (grp_ids[level], pvar->name);
        if (is_timeindex== adios_flag_no) {
            h5_dataspace_id = H5Screate_simple (rank, h5_localdims, NULL);
        }
        else {
            if (h5_dataset_id > 0) {
                h5_globaldims = (hsize_t *) malloc (rank * sizeof(hsize_t));
                h5_offsets = (hsize_t *) malloc (rank * sizeof(hsize_t));
                h5_strides = (hsize_t *) malloc (rank * sizeof(hsize_t));
                for (i=0;i<rank;i++) {
                    h5_offsets [i] = 0; 
                    h5_strides [i] = 1; 
                }
                h5_dataspace_id = H5Dget_space(h5_dataset_id);
                H5Sget_simple_extent_ndims (h5_dataspace_id);
                //fprintf(stderr, "var %s has time index %d %d \n"
                //       ,pvar->name, h5_offsets[1], h5_globaldims[1]); 
                H5Sget_simple_extent_dims (h5_dataspace_id, h5_globaldims, NULL);
                h5_offsets [dimindex] = h5_globaldims [dimindex];
                h5_globaldims [dimindex] = h5_globaldims [dimindex] + 1;
                H5Dextend (h5_dataset_id, h5_globaldims);
                h5_dataspace_id = H5Dget_space(h5_dataset_id);
                H5Sselect_hyperslab (h5_dataspace_id, H5S_SELECT_SET
                                ,h5_offsets, h5_strides, h5_localdims, 0
                                );
                h5_memspace_id = H5Screate_simple (rank, h5_localdims, NULL);
                //H5Sset_extent_simple (h5_dataspace_id, rank, h5_globaldims, NULL);
                //H5Sget_simple_extent_dims (h5_dataspace_id, h5_offsets, h5_globaldims);
                fprintf(stderr, "var %s has time index %llu %llu \n"
                       ,pvar->name, h5_offsets[1], h5_globaldims[1]); 
            }
            else {
                h5p_dset_id = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(h5p_dset_id, rank, h5_localdims);
                h5_dataspace_id = H5Screate_simple (rank, h5_localdims, NULL);
                h5_memspace_id = h5_dataspace_id;
            }
        }
       
        if ( h5_dataspace_id < 0) { 
            fprintf ( stderr, "PHDF5 ERROR: can not create memspace "
                              "for var: %s\n", pvar->name);
            return -1; 
        } 
        if ( h5_dataset_id < 0) {
            if (is_timeindex == adios_flag_yes) {
                h5_dataset_id = H5Dcreate (grp_ids[level]
                                      ,pvar->name
                                      ,h5_type_id
                                      ,h5_dataspace_id
                                      ,h5p_dset_id);
            } 
            else {
                h5_dataset_id = H5Dcreate (grp_ids[level]
                                      ,pvar->name
                                      ,h5_type_id
                                      ,h5_dataspace_id
                                      ,H5P_DEFAULT);
            } 
            if ( h5_dataset_id < 0) {
                fprintf ( stderr, "PHDF5 ERROR: can not create dataset: %s!\n", pvar->name); 
                return -2;
            } 
        } 
        if (myrank==0) 
        {
            if (is_timeindex == adios_flag_yes) {
               printf("dataspace: %d, memspace: %d\n",h5_memspace_id, h5_dataspace_id); 
               H5Dwrite (h5_dataset_id, h5_type_id,h5_memspace_id 
                        ,h5_dataspace_id, h5_plist_id, pvar->data
                        );
            }
            else {
               H5Dwrite (h5_dataset_id, h5_type_id, H5S_ALL
                        ,H5S_ALL, h5_plist_id, pvar->data
                        );
            }
        } 
        H5Dclose (h5_dataset_id);
        H5Sclose (h5_dataspace_id);
        H5Sclose (h5_memspace_id);
        free (h5_localdims);  
    }
    hw_gclose(grp_ids, level, adios_flag_yes);
    H5Tclose (h5_type_id);
    H5Pclose (h5_plist_id);
    return 0;
}

/*
 * Maps bp datatypes to h5 datatypes 
 *
 * The Mapping is according to HDF5 Reference Manual
 * (http://hdf.ncsa.uiuc.edu/HDF5/doc1.6/Datatypes.html)
 */
int getH5TypeId(enum ADIOS_DATATYPES type, hid_t* h5_type_id \
               ,enum ADIOS_FLAG fortran_flag) {
    int status=0;
    switch (type)
    {
        case adios_byte:
            *h5_type_id = H5Tcopy(H5T_NATIVE_CHAR);
            break;
        case adios_short:
            *h5_type_id = H5Tcopy(H5T_NATIVE_SHORT);
            break;
        case adios_integer:
            *h5_type_id = H5Tcopy(H5T_NATIVE_INT32);
            break;
        case adios_long:
                *h5_type_id = H5Tcopy(H5T_NATIVE_INT64);
            break;
        case adios_string:
            if (fortran_flag == adios_flag_yes)
                *h5_type_id = H5Tcopy(H5T_FORTRAN_S1);
            else if (fortran_flag == adios_flag_no)
                *h5_type_id = H5Tcopy(H5T_C_S1_g);
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
        case adios_complex:
        case adios_double_complex:
            fprintf(stderr, "ERROR in mapping ADIOS Data Types to HDF5: complex not supported yet.\n");
            status = -1;
            break;
        case adios_unsigned_byte:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UCHAR);
            break;
        case adios_unsigned_short:
            *h5_type_id = H5Tcopy(H5T_NATIVE_USHORT);
            break;
        case adios_unsigned_integer:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT32);
            break;
        case adios_unsigned_long:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT64);
            break;
        default:
            status = -1;
            break;
    }
    return status;
}

int adios_getsize(enum ADIOS_DATATYPES type, void *val) {

    switch (type)
    {
        case adios_byte:
            return 1;

        case adios_string:
            return strlen ((char *) val);

        case adios_short:
        case adios_unsigned_short:
            return 2;

        case adios_integer:
        case adios_unsigned_integer:
        case adios_long:
        case adios_unsigned_long:
            return 4;

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

// close group/dataset
// adios_flag_yes      group
// adios_flag_no       dataset

void hw_gclose (hid_t * grp_id, int level, enum ADIOS_FLAG flag) {
     int i;
     if (flag == adios_flag_unknown) {
         fprintf (stderr, "Unknown flag in hw_gclose!\n");
         return;
     }
     for ( i = 1; i <= level; i++) {
         if (i == level && flag == adios_flag_no)
             H5Dclose(grp_id[i]); 
         else
             H5Gclose(grp_id[i]); 

     }
}

// open group/dataset
// adios_flag_yes      group
// adios_flag_no       dataset
// adios_flag_unknown 

void hw_gopen (hid_t root_id, char * path, hid_t * grp_id, int * level, enum ADIOS_FLAG *flag) {

    //if (flag == adios_flag_unknown) {
    //   fprintf (stderr, "Unknown flag in hw_gclose!\n");
    //    return;
    //}
    int i, idx = 0, len = 0;
    char * pch, ** grp_name, * tmpstr;
    tmpstr= (char *)malloc(strlen(path)+1);
    strcpy (tmpstr, path);
    pch = strtok(tmpstr,"/");
    grp_name = (char **) malloc(NUM_GP);
    while ( pch!=NULL && *pch!=' ') {
        len = strlen(pch);
        grp_name[idx]  = (char *)malloc(len+1);
        grp_name[idx][0]='\0';
        strcat(grp_name[idx],pch);
        pch=strtok(NULL,"/");
        idx=idx+1;
    }
    *level = idx;
    grp_id [0] = root_id; 
    len = strlen(path);
    for ( i = 0; i < *level; i++) {
        //fprintf (stderr, "  -- open group [%s]\n", grp_name[i]);
        grp_id [i + 1] = H5Gopen (grp_id [i],grp_name [i]);
        if (grp_id [i + 1] < 0) {
            if ((i+1) == *level && (*flag == adios_flag_unknown)) {
                /* FIX: attribute may belong to a dataset or is a standalone thing */
                //MPI_Barrier(MPI_COMM_WORLD);
                grp_id [i + 1] = H5Dopen (grp_id [i],grp_name [i]);
                if (grp_id [i + 1] < 0) {
                    grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
                    *flag = adios_flag_yes; // it is a (new) group
                } else {
                    *flag = adios_flag_no; // it is an existing dataset
                }
                //MPI_Barrier(MPI_COMM_WORLD);
            }
            if ((i+1) == *level && (*flag == adios_flag_no)) {
                //MPI_Barrier(MPI_COMM_WORLD);
                grp_id [i + 1] = H5Dopen (grp_id [i],grp_name [i]);
                //MPI_Barrier(MPI_COMM_WORLD);
            }
            else {
                //MPI_Barrier(MPI_COMM_WORLD);
                grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
                //MPI_Barrier(MPI_COMM_WORLD);
//                printf("creat grp:name[%d]=%s id=%d level=%d\n", i,grp_name[i],grp_id[i+1],*level);
            }
        }
        if (grp_id [i + 1] < 0) { 
            fprintf (stderr, "PHDF5 ERROR: create group %s failed!\n", grp_name[i]);
            return; 
        }
    }
    for ( i = 0; i < *level; i++)
        if ( grp_name[i])
            free (grp_name[i]); 
    free (grp_name);
    free (tmpstr);
}
hsize_t parse_dimension(struct adios_var_struct *pvar_root,
                        struct adios_attribute_struct *patt_root, 
                        struct adios_dimension_item_struct * dim) {
    hsize_t dimsize;
    struct adios_var_struct *var_linked = NULL;
    struct adios_attribute_struct *attr_linked;
    if ( dim->var) {
        if ( dim->var->data){
            dimsize = *(int *)dim->var->data;
        }
    } else if ( dim->attr) {

        attr_linked = dim->attr;
        if (!attr_linked->var) {
            switch (attr_linked->type) {
                case adios_unsigned_byte:
                    dimsize = *(uint8_t *)attr_linked->value;
                    break;
                case adios_byte:
                    dimsize = *(int8_t *)attr_linked->value;
                    break;
                case adios_unsigned_short:
                    dimsize = *(uint16_t *)attr_linked->value;
                    break;
                case adios_short:
                    dimsize = *(int16_t *)attr_linked->value;
                    break;
                case adios_unsigned_integer:
                    dimsize = *(uint32_t *)attr_linked->value;
                    break;
                case adios_integer:
                    dimsize = *(int32_t *)attr_linked->value;
                    break;
                case adios_unsigned_long:
                    dimsize = *(uint64_t *)attr_linked->value;
                    break;
                case adios_long:
                    dimsize = *(int64_t *)attr_linked->value;
                    break;
                default:
                    fprintf (stderr, "Invalid datatype for array dimension on "
                            "var %s: %s\n"
                            ,attr_linked->name
                            ,adios_type_to_string_int (var_linked->type)
                            );
                    break;
            }
        }
        else {
            if (attr_linked->var->data) {
                dimsize = *(int *)attr_linked->var->data;
            }
        }
    }
    else {
        if (dim->is_time_index == adios_flag_yes)
            dimsize = 1;
        else
            dimsize = dim->rank;
    }
    return dimsize; 
}
#endif
