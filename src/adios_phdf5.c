#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mpi.h"
#include "hdf5.h"
#include "adios.h"
#include "adios_types.h"
#include "adios_bp_v1.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"

#define NUM_GP 24
///////////////////////////
// Function Declarations
///////////////////////////
int adios_getsize(enum ADIOS_DATATYPES type, void * val);
void hw_gopen (hid_t root_id, char * path, hid_t * grp_id, int * level, int flag) ;
void hw_gclose ( hid_t * grp_ids, int level, int flag);

int getH5TypeId(enum ADIOS_DATATYPES type, hid_t* h5_type_id, enum ADIOS_FLAG);

void *xrealloc (void *ptr, size_t size);

int hw_var ( hid_t root_id
             ,struct adios_var_struct *pvar_root
             ,struct adios_var_struct *pvar
             ,enum ADIOS_FLAG fortran_yes 
             ,int pid);

int hw_attribute ( hid_t root_id
             ,struct adios_var_struct *pvar_root
             ,struct adios_attribute_struct *pvar
             ,enum ADIOS_FLAG fortran_yes 
             ,int pid);

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
static void adios_var_to_comm_phdf5 (enum ADIOS_FLAG host_language_fortran
                                    ,void * data
                                    ,MPI_Comm * comm
                                    )
{
    if (data)
    {
        int t = *(int *) data;
        if (host_language_fortran == adios_flag_yes)
        {
            *comm = MPI_Comm_f2c (t);
        }
        else
        {
            *comm = *(MPI_Comm *) data;
        }
    }
    else
    {
        fprintf (stderr, "coordination-communication not provided. "
                         "Using MPI_COMM_WORLD instead\n"
                );
        *comm = MPI_COMM_WORLD;
    }
}
void adios_phdf5_init(const char *parameters
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
//    md->req = 0;
//    memset (&md->status, 0, sizeof (MPI_Status));
    md->rank = -1;
    md->size = 0;
    md->group_comm = MPI_COMM_NULL;
//    md->old_pg_root = 0;
//    md->old_vars_root = 0;

//    adios_buffer_struct_init (&md->b); 
}
enum ADIOS_FLAG adios_phdf5_should_buffer (struct adios_file_struct * fd
                            ,struct adios_method_struct * method
                            ,void * comm
                            )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct *)
                                                      method->method_data;
    char name[255];
    MPI_Info info = MPI_INFO_NULL;
    hid_t fapl_id;
    fapl_id = H5P_DEFAULT;
    adios_var_to_comm_phdf5 (fd->group->adios_host_language_fortran
                      ,comm
                      ,&md->group_comm
                      );
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
        fapl_id = H5Pcreate(H5P_FILE_ACCESS);
fprintf (stderr, "H5Pset_fapl_mpio not found on ewok so commented out\n");
//        H5Pset_fapl_mpio(fapl_id,md->group_comm,info);
    }
    else 
       md->group_comm=MPI_COMM_SELF;
    fd->group->process_id = md->rank;
    sprintf(name, "%s%s", method->base_path, fd->name);

    // create a new file. If file exists its contents will be overwritten. //
    md->fh = H5Fcreate (name, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    if (md->fh==-1)
    {
        md->fh = H5Fopen (name, H5F_ACC_RDWR, fapl_id);
    }
    else
    H5Pclose(fapl_id);
    md->root_id = H5Gopen(md->fh,"/");
    return 0;
}
int adios_phdf5_open(struct adios_file_struct *fd
                    ,struct adios_method_struct * method
                    )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct *)
                                                    method->method_data;
    return 1;
}

void adios_phdf5_write(struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,void * data
                     ,struct adios_method_struct * method
                     )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct *)
                                                       method->method_data;
    if(fd->mode != adios_mode_write && fd->mode != adios_mode_append)
    {
       printf("entering wrong mode!\n");
       return ;
    }
    v->data = data;
    hw_var (md->root_id, fd->group->vars, v, fd->group->adios_host_language_fortran,md->rank);
}

void adios_phdf5_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
}

void adios_phdf5_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffersize
                    ,struct adios_method_struct * method
                    )
{
    v->data = buffer;
    v->data_size = buffersize;
}
static void adios_phdf5_do_read (struct adios_file_struct * fd
                              ,struct adios_method_struct * method
                              )
{
}

void adios_phdf5_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_phdf5_data_struct * md = (struct adios_phdf5_data_struct*)
                                                  method->method_data;
      
    struct adios_attribute_struct * a = fd->group->attributes;

    while(a)
    {
//       if (a && md->rank==0) {
//            printf("\tatt name:%s\n",a->name);
            hw_attribute ( md->root_id, fd->group->vars, a
                          ,fd->group->adios_host_language_fortran
                          ,md->rank);
//         } 
         a = a->next;
    }

    if (md && md->fh && md->root_id) {
        H5Gclose (md->root_id);
        printf("rank:%d phdf5 file is closed;\n",md->rank);
    }
        H5Fclose (md->fh);
    //md->group_comm = -1;
    //md->fh = 0;
    //md->rank = -1;
    //md->size = 0;
}
void adios_phdf5_finalize (int mype, struct adios_method_struct * method)
{
    // nothing to do here
    if (adios_phdf5_initialized)
        adios_phdf5_initialized = 0;
}

void adios_phdf5_end_iteration (struct adios_method_struct * method)
{
}

void adios_phdf5_start_calculation (struct adios_method_struct * method)
{
}

void adios_phdf5_stop_calculation (struct adios_method_struct * method)
{
}

int hw_attribute ( hid_t root_id
		,struct adios_var_struct *pvar_root
		,struct adios_attribute_struct *patt
		,enum ADIOS_FLAG fortran_yes
		,int pid) {

    H5Eset_auto ( NULL, NULL);
    int i, rank = 0, level;
    char * name;
    herr_t  status; 
    hid_t   h5_plist_id, h5_type_id, h5_dataspace_id, h5_attribute_id, h5_memspace_id, grp_ids[NUM_GP];
    hsize_t * h5_globaldims, * h5_localdims, * h5_offsets, * h5_strides; 
    struct adios_dimension_struct * dims;
    struct adios_var_struct * var_linked;
     
    printf ("///////////////////////////\n");  
    h5_plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_plist_id=H5P_DEFAULT;   
fprintf (stderr, "H5Pset_dxpl_mpio not found on ewok so commented out\n");
    //H5Pset_dxpl_mpio(h5_plist_id, H5FD_MPIO_INDEPENDENT); 

    hw_gopen (root_id, patt->path,grp_ids, &level, 1);     

    if ( patt)
        printf ("write the attribute \t%s %d %d\n", patt->name, patt->type, patt->id);
    if (patt->type == -1) {
        var_linked = patt->var;
        if ( var_linked) {
            if (!(var_linked->data))
                printf ("\tERROR: invalid data in var_linked: %s %d\n", var_linked->name);  
            dims = var_linked->dimensions;
        }
        else 
            printf ("\tERROR: retrieving var_linked by var_id:  %s %d\n", var_linked->name);  
        getH5TypeId (var_linked->type, &h5_type_id, fortran_yes);

        if (!dims) {
             h5_dataspace_id = H5Screate ( H5S_SCALAR);
             if ( h5_dataspace_id < 0) {
                 printf (" ERROR in h5_dataspace_id in hw_attribute\n");  
                 return -2;
             }
             h5_attribute_id = H5Acreate ( grp_ids[level], patt->name
                                          ,h5_type_id,h5_dataspace_id,0);
             if ( h5_attribute_id < 0) {
                 printf (" ERROR: getting negative attribute_id in hw_attribute: %s %d\n",patt->name, h5_type_id);  
                 return -2;
             }
             if (pid == 0)	
                 H5Awrite ( h5_attribute_id, h5_type_id, var_linked->data);
             H5Aclose ( h5_attribute_id);
             H5Sclose ( h5_dataspace_id);
        }
        else {
             while ( dims) {
                 ++rank;
                 dims = dims->next;
             }
             h5_localdims = (hsize_t *) malloc (rank * sizeof(hsize_t));
             dims = var_linked->dimensions;
             for ( i = 0; i < rank; i++) {
                 if ( dims->dimension.rank == 0 && dims->dimension.id) { 
                     var_linked = adios_find_var_by_id (pvar_root , dims->dimension.id);
                     if ( var_linked) {
                         h5_localdims [i] = *(int *)var_linked->data;
                     }
                 }
                 else
                     h5_localdims [i] = dims->dimension.rank;
            }
#if  0
            h5_dataspace_id = H5Screate_simple(rank,h5_localdims, NULL);
            h5_attribute_id = H5Acreate ( grp_ids[level], patt->name
                                          ,h5_type_id,h5_dataspace_id,0);
            if ( h5_attribute_id < 0) {
                 printf (" ERROR: getting negative attribute_id in hw_attribute: %s %d\n",patt->name, h5_type_id);  
                 return -2;
            }
            if (pid == 0 && var_linked->data)
                 H5Awrite ( h5_attribute_id, h5_type_id, var_linked->data);
            H5Aclose ( h5_attribute_id);
            H5Sclose ( h5_dataspace_id);
#endif
            free (h5_localdims);
        } 
    }
    else {
        getH5TypeId (patt->type, &h5_type_id, fortran_yes);
        if ( h5_type_id < 0) {
            printf (" ERROR in getH5TypeId in hw_attribute\n");  
            return -2;
        }
        switch (patt->type) {
            case adios_string:
                h5_dataspace_id = H5Screate ( H5S_SCALAR);
                //printf("\tatt string size:%d\n", strlen(patt->value)+1);
                H5Tset_size ( h5_type_id, (strlen((char *)patt->value)+1));
                h5_attribute_id = H5Acreate ( grp_ids[level], patt->name 
                                        ,h5_type_id, h5_dataspace_id, 0);
                if (pid == 0)	
                    H5Awrite(h5_attribute_id, h5_type_id, patt->value);
                H5Aclose(h5_attribute_id);
                break; 
            case adios_integer: 
            case adios_long: 
                break; 
            default: 
                break;
        }
    } 
    H5Tclose (h5_type_id); 
    H5Pclose (h5_plist_id); 
    hw_gclose (grp_ids,level, 1);     
    return 0;
/*//////////////////////////////////////////////////////////////////////////
*/
} 

int hw_var ( hid_t root_id
             ,struct adios_var_struct *pvar_root
             ,struct adios_var_struct *pvar
             ,enum ADIOS_FLAG fortran_yes
             ,int pid) {

    H5Eset_auto ( NULL, NULL);

    int i, rank = 0, level;
    char * name;
    herr_t  status; 
    hid_t   h5_plist_id, h5_type_id, h5_dataspace_id, h5_dataset_id, h5_memspace_id, grp_ids[NUM_GP];
    hsize_t * h5_globaldims, * h5_localdims, * h5_offsets, * h5_strides; 
    struct adios_dimension_struct * dims = pvar->dimensions;
    struct adios_var_struct * var_linked;
     
    h5_plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_plist_id=H5P_DEFAULT;   
fprintf (stderr, "H5Pset_dxpl_mpio not found on ewok so commented out\n");
    //H5Pset_dxpl_mpio(h5_plist_id, H5FD_MPIO_INDEPENDENT); 

    getH5TypeId (pvar->type, &h5_type_id, fortran_yes);
    if ( h5_type_id <= 0) {
       printf (" ERROR in getH5TypeId in hw_var\n");  
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
    if (pvar->path)
        hw_gopen (root_id, pvar->path, grp_ids, &level, 1);

    if (!dims) {

        h5_dataspace_id = H5Screate(H5S_SCALAR);
        h5_dataset_id = H5Dcreate (grp_ids[level], pvar->name, h5_type_id
                                      ,h5_dataspace_id, H5P_DEFAULT
                                      );
        if ( h5_dataset_id < 0) { 
            printf ("ERROR in h5_dataset_id \n"); 
            return 0;
        }
        if ( pid==0)
            status = H5Dwrite (h5_dataset_id, h5_type_id, H5S_ALL
                                     ,H5S_ALL, h5_plist_id, pvar->data
                                     );
        status = H5Dclose (h5_dataset_id); 
        status = H5Sclose (h5_dataspace_id); 
        H5Tclose (h5_type_id);
        H5Pclose (h5_plist_id);
        hw_gclose (grp_ids,level, 1);
        return 0;
    }

    while ( dims) {
        ++rank;
        dims = dims->next;
    }

    h5_localdims = (hsize_t *) malloc (rank * sizeof(hsize_t));
    dims = pvar->dimensions;

    if ( dims->global_dimension.rank 
        || (dims->global_dimension.rank == 0  && dims->global_dimension.id)) {
        h5_strides = (hsize_t *) malloc (rank * sizeof(hsize_t));
        h5_globaldims = (hsize_t *) malloc (rank * sizeof(hsize_t));
        h5_offsets = (hsize_t *) malloc (rank * sizeof(hsize_t));
        // get the global/local/offset arrays       
        for ( i = 0; i < rank; i++) {

             h5_strides [i] = 1;
             //////////////////////////////////////////////////////////////////////
             if (dims->global_dimension.rank == 0 && dims->global_dimension.id) {
                  var_linked = adios_find_var_by_id (pvar_root , dims->global_dimension.id); 
                  h5_globaldims [i] = * ((int *)var_linked->data);
             }
             else
                  h5_globaldims [i] = dims->global_dimension.rank;

             //////////////////////////////////////////////////////////////////////
             if (dims->dimension.rank == 0 && dims->dimension.id) { 
                  var_linked = adios_find_var_by_id (pvar_root , dims->dimension.id);
                  if ( var_linked) {
                      h5_localdims [i] = *(int *)var_linked->data;
                  }
             }
             else
                  h5_localdims [i] = dims->dimension.rank;

             //////////////////////////////////////////////////////////////////////
             if (dims->local_offset.rank == 0 && dims->local_offset.id) { 
                  var_linked = adios_find_var_by_id (pvar_root , dims->local_offset.id); 
                  h5_offsets [i] = *(int *) var_linked->data;
             }
             else
                 h5_offsets [i] = dims->local_offset.rank;
  
             if ( dims) 
                 dims = dims -> next;
        }

        h5_dataspace_id = H5Screate_simple (rank, h5_globaldims, NULL);
        if ( h5_dataspace_id < 0) { 
            printf (" ERROR:  The system is out of memory, cannot create h5_dataset for var: %s %d\n"
                   ,pvar->name, rank);
            return -1; 
         } 

         status = H5Sselect_hyperslab (h5_dataspace_id, H5S_SELECT_SET
                                  ,h5_offsets, h5_strides, h5_localdims, 0 
                                  );

         h5_dataset_id  = H5Dopen ( grp_ids[level], pvar->name);
         if ( h5_dataset_id < 0) { 
             h5_dataset_id = H5Dcreate ( grp_ids[level], pvar->name, h5_type_id, h5_dataspace_id, H5P_DEFAULT);
             if ( h5_dataset_id < 0) {
                printf ("ERROR: The hdf5 is corrupted %d %s %d!\n", grp_ids[level], pvar->name, h5_type_id); 
                return -2; 
            } 
        } 
        h5_memspace_id = H5Screate_simple (rank, h5_localdims, NULL);
        if ( h5_memspace_id < 0) { 
            printf (" ERROR:  The system is out of memory, cannot create h5_dataset for var: %s\n", pvar->name);
            return -1; 
         } 
         status = H5Dwrite (h5_dataset_id, h5_type_id, h5_memspace_id 
                           ,h5_dataspace_id, h5_plist_id, pvar->data
                           );
        status = H5Dclose (h5_dataset_id);
        status = H5Sclose (h5_dataspace_id);
        status = H5Sclose (h5_memspace_id);
    }

    else {
        //////////////////////////////////////////////////////////////////////
        for ( i = 0; i < rank; i++) {
            if ( dims->dimension.rank == 0 && dims->dimension.id) { 
                var_linked = adios_find_var_by_id (pvar_root , dims->dimension.id);
                if ( var_linked) {
                    h5_localdims [i] = *(int *)var_linked->data;
                }
            }
            else
                h5_localdims [i] = dims->dimension.rank;
        }
        h5_dataspace_id = H5Screate_simple (rank, h5_localdims, NULL);
        if ( h5_dataspace_id < 0) { 
            printf (" ERROR:  The system is out of memory, cannot create h5_dataset for var: %s\n", pvar->name);
            return -1; 
         } 
        h5_dataset_id  = H5Dopen ( grp_ids[level], pvar->name);
        if ( h5_dataset_id < 0) { 
            h5_dataset_id = H5Dcreate ( grp_ids[level], pvar->name, H5T_NATIVE_FLOAT, h5_dataspace_id, H5P_DEFAULT);
            if ( h5_dataset_id < 0) {
                printf (" ERROR: The hdf5 is corrupted: %s!\n", pvar->name); 
                return -2;
            } 
        } 
        if ( pid==0)
           status = H5Dwrite (h5_dataset_id, h5_type_id, H5S_ALL
                              ,H5S_ALL, h5_plist_id, pvar->data
                              );
        
        status = H5Dclose (h5_dataset_id);
        status = H5Sclose (h5_dataspace_id);
    }
    hw_gclose(grp_ids, level, 1);
    H5Tclose (h5_type_id);
    H5Pclose (h5_plist_id);
    free (h5_localdims);  
}

/*
 * Maps bp datatypes to h5 datatypes 
 *
 * The Mapping is according to HDF5 Reference Manual
 * (http://hdf.ncsa.uiuc.edu/HDF5/doc1.6/Datatypes.html)
 */
int getH5TypeId(enum ADIOS_DATATYPES type, hid_t* h5_type_id \
               ,enum ADIOS_FLAG fortran_yes) {
    int size, status=0;
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
            if (fortran_yes == adios_flag_yes)
                *h5_type_id = H5Tcopy(H5T_FORTRAN_S1);
            else if (fortran_yes == adios_flag_no)
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
        case adios_unsigned_integer:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT32);
	    break;
        case adios_unsigned_long:
            *h5_type_id = H5Tcopy(H5T_NATIVE_UINT64);
	    break;
        default:
            status = -1;
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

void *xrealloc (void *ptr, size_t size) {
    register void *value = realloc (ptr, size);
    if (value == 0) {
        fprintf(stderr, "Virtual memory exhausted\n");
        exit(-1);
    }
    return value;
}
void hw_gclose (hid_t * grp_id, int level, int flag) {
     int i;
     for ( i = 1; i <= level; i++) {
         if (i == level && flag == 2)
             H5Dclose(grp_id[i]); 
         else
             H5Gclose(grp_id[i]); 

     }
}

// open group/dataset/

void hw_gopen (hid_t root_id, char * path, hid_t * grp_id, int * level, int flag) {

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
//    printf("parse group begin\n"); 
    for ( i = 0; i < *level; i++) {
        //printf("\t %s, %d\n", grp_name[i],grp_id[i]);
        grp_id [i + 1] = H5Gopen (grp_id [i],grp_name [i]);
        if (grp_id [i + 1] < 0) {
            if ((i+1) == *level && (flag == 2)) {
                grp_id [i + 1] = H5Dopen (grp_id [i],grp_name [i]);
                //if (grp_id [i + 1] < 0)  
                //    grp_id [i + 1] = H5Dcreate (grp_id [i], grp_name [i], 0);
            }
            else
                grp_id [i + 1] = H5Gcreate (grp_id [i], grp_name [i], 0);
        }
        if (grp_id [i + 1] < 0) { 
            printf ("ERROR: create group %s failed!\n", grp_name[i]);
            return; 
        }
    }
  //  printf("parse group end\n"); 
    for ( i = 0; i < *level; i++)
        if ( grp_name[i])
            free (grp_name[i]); 
    free (grp_name);
    free (tmpstr);
}
