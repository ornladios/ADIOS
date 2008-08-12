#include "hdf5.h"

int hw_makeh5(char *filename);
void hw_string_attr( hid_t parent_id, const char *name,const char *value);
void hw_scalar_attr(hid_t parent_id, const char *name, void *value,enum vartype_t type);

void hw_dset (hid_t root_id,struct adios_bp_element_struct* dims);
//void hw_dset(hid_t root_id, char *dirstr, char *name, void *data, enum vartype_t type,int rank, struct adios_bp_dimension_struct * dims, struct adios_bp_dimension_struct * global_dims);
void hw_scalar(hid_t root_id, char *dirstr, char *name, void *val, enum vartype_t type, int append);
void hw_attr_str_ds(hid_t root_id, char *dirstr, char *aname, char* aval);
void hw_attr_num_ds(hid_t root_id, char *dirstr, char *aname, void *, enum vartype_t type);
void hw_attr_str_gp(hid_t root_id, char *dirstr, char *aname, char* aval);
void hw_attr_num_gp(hid_t root_id, char *dirstr, char *aname, void *, enum vartype_t type);
void hw_dataset(hid_t parent_id, char* name, void* data,enum vartype_t type, int ndims, hsize_t* dims);
void hw_dataset1(hid_t parent_id, char* name, void* data,enum vartype_t type, int ndims, hsize_t* dims, int append);
int bp_getH5TypeId(enum vartype_t type, hid_t* h5_type_id, void * val);
