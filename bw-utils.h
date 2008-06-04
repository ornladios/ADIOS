#include "binpack-general.h"
#include "binpack-utils.h"

/***************************************/
/* public functions (Fortan interface) */
/***************************************/
// write a dataset (array)
void bw_dset_ (void *file_hd,int *startidx,int *endidx,char* path,char* name
              ,void* val, int* type, int* rank,int* dims
              );

// write a scalar (single value or string)
void bw_scalar_ (void *file_hd, int *startidx,int *endidx,char *str
                ,char *name, void *val, int *type
                );

// write a numeric attribute for a data item
void bw_an_d_ (void *file_hd, int*,int *,char *str,char *aname, void *aval
              ,int *type
              );

// write a numeric attribute for a group
void bw_an_g_ (void *file_hd, int*,int *,char *str,char *aname, void *aval
              ,int *type
              );

// write a string attribute for a data item
void bw_as_d_ (void *file_hd, int*,int *,char *str,char *aname, char *aval);

// write a string attribute for a group
void bw_as_g_ (void *file_hd, int*,int *,char *str,char *aname, char *aval);

/********************************/
/* public functions C interface */
/********************************/
// write a scalar
void bw_scalar (void*,int,int *, char *str, char *name, void *val
               ,enum vartype_t type
               );

// write a dataset (array)
void bw_dset (void* fbp,int start, int * end,char *path, char *name, void *val
             ,enum vartype_t type, int rank
             ,struct adios_bp_dimension_struct * dims
             );

///////////////////////////////////////////
// write a dataset attribute (string)
void bw_attr_str_ds (void*,int, int *,char *str, char *aname, char* aval);

// write a dataset attribute (numeric)
void bw_attr_num_ds (void*,int, int *,char *str, char *aname
                    ,void *, enum vartype_t
                    );

// write a group attribute (string)
void bw_attr_str_gp (void*,int, int *,char *str, char *aname, char* aval);

// write a group attribute (numeric)
void bw_attr_num_gp (void*,int, int *,char *str, char *aname
                    ,void *, enum vartype_t
                    );

///////////////////////////////////////////
// how big will this scalar element be?
int bcalsize_scalar (char *dir_name,char* name,enum vartype_t type, void * val);

// how big will this dataset (array) be?
int bcalsize_dset (char *dir_name,char* name,enum vartype_t type,int rank
                  ,struct adios_bp_dimension_struct * dims
                  );

// how big will this string attribute be?
int bcalsize_attr_str (char *dir_name,char*, char*);

// how big will this number attribute be?
int bcalsize_attr_num (char *dir_name,char* name,enum vartype_t type, void * val);

/**********************/
/* internal functions */
/**********************/
void bwrite (void *val,int size,int number, void *buf,int *idx);
void bw_stringtag (void*,int*,enum TAG_t tag, char *name);
void bw_scalartag (void*,int*,enum TAG_t tag, void *val, enum vartype_t type);
void bw_dsettag (void*,int*, enum TAG_t tag, void *val, enum vartype_t type
                ,int rank
                ,struct adios_bp_dimension_struct * dims
                );
void bw_set_write (int a);
void bw_fopen_ (char * filename, long long * file_hd);
void bw_fclose_ (long long * fptr);
