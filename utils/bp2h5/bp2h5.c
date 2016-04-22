/*  
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS bp2h5v2 utility 
 *  read all variables and attributes from 
 *    all groups in a BP file and output this to a hdf5 file
 *
 * This is a sequential program.
 */


#ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#endif


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <errno.h>
#include <limits.h>   // LONG_MAX
#include <math.h>     // NAN
#include <libgen.h>   // basename
#include <regex.h>    // regular expression matching
#include <fnmatch.h>  // shell pattern matching

#include "adios_read.h"
#include "adios_types.h"

#include "hdf5.h"
//#include "hdf5_hl.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

typedef int bool;
#define false 0
#define true  1

bool noindex = false;              // do no print array indices with data
bool printByteAsChar = false;      // print 8 bit integer arrays as string
bool formatgiven = false;           // true if format string is provided as argument
int  ncols = 6; // how many values to print in one row (only for -p)
char format[32];            // format string for one data element (e.g. %6.2f)

hid_t       HDF5_FILE;


//#define MAX_BUFFERSIZE 81 
#define MAX_BUFFERSIZE 10485760
#define MAX_DIMS 20
#define GMAX 100 
#define DEBUG 0
#define verbose 1

/* Support for complex types */
typedef struct {
    float re;   /*real part*/
    float im;   /*imaginary part*/
} complex_real_t;

typedef struct {
    double re;   /*real part*/
    double im;   /*imaginary part*/
} complex_double_t;

hid_t complex_real_id, complex_double_id;




int  istart[MAX_DIMS], icount[MAX_DIMS], ndimsspecified=0;


int bp_getH5TypeId(enum ADIOS_DATATYPES type, hid_t* h5_type_id);
int getTypeInfo( enum ADIOS_DATATYPES adiosvartype, int* elemsize);
int readVar(ADIOS_GROUP *gp, ADIOS_VARINFO *vi, const char * name);
const char * value_to_string (enum ADIOS_DATATYPES type, void * data, int idx);
char** bp_dirparser(char *str, int *nLevel);

int main (int argc, char ** argv)  
{
    int         gidx, i, j;
    MPI_Comm    comm_dummy = MPI_COMM_WORLD;  /* MPI_Comm is defined through adios_read.h */
    uint64_t    count[MAX_DIMS];
    herr_t      h5_err;
    char        h5name[256],aname[256],fname[256];
    int         h5i, level;
    hid_t       grp_id [GMAX+1], space_id;
    hid_t       att_id;
    char        ** grp_name;
    hid_t       h5_type_id;


    if (argc < 2) {
        printf("Usage: %s <BP-file> <HDF5-file>\n", argv[0]);
        return 1;
    }

    MPI_Init(&argc, &argv);
    h5_err = H5Eset_auto(NULL, NULL );
    ADIOS_FILE * f = adios_fopen (argv[1], comm_dummy);
    HDF5_FILE = H5Fcreate(argv[2],H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* create the complex types for HDF5 */
    complex_real_id = H5Tcreate (H5T_COMPOUND, sizeof (complex_real_t));
    H5Tinsert (complex_real_id, "real", HOFFSET(complex_real_t,re), H5T_NATIVE_FLOAT);
    H5Tinsert (complex_real_id, "imaginary", HOFFSET(complex_real_t,im), H5T_NATIVE_FLOAT);

    complex_double_id = H5Tcreate (H5T_COMPOUND, sizeof (complex_double_t));
    H5Tinsert (complex_double_id, "real", HOFFSET(complex_double_t,re), H5T_NATIVE_DOUBLE);
    H5Tinsert (complex_double_id, "imaginary", HOFFSET(complex_double_t,im), H5T_NATIVE_DOUBLE);

    if (f == NULL) {
        if (DEBUG) printf ("%s\n", adios_errmsg());
	return -1;
    }
    /* For all groups */
    for (gidx = 0; gidx < f->groups_count; gidx++) {
        if (DEBUG) printf("Group %s:\n", f->group_namelist[gidx]);
        ADIOS_GROUP * g = adios_gopen (f, f->group_namelist[gidx]);
        if (g == NULL) {
            if (DEBUG) printf ("%s\n", adios_errmsg());
            return -1;
        }
/* First create all of the groups */
        grp_id [0] = HDF5_FILE;
        for (i = 0; i < g->vars_count; i++) {
             strcpy(h5name,g->var_namelist[i]);
             grp_name = bp_dirparser (h5name, &level);
             for (j = 0; j < level-1; j++) {
                grp_id [j + 1] = H5Gopen (grp_id [j], grp_name [j]);
                if (grp_id [j + 1] < 0) {
                   grp_id [j + 1] = H5Gcreate (grp_id [j], grp_name [j], 0);
                }
             }
             for (j=1; j<level; j++) {
                  H5Gclose(grp_id[j]);
             }
        }
/* Now we can write data into these scalars */        
        /* For all variables */
        if (DEBUG) printf("  Variables=%d:\n", g->vars_count);
        for (i = 0; i < g->vars_count; i++) {
             ADIOS_VARINFO * v = adios_inq_var_byid (g, i);

            uint64_t total_size = adios_type_size (v->type, v->value);
            for (j = 0; j < v->ndim; j++)
                total_size *= v->dims[j];
            strcpy(h5name,g->var_namelist[i]);
            if (DEBUG) printf("    %-9s  %s", adios_type_to_string(v->type), g->var_namelist[i]);
            h5_err = bp_getH5TypeId (v->type, &h5_type_id);
            if (v->type==adios_string) H5Tset_size(h5_type_id,strlen(v->value)); 
            if (v->ndim == 0) {
                /* Scalars do not need to be read in, we get it from the metadata
                   when using adios_inq_var */
                if (DEBUG) printf(" = %s\n", value_to_string(v->type, v->value, 0));
                 // add the hdf5 dataset, these are scalars
                for (h5i = 0;h5i<MAX_DIMS;h5i++) 
                   count[0] = 0;
                count[0] = 1; // we are writing just 1 element, RANK=1
                h5_err = bp_getH5TypeId (v->type, &h5_type_id);
                H5LTmake_dataset(HDF5_FILE,h5name,1,count,h5_type_id,v->value);
            } else {

                    h5_err = readVar(g, v,  h5name);
            }
            adios_free_varinfo (v);
        } /* variables */

        /* For all attributes */
        if (DEBUG) printf("  Attributes=%d:\n", g->attrs_count);
        for (i = 0; i < g->attrs_count; i++) {
            enum ADIOS_DATATYPES atype;
            int  asize;
	    void *adata;
            adios_get_attr_byid (g, i, &atype, &asize, &adata);
            grp_name = bp_dirparser (g->attr_namelist[i], &level);
            strcpy(aname,grp_name[level-1]); 
// the name of the attribute is the last in the array
// we then need to concat the rest together
            strcpy(fname,"/");
            for (j=0;j<level-1;j++) {
              strcat(fname,grp_name[j]); 
            }
            h5_err = bp_getH5TypeId (atype, &h5_type_id);

            // let's create the attribute
            if (atype==adios_string) H5Tset_size(h5_type_id,strlen(adata)); 
            space_id = H5Screate(H5S_SCALAR); // just a scalar
            att_id = H5Acreate(HDF5_FILE, g->attr_namelist[i], h5_type_id, space_id,H5P_DEFAULT);
            h5_err = H5Awrite(att_id, h5_type_id, adata);
            h5_err = H5Aclose(att_id);
            h5_err = H5Sclose(space_id);

            if (DEBUG) printf("    %-9s  %s = %s\n", adios_type_to_string(atype), 
                    g->attr_namelist[i], value_to_string(atype, adata, 0));
            free(adata);
        } /* attributes */

        adios_gclose (g);
    } /* groups */

    adios_fclose (f);
    h5_err =  H5Fclose(HDF5_FILE);

    MPI_Finalize();
    return (int)h5_err;
}


const char * value_to_string (enum ADIOS_DATATYPES type, void * data, int idx)
{
    static char s [100];
    s [0] = 0;


    switch (type)
    {
        case adios_unsigned_byte:
            sprintf (s, "%u", ((uint8_t *) data)[idx]);
            break;

        case adios_byte:
            sprintf (s, "%d", ((int8_t *) data)[idx]);
            break;

        case adios_short:
            sprintf (s, "%hd", ((int16_t *) data)[idx]);
            break;

        case adios_unsigned_short:
            sprintf (s, "%hu", ((uint16_t *) data)[idx]);
            break;

        case adios_integer:
            sprintf (s, "%d", ((int32_t *) data)[idx]);
            break;

        case adios_unsigned_integer:
            sprintf (s, "%u", ((uint32_t *) data)[idx]);
            break;

        case adios_long:
            sprintf (s, "%" PRId64, ((int64_t *) data)[idx]);
            break;

        case adios_unsigned_long:
            sprintf (s, "%" PRIu64, ((uint64_t *) data)[idx]);
            break;

        case adios_real:
            sprintf (s, "%g", ((float *) data)[idx]);
            break;

        case adios_double:
            sprintf (s, "%lg", ((double *) data)[idx]);
            break;

        case adios_long_double:
            sprintf (s, "%Lg", ((long double *) data)[idx]);
            break;

        case adios_string:
            return (char*) ((char *)data+idx);
            break;

        case adios_string_array:
            return (char*) *((char **)data+idx);
            break;

        case adios_complex:
            sprintf (s, "(%g, %g)", 
                    ((float *) data)[2*idx], ((float *) data)[2*idx+1]);
            break;

        case adios_double_complex:
            sprintf (s, "(%lg, %lg)", 
                    ((double *) data)[2*idx], ((double *) data)[2*idx+1]);
            break;

        default:
        	break;
    }

    return s;
}

char** bp_dirparser(char *str, int *nLevel)
{
  char **grp_name;
  char *pch;
  int idx = 0, len=0;
  char *tmpstr;
  tmpstr= (char *)malloc(1*(strlen(str)+1));
  strcpy(tmpstr,str);
  pch = strtok(tmpstr,"/");
  grp_name = (char **)malloc(GMAX);
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

int readVar(ADIOS_GROUP *gp, ADIOS_VARINFO *vi, const char * name)
{
  int i,j;
  uint64_t start_t[MAX_DIMS], count_t[MAX_DIMS]; // processed <0 values in start/count
  uint64_t s[MAX_DIMS], c[MAX_DIMS]; // for block reading of smaller chunks
  hsize_t  h5_start[MAX_DIMS], h5_count[MAX_DIMS], h5_stride[MAX_DIMS];
  hsize_t  h5_dims[MAX_DIMS];
  uint64_t nelems;         // number of elements to read
  int      elemsize;            // size in bytes of one element
  uint64_t st, ct;
  void     *data;
  uint64_t sum;           // working var to sum up things
  int      maxreadn;     // max number of elements to read once up to a limit (10MB of data)
  int      actualreadn;       // our decision how much to read at once
  int      readn[MAX_DIMS];   // how big chunk to read in in each dimension?
  int64_t  bytes_read;     // retval from adios_get_var()
  int      incdim;            // used in incremental reading in
  hid_t    dataset, global_memspace;
  hid_t    local_memspace, h5_ndim ;
  hid_t    h5_err;
  hid_t    h5_type_id;


  if (getTypeInfo(vi->type, &elemsize)) {
    fprintf(stderr, "Adios type %d (%s) not supported in bpls. var=%s\n", 
	    vi->type, adios_type_to_string(vi->type), name);
    return 10;
  }

   h5_err = bp_getH5TypeId (vi->type, &h5_type_id);
   h5_ndim = (hsize_t) vi->ndim;
   for (j=0;j<h5_ndim;j++)
       h5_dims[j] = vi->dims[j];
// create the hdf5 dataspace.
 
  // create the counter arrays with the appropriate lengths
  // transfer start and count arrays to format dependent arrays
  for (j=0; j<vi->ndim; j++)  {
      icount[j]=-1;
      h5_stride[j]= (hsize_t) 1;
  }
  nelems = 1;
  for (j=0; j<vi->ndim; j++) {
    if (istart[j] < 0)  // negative index means last-|index|
      st = vi->dims[j]+istart[j];
    else
      st = istart[j];
    if (icount[j] < 0)  // negative index means last-|index|+1-start
      ct = vi->dims[j]+icount[j]+1-st;
    else
      ct = icount[j];
    if (verbose>2) 
      printf("    j=%d, st=%" PRIu64 " ct=%" PRIu64 "\n", j, st, ct);
    start_t[j] = st;
    count_t[j] = ct;
    nelems *= ct;
    if (verbose>1) 
      printf("    s[%d]=%" PRIu64 ", c[%d]=%" PRIu64 ", n=%" PRIu64 "\n", j, start_t[j], j, count_t[j], nelems);
  }
  if (verbose>1) {
    printf(" total size of data to read = %" PRIu64 "\n", nelems*elemsize);
  }
  maxreadn = MAX_BUFFERSIZE/elemsize;
  if (nelems < maxreadn)
	  maxreadn = nelems;
  // special case: string. Need to use different elemsize
  if (vi->type == adios_string) {
    if (vi->value)
      elemsize = strlen(vi->value)+1;
    maxreadn = elemsize;
  }
  // allocate data array
  data = (void *) malloc (maxreadn*elemsize); // SAK: don't want the +8.... 
  //+8 for just to be sure

  // determine strategy how to read in:
  //  - at once
  //  - loop over 1st dimension
  //  - loop over 1st & 2nd dimension
  //  - etc
  if (verbose>1) printf("Read size strategy:\n");
  sum = (uint64_t) 1;
  actualreadn = (uint64_t) 1;
  for (i=vi->ndim-1; i>=0; i--) {
    if (sum >= (uint64_t) maxreadn) {
      readn[i] = 1;
    } else {
      readn[i] = maxreadn / (int)sum; // sum is small for 4 bytes here
      // this may be over the max count for this dimension
      if (readn[i] > count_t[i]) 
	readn[i] = count_t[i];
    }
    if (verbose>1) printf("    dim %d: read %d elements\n", i, readn[i]);
    sum = sum * (uint64_t) count_t[i];
    actualreadn = actualreadn * readn[i];
  }
  if (verbose>1) printf("    read %d elements at once, %" PRId64 " in total (nelems=%" PRId64 ")\n", actualreadn, sum, nelems);


  // init s and c
  for (j=0; j<vi->ndim; j++) {
    s[j]=start_t[j];
    c[j]=readn[j];
    h5_count[j] = (hsize_t) c[j];
    h5_start[j] = (hsize_t) s[j];
  }

  // read until read all 'nelems' elements
  sum = 0;
  while (sum < nelems) {

    // how many elements do we read in next?
    actualreadn = 1;
    for (j=0; j<vi->ndim; j++) 
      actualreadn *= c[j];

    if (verbose>2) {
      printf("adios_read_var name=%s ", name);
      //PRINT_DIMS("  start", s, vi->ndim, j); 
      //PRINT_DIMS("  count", c, vi->ndim, j); 
      printf("  read %d elems\n", actualreadn);
    }

    // read a slice finally
    bytes_read = adios_read_var_byid (gp, vi->varid, s, c, data); 

    if (bytes_read < 0) {
      fprintf(stderr, "Error when reading variable %s. errno=%d : %s \n", name, adios_errno, adios_errmsg());
      free(data);
      return 11;
    }

// now we must place this inside the hdf5 file... we know the offset (s)
// we know the count, c
// we know the global rank v->ndim
// we know the global dimensions v->dims

// get the hdf5 calls for writing the hyperslab.....

    dataset = H5Dopen(HDF5_FILE,name);
    if (dataset> 0) {
       global_memspace = H5Dget_space(dataset);
       //hid_t rank_old = H5Sget_simple_extent_ndims(global_memspace);
       hsize_t *maxdims = (hsize_t *) malloc (h5_ndim * sizeof (hsize_t));
       h5_err = H5Sget_simple_extent_dims(global_memspace,maxdims,NULL);
       free(maxdims);
    } else {
       global_memspace = H5Screate_simple (h5_ndim, h5_dims, NULL);
       hid_t cparms = H5Pcreate(H5P_DATASET_CREATE);
       h5_err = H5Pset_chunk(cparms,h5_ndim,h5_count);
   
       dataset = H5Dcreate(HDF5_FILE,name,h5_type_id, global_memspace,cparms);
       H5Pclose(cparms);
    }

    local_memspace = H5Screate_simple (h5_ndim, h5_count, NULL);
     for (j=0;j<vi->ndim;j++) 

    h5_err =  H5Sselect_hyperslab (global_memspace, H5S_SELECT_SET
                                  ,h5_start, h5_stride, h5_count, NULL);

    h5_err = H5Dwrite(dataset,h5_type_id,local_memspace,global_memspace,H5P_DEFAULT, data);
    H5Sclose(local_memspace);
    H5Dclose(dataset);
    //H5Tclose(h5_type_id);

    if (verbose>2) printf("  read %" PRId64 " bytes\n", bytes_read);

    // print slice

    // prepare for next read
    sum += actualreadn;
    incdim=1; // largest dim should be increased 
    for (j=vi->ndim-1; j>=0; j--) {
      if (incdim==1) {
	if (s[j]+c[j] == start_t[j]+count_t[j]) {
	  // reached the end of this dimension
	  s[j] = start_t[j];
	  c[j] = readn[j];
          h5_count[j] = (hsize_t) c[j];
          h5_start[j] = (hsize_t) s[j];
           
	  incdim = 1; // next smaller dim can increase too
	} else {
	  // move up in this dimension up to total count
	  s[j] += readn[j];
          h5_start[j] = (hsize_t) s[j];
	  if (s[j]+c[j] > start_t[j]+count_t[j]) {
	    // do not reach over the limit
	    c[j] = start_t[j]+count_t[j]-s[j];
            h5_count[j] = (hsize_t) c[j];
	  }
	  incdim = 0;
	}
      }
    }
  } // end while sum < nelems
  H5Sclose(global_memspace);

  free(data);
  return (int)h5_err;
}

int getTypeInfo( enum ADIOS_DATATYPES adiosvartype, int* elemsize)
{
  switch(adiosvartype) {
  case adios_unsigned_byte:
    *elemsize = 1;
    break;
  case adios_byte:
    *elemsize = 1;
    break;
  case adios_string:
    *elemsize = 1;
    break;

  case adios_unsigned_short:  
    *elemsize = 2;
    break;
  case adios_short:
    *elemsize = 2;
    break;

  case adios_unsigned_integer:
    *elemsize = 4;
    break;
  case adios_integer:    
    *elemsize = 4;
    break;

  case adios_unsigned_long:
    *elemsize = 8;
    break;
  case adios_long:        
    *elemsize = 8;
    break;

  case adios_real:
    *elemsize = 4;
    break;

  case adios_double:
    *elemsize = 8;
    break;

  case adios_complex:  
    *elemsize = 8;
    break;

  case adios_double_complex:
    *elemsize = 16;
    break;

  case adios_long_double: // do not know how to print
    //*elemsize = 16;
  default:
    return 1;
  }
  return 0;
}

int bp_getH5TypeId(enum ADIOS_DATATYPES type, hid_t* h5_type_id)
{
    int status=0;

    switch (type)
    {
        case adios_byte:
            *h5_type_id = H5Tcopy(H5T_NATIVE_CHAR);
            break;
        case adios_string:
             *h5_type_id = H5Tcopy(H5T_C_S1);
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
            *h5_type_id = H5Tcopy(complex_real_id);
            break;
        case adios_double_complex:
            *h5_type_id = H5Tcopy(complex_double_id);
            //fprintf(stderr, "Error in mapping ADIOS Data Types to HDF5: complex not supported yet.\n");
            //status = -1;
            break;
        case adios_unknown:
        default:
            fprintf(stderr, "Error in mapping ADIOS Data Types to HDF5: unknown data type.\n");
            status = -1;
    }
    return status;
}

