#ifndef __BP_READ_H__
#define __BP_READ_H__

#ifdef NOMPI
#   include "mpidummy.h"
#else
#   include "mpi.h"
#endif

#include <stdint.h>
/*
header file for the subsetting read routines
March 2009, ORNL
*/

/* Types used in the API */

typedef struct {
        uint16_t namelist_true;     /* true: return namelists too in adios_inq_file()   */
        uint16_t vars_count;        /* Number of variables in this adios group          */
        char     ** var_namelist;   /* Variable names in a char* array                  */
        uint16_t attrs_count;       /* Number of attributes in this adios group         */
        char     ** attr_namelist;  /* Attribute names in a char* array                 */
}BP_GROUP_INFO;

typedef struct {
        uint16_t namelist_true;     /* true: return namelists too in adios_inq_file() */
        uint16_t groups_count;      /* Number of adios groups in file */
        uint16_t vars_count;        /* Number of variables in all groups */
        uint16_t attrs_count;       /* Number of attributes in all groups */
        uint32_t tidx_start;        /* First timestep in file, usually 1 */
        uint32_t tidx_stop;         /* Last timestep in file. There is always at least one timestep */
        uint32_t version;           /* ADIOS BP version of file format */
        uint32_t file_size;         /* Size of file in bytes */
        char     ** group_namelist; /* Names of the adios groups in the file (cf. groups_count) */
}BP_FILE_INFO;

/***************/
/* C interface */
/***************/

/*
        IN:  fname
             comm
        OUT: fh_p
*/
int adios_fopen ( int64_t * fh_p,
                  const char * fname,
                  MPI_Comm comm
                );

/*
        IN: fh 
*/
int adios_fclose ( int64_t fh);

/* Inquire information about a file:
   Usage in C code:
   
        int64_t fh;
        BP_FILE_INFO finfo;
        ...
        adios_fopen (&fh, <filename>, <MPI communicator>);
        adios_init_fileinfo ( &finfo, 1);
        adios_inq_file (fh, &finfo);
        ...
        adios_free_fileinfo ( &finfo );

 */ 
int adios_inq_file ( int64_t fh_p, BP_FILE_INFO * pfinfo); 
void adios_init_fileinfo  (BP_FILE_INFO * pfinfo, int flag);
void adios_free_fileinfo  (BP_FILE_INFO * pfinfo);
void adios_print_fileinfo (BP_FILE_INFO * pfinfo);

/*
        IN:  fh
             grpname 
        OUT: gh_p 
*/
int adios_gopen ( int64_t fh,
                  int64_t * gh_p,
                  const char * grpname);

/*
        IN:  fh
*/
int adios_gclose ( int64_t gh);

/*
        IN:  gh_p
*/
/* Inquire information about a group:
   Usage in C code:

        int64_t gh;
        BP_GROUP_INFO ginfo; 
        ...
        adios_gopen (fh, &gh, finfo.group_namelist[<index>]);
        adios_init_groupinfo ( &ginfo, 1);
        adios_inq_group (gh, &ginfo);
        ...
        adios_free_groupinfo (&ginfo);
   
*/
int  adios_inq_group (int64_t gh, BP_GROUP_INFO *);
void adios_init_groupinfo  (BP_GROUP_INFO * pginfo, int flag);
void adios_free_groupinfo  (BP_GROUP_INFO * pginfo);
void adios_print_groupinfo (BP_GROUP_INFO * pginfo);

/*
        IN:  gh
             varname 
             start     array of offsets to start reading in each dimension
             readsize  number of data elements to read in each dimensio
             timestep  timestep to read in
        OUT: var 
*/
int64_t adios_get_var (int64_t gh,
                       const char * varname,
                       void * var, 
                       const int  * start,
                       const int  * readsize, 
                       int    timestep);

int64_t adios_get_attr (int64_t gh_p
                       ,const char * attrname
                       ,void * data);

/*
        IN:  gh
             varname 
        OUT: type           enum ADIOS_DATATYPES, see adios_types.h; see bp_type_to_string() below
             ndim           Number of dimensions
             is_timebased   1: variable has timesteps in file, 0: read only with timestep=1
             dims           array of size ndim, dimensions of the variable
        RETURN: 0 OK
                -1 variable not found
                -2 gh is null
                -3 file handler to gh is invalid
*/
int adios_inq_var (int64_t gh, const char * varname,
                   int * type,
                   int * ndim,
                   int * is_timebased,
                   int * dims);

int adios_inq_attr (int64_t gh_p
                   ,const char * attrname
                   ,int * type
                   ,int * size);


const char * bp_type_to_string (int type);


/*********************/
/* Fortran interface */
/*********************/

void adios_fopen_ ( int64_t * fh,
                    char * fname,
                    void * comm,
                    int * err,
                    int fname_len
                  );

void adios_fclose_ ( int64_t * fh, int * err);

void adios_inq_file_ ( int64_t * fh_p,
                       int * groups_count,
                       int * vars_count,
                       int * attrs_count,
                       int * tstart,
                       int * tstop,
                       void * gnamelist,
                       int * err,
                       int gnamelist_len);

void adios_gopen_ (int64_t * fh, int64_t * gh_p, 
                   char * grpname, int *err, int grpname_len);

void adios_gclose_ ( int64_t * gh, int * err);

void adios_inq_group_ (int64_t * gh, int *vcnt, void *v, int * err, int);

void adios_get_var_ (int64_t * gh,
                     char * varname,
                     void * var, 
                     int  * start,
                     int  * readsize, 
                     int  * timestep,
                     int64_t  * err,
                     int  varname_len);

void adios_inq_var_ (int64_t * gh_p, char * varname,
                     int * type,
                     int * ndim,
                     int * is_timebased,
                     int * dims,
                     int * err,
                     int varname_len);

void adios_get_attr_ (int64_t * gh_p,
                      char * attrname,
                      void * attr,
                      int * err,
                      int attrname_len);

void adios_inq_attr_ (int64_t * gh_p,
                      char * attrname,
                      int * type,
                      int * size,
                      int * err,
                      int attrname_len);

#endif
