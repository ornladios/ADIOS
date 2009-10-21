/*
 *   Read API for ADIOS BP format files 
 *   Copyright 2009, Oak Ridge National Laboratory
 */
#ifndef __ADIOS_READ_H__
#define __ADIOS_READ_H__

#include "adios_types.h"

#ifdef NOMPI
    /* Sequential processes can use the library compiled with -DNOMPI */
#   include "mpidummy.h"
#else
    /* Parallel applications should use MPI to communicate file info and slices of data */
#   include "mpi.h"
#endif

#include <stdint.h>

/***************/
/* C interface */
/***************/

/*************************/
/* Types used in the API */
/*************************/

typedef struct {
        int64_t  fh;                /* File handler                                                       */
        uint16_t groups_count;      /* Number of adios groups in file                                     */
        uint16_t vars_count;        /* Number of variables in all groups                                  */
        uint16_t attrs_count;       /* Number of attributes in all groups                                 */
        uint32_t tidx_start;        /* First timestep in file, usually 1                                  */
        uint32_t ntimesteps;        /* Number of timesteps in file. There is always at least one timestep */
        uint32_t version;           /* ADIOS BP version of file format                                    */
        uint64_t file_size;         /* Size of file in bytes                                              */
        int      endianness;        /* 0: little endian, 1: big endian (you do not need to care)          */
        char     ** group_namelist; /* Names of the adios groups in the file (cf. groups_count)           */
} ADIOS_FILE;

typedef struct {
        int64_t  gh;                /* Group handler                                           */
        int      grpid;             /* group index (0..ADIOS_FILE.groups_count-1)              */
        uint16_t vars_count;        /* Number of variables in this adios group                 */
        char     ** var_namelist;   /* Variable names in a char* array                         */
        uint16_t attrs_count;       /* Number of attributes in this adios group                */
        char     ** attr_namelist;  /* Attribute names in a char* array                        */
        ADIOS_FILE * fp;            /* pointer to the parent ADIOS_FILE struct                 */
} ADIOS_GROUP;

typedef struct {
        int        grpid;           /* group index (0..ADIOS_FILE.groups_count-1)                     */
        int        varid;           /* variable index (0..ADIOS_GROUP.var_count-1)                    */
        enum ADIOS_DATATYPES type;  /* type of variable                                               */
        int        ndim;            /* number of dimensions, 0 for scalars                            */
        uint64_t * dims;            /* size of each dimension                                         */
        int        timedim;         /* -1: variable has no timesteps in file, >=0: which dimension is time     */
        void     * value;           /* value of a scalar variable, NULL for array.                    */
        void     * gmin;            /* minimum value in an array variable, = value for a scalar       */
        void     * gmax;            /* maximum value of an array variable (over all timesteps)        */
        void     * mins;            /* minimum per each timestep (array of timestep elements)         */
        void     * maxs;            /* maximum per each timestep (array of timestep elements)         */
} ADIOS_VARINFO;

/** Functions that return a pointer to some data structures (fopen, gopen), return NULL
    on error and set adios_errno to a non-zero value and writes an error string.  
    You have direct access to that string so you can print it. 
    Do not write anything into it, please. 
    Only the last error message is always available. 
*/
extern int adios_errno;
char *adios_errmsg();

/** Open an adios file.
 *  IN:  fname    pathname of file to be opened
 *       comm     the MPI communicator of all processes that want to read data from the file
 *                 if compile with -DNOMPI, pass any integer here. 
 *  RETURN:       pointer to an ADIOS_FILE struct, NULL on error (sets adios_errno)
 */
ADIOS_FILE * adios_fopen (const char * fname, MPI_Comm comm);

/** Close an adios file.
 *  It will free the content of the underlying data structures and the fp pointer itself.
 *  IN:   fp       pointer to an ADIOS_FILE struct 
 *  RETURN: 0 OK, !=0 on error (also sets adios_errno)
 */
int adios_fclose (ADIOS_FILE *fp);


/** Open an adios group. Usually there is one adios group in a file, 
 *  but there can be more than one. 
 *  IN:  fp       pointer to an (opened) ADIOS_FILE struct
 *       grpname  name of the group
 *  RETURN:       pointer to an ADIOS_GROUP struct, NULL on error (sets adios_errno)
 */
ADIOS_GROUP * adios_gopen (ADIOS_FILE *fp, const char * grpname);

/** Open a group by index 
 *       grpid    index of group (0..fp->groups_count-1)
 *                in fp->group_namelist of ADIOS_FILE struct
 */
ADIOS_GROUP * adios_gopen_byid (ADIOS_FILE *fp, int grpid);

/** Close an adios group. 
 *  To free the data structures allocated at gopen, you need to call this function
 *  IN:  gp       pointer to an (opened) ADIOS_GROUP struct
 *  RETURN: 0 OK, !=0 on error (also sets adios_errno)
 */
int adios_gclose (ADIOS_GROUP *gp);


/** Inquiry about one variable in a group.
 *  This function does not read anything from the file but processes info
 *  already in memory after fopen and gopen.
 *  It allocates memory for the ADIOS_VARINFO struct and content, so 
 *  you need to free resources later with adios_free_varinfo().
 *
 *  IN:  gp       pointer to an (opened) ADIOS_GROUP struct
 *       varname  name of the variable
 *  RETURN:       pointer to and ADIOS_VARINFO struct, NULL on error (sets adios_errno)
 */
ADIOS_VARINFO * adios_inq_var (ADIOS_GROUP *gp, const char * varname);

/** Inquiry a variable by index
 *       varid    index of variable (0..gp->vars_count-1)
 *                in gp->vars_namelist of ADIOS_GROUP struct
 */
ADIOS_VARINFO * adios_inq_var_byid (ADIOS_GROUP *gp, int varid);

/** Free memory used by an ADIOS_VARINFO struct */
void adios_free_varinfo (ADIOS_VARINFO *cp);

/** Read a variable (slice) from the file.
 *  You need to allocate the memory for the data. 
 *  IN:  gp        pointer to an (opened) ADIOS_GROUP struct
 *       varname   name of the variable
 *       start     array of offsets to start reading in each dimension
 *       readsize  number of data elements to read in each dimension
 *  OUT: data      data of the variable
 *  RETURN: the number of bytes read, <0 on error, sets adios_errno too
 */
int64_t adios_read_var (ADIOS_GROUP    * gp, 
                        const char     * varname, 
                        const uint64_t * start,
                        const uint64_t * readsize, 
                        void           * data);

/** Read a variable by index 
 *       varid    index of variable (0..gp->vars_count-1)
 *                in gp->vars_namelist of ADIOS_GROUP struct
 */
int64_t adios_read_var_byid (ADIOS_GROUP * gp, int varid,
                             const uint64_t * start, const uint64_t * readsize, 
                             void * data);

/** Get an attribute in a group.
 *  This function does not read anything from the file but processes info
 *  already in memory after fopen and gopen.
 *  The memory for the data is allocated within the library.
 *  You can use free() to free the memory after use. 
 *
 *  IN:  gp       pointer to an (opened) ADIOS_GROUP struct
 *       attrname name of the attribute
 *  OUT: type     adios type of attribute (see enum ADIOS_DATATYPES in adios_types.h)
 *       size     memory size of value (n+1 for a string of n characters)
 *       data     pointer to the value. You need to cast it afterward according to the type. 
 *  RETURN: 0 OK, error: set and return adios_errno
 */
int adios_get_attr (ADIOS_GROUP           * gp,
                    const char            * attrname,
                    enum ADIOS_DATATYPES  * type,
                    int                   * size,
                    void                 ** data);

/** Convenience function to get an attribute by name 
 *       attrid   index of attribute (0..gp->attrs_count-1)
 *                in gp->attr_namelist of ADIOS_GROUP struct
 */
int adios_get_attr_byid (ADIOS_GROUP * gp, int attrid, enum ADIOS_DATATYPES * type, 
                         int * size, void ** data); 

/** Return the name of an adios type */
const char * adios_type_to_string (enum ADIOS_DATATYPES type);

/** Return the memory size of one data element of an adios type.
 *  If the type is adios_string, and the second argument is
 *  the string itself, it returns strlen(data)+1.
 *  For other types, it does not care about data and returns
 *  the size occupied by one element. 
 */
int adios_type_size(enum ADIOS_DATATYPES type, void *data);


/*********************/
/* Fortran interface */
/*********************/

/* In fortran, you do not need to include this header file. 
   Just link the code with the read library.
   adios_errno is returned by each function as extra parameter.
*/

/*
void FC_FUNC_(adiosf_errmsg, ADIOSF_LASTERRMSG) (char *msg, int msg_len);

void FC_FUNC_(adiosf_fopen, ADIOSF_FOPEN) ( int64_t * fp,
                    char * fname,
                    void * comm,
                    int * err,
                    int fname_len
                  );

void FC_FUNC_(adiosf_fclose, ADIOSF_FCLOSE) ( int64_t * fp, int * err);

void FC_FUNC_(adiosf_inq_file, ADIOSF_INQ_FILE) ( int64_t * fp,
                       int * groups_count,
                       int * vars_count,
                       int * attrs_count,
                       int * tstart,
                       int * ntsteps,
                       void * gnamelist,
                       int * err,
                       int gnamelist_len);

void FC_FUNC_(adiosf_gopen, ADIOSF_GOPEN) (int64_t * fp, int64_t * gp, 
                   char * grpname, int *err, int grpname_len);

void FC_FUNC_(adiosf_gclose, ADIOSF_GCLOSE) ( int64_t * gp, int * err);

void FC_FUNC_(adiosf_inq_group, ADIOSF_INQ_GROUP) (int64_t * gp, 
                       int     * vcnt, 
                       void    * vnamelist, 
                       int     * acnt, 
                       void    * anamelist, 
                       int     * err, 
                       int       vnamelist_len, 
                       int       anamelist_len);

void FC_FUNC_(adiosf_inq_var, ADIOSF_INQ_VAR) (int64_t  * gp, 
                     char     * varname,
                     int      * type,
                     int      * ndim,
                     uint64_t * dims,
                     int      * timedim,
                     int      * err,
                     int varname_len);

void FC_FUNC_(adiosf_read_var, ADIOSF_READ_VAR) (int64_t   * gp,
                      char      * varname,
                      uint64_t  * start,
                      uint64_t  * readsize, 
                      void      * data, 
                      int64_t   * read_bytes,
                      int  varname_len);

void FC_FUNC_(adiosf_get_varminmax, ADIOSF_GET_VARMINMAX) (int64_t * gp,
                           char    * varname,
                           void    * value,
                           void    * gmin,
                           void    * gmax,
                           void    * mins,
                           void    * maxs,
                           int     * err,
                           int  varname_len);

void FC_FUNC_(adiosf_get_attr, ADIOSF_GET_ATTR) (int64_t * gp,
                      char * attrname,
                      void * attr,
                      int * err,
                      int attrname_len);

void FC_FUNC_(adiosf_inq_attr, ADIOSF_INQ_ATTR) (int64_t * gp,
                      char * attrname,
                      int * type,
                      int * size,
                      int * err,
                      int attrname_len);
*/
#endif
