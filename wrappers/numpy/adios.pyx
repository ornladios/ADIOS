"""
 ADIOS is freely available under the terms of the BSD license described
 in the COPYING file in the top level directory of this source distribution.

 Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
"""

import numpy as np
cimport numpy as np

import mpi4py.MPI as MPI 
cimport mpi4py.MPI as MPI

import cython
cimport cython

## ==========
## ADIOS Exported Functions
## ==========

from libc.stdint cimport uint32_t, int64_t, uint64_t
from libc.stdlib cimport malloc, free

cdef extern from "adios_types.h":
    ctypedef enum ADIOS_DATATYPES:
        adios_unknown
        adios_byte
        adios_short
        adios_integer
        adios_long
        adios_unsigned_byte
        adios_unsigned_short
        adios_unsigned_integer
        adios_unsigned_long
        adios_real
        adios_double
        adios_long_double
        adios_string
        adios_complex
        adios_double_complex

    ctypedef enum ADIOS_BUFFER_ALLOC_WHEN:
        ADIOS_BUFFER_ALLOC_UNKNOWN
        ADIOS_BUFFER_ALLOC_NOW
        ADIOS_BUFFER_ALLOC_LATER

    ctypedef enum ADIOS_FLAG:
        pass

cdef extern from "adios.h":
    ctypedef struct MPI_Comm:
        pass

    ctypedef char* const_char_ptr "const char*"

    cdef int adios_init (char * config, MPI_Comm)
    
    cdef int adios_finalize (int mype)
    
    cdef int adios_open (int64_t * fd,
                         char * group_name,
                         char * name, 
                         char * mode,
                         void * comm)
    
    cdef int adios_group_size (int64_t fd_p,
                               uint64_t data_size,
                               uint64_t * total_size)
    
    cdef int adios_write (int64_t fd_p,
                          char * name,
                          void * var)
    
    cdef int adios_read (int64_t fd_p,
                         char * name,
                         void * buffer,
                         uint64_t buffer_size)

    cdef int adios_close(int64_t fd_p)
    
    cdef int adios_init_noxml (MPI_Comm)
    
    cdef int adios_allocate_buffer (ADIOS_BUFFER_ALLOC_WHEN when,
                                    uint64_t buffer_size)
    
    cdef int adios_declare_group (int64_t * id,
                                  char * name,
                                  char * time_index,
                                  ADIOS_FLAG stats)
    
    cdef int adios_define_var (int64_t group_id,
                               char * name,
                               char * path,
                               ADIOS_DATATYPES type,
                               char * dimensions,
                               char * global_dimensions,
                               char * local_offsets)
    
    cdef int adios_define_attribute (int64_t group,
                                     char * name,
                                     char * path,
                                     ADIOS_DATATYPES type,
                                     char * value,
                                     char * var)
    
    cdef int adios_select_method (int64_t group,
                                  char * method,
                                  char * parameters,
                                  char * base_path)

cdef extern from "adios_read.h":
    ctypedef struct MPI_Comm:
        pass

    ctypedef struct ADIOS_FILE:
        uint64_t fh
        int      groups_count
        int      vars_count
        int      attrs_count
        int      tidx_start
        int      ntimesteps
        int      version
        uint64_t file_size
        int      endianness
        char     ** group_namelist
        void     * internal_data

    ctypedef struct ADIOS_GROUP:
        uint64_t gh
        int      grpid
        int      vars_count
        char     ** var_namelist
        int      attrs_count
        char     ** attr_namelist
        ADIOS_FILE * fp
        int      timestep
        int      lasttimestep

    ctypedef struct ADIOS_VARINFO:
        int        grpid
        int        varid
        ADIOS_DATATYPES   type
        int        ndim
        uint64_t * dims
        int        timedim
        int        characteristics_count
        void     * value
        void     * gmin
        void     * gmax
        double   * gavg
        double   * gstd_dev
        void     ** mins
        void     ** maxs
        double   ** avgs
        double   ** std_devs


    cdef ADIOS_FILE * adios_fopen (char *, MPI_Comm)
    cdef int adios_fclose (ADIOS_FILE *)

    cdef ADIOS_GROUP * adios_gopen (ADIOS_FILE *, char *)
    cdef ADIOS_GROUP * adios_gopen_byid (ADIOS_FILE *, int)
    cdef int adios_gclose (ADIOS_GROUP *)

    cdef ADIOS_VARINFO * adios_inq_var (ADIOS_GROUP *, char *)
    cdef ADIOS_VARINFO * adios_inq_var_byid (ADIOS_GROUP *, int)
    cdef void adios_free_varinfo(ADIOS_VARINFO *)
    cdef int64_t adios_read_var_byid (ADIOS_GROUP *, int,
                             uint64_t *, uint64_t *,
                             void *)

    cdef char * adios_errmsg()
    cdef int adios_errno


## ==========
## ADIOS Enum
## ==========

class DATATYPE(object):
    unknown = -1
    byte = 0
    short = 1
    integer = 2
    long = 4
    unsigned_byte = 50
    unsigned_short = 51
    unsigned_integer = 52
    unsigned_long = 54
    real = 5
    double = 6
    long_double = 7
    string = 9
    complex = 10
    double_complex = 11

class FLAG(object):
    UNKNOWN = 0
    YES = 1
    NO = 2
    
class BUFFER_ALLOC_WHEN(object):
    UNKNOWN = 0
    NOW = 1
    LATER = 2
    
## ==========
## ADIOS Write API
## ==========

cpdef init(char * config, MPI.Comm comm = MPI.COMM_WORLD):
    return adios_init(config, comm.ob_mpi)

cpdef int64_t open(char * group_name,
                   char * name,
                   char * mode,
                   MPI.Comm comm = MPI.COMM_WORLD):
    cdef int64_t fd
    cdef int result
    result = adios_open(&fd, group_name, name, mode, comm.ob_mpi)
    return fd

cpdef int64_t set_group_size(int64_t fd_p, uint64_t data_size):
    cdef uint64_t total_size
    cdef int result
    result = adios_group_size(fd_p, data_size, &total_size)
    return total_size

cpdef int write (int64_t fd_p, char * name, np.ndarray val):
    cdef np.ndarray val_
    if val.flags.contiguous:
        val_ = val
    else:
        val_ = np.array(val, copy=True)

    return adios_write (fd_p, name, <void *> val_.data)

cpdef int write_int (int64_t fd_p, char * name, int val):
    return adios_write (fd_p, name, &val)

cpdef int write_long (int64_t fd_p, char * name, long val):
    return adios_write (fd_p, name, &val)

cpdef int write_float (int64_t fd_p, char * name, float val):
    return adios_write (fd_p, name, &val)

cpdef int read(int64_t fd_p, char * name, np.ndarray val):
    assert val.flags.contiguous, 'Only contiguous arrays are supported.'
    print "Reading ... ", val.itemsize * val.size, "(bytes)"
    return adios_read(fd_p, name, <void *> val.data, val.itemsize * val.size)

cpdef int close(int64_t fd_p):
    return adios_close(fd_p)

cpdef finalize(int mype = 0):
    return adios_finalize(mype)

## ==========
## ADIOS No-XML API
## ==========
cpdef int init_noxml(MPI.Comm comm = MPI.COMM_WORLD):
    return adios_init_noxml(comm.ob_mpi)

cpdef int allocate_buffer(int when,
                          uint64_t buffer_size):
    return adios_allocate_buffer(<ADIOS_BUFFER_ALLOC_WHEN> when,
                                 buffer_size)

cpdef int64_t declare_group(char * name,
                            char * time_index,
                            int stats):
    cdef int64_t id = 0
    adios_declare_group (&id,
                         name,
                         time_index,
                         <ADIOS_FLAG> stats)
    return id

cpdef int define_var(int64_t group_id,
                     char * name,
                     char * path,
                     int type,
                     char * dimensions,
                     char * global_dimensions,
                     char * local_offsets):
    return adios_define_var(group_id,
                            name, path,
                            <ADIOS_DATATYPES> type,
                            dimensions,
                            global_dimensions,
                            local_offsets)

cpdef int define_attribute (int64_t group,
                            char * name,
                            char * path,
                            int type,
                            char * value,
                            char * var):
    return adios_define_attribute (group,
                                   name,
                                   path,
                                   <ADIOS_DATATYPES> type,
                                   value,
                                   var)

cpdef int select_method (int64_t group,
                         char * method,
                         char * parameters,
                         char * base_path):
    return adios_select_method (group,
                                method,
                                parameters,
                                base_path)


## ==========
## ADIOS Read API
## ==========

cpdef type adios2nptype(ADIOS_DATATYPES t):
    cdef type ntype = None
    if t == adios_byte:
        ntype = np.int8
    elif t == adios_short:
        ntype = np.int16
    elif t == adios_integer:
        ntype = np.int32
    elif t == adios_long:
        ntype = np.int64
    elif t == adios_unsigned_byte:
        ntype = np.uint8
    elif t == adios_unsigned_short:
        ntype = np.uint16
    elif t == adios_unsigned_integer:
        ntype = np.uint32
    elif t == adios_unsigned_long:
        ntype = np.uint64
    elif t == adios_real:
        ntype = np.float32
    elif t == adios_double:
        ntype = np.float64
    elif t == adios_long_double:
        ntype = np.float128
    elif t == adios_complex:
        ntype = np.complex64
    elif t == adios_double_complex:
        ntype = np.complex128
    else:
        ntype = None

    return ntype

"""
cpdef int np2adiostype(np.dtype t):
    cdef int atype = -1
    if t.type == np.float64:
        atype = DATATYPE.double
    else:
        atype = -1

    return atype
"""

cdef printAdiosFile(ADIOS_FILE * f):
    print '%15s : %lu' % ('fh', f.fh)
    print '%15s : %d' % ('groups_count', f.groups_count)
    print '%15s : %d' % ('vars_count', f.vars_count)
    print '%15s : %d' % ('attrs_count', f.attrs_count)      
    print '%15s : %d' % ('tidx_start', f.tidx_start)       
    print '%15s : %d' % ('ntimesteps', f.ntimesteps)       
    print '%15s : %d' % ('version', f.version)          
    print '%15s : %lu' % ('file_size', f.file_size)
    print '%15s : %d' % ('endianness', f.endianness)
    print '%15s : %s' % ('group_namelist', [f.group_namelist[i] for i in range(f.groups_count)])

cdef printAdiosGroup(ADIOS_GROUP * g):
    print '%15s : %lu' % ('gh', g.gh)
    print '%15s : %d' % ('grpid', g.grpid)
    print '%15s : %d' % ('vars_count', g.vars_count)
    print '%15s : %s' % ('var_namelist', [g.var_namelist[i] for i in range(g.vars_count)])
    print '%15s : %d' % ('attrs_count', g.attrs_count)      
    print '%15s : %s' % ('attr_namelist', [g.attr_namelist[i] for i in range(g.attrs_count)])
    print '%15s : %lu' % ('fp', <unsigned long> g.fp)

cdef printAdiosVariable(ADIOS_VARINFO * v):
    print '%15s : %d' % ('grpid', v.grpid)
    print '%15s : %d' % ('varid', v.varid)
    print '%15s : %s' % ('type', adios2nptype(v.type))
    print '%15s : %d' % ('ndim', v.ndim)
    print '%15s : %s' % ('dims', [v.dims[i] for i in range(v.ndim)])
    print '%15s : %d' % ('timedim', v.timedim)

cdef adios2scalar(ADIOS_DATATYPES t, void * val):
    if t == adios_byte :
        return (<char *> val)[0]
    elif t == adios_short:
        return (<short *> val)[0]
    elif t == adios_integer:
        return (<int *> val)[0]
    elif t == adios_long:
        return (<long *> val)[0]
    elif t == adios_unsigned_byte:
        return (<unsigned char *> val)[0]
    elif t == adios_unsigned_short:
        return (<unsigned short *> val)[0]
    elif t == adios_unsigned_integer:
        return (<unsigned int *> val)[0]
    elif t == adios_unsigned_long:
        return (<unsigned long *> val)[0]
    elif t == adios_real:
        return (<float *> val)[0]
    elif t == adios_double:
        return (<double *> val)[0]
    elif t == adios_long_double:
        return (<long double *> val)[0]
    else:
        return None

## ==========
## ADIOS Class Definition
## ==========
    
cdef class AdiosFile:
    """ Private Memeber """
    cpdef ADIOS_FILE * fp

    """ Public Memeber """
    cpdef public bytes name
    cpdef public int groups_count
    cpdef public int vars_count
    cpdef public int attrs_count
    cpdef public int tidx_start
    cpdef public int ntimesteps
    cpdef public int version
    cpdef public int file_size
    cpdef public int endianness
    
    cpdef public dict group
    
    def __init__(self, char * fname, MPI.Comm comm = MPI.COMM_WORLD):
        self.fp = NULL
        self.group = {}
        
        self.fp = adios_fopen(fname, comm.ob_mpi)
        assert self.fp != NULL, 'Not an open file'

        self.name         = fname.split('/')[-1]  ## basename
        self.groups_count = self.fp.groups_count
        self.vars_count   = self.fp.vars_count  
        self.attrs_count  = self.fp.attrs_count 
        self.tidx_start   = self.fp.tidx_start  
        self.ntimesteps   = self.fp.ntimesteps  
        self.version      = self.fp.version     
        self.file_size    = self.fp.file_size   
        self.endianness   = self.fp.endianness  
    
        cdef AdiosGroup g
        for grpname in [self.fp.group_namelist[i] for i in range(self.groups_count)]:
            g = AdiosGroup(self, grpname)
            self.group[grpname] = g

    def __del__(self):
        self.close()
            
    cpdef close(self):
        assert self.fp != NULL, 'Not an open file'
        for g in self.group.values():
            g.close()
        adios_fclose(self.fp)
        self.fp = NULL
        
    cpdef printself(self):
        assert self.fp != NULL, 'Not an open file'
        print '=== AdiosFile ==='
        print '%15s : %lu' % ('fp', <unsigned long> self.fp)
        printAdiosFile(self.fp)


cdef class AdiosGroup:
    """ Private Memeber """
    cdef AdiosFile file
    cdef ADIOS_GROUP * gp

    """ Public Memeber """
    cpdef public bytes name
    cpdef public int grpid
    cpdef public int vars_count
    cpdef public int attrs_count
    cpdef public int timestep
    cpdef public int lasttimestep
    
    cpdef public dict var
    
    def __init__(self, AdiosFile file, char * name):
        self.file = file
        self.var = {}
        
        self.gp = adios_gopen(self.file.fp, name)
        assert self.gp != NULL, 'Not an open group'

        self.name         = name
        self.grpid        = self.gp.grpid        
        self.vars_count   = self.gp.vars_count   
        self.attrs_count  = self.gp.attrs_count  
        self.timestep     = self.gp.timestep     
        self.lasttimestep = self.gp.lasttimestep 
        
        cdef AdiosVariable v
        for varname in [self.gp.var_namelist[i] for i in range(self.vars_count)]:
            v = AdiosVariable(self, varname)
            self.var[varname] = v

    def __del__(self):
        self.close()

    cpdef close(self):
        assert self.gp != NULL, 'Not an open file'
        for v in self.var.values():
            v.close()
        adios_gclose(self.gp)
        self.gp = NULL
        
    cpdef printself(self):
        assert self.gp != NULL, 'Not an open file'
        print '=== AdiosGroup ==='
        print '%15s : %lu' % ('gp', <unsigned long> self.gp)
        printAdiosGroup(self.gp)
        
cdef class AdiosVariable:
    """ Private Memeber """
    cdef AdiosGroup group
    cdef ADIOS_VARINFO * vp

    """ Public Memeber """
    cpdef public bytes name
    cpdef public int varid
    cpdef public type type
    cpdef public int ndim
    cpdef public tuple dims
    cpdef public int timedim
    cpdef public int characteristics_count
    
    def __init__(self, AdiosGroup group, char * name):
        self.group = group
        self.vp = NULL

        self.vp = adios_inq_var(self.group.gp, name)
        assert self.group.gp != NULL, 'Not an open group'

        self.name                  = name
        self.varid                 = self.vp.varid                
        self.type                  = adios2nptype(self.vp.type)
        self.ndim                  = self.vp.ndim                 
        self.timedim               = self.vp.timedim              
        self.characteristics_count = self.vp.characteristics_count
        
        self.dims = tuple([self.vp.dims[i] for i in range(self.vp.ndim)])
        
    def __del__(self):
        self.close()

    cpdef close(self):
        assert self.vp != NULL, 'Not an open file'
        adios_free_varinfo(self.vp)
        self.vp = NULL
        
    cpdef read(self, tuple offset = (), tuple count = ()):
        cdef type ntype = adios2nptype(self.vp.type)
        assert ntype is not None, 'Data type is not supported yet'
        
        cdef list lshape = [self.vp.dims[i] for i in range(self.vp.ndim)]
        cdef np.ndarray npshape = np.array(lshape, dtype=np.int64)
        
        cdef np.ndarray npoffset
        if len(offset) == 0:
            npoffset = npshape.copy()
            npoffset.fill(0)
        else:
            npoffset = np.array(offset, dtype=np.int64)
        
        cdef np.ndarray npcount
        if len(count) == 0:
            npcount = npshape - npoffset
        else:
            npcount = np.array(count, dtype=np.int64)

        assert npshape.ndim == npoffset.ndim, 'Offset dimension mismatch'
        assert npshape.ndim == npcount.ndim, 'Count dimension mismatch.'
        assert (npshape - npoffset > npcount).all(), 'Count is larger than shape.'

        cdef np.ndarray var = np.zeros(npcount, dtype=ntype)
        cdef int64_t nbytes = adios_read_var_byid(
            self.group.gp, 
            self.vp.varid, 
            <uint64_t *> npoffset.data, 
            <uint64_t *> npcount.data, 
            <void *> var.data
            )

        if nbytes < 0:
            print "[WARNING] bytes read :", nbytes

        return var

    """ Print self """
    cpdef printself(self):
        assert self.vp != NULL, 'Not an open file'
        print '=== AdiosVariable ==='
        print '%15s : %lu' % ('vp', <unsigned long> self.vp)
        printAdiosVariable(self.vp)
