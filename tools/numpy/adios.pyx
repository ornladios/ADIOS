"""
 ADIOS is freely available under the terms of the BSD license described
 in the COPYING file in the top level directory of this source distribution.

 Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
"""

import numpy as np
cimport numpy as np

import mpi4py.MPI as MPI 
cimport mpi4py.MPI as MPI 
##from mpi4py.mpi_c cimport * 
##from mpi4py.mpi_c import * 

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
    ctypedef char* const_char_ptr "const char*"
    cdef int adios_init (char * config)
    cdef int adios_finalize (int mype)
    cdef int adios_open (int64_t * fd, char * group_name, char * name, 
                         char * mode, void * comm)
    cdef int adios_group_size (int64_t fd_p, uint64_t data_size,
                               uint64_t * total_size)
    cdef int adios_write (int64_t fd_p, char * name, void * var)
    cdef int adios_read (int64_t fd_p, char * name, void * buffer, uint64_t buffer_size)

    cdef int adios_close (int64_t fd_p)
    
    cdef int adios_init_noxml ()
    cdef int adios_allocate_buffer (ADIOS_BUFFER_ALLOC_WHEN when, uint64_t buffer_size)
    cdef int adios_declare_group (int64_t * id, char * name, char * time_index, ADIOS_FLAG stats)
    cdef int adios_define_var (int64_t group_id, char * name, char * path, int type, char * dimensions, char * global_dimensions, char * local_offsets)
    cdef int adios_define_attribute (int64_t group, char * name, char * path, ADIOS_DATATYPES type, char * value, char * var)
    cdef int adios_select_method (int64_t group, char * method, char * parameters, char * base_path)

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

cpdef init(char * config):
    return adios_init(config)

cpdef int64_t open(char * group_name, char * name, char * mode, MPI.Comm comm = MPI.COMM_WORLD):
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
    ##assert val.flags.contiguous, 'Only contiguous arrays are supported.'
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
# adios_init_noxml
# adios_allocate_buffer
# adios_declare_group
# adios_define_var
# adios_define_attribute
# adios_select_method

cpdef int init_noxml():
    return adios_init_noxml()

cpdef int allocate_buffer(int when, uint64_t buffer_size):
    return adios_allocate_buffer(<ADIOS_BUFFER_ALLOC_WHEN> when, buffer_size)

cpdef int64_t declare_group(char * name, char * time_index, int stats):
    cdef int64_t id = 0
    adios_declare_group (&id, name, time_index, <ADIOS_FLAG> stats)
    return id

cpdef int define_var(int64_t group_id, char * name, char * path, int type,
                     char * dimensions,
                     char * global_dimensions,
                     char * local_offsets):
    return adios_define_var(group_id, name, path, type,
                            dimensions, global_dimensions, local_offsets)

cpdef int define_attribute (int64_t group, char * name, char * path,
                            int type, char * value, char * var):
    return adios_define_attribute (group, name, path,
                                   <ADIOS_DATATYPES> type, value, var)

cpdef int select_method (int64_t group,
                         char * method,
                         char * parameters,
                         char * base_path):
    return adios_select_method (group, method, parameters, base_path)


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
    cdef ADIOS_FILE * fp
    def __init__(self):
        self.fp = NULL

    def __init__(self, char * fname, MPI.Comm comm):
        self.open(fname, comm)
        
    def __init__(self, char * fname):
        self.open(fname)

    def __del__(self):
        self.close()

    cpdef open(self, char * fname, MPI.Comm comm = MPI.COMM_WORLD):
        self.fp = adios_fopen(fname, <MPI_Comm> comm.ob_mpi)

    cpdef close(self):
        assert self.fp != NULL, 'Not an open file'
        adios_fclose(self.fp)
        self.fp = NULL
        
    cpdef int groups_count(self):
        assert self.fp != NULL, 'Not an open file'
        return self.fp.groups_count;

    cpdef int vars_count(self):
        assert self.fp != NULL, 'Not an open file'
        return self.fp.vars_count;

    cpdef int attrs_count(self):
        assert self.fp != NULL, 'Not an open file'
        return self.fp.attrs_count;

    cpdef int tidx_start(self):
        assert self.fp != NULL, 'Not an open file'
        return self.fp.tidx_start

    cpdef int ntimesteps(self):
        assert self.fp != NULL, 'Not an open file'
        return self.fp.ntimesteps

    cpdef int version(self):
        assert self.fp != NULL, 'Not an open file'
        return self.fp.version

    cpdef uint64_t file_size(self):
        assert self.fp != NULL, 'Not an open file'
        return self.fp.file_size

    cpdef int endianness(self):
        assert self.fp != NULL, 'Not an open file'
        return self.fp.endianness

    cpdef group_namelist(self):
        return [self.fp.group_namelist[i] for i in range(self.fp.groups_count)]

    cpdef AdiosGroup group(self, char * grpname):
        assert self.fp != NULL, 'Not an open file'
        cdef AdiosGroup g = AdiosGroup(self)
        g.open(grpname)
        return g

    cpdef AdiosGroup group_byid(self, int grpid):
        assert self.fp != NULL, 'Not an open file'
        cdef AdiosGroup g = AdiosGroup(self)
        g.open_byid(grpid)
        return g

    cpdef printself(self):
        assert self.fp != NULL, 'Not an open file'
        print '=== AdiosFile ==='
        print '%15s : %lu' % ('fp', <unsigned long> self.fp)
        printAdiosFile(self.fp)


cdef class AdiosGroup:
    cdef AdiosFile file
    cdef ADIOS_GROUP * gp

    def __init__(self, AdiosFile file):
        self.file = file
        self.gp = NULL

    def __del__(self):
        self.close()

    cpdef open(self, char * name):
        assert self.file is not None, 'AdiosFile is not set'
        assert self.file.fp != NULL, 'Not an open file'
        self.gp = adios_gopen(self.file.fp, name)

    cpdef open_byid(self, int id):
        assert self.file is not None, 'AdiosFile is not set'
        assert self.file.fp != NULL, 'Not an open file'
        self.gp = adios_gopen_byid(self.file.fp, id)

    cpdef close(self):
        assert self.gp != NULL, 'Not an open group'
        adios_gclose(self.gp)
        self.file = None
        self.gp = NULL
            
    cpdef int grpid(self):
        assert self.gp != NULL, 'Not an open group'
        return self.gp.grpid

    cpdef int vars_count(self):
        assert self.gp != NULL, 'Not an open group'
        return self.gp.vars_count;

    cpdef var_namelist(self):
        assert self.gp != NULL, 'Not an open group'
        return [self.gp.var_namelist[i] for i in range(self.gp.vars_count)]

    cpdef int attrs_count(self):
        assert self.gp != NULL, 'Not an open group'
        return self.gp.attrs_count;

    cpdef attr_namelist(self):
        assert self.gp != NULL, 'Not an open group'
        return [self.gp.attr_namelist[i] for i in range(self.gp.attrs_count)]

    cpdef int timestep(self):
        assert self.gp != NULL, 'Not an open group'
        return self.gp.timestep

    cpdef int lasttimestep(self):
        assert self.gp != NULL, 'Not an open group'
        return self.gp.lasttimestep

    cpdef AdiosVariable variable_byid(self, int varid):
        assert self.gp != NULL, 'Not an open group'
        cdef AdiosVariable var = AdiosVariable(self)
        var.open_byid(varid)
        return var

    cpdef AdiosVariable variable(self, char * varname):
        assert self.gp != NULL, 'Not an open group'
        cdef AdiosVariable var = AdiosVariable(self)
        var.open(varname)
        return var

    cpdef np.ndarray read_var_byid(self, int varid, tuple offset = (), tuple count = ()):
        cdef AdiosVariable var = AdiosVariable(self)
        var.open_byid(varid)
        return var.read(offset, count)

    cpdef np.ndarray read_var(self, char * varname, tuple offset = (), tuple count = ()):
        cdef AdiosVariable var = AdiosVariable(self)
        var.open(varname)
        return var.read(offset, count)

    cpdef printself(self):
        assert self.gp != NULL, 'Not an open file'
        print '=== AdiosGroup ==='
        print '%15s : %lu' % ('gp', <unsigned long> self.gp)
        printAdiosGroup(self.gp)

cdef class AdiosVariable:
    cdef AdiosGroup group
    cdef ADIOS_VARINFO * vp

    def __init__(self, AdiosGroup group):
        self.group = group
        self.vp = NULL

    def __del__(self):
        self.close()

    def open(self, char * name):
        assert self.group is not None, 'AdiosGroup is not set'
        assert self.group.gp != NULL, 'Not an open group'
        self.vp = adios_inq_var(self.group.gp, name)

    def open_byid(self, int id):
        assert self.group is not None, 'AdiosGroup is not set'
        assert self.group.gp != NULL, 'Not an open group'
        self.vp = adios_inq_var_byid(self.group.gp, id)

    cpdef close(self):
        assert self.vp != NULL, 'Not an open variable'
        adios_free_varinfo(self.vp)
        self.group = None
        self.vp = NULL
    
    cpdef int grpid(self):
        assert self.vp != NULL, 'Not an open variable'
        return self.vp.grpid

    cpdef int varid(self):
        assert self.vp != NULL, 'Not an open variable'
        return self.vp.varid

    cpdef ADIOS_DATATYPES type(self):
        assert self.vp != NULL, 'Not an open variable'
        return self.vp.type

    cpdef int ndim(self):
        assert self.vp != NULL, 'Not an open variable'
        return self.vp.ndim

    cpdef dims(self):
        assert self.vp != NULL, 'Not an open variable'
        return [ self.vp.dims[i] for i in range(self.vp.ndim)]        

    cpdef int timedim(self):
        assert self.vp != NULL, 'Not an open variable'
        return self.vp.timedim

    cpdef int characteristics_count(self):
        assert self.vp != NULL, 'Not an open variable'
        return self.vp.characteristics_count

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
            assert npshape.ndim == npoffset.ndim, 'Offset dimension mismatch'
        
        cdef np.ndarray npcount
        if len(count) == 0:
            npcount = npshape - npoffset
        else:
            npcount = np.array(count, dtype=np.int64)
            npcount = npcount - npoffset

        assert npshape.ndim == npcount.ndim, 'Shape dimension mismatch.'
        assert (npshape - npoffset >= npcount).all(), 'Count is larger than shape.'
            
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

    cpdef printself(self):
        assert self.vp != NULL, 'Not an open file'
        print '=== AdiosVariable ==='
        print '%15s : %lu' % ('vp', <unsigned long> self.vp)
        printAdiosVariable(self.vp)
