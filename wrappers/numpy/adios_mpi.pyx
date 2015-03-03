"""
 ADIOS is freely available under the terms of the BSD license described
 in the COPYING file in the top level directory of this source distribution.

 Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
"""
"""
 This is a cython file. To generate a CPP file, use the following command:
 $ cython --cplus adios.pyx
"""

import numpy as np
cimport numpy as np

import mpi4py.MPI as MPI 
cimport mpi4py.MPI as MPI

import cython
cimport cython

## ====================
## ADIOS Exported Functions
## ====================

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
                         MPI_Comm comm)
    
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

cdef extern from "adios_selection.h":
    ctypedef enum ADIOS_SELECTION_TYPE:
        ADIOS_SELECTION_BOUNDINGBOX
        ADIOS_SELECTION_POINTS
        ADIOS_SELECTION_WRITEBLOCK
        ADIOS_SELECTION_AUTO

    ctypedef struct ADIOS_SELECTION_BOUNDINGBOX_STRUCT:
        int       ndim
        uint64_t *start
        uint64_t *count

    ctypedef struct ADIOS_SELECTION_POINTS_STRUCT:
        pass

    ctypedef struct ADIOS_SELECTION_WRITEBLOCK_STRUCT:
        pass
    
    ctypedef struct ADIOS_SELECTION_AUTO_STRUCT:
        pass

    cdef union ADIOS_SELECTION_UNION:
        ADIOS_SELECTION_BOUNDINGBOX_STRUCT bb
        ADIOS_SELECTION_POINTS_STRUCT points
        ADIOS_SELECTION_WRITEBLOCK_STRUCT block
        ADIOS_SELECTION_AUTO_STRUCT autosel

    ctypedef struct ADIOS_SELECTION:
        ADIOS_SELECTION_TYPE    type
        ADIOS_SELECTION_UNION   u

    cdef ADIOS_SELECTION * adios_selection_boundingbox (int ndim,
                                                        const uint64_t *start,
                                                        const uint64_t *count)

cdef extern from "adios_read.h":
    ctypedef enum ADIOS_READ_METHOD:
        pass

    ctypedef enum ADIOS_LOCKMODE:
        ADIOS_LOCKMODE_NONE
        ADIOS_LOCKMODE_CURRENT
        ADIOS_LOCKMODE_ALL

    ctypedef struct ADIOS_FILE:
        uint64_t fh               
        int      nvars            
        char     ** var_namelist  
        int      nattrs           
        char     ** attr_namelist 
        int      nmeshes          
        char     ** mesh_namelist 
        int      current_step     
        int      last_step        
        char     *path            
        int      endianness       
        int      version          
        uint64_t file_size

    ctypedef struct ADIOS_VARINFO:
        int        varid
        ADIOS_DATATYPES type  
        int        ndim       
        uint64_t * dims       
        int        nsteps     
        void     * value      
        int      * nblocks    
        int        sum_nblocks

    cdef int adios_read_init_method (ADIOS_READ_METHOD method, 
                                     MPI_Comm comm, 
                                     char * parameters)
    cdef int adios_read_finalize_method(ADIOS_READ_METHOD method)
    cdef ADIOS_FILE * adios_read_open (const char * fname, 
                                       ADIOS_READ_METHOD method, 
                                       MPI_Comm comm, 
                                       ADIOS_LOCKMODE lock_mode,
                                       float timeout_sec)
    cdef ADIOS_FILE * adios_read_open_file (const char * fname, 
                                            ADIOS_READ_METHOD method, 
                                            MPI_Comm comm)
    cdef int adios_read_close (ADIOS_FILE *fp)
    cdef int adios_advance_step (ADIOS_FILE *fp, int last, float timeout_sec)
    cdef void adios_release_step (ADIOS_FILE *fp)
    cdef ADIOS_VARINFO * adios_inq_var (ADIOS_FILE *fp, const char * varname)
    cdef ADIOS_VARINFO * adios_inq_var_byid (ADIOS_FILE *fp, int varid)
    cdef void adios_free_varinfo (ADIOS_VARINFO *cp)
    cdef int adios_schedule_read (const ADIOS_FILE * fp,
                                  const ADIOS_SELECTION * sel,
                                  const char * varname,
                                  int from_steps,
                                  int nsteps,
                                  void * data)
    cdef int adios_schedule_read_byid (const ADIOS_FILE * fp, 
                                       const ADIOS_SELECTION * sel,
                                       int varid,
                                       int from_steps,
                                       int nsteps,
                                       void * data)
    cdef int adios_perform_reads (const ADIOS_FILE *fp, int blocking)

    cdef int adios_get_attr (ADIOS_FILE *fp,
                             const char * attrname,
                             ADIOS_DATATYPES  * type,
                             int * size,
                             void ** data)
    
    cdef char * adios_type_to_string (ADIOS_DATATYPES type)

## ====================
## ADIOS Enum (public)
## ====================

class DATATYPE:
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

class FLAG:
    UNKNOWN = 0
    YES = 1
    NO = 2
    
class BUFFER_ALLOC_WHEN:
    UNKNOWN = 0
    NOW = 1
    LATER = 2

class READ_METHOD:
    BP            = 0 
    BP_AGGREGATE  = 1
    DATASPACES    = 3
    DIMES         = 4
    FLEXPATH      = 5
    ICEE          = 6


cpdef __parse_index(index, ndim):
    # Fix index, handling ellipsis and incomplete slices.
    if not isinstance(index, tuple):
        index = (index,)

    fixed = []
    length = len(index)
    
    for slice_ in index:
        if slice_ is Ellipsis:
            fixed.extend([slice(None)] * (ndim-length-len(fixed)+1))
        elif isinstance(slice_, (int, long)):
            fixed.append(slice(slice_, slice_+1, None))
        else:
            fixed.append(slice_)
        length -= 1

    index = tuple(fixed)
    if len(index) < ndim:
        index += (slice(None),) * (ndim-len(index))

    return index

## ====================
## ADIOS Write API
## ====================

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

cpdef int write (int64_t fd_p, char * name, val, dtype=None):
    cdef np.ndarray val_
    if isinstance(val, (np.ndarray)):
        if val.flags.contiguous:
            val_ = val
        else:
            val_ = np.array(val, copy=True)
    else:
        val_ = np.array(val, dtype=dtype)

    return adios_write (fd_p, name, <void *> val_.data)

cpdef int write_int (int64_t fd_p, char * name, int val):
    return adios_write (fd_p, name, &val)

cpdef int write_long (int64_t fd_p, char * name, long val):
    return adios_write (fd_p, name, &val)

cpdef int write_float (int64_t fd_p, char * name, float val):
    return adios_write (fd_p, name, &val)

cpdef int write_double (int64_t fd_p, char * name, double val):
    return adios_write (fd_p, name, &val)


cpdef int read(int64_t fd_p, char * name, np.ndarray val):
    assert val.flags.contiguous, 'Only contiguous arrays are supported.'
    print "Reading ... ", val.itemsize * val.size, "(bytes)"
    return adios_read(fd_p, name, <void *> val.data, val.itemsize * val.size)

cpdef int close(int64_t fd_p):
    return adios_close(fd_p)

cpdef int finalize(int mype = 0):
    return adios_finalize(mype)

## ====================
## ADIOS No-XML API
## ====================
cpdef int init_noxml(MPI.Comm comm = MPI.COMM_WORLD):
    return adios_init_noxml(comm.ob_mpi)

cpdef int allocate_buffer(int when,
                          uint64_t buffer_size):
    return adios_allocate_buffer(<ADIOS_BUFFER_ALLOC_WHEN> when,
                                 buffer_size)

cpdef int64_t declare_group(char * name,
                            char * time_index = "",
                            int stats = 1):
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
                     char * dimensions = "",
                     char * global_dimensions = "",
                     char * local_offsets = ""):
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
                         char * parameters = "",
                         char * base_path = ""):
    return adios_select_method (group,
                                method,
                                parameters,
                                base_path)


## ====================
## ADIOS Read API (V2)
## ====================

cpdef np.dtype adios2npdtype(ADIOS_DATATYPES t, int strlen = 1):
    """ strlen apply only to string type """
    cdef np.dtype ntype = None
    if t == adios_byte:
        ntype = np.dtype(np.int8)
    elif t == adios_short:
        ntype = np.dtype(np.int16)
    elif t == adios_integer:
        ntype = np.dtype(np.int32)
    elif t == adios_long:
        ntype = np.dtype(np.int64)
    elif t == adios_unsigned_byte:
        ntype = np.dtype(np.uint8)
    elif t == adios_unsigned_short:
        ntype = np.dtype(np.uint16)
    elif t == adios_unsigned_integer:
        ntype = np.dtype(np.uint32)
    elif t == adios_unsigned_long:
        ntype = np.dtype(np.uint64)
    elif t == adios_real:
        ntype = np.dtype(np.float32)
    elif t == adios_double:
        ntype = np.dtype(np.float64)
    elif t == adios_long_double:
        ntype = np.dtype(np.float128)
    elif t == adios_complex:
        ntype = np.dtype(np.complex64)
    elif t == adios_double_complex:
        ntype = np.dtype(np.complex128)
    elif t == adios_string:
        ntype = np.dtype((np.str_, strlen))
    else:
        ntype = None

    return ntype

cdef printfile(ADIOS_FILE * f):
    print '%15s : %lu' % ('fh', f.fh)
    print '%15s : %d' % ('nvars', f.nvars)
    print '%15s : %s' % ('var_namelist', [f.var_namelist[i] for i in range(f.nvars)])
    print '%15s : %d' % ('nattrs', f.nattrs)
    print '%15s : %s' % ('attr_namelist', [f.attr_namelist[i] for i in range(f.nattrs)])
    print '%15s : %d' % ('current_step', f.current_step)       
    print '%15s : %d' % ('last_step', f.last_step)       
    print '%15s : %s' % ('path', f.path)
    print '%15s : %d' % ('endianness', f.endianness)       
    print '%15s : %d' % ('version', f.version)       
    print '%15s : %lu' % ('file_size', f.file_size)

cdef printvar(ADIOS_VARINFO * v):
    print '%15s : %d' % ('varid', v.varid)
    print '%15s : %s' % ('type', adios2npdtype(v.type))
    print '%15s : %d' % ('ndim', v.ndim)
    print '%15s : %s' % ('dims', [v.dims[i] for i in range(v.ndim)])
    print '%15s : %d' % ('nsteps', v.nsteps)

cdef ADIOS_READ_METHOD str2adiosreadmethod(bytes name):
    if (name == "BP"):
        method = READ_METHOD.BP
    elif (name == "BP_AGGREGATE"):
        method = READ_METHOD.BP_AGGREGATE
    elif (name == "DATASPACES"):
        method = READ_METHOD.DATASPACES
    elif (name == "DIMES"):
        method = READ_METHOD.DIMES
    elif (name == "FLEXPATH"):
        method = READ_METHOD.FLEXPATH
    elif (name == "ICEE"):
        method = READ_METHOD.ICEE
    else:
        print '[WARN] Invalid read method name:', name, '. Use default BP method'
        method = READ_METHOD.BP
        
    return method

cpdef np2adiostype(type nptype):
    """ Ignored: int_, intc, intp """

    cdef atype = DATATYPE.unknown

    if (nptype == np.bool_):
        atype = DATATYPE.integer
    elif (nptype == np.int8):
        atype = DATATYPE.byte
    elif (nptype == np.int16):
        atype = DATATYPE.short
    elif (nptype == np.int32):
        atype = DATATYPE.integer
    elif (nptype == np.int64):
        atype = DATATYPE.long
    elif (nptype == np.uint8):
        atype = DATATYPE.unsigned_byte
    elif (nptype == np.uint16):
        atype = DATATYPE.unsigned_short
    elif (nptype == np.uint32):
        atype = DATATYPE.unsigned_integer
    elif (nptype == np.uint64):
        atype = DATATYPE.unsigned_long
    elif (nptype == np.float_):
        atype = DATATYPE.double
    elif (nptype == np.float16):
        atype = DATATYPE.real
    elif (nptype == np.float32):
        atype = DATATYPE.real
    elif (nptype == np.float64):
        atype = DATATYPE.double
    elif (nptype == np.complex_):
        atype = DATATYPE.double_complex
    elif (nptype == np.complex64):
        atype = DATATYPE.complex
    elif (nptype == np.complex128):
        atype = DATATYPE.double_complex
    elif (nptype == np.str_):
        atype = DATATYPE.string

    return atype

cpdef str adiostype2string (ADIOS_DATATYPES type):
    return str(adios_type_to_string(<ADIOS_DATATYPES> type))

## ====================
## ADIOS Class Definitions for Read
## ====================

""" Call adios_read_init_method """
cpdef int read_init(char * method_name = "BP",
                    MPI.Comm comm = MPI.COMM_WORLD,
                    char * parameters = ""):
    cdef method = str2adiosreadmethod(method_name)
    return adios_read_init_method (method, comm.ob_mpi, parameters)


""" Call adios_read_finalize_method """
cpdef int read_finalize(char * method_name = "BP"):
    cdef method = str2adiosreadmethod(method_name)
    return adios_read_finalize_method (method)

""" Python class for ADIOS_FILE structure """
cdef class file:
    """ Private Memeber """
    cpdef ADIOS_FILE * fp

    """ Public Memeber """
    cpdef public bytes name
    cpdef public int nvars
    cpdef public int nattrs
    cpdef public int current_step
    cpdef public int last_step
    cpdef public int endianness
    cpdef public int version
    cpdef public int file_size
    
    cpdef public dict var
    cpdef public dict attr

    cpdef public bint is_stream

    """ Initialization. Call adios_read_open and populate public members """
    def __init__(self, char * fname,
                 char * method_name = "BP",
                 MPI.Comm comm = MPI.COMM_WORLD,
                 is_stream = False,
                 ADIOS_LOCKMODE lock_mode = ADIOS_LOCKMODE_ALL,
                 float timeout_sec = 0.0):
        self.fp = NULL
        self.var = {}
        self.attr = {}
        self.is_stream = is_stream
        cdef method = str2adiosreadmethod(method_name)

        if (is_stream):
            self.fp = adios_read_open(fname, method, comm.ob_mpi,
                                      lock_mode, timeout_sec)
        else:
            self.fp = adios_read_open_file(fname, method, comm.ob_mpi)
            
        assert self.fp != NULL, 'Not an open file'

        self.name = fname.split('/')[-1]  ## basename
        self.nvars = self.fp.nvars
        self.nattrs = self.fp.nattrs
        self.current_step = self.fp.current_step
        self.last_step = self.fp.last_step
        self.endianness = self.fp.endianness
        self.version = self.fp.version
        self.file_size = self.fp.file_size
    
        for name in [self.fp.var_namelist[i] for i in range(self.nvars)]:
            self.var[name] = var(self, name)

        for name in [self.fp.attr_namelist[i] for i in range(self.nattrs)]:
            self.attr[name] = self.get_attr(name)

    def __del__(self):
            self.close()
            
    """ Call adios_read_close """
    cpdef close(self):
        assert self.fp != NULL, 'Not an open file'
        adios_read_close(self.fp)
        self.fp = NULL
        
    """ Print self """
    cpdef printself(self):
        assert self.fp != NULL, 'Not an open file'
        print '=== AdiosFile ==='
        print '%15s : %lu' % ('fp', <unsigned long> self.fp)
        printfile(self.fp)

    cpdef advance(self, int last = 0, float timeout_sec = 0.0):
        val = adios_advance_step(self.fp, last, timeout_sec)
        if (val >= 0):
            self.current_step = self.fp.current_step
            self.last_step = self.fp.last_step

            for v in self.var.values():
                v.advance()
                
        return val

    cpdef get_attr(self, char * attr_name):
        cdef ADIOS_DATATYPES type
        cdef int size
        cdef int64_t p
        adios_get_attr(self.fp, attr_name, &type, &size, <void **> &p)
        ##print '=== Attribute ==='
        ##print '%15s : %d (%s)' % ('type', type, adios2nptype(type))
        ##print '%15s : %d' % ('size', size)
        ##print '%15s : %lu' % ('data', <unsigned long> p)

        cdef DTYPE = adios2npdtype(type, size)
        cdef np.ndarray var = np.zeros(1, dtype=DTYPE)
        adios_get_attr(self.fp, attr_name, &type, &size, <void **> &var.data)
        
        return np.asscalar(var)
        
    def __getitem__(self, index):
        if not isinstance(index, tuple):
            index = (index,)

        if len(index) > 1:
            raise KeyError(index)
        
        for key_ in index:
            if not isinstance(key_, str):
                raise TypeError("Unhashable type")

            if key_ in self.var.keys():
                return self.var.get(key_)
            elif key_ in self.attr.keys():
                return self.attr.get(key_)
            else:
                raise KeyError(key_)
        
    def __repr__(self):
        return ("AdiosFile (path=%r, nvars=%r, var=%r, nattrs=%r, attr=%r, "
                "current_step=%r, last_step=%r, file_size=%r)") % \
                (self.fp.path,
                 self.nvars,
                 self.var.keys(),
                 self.nattrs,
                 self.attr.keys(),
                 self.current_step,
                 self.last_step,
                 self.file_size)

""" Python class for ADIOS_VARINFO structure """
cdef class var:
    """ Private Memeber """
    cdef file file
    cdef ADIOS_VARINFO * vp

    """ Public Memeber """
    cpdef public bytes name
    cpdef public int varid
    cpdef public np.dtype type
    cpdef public int ndim
    cpdef public tuple dims
    cpdef public int nsteps

    """ Initialization. Call adios_inq_var and populate public members """
    def __init__(self, file file, char * name):
        self.file = file
        self.vp = NULL

        assert self.file.fp != NULL, 'Not an open file'
        self.vp = adios_inq_var(self.file.fp, name)
        assert self.vp != NULL, 'Not a valid var'

        self.name = name
        self.varid = self.vp.varid
        self.ndim = self.vp.ndim                 
        self.dims = tuple([self.vp.dims[i] for i in range(self.vp.ndim)])
        self.nsteps = self.vp.nsteps

        if self.vp.type == DATATYPE.string:
            self.type = adios2npdtype(self.vp.type, len(<char*> self.vp.value))
        else:
            self.type = adios2npdtype(self.vp.type)

    def __del__(self):
        self.close()

    """ Call adios_free_varinfo """
    cpdef close(self):
        assert self.vp != NULL, 'Not an open var'
        adios_free_varinfo(self.vp)
        self.vp = NULL

    cpdef advance(self):
        self.vp = adios_inq_var(self.file.fp, self.name)
        assert self.vp != NULL, 'Not a valid var'
        self.nsteps = self.vp.nsteps

    """ Call adios_schedule_read and adios_perform_reads """
    cpdef read(self, tuple offset = (), tuple count = (), from_steps = None, nsteps = None, fill = 0):
        if from_steps is None:
            from_steps = 0 ##self.file.current_step

        if nsteps is None:
            nsteps = self.file.last_step - self.file.current_step + 1

        assert self.type is not None, 'Data type is not supported yet'
        if (self.nsteps > 0) and (from_steps + nsteps > self.nsteps):
            raise IndexError('Step index is out of range: from_steps=%r, nsteps=%r' % (from_steps, nsteps))
        
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

        if npshape.ndim != npoffset.ndim:
            raise IndexError('Offset dimension mismatch (offset dim: %r)' % (npoffset.ndim))

        if npshape.ndim != npcount.ndim:
            raise IndexError('Count dimension mismatch (count dim: %r)' % (npcount.ndim))

        if (npshape < npcount + npoffset).any():
            raise IndexError('Requested is larger than the shape.')

        shape = list(npcount)
        if (nsteps > 1):
            shape.insert(0, nsteps)
        cdef np.ndarray var = np.full(shape, fill, dtype=self.type)

        cdef ADIOS_SELECTION * sel
        sel = adios_selection_boundingbox (self.vp.ndim, <uint64_t *> npoffset.data, <uint64_t *> npcount.data)

        ##print 'npoffset', npoffset
        ##print 'npcount', npcount

        adios_schedule_read_byid (self.file.fp, sel, self.vp.varid, from_steps, nsteps, <void *> var.data)
        adios_perform_reads(self.file.fp, 1)

        ## Try not to return as scalar to be consistent
        ##if (var.ndim == 0):
        ##    return np.asscalar(var)
        ##else:
        ##    return var
        return var

    """ Print self """
    cpdef printself(self):
        assert self.vp != NULL, 'Not an open variable'
        print '=== AdiosVariable ==='
        print '%15s : %lu' % ('vp', <unsigned long> self.vp)
        print '%15s : %lu' % ('fp', <unsigned long> self.file.fp)
        printvar(self.vp)
        
    def __repr__(self):
        return "AdiosVar (varid=%r, type=%r, ndim=%r, dims=%r, nsteps=%r)" % \
               (self.varid,
                self.type,
                self.ndim,
                self.dims,
                self.nsteps)

    def __getitem__(self, index):
        ndim_ = self.ndim
        if (self.nsteps) > 1: ndim_ += 1

        index_ = __parse_index(index, ndim_)

        if (ndim_ > 0) and (len(index_) > ndim_):
            raise IndexError("Too many indices for data")

        if (ndim_ == 0) and (len(index_) > 1):
            raise IndexError("Too many indices for data")
        
        for slice_ in index_:
            if isinstance(slice_.step, (int, long)) and (slice_.step != 1):
                raise IndexError("Step size (%d) is not supported." % (slice_.step))
            if isinstance(slice_, str):
                raise IndexError("Name index (%r) is not supported." % (slice_))
        
        if (self.nsteps) > 1:
            dims_ = list(self.dims)
            dims_.insert(0, self.nsteps)
            indices = tuple(x[0].indices(x[1]) for x in zip(index_, dims_))
            z = zip(*indices)

            from_steps_ = z[0][0]
            nsteps_ = (z[1][0] - z[0][0]-1)%self.nsteps+1
            offset_ = z[0][1:]
            count_ = tuple((np.subtract(z[1][1:], z[0][1:])-1)%dims_[1:]+1)
        else:
            indices = tuple(x[0].indices(x[1]) for x in zip(index_, self.dims))
            z = zip(*indices)

            if len(z) == 0:
                from_steps_ = 0
                nsteps_ = self.nsteps
                offset_ = ()
                count_ = ()
            else:
                from_steps_ = 0
                nsteps_ = self.nsteps
                offset_ = z[0]
                count_ = tuple((np.subtract(z[1], z[0])-1)%self.dims+1)

        ##print "from_steps", from_steps_
        ##print "nsteps", nsteps_
        ##print "offset", offset_
        ##print "count", count_
        
        return self.read(offset=offset_,
                         count=count_,
                         from_steps=from_steps_,
                         nsteps=nsteps_)
        
## ====================
## ADIOS Global functions
## ====================

def readvar(fname, varname):
    f = file(fname, comm=MPI.COMM_SELF)
    if not f.var.has_key(varname):
        print "No valid variable"
        return

    v = f.var[varname]
    return v.read(from_steps=0, nsteps=v.nsteps)

def bpls(fname):
    f = file(fname, comm=MPI.COMM_SELF)
    return {'nvars': f.nvars,
            'nattrs': f.nattrs,
            'vars': tuple([ k for k in f.var.iterkeys() ]),
            'attrs': tuple([ k for k in f.attr.iterkeys() ]),
            'time_steps': (f.current_step, f.last_step),
            'file_size': f.file_size}
