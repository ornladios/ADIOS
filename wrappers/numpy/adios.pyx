# -*- coding: utf-8 -*-
"""ADIOS: ADIOS python module

.. moduleauthor:: Jong Choi <choij@ornl.gov>
"""

#cdef extern from "mpi-compat.h": pass

import numpy as np
cimport numpy as np

#import mpi4py.MPI_as MPI
#cimport mpi4py.MPI_as MPI

import cython
cimport cython

from libc.stdlib cimport malloc, free
from cpython.string cimport PyString_AsString

import os
from ._hl import selections as sel

cdef char ** to_cstring_array(list_str):
    cdef char **ret = <char **>malloc(len(list_str) * sizeof(char *))
    for i in xrange(len(list_str)):
        ret[i] = PyString_AsString(list_str[i])
    return ret

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
        adios_string_array

    ctypedef enum ADIOS_BUFFER_ALLOC_WHEN:
        ADIOS_BUFFER_ALLOC_UNKNOWN
        ADIOS_BUFFER_ALLOC_NOW
        ADIOS_BUFFER_ALLOC_LATER

    ctypedef enum ADIOS_FLAG:
        pass

cdef extern from "adios.h":
    ctypedef int MPI_Comm
    int MPI_COMM_WORLD
    int MPI_COMM_SELF

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

    cdef int adios_define_attribute_byvalue (int64_t group,
                                             char * name,
                                             char * path,
                                             ADIOS_DATATYPES type,
                                             int nelems,
                                             void * values)

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

    cdef ADIOS_SELECTION * adios_selection_boundingbox (uint64_t ndim,
                                                        const uint64_t *start,
                                                        const uint64_t *count)

    cdef ADIOS_SELECTION * adios_selection_points (uint64_t ndim,
                                                   uint64_t npoints,
                                                   const uint64_t *points)

    cdef void adios_selection_delete (ADIOS_SELECTION * sel)

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
    string_array = 12

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

class LOCKMODE:
    NONE = 0
    CURRENT = 1
    ALL =2

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
        elif isinstance(slice_, (float)):
            fixed.append(slice(int(slice_), int(slice_)+1, None))
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

cdef MPI_Comm init_comm
cdef MPI_Comm read_init_comm

cpdef init(char * config, MPI_Comm comm = MPI_COMM_WORLD):
    global init_comm
    init_comm = comm
    return adios_init(config, init_comm)

cpdef int64_t open(char * group_name,
                   char * name,
                   char * mode,
                   MPI_Comm comm = MPI_COMM_WORLD):
    cdef int64_t fd
    cdef int result
    result = adios_open(&fd, group_name, name, mode, comm)
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

    cdef void * ptr
    if (val_.dtype.char == 'S'):
        ptr = <void *> PyString_AsString(val_.data)
    else:
        ptr = <void *> val_.data

    return adios_write (fd_p, name, ptr)

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

cpdef int init_noxml(MPI_Comm comm = MPI_COMM_WORLD):
    global init_comm
    init_comm = comm
    return adios_init_noxml(init_comm)

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
                     int atype,
                     char * dimensions = "",
                     char * global_dimensions = "",
                     char * local_offsets = ""):
    return adios_define_var(group_id,
                            name, path,
                            <ADIOS_DATATYPES> atype,
                            dimensions,
                            global_dimensions,
                            local_offsets)

cpdef int define_attribute (int64_t group,
                            char * name,
                            char * path,
                            int atype,
                            char * value,
                            char * var):
    return adios_define_attribute (group,
                                   name,
                                   path,
                                   <ADIOS_DATATYPES> atype,
                                   value,
                                   var)

cpdef int define_attribute_byvalue (int64_t group,
                                    char * name,
                                    char * path,
                                    val):
    cdef np.ndarray val_
    if isinstance(val, (np.ndarray)):
        if val.flags.contiguous:
            val_ = val
        else:
            val_ = np.array(val, copy=True)
    else:
        val_ = np.array(val)

    atype = np2adiostype(val_.dtype)

    cdef char * pt1
    cdef char ** pt2
    if (val_.dtype.char == 'S'):
        if (val_.size == 1):
            pt1 = PyString_AsString(val)
            adios_define_attribute_byvalue (group,
                                            name,
                                            path,
                                            DATATYPE.string,
                                            1,
                                            <void *> pt1)
        else:
            pt2 = to_cstring_array(val)
            adios_define_attribute_byvalue (group,
                                            name,
                                            path,
                                            DATATYPE.string_array,
                                            len(val),
                                            <void *> pt2)
            free(pt2)
    else:
        adios_define_attribute_byvalue (group,
                                        name,
                                        path,
                                        <ADIOS_DATATYPES> atype,
                                        val_.size,
                                        <void *> val_.data)

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

cpdef np2adiostype(np.dtype nptype):
    """ Convert Numpy.dtype to Adios Datatype
    """

    cdef atype = DATATYPE.unknown

    if (nptype == np.bool_):
        atype = DATATYPE.byte
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
    elif (nptype.char == 'S'):
        atype = DATATYPE.string
    else:
        atype = DATATYPE.unknown

    return atype

cpdef str adiostype2string (ADIOS_DATATYPES type):
    return str(adios_type_to_string(<ADIOS_DATATYPES> type))

## ====================
## ADIOS Class Definitions for Read
## ====================

""" Call adios_read_init_method """
cpdef int read_init(char * method_name = "BP",
                    MPI_Comm comm = MPI_COMM_WORLD,
                    char * parameters = ""):
    global read_init_comm
    read_init_comm = comm
    cdef method = str2adiosreadmethod(method_name)
    return adios_read_init_method (method, read_init_comm, parameters)


""" Call adios_read_finalize_method """
cpdef int read_finalize(char * method_name = "BP"):
    cdef method = str2adiosreadmethod(method_name)
    return adios_read_finalize_method (method)

## Python class for ADIOS_FILE structure
cdef class file(object):
    """
    file class for Adios file read and write.

    Args:
        fname (str): filename.
        method_name (str, optional): Adios read method (default: 'BP').
        comm (MPI_Comm, optional): MPI_comm for parallel read/write (default: MPI_COMM_WORLD).
        is_stream (bool, optional): Set True if use stream reader (default: False).
        lock_mode (int, optional): ADIOS_LOCKMODE for stream reader (default: ADIOS_LOCKMODE_ALL).
        timeout_sec (float, optional): Timeout seconds for stream reader (default: 0.0).

    Example:

    >>> import adios as ad
    >>> f = ad.file('adiosfile.bp')

    """

    cpdef ADIOS_FILE * fp
    cpdef bytes name
    cpdef int nvars
    cpdef int nattrs
    cpdef int current_step
    cpdef int last_step
    cpdef int endianness
    cpdef int version
    cpdef int file_size
    cpdef bint is_stream

    ## Public Memeber
    cpdef public dict vars
    cpdef public dict attrs
    cpdef public var
    cpdef public attr

    property name:
        """ The filename (or stream name) associated with. """
        def __get__(self):
            return self.name

    property nvars:
        """ The number of variables. """
        def __get__(self):
            return self.nvars

    property nattrs:
        """ The number of attributes. """
        def __get__(self):
            return self.nattrs

    property current_step:
        """ The current timestep index. """
        def __get__(self):
            return self.current_step

    property last_step:
        """ The last timestep index. """
        def __get__(self):
            return self.last_step

    property endianness:
        """ The endianness of the stored data. """
        def __get__(self):
            return self.endianness

    property version:
        """ The version of Adios. """
        def __get__(self):
            return self.version

    property file_sizec:
        """ The size of Adios file. """
        def __get__(self):
            return self.file_size

    property is_stream:
        """ Indicating reader type; file reader or stream reader """
        def __get__(self):
            return self.is_stream

    def __init__(self, char * fname,
                 char * method_name = 'BP',
                 MPI_Comm comm = MPI_COMM_WORLD,
                 is_stream = False,
                 ADIOS_LOCKMODE lock_mode = ADIOS_LOCKMODE_ALL,
                 float timeout_sec = 0.0):
        self.fp = NULL
        self.vars = {}
        self.attrs = {}
        self.is_stream = is_stream
        cdef method = str2adiosreadmethod(method_name)

        if (is_stream):
            self.fp = adios_read_open(fname, method, comm,
                                      lock_mode, timeout_sec)
        else:
            self.fp = adios_read_open_file(fname, method, comm)

        assert self.fp != NULL, 'Not an open file'

        self.name = fname.split('/')[-1]  ## basename
        self.nvars = self.fp.nvars
        self.nattrs = self.fp.nattrs
        self.current_step = self.fp.current_step
        self.last_step = self.fp.last_step
        self.endianness = self.fp.endianness
        self.version = self.fp.version
        self.file_size = self.fp.file_size

        for name in [self.fp.attr_namelist[i] for i in range(self.nattrs)]:
            self.attrs[name] = attr(self, name)

        for name in [self.fp.var_namelist[i] for i in range(self.nvars)]:
            self.vars[name] = var(self, name)

        self.var = self.vars
        self.attr = self.attrs

    def __del__(self):
        """ Close file on destruction. """
        self.close()

    cpdef close(self):
        """ Close the open file. """
        assert self.fp != NULL, 'Not an open file'
        adios_read_close(self.fp)
        self.fp = NULL

    cpdef printself(self):
        """ Print native ADIOS_FILE structure. """
        assert self.fp != NULL, 'Not an open file'
        print '=== AdiosFile ==='
        print '%15s : %lu' % ('fp', <unsigned long> self.fp)
        printfile(self.fp)

    cpdef release_step(self):
        """ Release the current step lock and let the writer code to progress """
        adios_release_step(self.fp)

    cpdef advance(self, int last = 0, float timeout_sec = 0.0):
        """
        Advance a timestep for stream reader.

        Args:
            last (int, optional): last timestep index (default: 0).
            timeout_sec (float, optional): timeout seconds (default: 0.0).

        Returns:
            int: 0 if successful, non-zero otherwise.
        """
        val = adios_advance_step(self.fp, last, timeout_sec)
        if (val >= 0):
            self.current_step = self.fp.current_step
            self.last_step = self.fp.last_step

            for v in self.vars.values():
                v.advance()

        return val

    def __getitem__(self, varname):
        """
        Return Adios variable.

        Args:
            varname (str): variable name.

        Raises:
            KeyError: If no varname exists.

        """
        if not isinstance(varname, tuple):
            varname = (varname,)

        if len(varname) > 1:
            raise KeyError(varname)

        for key_ in varname:
            if not isinstance(key_, str):
                raise TypeError("Unhashable type")

            if key_ in self.vars.keys():
                return self.vars.get(key_)
            elif key_ in self.attrs.keys():
                return self.attrs.get(key_)

            #TODO: return group (self, groupname)
            for name in self.vars.keys():
                if key_ == os.path.dirname(name):
                    return group(self, key_)

            for name in self.attrs.keys():
                if key_ == os.path.dirname(name):
                    return group(self, key_)

        raise KeyError(key_)


    def __repr__(self):
        """ Return string representation. """
        return ("AdiosFile (path=%r, nvars=%r, vars=%r, nattrs=%r, attrs=%r, "
                "current_step=%r, last_step=%r, file_size=%r)") % \
                (self.fp.path,
                 self.nvars,
                 self.vars.keys(),
                 self.nattrs,
                 self.attrs.keys(),
                 self.current_step,
                 self.last_step,
                 self.file_size)

cdef class var(object):
    """
    Adios variable class.

    Unlike attributes whose values are populated on initialization,
    variable's values will be returned by explicitly calling read() or
    array access interface ([]).

    Args:
        file (file): Associated file class
        name (str): variable name

    Note:
        Users do not need to create this class manually.
    """

    cdef file file
    cdef ADIOS_VARINFO * vp

    cpdef bytes name
    cpdef int varid
    cpdef np.dtype dtype
    cpdef int ndim
    cpdef tuple dims
    cpdef int nsteps
    cpdef dict attrs

    property name:
        """ The variable name. """
        def __get__(self):
            return self.name

    property varid:
        """ Internal variable id. """
        def __get__(self):
            return self.varid

    property dtype:
        """ Variable type as in numpy.dtype. """
        def __get__(self):
            return self.dtype

    property ndim:
        """ The number of dimensions of the variable. """
        def __get__(self):
            return self.ndim

    property dims:
        """ The shape of the variable. """
        def __get__(self):
            return self.dims

    property nsteps:
        """ The number of time steps of the variable. """
        def __get__(self):
            return self.nsteps

    property attrs:
        def __get__(self):
            return self.attrs

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
            self.dtype = adios2npdtype(self.vp.type, len(<char*> self.vp.value))
        else:
            self.dtype = adios2npdtype(self.vp.type)

        self.attrs = {}
        for name in self.file.attrs.keys():
            if name.startswith(self.name + '/'):
                self.attrs[name.replace(self.name + '/', '')] = self.file.attrs[name]
            if name.startswith('/' + self.name + '/'):
                self.attrs[name.replace('/' + self.name + '/', '')] = self.file.attrs[name]

    def __del__(self):
        self.close()

    cpdef close(self):
        """ Close and free variable information """
        assert self.vp != NULL, 'Not an open var'
        adios_free_varinfo(self.vp)
        self.vp = NULL

    cpdef advance(self):
        """ Update variable information after the stream advanced """
        self.vp = adios_inq_var(self.file.fp, self.name)
        assert self.vp != NULL, 'Not a valid var'
        self.nsteps = self.vp.nsteps

    cpdef read_points(self, tuple points = (), from_steps = None, nsteps = None):
        """ Perform points read.

        Read data from an ADIOS BP file based on the given list of point index.

        Args:
            points (tuple of int, optional): points index defined by ((o1,o2,...,oN),...) (default: ())
            from_steps (int, optional): starting step index (default: None)
            nsteps (int, optional): number of time dimensions (default: None)

        Returns:
            NumPy 1-D ndarray

        Raises:
            IndexError: If dimension is mismatched or out of the boundary.
        """
        if from_steps is None:
            from_steps = 0 ##self.file.current_step

        if nsteps is None:
            nsteps = self.file.last_step - from_steps + 1

        assert self.dtype is not None, 'Data type is not supported yet'

        if (self.nsteps > 0) and (from_steps + nsteps > self.nsteps):
            raise IndexError('Step index is out of range: from_steps=%r, nsteps=%r' % (from_steps, nsteps))

        if not isinstance(points, tuple):
            points = (points,)

        if len(points) > 1:
            plen = len(points[0])
            if not all([len(x) == plen for x in points]):
                raise IndexError('All points must have the same length %r' % (points,))

        cpdef uint64_t ndim = self.ndim
        cpdef uint64_t npoints = len(points)
        ##print 'ndim, npoints = %r, %r' % (ndim, npoints)

        cdef np.ndarray nppoints = np.array(points, dtype=np.int64, order='C')
        ##print 'nppoints.ndim = %r' % (nppoints.ndim)
        ##print 'nppoints.shape = (%r, %r)' % (nppoints.shape[0], nppoints.shape[1])

        cdef np.ndarray var = np.zeros((npoints * nsteps,), dtype=self.dtype)
        ##print 'nppoints.ndim = %r' % (nppoints.ndim)
        ##print 'nppoints.shape = (%r, %r)' % (nppoints.shape[0], nppoints.shape[1])

        cdef ADIOS_SELECTION * sel
        sel = adios_selection_points (ndim, npoints, <uint64_t *> nppoints.data)

        adios_schedule_read_byid (self.file.fp, sel, self.vp.varid, from_steps, nsteps, <void *> var.data)
        adios_perform_reads(self.file.fp, 1)
        adios_selection_delete(sel)

        return var

    cpdef read(self, tuple offset = (), tuple count = (), tuple scalar = (),
               from_steps = None, nsteps = None, fill = 0, step_scalar = True):
        """ Perform read.

        Read data from an ADIOS BP file. Subset reading is
        supported. Without any options, this will read out a whole
        data.

        Args:
            offset (tuple of int, optional): offset (default: ())
            count (tuple of int, optional): count (default: ())
            scalar (tuple of bool, optional): scalar (default: ())
            from_steps (int, optional): starting step index (default: None)
            nsteps (int, optional): number of time dimensions (default: None)
            fill (value, optional): default fill value (default: 0)
            step_scalar (bool, optional): add time dim or not (default: True)

        Returns:
            NumPy ndarray

        Raises:
            IndexError: If dimension is mismatched or out of the boundary.


        Example:

        The following command will read the full data:

        >>> var.read()

        which is equvalent to

        >>> var[]


        The following command is for subset reading:

        >>> var.read(offset=(1,2), count=(3,4))

        which will return an 3x4 array offset by (1,2) in the original
        data. With Numpy's array notation, the following command does the same job:

        >>> var[1:4, 2:6]

        Similarly, the following two commands are same:

        >>> var.read(count=(5,6))
        >>> var[:5, :6]

        """
        if from_steps is None:
            from_steps = 0 ##self.file.current_step

        if nsteps is None:
            nsteps = self.file.last_step - from_steps + 1

        assert self.dtype is not None, 'Data type is not supported yet'

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

        if len(scalar) == 0:
            scalar = tuple((False,) * len(npshape))

        if len(npshape) != len(npoffset):
            raise IndexError('Offset dimension mismatch (offset dim: %r)' % len(npoffset))

        if len(npshape) != len(npcount):
            raise IndexError('Count dimension mismatch (count dim: %r)' % len(npcount))

        if len(npshape) != len(scalar):
            raise IndexError('Scalar dimension mismatch (scalar dim: %r)' % len(scalar))

        if (npshape < npcount + npoffset).any():
            raise IndexError('Requested is larger than the shape.')

        shape = [x for x, y in zip(npcount, scalar) if not y]
        if (nsteps > 1) or (self.nsteps>1 and not step_scalar):
            shape.insert(0, nsteps)
        cdef np.ndarray var = np.zeros(shape, dtype=self.dtype)

        if len(shape) > 0:
            var[:] = fill

        cdef ADIOS_SELECTION * sel
        sel = adios_selection_boundingbox (self.vp.ndim, <uint64_t *> npoffset.data, <uint64_t *> npcount.data)

        ##print 'npoffset', npoffset
        ##print 'npcount', npcount

        adios_schedule_read_byid (self.file.fp, sel, self.vp.varid, from_steps, nsteps, <void *> var.data)
        adios_perform_reads(self.file.fp, 1)
        adios_selection_delete(sel)

        if (var.ndim == 0):
            return np.asscalar(var)
        else:
            return var

    cpdef printself(self):
        """ Print native ADIOS_VARINFO structure. """
        assert self.vp != NULL, 'Not an open variable'
        print '=== AdiosVariable ==='
        print '%15s : %lu' % ('vp', <unsigned long> self.vp)
        print '%15s : %lu' % ('fp', <unsigned long> self.file.fp)
        printvar(self.vp)

    def __repr__(self):
        return "AdiosVar (varid=%r, name=%r, dtype=%r, ndim=%r, dims=%r, nsteps=%r)" % \
               (self.varid,
                self.name,
                self.dtype,
                self.ndim,
                self.dims,
                self.nsteps)

    def __getitem__(self, args):
        shape = list(self.dims)
        if self.nsteps > 1:
            shape.insert(0, self.nsteps)
        asel = sel.select(shape, args)

        if isinstance(asel, sel.SimpleSelection):
            if (self.nsteps) > 1:
                return self.read(offset=asel.sel[0][1:],
                                 count=asel.sel[1][1:],
                                 scalar=asel.sel[3][1:],
                                 from_steps=asel.sel[0][0],
                                 nsteps=asel.sel[1][0],
                                 step_scalar=asel.sel[3][0])
            else:
                return self.read(offset=asel.sel[0],
                                 count=asel.sel[1],
                                 scalar=asel.sel[3],
                                 from_steps=0,
                                 nsteps=1)

        elif isinstance(asel, sel.FancySelection):
            shape = list(asel.sel[0][1])
            shape[asel.morder[0]] = 0
            var = np.ndarray(shape, dtype=self.dtype)
            for idx, obj in enumerate(asel.sel):
                if (self.nsteps) > 1:
                    v = self.read(offset=obj[0][1:],
                                  count=obj[1][1:],
                                  scalar=obj[3][1:],
                                  from_steps=obj[0][0],
                                  nsteps=obj[1][0],
                                  step_scalar=obj[3][0])
                else:
                    v = self.read(offset=obj[0],
                                  count=obj[1],
                                  scalar=obj[3],
                                  from_steps=0,
                                  nsteps=1)

                var = np.concatenate((var, v), axis=asel.morder[idx])

            var = np.reshape(var, asel.mshape)

            return var

        else:
            raise NotImplementedError("Not implemented yet")

cdef class attr(object):
    """
    Adios attribute class.

    Attribute values are loaded on initialization.

    Args:
        attr_name (str): attribute name

    Raises:
        KeyError: If no attribute name exists.

    Note:
        Users do not need to create this class manually.
    """
    cdef file file
    cpdef bytes name
    cpdef np.dtype dtype
    cdef np.ndarray value

    property name:
        """ The attribute name """
        def __get__(self):
            return self.name

    property dtype:
        """ The attribute type as in numpy.dtype """
        def __get__(self):
            return self.dtype

    property value:
        """ The attribute's value """
        def __get__(self):
            if (self.value.ndim == 0):
                return np.asscalar(self.value)
            else:
                return self.value

    def __init__(self, file file, char * name):
        self.file = file
        self.name = name

        cdef int64_t p
        cdef ADIOS_DATATYPES atype
        cdef int bytes
        cdef list strlist
        cdef int len

        err = adios_get_attr(self.file.fp, self.name, &atype, &bytes, <void **> &p)

        if err == 0:
            if atype == DATATYPE.string:
                bytes = bytes - 1 ## Remove the NULL terminal
            self.dtype = adios2npdtype(atype, bytes)
            if atype == DATATYPE.string_array:
                strlist = list()
                len = bytes/sizeof(p)
                for i in range(len):
                    strlist.append((<char **>p)[i])
                self.value = np.array(strlist)
                self.dtype = self.value.dtype

            elif self.dtype is None:
                print 'Warning: No support yet: %s (type=%d, bytes=%d)' % \
                      (self.name, atype, bytes)
            else:
                len = bytes/self.dtype.itemsize
                if len == 1:
                    self.value = np.array(len, dtype=self.dtype)
                else:
                    self.value = np.zeros(len, dtype=self.dtype)
                self.value.data = <char *> p
        else:
            raise KeyError(name)

    def __getitem__(self, args):
        val = self.value[args]
        if (val.ndim == 0):
            return np.asscalar(val)
        else:
            return val

    def __repr__(self):
        return "AdiosAttr (name=%r, dtype=%r, value=%r)" % \
               (self.name, self.dtype, self.value)

cdef class group(object):
    """
    Adios group class.

    Note:
        Users do not need to create this class manually.
    """
    cdef file file
    cpdef bytes name

    ## Public Memeber
    cpdef public dict vars
    cpdef public dict attrs

    def __init__(self, file file, char * name):
        self.file = file
        self.name = name.rstrip('/')

        self.vars = {}
        for name in self.file.vars.keys():
            if name.startswith(self.name + '/'):
                self.vars[name.replace(self.name + '/', '', 1)] = self.file.vars[name]
            if name.startswith('/' + self.name + '/'):
                self.vars[name.replace('/' + self.name + '/', '', 1)] = self.file.vars[name]

        self.attrs = {}
        for name in self.file.attrs.keys():
            if name.startswith(self.name + '/'):
                self.attrs[name.replace(self.name + '/', '', 1)] = self.file.attrs[name]
            if name.startswith('/' + self.name + '/'):
                self.attrs[name.replace('/' + self.name + '/', '', 1)] = self.file.attrs[name]

    def __getitem__(self, varname):
        """
        Return Adios variable or attribute.

        Args:
            varname (str): variable or attribute name.

        Raises:
            KeyError: If no name exists.

        """
        if not isinstance(varname, tuple):
            varname = (varname,)

        if len(varname) > 1:
            raise KeyError(varname)

        for key_ in varname:
            if not isinstance(key_, str):
                raise TypeError("Unhashable type")

            if key_ in self.vars.keys():
                return self.vars.get(key_)
            elif key_ in self.attrs.keys():
                return self.attrs.get(key_)

        raise KeyError(key_)

    def __repr__(self):
        """ Return string representation. """
        return ("AdiosGroup (vars=%r, attrs=%r)") % \
                (self.vars.keys(),
                 self.attrs.keys())

## Helper dict
cdef class smartdict(dict):
    cdef factory
    def __init__(self, factory):
        dict.__init__(self)
        self.factory = factory

    def __setitem__(self, key, value):
        if key in dict.keys(self):
            dict.__setitem__(self, key, value)
        else:
            self.factory(key, value)

cdef class writer(object):
    """
    writer class for Adios write.

    Args:
        fname (str): filename.
        is_noxml (bool, optional): Set True if use noxml APIs (default: True).
        comm (MPI_Comm, optional): MPI_comm for parallel read/write (default: MPI_COMM_WORLD).

    Example:

    >>> import adios as ad
    >>> f = ad.writer('adiosfile.bp')

    """

    cdef int64_t gid
    cpdef bytes fname
    cpdef bytes gname
    cpdef bytes method
    cpdef bytes method_params
    cpdef bint is_noxml
    cpdef bytes mode
    cpdef MPI_Comm comm

    cpdef dict vars
    cpdef dict attrs

    property fname:
        """ The filename to write. """
        def __get__(self):
            return self.fname

    property gname:
        """ The groupname associated with the file. """
        def __get__(self):
            return self.gname

    property is_noxml:
        """ Boolean to indicate using No-XML or not. """
        def __get__(self):
            return self.is_noxml

    property mode:
        """ Writing mode: overwrite or append. """
        def __get__(self):
            return self.mode

    property vars:
        """ Dictionary of variables to write. """
        def __get__(self):
            return self.vars

    property attrs:
        """ Dictionary of attributes to write. """
        def __get__(self):
            return self.attrs

    def __init__(self,char * fname,
                 bint is_noxml = True,
                 char * mode = "w",
                 MPI_Comm comm = MPI_COMM_WORLD):
        self.fname = fname
        self.method = <bytes>""
        self.method_params = <bytes>""
        self.is_noxml = is_noxml
        self.mode = mode
        self.comm = comm
        self.vars = dict()
        self.attrs = dict()

    ##def __var_factory__(self, name, value):
    ##    print "var_factory:", name, value
    ##
    ##def __attr_factory__(self, name, value):
    ##    print "attr_factory:", name, value

    def declare_group(self, char * gname,
                      char * method = "POSIX1",
                      char * method_params = ""):
        """
        Define a group associated with the file.

        Args:
            gname (str): group name.
            method (str, optional): Adios write method (default: 'POSIX1')
            method_params (str, optional): parameters for the write method (default: '')

        Example:

        >>>  fw.declare_group('group', method='MPI_, method_params='verbose=3')

        """
        self.gid = declare_group(gname, "", 1)
        self.gname = gname
        self.method = method
        self.method_params = method_params
        select_method(self.gid, self.method, self.method_params, "")

    def define_var(self, char * varname,
                   ldim = tuple(),
                   gdim = tuple(),
                   offset = tuple()):
        """
        Define a variable associated with the file.

        Args:
            varname (str): variable name.
            ldim (tuple, optional): local dimension (default: tuple())
            gdim (tuple, optional): global dimension (default: tuple())
            offset (tuple, optional): offset (default: tuple())

        Example:

        Write 'temperature' variable of size of 2x3 array.

        >>>  fw.define_var ('temperature', (2,3))

        """
        self.vars[varname] = varinfo(varname, ldim, gdim, offset)

    def define_attr(self, char * attrname):
        """
        Define attribute in the file.

        Args:
            attrname (str): attribute name.
        """

        self.attrs[attrname] = attrinfo(attrname, is_static=True)

    def define_dynamic_attr(self, char * attrname,
                            char * varname,
                            dtype):
        self.attrs[attrname] = attrinfo(attrname, varname, dtype, is_static=False)
    def __setitem__(self, name, val):
        if self.vars.has_key(name):
            self.vars[name] = val
        elif self.attrs.has_key(name):
            self.attrs[name] = val
        else:
            self.vars[name] = val

    def __getitem__(self, name):
        if self.vars.has_key(name):
            return self.vars[name]
        elif self.attrs.has_key(name):
            return self.attrs[name]
        else:
            raise KeyError(name)

    def close(self):
        """
        Write variables and attributes to a file and close the writer.
        """
        fd = open(self.gname, self.fname, self.mode)

        extra_vars = dict()
        extra_attrs = dict()

        for key, val in self.vars.iteritems():
            if not isinstance(val, varinfo):
                n = np.array(val)
                extra_vars[key] = varinfo(key, n.shape)
                extra_vars[key].value = val
            else:
                if self.is_noxml: val.define(self.gid)

        for key, val in extra_vars.iteritems():
            if self.is_noxml: val.define(self.gid)
            self.vars[key] = val

        for key, val in self.attrs.iteritems():
            if not isinstance(val, attrinfo):
                extra_attrs[key] = attrinfo(key, val, np.array(val).dtype)
            else:
                if self.is_noxml: val.define(self.gid)

        for key, val in extra_attrs.iteritems():
            if self.is_noxml: val.define(self.gid)

        groupsize = 0
        for var in self.vars.values():
            groupsize = groupsize + var.bytes()

        set_group_size(fd, groupsize)

        for var in self.vars.values():
            var.write(fd)

        close(fd)

    def __repr__(self):
        return ("AdiosWriter (fname=%r, gname=%r, "
                "method=%r, method_params=%r, vars=%r, attrs=%r, mode=%r)") % \
                (self.fname,
                 self.gname,
                 self.method,
                 self.method_params,
                 self.vars.keys(),
                 self.attrs.keys(),
                 self.mode)

cdef class attrinfo(object):
    cdef bytes name
    cdef bint is_static # Use define_byvalue, if True
    cdef dtype
    cdef value # Either varname or nparray

    property name:
        def __get__(self):
            return self.name

    property is_static:
        def __get__(self):
            return self.is_static

    property dtype:
        def __get__(self):
            return self.dtype

    property value:
        def __get__(self):
            return self.value

        def __set__(self, value):
            self.value = value

    def __init__(self, char * name,
                 value = None,
                 dtype = None,
                 bint is_static = True):
        self.name = name
        self.value = value
        self.dtype = dtype
        self.is_static = is_static

    def define(self, int64_t gid):
        if self.is_static:
            if self.value is None:
                raise TypeError("Value is none")

            define_attribute_byvalue(gid, self.name, "", self.value)
        else:
            ##atype = np2adiostype(np.dtype(self.dtype))
            ##define_attribute(gid, self.name, "",
            ##                 atype, "", str(self.value))
            raise NotImplementedError

    def __repr__(self):
        return ("AdiosAttrinfo (name=%r, is_static=%r, value=%r, dtype=%r)") % \
                (self.name,
                 self.is_static,
                 self.value,
                 self.dtype)

cdef class varinfo(object):
    cdef bytes name
    cdef public ldim
    cdef public gdim
    cdef public offset
    cdef public value

    def __init__(self, char * name,
                 ldim = tuple(),
                 gdim = tuple(),
                 offset = tuple(),
                 value = None):
        self.name = name
        self.ldim = ldim
        self.gdim = gdim
        self.offset = offset
        self.value = value

    def define(self, int64_t gid):
        if self.value is None:
            raise TypeError("Value is none")

        ldim_ = self.ldim
        if isinstance(self.ldim, (tuple, list)):
            ldim_ = tuple(self.ldim)

        gdim_ = self.gdim
        if isinstance(self.gdim, (tuple, list)):
            gdim_ = tuple(self.gdim)

        offset_ = self.offset
        if isinstance(self.offset, (tuple, list)):
            offset_ = tuple(self.offset)

        val_ = self.value
        if not isinstance(self.value, (np.ndarray)):
            val_ = np.array(self.value)

        atype = np2adiostype(val_.dtype)
        ## No space allowed
        define_var(gid, self.name, "", atype,
                   str(ldim_).replace(' ', '').strip('(,)'),
                   str(gdim_).replace(' ', '').strip('(,)'),
                   str(offset_).replace(' ', '').strip('(,)'))

    def bytes(self):
        val_ = self.value
        if not isinstance(self.value, (np.ndarray)):
            val_ = np.array(self.value)

        return val_.size * val_.itemsize

    def write(self, int64_t fd):
        val_ = self.value
        if not isinstance(self.value, (np.ndarray)):
            val_ = np.array(self.value)

        write(fd, self.name, val_)

    def __repr__(self):
        return ("AdiosVarinfo (name=%r, ldim=%r, gdim=%r, offset=%r, value=%r)") % \
                (self.name, self.ldim, self.gdim, self.offset, self.value)

## Aliases
File = file
Var = var
Attr = attr
Writer = writer
Attrinfo = attrinfo
Varinfo = varinfo

## ====================
## ADIOS Global functions
## ====================

def readvar(fname, varname):
    """ Retrieve a variable value from an Adios file.

    Args:
        fname (str): Adios file name
        varname (str): Variable name to retrieve

    Returns:
        NumPy ndarray: variable value
    """
    f = file(fname, comm=MPI_COMM_SELF)
    if not f.vars.has_key(varname):
        raise KeyError(varname)

    v = f.vars[varname]
    return v.read(from_steps=0, nsteps=v.nsteps)

def bpls(fname):
    """ Return meta data of an Adios file as a Python dictionary object.

    Args:
        fname (str): Adios file name

    Returns:
        dict: Adios file meta data
    """
    f = file(fname, comm=MPI_COMM_SELF)
    return {'nvars': f.nvars,
            'nattrs': f.nattrs,
            'vars': tuple([ k for k in f.vars.iterkeys() ]),
            'attrs': tuple([ k for k in f.attrs.iterkeys() ]),
            'time_steps': (f.current_step, f.last_step),
            'file_size': f.file_size}
