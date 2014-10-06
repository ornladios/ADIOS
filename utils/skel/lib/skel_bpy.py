#!/usr/bin/env python

import yaml

import sys


# A representation of a bp file that is stored in a yaml document
class skel_bpy:
    def __init__ (self, filename):
        stream = file (filename, 'r')
        self.doc = yaml.load(stream)

        #print self.doc

        if self.get_language().lower() == "fortran":
            flip = True
        else:
            flip = False

        self.vars = {}
        for v in self.doc['variables']:
            name = v ['name']
            self.vars [name] = var (name, v, flip)
            #print "Added variable %s\n" % name

    def get_num_procs (self):
        return self.doc['procs']

    def get_steps (self):
        return self.doc.get('steps', 1)

    def get_vars (self):
        return filter (lambda x: not x.get_type() ==  'string', self.vars.values() )

    def get_var_dict (self):
        return self.vars

    def get_var (self, name):
        return self.vars[name]

    def get_var_names (self):
        return self.vars.keys()

    def get_language (self):
        return self.doc['lang']

    def get_group_name (self):
        return self.doc['group']

    def get_method (self):
        ret_val = self.doc.get ('method', None)
        if ret_val:
            return ret_val
        return "POSIX" # default method

    def get_parameters (self):
        ret_val = self.doc.get ('parameters', None)
        if ret_val:
            return ret_val
        return "" # default parameters 

    def get_buf_size_mb (self):
        ret_val = self.doc.get ('buf_size_mb', None)
        if ret_val:
            return ret_val
        return "50"

class var:
    def __init__ (self, name, vardict, flip):
        self.name = name
        self.vardict = vardict
        self.global_dims = None

        if self.is_scalar():
            # This is a scalar
            #print "%s is a scalar" % name
            return

        # If we're writing fortran, we need to flip the order of the dimensions, as they were reported by the C API
        if flip:
            if self.vardict.get('dims', None) != None and self.vardict['dims'] != 'scalar':
                self.vardict['dims'].reverse()
            if self.vardict.get('decomposition', None) != None:
                for proc in self.vardict['decomposition']:
                    proc.reverse()
            if self.vardict.get('decomp-procs', None) != None:
                self.vardict['decomp-procs'].reverse()
            if self.vardict.get('global-dims', None) != None:
                self.vardict['global-dims'].reverse()

        # We want to accept either the individual blocks of the decomposition
        # or do a regular decomposition based on given global dims and the
        # number of procs to split in each dimension

        # Case 1, we are given the individual blocks in vardict['decomposition']
        if self.get_decomposition() != None:
            # Calculate global dimensions
            # Find the global dims by looking for the largest end of block in each dimension 
            self.global_dims = []
            for i in range (self.get_ndims() ):
                self.global_dims.append (max ([x[i][1] for x in self.get_decomposition() ] ) + 1 )

        # Case 2, we are given the global dimensions and decomposition parameters
        elif self.get_decomp_procs() != None and self.get_global_dims() != None:
            dp = self.get_decomp_procs()
            gd = self.get_global_dims()

            # Sanity check
            if len(dp) != len(gd):
                print 'Dimensional mismatch, exiting'
                exit(13)
            for i in range(len(dp)):
                #divides evenly
                if gd[i] % dp[i] != 0:
                    print 'irregular decomposition not supported, exiting'
                    exit(17)

            # Calculate individual blocks
            self.vardict['dims'] = []
            procs = 1
            for i in range(len(dp)):
                self.vardict['dims'].append (gd[i] / dp[i])
                procs = procs * dp[i]


            self.vardict['decomposition'] = []
            #Assume there is at least one dimension specified
            for i in range(dp[0]):
                if len(dp) == 1: #1D
                    start0 = i * self.vardict['dims'][0]
                    end0 = (i+1) * self.vardict['dims'][0] - 1
                    self.vardict['decomposition'].append([[start0,end0]])
                else:
                    for j in range(dp[1]):
                        if len(dp) == 2: #2D
                            start0 = i * self.vardict['dims'][0]
                            end0 = (i+1) * self.vardict['dims'][0] - 1
                            start1 = j * self.vardict['dims'][1]
                            end1 = (j+1) * self.vardict['dims'][1] - 1
                            self.vardict['decomposition'].append([[start0,end0],[start1,end1]])
                        else:
                            for k in range(dp[2]):
                                if len(dp) == 2: #3D
                                    start0 = i * self.vardict['dims'][0]
                                    end0 = (i+1) * self.vardict['dims'][0] - 1
                                    start1 = j * self.vardict['dims'][1]
                                    end1 = (j+1) * self.vardict['dims'][1] - 1
                                    start2 = k * self.vardict['dims'][2]
                                    end2 = (k+1) * self.vardict['dims'][2] - 1
                                    self.vardict['decomposition'].append([[start0,end0],[start1,end1],[start2,end2]])
                                else:
                                    print "Failed sanity check, more than 3 dimensions specified"
                                    exit (13) 
                        
        # Case 3, no decomposition given
        else:
            #What to do here? PANIC!
            print "No decomposition found, exiting."
            exit(11)


    def get_name (self):
        return self.name

    def get_safe_name (self):
        # Remove slashes
        return 'SLASH'.join(self.get_name().split('/'))

    def get_type (self):
        return self.vardict['type']

    def get_lang_type (self, lang):
        #print "getting type for lang"
        if lang == 'C' or lang == 'c':
            return self.get_c_type()
        else:
            return self.get_fortran_type()
          

    def get_adios_type (self):

        self.ftypes = {
            'string' : 'adios_string',
            'byte' : 'adios_byte',
            'integer*1' : 'adios_byte',
            'short' : 'adios_short',
            'integer*2' : 'adios_short',
            'integer' : 'adios_integer',
            'integer*4' : 'adios_integer',
            'long' : 'adios_long',
            'long long' : 'adios_long',
            'integer*8' : 'adios_long',
            'unsigned byte' : 'adios_unsigned_byte',
            'unsigned integer*1' : 'adios_unsigned_byte',
            'unsigned short' : 'adios_unsigned_short',
            'unsigned integer*2' : 'adios_unsigned_short',
            'unsigned integer' : 'adios_unsigned_integer',
            'unsigned integer*4' : 'adios_unsigned_integer',
            'unsigned long' : 'adios_unsigned_long',
            'unsigned integer*8' : 'adios_unsigned_long',
            'float' : 'adios_real',
            'real' : 'adios_real',
            'real*4' : 'adios_real',
            'double' : 'adios_double',
            'real*8' : 'adios_double',
            'complex' : 'adios_complex',
            'double complex' : 'adios_double_complex'
        }

        return self.ftypes.get(self.get_type(), "UNKNOWN_TYPE")


    def get_fortran_type (self):
        self.ftypes = {
            'string' : 'string',
            'byte' : 'integer*1',
            'integer*1' : 'integer*1',
            'short' : 'integer*2',
            'integer*2' : 'integer*2',
            'integer' : 'integer*4',
            'integer*4' : 'integer*4',
            'long' : 'integer*8',
            'long long' : 'integer*8',
            'integer*8' : 'integer*8',
            'unsigned byte' : 'unsigned integer*1',
            'unsigned integer*1' : 'unsigned integer*1',
            'unsigned short' : 'unsigned integer*2',
            'unsigned integer*2' : 'unsigned integer*2',
            'unsigned integer' : 'unsigned integer*4',
            'unsigned integer*4' : 'unsigned integer*4',
            'unsigned long' : 'unsigned integer*8',
            'unsigned integer*8' : 'unsigned integer*8',
            'float' : 'real*4',
            'real' : 'real*4',
            'real*4' : 'real*4',
            'unsigned float' : 'unsigned real*4',
            'unsigned real' : 'unsigned real*4',
            'unsigned real*4' : 'unsigned real*4',
            'double' : 'real*8',
            'real*8' : 'real*8',
            'unsigned double' : 'unsigned real*8',
            'unsigned real*8' : 'unsigned real*8',
            'complex' : 'complex',
            'double complex' : 'double complex'
        }

        return self.ftypes.get(self.get_type(), "UNKNOWN_TYPE")

    def get_c_type (self):
        self.ftypes = {
            'string' : 'string',
            'byte' : 'unsigned char',
            'integer*1' : 'char',
            'short' : 'short',
            'integer*2' : 'short',
            'integer' : 'int',
            'integer*4' : 'int',
            'long' : 'long',
            'long long' : 'long',
            'integer*8' : 'long',
            'unsigned byte' : 'unsigned byte',
            'unsigned integer*1' : 'unsigned byte',
            'unsigned short' : 'unsigned short',
            'unsigned integer*2' : 'unsigned short',
            'unsigned integer' : 'unsigned integer',
            'unsigned integer*4' : 'unsigned integer',
            'unsigned long' : 'unsigned long',
            'unsigned integer*8' : 'unsigned long',
            'float' : 'float',
            'real' : 'float',
            'real*4' : 'float',
            'unsigned float' : 'unsigned float',
            'unsigned real' : 'unsigned float',
            'unsigned real*4' : 'unsigned float',
            'double' : 'double',
            'real*8' : 'double',
            'unsigned double' : 'unsigned double',
            'unsigned real*8' : 'unsigned double',
            'complex' : 'complex',
            'double complex' : 'double complex'
        }

        return self.ftypes.get(self.get_type(), "UNKNOWN_TYPE")

    def is_scalar (self):
        d = self.vardict.get('dims', None)
        return d == 'scalar'


    def get_dims (self):
        d = self.vardict.get('dims', None)
        if d == None or d == 'scalar':
            return None
        else:
            return self.vardict['dims']

    def get_dims_str (self):
        dims = self.get_dims()
        if dims == None:
            return ''
        else:
            return ','.join(str (d) for d in dims)

    def get_ndims (self):
        if self.get_dims() is None:
            return 0
        else:
            return len (self.get_dims()) 
    
    # For strings...
    def get_len (self):
        return self.vardict ['len']


    # This gives the size of one element of this type
    def get_unit_size (self):

        #print "Checking size of %s\n" % self.get_type()

        type = self.get_type()

        type_sizes = {

            'string' : 1,
            'byte' : 1,
            'integer*1' : 1,
            'short' : 2,
            'integer*2' : 2,
            'integer' : 4,
            'integer*4' : 4,
            'long' : 8,
            'long long' : 8,
            'integer*8' : 8,
            'unsigned byte' : 1,
            'unsigned integer*1' : 1,
            'unsigned short' : 2,
            'unsigned integer*2' : 2,
            'unsigned integer' : 4,
            'unsigned integer*4' : 4,
            'unsigned long' : 8,
            'unsigned integer*8' : 8,
            'float' : 4,
            'real' : 4,
            'real*4' : 4,
            'unsigned float' : 4,
            'unsigned real' : 4,
            'unsigned real*4' : 4,
            'double' : 8,
            'real*8' : 8,
            'unsigned double' : 8,
            'unsigned real*8' : 8,
            'complex' : 8,
            'double complex' : 16

        }
        size = type_sizes.get (type, None)
        if size is not None:
            return "%i" % size
        else:
            print "Unknown type: %s in get_unit_size()" % self.get_type()
            sys.exit()


    # This gives the size of a scalar or an array
    def get_size (self):
        if self.vardict.get ('dims', None) == 'scalar':
            return self.get_unit_size()
        else:
            return "%s * %s" % (self.get_unit_size(), '*'.join (str(x) for x in self.get_dims() ) ) 

    # Okay, so this is pinned to the decomposition data from the yaml file. We'll just take the max value
    # in each dimension...
    def get_global_dims (self):
        if self.global_dims:
            return self.global_dims

        return self.vardict.get ('global-dims', None)


    def get_global_dims_str (self):
        gd = self.get_global_dims()
        if gd is None:
            return ""
        else:
            return ",".join (str(glob_dim) for glob_dim in gd)

    def get_offset_values_str (self, var_name):
        dim = int (var_name[-1:] )
        return ",".join ( [str(x[dim][0]) for x in self.get_decomposition()] )

    # This is a list of var names to represent the global offsets in each dimension
    def get_offset_vars (self):
        return ["skel__offset_%s_%i" % (self.get_safe_name(),i) for i in range (self.get_ndims() )]

    def get_offset_vars_str (self):
        o = self.get_offset_vars()
        if o is None:
            return ""
        else:
            return ",".join (str(off) for off in o)

    def get_offsets (self):
        return self.vardict.get ('offsets', None)

    def get_offsets_str (self):
        o = self.get_offsets()
        if o is None:
            return ""
        else:
            return ",".join (str(off) for off in o)

    def has_global_bounds (self):
        return self.get_decomposition() is not None 

    def get_value (self):

        # Look for a single value first
        val = self.vardict.get ('value', None)
        if not val is None:
            return val

        # Multiple values not currently handled by template, just return first one.
        vals = self.vardict.get('values', None)
        if vals is None or length (vals) < 1:
            return None
        return vals[0]['val']

    def get_decomposition (self):
        return self.vardict.get ('decomposition', None)

    def get_decomp_procs (self):
        return self.vardict.get ('decomp-procs', None)


def main(argv=None):
    b = skel_bpy ("test.yaml")
    print "Num Procs is %d\n" % b.get_num_procs()

    vardict = b.get_vars()
    print vardict

    for vname in vardict:
        v = vardict [vname]
        print "%s: %s, %s\n" % (v.get_name(), v.get_type(), v.get_dims() )


if __name__ == "__main__":
    main()
