#!/usr/bin/env python

import yaml

import sys


# A representation of a bp file that is stored in a yaml document
class skel_bpy:
    def __init__ (self, filename):
        stream = file (filename, 'r')
        self.doc = yaml.load(stream)

        print self.doc

        self.vars = {}
        for v in self.doc['variables']:
            print "found variable %s" % v['name']
            name = v ['name']
            self.vars [name] = var (name, v)

    def get_num_procs (self):
        return self.doc['procs']

    def get_steps (self):
        return self.doc.get('steps', 1)

    def get_vars (self):
        return self.vars.values()

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

class var:
    def __init__ (self, name, vardict):
        self.name = name
        self.vardict = vardict

    def get_name (self):
        return self.name

    def get_safe_name (self):
        # Remove slashes
        return '_slash_'.join(self.get_name().split('/'))

    def get_type (self):
        return self.vardict['type']

    def get_fortran_type (self):
        self.ftypes = {
            "double": "real*8",
            "integer": "integer*4",
            "double complex": "double complex",
            "byte": "byte"
        }

        return self.ftypes.get(self.get_type(), "UNKNOWN_TYPE")

    def get_c_type (self):
        self.ftypes = {
            "double": "double",
            "integer": "int",
            "long": "long",
            "long long": "long long"
        }

        return self.ftypes.get(self.get_type(), "UNKNOWN_TYPE")

    def get_dims (self):
        if self.vardict['dims'] == 'scalar':
            return None
        else:
            return self.vardict['dims']

    def get_ndims (self):
        if self.get_dims() is None:
            return 0
        else:
            return len (self.get_dims()) 
         
    # This gives the size of one element of this type
    def get_unit_size (self):

        #print "Checking size of %s\n" % self.get_type()

        type = self.get_type()

        if type == "int" or type == "integer":
            return "4"
        elif type == "double":
            return "8"
        elif type == "double complex":
            return "16"
        elif type == "long long":
            return "8"
        elif type == "byte":
            return "1"
        else:
            print "Unknown type: %s in get_unit_size()" % self.get_type()
            sys.exit()


    # This gives the size of a scalar or an array
    def get_size (self):
        if self.vardict['dims'] == 'scalar':
            return self.get_unit_size()
        else:
            return "%s * %s" % (self.get_unit_size(), '*'.join (str(x) for x in self.get_dims() ) ) 

    def get_global_dims (self):
        return self.vardict.get ('global_dims', None)

    def get_global_dims_str (self):
        gd = self.get_global_dims()
        if gd is None:
            return ""
        else:
            return ",".join (str(glob_dim) for glob_dim in gd)

    def get_offsets (self):
        return self.vardict.get ('offsets', None)

    def get_offsets_str (self):
        o = self.get_offsets()
        if o is None:
            return ""
        else:
            return ",".join (str(off) for off in o)

    def has_global_bounds (self):
        return self.get_global_dims() is not None 

    def get_value (self):
        vals = self.vardict.get('values', None)
        if vals is None or length (vals) < 1:
            return None
        return vals[0]['val']

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
