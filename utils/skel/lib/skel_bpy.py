#!/usr/bin/env python

import yaml

# A representation of a bp file that is stored in a yaml document
class skel_bpy:
    def __init__ (self, filename):
        stream = file (filename, 'r')
        self.doc = yaml.load(stream)

        self.vars = {}
        for v in self.doc['variables']:
            name = v ['name']
            self.vars [name] = var (name, v)

    def get_num_procs (self):
        return self.doc['procs']

    def get_vars (self):
        return self.vars

class var:
    def __init__ (self, name, vardict):
        self.name = name
        self.vardict = vardict

    def get_name (self):
        return self.name

    def get_type (self):
        return self.vardict['type']

    def get_dims (self):
        if self.vardict['dims'] == 'scalar':
            return None
        else:
            return self.vardict['dims']



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
