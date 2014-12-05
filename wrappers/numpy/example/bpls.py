#!/usr/bin/env python
"""
Example:

$ python ./bpls.py bpfile
"""

import adios as ad
import numpy as np
import getopt, sys

def usage():
    print "USAGE: %s filename" % sys.argv[0]
    
def main():
    fname = ""
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    else:
        fname = sys.argv[1]

    f = ad.file(fname)

    print "File info:"
    print "  %-18s %d" % ("of variables:", f.nvars)
    print "  %-18s %d - %d" % ("time steps:", f.current_step, f.last_step)
    print "  %-18s %d" % ("file size:", f.file_size)
    print "  %-18s %d" % ("bp version:", f.version)
    print ""
    
    for k in sorted(f.var.keys()):
        v = f.var[k]
        print "  %-17s  %-12s  %d*%s" % (np.typename(np.sctype2char(v.type)), v.name, v.nsteps, v.dims)
            

if __name__ == "__main__":
    main()
