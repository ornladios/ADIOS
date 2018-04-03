#!/usr/bin/env python

import numpy as np
import adios as ad
from scipy.io import netcdf
import argparse
import logging

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('BPFILE', help='Input Adios file')
    parser.add_argument('OUTFILE', help='Output NCDF file')
    parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

    fn = netcdf.netcdf_file(args.OUTFILE, 'w')
    with ad.file(args.BPFILE) as f:
        for v in f.var.values():
            logging.debug('Reading: %s (shape: %r)'%(v.name, v.shape))
            dnames = list()
            if v.nsteps > 0:
                dname = '_'.join([v.name, 'timedim'])
                dnames.append(dname)
                fn.createDimension(dname, v.nsteps)

            for d in range(v.ndim):
                dname = '_'.join([v.name, 'dim%d'%(d)])
                dnames.append(dname)
                fn.createDimension(dname, v.dims[d])
            
            nv = fn.createVariable(v.name, v.dtype, dnames)
            nv[...] = v[...]
    fn.close()

if __name__ == '__main__':
    main()
