#!/usr/bin/env python
# Author:  Jong Choi
# Contact: choij@ornl.gov

from distutils.extension import Extension
import numpy as np

# Use mpi4py dist utils: https://bitbucket.org/mpi4py/mpi4py
#from mpidistutils import setup
from distutils.core import setup
from distutils.spawn import find_executable

import subprocess
import sys

m1 = Extension('adios', 
               sources=['adios.cpp'], 
               define_macros=[('_NOMPI', None)],
               include_dirs = [np.get_include()],
               library_dirs = [],
               libraries = [],
               extra_objects = [])

cmd = find_executable("adios_config")
if cmd == None:
    sys.stderr.write(
        "adios_config is not installed nor found. "
        "Please install Adios or check PATH.\n")
    sys.exit(-1)

p = subprocess.Popen(["adios_config", "-c", "-s"], stdout=subprocess.PIPE)
for path in p.communicate()[0].strip().split(" "):
    if path.startswith('-I'):
        m1.include_dirs.append(path.replace('-I', '', 1))

p = subprocess.Popen(["adios_config", "-l", "-s"], stdout=subprocess.PIPE)
for path in p.communicate()[0].strip().split(" "):
    if path.startswith('-L'):
        m1.library_dirs.append(path.replace('-L', '', 1))
    if path.startswith('-l'):
        m1.libraries.append(path.replace('-l', '', 1))

setup(name = 'Adios',
      version = '1.0.1',
      description = 'Python Module for Adios',
      author = 'Jong Choi',
      author_email = 'yyalli@gmail.com',
      url = 'http://www.olcf.ornl.gov/center-projects/adios/',
      ext_modules = [m1])
