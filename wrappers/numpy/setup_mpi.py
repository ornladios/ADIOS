#!/usr/bin/env python
# Author:  Jong Choi
# Contact: choij@ornl.gov

from distutils.extension import Extension
import numpy as np

# Use mpi4py dist utils: https://bitbucket.org/mpi4py/mpi4py
from conf.mpidistutils import setup
#from distutils.core import setup
from distutils.spawn import find_executable
from distutils.core import Command

import subprocess

m1 = Extension('adios_mpi', 
               sources=['adios_mpi.cpp'], 
               define_macros=[],
               include_dirs = [np.get_include()],
               library_dirs = [],
               libraries = [],
               extra_objects = [],
               extra_compile_args = ['-Wno-#warnings', '-Wno-uninitialized', '-Wno-unused-function'])

cmd = find_executable("adios_config")
if cmd == None:
    sys.stderr.write(
        "adios_config is not installed nor found. "
        "Please install Adios or check PATH.\n")
    sys.exit(-1)

p = subprocess.Popen(["adios_config", "-c"], stdout=subprocess.PIPE)
for path in p.communicate()[0].strip().split(" "):
    if path.startswith('-I'):
        m1.include_dirs.append(path.replace('-I', '', 1))

p = subprocess.Popen(["adios_config", "-l"], stdout=subprocess.PIPE)
for path in p.communicate()[0].strip().split(" "):
    if path.startswith('-L'):
        m1.library_dirs.append(path.replace('-L', '', 1))
    if path.startswith('-l'):
        m1.libraries.append(path.replace('-l', '', 1))

class adios_test(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        import sys
        errno = subprocess.call([sys.executable, 'tests/test_adios_mpi.py', 'tests/config_mpi.xml'])
        raise SystemExit(errno)
    
setup(name = 'adios_mpi',
      version = '1.0.3',
      description = 'Python Module for Adios MPI',
      author = 'Jong Choi',
      author_email = 'yyalli@gmail.com',
      url = 'http://www.olcf.ornl.gov/center-projects/adios/',
      cmdclass={'test': adios_test},
      executables = [],
      ext_modules = [m1])
