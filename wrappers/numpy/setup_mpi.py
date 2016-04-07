#!/usr/bin/env python
# Author:  Jong Choi
# Contact: choij@ornl.gov


import os
import sys
import getopt

## Credit: http://svn.apache.org/repos/asf/subversion/tags/0.20.1/subversion/bindings/swig/python/setup.py
def _do_usage():
  print "Usage: setup.py [OPTIONS] build"
  print "       setup.py install [--prefix PREFIX]"
  print "       setup.py install_lib [--install-dir DIR]"
  print ""
  print "Options:"
  print "   -I dir      " + \
        "search DIR for includes (multiple instances allowed)"
  print "   -L dir      " + \
        "search DIR for libraries (multiple instances allowed)"
  print "   -C option   " + \
        "pass OPTION to the compiler at compile time (multiple instances " + \
        "allowed)"
  print "   -R option   " + \
        "pass OPTION to the compiler at link time (multiple instances " + \
        "allowed)"
  sys.exit(0)

# Default option values
include_dirs = []
library_dirs = []
extra_compile_args = []
extra_link_args = []

# No args?  Give usage.
if len(sys.argv) < 2:
  _do_usage()

# Parse the command-line arguments, keeping what we want and letting
# distutils have the rest.  Distutils parameters should come after
# the target as in 'python setup.py build --prefix=/usr/local' and
# parameters for us should appear before the target as in
# 'python setup.py -I/usr/include build'.
options, leftovers = getopt.getopt(sys.argv[1:], "I:L:C:R:h",
                                   ["help"])
for option in options:
  if option[0] == '-I':
    include_dirs.append(option[1])
  if option[0] == '-L':
    library_dirs.append(option[1])
  if option[0] == '-C':
    extra_compile_args.append(option[1])
  if option[0] == '-R':
    extra_link_args.append(option[1])
  if option[0] == '-h':
    _do_usage()

  if option[0] == '--help':
    _do_usage()

  # All long options just get passed through
  if option[0][:2] == '--':
    leftovers.append(option[0])
    leftovers.append(option[1])
sys.argv[1:] = leftovers

from distutils.extension import Extension
import numpy as np

# Use mpi4py dist utils: https://bitbucket.org/mpi4py/mpi4py
from conf.mpidistutils import setup
#from distutils.core import setup
from distutils.spawn import find_executable
from distutils.core import Command

import subprocess

include_dirs.insert(0, np.get_include())
extra_compile_args.insert(0, '-Wno-uninitialized')
extra_compile_args.insert(0, '-Wno-unused-function')

m1 = Extension('adios_mpi.adios_mpi',
               sources=['adios_mpi.cpp'],
               define_macros=[],
               include_dirs = include_dirs,
               library_dirs = library_dirs,
               libraries = [],
               extra_objects = [],
               extra_compile_args = extra_compile_args,
               extra_link_args = extra_link_args)

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

NAME = 'adios_mpi'
DESCRIPTION = 'Python Module for Adios MPI'
AUTHOR = 'Jong Choi'
AUTHOR_EMAIL = 'choij@ornl.gov'
URL = 'http://www.olcf.ornl.gov/center-projects/adios/'

import re
module_file = open("src/__init__.py").read()
metadata = dict(re.findall("__([a-z]+)__\s*=\s*'([^']+)'", module_file))
VERSION = metadata['version']

setup(name = NAME,
      version = VERSION,
      description = DESCRIPTION,
      author = AUTHOR,
      author_email = AUTHOR_EMAIL,
      url = URL,
      cmdclass={'test': adios_test},
      executables = [],
      ext_modules = [m1],
      packages=['adios_mpi', 'adios_mpi._hl'],
      package_dir = {'adios_mpi': 'src_mpi', 'adios_mpi._hl': '_hl'},
      )
