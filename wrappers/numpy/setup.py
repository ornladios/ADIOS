#!/usr/bin/env python
# Author:  Jong Choi
# Contact: choij@ornl.gov

from distutils.extension import Extension
import numpy as np

# Use mpi4py dist utils: https://bitbucket.org/mpi4py/mpi4py
#from mpidistutils import setup
from distutils.core import setup
from distutils.spawn import find_executable
from distutils.core import Command

import subprocess
import sys

m1 = Extension('adios.adios',
               sources=['adios.cpp'],
               define_macros=[('_NOMPI', None)],
               include_dirs = [np.get_include()],
               library_dirs = [],
               libraries = [],
               extra_objects = [],
               extra_compile_args = ['-Wno-uninitialized',
                                     '-Wno-unused-function'])

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

class adios_test(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        ##import subprocess
        ##import sys
        ##errno = subprocess.call([sys.executable, 'tests/test_adios.py', 'tests/config.xml'])
        ##raise SystemExit(errno)
        import os
        import sys
        import unittest
        setup_file = sys.modules['__main__'].__file__
        setup_dir = os.path.abspath(os.path.dirname(setup_file))
        test_loader = unittest.defaultTestLoader
        test_runner = unittest.TextTestRunner()
        test_suite = test_loader.discover(os.path.join(setup_dir, 'test'))
        test_runner.run(test_suite)

NAME = 'adios'
DESCRIPTION = 'Python Module for Adios'
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
      ext_modules = [m1],
      packages=['adios', 'adios._hl'],
      package_dir = {'adios': 'src', 'adios._hl': '_hl'},
      )
