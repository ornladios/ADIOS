#!/bin/bash -v
rm -rf ./build
sed -f adios_mpi2serial.sed adios_mpi.pyx > adios.pyx
cython -X embedsignature=True --cplus adios.pyx    
cython -X embedsignature=True --cplus adios_mpi.pyx    
[ $? -ne 0 ] && exit $?
CC=clang CXX=clang++ python setup.py install --user
python setup_mpi.py install --user

#cd doc
#touch source/index.rst
#make html
#cd ..

