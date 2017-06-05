.. _build:

Installation
============

Adios Python wrapper requires Adios built with the GNU C compiler with
relocatable codes. Add -fPIC flag to CFLAGS before configuring Adios.

::

  $ ./configure --prefix=/dir/tp/install \
            MPICC=mpicc MPICXX=mpicxx MPIFC=mpif90 \
            CC=gcc CXX=g++ FC=gfortran \
            CFLAGS="-fPIC"

.. note:: Adios provides various functions (such as, staging, compression,
          hdf5 conversion, etc) which can be turned with proper configure options.
          Please check 'configure --help' for more options.

::

  $ ./configure --help


Quick install with pip
----------------------

ADIOS Python wrapper can be installed with pip. Check if pip is
installed already. Otherwise, install pip first
::

  $ wget https://bootstrap.pypa.io/get-pip.py
  $ sudo python get-pip.py



Or,
::

  $ python get-pip.py --user

to install in a local directory, $HOME/.local

Then, install Adios and Adios-MPI wrapper as follows:
::

  $ pip install adios
  $ pip install adios_mpi

If you want to install in a custom directory, use the following:
::

  $ pip install --install-option="--prefix=$PREFIX" adios
  $ pip install --install-option="--prefix=$PREFIX" adios_mpi


Trouble Shooting
----------------

Custom MPICC and MPICXX
^^^^^^^^^^^^^^^^^^^^^^^

If one needs to use a custom MPICC and MPICXX command (e.g., Titan),
then use the following command:
::

  $ pip --global-option=build_ext \
        --global-option=--mpicc=cc --global-option=--mpicxx=CC adios


"Could not find any downloads that satisfy the requirement" with pip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the error is caused by a certificate error, then, try
::

  $ wget http://curl.haxx.se/ca/cacert.pem
  $ pip --cert cacert.pem search adios
  $ pip --cert cacert.pem install adios

Custom MPI4py build
^^^^^^^^^^^^^^^^^^^

You may want to build mpi4py by yourself. The source code is available here[https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-2.0.0.tar.gz]

Then, run setup.py as follows:
::
  $ python setup.py build_ext --mpicc=cc --mpicxx=CC

Install in a cutom location and set PYTHONPATH env:
::
  $ PREFIX=/dir/to/install
  $ python setup.py install --prefix=$PREFIX
  $ export PYTHONPATH=$PREFIX/lib/python2.7/site-packages:$PYTHONPATH


Build script on Titan
^^^^^^^^^^^^^^^^^^^^^

The following script is a reference for building Adios python module on Titan:
::
  #!/bin/bash -ex
  module unload PrgEnv-cray
  module unload PrgEnv-pgi
  module unload PrgEnv-intel
  module unload PrgEnv-gnu
  module unload python
  module unload python_anaconda
  module unload python_anaconda_mpi4py
  module unload python_anaconda_adios

  module load PrgEnv-gnu
  module load python_anaconda
  module load python_anaconda_mpi4py
  module load adios/1.12-devel

  git fetch origin
  VER=`git describe --always`
  PVER=1.11.1.post2

  PREFIX=$WORLDWORK_CSC143/sw/python_anaconda_adios/$VER
  PIP=`which pip`

  [ ! -f cacert.pem ] && wget http://curl.haxx.se/ca/cacert.pem

  $PIP install -I -U \
      --global-option build_ext \
      --global-option -lrt --install-option="--prefix=$PREFIX" \
      --ignore-installed \
      --cert cacert.pem \
      adios==$PVER -v

  $PIP install -I -U \
      --global-option build_ext \
      --global-option=--mpicc=cc \
      --global-option=--mpicxx=CC \
      --global-option -lrt --install-option="--prefix=$PREFIX" \
      --ignore-installed \
      --cert cacert.pem \
      adios_mpi==$PVER -v
