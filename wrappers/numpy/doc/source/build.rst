.. _build:

Installation
============

Adios Python wrapper requires Adios built with the GNU C compiler with
relocatable codes. Add -fPIC flag to CFLAGS before configuring Adios.

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
