# readme.txt
# Created on: Aug 21, 2013
# Author: Magda S. aka Magic Magg magg dot gatech at gmail.com
#
DESCRIPTION
===========
The idea is to test how to write in an appended mode and read from
that mode. This test is based on the Maya project.

It is designed for two methods:

1. MPI/ADIOS_READ_METHOD_BP
2. FLEXPATH/ADIOS_READ_METHOD_FLEXPATH

To switch between those two modes you need to run the make without or
with the CFLAGS set to -DFLEXPATH_METHOD. See the build section.

BUILD
=======
# you need to set the environment variables as Makefile uses those locations 
# to locate libraries and headers

export ADIOS_ROOT=/rock/opt/adios/git-dbg
export MXML_ROOT=/rock/opt/mxml/2.7
export MPI_ROOT=/rock/opt/openmpi/1.6.3
export EVPATH_ROOT=/rock/opt/evpath

# in certain cases you might need the lustre directory (e.g., on kraken)
export LUSTRE_ROOT=/opt/cray/lustre-cray_ss_s/default

# build the MPI/ADIOS_READ_METHOD_BP
$ make -f Makefile.generic

# build FLEXPATH/ADIOS_READ_METHOD_FLEXPATH
$ make -f Makefile.generic CFLAGS="-DFLEXPATH_METHOD"

# should remove all unnecessary exec files 
$ make -f Makefile.generic clean

# cleans files hanging around after previous runs
$ make -f Makefile.generic clean_test

RUN
===== 
# to clean files hanging around after previous runs
$ make -f Makefile.generic clean_test

# you can run as many writers as you want and as many readers as you want
# they write and read independently; the default is to run one writer
# and one reader
$ ./writer
$ ./reader 

See Makefile for other options of running the test.

PLAYING WITH TEST CONFIGURATION
===============================
To play with the test configuration, you can modify macros in the cfg.h file

TIMESTEP_COUNT seems to be the only reasonable setting to play so far

NOTES
======
2013-08-28 Test passes with the MPI method on my laptop; it fails with the FLEXPATH method
enabled (branch v1.5.1)

$ ./reader
ERROR: FLEXPATH staging method does not support file mode for reading. Use adios_read_open() to open a staged dataset.



TROUBLESHOOTING
================

2013-08-06, ERROR adios_allocate_buffer(): insufficient memory

ERROR: adios_allocate_buffer (): insufficient memory: 5242880000 requested, 860221440 available.  Using available.

$ grep ADS_BUFFER_SIZE cfg.h
#define ADS_BUFFER_SIZE 50

Try changing the ADS_BUFFER_SIZE in cfg.h to a smaller value.

2013-0x-0y, txt files left

There might be text files left; they should be removed for the next run.

# EOF

# EOF

