#!/bin/bash

module unload mxml PE-pgi PE-gnu PE-intel
module load PE-gnu
module load mxml/2.7
export LD_LIBRARY_PATH=/sw/redhat6/fastbit/svn/rhel6_gnu4.7.2/lib:$LD_LIBRARY_PATH
