# Copied from FindNumpy.cmake
# ---
#
# $Id: $
#
# Author(s):  Anton Deguet
# Created on: 2010-01-20
#
# (C) Copyright 2010 Johns Hopkins University (JHU), All Rights
# Reserved.
#
# --- begin cisst license - do not edit ---
# 
# This software is provided "as is" under an open source license, with
# no warranty.  The complete license can be found in license.txt and
# http://www.cisst.org/cisst/license.txt.
# 
# --- end cisst license ---
#
# File based on FindNUMARRAY distributed with ITK 3.4 (see itk.org)
#
# Main modifications:
# - use Numpy instead of Numarray for all naming
# - added path for Python 2.5 and 2.6
# - renamed python script generated (det_npp became determineNumpyPath)
# - use lower case for CMake commands and keywords
# - updated python script to use get_include, not get_numpy_include which is now deprecated
#
# ---
#
# Try to find mpi4py python package
# Once done this will define
#
# PYTHON_MPI4PY_FOUND        - system has numpy development package and it should be used
# PYTHON_MPI4PY_INCLUDE_DIR  - directory where the arrayobject.h header file can be found
#
#

# include this to handle the QUIETLY and REQUIRED arguments
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
include(GetPrerequisites)

if(PYTHON_EXECUTABLE)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/determineMPI4PyPath.py "try: import mpi4py; print mpi4py.get_include()\nexcept: pass\n")
    exec_program("${PYTHON_EXECUTABLE}"
                 ARGS "\"${CMAKE_CURRENT_BINARY_DIR}/determineMPI4PyPath.py\""
                 OUTPUT_VARIABLE MPI4PY_PATH
                 )
endif(PYTHON_EXECUTABLE)

find_path(PYTHON_MPI4PY_INCLUDE_DIR mpi4py/mpi4py.h
          "${MPI4PY_PATH}"
          "${PYTHON_INCLUDE_PATH}/mpi4py/"
          /usr/include/python2.6/mpi4py/
          /usr/include/python2.5/mpi4py/
          /usr/include/python2.4/mpi4py/
          /usr/include/python2.3/mpi4py/
          DOC "Directory where the mpi4py.h header file can be found. This file is part of the mpi4py package"
    )

##if(PYTHON_MPI4PY_INCLUDE_DIR)
##    set(PYTHON_MPI4PY_FOUND 1 CACHE INTERNAL "Python mpi4py development package is available")
##endif(PYTHON_MPI4PY_INCLUDE_DIR)

find_package_handle_standard_args(PYTHON_MPI4PY DEFAULT_MSG PYTHON_MPI4PY_INCLUDE_DIR)
