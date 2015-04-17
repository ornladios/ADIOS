** This cray folder is cloned from https://github.com/chuckatkins/miscelaneous-scripts

Cray-CMake-Modules
==================

A collection CMake code for use on various Cray supercomputing systems

* Platform: Cross compilation platform files not yet upstream
  * ComputeNodeLinux.cmake: Cray Compute Node Linux based on the Cray environment modules.
* ToolChain: CMake cross-compiling toolchain files
  * CrayPrgEnv-ToolChain.cmake: The Cray Programming Environment

The toolchain was tested on Titan, eos and chester.

How to use it:
in ADIOS top source file
$mkdir build
then
$cd build
$cmake .. -DCMAKE_TOOLCHAIN_FILE=../toolchain/cray/ToolChain/CrayPrgEnv-ToolChain.cmake
