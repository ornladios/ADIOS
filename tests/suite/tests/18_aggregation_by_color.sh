#!/bin/bash
#
# Test if adios aggregates data by color properly.
# Uses the example code from programs/examples/global_array/global_array_aggregate_by_color_C
#
# Environment variables set by caller:
# MPIRUN        Run command
# NP_MPIRUN     Run commands option to set number of processes
# MAXPROCS      Max number of processes allowed
# HAVE_FORTRAN  yes or no
# SRCDIR        Test source dir (.. of this script)
# TRUNKDIR      ADIOS trunk dir

PROCS=7
READPROCS=3

if [ $MAXPROCS -lt $PROCS ]; then
    echo "WARNING: Needs $PROCS processes at least"
    exit 77  # not failure, just skip
fi

# copy codes and inputs to . 
cp $SRCDIR/programs/examples/global_array/global_array_aggregate_by_color_C .

echo "Run C adios_global_aggregate_by_color"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./global_array_aggregate_by_color_C
EX=$?
if [ ! -f global_array_aggregate_by_color_C.bp ] \
|| [! -f global_array_aggregate_by_color_C.bp.dir/global_array_aggregate_by_color_C.bp.0 ] \
|| [! -f global_array_aggregate_by_color_C.bp.dir/global_array_aggregate_by_color_C.bp.1 ]; then
    echo "ERROR: adios_global_aggregate_by_color failed. No BP file or subfiles are created. Exit code=$EX"
    exit 1
fi
