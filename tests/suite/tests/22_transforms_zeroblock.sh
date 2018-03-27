#!/bin/bash
#
# Test if adios can write and read back a transformed variable which contains a zero block
# Uses ../programs/zerolength
#
# Environment variables set by caller:
# MPIRUN        Run command
# NP_MPIRUN     Run commands option to set number of processes
# MAXPROCS      Max number of processes allowed
# HAVE_FORTRAN  yes or no
# SRCDIR        Test source dir (.. of this script)
# TRUNKDIR      ADIOS trunk dir

PROCS=2

if [ $MAXPROCS -lt $PROCS ]; then
    echo "WARNING: Needs $PROCS processes at least"
    exit 77  # not failure, just skip
fi

# copy codes and inputs to . 
cp $SRCDIR/programs/zerolength .

echo "Run test_singlevalue"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./zerolength
EX=$?

if [ $EX != 0 ]; then
    echo "ERROR: test_singlevalue failed with exit code=$EX"
    exit 1
fi

