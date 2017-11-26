#!/bin/bash
#
# Test if adios can stream back a file where the first variable
# is written only in the first step. This bug was fixed after 1.11.1
# Uses ../programs/test_singlevalue
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
cp $SRCDIR/programs/test_singlevalue .

echo "Run test_singlevalue"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./test_singlevalue
EX=$?
if [ ! -f test_singlevalue.bp ]; then
    echo "ERROR: test_singlevalue failed at creating the BP file, test_singlevalue.bp. Exit code=$EX"
    exit 1
fi

if [ $EX != 0 ]; then
    echo "ERROR: test_singlevalue failed with exit code=$EX"
    exit 1
fi

