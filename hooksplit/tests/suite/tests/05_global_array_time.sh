#!/bin/bash
#
# Test if adios can write and read global arrays over time correctly
# Uses codes from examples/C/global-array-time and examples/Fortran/global-array-time
#
# Environment variables set by caller:
# MPIRUN        Run command
# NP_MPIRUN     Run commands option to set number of processes
# MAXPROCS      Max number of processes allowed
# HAVE_FORTRAN  yes or no
# SRCDIR        Test source dir (.. of this script)
# TRUNKDIR      ADIOS trunk dir

PROCS=9

if [ $MAXPROCS -lt $PROCS ]; then
    echo "WARNING: Needs $PROCS processes at least"
    exit 77  # not failure, just skip
fi

# copy codes and inputs to . 
cp $TRUNKDIR/examples/C/global-array-time/adios_globaltime .
cp $TRUNKDIR/examples/C/global-array-time/adios_read_globaltime .
cp $TRUNKDIR/examples/C/global-array-time/adios_globaltime.xml .

echo "Run C adios_globaltime"
$MPIRUN $NP_MPIRUN $PROCS ./adios_globaltime
EX=$?
if [ ! -f adios_globaltime.bp ]; then
    echo "ERROR: C version of adios_globaltime failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav adios_globaltime.bp | grep -v endianness > c_bpls.txt
diff -q c_bpls.txt $SRCDIR/reference/global_array_time_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: C version of adios_globaltime produced a file different from the reference."
    echo "Compare \"bpls -lav $PWD/adios_globaltime.bp | grep -v endianness\" to reference $SRCDIR/reference/global_array_time_bpls.txt"
    exit 1
fi

echo "Run C adios_read_globaltime"
$MPIRUN $NP_MPIRUN $PROCS ./adios_read_globaltime | sort > c_read.txt
EX=$?
if [ $? != 0 ]; then
    echo "ERROR: C version of adios_read_globaltime exited with $EX"
    echo "Check $PWD/c_read.txt"
    exit 1
fi

echo "Check output"
diff -q c_read.txt $SRCDIR/reference/global_array_time_read.txt
if [ $? != 0 ]; then
    echo "ERROR: C version of adios_read_globaltime produced a file different from the reference."
    echo "$PWD/c_read.txt to reference $SRCDIR/reference/global_array_time_read.txt"
    exit 1
fi


if [ $HAVE_FORTRAN != yes ]; then
    exit 0
fi
# run the Fortran tests too if available

mv adios_globaltime.xml adios_globaltime_c.xml
mv adios_globaltime.bp adios_globaltime_c.bp
cp $TRUNKDIR/examples/Fortran/global-array-time/adios_globaltime adios_globaltime_f
cp $TRUNKDIR/examples/Fortran/global-array-time/adios_globaltime.xml .

echo "Run Fortran adios_globaltime_f"
$MPIRUN $NP_MPIRUN $PROCS ./adios_globaltime_f
EX=$?
if [ ! -f adios_globaltime.bp ]; then
    echo "ERROR: Fortran version of adios_globaltime failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav adios_globaltime.bp | grep -v endianness > f_bpls.txt
diff -q f_bpls.txt $SRCDIR/reference/global_array_time_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: Fortran version of adios_globaltime produced a file different from the reference."
    echo "Compare \"bpls -lav $PWD/adios_globaltime.bp | grep -v endianness\" to reference $SRCDIR/reference/global_array_time_bpls.txt"
    exit 1
fi


