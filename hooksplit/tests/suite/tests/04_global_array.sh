#!/bin/bash
#
# Test if adios can write and read global arrays over time correctly
# Uses codes from examples/C/global-array and examples/Fortran/global-array
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
cp $TRUNKDIR/examples/C/global-array/adios_global .
cp $TRUNKDIR/examples/C/global-array/adios_read_global .
cp $TRUNKDIR/examples/C/global-array/adios_global.xml .

echo "Run C adios_global"
$MPIRUN $NP_MPIRUN $PROCS ./adios_global
EX=$?
if [ ! -f adios_global.bp ]; then
    echo "ERROR: C version of adios_global failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav adios_global.bp -d -n 10  | grep -v endianness > c_bpls.txt
diff -q c_bpls.txt $SRCDIR/reference/global_array_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: C version of adios_global produced a file different from the reference."
    echo "Compare \"bpls -lav $PWD/adios_global.bp -d -n 10 | grep -v endianness\" to reference $SRCDIR/reference/global_array_bpls.txt"
    exit 1
fi

echo "Run C adios_read_global"
$MPIRUN $NP_MPIRUN $READPROCS ./adios_read_global | sort > c_read.txt
EX=$?
if [ $? != 0 ]; then
    echo "ERROR: C version of adios_read_global exited with $EX"
    echo "Check $PWD/c_read.txt"
    exit 1
fi

echo "Check output"
diff -q c_read.txt $SRCDIR/reference/global_array_read.txt
if [ $? != 0 ]; then
    echo "ERROR: C version of adios_read_global produced a file different from the reference."
    echo "$PWD/c_read.txt to reference $SRCDIR/reference/global_array_read.txt"
    exit 1
fi


if [ $HAVE_FORTRAN != yes ]; then
    exit 0
fi
# run the Fortran tests too if available

mv adios_global.xml adios_global_c.xml
mv adios_global.bp adios_global_c.bp
cp $TRUNKDIR/examples/Fortran/global-array/adios_global adios_global_f
cp $TRUNKDIR/examples/Fortran/global-array/adios_global.xml .

echo "Run Fortran adios_global_f"
$MPIRUN $NP_MPIRUN $PROCS ./adios_global_f
EX=$?
if [ ! -f adios_global.bp ]; then
    echo "ERROR: Fortran version of adios_global failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav adios_global.bp -d -n 10 | grep -v endianness > f_bpls.txt
diff -q f_bpls.txt $SRCDIR/reference/global_array_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: Fortran version of adios_global produced a file different from the reference."
    echo "Compare \"bpls -lav $PWD/adios_global.bp -d -n 10 | grep -v endianness\" to reference $SRCDIR/reference/global_array_bpls.txt"
    exit 1
fi


