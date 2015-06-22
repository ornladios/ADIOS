#!/bin/bash
#
# Test if adios can write and read local_array_C correctly
# Uses codes from examples/C/local_array_C 
#
# Environment variables set by caller:
# MPIRUN        Run command
# NP_MPIRUN     Run commands option to set number of processes
# MAXPROCS      Max number of processes allowed
# HAVE_FORTRAN  yes or no
# SRCDIR        Test source dir (.. of this script)
# TRUNKDIR      ADIOS trunk dir

PROCS=5

if [ $MAXPROCS -lt $PROCS ]; then
    echo "WARNING: Needs $PROCS processes at least"
    exit 77  # not failure, just skip
fi

# copy codes and inputs to . 
cp $SRCDIR/programs/examples/local_array/local_array_write_C .
cp $SRCDIR/programs/examples/local_array/local_array_read_C .
cp $SRCDIR/programs/examples/local_array/local_array_C.xml .

# Insert transform=X if requested by user
add_transform_to_xmls

echo "Run C local_array_write_C"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./local_array_write_C
EX=$?
if [ ! -f local_array_C.bp ]; then
    echo "ERROR: C version of local_array_write_C failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav local_array_C.bp | grep -v -e endianness -e 'file size' > c_bpls.txt
diff -q c_bpls.txt $SRCDIR/reference/local_array_C_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: C version of local_array_write_C produced a file different from the reference."
    echo "Compare \"bpls -lav $PWD/local_array_C.bp | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/local_array_C_bpls.txt"
    exit 1
fi

echo "Run C local_array_read_C"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./local_array_read_C 
EX=$?
if [ $? != 0 ]; then
    echo "ERROR: C version of local_array_read_C failed with exit code $EX"
    exit 1
fi

echo "Check output"
cat log_read_C.[0-9]* | grep -F -v -e "DEBUG:" -e "INFO:" -e "WARN:" > c_read.txt
diff -q c_read.txt $SRCDIR/reference/local_array_C_read.txt
if [ $? != 0 ]; then
    echo "ERROR: C version of local_array_read_C produced an output different from the reference."
    echo "Compare $PWD/c_read.txt reference $SRCDIR/reference/local_array_C_read.txt"
    exit 1
fi


