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
cp $SRCDIR/programs/examples/global_array_time/global_array_time_write_C .
cp $SRCDIR/programs/examples/global_array_time/global_array_time_read_as_file_C .
cp $SRCDIR/programs/examples/global_array_time/global_array_time_read_as_stream_C .
cp $SRCDIR/programs/examples/global_array_time/global_array_time_C.xml .

# Insert transform=X if requested by user
add_transform_to_xmls

echo "Run C global_array_time_write_C"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./global_array_time_write_C
EX=$?
if [ ! -f global_array_time_C.bp ]; then
    echo "ERROR: global_array_time_write_C failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -la global_array_time_C.bp | grep -v -e endianness -e 'file size' > c_bpls.txt
diff -q c_bpls.txt $SRCDIR/reference/global_array_time_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: global_array_time_write_C produced a file different from the reference."
    echo "Compare \"bpls -la $PWD/global_array_time_C.bp | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/global_array_time_bpls.txt"
fi

###################################################
echo "Run C global_array_time_read_as_file_C"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./global_array_time_read_as_file_C 
EX=$?
if [ $? != 0 ]; then
    echo "ERROR: global_array_time_read_as_file_C exited with $EX"
    echo "Check $PWD/c_read_as_file.txt"
    exit 1
fi

echo "Check output"
cat log_read_as_file_C.[0-9]* | grep -F -v -e "DEBUG:" -e "INFO:" -e "WARN:" > c_read_as_file.txt
diff -q c_read_as_file.txt $SRCDIR/reference/global_array_time_read_as_file.txt
if [ $? != 0 ]; then
    echo "ERROR: global_array_time_read_as_file_C produced a file different from the reference."
    echo "$PWD/c_read_as_file.txt to reference $SRCDIR/reference/global_array_time_read_as_file.txt"
    exit 1
fi

###################################################
echo "Run C global_array_time_read_as_stream_C"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./global_array_time_read_as_stream_C > c_read_as_stream.txt
EX=$?
if [ $? != 0 ]; then
    echo "ERROR: global_array_time_read_as_stream_C exited with $EX"
    echo "Check $PWD/c_read_as_stream.txt"
    exit 1
fi

echo "Check output"
cat log_read_as_stream_C.[0-9]* | grep -F -v -e "DEBUG:" -e "INFO:" -e "WARN:" > c_read_as_stream.txt
diff -q c_read_as_stream.txt $SRCDIR/reference/global_array_time_read_as_stream.txt
if [ $? != 0 ]; then
    echo "ERROR: global_array_time_read_as_stream_C produced a stream different from the reference."
    echo "$PWD/c_read_as_stream.txt to reference $SRCDIR/reference/global_array_time_read_as_stream.txt"
    exit 1
fi


if [ $HAVE_FORTRAN != yes ]; then
    exit 0
fi
###################################################
# run the Fortran tests too if available

cp $SRCDIR/programs/examples/global_array_time/global_array_time_write_F .
cp $SRCDIR/programs/examples/global_array_time/global_array_time_F.xml .

# Insert transform=X if requested by user
add_transform_to_xmls

echo "Run Fortran global_array_time_write_F"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./global_array_time_write_F
EX=$?
if [ ! -f global_array_time_F.bp ]; then
    echo "ERROR: Fortran version of global_array_time_write_C failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -la global_array_time_F.bp | grep -v -e endianness -e 'file size' > f_bpls.txt
diff -q f_bpls.txt $SRCDIR/reference/global_array_time_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: global_array_time_write_F produced a file different from the reference."
    echo "Compare \"bpls -la $PWD/global_array_time_F.bp | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/global_array_time_bpls.txt"
    exit 1
fi


