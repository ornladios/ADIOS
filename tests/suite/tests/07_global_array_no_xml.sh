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
cp $SRCDIR/programs/examples/global_array/global_array_write_noxml_C .
cp $SRCDIR/programs/examples/global_array/global_array_read_noxml_C .

echo "Run C global_array_write_noxml_C"
echo $MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./global_array_write_noxml_C
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./global_array_write_noxml_C
EX=$?
if [ ! -f global_array_noxml_C.bp ]; then
    echo "ERROR: global_array_write_noxml_C failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav global_array_noxml_C.bp -d -n 10 | grep -v -e endianness -e 'file size' > c_bpls.txt
diff -q c_bpls.txt $SRCDIR/reference/global_array_no_xml_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: global_array_write_noxml_C produced a file different from the reference."
    echo "Compare \"bpls -lav $PWD/global_array_noxml_C.bp -d -n 10 | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/global_array_no_xml_bpls.txt"
    exit 1
fi

echo "Run C global_array_read_noxml_C"
echo $MPIRUN $NP_MPIRUN $READPROCS $EXEOPT ./global_array_read_noxml_C 
$MPIRUN $NP_MPIRUN $READPROCS $EXEOPT ./global_array_read_noxml_C 
EX=$?
if [ $? != 0 ]; then
    echo "ERROR: global_array_read_noxml_C exited with $EX"
    echo "Check $PWD/log_read_noxml.* files."
    exit 1
fi

echo "Check output"
cat log_read_C.[0-9]* | grep -F -v -e "DEBUG:" -e "INFO:" -e "WARN:" > c_read.txt
diff -q c_read.txt $SRCDIR/reference/global_array_no_xml_read.txt
if [ $? != 0 ]; then
    echo "ERROR: global_array_read_noxml_C produced a file different from the reference."
    echo "$PWD/c_read.txt to reference $SRCDIR/reference/global_array_no_xml_read.txt"
    exit 1
fi


if [ $HAVE_FORTRAN != yes ]; then
    exit 0
fi
# run the Fortran tests too if available

cp $SRCDIR/programs/examples/global_array/global_array_write_noxml_F .

echo "Run Fortran adios_global_f"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./global_array_write_noxml_F
EX=$?
if [ ! -f global_array_noxml_F.bp ]; then
    echo "ERROR: global_array_write_noxml_F failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav global_array_noxml_F.bp -d -n 10 | grep -v -e endianness -e 'file size' > f_bpls.txt
diff -q f_bpls.txt $SRCDIR/reference/global_array_no_xml_f_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: global_array_write_noxml_F produced a file different from the reference."
    echo "Compare \"bpls -lav $PWD/global_array_F.bp -d -n 10 | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/global_array_no_xml_f_bpls.txt"
    exit 1
fi


