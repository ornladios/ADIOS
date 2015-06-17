#!/bin/bash
#
# Test if adios can write and read scalars correctly
# Uses codes from examples/C/scalars and examples/Fortran/scalars 
#
# Environment variables set by caller:
# MPIRUN        Run command
# NP_MPIRUN     Run command's option to set number of processes
# MAXPROCS      Max number of processes allowed
# HAVE_FORTRAN  yes or no
# SRCDIR        Test source dir (.. of this script)
# TRUNKDIR      ADIOS trunk dir
# EXEOPT        Run command's option to set executable

PROCS=3

if [ $MAXPROCS -lt $PROCS ]; then
    echo "WARNING: Needs $PROCS processes at least"
    exit 77  # not failure, just skip
fi

# copy codes and inputs to . 
cp $SRCDIR/programs/examples/scalars/scalars_read_C .
cp $SRCDIR/programs/examples/scalars/scalars_write_C .
cp $SRCDIR/programs/examples/scalars/scalars_C.xml .

# Insert transform=X if requested by user
add_transform_to_xmls

echo "Run C scalars_write"
echo $MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./scalars_write_C
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./scalars_write_C
EX=$?
if [ ! -f scalars_C.bp ]; then
    echo "ERROR: C version of scalars_write failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav scalars_C.bp | grep -v -e endianness -e 'file size' > c_bpls.txt
diff -q c_bpls.txt $SRCDIR/reference/scalars_write_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: C version of scalars_write produced a file different from the reference."
    echo "Compare \"bpls -lav $PWD/scalars_C.bp | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/scalars_write_bpls.txt"
    exit 1
fi

echo "Run C scalars_read"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./scalars_read_C > c_read.txt
EX=$?
if [ $? != 0 ]; then
    echo "ERROR: C version of scalars_read failed with exit code $EX"
    exit 1
fi
diff -q c_read.txt $SRCDIR/reference/scalars_read_c.txt
if [ $? != 0 ]; then
    echo "ERROR: C version of scalars_read produced an output different from the reference."
    echo "Compare $PWD/c_read.txt reference $SRCDIR/reference/scalars_read_c.txt"
    exit 1
fi


if [ $HAVE_FORTRAN != yes ]; then
    exit 0
fi
# run the Fortran tests too if available

cp $SRCDIR/programs/examples/scalars/scalars_read_F .
cp $SRCDIR/programs/examples/scalars/scalars_write_F .
cp $SRCDIR/programs/examples/scalars/scalars_F.xml .

# Insert transform=X if requested by user
add_transform_to_xmls

echo "Run Fortran scalar_write"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./scalars_write_F
EX=$?
if [ ! -f scalars_F.bp ]; then
    echo "ERROR: Fortran version of scalars_write failed. No BP file is created. Exit code=$EX"
    exit 1
fi

echo "Check output with bpls"
$TRUNKDIR/utils/bpls/bpls -lav scalars_F.bp | grep -v -e endianness -e 'file size' > f_bpls.txt
diff -q f_bpls.txt $SRCDIR/reference/scalars_write_bpls.txt
if [ $? != 0 ]; then
    echo "ERROR: Fortran version of scalars_write produced a file different from the reference"
    echo "Compare \"bpls -lav scalars_F.bp | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/scalars_write_bpls.txt"
    exit 1
fi

echo "Run Fortran scalars_read"
$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ./scalars_read_F > f_read.txt
EX=$?
if [ $? != 0 ]; then
    echo "ERROR: Fortran version of scalars_read failed with exit code $EX"
    exit 1
fi
diff -q f_read.txt $SRCDIR/reference/scalars_read_f.txt
if [ $? != 0 ]; then
    echo "ERROR: Fortran version of scalars_read produced an output different from the reference."
    echo "Compare $PWD/f_read.txt reference $SRCDIR/reference/scalars_read_f.txt"
    exit 1
fi


