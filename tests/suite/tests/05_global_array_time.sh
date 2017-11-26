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
cp $SRCDIR/programs/examples/global_array_time/global_array_time_aggr_C.xml .
cp $SRCDIR/programs/examples/global_array_time/global_array_time_write_multifile_C .

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
    exit 1
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

###################################################
# run the time-aggregation test

echo "Run C global_array_time_write_C with time-aggregation turned on in subdir time_aggr/ with various buffer sizes:"
mkdir -p time_aggr
pushd time_aggr >/dev/null
$TRUNKDIR/utils/bpls/bpls -la ../global_array_time_C.bp -D -d temperature -n 10 > c_bpls_non_aggr.txt

for BUFSIZE in 3000 300000 10000 1000 0 ; do
            
    echo "  Run time-aggregation with buffer size = $BUFSIZE"
    cat ../global_array_time_aggr_C.xml | sed -e "s/buffer-size=\"[0-9]*\"/buffer-size=\"$BUFSIZE\"/" > global_array_time_C.xml
    $MPIRUN $NP_MPIRUN $PROCS $EXEOPT ../global_array_time_write_C
    EX=$?
    if [ ! -f global_array_time_C.bp ]; then
        echo "ERROR: global_array_time_write_C failed. No BP file is created. Exit code=$EX"
        exit 1
    fi

    mv global_array_time_C.xml global_array_time_C_$BUFSIZE.xml
    mv global_array_time_C.bp  global_array_time_C_$BUFSIZE.bp

    echo "    Check output with bpls"
    $TRUNKDIR/utils/bpls/bpls -la global_array_time_C_$BUFSIZE.bp | grep -v -e endianness -e 'file size' > c_bpls_$BUFSIZE.txt
    diff -q c_bpls_$BUFSIZE.txt $SRCDIR/reference/global_array_time_bpls.txt
    if [ $? != 0 ]; then
        echo "ERROR: global_array_time_write_C with time aggregation with buffer size $BUFSIZE produced a file different from the reference."
        echo "Compare \"bpls -la $PWD/global_array_time_C_$BUFSIZE.bp | grep -v -e endianness -e 'file size'\" to reference $SRCDIR/reference/global_array_time_bpls.txt"
        exit 1
    fi

    echo "    Check output with bpls even more. Dump data and compare to non-time-aggregated version"
    $TRUNKDIR/utils/bpls/bpls -la    global_array_time_C_$BUFSIZE.bp -D -d temperature -n 10 > c_bpls_aggr_$BUFSIZE.txt
    diff -q c_bpls_non_aggr.txt c_bpls_aggr_$BUFSIZE.txt
    if [ $? != 0 ]; then
        echo "ERROR: global_array_time_write_C with time aggregation with buffer size $BUFSIZE produced a file different from the file without time aggregation."
        echo "Compare \"bpls -la $PWD/global_array_time_C.bp -D -d temperature -n 10\" "
        echo "to      \"bpls -la $PWD/../global_array_time_C_$BUFSIZE.bp -D -d temperature -n 10\" "
        exit 1
    fi
done


popd >/dev/null


##################################################################################
# run the time-aggregation test where multiple output files are written over time
# Check if we have the correct number of steps in each file

echo "Run C global_array_time_write_C with time-aggregation turned on in subdir time_aggr_multifile/ with various buffer sizes:"
mkdir -p time_aggr_multifile
pushd time_aggr_multifile >/dev/null

for BUFSIZE in 3000 300000 10000 1000 0 ; do
            
	echo "  Run multi-file time-aggregation with buffer size = $BUFSIZE"
    cat ../global_array_time_aggr_C.xml | sed -e "s/buffer-size=\"[0-9]*\"/buffer-size=\"$BUFSIZE\"/" > global_array_time_aggr_C.xml
	$MPIRUN $NP_MPIRUN $PROCS $EXEOPT ../global_array_time_write_multifile_C
    EX=$?
    if [ $EX != 0 ]; then
    	echo "ERROR: global_array_time_write_multifile_C failed. Exit code=$EX"
        exit 1
    fi

    if [ ! -f global_array_time_C_1.bp ]; then
        echo "ERROR: global_array_time_write_C failed. global_array_time_C_1.bp file is missing"
        exit 1
    fi
    if [ ! -f global_array_time_C_2.bp ]; then
        echo "ERROR: global_array_time_write_C failed. global_array_time_C_2.bp file is missing"
        exit 1
    fi
    if [ ! -f global_array_time_C_3.bp ]; then
        echo "ERROR: global_array_time_write_C failed. global_array_time_C_3.bp file is missing"
        exit 1
    fi

	ADIR=bufsize_$BUFSIZE
	mkdir $ADIR
    mv global_array_time_aggr_C.xml $ADIR
    mv global_array_time_C_*.bp  $ADIR
    pushd $ADIR >/dev/null 

    echo "    Check output with bpls"
    N1=`$TRUNKDIR/utils/bpls/bpls -l global_array_time_C_1.bp temperature | cut -d 'e' -f 5 | cut -d '*' -f 1`
    N2=`$TRUNKDIR/utils/bpls/bpls -l global_array_time_C_2.bp temperature | cut -d 'e' -f 5 | cut -d '*' -f 1`
    N3=`$TRUNKDIR/utils/bpls/bpls -l global_array_time_C_3.bp temperature | cut -d 'e' -f 5 | cut -d '*' -f 1`
	echo $N1 > steps.txt
	echo $N2 >> steps.txt
	echo $N3 >> steps.txt
			
    if [ "$N1" -ne "10" ]; then
    	echo "ERROR: Expected 10 steps in global_array_time_C_1.bp. Instead got $N1 steps."
        exit 1
    fi
    if [ "$N2" -ne "10" ]; then
    	echo "ERROR: Expected 10 steps in global_array_time_C_2.bp. Instead got $N2 steps."
        exit 1
    fi
    if [ "$N3" -ne "5" ]; then
    	echo "ERROR: Expected 5 steps in global_array_time_C_3.bp. Instead got $N3 steps."
        exit 1
    fi
    popd >/dev/null
    
done


popd >/dev/null



###################################################
# run the Fortran tests too if available
if [ $HAVE_FORTRAN != yes ]; then
    exit 0
fi


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


