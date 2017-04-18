#!/bin/bash
VERSION="1.0"

# bad arguments
readonly E_BADARGS=65


#  Currently attributes and maya_append do not work
all_tests="1D_arr_global 1D_arr_global_noxml scalar maya_append maya_noxml local_arrays global_range_select"

count=0
passed=0
failed=0
for t in $all_tests; do
    pushd $t > /dev/null
    let count++
    rm 	-rf *_writer_info.txt *_reader_info.txt test.bp *_reader_ready.txt *_writer_ready.txt
    ../run_test > /dev/null 2>&1
    rc=$?
    if [[ $rc != 0 ]]; then 
	echo Test $t FAILED
	let failed++
    else 
	echo Test $t PASSED
	let passed++
    fi
    popd > /dev/null
done
echo "Results:  $count tests run.  $passed passed,  $failed failed"

# eof