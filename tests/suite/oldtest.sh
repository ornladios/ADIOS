#!/bin/bash
#
# Run this script in an interactive-job environment, with 16 cores at least.
# Define MPIRUN and MPIRUN_NP environment variables for running parallel programs
#   like mpirun and -np  or  aprun and -n
#
# Run the script in the job directory with full path.
#
# Log output of each test can be found in ./log/
#

function Usage() {
    echo "Usage:  MPIRUN=mpirun MPIRUN_NP=-np  <path>/test.sh   [pattern [pattern2 pattern3 ...]]
  <path>    is used to find all the test codes, inputs and reference outputs. 
            codes and inputs will be copied to the job directory
  pattern   if given, only those tests are executed, that match tests/*<pattern>*.sh
            otherwise all tests are executed.
            e.g. arguments   1[0-3] attr global   will execute all tests that have any of the 
            strings 10,11,12,13,attr or global in their name (and end with .sh)"
}

if [ "x$1" == "x-h" -o "x$1" == "x--help" ]; then
    Usage
    exit 0
fi

# check environment variables
if [ -z "$MPIRUN" ]; then
    echo "WARNING: Environment variable MPIRUN is not set. Set to 'mpirun'"
    MPIRUN=mpirun
fi

if [ -z "$MPIRUN_NP" ]; then
    echo "WARNING: Environment variable MPIRUN_NP is not set. Set to '-np'"
    MPIRUN_NP=-np
fi

# get source path
SRCDIR=`dirname $0`
TRUNKDISTANCE=../..
echo "Source dir = $SRCDIR"

if [ "${SRCDIR:0:1}" != "/" ]; then
    echo "WARNING: Jobs on some systems do not have access to any directory but to "
    echo "  parallel file systems. If the tests fail because of this, run this script"
    echo "  from a directory on the parallel file system and the script will copy the"
    echo "  data automatically."
    echo " "
fi

# find and list tests to be executed
TESTS=
if [ ! -z "$1" ]; then
    while [ ! -z "$1" ]; do
        TESTS1=`ls $SRCDIR/tests/*$1*.sh`
        if [ -z "$TESTS1" ]; then
            echo "ERROR: Did not find any test with name \"$1\" in $SRCDIR/tests."
        else
            TESTS="$TESTS $TESTS1"
        fi
        shift 1
    done
    if [ -z "$TESTS" ]; then
        exit 1
    fi
    echo "Tests found: $TESTS"
else
    TESTS=`ls $SRCDIR/tests/*.sh`
    if [ -z "$TESTS" ]; then
        echo "ERROR: Did not find any test in $SRCDIR/tests."
        exit 1
    fi
fi

# check if Fortran codes were built
S=`grep BUILD_FORTRAN_TRUE $SRCDIR/$TRUNKDISTANCE/Makefile`
if [ -z ${S#BUILD_FORTRAN_TRUE = } ]; then
    HAVE_FORTRAN=yes
else
    echo "WARNING: Fortran binaries are not built, so test will not use them"
    HAVE_FORTRAN=no
fi

# run tests one by one
for TESTSCRIPT in $TESTS; do
    TEST=`basename ${TESTSCRIPT%%.sh}`
    echo -n "  Test $TEST ... "
    rm -rf log.$TEST work.$TEST
    mkdir -p work.$TEST
    pushd work.$TEST >/dev/null
    if [ "${SRCDIR:0:1}" == "/" ]; then
        TESTSRCDIR=$SRCDIR       
        TRUNKDIR=$SRCDIR/$TRUNKDISTANCE
    else
        # we are one level deeper now
        TESTSRCDIR=../$SRCDIR
        TRUNKDIR=../$SRCDIR/$TRUNKDISTANCE
        TESTSCRIPT=../$TESTSCRIPT
    fi
    #echo "MPIRUN=$MPIRUN MPIRUN_NP=$MPIRUN_NP $TESTSCRIPT $TESTSRCDIR $TRUNKDIR $HAVE_FORTRAN &> ../log.$TEST"
    MPIRUN="$MPIRUN" MPIRUN_NP="$MPIRUN_NP" $TESTSCRIPT $TESTSRCDIR $TRUNKDIR $HAVE_FORTRAN &> ../log.$TEST
    EX=$?
    if [ $EX == 0 ]; then
        echo "OK"
        #rm -rf work.$TEST log.$TEST
    else
        echo "FAILED with exit code $EX. Check log.$TEST for details and work.$TEST/ for outputs."
    fi
    popd >/dev/null
done




