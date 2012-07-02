#!/bin/bash

function usage {
    echo "USAGE : `basename $0`"
}

if [ $# -lt 1 ]  ; then
    usage
    exit 1
fi

SRCDIR=$1

PYTHONPATH=.:$PYTHONPATH python $SRCDIR/adios_noxml_test.py
exit $?
