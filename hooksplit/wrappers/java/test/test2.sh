#!/bin/bash

function usage {
    echo "USAGE : `basename $0` ADIOSJAR SRCDIR TARGET"
}

if [ $# -lt 3 ]  ; then
    usage
    exit 1
fi

ADIOSJAR=$1
SRCDIR=$2
TARGET=$3

rm -f adios_noxml.bp
javac -classpath $ADIOSJAR -d . $SRCDIR/$TARGET.java
java -Djava.library.path=. -classpath $ADIOSJAR:. $TARGET
exit $?
