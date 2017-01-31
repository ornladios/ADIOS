#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo DIR=$DIR
pushd $DIR >& /dev/null

rm -f ornl/adios/ext/*.c
rm -f ornl/adios/ext/*.java
CMD="swig -java -package ornl.adios.ext `adios_config -s -c` -o adiosJAVA_wrap.c -outdir ornl/adios/ext adios.i"
echo $CMD
$CMD
popd $DIR >& /dev/null
