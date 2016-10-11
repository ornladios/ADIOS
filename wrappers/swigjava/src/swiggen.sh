#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo DIR=$DIR
pushd $DIR >& /dev/null

rm -f ornl/adios/adioslib/*.c
rm -f ornl/adios/adioslib/*.java
CMD="swig -java -package ornl.adios.adioslib `adios_config -s -c` -o adiosJAVA_wrap.c -outdir ornl/adios/adioslib adios.i"
echo $CMD
$CMD
popd $DIR >& /dev/null
