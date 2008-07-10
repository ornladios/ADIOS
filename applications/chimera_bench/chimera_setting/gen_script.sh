#!/bin/sh

PBS_FILE=chimera_$1_$2_$3.pbs

cp wwww/chimera_setting/chimera.pbs $PBS_FILE

sed -i.bak -e "s/xxxx/$1/g" \
           -e "s/mmmm/$2/g" \
           -e "s/nnnn/$3/g" \
           $PBS_FILE

rm $PBS_FILE.bak

