#!/bin/bash

mkdir output/

rm -rf conf/
cp -r conf_genarray3D/ conf/

for xml_file in `ls conf | grep ".xml$"`
do
	conf_name=`echo $xml_file | awk -F"." '{printf $1}'`
	echo $conf_name
	#mpirun -np 4 ./adios_global `echo $xml_file | awk -F"." '{printf $1}'`
	#mpirun -np 4 ./adios_global $conf_name
	mpirun -np 8 ./genarray3D 64 32 $conf_name
done
