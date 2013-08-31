#!/bin/bash

all_bp_files=""

for bp_file in `ls output/ | grep ".bp$"`
do
	#echo $bp_file
	all_bp_files=$all_bp_files" output/"$bp_file
done

echo $all_bp_files

./adios_read_box $all_bp_files