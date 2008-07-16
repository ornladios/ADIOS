#!/bin/sh

if [ "X$WORKING_DIR" = "X" ]
then
    WORKING_DIR=/tmp/work/$LOGNAME/chimera_bench
fi

if [ "X$NUMBER_TEST" = "X" ]
then
    NUMBER_TEST=5
fi

if [ "X$LOW_BOUND_PROC" = "X" ]
then
    LOW_BOUND_PROC=512
fi

if [ "X$UP_BOUND_PROC" = "X" ]
then
    UP_BOUND_PROC=16384
fi

cd $WORKING_DIR/chimera_run

d=$LOW_BOUND_PROC 
while [ $d -le $UP_BOUND_PROC ]
do
    mkdir $WORKING_DIR/chimera_run/$d
    for m in MPI MPI_CIO POSIX NULL       
    do
        mkdir $WORKING_DIR/chimera_run/$d/$m 
        exe_file=$WORKING_DIR/src/Execute/build/RadHyd3D_mpi 
        if [ "X$m" = "XPOSIX" ]
	then
	    exe_file=$WORKING_DIR/src/Execute/build/RadHyd3D_posix
        fi 

        n=0
	while [ $n -lt $NUMBER_TEST ]
        do
            mkdir $WORKING_DIR/chimera_run/$d/$m/$n
	    cd $WORKING_DIR/chimera_run/$d/$m/$n
            ln -s $WORKING_DIR/chimera_run/working_setting/Equation_of_State Equation_of_State
            ln -s $WORKING_DIR/chimera_run/working_setting/FullNet FullNet
            ln -s $WORKING_DIR/chimera_run/working_setting/MGFLD MGFLD
            mkdir -p Execute/build/Data3
            cd Execute/build/
            ln -s $exe_file RadHyd3D
            cp $WORKING_DIR/chimera_setting/cleanup Data3
            ./Data3/cleanup
	    cd $WORKING_DIR/chimera_run/$d/$m
	    let n=n+1
        done
    done
    let d=d*2
done


