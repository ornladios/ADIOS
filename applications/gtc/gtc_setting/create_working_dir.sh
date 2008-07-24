#!/bin/sh

if [ "X$WORKING_DIR" = "X" ]
then
    WORKING_DIR=/tmp/work/$LOGNAME/gtc_bench
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

cd $WORKING_DIR/gtc_run

for p in P2 P4 P8 P16 P32 
do
    mkdir $WORKING_DIR/gtc_run/$p
    d=$LOW_BOUND_PROC 
    while [ $d -le $UP_BOUND_PROC ]
    do
        mkdir $WORKING_DIR/gtc_run/$p/$d
        for m in MPI MPI_CIO POSIX NULL       
        do
            mkdir $WORKING_DIR/gtc_run/$p/$d/$m 
            exe_file=$WORKING_DIR/src/gtc 
            if [ "X$m" = "XPOSIX" ]
	    then
	        exe_file=$WORKING_DIR/src/gtc_posix
            fi 

            n=0
            while [ $n -lt $NUMBER_TEST ]
            do
                mkdir $WORKING_DIR/gtc_run/$p/$d/$m/$n
	        cd $WORKING_DIR/gtc_run/$p/$d/$m/$n
                ln -s $exe_file gtc
                cp $WORKING_DIR/gtc_setting/cleanup . 
                ./cleanup
	        cd $WORKING_DIR/gtc_run/$p/$d/$m
	        let n=n+1
            done
        done
        let d=d*2
    done
done

