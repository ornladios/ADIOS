#!/bin/bash
# Submission script for the $$TARGET$$ target


#PBS -A $$ACCOUNT$$
#PBS -N $$JOB_NAME$$
#PBS -j oe
#PBS -l walltime=$$WALLTIME$$,size=$$CORES_TOTAL$$16$$
#PBS -l gres=widow2

cd $PBS_O_WORKDIR
date

export MPICH_UNEX_BUFFER_SIZE=251658240
export MPICH_MAX_SHORT_MSG_SIZE=64000


$$START_TEST$$
for i in {1..$$ITERATIONS$$}
do
echo $$METHOD$$
$$PRE_RM$$
set_method.sh $$METHOD$$ $$APP$$_skel.xml.in $$APP$$_skel.xml $$METHOD_PARAMS$$
aprun -n $$CORES_USED$$ ./$$EXEC$$
$$POST_RM$$
done

$$END_TEST$$
