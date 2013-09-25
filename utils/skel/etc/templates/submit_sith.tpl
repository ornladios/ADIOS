#!/bin/bash
# Submission script for the $$TARGET$$ target
# To make changes to this file, edit the appropriate submit_<target>.tpl file in ~/.skel/templates/

#PBS -A $$ACCOUNT$$
#PBS -N $$JOB_NAME$$
#PBS -j oe
#PBS -l walltime=$$WALLTIME$$,nodes=$$NODES_TOTAL$$32$$:ppn=32
#PBS -l gres=widow2

cd $PBS_O_WORKDIR

date

$$START_TEST$$
for i in {1..$$ITERATIONS$$}
do
echo $$METHOD$$
$$PRE_RM$$
./set_method.sh $$METHOD$$ $$APP$$_skel.xml.in $$APP$$_skel.xml 
mpirun -n $$CORES_USED$$ ./$$EXEC$$
mv skel_time.xml $$APP$$_skel_time_${PBS_JOBID}_${i}.xml
$$POST_RM$$
done

$$END_TEST$$



