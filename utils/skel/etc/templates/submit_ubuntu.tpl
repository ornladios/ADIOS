#!/bin/bash
# Submission script for the $$TARGET$$ target
# To make changes to this file, edit the appropriate submit_<target>.tpl file in ~/.skel/templates/

ulimit -c unlimited

mpd&

sleep 5   #Make sure mpd has a chance to start

$$START_TEST$$
echo $$METHOD$$
$$PRE_RM$$
set_method.sh $$METHOD$$ $$APP$$_skel.xml.in $$APP$$_skel.xml 
mpirun -n $$CORES_USED$$ ./$$EXEC$$
$$POST_RM$$

$$END_TEST$$

mpdallexit
