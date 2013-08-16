#!/bin/bash

export LIBGL_ALWAYS_INDIRECT=1
#export PATH=$PATH:/opt/plotter/bin
export PATH=$PATH:/ccs/proj/e2e/plotter/titan

plotter -f adios_global_no_xml.bp -v temperature \
   -o temperature -imgsize 720 600

plotter -f adios_global_no_xml.bp -v blocks \
   -o blocks -imgsize 720 600

plotter -f adios_global_no_xml.bp -v temperature -v blocks \
   -o multi -multiplot -imgsize 720 600 \
    -subtitle="" \
    -title="Temperature and Blocks" 

