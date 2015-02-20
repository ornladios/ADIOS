#!/bin/bash

FN=${1:-g1.bp}
echo Filename=$FN

/ccs/proj/e2e/plotter/sith/plotter2d -f $FN -v xy -o xy -zonal -min -3 -max 2
/ccs/proj/e2e/plotter/sith/plotter2d -f $FN -v yz -o yz -zonal -min -3 -max 2
/ccs/proj/e2e/plotter/sith/plotter2d -f $FN -v xz -o xz -zonal -min -3 -max 2

