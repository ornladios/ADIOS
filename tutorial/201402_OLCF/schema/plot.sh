#!/bin/bash 

plotter2d -f rectilinear2d.bp -v data -x X -y Y -rectilinear -o r

plotter2d -f tri2d.bp -v N -x points -y cells --triangle -o t -color hotdesaturated
