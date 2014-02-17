#!/bin/bash

#export LIBGL_ALWAYS_INDIRECT=1
#export PATH=$PATH:/opt/plotter/bin

plotter2d -f writer00.bp -v xy -o xy.00 -min -4.0 -max 4.0 -colormap HotDesaturated

plotter2d -f writer01.bp -v xy -o xy.01 -min -4.0 -max 4.0 -colormap HotDesaturated

