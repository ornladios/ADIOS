#!/bin/bash

export LIBGL_ALWAYS_INDIRECT=1
export PATH=$PATH:/ccs/proj/e2e/plotter/titan

for f in reader_0*.bp; do
    NNN=${f:7:3}
    echo
    echo Plot file reader_$NNN.bp
    plotter2d -f $f -v xy -o xy.$NNN -z -colormap HotDesaturated -min 0 -max 15
    plotter2d -f $f -v xy2 -o xy2.$NNN -z -colormap HotDesaturated -min 0 -max 2.75
done

