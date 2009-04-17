#!/bin/bash

function split() {
    local inpf=$1
    local outf=$2
    local n=$3
    local m=$4
    echo "../bpsplit -n $n -m $m $inpf $outf"
    ../bpsplit -n $n -m $m $inpf $outf 
    if [ $? != 0 ]; then
        echo "Command failed: ../bpsplit -v -n $n -m $m $inpf $outf" 1>&2
        exit 1
    fi
    if [ ! -f $outf ]; then
        echo "Output $outf is not created" 1>&2
        exit 2
    fi
}

function merge() {
    local inpf=$1
    local outf=$2
    echo "../bpappend $inpf $outf"
    ../bpappend $inpf $outf
    if [ $? != 0 ]; then
        echo "Command failed: ../bpappend $inpf $outf" 1>&2
        exit 1
    fi
    if [ ! -f $outf ]; then
        echo "Output $outf is not created" 1>&2
        exit 1
    fi
}

split flow11.bp flow_1_5.bp 1 5
split flow11.bp flow_6_8.bp 6 8
split flow11.bp flow_9_11.bp 9 -1

merge flow_1_5.bp flow_merged.bp
merge flow_6_8.bp flow_merged.bp
merge flow_9_11.bp flow_merged.bp

diff flow11.bp flow_merged.bp &>/dev/null
if [ $? != 0 ]; then
    echo "Merged BP file flow_merged.bp differs from original flow11.bp" 1>&2
    echo "Files to be examined:" 1>&2
    echo "  flow11.bp -> splitted into" 1>&2
    echo "  flow_1_5.bp  flow_6_8.bp  flow_9_11.bp  -> merged into" 1>&2
    echo "  flow_merged.bp" 1>&2
    exit 1
else
    echo "Split/merge test on flow11.bp was successful"
    rm -f flow_1_5.bp flow_6_8.bp flow_9_11.bp flow_merged.bp
fi


