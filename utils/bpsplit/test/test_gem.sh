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

split gem.bp gem_1_5.bp 1 5
split gem.bp gem_6_8.bp 6 8
split gem.bp gem_9_11.bp 9 -1

merge gem_1_5.bp gem_merged.bp
merge gem_6_8.bp gem_merged.bp
merge gem_9_11.bp gem_merged.bp

diff gem.bp gem_merged.bp &>/dev/null
if [ $? != 0 ]; then
    echo "Merged BP file gem_merged.bp differs from original gem.bp" 1>&2
    echo "Files to be examined:" 1>&2
    echo "  gem.bp -> splitted into" 1>&2
    echo "  gem_1_5.bp  gem_6_8.bp  gem_9_11.bp  -> merged into" 1>&2
    echo "  gem_merged.bp" 1>&2
    exit 1
else
    echo "Split/merge test on gem.bp was successful"
    rm -f gem_1_5.bp gem_6_8.bp gem_9_11.bp gem_merged.bp
fi


