#!/bin/bash

function usage {
    echo "USAGE : `basename $0`"
}

if [ $# -lt 1 ]  ; then
    usage
    exit 1
fi

SRCDIR=$1

cat << EOF > config.xml
<?xml version="1.0"?>
<adios-config host-language="C">
    <adios-group name="temperature">
    <var name="NX" type="integer"/>
    <var name="size" type="integer"/>
    <var name="rank" type="integer"/>
    <global-bounds dimensions="size,NX" offsets="rank,0">
       <var name="temperature" gwrite="t" type="double" dimensions="1,NX"/>
    </global-bounds>
    <attribute name="description" path="/temperature" value="Global array written from 'size' processes" type="string"/>
</adios-group>

<method group="temperature" method="POSIX"/>

<buffer size-MB="2" allocate-time="now"/>

</adios-config>
EOF

PYTHONPATH=.:$PYTHONPATH python $SRCDIR/adios_read_test.py
exit $?
