#!/bin/bash

function usage {
    echo "USAGE : `basename $0` ADIOSJAR SRCDIR TARGET"
}

if [ $# -lt 3 ]  ; then
    usage
    exit 1
fi

ADIOSJAR=$1
SRCDIR=$2
TARGET=$3

cat << EOF > config.xml
<?xml version="1.0"?>
<adios-config host-language="C">
    <adios-group name="temperature" coordination-communicator="comm" stats="On">
    <var name="NX" type="integer"/>
    <var name="size" type="integer"/>
    <var name="rank" type="integer"/>
    <global-bounds dimensions="size,NX" offsets="rank,0">
       <var name="temperature" gwrite="t" type="double" dimensions="1,NX"/>
       <var name="ascii" gwrite="t" type="byte" dimensions="1,NX"/>
    </global-bounds>
    <attribute name="description" path="/temperature" value="Global array written from 'size' processes" type="string"/>
</adios-group>

<method group="temperature" method="MPI">stripe_count=1;stripe_size=10485760;num_aggregators=2;merging_pgs=0</method>

<buffer max-size-MB="2"/>

</adios-config>
EOF

rm -f adios_global.bp
javac -classpath $ADIOSJAR -d . $SRCDIR/$TARGET.java
java -Djava.library.path=. -classpath $ADIOSJAR:. $TARGET
exit $?
