#!/bin/bash

CONFIG_FILE=$1 ; shift
SERVER_CONTACT_INFO=$1 ; shift

CONFIG_FRAG=/tmp/config.frag.$PBS_JOBID

while [ ! -s ${SERVER_CONTACT_INFO} ]; do
  echo "waiting for contact info.  retry in 1 second..."
  sleep 1
done
ls -l ${SERVER_CONTACT_INFO}


cat >$CONFIG_FRAG<<HERE
<?xml version="1.0"?>
<nssi-config>
  <staging-group write-type="WRITE_DIRECT">
HERE

SVC_LIST=`cat ${SERVER_CONTACT_INFO} | awk '{ ML=$0" "ML } END { print ML }'`
for SVC in $SVC_LIST; do
SVC_NID=`echo $SVC | sed -e 's/\(.*\)@\(.*\)@\(.*\)@\(.*\)/\1/'`
SVC_PID=`echo $SVC | sed -e 's/\(.*\)@\(.*\)@\(.*\)@\(.*\)/\2/'`
SVC_ADDR=`echo $SVC | sed -e 's/\(.*\)@\(.*\)@\(.*\)@\(.*\)/\3/'`
SVC_PORT=`echo $SVC | sed -e 's/\(.*\)@\(.*\)@\(.*\)@\(.*\)/\4/'`
cat>>$CONFIG_FRAG<<HERE
    <staging-service nid="$SVC_NID" pid="$SVC_PID" hostname="$SVC_ADDR" port="$SVC_PORT" />
HERE
done

cat >>$CONFIG_FRAG<<HERE
  </staging-group>
</nssi-config>
HERE

mv $CONFIG_FRAG $CONFIG_FILE
