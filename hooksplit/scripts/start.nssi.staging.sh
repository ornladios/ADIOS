#!/bin/sh

SERVICE_COUNT=$1 ; shift
SERVICE_CONFIG_FILE=$1 ; shift

rm $ADIOS_NSSI_CONTACT_INFO

aprun -n $SERVICE_COUNT -N 1 $ADIOS_DIR/bin/nssi-staging-server $SERVICE_CONFIG_FILE

while [ ! -s $ADIOS_NSSI_CONTACT_INFO ]
  do
  sleep 1
done
