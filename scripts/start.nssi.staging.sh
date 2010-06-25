#!/bin/sh

export ADIOS_SRC_PATH=$HOME/projects/adios/src/adios

SERVICE_COUNT=$1 ; shift
SERVICE_CONFIG_FILE=$1 ; shift

rm $ADIOS_NSSI_CONTACT_INFO

aprun -n $SERVICE_COUNT -N 1 $ADIOS_SRC_PATH/src/nssi-staging-server $SERVICE_CONFIG_FILE

while [ ! -s $ADIOS_NSSI_CONTACT_INFO ]
  do
#  printf "%10s \r" waiting
  sleep 1
done
echo "found $ADIOS_NSSI_CONTACT_INFO"
cat $$ADIOS_NSSI_CONTACT_INFO
