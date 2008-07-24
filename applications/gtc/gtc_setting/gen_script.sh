#!/bin/sh

PBS_FILE=$1_$2_$3_$4.pbs

cp wwww/gtc_setting/gtc.pbs $PBS_FILE

half_x=`expr $2 / 2 `

sed -i.bak -e "s/SSSS/$1/g" \
           -e "s/xxxx/$2/g" \
           -e "s/XXXX/$half_x/g" \
           -e "s/mmmm/$3/g" \
           -e "s/nnnn/$4/g" \
           $PBS_FILE

rm $PBS_FILE.bak

