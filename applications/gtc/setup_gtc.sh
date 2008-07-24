#!/bin/sh

# set working directory
if [ "X$1" = "X" ]
then
    export WORKING_DIR=/tmp/work/$LOGNAME/gtc
else
    export WORKING_DIR=$1
fi

if [ "X$2" = "X" ]
then
    export LOW_BOUND_PROC=512
else 
    export LOW_BOUND_PROC=$2
fi

if [ "X$3" = "X" ]
then
    export UP_BOUND_PROC=16384
else 
    export UP_BOUND_PROC=$3
fi

if [ "X$4" = "X" ]
then
    export NUMBER_TEST=5
else 
    export NUMBER_TEST=$4
fi

######################################
# Needed before any module commands
######################################
if [ -f /opt/modules/default/etc/modules.sh ]; then
. /opt/modules/default/etc/modules.sh
fi
if [ -f /etc/profile.d/modules.sh ]; then
. /etc/profile.d/modules.sh
fi

echo "GTC Benchmark is set to " $WORKING_DIR
echo "Number of processes from " $LOW_BOUND_PROC " to " $UP_BOUND_PROC
echo "For each (number of processes, ADIOS Method), run " $NUMBER_TEST " job(s)"
echo

sleep 2
export SVN_URL="https://svn.ccs.ornl.gov/svn-ewok/ADIOS/trunk/applications/gtc"
echo 
echo "Check out from SVN ("$SVN_URL") ..."
echo 
sleep 1

mkdir -p $WORKING_DIR
cd $WORKING_DIR
module load subversion
svn checkout $SVN_URL
mv gtc/* . 
rm gtc -rf

sleep 2
echo 
echo "Build executables in " $WORKING_DIR"/src/ ..." 
echo 
sleep 1

# build executables
cd src/
make gtc
make gtc_posix

sleep 2
echo
echo "Correct path in scripts ..."
echo
sleep 1

# correct paths in scripts and input files to refer to WORKING_DIR
cd $WORKING_DIR/gtc_setting

WD_STRING=`echo $WORKING_DIR|sed -e "s:\/:\\\\/:g"`
for i in gen_script.sh gtc.pbs post_script.sh
do
    sed -i.bak -e "s:wwww:$WD_STRING:g" ./$i
    rm $i.bak
done

sleep 2
echo
echo "Create sub-directories in $WORKING_DIR/gtc_run for job execution ..."
echo
sleep 1

# setup gtc_run
sh ./create_working_dir.sh

echo 
echo "Benchmark environment is setup in " $WORKING_DIR
echo 
