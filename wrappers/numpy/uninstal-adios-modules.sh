#/bin/bash

usage () {
    echo "USAGE: $0 [-h|-d]"
    exit
}

DRY_RUN=false

while getopts ":hd" arg; do
  case $arg in
    h) usage ;;
    d) DRY_RUN=true ;;
    *) usage ;; 
  esac
done

PRECMD=
[ $DRY_RUN = true ] && PRECMD=echo && echo "Dry-run mode is ON"

$PRECMD pip uninstall --yes adios adios_mpi
$PRECMD pip uninstall --yes adios adios_mpi

echo "Extra cleaning:"
DIR1=`python -c "import site; print(site.getsitepackages()[0])"`
DIR2=`python -m site --user-site`
for DIR in $DIR1 $DIR2; do
    echo "$DIR"
    if [ -d $DIR ]; then
        pushd $DIR >> /dev/null
        $PRECMD rm -rf adios*
        popd >> /dev/null
    fi
done
