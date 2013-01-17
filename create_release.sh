#!/bin/bash
#
# Create a adios-VERSION.tar.gz for release
#
# It performs a
#    make dist
#    extracts adios-VERSION.tar.gz
#    removes the research transport method files from adios-VERSION/src
#    remakes the tar.gz
#
# Requirements
#    directory should be configured (to have Makefiles)

if [ ! -f Makefile ]; then
    echo "Prepare your directory with autogen.sh and runconf/configure before running this script"
    exit 1
fi

VERSION=`grep "^VERSION =" Makefile |  sed -e "s/.*= //"`

if [ -z $VERSION ]; then
    echo "Could not get the VERSION = ... line from the Makefile."
    echo "Cannot proceed."
    exit 2
fi

echo "Run make dist..."
make dist &>create_release.log
EX=$?

if [ $EX != 0 ]; then
    echo "Make dist failed. See create_release.log for details"
    exit 3
fi

if [ ! -f adios-$VERSION.tar.gz ]; then
    echo "Strange: missing adios-$VERSION.tar.gz after make dist."
    echo "Give up"
    exit 4
fi

echo "Extract adios-$VERSION.tar.gz"
rm -rf adios-$VERSION
tar zxf adios-$VERSION.tar.gz

echo "Remove research transport methods"
rm -f adios-$VERSION/src/write/adios_adaptive.c
rm -f adios-$VERSION/src/write/adios_mpi_stripe.c
rm -f adios-$VERSION/src/write/adios_mpi_cio.c
rm -f adios-$VERSION/src/write/adios_mpi_stagger.c
rm -f adios-$VERSION/src/write/adios_mpi_aggregate.c
rm -f adios-$VERSION/src/write/adios_mpi_amr1.c
rm -f adios-$VERSION/src/write/adios_provenance.c

echo "Add language wrappers"
#cp -rpP wrappers adios-$VERSION
#tar rf adios-$VERSION.tar --exclude .svn wrappers
#find wrappers | grep -v \.svn | cp  adios-$VERSION 
tar -cO wrappers --exclude "\.svn" | tar -x -C adios-$VERSION

echo "Add skel examples"
tar -cO examples/skel --exclude "\.svn" | tar -x -C adios-$VERSION

echo "Clean staging/coupling examples"
(cd examples/coupling; make distclean)
for sd in examples/staging; do
    (cd $sd; make distclean)
done
(cd examples/staging/staging_write; make -f Makefile.genarray_stream clean)

echo "Add staging/coupling examples"
tar -cO examples/staging examples/coupling  --exclude "\.svn" | tar -x -C adios-$VERSION

echo "Repack adios-$VERSION.tar.gz"
rm -rf adios-$VERSION.tar.gz
tar zcf adios-$VERSION.tar.gz adios-$VERSION
rm -rf adios-$VERSION

echo "Done"
ls -l adios-$VERSION.tar.gz

