
etcdir=@prefix@/etc/skel

etc_DATA=\
    compiler_fragment.mk \
    settings


CC=@MPICC@
FC=@MPIFC@

all: compiler_fragment.mk

compiler_fragment.mk: Makefile
	grep -e'^prefix\ =' -e'^CC\ =' -e'^FC\ =' -e'^CFLAGS\ =' -e'^LDFLAGS\ =' -e'^FCFLAGS\ =' -e'^FCLIBS\ =' < Makefile > compiler_fragment.mk
