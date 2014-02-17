
## MPIF90/MPICC SETTINGS
#FC=mpif90
#FFLAGS=-g
#CC=mpicc
#CFLAGS=-g 
#LDFLAGS=-g

## FTN/CC SETTINGS
FC=ftn
FFLAGS = -O3
CC=cc
CFLAGS = -g -O0
LDFLAGS =

export CC CFLAGS LDFLAGS
export FC FCFLAGS

## Set ADIOS_DIR here or before doing make
#override ADIOS_DIR=/opt/adios/1.6.0
override ADIOS_INC:=$(shell ${ADIOS_DIR}/bin/adios_config -c)
override ADIOS_CLIB:=$(shell ${ADIOS_DIR}/bin/adios_config -l)
override ADIOS_FLIB:=$(shell ${ADIOS_DIR}/bin/adios_config -l -f)
override GPP := ${ADIOS_DIR}/bin/gpp.py

export ADIOS_INC 
export ADIOS_CLIB 
export ADIOS_FLIB 
export GPP 


default: writer stage_write 
all: default

stage_write: stage_write.c
	make -f Makefile.stage_write 

writer: writer.F90
	make -f Makefile.writer 

clean:
	make -f Makefile.writer clean
	make -f Makefile.stage_write clean
	rm -f staging.o[0-9]*


