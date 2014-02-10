## MPICC SETTINGS
#CC=mpicc
#CFLAGS=-g

## CC SETTINGS
CC=cc
CFLAGS = -g -O0

## Set ADIOS_DIR here or before doing make
#override ADIOS_DIR=/opt/adios/1.6.0
override ADIOS_INC:=$(shell ${ADIOS_DIR}/bin/adios_config -c)
override ADIOS_LIB:=$(shell ${ADIOS_DIR}/bin/adios_config -l)
override GPP := ${ADIOS_DIR}/bin/gpp.py

default: buginxml
all: default

buginxml:  buginxml.c
	${CC} -g -o buginxml ${ADIOS_INC} buginxml.c ${ADIOS_LIB} 

clean:
	rm -f *.o *.ch core.* 
	rm -f buginxml buginxml.bp
        
