## MPICC SETTINGS
#CC=mpicc
#CFLAGS=-g

## CC SETTINGS
CC=cc
CFLAGS = -g -O0

## Set ADIOS_DIR here or before doing make
#override ADIOS_DIR=/opt/adios/1.6.0
override ADIOS_INC=`${ADIOS_DIR}/bin/adios_config -c`
override ADIOS_LIB=`${ADIOS_DIR}/bin/adios_config -l`
override GPP = ${ADIOS_DIR}/bin/gpp.py

default: uniform2d rectilinear2d structured2d tri2d
all: default

uniform2d:  uniform2d.c
	${CC} -g -o uniform2d ${ADIOS_INC} uniform2d.c ${ADIOS_LIB} 

rectilinear2d:  rectilinear2d.c
	${CC} -g -o rectilinear2d ${ADIOS_INC} rectilinear2d.c ${ADIOS_LIB} 

structured2d:  structured2d.c
	${CC} -g -o structured2d ${ADIOS_INC} structured2d.c ${ADIOS_LIB} 

tri2d:  tri2d.c
	${CC} -g -o tri2d ${ADIOS_INC} tri2d.c ${ADIOS_LIB} 

clean:
	rm -f *.o *.ch core.* 
	rm -f uniform2d.bp rectilinear2d.bp structured2d.bp tri2d.bp 
	rm -f uniform2d rectilinear2d structured2d tri2d
        
