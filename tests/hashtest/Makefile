## mpicc SETTINGS
CC=mpicc
CFLAGS=-g 
LDFLAGS=-g

## cc SETTINGS
#CC=cc
#CFLAGS = -O3
#LDFLAGS =

#TARGET=gnu
TARGET=pgi

## Set ADIOS_DIR here or before doing make
#override ADIOS_DIR=/ccs/proj/e2e/pnorbert/ADIOS/sith.${TARGET}
#override ADIOS_INC=`${ADIOS_DIR}/bin/adios_config -c`
#override ADIOS_CLIB=`${ADIOS_DIR}/bin/adios_config -l`

INC = -I../../src/core
LIB = -L../../src/core

default: hashtest 
all: default

OBJS = hashtest.o ../../src/libadios_a-qhashtbl.o

%.o: %.c
	${CC} ${CFLAGS} -c ${INC}  $<

hashtest: ${OBJS} 
	${CC} ${LDFLAGS} -o hashtest ${OBJS} # ${ADIOS_CLIB}

clean:
	rm -f *.o *.mod *.fh core.*
	rm -f hashtest


distclean: clean
	rm -f log.* *.bp

