CC=gcc
CFLAGS=-O3 

all: libmpinull.a

clean:
	rm -f mpi-dummy.o libmpinull.a

mpi-dummy.o: mpi-dummy.c mpi-dummy.h
	${CC} -c ${CFLAGS} mpi-dummy.c

libmpinull.a: mpi-dummy.o
	ar qcv libmpinull.a mpi-dummy.o
	ranlib libmpinull.a
