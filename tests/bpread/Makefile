# module load lFLAGS = -g 

CC = mpicc
#INCLUDES = -I${ADIOS_DIR}/include -I${ADIOS_DIR}/adios-0.9.8/ -I${ADIOS_DIR}/adios-0.9.8/src/
#LIBS = -L${ADIOS_DIR}/lib -L${ADIOS_DIR}/adios-0.9.8/src -ladios -ladios_int -lmxml
SRCS = bm_read.c
INCLUDES = -I../../src -I../../
LIBS = -L../../src -lbpread
all: bm_read

bm_read: bm_read.o 
	${CC} -o bm_read bm_read.o ${LIBS}
bm_read.o: bm_read.c
	${CC} -c bm_read.o bm_read.c ${INCLUDES}
depend:
	makedepend ${SRCS}	
clean:
	rm *.o bm_read 

