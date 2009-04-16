# module load lFLAGS = -g 

CC = mpicc
INCLUDES = -I${ADIOS_DIR}/include -I${ADIOS_DIR}/adios-0.9.8/ -I${ADIOS_DIR}/adios-0.9.8/src/
LIBS = -L${ADIOS_DIR}/lib -L${ADIOS_DIR}/adios-0.9.8/src -ladios -ladios_int -lmxml
SRCS = bm_read.c bm_config.xml
INCLUDES = -I../../src -I../../
LIBS = -L../../src -lbpread
all: bm_read

#bp_utils.o: bp_utils.c
#	${CC} -c bp_utils.o bp_utils.c ${INCLUDES}
#bp_apis.o: bp_apis.c
#	${CC} -c bp_apis.o bp_apis.c ${INCLUDES}
bm_read.o: bm_read.c
	${CC} -c bm_read.o bm_read.c ${INCLUDES}
bm_read: bm_read.o 
	${CC} -o bm_read bm_read.o ${LIBS}
#bm_read.o: bm_read.c
#	${CC} -c bm_read.o bm_read.c ${INCLUDES}
#bm_read: bm_read.o bp_apis.o bp_utils.o 
#	${CC} -o bm_read bm_read.o bp_apis.o bp_utils.o ${LIBS}
depend:
	makedepend ${SRCS}	
clean:
	rm *.o bm_read 

