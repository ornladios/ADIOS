# module load lFLAGS = -g 

CC = mpicc
FC = mpif90
#INCLUDES = -I${ADIOS_DIR}/include -I${ADIOS_DIR}/adios-0.9.8/ -I${ADIOS_DIR}/adios-0.9.8/src/
#LIBS = -L${ADIOS_DIR}/lib -L${ADIOS_DIR}/adios-0.9.8/src -ladios -ladios_int -lmxml
SRCS = bp_read_c.c
INCLUDES = -I../../src -I../../
LIBS = -L../../src -lbpread
all: bp_read_c bp_read_f genbp


bp_read_f: bp_read_f.o 
	${FC} -o bp_read_f bp_read_f.o ${LIBS}
bp_read_f.o: bp_read_f.F90
	${FC} -c bp_read_f.o bp_read_f.F90 ${INCLUDES}
bp_read_c: bp_read_c.o 
	${CC} -o bp_read_c bp_read_c.o ${LIBS}
bp_read_c.o: bp_read_c.c
	${CC} -c bp_read_c.o bp_read_c.c ${INCLUDES}
genbp: genbp.o
	${CC} -o genbp genbp.o ${LIBS} -ladios -L/autofs/na1_home/cgj/work/ewok/lib -lmxml 
genbp.o: genbp.c
	${CC} -c genbp.o genbp.c ${INCLUDES}

depend:
	makedepend ${SRCS}	
clean:
	rm *.o bp_read_c bp_read_f genbp

