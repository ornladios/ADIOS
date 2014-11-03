## Chester, Titan
FC=ftn
## Rhea, Sith
#FC=mpif90
CC=cc
LINKER=${CC}
HOME=/home/hobbes/hobbes/ADIOS

# required dependencies
MXML_INC=-I/home/hobbes/hobbes/nvl/src/nvl/pisces/xpmem/include
MXML_LDFLAGS=-L/home/hobbes/hobbes/nvl/src/nvl/pisces/xpmem/lib
MXML_LIBS= -lmxml -lxpmem -lpthread -lm

ADIOS_DIR = ${HOME}
ADIOS_INC := -I${HOME}/src/ -I${HOME}/src/public -I${HOME}/vanilla_build ${MXML_INC}
ADIOS_FLIB := -L{HOME}/vanilla_build ${MXML_LDFLAGS} -ladiosf ${MXML_LIBS}
ADIOS_LIB := -L${HOME}/vanilla_build ${MXML_LDFLAGS} -ladios ${MXML_LIBS} 
ADIOSREAD_FLIB := -L${HOME} ${MXML_LDFLAGS} -ladiosf ${MXML_LIBS}
GPP = ${ADIOS_DIR}/bin/gpp.py


default: all
all: write

write:
	${CC} -g -c -o writer.o ${ADIOS_INC} writer.c
	${LINKER} -static -g -o writer writer.o ${ADIOS_LIB}
	${CC} -g -c -o read.o ${ADIOS_INC} read.c
	${LINKER} -static -g -o reader read.o ${ADIOS_LIB}

writer: writer.F90 gwrite_writer.fh
	${FC}   -g -c -o writer.o ${ADIOS_INC} writer.F90
	${LINKER}   -g -o writer writer.o ${ADIOS_FLIB}

gwrite_writer.fh: writer.xml
	${GPP} writer.xml

reader: reader.F90 
	${FC}   -g -c -o reader.o ${ADIOS_INC} reader.F90
	${LINKER}   -g -o reader reader.o ${ADIOSREAD_FLIB} 

clean:
	rm -f *.o *.mod *.fh core.* fort.*
	rm -f writer reader  writer*.bp writer*.h5 writer*.nc
	rm -f writer_base reader_base

