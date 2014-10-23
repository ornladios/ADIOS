## Chester, Titan
FC=ftn
## Rhea, Sith
#FC=mpif90
CC=cc
LINKER=${CC}
#ADIOS_DIR = /opt/adios/1.6.0
ADIOS_INC := -I/ccs/home/habbasi/work/ADIOS/src/public -I/ccs/home/habbasi/work/ADIOS/vanilla_build/
ADIOS_FLIB := -L/ccs/home/habbasi/work/ADIOS/vanilla_build -ladiosf ${MXML_LIB} 
ADIOS_LIB := -L/ccs/home/habbasi/work/ADIOS/vanilla_build -ladios ${MXML_LIB} 
ADIOSREAD_FLIB := -L/ccs/home/habbasi/work/ADIOS/vanilla_build -ladiosf ${MXML_LIB}
GPP = ${ADIOS_DIR}/bin/gpp.py


default: all
all: write

write: writer.c
	${CC} -g -c -o writer.o ${ADIOS_INC} writer.c
	${LINKER} -g -o writer writer.o ${ADIOS_LIB}

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

