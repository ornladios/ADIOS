GPP ?= gpp.py 
COMP=ftn
LINKER=${COMP}

default: all
all: coupling2D 

coupling2D: coupling_writer_2D coupling_writer_2D_base

coupling_writer_2D: coupling_writer_2D.F90 gwrite_writer2D.fh
	${COMP}   -g -c -o coupling_writer_2D.o ${ADIOS_INC} coupling_writer_2D.F90
	${LINKER}   -g -o coupling_writer_2D coupling_writer_2D.o ${ADIOS_FLIB} ${DART_FLIB}

gwrite_writer2D.fh: coupling2D_writer.xml
	${GPP} coupling2D_writer.xml

coupling_writer_2D_base: coupling_writer_2D_base.F90 
	${COMP}   -g -c -o coupling_writer_2D_base.o coupling_writer_2D_base.F90
	${LINKER}   -g -o coupling_writer_2D_base coupling_writer_2D_base.o 

clean:
	rm -f *.o *.mod *.fh core.* fort.*
	rm -f coupling_writer_2D coupling_reader_2D  cpes*
	rm -f coupling_writer_2D_base coupling_reader_2D_base
