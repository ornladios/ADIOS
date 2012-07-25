GPP ?= gpp.py 
COMP=ftn
LINKER=${COMP}

default: all
all: coupling2D 

coupling2D: coupling_reader_2D coupling_reader_2D_base

coupling_reader_2D: coupling_reader_2D.F90 
	${COMP}   -g -c -o coupling_reader_2D.o coupling_reader_2D.F90
	${LINKER}   -g -o coupling_reader_2D coupling_reader_2D.o ${ADIOSREAD_FLIB} 

coupling_reader_2D_base: coupling_reader_2D_base.F90 
	${COMP}   -g -c -o coupling_reader_2D_base.o coupling_reader_2D_base.F90
	${LINKER}   -g -o coupling_reader_2D_base coupling_reader_2D_base.o 

clean:
	rm -f *.o *.mod *.fh core.* fort.*
	rm -f coupling_writer_2D coupling_reader_2D  cpes*
	rm -f coupling_writer_2D_base coupling_reader_2D_base
