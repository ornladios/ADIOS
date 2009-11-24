ADIOS_INC ?= -I../../src
ADIOS_FLIB ?= -L../../src -ladiosf -L/ccs/proj/e2e/pnorbert/mxml/lib -lmxml
ADIOSREAD_FLIB ?= -L../../src -ladiosreadf
GPP ?= ../../utils/gpp/gpp.py

default: coupling2D
all: coupling2D

coupling2D: coupling_writer_2D coupling_reader_2D

coupling_writer_2D: coupling_writer_2D.F90 gwrite_writer2D.fh
	mpif90 -g -c -o coupling_writer_2D.o ${ADIOS_INC} coupling_writer_2D.F90
	mpif90 -g -o coupling_writer_2D coupling_writer_2D.o ${ADIOS_FLIB} 

coupling_reader_2D: coupling_reader_2D.F90 gread_writer2D.fh gwrite_reader2D.fh
	mpif90 -g -c -o coupling_reader_2D.o ${ADIOS_INC} coupling_reader_2D.F90
	mpif90 -g -o coupling_reader_2D coupling_reader_2D.o ${ADIOS_FLIB} ${ADIOSREAD_FLIB}

gwrite_writer2D.fh: coupling2D.xml
	${GPP} coupling2D.xml

gread_writer2D.fh: coupling2D.xml
	${GPP} coupling2D.xml

gwrite_reader2D.fh: coupling2D.xml
	${GPP}  coupling2D.xml

clean:
	rm -f coupling_writer_2D coupling_reader_2D *.o *.mod *.fh core.*
        
