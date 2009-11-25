ADIOS_INC ?= -I../../src
ADIOS_FLIB ?= -L../../src -ladiosf -L/ccs/proj/e2e/pnorbert/mxml/lib -lmxml
ADIOSREAD_FLIB ?= -L../../src -ladiosreadf
#GPP ?= ../../utils/gpp/gpp.py
GPP ?= gpp.py
#COMP=mpif90
COMP=ftn
LINKER=${COMP}
DART_INC ?= 
DART_FLIB ?= -L/ccs/proj/e2e/pnorbert/coupling/src -ldart2

default:  coupling2Ddart
all: coupling2D coupling2Ddart

coupling2D: coupling_writer_2D coupling_reader_2D
coupling2Ddart: coupling_writer_2D_dart coupling_reader_2D_dart

# COUPLING2D DART 
coupling_writer_2D_dart: coupling_writer_2D_dart.F90 
	${COMP} -Mprof=func -g -c -o coupling_writer_2D_dart.o ${DART_INC} coupling_writer_2D_dart.F90
	${LINKER} -Mprof=func -g -o coupling_writer_2D_dart coupling_writer_2D_dart.o ${DART_FLIB}

coupling_reader_2D_dart: coupling_reader_2D_dart.F90 gwrite_reader2D.fh
	${COMP} -Mprof=func -g -c -o coupling_reader_2D_dart.o ${ADIOS_INC} ${DART_INC} coupling_reader_2D_dart.F90
	${LINKER} -Mprof=func -g -o coupling_reader_2D_dart coupling_reader_2D_dart.o ${ADIOS_FLIB} ${DART_FLIB} 


# COUPLING2D ADIOS FILE-BASED 
coupling_writer_2D: coupling_writer_2D.F90 gwrite_writer2D.fh
	${COMP} -g -c -o coupling_writer_2D.o ${ADIOS_INC} coupling_writer_2D.F90
	${LINKER} -g -o coupling_writer_2D coupling_writer_2D.o ${ADIOS_FLIB} 

coupling_reader_2D: coupling_reader_2D.F90 gread_writer2D.fh gwrite_reader2D.fh
	${COMP} -g -c -o coupling_reader_2D.o ${ADIOS_INC} coupling_reader_2D.F90
	${LINKER} -g -o coupling_reader_2D coupling_reader_2D.o ${ADIOS_FLIB} ${ADIOSREAD_FLIB}

gwrite_writer2D.fh: coupling2D.xml
	${GPP} coupling2D.xml

gread_writer2D.fh: coupling2D.xml
	${GPP} coupling2D.xml

gwrite_reader2D.fh: coupling2D.xml
	${GPP}  coupling2D.xml

clean:
	rm -f *.o *.mod *.fh core.*
	rm -f coupling_writer_2D coupling_reader_2D *.o *.mod *.fh core.*
	rm -f coupling_writer_2D_dart coupling_reader_2D_dart
        
