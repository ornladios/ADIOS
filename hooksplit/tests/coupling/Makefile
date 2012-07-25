ADIOS_INC ?= -I../../src
ADIOS_FLIB ?= -L../../src -ladiosf -L/ccs/proj/e2e/pnorbert/mxml/lib -lmxml
ADIOSREAD_FLIB ?= -L../../src -ladiosreadf
GPP ?= ../../utils/gpp/gpp.py
#GPP ?= gpp.py
#COMP=mpif90
#COMP=ftn -Mprof=func
COMP=ftn
LINKER=${COMP}
DART_INC ?= 
DART_FLIB ?= -LcompareMPI -ldart2
#DART_FLIB ?= -L/ccs/proj/e2e/dart/lib -ldart2

default: coupling2D
all: coupling2D coupling2Ddart 

coupling2D: coupling_writer_2D coupling_reader_2D
coupling3D: coupling_writer_3D coupling_reader_3D
coupling2Ddart: coupling_writer_2D_dart coupling_reader_2D_dart
coupling3Ddart: coupling_writer_3D_dart coupling_reader_3D_dart

# COUPLING3D NATIVE DART
coupling_writer_3D_dart: coupling_writer_3D_dart.F90 
	${COMP}  -g -c -o coupling_writer_3D_dart.o ${DART_INC} coupling_writer_3D_dart.F90
	${LINKER}  -g -o coupling_writer_3D_dart coupling_writer_3D_dart.o ${DART_FLIB}

coupling_reader_3D_dart: coupling_reader_3D_dart.F90 gwrite_reader3D.fh
	${COMP}  -g -c -o coupling_reader_3D_dart.o ${ADIOS_INC} ${DART_INC} coupling_reader_3D_dart.F90
	${LINKER}  -g -o coupling_reader_3D_dart coupling_reader_3D_dart.o ${ADIOS_FLIB} ${DART_FLIB} 


# COUPLING2D NATIVE DART
coupling_writer_2D_dart: coupling_writer_2D_dart.F90 
	${COMP}  -g -c -o coupling_writer_2D_dart.o ${DART_INC} coupling_writer_2D_dart.F90
	${LINKER}  -g -o coupling_writer_2D_dart coupling_writer_2D_dart.o ${DART_FLIB}

coupling_reader_2D_dart: coupling_reader_2D_dart_noadios.F90 gwrite_reader2D.fh
	${COMP}  -g -c -o coupling_reader_2D_dart.o ${ADIOS_INC} ${DART_INC} coupling_reader_2D_dart_noadios.F90
	${LINKER}  -g -o coupling_reader_2D_dart coupling_reader_2D_dart.o ${ADIOS_FLIB} ${DART_FLIB} 


# COUPLING3D ADIOS FILE-BASED 
coupling_writer_3D: coupling_writer_3D.F90 gwrite_writer3D.fh
	${COMP}  -g -c -o coupling_writer_3D.o ${ADIOS_INC} coupling_writer_3D.F90
	${LINKER}  -g -o coupling_writer_3D coupling_writer_3D.o ${ADIOS_FLIB} ${DART_FLIB}

coupling_reader_3D: coupling_reader_3D.F90 gread_writer3D.fh gwrite_reader3D.fh
	${COMP}  -g -c -o coupling_reader_3D.o ${ADIOS_INC} coupling_reader_3D.F90
	${LINKER}  -g -o coupling_reader_3D coupling_reader_3D.o ${ADIOS_FLIB} ${ADIOSREAD_FLIB} ${DART_FLIB}

gwrite_writer3D.fh: coupling3D.xml
	${GPP} coupling3D.xml

gread_writer3D.fh: coupling3D.xml
	${GPP} coupling3D.xml

gwrite_reader3D.fh: coupling3D.xml
	${GPP}  coupling3D.xml

# COUPLING2D ADIOS BOTH FILE-BASED AND DART
coupling_writer_2D: coupling_writer_2D.F90 gwrite_writer2D.fh
	${COMP}   -g -c -o coupling_writer_2D.o ${ADIOS_INC} coupling_writer_2D.F90
	${LINKER}   -g -o coupling_writer_2D coupling_writer_2D.o ${ADIOS_FLIB} ${DART_FLIB}

coupling_reader_2D: coupling_reader_2D.F90 gread_writer2D.fh gwrite_reader2D.fh
	${COMP}   -g -c -o coupling_reader_2D.o ${ADIOS_INC} coupling_reader_2D.F90
	${LINKER}   -g -o coupling_reader_2D coupling_reader_2D.o ${ADIOS_FLIB} ${ADIOSREAD_FLIB} ${DART_FLIB}

gwrite_writer2D.fh: coupling2D_writer.xml
	${GPP} coupling2D_writer.xml

gwrite_reader2D.fh: coupling2D_reader.xml
	${GPP}  coupling2D_reader.xml

clean:
	rm -f *.o *.mod *.fh core.*
	rm -f coupling_writer_2D coupling_reader_2D 
	rm -f coupling_writer_3D coupling_reader_3D 
	rm -f coupling_writer_2D_dart coupling_reader_2D_dart
	rm -f coupling_writer_3D_dart coupling_reader_3D_dart
	rm -f coupling_writer_2D_combined coupling_reader_2D_combined
        
