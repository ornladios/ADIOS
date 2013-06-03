## MPIF90 SETTINGS
FC=mpif90
FFLAGS=-g

## FTN SETTINGS
#FC=ftn
#FFLAGS = -O3

#BG/P IBM compiler
#FC=mpixlf90
#FFLAGS=-O0 -qarch=450 -qtune=450

#TARGET=gnu
TARGET=pgi

## Set ADIOS_DIR here or before doing make
override ADIOS_DIR:=/ccs/proj/e2e/pnorbert/ADIOS/sith.pgi
override ADIOS_INC:=` ${ADIOS_DIR}/bin/adios_config -c -f`
override ADIOS_FLIB:=`${ADIOS_DIR}/bin/adios_config -l -f`
override GPP = ${ADIOS_DIR}/bin/gpp.py

default: coupling_writer_2D coupling_reader_2D 
all: default

coupling_writer_2D.o : coupling_writer_2D.F90 gwrite_writer2D.fh
	${FC} -g -c ${ADIOS_INC} coupling_writer_2D.F90  $<

coupling_writer_2D:  coupling_writer_2D.o gwrite_writer2D.fh
	${FC} -g -o coupling_writer_2D coupling_writer_2D.o ${ADIOS_FLIB} 

gwrite_writer2D.fh: coupling_writer_2D.xml
	${GPP} coupling_writer_2D.xml

coupling_reader_2D.o : coupling_reader_2D.F90 gwrite_reader2D.fh
	${FC} -g -c ${ADIOS_INC} coupling_reader_2D.F90  $<

coupling_reader_2D:  coupling_reader_2D.o gwrite_reader2D.fh
	${FC} -g -o coupling_reader_2D coupling_reader_2D.o ${ADIOS_FLIB} 

gwrite_reader2D.fh: coupling_reader_2D.xml
	${GPP} coupling_reader_2D.xml


clean:
	rm -f *.o *.mod *.fh core.* 
	rm -f coupling_writer_2D coupling_reader_2D 
        
distclean: clean
	rm -f log.* dataspaces.conf conf srv.lck *.bp
