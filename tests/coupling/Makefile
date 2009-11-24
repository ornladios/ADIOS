MXML_LIB = -L/ccs/proj/e2e/pnorbert/mxml/lib -lmxml
default: coupling2D
all: coupling2D

coupling2D: coupling_writer_2D coupling_reader_2D

coupling_writer_2D: coupling_writer_2D.F90 gwrite_writer2D.fh
	mpif90 -g -c -o coupling_writer_2D.o -I../../src coupling_writer_2D.F90
	mpif90 -g -o coupling_writer_2D coupling_writer_2D.o -L../../src -ladiosf ${MXML_LIB}

coupling_reader_2D: coupling_reader_2D.F90 gread_writer2D.fh gwrite_reader2D.fh
	mpif90 -g -c -o coupling_reader_2D.o -I../../src coupling_reader_2D.F90
	mpif90 -g -o coupling_reader_2D coupling_reader_2D.o -L../../src -ladiosf -ladiosreadf ${MXML_LIB}

gwrite_writer2D.fh: coupling2D.xml
	../../utils/gpp/gpp.py coupling2D.xml

gread_writer2D.fh: coupling2D.xml
	../../utils/gpp/gpp.py coupling2D.xml

gwrite_reader2D.fh: coupling2D.xml
	../../utils/gpp/gpp.py coupling2D.xml

clean:
	rm -f coupling_writer_2D coupling_reader_2D *.o *.mod *.fh core.*
        
