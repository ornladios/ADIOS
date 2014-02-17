## Chester, Titan
FC=ftn
CC=cc
## Rhea, Sith
#FC=mpif90
#CC=mpicc

#ADIOS_DIR = /opt/adios/1.6.0
ADIOS_INC = $(shell ${ADIOS_DIR}/bin/adios_config -c)
ADIOS_CLIB = $(shell ${ADIOS_DIR}/bin/adios_config -l)
ADIOS_FLIB = $(shell ${ADIOS_DIR}/bin/adios_config -lf)
ADIOSREAD_SEQ_CLIB = $(shell ${ADIOS_DIR}/bin/adios_config -l -r -s)
GPP = ${ADIOS_DIR}/bin/gpp.py


default: all
all: adios_global_no_xml writer_noxml

writer_noxml: writer_noxml.F90 
	${FC} -g -c -o writer_noxml.o ${ADIOS_INC} writer_noxml.F90
	${FC} -g -o writer_noxml writer_noxml.o ${ADIOS_FLIB}

adios_global_no_xml: adios_global_no_xml.c
	${CC}  -g -o adios_global_no_xml ${ADIOS_INC} adios_global_no_xml.c ${ADIOS_CLIB}

clean:
	rm -f *.o *.mod *.fh core.* fort.* *.png *.bp 
	rm -f adios_global_no_xml writer_noxml

