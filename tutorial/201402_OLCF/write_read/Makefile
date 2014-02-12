## Chester, Titan
FC=ftn
## Rhea, Sith
#FC=mpif90

LINKER=${FC}
#ADIOS_DIR = /opt/adios/1.6.0
ADIOS_INC := $(shell ${ADIOS_DIR}/bin/adios_config -c -f)
ADIOS_FLIB := $(shell ${ADIOS_DIR}/bin/adios_config -l -f)
ADIOSREAD_FLIB := $(shell ${ADIOS_DIR}/bin/adios_config -l -f -r)
GPP = ${ADIOS_DIR}/bin/gpp.py


default: all
all: writer reader

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

