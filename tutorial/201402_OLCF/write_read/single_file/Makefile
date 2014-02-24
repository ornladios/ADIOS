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
all: writer_single_file reader_all_at_once reader_step_by_step

writer_single_file: writer_single_file.F90 gwrite_writer.fh
	${FC} -g -o writer_single_file ${ADIOS_INC} writer_single_file.F90 ${ADIOS_FLIB}

gwrite_writer.fh: writer.xml
	${GPP} writer.xml

reader_all_at_once: reader_all_at_once.F90 
	${FC}  -g -o reader_all_at_once ${ADIOS_INC} reader_all_at_once.F90 ${ADIOSREAD_FLIB} 

reader_step_by_step: reader_step_by_step.F90 
	${FC}  -g -o reader_step_by_step ${ADIOS_INC} reader_step_by_step.F90 ${ADIOSREAD_FLIB} 

clean:
	rm -f *.o *.mod *.fh core.* fort.*
	rm -f writer*.bp writer*.h5 writer*.nc
	rm -f writer_single_file reader_all_at_once reader_step_by_step

