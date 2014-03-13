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
all: nc4_write nc4_read adios_write adios_read

nc4_write: nc4_write.F90 
	${FC}   -g -o nc4_write -I${NETCDF_DIR}/include  nc4_write.F90 -L${NETCDF_DIR}/lib/libnetcdff

nc4_read: nc4_read.F90 
	${FC}   -g -o nc4_read -I${NETCDF_DIR}/include  nc4_read.F90 -L${NETCDF_DIR}/lib/libnetcdff

adios_write: adios_write.F90 gwrite_writer.fh
	${FC}   -g -o adios_write ${ADIOS_INC} adios_write.F90 ${ADIOS_FLIB}

gwrite_writer.fh: adios_write.xml
	${GPP} adios_write.xml

adios_read: adios_read.F90 
	${FC}   -g -o adios_read ${ADIOS_INC} adios_read.F90 ${ADIOSREAD_FLIB}

clean:
	rm -f *.o *.mod *.fh core.* fort.*
	rm -f *.nc *.bp
	rm -f nc4_write adios_write
	rm -f nc4_read adios_read


