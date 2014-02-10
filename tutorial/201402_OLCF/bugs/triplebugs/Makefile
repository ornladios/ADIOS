## Chester, Titan
FC=ftn
## Rhea, Sith
#FC=mpif90

LINKER=${FC}
#ADIOS_DIR = /opt/adios/1.6.0
ADIOS_INC := -I ${ADIOS_DIR}/include
ADIOS_FLIB := -L ${ADIOS_DIR}/lib -ladiosf -lrdmacm -libverbs -llustreapi -lz 


default: all
all: write

write: write.F90 
	${FC}   -g -c -o write.o ${ADIOS_INC} write.F90
	${LINKER}   -g -o write write.o ${ADIOS_FLIB}

clean:
	rm -f *.o *.mod *.fh core.* fort.*
	rm -f write write*.bp 

