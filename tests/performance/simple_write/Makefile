# 
# Load adios and cray-hdf5-parallel modules before compilation
#

####################
## Chester, Titan ##
####################
CC=cc
H5PCC=cc 
## ADIOS with VAMPIR
#ADIOS_DIR = /ccs/proj/e2e/pnorbert/ADIOS/xk6.pgi.vampir
#ADIOS_INC = -I${ADIOS_DIR}/include
#ADIOS_CLIB = -L${ADIOS_DIR}/lib -ladios -L/sw/xk6/mxml/2.6/cle4.0_pgi11.8.0/lib -lmxml /usr/lib64/liblustreapi.a
#override LINKER=vtcc -vt:mpi -vt:cc cc
## ADIOS vanilla
ADIOS_DIR=/ccs/proj/e2e/pnorbert/ADIOS/xk6.pgi
ADIOS_INC = $(shell ${ADIOS_DIR}/bin/adios_config -c)
ADIOS_CLIB = $(shell ${ADIOS_DIR}/bin/adios_config -l)
override LINKER=cc

################
## Rhea, Sith ##
################
#CC=mpicc -g 
#H5PCC=h5pcc -g
#ADIOS_DIR = /ccs/proj/e2e/pnorbert/ADIOS/rhea.vampir.pgi
#ADIOS_INC = $(shell ${ADIOS_DIR}/bin/adios_config -c)
#ADIOS_CLIB = $(shell ${ADIOS_DIR}/bin/adios_config -l)
#override LINKER=vtcc -vt:mpi -vt:cc mpicc


default: all
all: writer_adios writer_posix writer_none writer_hdf5_nto1 writer_hdf5_seq

%.o:   %.c
	${CC} -c -o $@ $<

writer.o:   writer.c
	${CC} -c -o writer.o writer.c 

writer_adios: writer.o timing.o output_adios.c 
	${CC} -c -o output_adios.o ${ADIOS_INC} output_adios.c
	${LINKER} -o writer_adios writer.o timing.o output_adios.o ${ADIOS_CLIB}

writer_posix: writer.o timing.o output_posix.c 
	${CC} -c -o output_posix.o output_posix.c
	${CC} -o writer_posix writer.o timing.o output_posix.o 

writer_none: writer.o timing.o output_none.c 
	${CC} -c -o output_none.o output_none.c
	${CC} -o writer_none writer.o timing.o output_none.o 

writer_hdf5_nto1: writer.o timing.o output_hdf5_nto1.c 
	${H5PCC} -c -o output_hdf5_nto1.o output_hdf5_nto1.c
	${H5PCC} -o writer_hdf5_nto1 writer.o timing.o output_hdf5_nto1.o 

writer_hdf5_seq: writer.o timing.o output_hdf5_seq.c 
	${H5PCC} -c -o output_hdf5_seq.o output_hdf5_seq.c
	${H5PCC} -o writer_hdf5_seq writer.o timing.o output_hdf5_seq.o 

clean:
	rm -f *.o core*  
	rm -f writer_adios writer_posix writer_none writer_hdf5_nto1 writer_hdf5_seq
	rm -rf data[0-9]*

