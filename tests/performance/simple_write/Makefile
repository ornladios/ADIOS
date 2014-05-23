# 
# Load adios and cray-hdf5-parallel modules before compilation
#

####################
## Chester, Titan ##
####################
CC=cc
#H5PCC=cc 
## ADIOS with VAMPIR
#ADIOS_DIR = /ccs/proj/e2e/pnorbert/ADIOS/xk6.pgi.vampir
ADIOS_DIR = /ccs/proj/e2e/pnorbert/ADIOS/chester.gnu.vampir
#ADIOS_INC = -I${ADIOS_DIR}/include
#ADIOS_CLIB = -L${ADIOS_DIR}/lib -ladios -L/sw/xk6/mxml/2.6/cle4.0_pgi11.8.0/lib -lmxml #/usr/lib64/liblustreapi.a
override LINKER=vtcc -vt:mpi -vt:cc cc
## ADIOS vanilla
#ADIOS_DIR=/ccs/proj/e2e/pnorbert/ADIOS/xk6.pgi
ADIOS_INC = $(shell ${ADIOS_DIR}/bin/adios_config -c)
ADIOS_CLIB = $(shell ${ADIOS_DIR}/bin/adios_config -l)
#override LINKER=${CC}

################
## Rhea, Sith ##
################
#CC=mpicc -g 
#H5PCC=h5pcc -g
##ADIOS_DIR = /ccs/proj/e2e/pnorbert/ADIOS/rhea.vampir.pgi
#ADIOS_DIR = /ccs/proj/e2e/pnorbert/ADIOS/sith.vampir.gnu
#ADIOS_INC = $(shell ${ADIOS_DIR}/bin/adios_config -c)
#ADIOS_CLIB = $(shell ${ADIOS_DIR}/bin/adios_config -l)
#override LINKER=vtcc -vt:mpi -vt:cc mpicc
#override LINKER=${CC}

################
## AdiosVM    ##
################
#CC=mpicc -g 
#H5PCC=/opt/hdf5-1.8.12-parallel/bin/h5pcc -g
#ADIOS_DIR = /opt/adios/1.6.0
#ADIOS_INC = $(shell ${ADIOS_DIR}/bin/adios_config -c)
#ADIOS_CLIB = $(shell ${ADIOS_DIR}/bin/adios_config -l)
#override LINKER=${CC}

####################
## EOS with TAU   ##
####################
# call
#     source /ccs/proj/super/TOOLS/tau/tau.bashrc
# before doing make so we have
# TAU=/ccs/proj/super/TOOLS/tau/tau_latest/craycnl/lib
#
#CC=TAU_OPTIONS="-optTrackIO -optVerbose" TAU_MAKEFILE=/${TAU}/Makefile.tau-intel-mpi-pdt  tau_cc.sh 
#H5PCC=cc 
## ADIOS with TAU
#ADIOS_DIR = /ccs/proj/e2e/pnorbert/ADIOS/eos.intel
##ADIOS_DIR = /ccs/proj/e2e/pnorbert/ADIOS/eos_tau.intel
#ADIOS_INC = $(shell ${ADIOS_DIR}/bin/adios_config -c)
#ADIOS_CLIB = $(shell ${ADIOS_DIR}/bin/adios_config -l)
#override LINKER=${CC}

default: all
all: setup_comp_comm_ratio writer_adios writer_posix writer_none writer_hdf5_nto1 writer_hdf5_seq

OBJS = writer.o timing.o comp_comm.o

%.o:   %.c
	${CC} -c -o $@ $<

setup_comp_comm_ratio: comp_comm.o setup_comp_comm_ratio.c
	${CC} -o setup_comp_comm_ratio setup_comp_comm_ratio.c comp_comm.o

writer_adios: ${OBJS} output_adios.c 
	${CC} -c -o output_adios.o ${ADIOS_INC} output_adios.c
	${LINKER} -o writer_adios writer.o timing.o comp_comm.o output_adios.o ${ADIOS_CLIB}

writer_posix: ${OBJS} output_posix.c 
	${CC} -c -o output_posix.o output_posix.c
	${CC} -o writer_posix writer.o timing.o comp_comm.o output_posix.o 

writer_none: ${OBJS} output_none.c 
	${CC} -c -o output_none.o output_none.c
	${CC} -o writer_none writer.o timing.o comp_comm.o output_none.o 

writer_hdf5_nto1: ${OBJS} output_hdf5_nto1.c 
	${H5PCC} -c -o output_hdf5_nto1.o output_hdf5_nto1.c
	${H5PCC} -o writer_hdf5_nto1 writer.o timing.o comp_comm.o output_hdf5_nto1.o 

writer_hdf5_seq: ${OBJS} output_hdf5_seq.c 
	${H5PCC} -c -o output_hdf5_seq.o output_hdf5_seq.c
	${H5PCC} -o writer_hdf5_seq writer.o timing.o comp_comm.o output_hdf5_seq.o 

clean:
	rm -f *.o core*  
	rm -f writer_adios writer_posix writer_none writer_hdf5_nto1 writer_hdf5_seq
	rm -f setup_comp_comm_ratio
	rm -rf data[0-9]*

