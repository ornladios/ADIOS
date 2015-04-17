include ../../adios_config.flags

MEXOPTS=-g -DDEBUG
default:
	@echo "ADIOS reader methods for MATLAB"
	@echo "Have 'mex' compiler in the path and run"
	@echo "> make matlab"


matlab: adiosopenc adiosclosec adiosreadc

INC:=$(shell echo ${ADIOSREAD_SEQ_INC})
LIB:=$(shell echo ${ADIOSREAD_SEQ_CLIB})

adiosopenc:  
	mex adiosopenc.c  ${MEXOPTS} -D_NOMPI -DADIOS_USE_READ_API_1 -I. -I../../src/public ${INC} -L. -L../../src ${LIB}

adiosreadc: 
	mex adiosreadc.c  ${MEXOPTS} -D_NOMPI -DADIOS_USE_READ_API_1 -I. -I../../src/public ${INC} -L. -L../../src ${LIB}

adiosclosec:
	mex adiosclosec.c ${MEXOPTS} -D_NOMPI -DADIOS_USE_READ_API_1 -I. -I../../src/public ${INC} -L. -L../../src ${LIB}

#adiosdbgc:
#       mex adiosdbgc.c ${MEXOPTS} -D_NOMPI -I. -I../../src/public ${INC} -L. -L../../src ${LIB} 


clean:
	rm -f *.mexa64 core* *.o 

install:
	@echo "Install *.m and *.mexa64 into MATLABPATH"
	@echo "Do not have 'adios' anywhere in the path where you install the files"
	@echo "otherwise 'help adios' will not work in Matlab"
