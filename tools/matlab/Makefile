include ../../adios_config.flags

MEXOPTS="-g -DDEBUG"
default:
	@echo "ADIOS reader methods for MATLAB"
	@echo "Have 'mex' compiler in the path and run"
	@echo "> make matlab"


matlab: adiosopenc adiosclosec adiosreadc

adiosopenc:  
	mex adiosopenc.c ${MEXOPTS} -D_NOMPI  -I. -I../../src -L. -L../../src ${ADIOSREAD_SEQ_INC} ${ADIOSREAD_SEQ_CLIB}

adiosreadc: 
	mex adiosreadc.c ${MEXOPTS} -D_NOMPI -I. -I../../src -L. -L../../src ${ADIOSREAD_SEQ_INC} ${ADIOSREAD_SEQ_CLIB}

adiosclosec:
	mex adiosclosec.c ${MEXOPTS} -D_NOMPI -I. -I../../src -L. -L../../src ${ADIOSREAD_SEQ_INC} ${ADIOSREAD_SEQ_CLIB}

#adiosdbgc:
#	mex adiosdbgc.c ${MEXOPTS} -D_NOMPI -I. -I../../src -L. -L../../src ${ADIOSREAD_SEQ_INC} ${ADIOSREAD_SEQ_CLIB}

clean:
	rm -f *.mexa64 core* *.o 

install:
	@echo "Install *.m and *.mexa64 into MATLABPATH"
	@echo "Do not have 'adios' anywhere in the path where you install the files"
	@echo "otherwise 'help adios' will not work in Matlab"
