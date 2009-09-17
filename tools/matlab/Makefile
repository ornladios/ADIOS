
default:
	@echo "ADIOS reader methods for MATLAB"
	@echo "Have 'mex' compiler in the path and run"
	@echo "> make matlab"


matlab: adiosopenc adiosclosec adiosreadc 

adiosopenc:  
	mex adiosopenc.c -DNOMPI  -I. -I../../src -L. -L../../src -ladiosread_nompi

adiosreadc: 
	mex adiosreadc.c -DNOMPI -I. -I../../src -L. -L../../src -ladiosread_nompi

adiosclosec:
	mex adiosclosec.c -DNOMPI -I. -I../../src -L. -L../../src -ladiosread_nompi

clean:
	rm -f *.mexa64 core* *.o 

install:
	@echo "Install *.m and *.mexa64 into MATLABPATH"
	@echo "Do not have 'adios' anywhere in the path where you install the files"
	@echo "otherwise 'help adios' will not work in Matlab"
