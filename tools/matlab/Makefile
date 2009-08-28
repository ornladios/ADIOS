
default:
	@echo "ADIOS reader methods for MATLAB"
	@echo "Have 'mex' compiler in the path and run"
	@echo "> make matlab"


matlab: adiosopenc adiosreadc adiosclosec

adiosopenc:  adios_readutil
	mex adiosopenc.c  -I. -I../../src -L. -ladios_readutil -L../../src -ladiosread_nompi

adiosreadc:  adios_readutil
	mex adiosreadc.c  -I. -I../../src -L. -ladios_readutil -L../../src -ladiosread_nompi

adiosclosec:  adios_readutil
	mex adiosclosec.c -I. -I../../src -L. -ladios_readutil -L../../src -ladiosread_nompi

adiosreadc:  adios_readutil
adios_readutil:
	gcc -c -fPIC -DNOMPI -I../../src adios_readutil.c
	rm -rf libadios_readutil.a
	ar cru libadios_readutil.a adios_readutil.o
	ranlib libadios_readutil.a
	rm -rf adios_readutil.o

clean:
	rm -f libadios_readutil.a *.mexa64 core* *.o 

install:
	@echo "Install *.m and *.mexa64 into MATLABPATH"
