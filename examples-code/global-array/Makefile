adios_global: adios_global.o
	mpicc -o adios_global adios_global.o -L../../src  -L/ccs/home/cgj/work/ewok/lib -ladios -lmxml ${HDF5_CLIB}
adios_global.o:adios_global.c
	mpicc -c ${HDF5_CLIB} adios_global.c -I../../src -I/ccs/home/cgj/work/ewok
