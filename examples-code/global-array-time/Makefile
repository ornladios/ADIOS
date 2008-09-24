adios_globaltime: adios_globaltime.o
	mpicc -o adios_globaltime adios_globaltime.o -L../../src  -L/ccs/home/cgj/work/ewok/lib -ladios -lmxml ${HDF5_CLIB}
adios_globaltime.o:adios_globaltime.c
	mpicc -c ${HDF5_CLIB} adios_globaltime.c -I../../src -I/ccs/home/cgj/work/ewok
