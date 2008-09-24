adios_timeindex: adios_timeindex.o
	mpicc -o adios_timeindex adios_timeindex.o -L../../src  -L/ccs/home/cgj/work/ewok/lib -ladios -lmxml ${HDF5_CLIB}
adios_timeindex.o:adios_timeindex.c
	mpicc -c ${HDF5_CLIB} adios_timeindex.c -I../../src -I/ccs/home/cgj/work/ewok
