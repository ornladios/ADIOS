CCMP = mpicc -g
EWOK=y
JAGUAR=n

ifeq  ($(JAGUAR),y)
  ARCH=jaguar
endif
ifeq ($(EWOK),y)
  ARCH=ewok
endif
#WORK=/apps/adios/
INC=/autofs/na1_home/ge1/work/$(ARCH)/include
LIB=/autofs/na1_home/ge1/work/$(ARCH)/lib
HDF5_CLIB= -I/apps/HDF5/1.6.5/ewok_par/pgi625/include -L/apps/HDF5/1.6.5/ewok_par/pgi625/lib  -lhdf5  -lz
prog:  hw-utils.c bp2ascii.c
	$(CCMP)  $(HDF5_CLIB) -o bp2ascii -L$(LIB) -I$(INC) -ladios hw-utils.c bp2ascii.c

clean:
	rm -f *.o bp2h5
