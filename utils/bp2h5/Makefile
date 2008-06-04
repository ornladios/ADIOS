EWOK=n
JAGUAR=y

ifeq ($(JAGUAR),y)
   CCMP = cc -g
   ARCH=jaguar
   INCDIR1=
   LIBDIR1=
endif

ifeq ($(EWOK),y)
   CCMP = mpicc -g
   ARCH=ewok
   INCDIR1=
   LIBDIR1=
   HDF5_CLIB=-L/apps/HDF5/1.6.5/ewok_par/pgi625/lib -I/apps/HDF5/1.6.5/ewok_par/pgi625/include -lhdf5 -lz
endif

INCDIR=$(INCDIR1) -I. -I$(HOME)/work/$(ARCH)/include
LIBDIR=$(LIBDIR1) -L$(HOME)/work/$(ARCH)/lib
LIBS=-ladios -lmxml

prog:  hw-utils.c bp2h5.c
	$(CCMP) $(INCDIR) $(LIBDIR) -o bp2h5 hw-utils.c bp2h5.c $(HDF5_CLIB) $(LIBS)

clean:
	rm -f *.o bp2h5
