JAGUAR=n
EWOK=y
DRAGON=n
CCMP = cc -g
DEBUG=y

ifeq ($(DEBUG),y)
   DOPT=-DDEBUG
else
   DOPT=
endif
ifeq ($(JAGUAR),y)
   ARCH=jaguar
   #NETCDF_CLIB= $(NETCDF_LIB)#-I/usr/lib/netcdf-3.6.1/incude -L/usr/lib/netcdf-3.6.1/lib -lnetcdf
endif
ifeq ($(DRAGON),y)
   ARCH=ewok
   NETCDF_CLIB=-I/usr/lib/netcdf-3.6.1/incude -L/usr/lib/netcdf-3.6.1/lib -lnetcdf
endif
ifeq ($(EWOK),y)
   CCMP = pgcc -g
   ARCH=ewok
   NETCDF_CLIB=${NETCDF_LIB}
endif

INC=-I. -I$(HOME)/work/$(ARCH)/include
LIBDIR=$(HOME)/work/$(ARCH)/lib

all: bp2ncd bp2ncdphi bp2ncdxgc

bp2ncd:  bp2ncd.o Makefile
	$(CCMP) $(INC) $(DOPT) -o bp2ncd bp2ncd.o ${NETCDF_CLIB} -L$(LIBDIR) -ladios -lmxml
bp2ncdphi:  bp2ncdphi.o Makefile
	$(CCMP) $(INC) $(DOPT) -o bp2ncdphi bp2ncdphi.o ${NETCDF_CLIB} -L$(LIBDIR) -ladios -lmxml
bp2ncdxgc:  bp2ncdxgc.o Makefile
	$(CCMP) $(INC) $(DOPT) -o bp2ncdxgc bp2ncdxgc.o ${NETCDF_CLIB} -L$(LIBDIR) -ladios -lmxml

bp2ncdxgc.o: bp2ncdxgc.c
	$(CCMP) $(INC) -c bp2ncdxgc.c ${NETCDF_CLIB}
bp2ncd.o: bp2ncd.c
	$(CCMP) $(INC) -c bp2ncd.c ${NETCDF_CLIB}
bp2ncdphi.o: bp2ncdphi.c
	$(CCMP) $(INC) -c bp2ncdphi.c ${NETCDF_CLIB}

clean:
	rm -f *.o bp2ncd
