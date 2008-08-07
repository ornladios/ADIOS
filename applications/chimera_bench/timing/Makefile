#use this for ewok
EWOK=y
JAGUAR=n

ifeq ($(JAGUAR),y)
    #use this for cray compute nodes
    CC=cc -c -g -target=linux

    ADIOS_LIB_HOME=$(HOME)/work/jaguar
    ADIOS_LIB_I=$(ADIOS_LIB_HOME)/include

    I_PATH=-I$(ADIOS_LIB_I)

    INSTALL_PREFIX=$(HOME)/work/jaguar
endif

ifeq ($(EWOK),y)
    CC=mpicc -c -g
    ADIOS_LIB_HOME=$(HOME)/work/ewok
    ADIOS_LIB_I=$(ADIOS_LIB_HOME)/include

    I_PATH=-I$(ADIOS_LIB_I)

    INSTALL_PREFIX=$(HOME)/work/ewok
endif

all: libadios_timing.a

libadios_timing.a: timing.o
	rm -f libadios_timing.a
	ar crvs libadios_timing.a timing.o
	ranlib libadios_timing.a

.c.o :
	$(CC) $(I_PATH) $<

clean:
	rm -f *.o libadios_timing.a

install:
	cp adios-timing.h $(INSTALL_PREFIX)/include
	cp libadios_timing.a $(INSTALL_PREFIX)/lib

