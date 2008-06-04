CCMP = cc -g
EWOK=n
JAGUAR=y

ifeq ($(JAGUAR),y)
   ARCH=jaguar
endif
ifeq ($(EWOK),y)
   ARCH=ewok
endif

INC=-I. -I$(HOME)/work/$(ARCH)/include
LIBDIR=$(HOME)/work/$(ARCH)/lib

bpdump:  bpdump.o Makefile
	$(CCMP) $(INC) -o bpdump bpdump.o -L$(LIBDIR) -ladios -lmxml

bpdump.o: bpdump.c
	$(CCMP) $(INC) -c bpdump.c

clean:
	rm -f *.o bpdump
