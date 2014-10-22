HOST=titan

# Do not modify below. Add your settings in make.settings.XXXX

include make.settings.$(HOST)

default: makeobjdirs libadios.a libadiosf.a bpls bpdump
all: default
include Makefile.srcdef

INC:=$(INCLUDES) $(MXML_INC) -I.

libadios_nompi_a_CSRC := $(filter %.c,$(libadios_nompi_a_SOURCES))
libadios_nompi_a_FSRC := $(filter %.f90,$(libadios_nompi_a_SOURCES))
libadios_nompi_a_COBJS := $(libadios_nompi_a_CSRC:.c=.o)
libadios_nompi_a_FOBJS := $(libadios_nompi_a_FSRC:.f90=.o)

libadiosf_nompi_a_CSRC := $(filter %.c,$(libadiosf_nompi_a_SOURCES))
libadiosf_nompi_a_FSRC := $(filter %.f90,$(libadiosf_nompi_a_SOURCES))
libadiosf_nompi_a_COBJS := $(libadiosf_nompi_a_CSRC:.c=.o)
libadiosf_nompi_a_FOBJS := $(libadiosf_nompi_a_FSRC:.f90=.o)


OBJDIRS = core core/transforms transforms read write 

help:
	@echo "Settings imported from: make.settings.$(HOST)"
	@echo "CC=$(CC)"
	@echo "CPPFLAGS=$(CPPFLAGS)"
	@echo "CFLAGS=$(CFLAGS)"
	@echo "FC=$(FC)"
	@echo "FFLAGS=$(FFLAGS)"
	@echo "INC=$(INC)"
	@echo " "
	@echo "Install path PREFIX=$(PREFIX)"
	@echo " "
	@echo "Build steps: "
	@echo "   make "
	@echo "   make test"
	@echo "   ./bpls -la write_read_1.bp"
	@echo "   ./bpdump -la write_read_1.bp"
	@echo "   make install"
	@echo " "
	@echo "Application building flags" 
	@echo "   CFLAGS=-I$(PREFIX)/include"
	@echo "   LDFLAGS=-L$(PREFIX)/lib $(MXML_LDFLAGS)"
	@echo "   CLIBS=-ladios  $(MXML_LIBS) -lm"
	@echo "   FLIBS=-ladiosf $(MXML_LIBS) -lm"



%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -D_NOMPI $(INC) $< -o $@


%.o: %.f90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c -D_NOMPI $(INC) $< -o $@


makeobjdirs:
	@echo "================================================================="
	@echo "Create subdirectories to store object files"
	mkdir -p $(OBJDIRS)
	

libadios.a: $(libadios_nompi_a_COBJS) $(libadios_nompi_a_FOBJS)
	@echo "================================================================="
	@echo "Create libadios.a"
	ar rc libadios.a $(libadios_nompi_a_COBJS) $(libadios_nompi_a_FOBJS)


libadiosf.a: $(libadiosf_nompi_a_COBJS) $(libadiosf_nompi_a_FOBJS)
	@echo "================================================================="
	@echo "Create libadiosf.a"
	ar rc libadiosf.a $(libadiosf_nompi_a_COBJS) $(libadiosf_nompi_a_FOBJS)


bpls:   libadios.a
	@echo "================================================================="
	@echo "Build bpls"
	$(CC)  $(CPPFLAGS) \
	$(CFLAGS) -D_NOMPI \
	-I../src/public -I../utils/bpls \
	../utils/bpls/bpls.c -o bpls \
	-L. $(MXML_LDFLAGS) -ladios $(MXML_LIBS) -lm


bpdump: libadios.a
	@echo "================================================================="
	@echo "Build bpdump"
	$(CC)  $(CPPFLAGS) \
	$(CFLAGS) -D_NOMPI \
	$(INCLUDES) -I../src/public -I. \
	../utils/bpdump/bpdump.c -o bpdump \
	-L. $(MXML_LDFLAGS) -ladios $(MXML_LIBS) -lm


write_read.xml: ../tests/suite/programs/write_read.xml
	cat ../tests/suite/programs/write_read.xml | sed -s "s/method.*group=\"alltypes\".*method=\".*\"/method group=\"alltypes\" method=\"POSIX1\"/" > write_read.xml

write_read: libadios.a ../tests/suite/programs/write_read.c write_read.xml
	$(CC)  $(CPPFLAGS) \
	$(CFLAGS) -D_NOMPI \
	-I../src/public \
	../tests/suite/programs/write_read.c -o write_read \
	-L. $(MXML_LDFLAGS) -ladios $(MXML_LIBS) -lm

test:   write_read
	./write_read
	

check:  test

install: libadios.a libadiosf.a bpls bpdump
	mkdir -p $(PREFIX) $(PREFIX)/include $(PREFIX)/lib 
	(cd ../src; $(INSTALLCMD) $(include_HEADERS) $(PREFIX)/include)
	$(INSTALLCMD) adios_*.mod $(PREFIX)/include
	$(INSTALLCMD) libadios.a libadiosf.a $(PREFIX)/lib


clean: 
	rm -f bpls
	rm -f bpdump
	rm -f libadios.a 
	rm -f libadiosf.a
	rm -f adios_*.mod
	rm -f $(libadiosf_nompi_a_COBJS) 
	rm -f $(libadiosf_nompi_a_FOBJS)
	rm -rf $(OBJDIRS)
	rm -rf write_read*

