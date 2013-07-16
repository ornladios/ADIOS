SHELL = /bin/sh

# clear out the suffix list
.SUFFIXES:
# define new suffixes
.SUFFIXES: .c .o
CC=gcc

ADIOS_ROOT=/rock/opt/adios/git-dbg
MXML_ROOT=/rock/opt/mxml/2.7
MPI_ROOT=/rock/opt/openmpi/1.6.3
EVPATH_ROOT=/rock/opt/evpath

#LIBS=`/rock/opt/adios/1.5.0-dbg/bin/adios_config -l`
LIBS_PATHS=-L$(ADIOS_ROOT)/lib -L$(MXML_ROOT)/lib -L$(MPI_ROOT)/lib -L$(EVPATH_ROOT)/lib
LIBS=-ladios -ladiosread -lm -lmpi -levpath -latl -lcercs_env -ldill -lffs -lgen_thread  -lmxml  
#INCLUDES=`/rock/opt/adios/1.5.0-dbg/bin/adios_config -c`
INCLUDES=-I$(ADIOS_ROOT)/include -I$(MXML_ROOT)/include -I$(MPI_ROOT)/include -I$(EVPATH_ROOT)/include

# $^ automatic var - source file whatever it is
# $@ output automatic var

all: reader writer

reader: reader.o  utils.o
	$(CC) -g -o $@ $^ $(INCLUDES) $(LIBS_PATHS) $(LIBS) 
	 
writer: writer.o  utils.o
	$(CC) -g -o $@ $^ $(INCLUDES) $(LIBS_PATHS) $(LIBS) 

.c.o: utils.h misc.h
	$(CC) -g -c $< -o $@ $(INCLUDES) $(LIBS_PATHS) $(LIBS)
	
	
run-w:
	mpirun -np 2 ./writer
	
run-r:
	mpirun -np 2 ./reader
	
clean_test:
	rm 	-rf *_writer_info.txt *_reader_info.txt test.bp *_reader_ready.txt *_writer_ready.txt
	
clean:
	rm -rf reader writer  *.o 
