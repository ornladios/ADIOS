
# Change setup in make.settings

default: genarray_stream stage_write 
all: default

stage_write: stage_write.c
	make -f Makefile.stage_write 

genarray_stream: genarray_stream.F90
	make -f Makefile.genarray_stream

clean:
	make -f Makefile.stage_write clean
	make -f Makefile.genarray_stream clean
	rm -f staging.o[0-9]*

distclean: 
	make -f Makefile.stage_write distclean
	make -f Makefile.genarray_stream distclean
	rm -f staging.o[0-9]*

