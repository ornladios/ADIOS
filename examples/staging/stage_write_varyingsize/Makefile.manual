
# Change setup in make.settings

default: genarray_varying stage_write 
all: default

stage_write: stage_write.c
	make -f Makefile.stage_write 

genarray_varying: genarray_varying.F90
	make -f Makefile.genarray_varying

clean:
	make -f Makefile.stage_write clean
	make -f Makefile.genarray_varying clean
	rm -f staging.o[0-9]*

distclean: 
	make -f Makefile.stage_write distclean
	make -f Makefile.genarray_varying distclean
	rm -f staging.o[0-9]*

