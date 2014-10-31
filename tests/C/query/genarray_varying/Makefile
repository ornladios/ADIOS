
# Change setup in make.settings

default: genarray_varying tests
all: default

tests:
	make -f Makefile.tests

genarray_varying: genarray_varying.F90 test_query.F90
	make -f Makefile.genarray_varying

clean:
	make -f Makefile.tests clean
	make -f Makefile.genarray_varying clean

distclean: 
	make -f Makefile.tests distclean
	make -f Makefile.genarray_varying distclean
	rm -f *.png

