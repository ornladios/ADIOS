
CC=$(MPICC)


lib_LIBRARIES = libskel.a
libskel_a_SOURCES = skel_xml_output.c

if BUILD_FORTRAN
libskel_a_SOURCES += skel_xml_output_f.c
endif

library_includedir=$(includedir)/skel
library_include_HEADERS = skel_xml_output.h

