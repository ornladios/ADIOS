AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src -I$(top_builddir)/src/public  -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp -p $(srcdir)/*.xml $(builddir)

check_PROGRAMS = local_array_write_C local_array_read_C

local_array_write_C_SOURCES = local_array_write_C.c
local_array_write_C_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
local_array_write_C_LDFLAGS = $(ADIOSLIB_LDFLAGS)
local_array_write_C.o: gwrite_arrays.ch 
gwrite_arrays.ch: local_array_C.xml
	$(top_builddir)/utils/gpp/gpp.py $(srcdir)/local_array_C.xml

local_array_read_C_SOURCES = local_array_read_C.c
local_array_read_C_LDADD = $(top_builddir)/src/libadiosread.a $(ADIOSREADLIB_LDADD)
local_array_read_C_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)

CLEANFILES = *.bp *.ch *.fh
CC=$(MPICC)
EXTRA_DIST = local_array_C.xml 
