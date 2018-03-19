AM_CPPFLAGS = $(all_includes)
AM_CPPFLAGS += -I$(top_builddir)/src/public -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

#all-local:
#	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)
	
noinst_PROGRAMS = strings_write 
#strings_read 

strings_write_SOURCES = strings_write.c
strings_write_LDADD = $(top_builddir)/src/libadios.a 
strings_write_LDFLAGS = $(ADIOSLIB_LDFLAGS) $(ADIOSLIB_EXTRA_LDFLAGS)
strings_write_LDADD += $(ADIOSLIB_LDADD)
#strings_read_SOURCES = strings_read.c
#strings_read_LDADD = $(top_builddir)/src/libadiosread.a 
#strings_read_LDFLAGS = $(ADIOSREADLIB_LDFLAGS)
#strings_read_LDADD += $(ADIOSREADLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
