INCLUDES = $(all_includes)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.xml $(builddir)

noinst_PROGRAMS=posix_fgr

posix_fgr_SOURCES = posix_fgr.c
posix_fgr_CPPFLAGS = $(ADIOSLIB_CPPFLAGS)
posix_fgr_CFLAGS = $(ADIOSLIB_CFLAGS)
posix_fgr_LDADD = $(top_builddir)/src/libadios.a $(ADIOSLIB_LDADD)
posix_fgr_LDFLAGS = $(ADIOSLIB_LDFLAGS) 

CLEANFILES = posix_fgr.bp
CC=$(MPICC)

EXTRA_DIST = posix_fgr.xml
