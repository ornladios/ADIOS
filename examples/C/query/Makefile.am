AM_CPPFLAGS = $(all_includes) -I$(top_builddir)/src/public  -I$(top_srcdir)/src/public -I$(top_srcdir)/src/query

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = write_vars query_vars write_table query_table

write_vars_SOURCES = write_vars.c
write_vars_LDADD = $(top_builddir)/src/libadios_nompi.a 
write_vars_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS)
write_vars_LDADD += $(ADIOSLIB_SEQ_LDADD)

query_vars_SOURCES = query_vars.c
query_vars_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
query_vars_LDFLAGS = $(ADIOSREADLIB_SEQ_LDFLAGS)
query_vars_LDADD += $(ADIOSREADLIB_SEQ_LDADD)

write_table_SOURCES = write_table.c
write_table_LDADD = $(top_builddir)/src/libadios_nompi.a 
write_table_LDFLAGS = $(ADIOSLIB_SEQ_LDFLAGS)
write_table_LDADD += $(ADIOSLIB_SEQ_LDADD)

query_table_SOURCES = query_table.c
query_table_LDADD = $(top_builddir)/src/libadiosread_nompi.a 
query_table_LDFLAGS = $(ADIOSREADLIB_SEQ_LDFLAGS)
query_table_LDADD += $(ADIOSREADLIB_SEQ_LDADD)

EXTRA_DIST = README_vars.txt README_table.txt

CLEANFILES = *.bp *.idx

CC=$(MPICC)


