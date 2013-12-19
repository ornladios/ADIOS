INCLUDES = $(all_includes)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/public

AUTOMAKE_OPTIONS = no-dependencies

all-local:
	test "$(srcdir)" = "$(builddir)" || cp $(srcdir)/*.sh $(builddir)
	test "$(srcdir)" = "$(builddir)" || cp -r $(srcdir)/conf_genarray3D $(builddir)
	test "$(srcdir)" = "$(builddir)" || cp -r $(srcdir)/conf_simple $(builddir)

noinst_PROGRAMS = adios_global adios_read_box adios_read_points genarray3D 

adios_global_SOURCES = adios_global.c
adios_global_LDADD = $(top_builddir)/src/libadios.a 
adios_global_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_LDADD += $(ADIOSLIB_LDADD)
adios_read_box_SOURCES = adios_read_box.c
adios_read_box_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_box_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_box_LDADD += $(ADIOSLIB_LDADD)
adios_read_points_SOURCES = adios_read_points.c
adios_read_points_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_points_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_points_LDADD += $(ADIOSLIB_LDADD)
genarray3D_SOURCES = genarray3D.c
genarray3D_LDADD = $(top_builddir)/src/libadios.a 
genarray3D_LDFLAGS = $(ADIOSLIB_LDFLAGS)
genarray3D_LDADD += $(ADIOSLIB_LDADD)


CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = gread_temperature.ch gwrite_temperature.ch \
             read_test.h \
	     test_read_box.sh  test_read_points.sh  write_all.sh \
	     conf_genarray3D conf_simple 

