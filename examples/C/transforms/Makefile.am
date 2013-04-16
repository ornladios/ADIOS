INCLUDES = $(all_includes)
INCLUDES += -I$(top_srcdir)/src -I$(top_srcdir)/src/public $(TIMER_CPPFLAGS)

AUTOMAKE_OPTIONS = no-dependencies

noinst_PROGRAMS = adios_global adios_read_subv adios_read_points adios_read_all_3D adios_write_all_3D

adios_global_SOURCES = adios_global.c
adios_global_LDADD = $(top_builddir)/src/libadios.a 
adios_global_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_global_LDADD += $(ADIOSLIB_LDADD)
adios_read_subv_SOURCES = adios_read_subv.c
adios_read_subv_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_subv_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_subv_LDADD += $(ADIOSLIB_LDADD)
adios_read_points_SOURCES = adios_read_points.c
adios_read_points_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_points_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_points_LDADD += $(ADIOSLIB_LDADD)
adios_read_all_3D_SOURCES = adios_read_all_3D.c
adios_read_all_3D_LDADD = $(top_builddir)/src/libadiosread.a 
adios_read_all_3D_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_read_all_3D_LDADD += $(ADIOSLIB_LDADD)
adios_write_all_3D_SOURCES = adios_write_all_3D.c
adios_write_all_3D_LDADD = $(top_builddir)/src/libadios.a 
adios_write_all_3D_LDFLAGS = $(ADIOSLIB_LDFLAGS)
adios_write_all_3D_LDADD += $(ADIOSLIB_LDADD)

CLEANFILES = *.bp
CC=$(MPICC)
EXTRA_DIST = adios_global.xml gread_temperature.ch gwrite_temperature.ch job.pbs

 
