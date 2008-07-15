INCLUDES = $(all_includes)
INCLUDES += -I$(top_builddir)/src

AUTOMAKE_OPTIONS = no-dependencies

bin_PROGRAMS = bp2ascii

bp2ascii_SOURCES = bp2ascii.c \
                 hw-utils.c
bp2ascii_LDFLAGS = $(MXML_LDFLAGS)
bp2ascii_CPPFLAGS = $(MXML_CPPFLAGS)
bp2ascii_LDADD = $(top_builddir)/src/libadios.a
bp2ascii_LDADD += $(top_builddir)/src/libadios_int.a
bp2ascii_LDADD += $(MXML_LIBS)
