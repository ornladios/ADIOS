MXML_SRCDIR=mxml-2.9

AM_CPPFLAGS=-I$(top_srcdir)/mxml/${MXML_SRCDIR}

AUTOMAKE_OPTIONS = no-dependencies subdir-objects
METASOURCES=AUTO

CLEANFILES= 


noinst_LIBRARIES = libmxml.a
libmxml_a_SOURCES = ${MXML_SRCDIR}/mxml-attr.c     \
                     ${MXML_SRCDIR}/mxml-get.c      \
                     ${MXML_SRCDIR}/mxml-private.c  \
                     ${MXML_SRCDIR}/mxml-string.c   \
                     ${MXML_SRCDIR}/mxml-entity.c   \
                     ${MXML_SRCDIR}/mxml-index.c    \
                     ${MXML_SRCDIR}/mxml-search.c   \
                     ${MXML_SRCDIR}/mxml-file.c     \
                     ${MXML_SRCDIR}/mxml-node.c     \
                     ${MXML_SRCDIR}/mxml-set.c 

libmxml_a_CPPFLAGS=$(AM_CPPFLAGS)
libmxml_a_CFLAGS= 



EXTRA_DIST = config.h CMakeLists.txt                                    \
             ${MXML_SRCDIR}/ANNOUNCEMENT    ${MXML_SRCDIR}/CHANGES      \
             ${MXML_SRCDIR}/COPYING         ${MXML_SRCDIR}/Makefile.in  \
             ${MXML_SRCDIR}/README          ${MXML_SRCDIR}/config.h.in  \
             ${MXML_SRCDIR}/configure       ${MXML_SRCDIR}/configure.in \
             ${MXML_SRCDIR}/install-sh      ${MXML_SRCDIR}/mxml.list.in \
             ${MXML_SRCDIR}/mxml.pc.in      ${MXML_SRCDIR}/mxml.spec    \
             ${MXML_SRCDIR}/mxml.xml        ${MXML_SRCDIR}/mxmldoc.c    \
             ${MXML_SRCDIR}/test.xml        ${MXML_SRCDIR}/testmxml.c   \
             ${MXML_SRCDIR}/mxml-private.h  ${MXML_SRCDIR}/mxml.h       \
             ${MXML_SRCDIR}/doc                                         \
             ${MXML_SRCDIR}/test                                        \
             ${MXML_SRCDIR}/vcnet

