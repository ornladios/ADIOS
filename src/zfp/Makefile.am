ZFP_SRCDIR=zfp-0.5.0

AM_CPPFLAGS=-I${srcdir}/${ZFP_SRCDIR}/inc -I${srcdir}/${ZFP_SRCDIR}/src

AUTOMAKE_OPTIONS = no-dependencies subdir-objects
METASOURCES=AUTO

CLEANFILES= 

noinst_LIBRARIES = libzfp.a
libzfp_a_SOURCES = ${ZFP_SRCDIR}/src/bitstream.c \
                   ${ZFP_SRCDIR}/src/decode1f.c \
                   ${ZFP_SRCDIR}/src/decode1d.c \
                   ${ZFP_SRCDIR}/src/encode1f.c \
                   ${ZFP_SRCDIR}/src/encode1d.c \
                   ${ZFP_SRCDIR}/src/decode2f.c \
                   ${ZFP_SRCDIR}/src/decode2d.c \
                   ${ZFP_SRCDIR}/src/encode2f.c \
                   ${ZFP_SRCDIR}/src/encode2d.c \
                   ${ZFP_SRCDIR}/src/decode3f.c \
                   ${ZFP_SRCDIR}/src/decode3d.c \
                   ${ZFP_SRCDIR}/src/encode3f.c \
                   ${ZFP_SRCDIR}/src/encode3d.c \
                   ${ZFP_SRCDIR}/src/zfp.c 


libzfp_a_CPPFLAGS=$(AM_CPPFLAGS)
libzfp_a_CFLAGS= 



EXTRA_DIST = ${ZFP_SRCDIR}/src/block1.h \
	     ${ZFP_SRCDIR}/src/block2.h \
             ${ZFP_SRCDIR}/src/block3.h \
	     ${ZFP_SRCDIR}/src/traitsd.h \
             ${ZFP_SRCDIR}/src/traitsf.h \
             ${ZFP_SRCDIR}/src/Makefile \
             ${ZFP_SRCDIR}/src/inline/inline.h \
             ${ZFP_SRCDIR}/src/template/codec.h \
             ${ZFP_SRCDIR}/src/template/codec1.c \
             ${ZFP_SRCDIR}/src/template/codec2.c \
             ${ZFP_SRCDIR}/src/template/codec3.c \
             ${ZFP_SRCDIR}/src/template/compress.c \
             ${ZFP_SRCDIR}/src/template/decode.c \
             ${ZFP_SRCDIR}/src/template/decode1.c \
             ${ZFP_SRCDIR}/src/template/decode2.c \
             ${ZFP_SRCDIR}/src/template/decode3.c \
             ${ZFP_SRCDIR}/src/template/decompress.c \
             ${ZFP_SRCDIR}/src/template/encode.c \
             ${ZFP_SRCDIR}/src/template/encode1.c \
             ${ZFP_SRCDIR}/src/template/encode2.c \
             ${ZFP_SRCDIR}/src/template/encode3.c \
             ${ZFP_SRCDIR}/src/template/template.h \
             ${ZFP_SRCDIR}/inc/bitstream.h \
             ${ZFP_SRCDIR}/inc/macros.h   \
             ${ZFP_SRCDIR}/inc/system.h   \
             ${ZFP_SRCDIR}/inc/types.h  \
             ${ZFP_SRCDIR}/inc/zfp.h \
             ${ZFP_SRCDIR}/API ${ZFP_SRCDIR}/Config ${ZFP_SRCDIR}/FAQ ${ZFP_SRCDIR}/LICENSE \
             ${ZFP_SRCDIR}/Makefile ${ZFP_SRCDIR}/README ${ZFP_SRCDIR}/VERSIONS \
             ${ZFP_SRCDIR}/array ${ZFP_SRCDIR}/examples 
