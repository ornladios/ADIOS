ZFP_SRCDIR=zfp-0.5.5

AM_CPPFLAGS=-I${srcdir}/${ZFP_SRCDIR}/include -I${srcdir}/${ZFP_SRCDIR}/src

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
                   ${ZFP_SRCDIR}/src/decode1i.c \
                   ${ZFP_SRCDIR}/src/decode1l.c \
                   ${ZFP_SRCDIR}/src/decode2i.c \
                   ${ZFP_SRCDIR}/src/decode2l.c \
                   ${ZFP_SRCDIR}/src/decode3i.c \
                   ${ZFP_SRCDIR}/src/decode3l.c \
                   ${ZFP_SRCDIR}/src/decode4d.c \
                   ${ZFP_SRCDIR}/src/decode4f.c \
                   ${ZFP_SRCDIR}/src/decode4i.c \
                   ${ZFP_SRCDIR}/src/decode4l.c \
                   ${ZFP_SRCDIR}/src/encode1i.c \
                   ${ZFP_SRCDIR}/src/encode1l.c \
                   ${ZFP_SRCDIR}/src/encode2i.c \
                   ${ZFP_SRCDIR}/src/encode2l.c \
                   ${ZFP_SRCDIR}/src/encode3i.c \
                   ${ZFP_SRCDIR}/src/encode3l.c \
                   ${ZFP_SRCDIR}/src/encode4d.c \
                   ${ZFP_SRCDIR}/src/encode4f.c \
                   ${ZFP_SRCDIR}/src/encode4i.c \
                   ${ZFP_SRCDIR}/src/encode4l.c \
                   ${ZFP_SRCDIR}/src/zfp.c 


libzfp_a_CPPFLAGS=$(AM_CPPFLAGS)
libzfp_a_CFLAGS= 



EXTRA_DIST = ${ZFP_SRCDIR}/src/block1.h \
             ${ZFP_SRCDIR}/src/block2.h \
             ${ZFP_SRCDIR}/src/block3.h \
             ${ZFP_SRCDIR}/src/block4.h \
             ${ZFP_SRCDIR}/src/traitsd.h \
             ${ZFP_SRCDIR}/src/traitsi.h \
             ${ZFP_SRCDIR}/src/traitsf.h \
             ${ZFP_SRCDIR}/src/traitsl.h \
             ${ZFP_SRCDIR}/src/Makefile \
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
             ${ZFP_SRCDIR}/src/template/codec4.c \
             ${ZFP_SRCDIR}/src/template/codecf.c \
             ${ZFP_SRCDIR}/src/template/cudacompress.c \
             ${ZFP_SRCDIR}/src/template/cudadecompress.c \
             ${ZFP_SRCDIR}/src/template/decode4.c \
             ${ZFP_SRCDIR}/src/template/decodef.c \
             ${ZFP_SRCDIR}/src/template/decodei.c \
             ${ZFP_SRCDIR}/src/template/encode4.c \
             ${ZFP_SRCDIR}/src/template/encodef.c \
             ${ZFP_SRCDIR}/src/template/encodei.c \
             ${ZFP_SRCDIR}/src/template/ompcompress.c \
             ${ZFP_SRCDIR}/src/template/revcodecf.c \
             ${ZFP_SRCDIR}/src/template/revdecode1.c \
             ${ZFP_SRCDIR}/src/template/revdecode2.c \
             ${ZFP_SRCDIR}/src/template/revdecode3.c \
             ${ZFP_SRCDIR}/src/template/revdecode4.c \
             ${ZFP_SRCDIR}/src/template/revdecode.c \
             ${ZFP_SRCDIR}/src/template/revdecodef.c \
             ${ZFP_SRCDIR}/src/template/revencode1.c \
             ${ZFP_SRCDIR}/src/template/revencode2.c \
             ${ZFP_SRCDIR}/src/template/revencode3.c \
             ${ZFP_SRCDIR}/src/template/revencode4.c \
             ${ZFP_SRCDIR}/src/template/revencode.c \
             ${ZFP_SRCDIR}/src/template/revencodef.c \
             ${ZFP_SRCDIR}/src/template/template.h \
             ${ZFP_SRCDIR}/include/bitstream.h \
             ${ZFP_SRCDIR}/include/zfp/macros.h   \
             ${ZFP_SRCDIR}/include/zfp/system.h   \
             ${ZFP_SRCDIR}/include/zfp/types.h  \
             ${ZFP_SRCDIR}/include/zfp.h \
             ${ZFP_SRCDIR}/Config ${ZFP_SRCDIR}/LICENSE ${ZFP_SRCDIR}/CMakeLists.txt \
             ${ZFP_SRCDIR}/cmake ${ZFP_SRCDIR}/cfp ${ZFP_SRCDIR}/fortran \
             ${ZFP_SRCDIR}/Makefile ${ZFP_SRCDIR}/README.md ${ZFP_SRCDIR}/VERSIONS.md \
             ${ZFP_SRCDIR}/array ${ZFP_SRCDIR}/examples \
             ${ZFP_SRCDIR}/python ${ZFP_SRCDIR}/tests ${ZFP_SRCDIR}/utils 
