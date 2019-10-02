set(ZFP_SRCDIR "zfp-0.5.5")
include_directories(${PROJECT_SOURCE_DIR}/src/zfp/${ZFP_SRCDIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src/zfp/${ZFP_SRCDIR})
include_directories(${PROJECT_BINARY_DIR}/src/zfp)
link_directories(${PROJECT_BINARY_DIR}/src/zfp)


set(CLEANFILES "")
set(include_HEADERS "")
set(nodist_include_HEADERS "")
set(lib_LIBRARIES "")

if(BUILD_ZFP)
    set(libzfp_a_SOURCES 
        ${ZFP_SRCDIR}/src/bitstream.c
        ${ZFP_SRCDIR}/src/decode1f.c
        ${ZFP_SRCDIR}/src/decode1d.c
        ${ZFP_SRCDIR}/src/encode1f.c
        ${ZFP_SRCDIR}/src/encode1d.c
        ${ZFP_SRCDIR}/src/decode2f.c
        ${ZFP_SRCDIR}/src/decode2d.c
        ${ZFP_SRCDIR}/src/encode2f.c
        ${ZFP_SRCDIR}/src/encode2d.c
        ${ZFP_SRCDIR}/src/decode3f.c
        ${ZFP_SRCDIR}/src/decode3d.c
        ${ZFP_SRCDIR}/src/encode3f.c
        ${ZFP_SRCDIR}/src/encode3d.c
        ${ZFP_SRCDIR}/src/decode1i.c
        ${ZFP_SRCDIR}/src/decode1l.c
        ${ZFP_SRCDIR}/src/decode2i.c
        ${ZFP_SRCDIR}/src/decode2l.c
        ${ZFP_SRCDIR}/src/decode3i.c
        ${ZFP_SRCDIR}/src/decode3l.c
        ${ZFP_SRCDIR}/src/decode4d.c
        ${ZFP_SRCDIR}/src/decode4f.c
        ${ZFP_SRCDIR}/src/decode4i.c
        ${ZFP_SRCDIR}/src/decode4l.c
        ${ZFP_SRCDIR}/src/encode1i.c
        ${ZFP_SRCDIR}/src/encode1l.c
        ${ZFP_SRCDIR}/src/encode2i.c
        ${ZFP_SRCDIR}/src/encode2l.c
        ${ZFP_SRCDIR}/src/encode3i.c
        ${ZFP_SRCDIR}/src/encode3l.c
        ${ZFP_SRCDIR}/src/encode4d.c
        ${ZFP_SRCDIR}/src/encode4f.c
        ${ZFP_SRCDIR}/src/encode4i.c
        ${ZFP_SRCDIR}/src/encode4l.c
        ${ZFP_SRCDIR}/src/zfp.c 
       )
    set(libzfp_a_CPPFLAGS "${AM_CPPFLAGS}")
    set(libzfp_a_CFLAGS "${AM_CFLAGS}")

    add_library(zfp OBJECT ${libzfp_a_SOURCES})
    #set_target_properties(zfp PROPERTIES COMPILE_FLAGS "")

## Do not install libzfp.a, adios does not need it
#install(FILES ${PROJECT_BINARY_DIR}/zfp/libzfp.a
#        DESTINATION ${libdir})

endif(BUILD_ZFP)
