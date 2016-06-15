set(MXML_SRCDIR "mxml-2.9")
include_directories(${PROJECT_BINARY_DIR}/mxml)
include_directories(${PROJECT_SOURCE_DIR}/src/mxml/${MXML_SRCDIR})
link_directories(${PROJECT_BINARY_DIR}/mxml)


set(CLEANFILES "")
set(include_HEADERS "")
set(nodist_include_HEADERS "")
set(lib_LIBRARIES "")

if(BUILD_MXML)
    set(libmxml_la_SOURCES 
            ${MXML_SRCDIR}/mxml-attr.c
            ${MXML_SRCDIR}/mxml-get.c 
            ${MXML_SRCDIR}/mxml-private.c
            ${MXML_SRCDIR}/mxml-string.c
            ${MXML_SRCDIR}/mxml-entity.c
            ${MXML_SRCDIR}/mxml-index.c
            ${MXML_SRCDIR}/mxml-search.c
            ${MXML_SRCDIR}/mxml-file.c 
            ${MXML_SRCDIR}/mxml-node.c
            ${MXML_SRCDIR}/mxml-set.c
       )
    set(libmxml_la_CPPFLAGS "${AM_CPPFLAGS}")
    set(libmxml_la_CFLAGS "")

    add_library(mxml STATIC ${libmxml_la_SOURCES})
    set_target_properties(mxml PROPERTIES COMPILE_FLAGS "${MXML_CPPFLAGS}")

#add_subdirectory(mxml-2.9)

install(FILES ${PROJECT_BINARY_DIR}/mxml/libmxml.la
        DESTINATION ${libdir})

endif(BUILD_MXML)
