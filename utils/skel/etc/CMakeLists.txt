include_directories(${PROJECT_SOURCE_DIR}/utils/skel/etc)
include_directories(${PROJECT_BINARY_DIR}/utils/skel/etc)
link_directories(${PROJECT_BINARY_DIR}/utils/skel/etc)

set(etcdir ${prefix}/etc/skel)

install(FILES ${PROJECT_BINARY_DIR}/utils/skel/etc/compiler_fragment.mk settings DESTINATION ${etcdir})

#ADD_EXECUTABLE(all compiler_fragment.mk)

ADD_CUSTOM_TARGET(Dummy ALL)

add_custom_command(TARGET Dummy
        DEPENDS compiler_fragment.mk)


#add_custom_target(dummy
#    DEPENDS compiler_fragment.mk)

add_custom_command(
	OUTPUT compiler_fragment.mk
	COMMAND grep -e'^prefix\ =' -e'^CC\ =' -e'^FC\ =' -e'^CFLAGS\ =' -e'^LDFLAGS\ =' -e'^FCFLAGS\ =' -e'^FCLIBS\ =' < Makefile > compiler_fragment.mk
	DEPENDS ${PROJECT_BINARY_DIR}/utils/skel/Makefile)
