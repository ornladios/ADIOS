include_directories(${PROJECT_SOURCE_DIR}/utils/skel/etc)
include_directories(${PROJECT_BINARY_DIR}/utils/skel/etc)
link_directories(${PROJECT_BINARY_DIR}/utils/skel/etc)

set(etcdir ${prefix}/etc/skel)

#ADD_EXECUTABLE(all compiler_fragment.mk)

#ADD_CUSTOM_TARGET(Dummy ALL)

#add_custom_command(TARGET Dummy
#        DEPENDS compiler_fragment.mk)


add_custom_target(dummy ALL
    DEPENDS compiler_fragment.mk)

add_custom_command(
	OUTPUT compiler_fragment.mk
	COMMAND cat compiler_fragment.tmp | sed -e 's/CMAKE_C_COMPILER:FILEPATH=/CC=/' -e 's/CMAKE_C_FLAGS:STRING=/CFLAGS=/' -e 's/CMAKE_Fortran_COMPILER:FILEPATH=/FC=/' -e 's/CMAKE_Fortran_FLAGS:STRING=/FCFLAGS=/' > compiler_fragment.mk
	DEPENDS compiler_fragment.tmp)

add_custom_command(
	OUTPUT compiler_fragment.tmp
	COMMAND grep -e'^prefix\ =' -e '^CMAKE_C_COMPILER:FILEPATH=' -e'^CMAKE_Fortran_COMPILER:FILEPATH=' -e'^CMAKE_C_FLAGS:STRING=' -e'^CMAKE_Fortran_FLAGS:STRING='  < ${PROJECT_BINARY_DIR}/CMakeCache.txt > compiler_fragment.tmp
	DEPENDS ${PROJECT_BINARY_DIR}/CMakeCache.txt)

install(FILES ${PROJECT_BINARY_DIR}/utils/skel/etc/compiler_fragment.mk settings DESTINATION ${etcdir})
