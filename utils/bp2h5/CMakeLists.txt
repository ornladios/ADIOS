include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/public)
include_directories(${PROJECT_SOURCE_DIR}/utils/bp2h5)
link_directories(${PROJECT_BINARY_DIR}/utils/bp2h5)

add_executable(bp2h5 bp2h5.c)
target_link_libraries(bp2h5 adiosread_nompi ${ADIOSREADLIB_SEQ_LDADD} ${HDF5_LIBS})
set_target_properties(bp2h5 PROPERTIES COMPILE_FLAGS "-DH5_USE_16_API -DADIOS_USE_READ_API_1 ${ADIOSLIB_EXTRA_CPPFLAGS} ${ADIOSREADLIB_SEQ_CPPFLAGS} ${HDF5_CPPFLAGS}")
