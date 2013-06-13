include_directories(${PROJECT_SOURCE_DIR}/utils/skel/src)
include_directories(${PROJECT_BINARY_DIR})
link_directories(${PROJECT_BINARY_DIR}/utils/skel/src)

if(BUILD_FORTRAN)
  add_library(skel STATIC skel_xml_output.c skel_xml_output_f.c)
else()
  add_library(skel STATIC skel_xml_output.c)
endif(BUILD_FORTRAN)

install(FILES skel_xml_output.h DESTINATION ${includedir}/skel)
install(FILES ${PROJECT_BINARY_DIR}/utils/skel/src/libskel.a DESTINATION ${libdir})
