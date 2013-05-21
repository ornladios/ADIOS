include_directories(${PROJECT_SOURCE_DIR}/utils/skel/src)
include_directories(${PROJECT_BINARY_DIR})
link_directories(${PROJECT_BINARY_DIR}/utils/skel/src)

add_library(skel STATIC skel_xml_output.c)
install(FILES skel_xml_output.h DESTINATION ${includedir}/skel)
install(FILES ${PROJECT_BINARY_DIR}/utils/skel/src/libskel.a DESTINATION ${libdir})
