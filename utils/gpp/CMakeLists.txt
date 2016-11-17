CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/gpp.py.in ${CMAKE_CURRENT_BINARY_DIR}/gpp.py @ONLY)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gpp.py DESTINATION ${bindir})
install(FILES ad_config.py type_mapper.py DESTINATION ${libdir}/python)

file(COPY ad_config.py type_mapper.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/) 
