add_subdirectory(etc)
add_subdirectory(src)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/lib/skel_have_adios_timing.py.in ${CMAKE_CURRENT_BINARY_DIR}/lib/skel_have_adios_timing.py)

install(PROGRAMS bin/skel bin/set_method.sh bin/skel_extract.py bin/skel_cat.py DESTINATION ${bindir})

install(FILES lib/adios.py 
	lib/argparse.py 
	lib/skel_bpls.py
	lib/skel_bpy.py
	lib/skel_install.py 
	lib/skel_makefile.py
	lib/skel_params.py
	lib/skel_replay.py
	lib/skel_settings.py 
	lib/skel_source.py 
	lib/skel_submit.py 
	lib/skel_xml.py 
	lib/skelconf.py
	${CMAKE_CURRENT_BINARY_DIR}/lib/skel_have_adios_timing.py 
	lib/typeMapper.py 
	DESTINATION ${libdir}/python)


install(FILES etc/templates/Makefile.default.tpl
	etc/templates/submit_sith.tpl
	etc/templates/submit_jaguar.tpl
	etc/templates/Makefile.tmpl
	etc/templates/source_write_c.tmpl
	etc/templates/source_write_fortran.tmpl
        etc/templates/submit_sith.tmpl
	etc/templates/xml.tmpl
	etc/templates/replay_bp.tmpl
	etc/templates/replay_yaml.tmpl
	DESTINATION ${prefix}/etc/skel/templates)




#install(FILES etc/settings etc/compiler_fragment.mk
#	etc/templates/Makefile.default.tpl
#	etc/templates/submit_jaguar.tpl
#	etc/templates/submit_sith.tpl
#	DESTINATION ${PROJECT_BINARY_DIR}/utils/skel)

#install(FILES bin/skel bin/set_method.sh bin/skel_extract.py DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#install(FILES lib/adios.py lib/argparse.py lib/skel_cat.py lib/skel_install.py lib/skel_makefile.py lib/skel_params.py lib/skel_settings.py lib/skel_source.py lib/skel_submit.py lib/skel_xml.py lib/skelconf.py lib/typeMapper.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python)
