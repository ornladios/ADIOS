SUBDIRS = etc src

dist_bin_SCRIPTS = bin/skel bin/set_method.sh bin/skel_extract.py bin/skel_cat.py

pythonlibdir=@libdir@/python

dist_pythonlib_DATA=\
    lib/adios.py \
	lib/argparse.py \
    lib/skel_bpls.py \
    lib/skel_bpy.py \
    lib/skel_install.py \
    lib/skel_makefile.py \
    lib/skel_params.py \
    lib/skel_settings.py \
    lib/skel_source.py \
    lib/skel_submit.py \
    lib/skel_xml.py \
    lib/skelconf.py \
    lib/skel_have_adios_timing.py \
    lib/typeMapper.py 

templatedir=@prefix@/etc/skel/templates

dist_template_DATA=\
    etc/templates/Makefile.default.tpl\
	etc/templates/submit_sith.tpl\
	etc/templates/submit_nautilus.tpl\
	etc/templates/submit_jaguar.tpl

EXTRA_DIST=etc/settings \
           etc/templates/Makefile.default.tpl \
           etc/templates/submit_jaguar.tpl \
           etc/templates/submit_sith.tpl \
           etc/templates/submit_nautilus.tpl \
           lib/skel_have_adios_timing.py.in
           ##etc/compiler_fragment.mk 

