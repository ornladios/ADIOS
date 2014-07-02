SUBDIRS = etc src

dist_bin_SCRIPTS = bin/skel bin/set_method.sh bin/skel_extract.py bin/skel_cat.py

pythonlibdir=@libdir@/python

all-local:
	test "$(srcdir)" = "$(builddir)" || { \
	 cp -r $(srcdir)/bin $(builddir); \
	 for f in $(dist_pythonlib_DATA); do \
	     test ! -f "$(srcdir)/$$f" || cp $(srcdir)/"$$f" $(builddir)/lib; \
	 done; \
	 for f in $(dist_template_DATA); do \
	     test ! -f "$(srcdir)/$$f" || cp $(srcdir)/"$$f" $(builddir)/etc; \
	 done; \
	 cp -r $(srcdir)/etc/templates $(builddir)/etc; }

dist_pythonlib_DATA=\
    lib/adios.py \
	lib/argparse.py \
    lib/skel_bpls.py \
    lib/skel_bpy.py \
    lib/skel_install.py \
    lib/skel_makefile.py \
    lib/skel_params.py \
    lib/skel_replay.py \
    lib/skel_settings.py \
    lib/skel_source.py \
    lib/skel_submit.py \
    lib/skel_suite.py \
    lib/skel_template.py \
    lib/skel_test_plan.py \
    lib/skel_xml.py \
    lib/skelconf.py \
    lib/skel_have_adios_timing.py \
    lib/typeMapper.py 

templatedir=@prefix@/etc/skel/templates

dist_template_DATA=\
    etc/templates/Makefile.default.tpl\
	etc/templates/submit_sith.tpl\
	etc/templates/submit_nautilus.tpl\
	etc/templates/submit_jaguar.tpl \
        etc/templates/create_suite.tmpl \
        etc/templates/Makefile.tmpl \
        etc/templates/replay_bp.tmpl \
        etc/templates/replay.tmpl \
        etc/templates/replay_yaml.tmpl \
        etc/templates/source_write_c.tmpl \
        etc/templates/source_write_fortran.tmpl \
        etc/templates/submit_nautilus.tmpl \
        etc/templates/submit_sith.tmpl \
        etc/templates/xml.tmpl


EXTRA_DIST=etc/settings \
           etc/templates/Makefile.default.tpl \
           etc/templates/submit_jaguar.tpl \
           etc/templates/submit_sith.tpl \
           etc/templates/submit_nautilus.tpl \
           etc/templates/Makefile.tmpl \
           etc/templates/create_suite.tmpl \
           etc/templates/replay_bp.tmpl \
           etc/templates/replay.tmpl \
           etc/templates/replay_yaml.tmpl \
           etc/templates/source_write_c.tmpl \
           etc/templates/source_write_fortran.tmpl \
           etc/templates/submit_nautilus.tmpl \
           etc/templates/submit_sith.tmpl \
           etc/templates/xml.tmpl \
           lib/skel_have_adios_timing.py.in
           ##etc/compiler_fragment.mk 

