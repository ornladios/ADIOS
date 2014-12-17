ACLOCAL_AMFLAGS=-I config

SUBDIRS=src utils tests examples

EXTRA_DIST = runconf \
             KNOWN_BUGS \
             adios_config.in \
             scripts/create.nssi.config.sh \
             scripts/kill.nssi.staging.sh \
             scripts/start.nssi.staging.sh \
	     scripts/FindADIOS.cmake

#             wrappers/matlab/Makefile \
#             wrappers/matlab/adios.m \
#             wrappers/matlab/adiosclose.m \
#             wrappers/matlab/adiosopen.m \
#             wrappers/matlab/adiosread.m \
#             wrappers/matlab/adiosclosec.c \
#             wrappers/matlab/adiosopenc.c \
#             wrappers/matlab/adiosreadc.c 

bin_SCRIPTS = adios_config 

if HAVE_NSSI
    bin_SCRIPTS += scripts/create.nssi.config.sh \
                   scripts/kill.nssi.staging.sh  \
                   scripts/start.nssi.staging.sh
endif

sysconf_DATA = scripts/FindADIOS.cmake

#install-data-hook:
#	echo -n "VERSIONSTRING=\"$(PACKAGE) version $(VERSION) svn revision " >> adios_config.flags
#	if type svnversion &>/dev/null; then \
#	    echo `svnversion`\" >> adios_config.flags; \
#	else \
#	    echo "unknown\"" >> adios_config.flags; \
#	    echo "Error: svnversion command not found"; \
#	fi 
#	cp adios_config.flags $(DESTDIR)$(bindir)
        
install-data-hook:
	echo "SRCDIR=\"$(PWD)\"" >> adios_config.flags; \
	if git rev-parse --git-dir &>/dev/null; then \
	    echo "GITLOG=\"`git log -n 1 --format=\"Last commit on %ad hash %h\"`\"" > adios_git_status; \
	    echo "read -d '' GITSTAT <<'_EOF_'" >> adios_git_status; \
	    git status -uno >> adios_git_status; \
	    echo "_EOF_" >> adios_git_status; \
	else \
	    echo "GITLOG=\"git log not available\"" >> adios_git_status; \
	    echo "GITSTAT=\"git status not available\"" >> adios_git_status; \
	    echo "Note: git command not found or $(PWD) is not a git repository."; \
	fi 
	$(INSTALL_DATA) adios_config.flags $(DESTDIR)$(sysconfdir)
	$(INSTALL_DATA) adios_git_status $(DESTDIR)$(sysconfdir)
        
CLEANFILES = adios_config adios_config.flags
