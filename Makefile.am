ACLOCAL_AMFLAGS=-I config

SUBDIRS=src utils tests examples

EXTRA_DIST = runconf \
             KNOWN_BUGS \
             adios_config \
             scripts/create.nssi.config.sh \
             scripts/kill.nssi.staging.sh \
             scripts/start.nssi.staging.sh

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
	    echo "GITLOG=\"`git log -n 1 --format=\"Last commit on %ad hash %h\"`\"" >> adios_config.flags; \
	    echo "read -d '' GITSTAT <<'_EOF_'" >> adios_config.flags; \
	    git status -uno >> adios_config.flags; \
	    echo "_EOF_" >> adios_config.flags; \
	else \
	    echo "GITLOG=\"\"" >> adios_config.flags; \
	    echo "GITSTAT=\"\"" >> adios_config.flags; \
	    echo "Error: git command not found"; \
	fi 
	cp adios_config.flags $(DESTDIR)$(bindir)
        
