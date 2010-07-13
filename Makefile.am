SUBDIRS=src utils tests examples

EXTRA_DIST = runconf \
             KNOWN_BUGS \
             adios_config \
             tools/matlab/Makefile \
             tools/matlab/adios.m \
             tools/matlab/adiosclose.m \
             tools/matlab/adiosopen.m \
             tools/matlab/adiosread.m \
             tools/matlab/adiosclosec.c \
             tools/matlab/adiosopenc.c \
             tools/matlab/adiosreadc.c \
             scripts/create.nssi.config.sh \
             scripts/kill.nssi.staging.sh \
             scripts/start.nssi.staging.sh

bin_SCRIPTS = adios_config 

if HAVE_NSSI
    bin_SCRIPTS += scripts/create.nssi.config.sh \
                   scripts/kill.nssi.staging.sh  \
                   scripts/start.nssi.staging.sh
endif

install-data-hook:
	cp adios_config.flags $(DESTDIR)$(bindir)
        
