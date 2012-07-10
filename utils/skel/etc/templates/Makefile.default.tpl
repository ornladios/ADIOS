
# To modify this Makefile, edit the appropriate Makefile_<target>.tpl file in ~/.skel/templates/
# and run skel makefile  

include $$INCLUDE$$

ADIOS_DIR=`adios_config -d`
LDFLAGS += `adios_config -l`
CFLAGS += `adios_config -c`
FCFLAGS += `adios_config -fc`          
FCLIBS += `adios_config -fl`
               
APP=$$APP$$
CTESTS=$$CTESTS$$
FTESTS=$$FTESTS$$

all: $(FTESTS) $(CTESTS)


$(CTESTS): $(CTESTS:=.c)
	$(CC) $(CFLAGS) -o $@ ${@}.c $(LDFLAGS)

$(FTESTS): $(FTESTS:=.f90)
	$(FC) $(FCFLAGS) -o $@ ${@}.f90 $(FCLIBS)

SKEL_DIR=$$SKEL_HOME$$

SKEL_SOURCE=$(SKEL_DIR)/bin/skel_source.py
SKEL_MAKEFILE=$(SKEL_DIR)/bin/skel_makefile.py
SKEL_SUBMIT=$(SKEL_DIR)/bin/skel_submit.py
SKEL_PARAMS=$(SKEL_DIR)/bin/skel_params.py
SKEL_XML=$(SKEL_DIR)/bin/skel_xml.py



DEST_DIR=$$DEPLOY_DIR$$/$(APP)/$$CORES_USED$$



#$(TESTS:%=%.c): $(APP).xml $(APP)_params.xml 
#	$(SKEL_SOURCE) $(APP)
#	$(SKEL_SUBMIT) $(APP)

 
deploy:
	#Make sure this exists
	mkdir -p $(DEST_DIR)

	rm -f $(DEST_DIR)/submit*
	cp $(APP)_skel.xml $(DEST_DIR)/$(APP)_skel.xml.in
	cp $(CTESTS) submit* $(DEST_DIR)
	cp $(FTESTS) submit* $(DEST_DIR)
	cp $(SKEL_DIR)/bin/set_method.sh $(DEST_DIR)

clean:
	rm -f *.c
	rm -f *.f90
	#rm -f submit*
	rm -f $(CTESTS) $(FTESTS)


