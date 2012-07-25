
# To modify this Makefile, edit the appropriate Makefile_<target>.tpl file in ~/.skel/
# and run skel makefile  
                    
ADIOS_DIR=/ccs/home/lot/adios/adios-install/jaguar

APP=$$APP$$
TESTS=$$TESTS$$
     
SKEL_DIR=$$SKEL_HOME$$

SKEL_SOURCE=$(SKEL_DIR)/bin/skel_source.py
SKEL_SUBMIT=$(SKEL_DIR)/bin/skel_submit.py
SKEL_PARAMS=$(SKEL_DIR)/bin/skel_params.py
SKEL_MAKEFILE=$(SKEL_DIR)/bin/skel_makefile.py
SKEL_XML=$(SKEL_DIR)/bin/skel_xml.py

CC=cc
FC=ftn -fast 
LDFLAGS=`$(ADIOS_DIR)/bin/adios_config -f -l`
CFLAGS=`$(ADIOS_DIR)/bin/adios_config -c`


DEST_DIR=$$DEPLOY_DIR$$/$(APP)/$$TARGET$$/$$CORES_USED$$


#$(TESTS): $(TESTS:%=%.c)
#	$(CC) $(CFLAGS) $@.c -o $@ $(LDFLAGS)

$(TESTS): $(TESTS:%=%.f90)
	$(FC) $(CFLAGS) $@.f90 -o $@ $(LDFLAGS)


#$(TESTS:%=%.c): $(APP).xml $(APP)_params.xml 
#	$(SKEL_SOURCE) $(APP)
#	$(SKEL_SUBMIT) $(APP)

$(TESTS:%=%.f90): $(APP).xml $(APP)_params.xml 
	$(SKEL_SOURCE) $(APP)
	$(SKEL_SUBMIT) $(APP)


makefile:
	$(SKEL_MAKEFILE) $(APP)

params: $(APP)_skel.xml
	$(SKEL_PARAMS) $(APP)

xml: $(APP).xml
	$(SKEL_XML) $(APP)

deploy:
	# Make sure this directory exists
	mkdir -p $(DEST_DIR)

	rm -f $(DEST_DIR)/submit*
	cp $(APP)_skel.xml $(DEST_DIR)/$(APP)_skel.xml.in
	cp $(TESTS) submit_$$TARGET$$* $(DEST_DIR)
	cp $(SKEL_DIR)/bin/set_method.sh $(DEST_DIR)
	cp $(SKEL_DIR)/bin/skel_cat.py $(DEST_DIR)

clean:
	rm -f *.c
	rm -f *.f90
	#rm -f submit*
	rm -f $(TESTS)





