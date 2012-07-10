
# To modify this Makefile, edit the appropriate Makefile_<target>.tpl file in ~/.skel/templates/
# and run skel makefile  

ADIOS_DIR=/opt/adios
                         
APP=$$APP$$
TESTS=$$TESTS$$

SKEL_DIR=$$SKEL_HOME$$

SKEL_SOURCE=$(SKEL_DIR)/bin/skel_source.py
SKEL_MAKEFILE=$(SKEL_DIR)/bin/skel_makefile.py
SKEL_SUBMIT=$(SKEL_DIR)/bin/skel_submit.py
SKEL_PARAMS=$(SKEL_DIR)/bin/skel_params.py
SKEL_XML=$(SKEL_DIR)/bin/skel_xml.py

FC=mpif90 -ffree-line-length-256
LDFLAGS=-L/opt/adios/lib -ladiosf -L/opt/mxml/lib -L/opt/nc4par/lib -L/opt/phdf5/lib -lm -lmxml -lpthread -lnetcdf -lhdf5_hl -lhdf5 -lz
CFLAGS=-I/opt/adios/include -I/opt/mxml/include -I/opt/nc4par/include -I/opt/phdf5/include

FCFLAGS=-I/opt/adios/include -I/opt/mxml/include -I/opt/nc4par/include -I/opt/phdf5/include

DEST_DIR=$$DEPLOY_DIR$$/$(APP)/$$TARGET$$/$$CORES_USED$$


$(TESTS): $(TESTS:%=%.f90)
	$(FC) $(CFLAGS) $@.f90 -o $@ $(LDFLAGS)


$(TESTS:%=%.f90): $(APP).xml $(APP)_params.xml 
	$(SKEL_SOURCE) $(APP)
	$(SKEL_SUBMIT) $(APP)
	chmod u+x submit*


makefile:
	$(SKEL_MAKEFILE) $(APP)

params: $(APP)_skel.xml
	$(SKEL_PARAMS) $(APP)

xml: $(APP).xml
	$(SKEL_XML) $(APP)
	cp $(APP)_skel.xml $(APP)_skel.xml.in
 
deploy:
	#Make sure this exists
	mkdir -p $(DEST_DIR)

	rm -f $(DEST_DIR)/submit*
	cp $(APP)_skel.xml $(DEST_DIR)/$(APP)_skel.xml.in
	cp $(TESTS) submit_sith* $(DEST_DIR)
	cp $(SKEL_DIR)/bin/set_method.sh $(DEST_DIR)

clean:
	rm -f *.c *.f90
	#rm -f submit*
	rm -f $(TESTS)

