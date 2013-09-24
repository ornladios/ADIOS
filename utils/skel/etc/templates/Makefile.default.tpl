
# To modify this Makefile, edit Makefile_default.tpl file in ~/.skel/templates/
# and run skel makefile  

include $$INCLUDE$$

ADIOS_DIR=`adios_config -d`
LDFLAGS += `adios_config -l` -L${ADIOS_DIR}/lib/skel -lskel 
CFLAGS += `adios_config -c`
FCFLAGS += `adios_config -fc`          
FCLIBS += `adios_config -fl` -L${ADIOS_DIR}/lib/skel -lskel
               
APP=$$APP$$
CTESTS=$$CTESTS$$
FTESTS=$$FTESTS$$

DEST_DIR=$$DEPLOY_DIR$$/$(APP)/$$CORES_USED$$

all: $(FTESTS) $(CTESTS)


$(CTESTS): $(CTESTS:=.c)
	$(CC) $(CFLAGS) -o $@ ${@}.c $(LDFLAGS)

$(FTESTS): $(FTESTS:=.f90)
	$(FC) $(FCFLAGS) -o $@ ${@}.f90 $(FCLIBS)


deploy:
	#Make sure this exists
	mkdir -p $(DEST_DIR)

	rm -f $(DEST_DIR)/submit*
	cp $(APP)_skel.xml $(DEST_DIR)/$(APP)_skel.xml.in
	cp $(CTESTS) submit* $(DEST_DIR)
	cp $(FTESTS) submit* $(DEST_DIR)
	cp $(prefix)/bin/set_method.sh $(DEST_DIR)
	cp $(prefix)/bin/skel_cat.py $(DEST_DIR)

clean:
	rm -f *.c
	rm -f *.f90
	#rm -f submit*
	rm -f $(CTESTS) $(FTESTS)


