JAGUAR=n
EWOK=y
ifeq ($(JAGUAR),y)
     HOME=/ccs/proj/e2e/adios/jaguar
     FC=ftn
     CC=cc
endif

ifeq ($(EWOK),y)
     HOME=/ccs/proj/e2e/adios/ewok.gcc
     HOME1=/ccs/proj/e2e/adios/ewok.gcc
#     HOME=/ccs/home/ge1/work/ADIOS/trunk/src
#     HOME1=/ccs/home/ge1/work/ewok
     FC=mpif90
     CC=mpicc
endif

LIBD=-L${HOME} -L${HOME1}/lib -ladios -lmxml
INCD=-I${HOME} -I${HOME1}/include

all: heat_adios 

gwrite_output.ch: heat.xml
	python gpp.py heat.xml

heat_adios: heat_adios.c gwrite_output.ch
	$(CC) -g $(INCD) heat_adios.c -o heat_adios $(LIBD)

clean:
	rm -rf *demo* *.o *.dat *.bp *.h5 heat heat_adios heat_adios_global
