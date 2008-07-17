JAGUAR=n
EWOK=y
ifeq ($(JAGUAR),y)
     HOME=/ccs/proj/e2e/adios/jaguar
     FC=ftn
     CC=cc
endif

ifeq ($(EWOK),y)
     HOME=/ccs/proj/e2e/adios/ewok.gcc
     FC=mpif90
     CC=mpicc
endif

LIBD=${HOME}/lib -ladios -lmxml
INCD=${HOME}/include

all: heat_adios 

heat_adios: heat_adios.c heat.xml
	python gwrite.py heat.xml temp
	$(CC) -I$(INCD)  temp/heat_adios.c -o heat_adios -L$(LIBD)

clean:
	rm -rf temp *demo* *.o *.dat *.bp *.h5 heat heat_adios heat_adios_global
