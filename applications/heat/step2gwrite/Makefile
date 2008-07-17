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

heat_adios: heat_adios.c
	python gwrite.py heat.xml temp
	$(CC) -I$(INCD) temp/heat_adios.c -o heat_adios -L$(LIBD)

clean:
	rm -rf *.h5 *.bp temp heat_adios_demo.o* *.o *.dat heat_adios 
