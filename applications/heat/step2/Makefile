JAGUAR=y
EWOK=n
ifeq ($(JAGUAR),y)
     HOME=/ccs/proj/e2e/adios/jaguar
     FC=ftn
     CC=cc
endif

ifeq ($(ewok),y)
     HOME=/ccs/proj/e2e/adios/ewok.gcc
     FC=mpif90
     CC=mpicc
endif

LIBD=${HOME}/lib -ladios -lmxml
INCD=${HOME}/include

all: heat_adios 

heat_adios: heat_adios.c
	$(CC) -I$(INCD) -o heat_adios heat_adios.c -L$(LIBD) 

clean:
	rm -f *.nc *demo* *.txt *.h5 *.bp *.o *.dat heat_adios 
