
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

all: heat 

heat: heat.c
	$(CC) heat.c -o heat

clean:
	rm -f *demo* *.o *.dat heat
