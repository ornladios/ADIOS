.SUFFIXES: .o .F90

CMD=xgc1
ALL:$(CMD)

include ${PETSC_DIR}/bmake/common/variables

OBJ=module.o $(EXTRA_OBJ) search.o charge.o main.o read.o push.o load.o \
	setup.o efield.o interpolation.o $(MPI_OBJ) diagnosis.o limiter.o \
	bounce.o diagnosis2.o collision.o collision2.o diagnosis-f.o heat.o \
	turbulence.o neutral.o neutral2.o linearsolver.o therm2d.o poisson.o
IMSL_OBJ=my_imsl.o
PORT_OBJ=bspline90_22.o taus88.o derf.o datanh.o pppack.o fmin.o
SER_OBJ=mpisingle.o
PAR_OBJ=mpi.o

