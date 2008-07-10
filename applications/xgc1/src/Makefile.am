# -*- Automake -*-
# search directory for m4 macros
ACLOCAL_AMFLAGS = -I ../acfiles

# main application
bin_PROGRAMS = xgc1

# conditionals
if COND_IMSL
IMSL_SRCS = 
else
IMSL_SRCS = taus88.F90 derf.F90 datanh.F90 fmin.F90 bspline90_22.F90
endif

if COND_PSPLINE
INTERP_SRCS = 
else
INTERP_SRCS = pppack.F90
endif

if COND_MPI
MPI_SRCS = mpi.F90
else
MPI_SRCS = mpisingle.F90
endif

# source code files
xgc1_SOURCES = module.F90 search.F90 charge.F90 read.F90 push.F90 \
	setup.F90 efield.F90 diagnosis.F90 limiter.F90 \
	bounce.F90 diagnosis2.F90 collision.F90 collision2.F90 \
	diagnosis-f.F90 heat.F90 turbulence.F90 neutral.F90 neutral2.F90 \
	linearsolver.F90 therm2d.F90 poisson.F90 \
	$(INTERP_SRCS) $(IMSL_SRCS) $(MPI_SRCS) \
	interpolation.F90 load.F90 main.F90

# conditional sources
EXTRA_xgc1_SOURCES = my_imsl.F90 bspline90_22.F90 taus88.F90 derf.F90 \
	datanh.F90 pppack.F90 fmin.F90 mpisingle.F90 mpi.F90

# compiler flags 
AM_FCFLAGS = $(FCFLAGS_f90)
xgc1_FCFLAGS = $(MPIINCLUDES) $(BINPACK_FLAGS) $(CAMTIMERS_FLAGS) \
	$(PSPLINE_FLAGS) $(IMSL_FLAGS) $(NETCDF_FLAGS) $(PETSC_FLAGS)

# loader flags
xgc1_LDADD = $(INTERP_OBJ) $(IMSL_OBJ) $(MPI_OBJ) \
	$(BINPACK_LIBS) $(CAMTIMERS_LIBS) \
	$(PSPLINE_LIBS) $(NETCDF_LIBS) $(PETSC_LIBS) $(MPILIBS)

# dependencies
xgc1_DEPENDENCIES = $(INTERP_OBJ) $(IMSL_OBJ) $(MPI_OBJ) \
	$(BINPACK_DEPS) $(CAMTIMERS_DEPS)

noinst_LIBRARIES = 

# binpack pack library
BPPACK_DIR = $(srcdir)/../binpack/pack
BPSHARE_DIR = $(srcdir)/../binpack/share
if COND_BINPACK
  noinst_LIBRARIES += libbppack.a
  libbppack_a_SOURCES = \
	$(BPPACK_DIR)/binpack-pack.c $(BPPACK_DIR)/binpack-read.c \
	$(BPSHARE_DIR)/binpack-utils.c $(BPPACK_DIR)/binpack-wrapper.c 
  libbppack_a_CFLAGS = -I$(BPPACK_DIR) -I$(BPSHARE_DIR)
endif

# CAM timing library
CAMTIMERS_DIR = $(srcdir)/../camtimers
CAMTIMERS_SRCS = \
	$(CAMTIMERS_DIR)/GPTLget_memusage.c \
	$(CAMTIMERS_DIR)/GPTLprint_memusage.c \
	$(CAMTIMERS_DIR)/GPTLutil.c \
	$(CAMTIMERS_DIR)/f_wrappers.c \
	$(CAMTIMERS_DIR)/gptl.c $(CAMTIMERS_DIR)/gptl_papi.c \
	$(CAMTIMERS_DIR)/threadutil.c \
	$(CAMTIMERS_DIR)/perf_utils.F90 $(CAMTIMERS_DIR)/perf_mod.F90
CAMTIMERS_FCFLAGS = -Mfree -i4 -Mdalign -Mextend -byteswapio \
	-DLINUX -I$(CAMTIMERS_DIR)
CAMTIMERS_CFLAGS = -DFORTRANUNDERSCORE -DLINUX -DSPMD -I$(CAMTIMERS_DIR)
if COND_PAPI
  CAMTIMERS_LIBNAME = libtimers_wpapi.a
else
  CAMTIMERS_LIBNAME = libtimers.a
endif
if COND_CAMTIMERS
  noinst_LIBRARIES += $(CAMTIMERS_LIBNAME)
  libtimers_wpapi_a_SOURCES = $(CAMTIMERS_SRCS)
  libtimers_wpapi_a_FCFLAGS = $(CAMTIMERS_FCFLAGS) $(PAPI_FLAGS)
  libtimers_wpapi_a_CFLAGS = $(CAMTIMERS_CFLAGS) $(PAPI_FLAGS)
  libtimers_a_SOURCES = $(CAMTIMERS_SRCS)
  libtimers_a_FCFLAGS = $(CAMTIMERS_FCFLAGS)
  libtimers_a_CFLAGS = $(CAMTIMERS_CFLAGS)
endif

# PETSc variable definitions
if COND_PETSC
PETSC_VARFILE = ${PETSC_DIR}/bmake/common/variables
else
PETSC_VARFILE = dummy.mk
endif
-include $(PETSC_VARFILE)

clean-local:
	-rm -f *.$(MODEXT)
