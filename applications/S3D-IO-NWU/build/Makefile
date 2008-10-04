#  $Id: Makefile,v 1.6.6.1 2006/04/04 02:01:25 rsankar Exp $
#----------------------------------------------------------------------
# makefile for program S3D
#
# usage:  make MACH=[mach]
#         where [mach] corresponds to an extension on a "make.[mach]" file.
#
# For example, if I was building on an SGI, I would create a file called make.SGI
# (see below) and then issue the command "make MACH=SGI"
#
# Alternatively (and perhaps a better approach) you may set an environment variable
# "MACH" to the appropriate platform...  eg "setenv MACH SGI"
#
#----------------------------------------------------------------------

# the include file provides specific rules for building on your platform.
# You need to provide a compiler for F90 and F77 files, as well as a link 
# statement that will link with the appropriate libraries (MPI, math, etc).
#----------------------------------------------------------------------
# Revisions - 
# 01/03/05 - Ramanan Sankaran
# New version of makefile with dependencies auto-inserted using sfmakedepend
#----------------------------------------------------------------------

MACH = CNL
MACH = ewok

include make.$(MACH)

MYERRMSG = 
MAKEDEPEND = ./sfmakedepend

#-----------------------------------------------------------
# Set paths for various source files.

MODULES	= ../source/modules/
MAIN	= ../source/drivers/
INIT	= ../source/f90_files/init/
SOLVE	= ../source/f90_files/solve/
POST	= ../source/f90_files/post/
FSEVEN	= ../source/f77_files/
FSEVEN_HDR	= ../source/include/
MORPH	= ../source/utils/newmorph/
OTHER   = ../source/other_files/

#----------------------------------------------------------------------
# Default compilation rules based on file extension

.f.o:
	$(COMPILE_F77) -I$(FSEVEN_HDR) $<

%.o:%.f90
	$(COMPILE_F90) $<

.c.o:
	$(COMPILE_C) $<

#----------------------------------------------------------------------
#List of source files
SRC_MODULES = \
  $(MODULES)param_m.f90             \
	$(MODULES)chemkin_m.f90           \
	$(MODULES)topology_m.f90          \
	$(MODULES)work_m.f90              \
  $(MODULES)bc_m.f90                \
  $(MODULES)derivative_m.f90        \
  $(MODULES)grid_m.f90              \
  $(MODULES)reference_m.f90         \
  $(MODULES)runtime_m.f90           \
  $(MODULES)thermchem_m.f90         \
  $(MODULES)variables_m.f90         \
  $(MODULES)erk_m.f90               \
  $(MODULES)mpi_io_m.f90            \
  $(MODULES)pnetcdf_m.f90           \
  $(MODULES)hdf5_m.f90              \
  $(MODULES)adios_m.f90


SRC_DRIVERS =  \
  $(MAIN)main.f90                   \
  $(MAIN)solve_driver.f90           \
  $(MAIN)post_driver.f90

SRC_INIT = \
	$(INIT)init_pressure_wave_test.f90

SRC_SOLVE = \
	$(SOLVE)calc_cfl_limit.f90       \
	$(SOLVE)init_field.f90	         \
	$(SOLVE)io.f90		               \
	$(SOLVE)terminate_run.f90        

SRC_FSEVEN = \
	$(FSEVEN)ck2f90_tools.f		       \
	$(FSEVEN)ckheadersetup.f		     \
	$(FSEVEN)cklib.f			           \
	$(FSEVEN)getrates.f		           \
	$(FSEVEN)ckerr.f                 \
  $(FSEVEN)xerror.f                \
  $(FSEVEN)tranlib.f

#----------------------------------------------------------------------
# Add mkdir.c routine for machines on which SYSTEMCALLWONTWORK
ifeq ($(MACH), $(filter $(MACH), JETPTH JET CNL ewok))
  SRC_SOLVE += $(SOLVE)mkdir.c
endif 

#----------------------------------------------------------------------
# Choose the transport module to use
#ifndef $(TRANSPORT)
#  TRANSPORT = MIXAVG
#endif
#
#ifeq ($(TRANSPORT), LEWIS)
#  SRC_MODULES += $(MODULES)lewis_transport_m.f90
#else 
#  ifeq ($(TRANSPORT), MIXAVG)
#    SRC_MODULES += $(MODULES)mixavg_transport_m.f90
#  else
#    #Unknown option for INTMETHOD - error
#    MYERRMSG += Unknown option $(TRANSPORT) for TRANSPORT.
#  endif
#endif

#----------------------------------------------------------------------
# Choose the integration method

#SRC_ARK = \
#  $(MODULES)ark_m.f90              \
#  $(FSEVEN)g_getrates.f            \
#  $(FSEVEN)g_cklib.f               \
#  $(FSEVEN)eshfid.f 
#
#SRC_ERK = \
#  $(MODULES)erk_m.f90              \
#
#ifndef $(INTMETHOD)
#  INTMETHOD = ERK
#endif
#
#ifeq ($(INTMETHOD), ARK)
#  SRC_INT = $(SRC_ARK)
#else
#  ifeq ($(INTMETHOD), ERK)
#    SRC_INT = $(SRC_ERK)
#  else
#    #Unknown option for INTMETHOD - error
#		MYERRMSG += Unknown option $(INTMETHOD) for INTMETHOD.
#  endif
#endif

#----------------------------------------------------------------------
# Choose the math library to use
# Whether to use the math.f source file or the system lapack library

ifndef ($(MATHLIB))
  MATHLIB = FILE
endif

ifeq ($(MATHLIB), FILE)
  SRC_FSEVEN += $(FSEVEN)math.f
else
  ifeq ($(MATHLIB), SYSMATHLIB)
    LIBS += $(SYSMATHLIB)
  else
    #Unknown option for MATHLIB - error
		MYERRMSG += Unknown option $(MATHLIB) for MATHLIB.
  endif
endif

#----------------------------------------------------------------------
# Sources for the main program

SRC_S3D = \
  $(SRC_MODULES) $(SRC_DRIVERS) $(SRC_INT) \
  $(SRC_INIT) $(SRC_SOLVE) \
  $(SRC_BC) $(SRC_FSEVEN) 


# Extract the list of object files
OBJ_S3D2 = $(SRC_S3D:.f90=.o)
#One more pass for converting the .f extensions of f77 files
#Is there a smarter way?
OBJ_S3D1 = $(OBJ_S3D2:.f=.o)
OBJ_S3D  = $(OBJ_S3D1:.c=.o)

#----------------------------------------------------------------------
# TARGETS:

S3D = ../run/s3d_f90.x

morph = ../run/morph.x

# build the s3d_f90 executable by default.

default: errchk $(S3D)

errchk:
ifneq ($(MYERRMSG), )
	$(error Error:!! $(MYERRMSG))
endif

none: 
	echo "Doing nothing"

$(S3D): $(OBJ_S3D)
	$(LINK) $(OBJ_S3D) -o $(S3D) $(LIBS) $(IO_LIB)

morph: $(MORPH)morph.o
	$(LINK) $(MORPH)morph.o $(LIBS) -o ../run/morph.x

clean: 
	find .. -name "*.mod" -exec rm {} \; ; find .. -name "*.o" -exec rm {} \;

new: clean default

depend: $(SRC_S3D)
	$(MAKEDEPEND) $(DIRECTIVES) -I $(FSEVEN_HDR) $(SRC_S3D)


#----------------------------------------------------------------------
#Auto create the dependencies using make depend
# Do not add anything below the following line
# They will be removed when you do 'make depend'
# DO NOT DELETE THIS LINE - used by make depend
bc_m.xyz = ../source/modules/bc_m.o
chemkin_m.xyz = ../source/modules/chemkin_m.o
derivative_m.xyz = ../source/modules/derivative_m.o
grid_m.xyz = ../source/modules/grid_m.o
param_m.xyz = ../source/modules/param_m.o
reference_m.xyz = ../source/modules/reference_m.o
runtime_m.xyz = ../source/modules/runtime_m.o
thermchem_m.xyz = ../source/modules/thermchem_m.o
topology_m.xyz = ../source/modules/topology_m.o
variables_m.xyz = ../source/modules/variables_m.o
work_m.xyz = ../source/modules/work_m.o




../source/drivers/main.o: \
         $(chemkin_m.xyz) $(param_m.xyz) $(runtime_m.xyz) $(topology_m.xyz)

../source/drivers/post_driver.o: \
         $(chemkin_m.xyz) $(derivative_m.xyz) $(grid_m.xyz) $(param_m.xyz) \
         $(reference_m.xyz) $(runtime_m.xyz) $(thermchem_m.xyz) \
         $(topology_m.xyz) $(variables_m.xyz) $(work_m.xyz)

../source/drivers/solve_driver.o: \
         $(bc_m.xyz) $(derivative_m.xyz) $(grid_m.xyz) $(param_m.xyz) \
         $(reference_m.xyz) $(rk_m.xyz) $(runtime_m.xyz) $(thermchem_m.xyz) \
         $(topology_m.xyz) $(variables_m.xyz) $(work_m.xyz)
../source/f77_files/ck2f90_tools.o: \
         ../source/include//ckstrt.h

../source/f77_files/ck2f90_tools.o: \
        
../source/f77_files/ckheadersetup.o: \
         ../source/include//ckstrt.h

../source/f77_files/ckheadersetup.o: \
         $(topology_m.xyz)
../source/f77_files/cklib.o: \
         ../source/include//ckstrt.h

../source/f77_files/cklib.o: \
        

../source/f77_files/getrates.o: \
         $(thermchem_m.xyz)

../source/f90_files/init/init_pressure_wave_test.o: \
         $(chemkin_m.xyz) $(grid_m.xyz) $(param_m.xyz) $(reference_m.xyz) \
         $(topology_m.xyz) $(work_m.xyz)

../source/f90_files/solve/calc_cfl_limit.o: \
         $(grid_m.xyz) $(param_m.xyz) $(topology_m.xyz) $(work_m.xyz)

../source/f90_files/solve/init_field.o: \
         $(bc_m.xyz) $(grid_m.xyz) $(param_m.xyz) $(reference_m.xyz) \
         $(rk_m.xyz) $(runtime_m.xyz) $(thermchem_m.xyz) $(topology_m.xyz) \
         $(variables_m.xyz)

../source/f90_files/solve/io.o: \
         $(bc_m.xyz) $(grid_m.xyz) $(param_m.xyz) $(reference_m.xyz) \
         $(rk_m.xyz) $(runtime_m.xyz) $(thermchem_m.xyz) $(topology_m.xyz) \
         $(variables_m.xyz)

../source/f90_files/solve/terminate_run.o: \
         $(param_m.xyz) $(topology_m.xyz)

../source/modules/bc_m.o: \
         $(param_m.xyz) $(topology_m.xyz)

../source/modules/chemkin_m.o: \
         $(param_m.xyz) $(variables_m.xyz)

../source/modules/derivative_m.o: \
         $(param_m.xyz) $(topology_m.xyz)

../source/modules/grid_m.o: \
         $(param_m.xyz) $(reference_m.xyz) $(topology_m.xyz)

../source/modules/reference_m.o: \
         $(chemkin_m.xyz) $(topology_m.xyz)

../source/modules/thermchem_m.o: \
         $(chemkin_m.xyz) $(param_m.xyz) $(reference_m.xyz) $(topology_m.xyz)

../source/modules/variables_m.o: \
         $(param_m.xyz)

../source/modules/work_m.o: \
         $(param_m.xyz)

../source/modules/mpi_io_m.o: $(param_m.xyz) $(topology_m.xyz) $(variables_m.xyz) $(bc_m.xyz) $(runtime_m.xyz)

../source/modules/pnetcdf_m.o: ../source/modules/pnetcdf_m.f90 $(param_m.xyz) $(topology_m.xyz) $(variables_m.xyz) $(bc_m.xyz) $(runtime_m.xyz)
	$(COMPILE_F90) $(BUILD_IO_METHODS) -I$(FSEVEN_HDR) -I$(PNETCDF_DIR)/include $<

../source/modules/hdf5_m.o: ../source/modules/hdf5_m.f90 $(param_m.xyz) $(topology_m.xyz) $(variables_m.xyz) $(bc_m.xyz) $(runtime_m.xyz)
	$(COMPILE_F90) $(BUILD_IO_METHODS) -I$(FSEVEN_HDR) $(HDF5_FLIB) $<

../source/modules/adios_m.o: ../source/modules/adios_m.f90 $(param_m.xyz) $(topology_m.xyz) $(variables_m.xyz) $(bc_m.xyz) $(runtime_m.xyz)
	$(COMPILE_F90) $(BUILD_IO_METHODS) -I$(FSEVEN_HDR) $(HDF5_FLIB) $<
