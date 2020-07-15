#############################################################################
#             Makefile to build the GTS code
#           ==================================
#
# You only need to type "gmake" to build the code on the platforms
# defined below. The makefile runs the "uname -s" command to detect
# the operating system automatically.
#
# Other options are:
#
#  % gmake OPENMP=y       Builds the code with OpenMP support
#
#  % gmake DEBUG=y        Compiles the files with debug option (-g)
#
#  % gmake ESSL=y         On AIX, uses the FFT routine from ESSL library
#                         The default is to use the NAG library routine
#
# You can combine more than one option on the command line:
#
#  % gmake OPENMP=y ESSL=y
#
# Special targets:
#
#  % gmake clean      Removes the executable and all object files (*.o *.mod *.lst)
#
#  % gmake cleanomp   Removes the executable and the object files
#                     containing OpenMP directives
#
#############################################################################

ifeq (,$(filter _%,$(notdir $(CURDIR))))
OPENMP:=n
DEBUG:=n
ESSL:=n
ACC:=n
CUDA:=n
ifeq ($(ACC), y)
  ifeq ($(CUDA), y)
    BIN_POSTFIX := $(BIN_POSTFIX)-ACC-CUDA
  else
	 BIN_POSTFIX := $(BIN_POSTFIX)-ACC
  endif
else
endif

ifeq ($(OPENMP), y)
  BIN_POSTFIX := $(BIN_POSTFIX)-OMP
else
endif

ifeq ($(DEBUG), y)
  BIN_POSTFIX := $(BIN_POSTFIX)-G
else
endif

include target.mk

else
#----- End Boilerplate

VPATH = $(SRCDIR)

include $(SRCDIR)/$(_ARCH).mk

CMD := GTS

# Common file for fft routine using NAG except on the CRAY T3E
FFT:=fft_gl.o

# List of all the object files needed to build the code
SOLVER:=solver_petsc.o
OBJ:=module.o spl_profile.o function_gg.o Structured_Grid.o ${SOLVER} ptl_src_loss.o specfun.o mgy_petsc.o setup_v2.o Er_provider.o Theory_jb.o \
     load.o restart.o diagnosis.o snapshot.o charge_ion.o current_ion.o \
     charge_eon.o current_eon.o solve_dphidt.o poisson.o smooth.o \
     field.o pushi_new.o pushe_new.o shifti.o shifte.o collision.o collision_eon.o CTPandCELL.o \
     $(FFT) escZ_march03.o check_eq.o view_distribution.o view_orbit.o $(NAG_FILES) \
     check_temperature.o get_fluxes.o get_3D_flow.o get_Te_fluctuation.o get_fluxes_e.o \
     get_ZF_shear.o open_files.o csr2css.o check_tracers.o create_directory.o\
     weight_damping_new.o weight_damping_eon_new2.o sobol.o main.o git_info.o
#output_profiles.o

#hypre openmp solver
ifeq ($(OPENMP), y)
OBJ:=$(OBJ) hypresolve.o
endif

all : libcm libfem $(CMD)

$(CMD): $(OBJ)
	$(LOADER) $(OBJ) $(OPT) $(LIB) -o $@

libcm:
	(cd $(VPATH)/cmlib; $(MAKE) FC=$(F77) FFLAGS="$(OPT)" ; cd ..)

#finite element and mesh generator
libfem: module.o
	(cd $(VPATH)/fem; $(MAKE) FC=$(F90) FFLAGS="$(OPT)" CC=$(CC) CFLAGS="$(COPT)" DEBUG="$(GOPT)"; cd ..)

$(OB): module.o

$(SOLVER): solver_petsc.F90
	$(F90) $(OPT) -DXUSERZ -c $< -o $@

hypresolve:
	-${CC} -o hypresolve -DXEX2_TEST $(VPATH)/hypresolve.c $(COPT) ${PETSC_FLIB}
	${RM} hypresolve.o

#----- Begin Boilerplate
endif
