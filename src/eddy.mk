ADIOS := n
CPP := mpiCC
CC := mpicc
F90 := mpif90
F77 := mpif90
LOADER := mpif90 -Wl,-z,muldefs -Wl,-z,noexecstack

include  $(SRCDIR)/git_commands.mk

# PE_ENV has to be set as an environment variable before running make.
# It does not get picked up since target.mk is called before eddy.mk
#PE_ENV := INTEL

ifeq ($(DEBUG),y)
  #LIB_DIR:=/home/ethier/software/INTEL18_DEBUG/lib
  #INC_DIR:=/home/ethier/software/INTEL18_DEBUG/include
  #MOD_DIR:=/home/ethier/software/INTEL18_DEBUG/mod
   PSPLINE_LIB:=/home/jinchen/LIB/INTEL18/lib
   PSPLINE_INC:=/home/jinchen/LIB/INTEL18/include
   PSPLINE_MOD:=/home/jinchen/LIB/INTEL18/mod
   SPRNG_LIB:=/home/jinchen/LIB/INTEL18/lib
   SPRNG_INC:=/home/jinchen/LIB/INTEL18/include
   SPRNG_MOD:=/home/jinchen/LIB/INTEL18/mod
   PETSC_DIR=/home/jinchen/LIB/INTEL18/petsc-3.12.0
   PETSC_ARCH=eddy-intel-openmpi-gts-2018-debug
else
  #LIB_DIR:=/home/ethier/software/INTEL18/lib
  #INC_DIR:=/home/ethier/software/INTEL18/include
  #MOD_DIR:=/home/ethier/software/INTEL18/mod
   PSPLINE_LIB:=/home/jinchen/LIB/INTEL18/intel-mpi/lib
   PSPLINE_INC:=/home/jinchen/LIB/INTEL18/intel-mpi/include
   PSPLINE_MOD:=/home/jinchen/LIB/INTEL18/intel-mpi/mod
   SPRNG_LIB:=/home/jinchen/LIB/INTEL18/intel-mpi/lib
   SPRNG_INC:=/home/jinchen/LIB/INTEL18/intel-mpi/include
   SPRNG_MOD:=/home/jinchen/LIB/INTEL18/intel-mpi/mod
 #2PETSC_DIR=/home/jinchen/LIB/petsc-3.8.3
 #2PETSC_ARCH=eddy-intel-openmpi-gts-2018
   PETSC_DIR=/home/jinchen/LIB/INTEL18/petsc-3.12.0
   PETSC_ARCH=eddy-intelmpi-gts-2019
 #1PETSC_DIR=/home/jinchen/LIB/petsc-3.7.7
 #1PETSC_ARCH=intel-openmpi-real
endif

FFLAGS:=-r8 -Dflush_=flush -assume no2underscore -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR)
OPT:=-O3
CFLAGS:=
COPT:=-O3
OMPFLAG:=-qopenmp
NOOMPFLAG:=
###LOCAL_LIB:=-L$(LIB_DIR) -lsprng -lpspline -lezcdf -lnetcdf
HDF5_DIR:=/home/myoo/software/eddy_GTS
HDF5_INC:=-I${HDF5_DIR}/include
HDF5_LIB:=-L${HDF5_DIR}/lib -lhdf5_hl -lhdf5 -lfftw3 -lhdf5_fortran -lhdf5hl_fortran
LOCAL_LIB:=-L${SPRNG_LIB} -lsprng -lgmp -L${PSPLINE_LIB} -lpspline
ifeq ($(DEBUG),y)
   PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lHYPRE -lX11 -lflapack -lfblas
else
   PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lHYPRE -lX11 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
endif
PETSC_INC:=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include #-I/usr/local/openmpi/3.0.0/intel180/x86_64/include

FFTLIB:=
NAG_FILES:=fake_nag.o nag.o
NAG:=-L$(VPATH)/cmlib/$(OBJDIR) -lcmlib
NAG_REPLACEMENT:=libcm
FEM:=-L$(VPATH)/fem/$(OBJDIR) -lfemlib
CFLAGS:=$(CFLAGS) -DAdd_ -I${PSPLINE_INC} $(SPRNG_INC) $(PETSC_INC)
FFLAGS:=$(FFLAGS) $(HDF5_INC) -I${PSPLINE_INC} -I${SPRNG_INC} -I${PSPLINE_MOD} -I${SPRNG_MOD} $(PETSC_INC)

   ifeq ($(ADIOS),y)   # Using ADIOS
      FFLAGS:=-DADIOS $(FFLAGS)
      SNAPSHOT:=snapshot_adios.o
      ADIOS_FLIB:=-ladiosf -lmxml
      CMD:=GTS_adios
   else
      SNAPSHOT:=snapshot.o
   endif

   ifeq ($(DEBUG),y)
     OPT:=-g  $(FFLAGS)
     COPT:= -g  $(CFLAGS)
     GOPT:= -g  $(CFLAGS)
   else
     OPT:=$(OPT) $(FFLAGS)
     COPT:=$(COPT) $(CFLAGS)
     GOPT:= -g  $(CFLAGS)
   endif

   ifeq ($(OPENMP),y)
     CMD:=GTS_OMP
     OPT:=$(OMPFLAG) $(OPT)
     COPT:=$(OMPFLAG) $(COPT)
     GOPT:= $(OMPFLAG) $(GOPT)
   else
     CMD:=GTS
     OPT:=$(NOOMPFLAG) $(OPT)
     COPT:=$(NOOMPFLAG) $(COPT)
     GOPT:= $(NOOMPFLAG) $(GOPT)
   endif

LIB:=$(NAG) $(HDF5_LIB) $(FEM) $(LOCAL_LIB) $(FFTLIB) $(PETSC_FLIB) $(ADIOS_FLIB)

%.o : %.cpp
	$(CPP) -c $(COPT) -c $< -o $@

%.o : %.c
	$(CC)  -c $(COPT) -c $< -o $@

%.o: %.f
	$(F77) -c $(OPT)  $< -o $@

%.o: %.F
	$(F77) -c $(OPT)  $< -o $@

%.o: %.f90
	$(F90) -c $(OPT)  $< -o $@

%.o: %.F90
	$(F90) -c $(OPT)  $< -o $@
