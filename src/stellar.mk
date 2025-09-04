ADIOS := n
CPP := mpiCC
CC := mpicc
F90 := mpif90
F77 := mpif90
LOADER := mpif90

include  $(SRCDIR)/git_commands.mk

# PE_ENV has to be set as an environment variable before running make.
# It does not get picked up since target.mk is called before eddy.mk
#PE_ENV := INTEL

ifeq ($(DEBUG),y)
  #LIB_DIR:=/home/ethier/software/INTEL20201_DEBUG/lib
  #INC_DIR:=/home/ethier/software/INTEL20201_DEBUG/include
  #MOD_DIR:=/home/ethier/software/INTEL20201_DEBUG/mod
   PSPLINE_LIB:=/home/jinchen/LIB/INTEL20201/lib
   PSPLINE_INC:=/home/jinchen/LIB/INTEL20201/include
   PSPLINE_MOD:=/home/jinchen/LIB/INTEL20201/mod
   SPRNG_LIB:=/home/jinchen/LIB/INTEL20201/lib
   SPRNG_INC:=/home/jinchen/LIB/INTEL20201/include
   SPRNG_MOD:=/home/jinchen/LIB/INTEL20201/mod
   PETSC_DIR=/home/jinchen/LIB/PETSC
   PETSC_ARCH=
else
  #LIB_DIR:=/home/ethier/software/INTEL20201/lib
  #INC_DIR:=/home/ethier/software/INTEL20201/include
  #MOD_DIR:=/home/ethier/software/INTEL20201/mod
   PSPLINE_LIB:=/home/jinchen/LIB/INTEL20201/lib
   PSPLINE_INC:=/home/jinchen/LIB/INTEL20201/include
   PSPLINE_MOD:=/home/jinchen/LIB/INTEL20201/mod
   SPRNG_LIB:=/home/jinchen/LIB/INTEL20201/lib
   SPRNG_INC:=/home/jinchen/LIB/INTEL20201/include
   SPRNG_MOD:=/home/jinchen/LIB/INTEL20201/mod
   PETSC_DIR=/home/jinchen/LIB/PETSC
   PETSC_ARCH=
endif

FFLAGS:=-r8 -Dflush_=flush -assume no2underscore -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR)
OPT:=-O3
CFLAGS:=
COPT:=-O3
OMPFLAG:=-qopenmp
NOOMPFLAG:=
###LOCAL_LIB:=-L$(LIB_DIR) -lsprng -lpspline -lezcdf -lnetcdf
LOCAL_LIB:=-L${SPRNG_LIB} -lsprng -lgmp -L${PSPLINE_LIB} -lpspline -L$(NETCDFDIR)/lib64/ -lnetcdf -L$(HDF5DIR)/lib64 -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran
ifeq ($(DEBUG),y)
   #PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lHYPRE -lX11 -lflapack -lfblas
   PETSC_FLIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/home/jinchen/LIB/PETSC/lib -L/home/jinchen/LIB/PETSC/lib -L/usr/local/openmpi/4.1.0/intel20211/lib64 -L/usr/local/hdf5/intel-2021.1/openmpi-4.1.0/1.10.6/lib64 -L/usr/local/fftw/intel-2021.1/openmpi-4.1.0/3.3.9/lib64 -L/opt/intel/oneapi/mkl/2021.1.1/lib/intel64 -L/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8 -L/opt/intel/oneapi/compiler/2021.1.2/linux/lib -L/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/8 -Wl,-rpath,/usr/local/openmpi/4.1.0/intel20211/lib64 -lpetsc -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lflapack -lfblas -lparmetis -lmetis -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -lquadmath -lstdc++ -ldl
else
   #PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lHYPRE -lX11 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
   PETSC_FLIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/home/jinchen/LIB/PETSC/lib -L/home/jinchen/LIB/PETSC/lib -L/usr/local/openmpi/4.1.0/intel20211/lib64 -L/usr/local/hdf5/intel-2021.1/openmpi-4.1.0/1.10.6/lib64 -L/usr/local/fftw/intel-2021.1/openmpi-4.1.0/3.3.9/lib64 -L/opt/intel/oneapi/mkl/2021.1.1/lib/intel64 -L/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8 -L/opt/intel/oneapi/compiler/2021.1.2/linux/lib -L/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/8 -Wl,-rpath,/usr/local/openmpi/4.1.0/intel20211/lib64 -lpetsc -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu -lsuperlu_dist -lflapack -lfblas -lparmetis -lmetis -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -lquadmath -lstdc++ -ldl
   PETSC_INC:=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include #-I/usr/local/openmpi/3.0.0/intel180/x86_64/include
endif

   FFTLIB:=
   ifeq ($(NAG),y)   # Use nag libraries
	NAG_FILES:=fake_nag.o nag.o
	NAG:=-L$(VPATH)/cmlib/$(OBJDIR) -lcmlib
   else   # Use nag replacement
	NAG_FILES:=fake_nag.o
	NAG:=-L$(VPATH)/cmlib/$(OBJDIR) -lcmlib -L$(VPATH)/quadpack -lquadpack
   endif
   FEM:=-L$(VPATH)/fem/$(OBJDIR) -lfemlib

CFLAGS:=$(CFLAGS) -DAdd_ -I${PSPLINE_INC} $(SPRNG_INC) $(PETSC_INC)
FFLAGS:=$(FFLAGS) -I${PSPLINE_INC} -I${SPRNG_INC} -I${PSPLINE_MOD} -I${SPRNG_MOD} $(PETSC_INC)

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

LIB:=$(NAG) $(FEM) $(LOCAL_LIB) $(FFTLIB) $(PETSC_FLIB) $(ADIOS_FLIB)

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

