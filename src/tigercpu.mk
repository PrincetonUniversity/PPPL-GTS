
   CPP := mpicxx
  CC := mpicc
  F90 := mpif90
  F77 := mpif77

  LOADER := mpifc -Wl,-z,muldefs -Wl,-z,noexecstack

include  $(SRCDIR)/git_commands.mk


 # If working tree is dirty, append dirty flag (D)
ifneq ($(strip $(shell git status --porcelain 2>/dev/null)),)
  GIT_VERSION := $(GIT_VERSION)-D
  GIT_WTREE := $(shell git status --porcelain)
endif

 GIT_FLAGS := 	-DGIT_VERSION="\"$(GIT_VERSION)\""\
 					-DGIT_BRANCH="\"$(GIT_BRANCH)\""\
					-DGIT_COM_MSG="\"$(GIT_COM_MSG)\""\
 					-DGIT_COM_SHA="\"$(GIT_COM_SHA)\""\
					-DGIT_WTREE="\"$(GIT_WTREE)\""\
 					-DGIT_DATE="\"$(GIT_DATE)\""\
 					-DBUILD_DATE="\"$(BUILD_DATE)\""\
 					-DBUILD_DIR="\"$(BUILD_DIR)\""

 ifeq ($(F90),mpif90)
       ifeq ($(PE_ENV),PGI)   # Using PGI compiler
          FFLAGS:= -r8 -Kieee -Dflush_=flush -module .
          OPT:=-O2
          COPT:=-O2
          CFLAGS:=
          OMPFLAG:=-mp=nonuma
          NOOMPFLAG:=
          SPRNG_LIB:=-L${SPRNG_DIR}/lib -lsprng
          SPRNG_INC:=-I${SPRNG_DIR}/include
          PSPLINE_INC:=-I${PSPLINE_DIR}/mod
          PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline -lezcdf -L${NETCDF_DIR}/lib -lnetcdf
       else   # Using INTEL compiler
          FFLAGS:=-r8 -Dflush_=flush -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR)
          OPT:=-O2 -no-ipo
          CFLAGS:=
          COPT:=-O2 -no-ipo
          OMPFLAG:=-qopenmp -DUSE_OPENMP
          NOOMPFLAG:=
          NETCDF_INC:=-I${NETCDFDIR}/include
          SPRNG_DIR=/home/ethier/software/tigercpu_intel19_GTS
          SPRNG_LIB:=-L${SPRNG_DIR}/lib -lsprng
          SPRNG_INC:=-I${SPRNG_DIR}/include
          PSPLINE_DIR=/home/ethier/software/tigercpu_intel19_GTS
          PSPLINE_INC:=-I${PSPLINE_DIR}/mod
          PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline
          #PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline -L${NETCDF_DIR}/lib -lnetcdf
          HDF5_DIR=/home/myoo/software/tigercpu_GTS
          HDF5_INC:=-I${HDF5_DIR}/include
          HDF5_LIB:=-L${HDF5_DIR}/lib -lhdf5_hl -lhdf5 -lfftw3 -lhdf5_fortran -lhdf5hl_fortran
          PETSC_DIR=/home/ethier/software/tigercpu_intel19_GTS
          PETSC_ARCH=
          PETSC_FINC:= -I${PETSC_DIR}/include
          #PETSC_FINC:= -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include
          ###PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lcraypetsc_intel_real
          PETSC_FLIB:=-L${PETSC_DIR}/lib -lpetsc -lHYPRE -lhdf5 -lz  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
          #PETSC_FLIB:=-L${PETSC_DIR}/lib -lpetsc -lHYPRE -llapack -lblas -lhdf5 -lz
          #PETSC_FLIB:=-L${PETSC_DIR}/$(PETSC_ARCH)/lib -lpetsc -lHYPRE -lflapack -lfblas -lscalapack -lsuperlu -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -ldl -lstdc++
       endif

    FFTLIB:=
   NAG_FILES:=fake_nag.o nag.o
   NAG:=-L$(VPATH)/cmlib/$(OBJDIR) -lcmlib
   NAG_REPLACEMENT:=libcm
   FEM:=-L$(VPATH)/fem/$(OBJDIR) -lfemlib
   CFLAGS:=$(CFLAGS) -DAdd_ ${SPRNG_INC} $(PETSC_FINC)
   FFLAGS:=$(FFLAGS) $(HDF5_INC) $(SPRNG_INC) $(NETCDF_INC) $(PSPLINE_INC) $(PETSC_FINC)

endif

    ifeq ($(DEBUG),y)
     OPT:=-g  $(FFLAGS)   -dynamic -debug inline-debug-info -parallel-source-info=2 -O3
     COPT:= -g  $(CFLAGS) -dynamic -debug inline-debug-info -parallel-source-info=2 -O3
     GOPT:= -g  $(CFLAGS) -dynamic -debug inline-debug-info -parallel-source-info=2 -O3
   else
     OPT:=$(OPT) $(FFLAGS)
     COPT:=$(COPT) $(CFLAGS)
     GOPT:= -g  $(CFLAGS)
   endif

    ifeq ($(OPENMP),y)
     OPT:=$(OMPFLAG) $(OPT)
     COPT:=$(OMPFLAG) $(COPT)
     GOPT:= $(OMPFLAG) $(GOPT)
   else
     OPT:=$(NOOMPFLAG) $(OPT)
     COPT:=$(NOOMPFLAG) $(COPT)
     GOPT:= $(NOOMPFLAG) $(GOPT)
   endif

 LIB:=$(NAG) $(HDF5_LIB) $(FEM) $(SPRNG_LIB) $(PSPLINE_LIB) $(FFTLIB) $(PETSC_FLIB)

 %.o : %.cpp
	$(CPP) $(COPT) -c $< -o $@

 %.o : %.c
	$(CC) $(COPT) -I$(VPATH)/$(OBJDIR) $(GIT_FLAGS) -c $< -o $@

 %.o: %.f
	$(F77) $(OPT) -c $< -o $@

 %.o: %.F
	$(F77) $(OPT) -c $< -o $@

 %.o: %.f90
	$(F90) $(OPT) -c $< -o $@

 %.o: %.F90
	$(F90) $(OPT) -c $< -o $@
