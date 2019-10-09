
  CPP := CC
  CC := cc
  F90 := ftn
  F77 := ftn
  LOADER := ftn -Wl,-z,muldefs

include  $(SRCDIR)/git_commands.mk

ifeq ($(F90),ftn)
   # Flags for Cray XT system
   #####ifeq ($(shell uname -n | cut -c-3),hop)
   #####ifeq ($(NERSC_HOST),hopper)
       ifeq ($(PE_ENV),PGI)   # Using PGI compiler on Hopper
          FFLAGS:= -r8 -Kieee -Dflush_=flush -module .
          OPT:=-fast
          COPT:=-fast
          CFLAGS:=
          OMPFLAG:=-mp=nonuma
          NOOMPFLAG:=
          SPRNG_LIB:=-L${SPRNG_DIR}/lib -lsprng
          SPRNG_INC:=-I${SPRNG_DIR}/include
          PSPLINE_INC:=-I${PSPLINE_DIR}/mod
          PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline -lezcdf -L${NETCDF_DIR}/lib -lnetcdf
       else   # Using CRAY compiler
          ifeq ($(PE_ENV),INTEL)   # Using the INTEL compiler on Hopper or Edison
             FFLAGS:=-r8 -Dflush_=flush -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR)
             OPT:=-fast -no-ipo
             CFLAGS:=
             COPT:=-fast -no-ipo
             OMPFLAG:=-qopenmp -DUSE_OPENMP
             NOOMPFLAG:=
             SPRNG_DIR=/project/projectdirs/mp288/jinchen/LIB/INTEL_cori_haswell
             SPRNG_LIB:=-L${SPRNG_DIR}/lib -lsprng -lgmp          # for Edison
             SPRNG_INC:=-I${SPRNG_DIR}/include
             PSPLINE_DIR=/project/projectdirs/mp288/jinchen/LIB/INTEL_cori_haswell
             PSPLINE_INC:=-I${PSPLINE_DIR}/mod
             PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline -L${NETCDF_DIR}/lib -lnetcdf
             PETSC_DIR=/global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.7.7
             PETSC_ARCH=cori-hsw-intel1903199-mpich776
             PETSC_FINC:= -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include
             ###PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lcraypetsc_intel_real
             PETSC_FLIB:=-L${PETSC_DIR}/$(PETSC_ARCH)/lib -lpetsc -lHYPRE -lflapack -lfblas -lscalapack -lsuperlu -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -ldl -lstdc++
          else  # Using the CRAY compiler (cce)
             #initialize array to zero: -e0 -ez
             #check all array bound:-m2 -R bcsp, but too slow, so not used unless debug code
             #FFLAGS:=-s real64 -Dflush_=flush -I. -static -e F -e m -rm
             #FFLAGS:= -O 2 -s real64 -Dflush_=flush -I. -static -e F -e m -rm -m2 -R bcsp -e0 -ez
             #CFLAGS:= -m2 -R bcsp -e0 -ez
             FFLAGS:= -s real64 -Dflush_=flush -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR) -static -e F -e m -rm -e0 -ez
             OPT:= -O2
             CFLAGS:= -e0 -ez
             COPT:= -O2
             OMPFLAG:= -DUSE_OPENMP
             NOOMPFLAG:=-hnoomp
             #PSPLINE_DIR:=/project/projectdirs/mp19/software/cray
             #SPRNG_DIR:=/project/projectdirs/mp19/software/cray_hopper
             SPRNG_DIR=/project/projectdirs/mp288/jinchen/LIB/CRAY_cori_haswell
             SPRNG_INC:=-I${SPRNG_DIR}/include -I${SPRNG_DIR}/mod
             SPRNG_LIB:=-L${SPRNG_DIR}/lib -lsprng -lgmp
             #PSPLINE_DIR:=/usr/common/usg/pspline/nersc1.1/cray
             #PSPLINE_DIR:=/global/homes/j/jinchen/project/LIB/pspline
             PSPLINE_DIR=/project/projectdirs/mp288/jinchen/LIB/CRAY_cori_haswell
             PSPLINE_INC:=-I${PSPLINE_DIR}/include -I${PSPLINE_DIR}/mod -I$(NETCDF_DIR)/include
             PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline -L${NETCDF_DIR}/lib -lnetcdf
             PETSC_DIR=/global/project/projectdirs/mp288/jinchen/PETSC/petsc-3.7.7
             #PETSC_ARCH=cori-hsw-cray604-mpich773
             PETSC_ARCH=cori-hsw-cray605-mpich776-noomp
             PETSC_FINC:= -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include
             #PETSC_FLIB:=-L${PETSC_DIR}/$(PETSC_ARCH)/lib -lcraypetsc_cray_real -ldl -lstdc++ -L$(HYPRE_DIR)/lib -lHYPRE
             PETSC_FLIB:=-L${PETSC_DIR}/$(PETSC_ARCH)/lib -lpetsc -lHYPRE -lflapack -lfblas -lscalapack -lsuperlu -lsuperlu_dist -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -ldl -lstdc++
          endif
       endif

   FFTLIB:=
   NAG_FILES:=fake_nag.o nag.o
   NAG:=-L$(VPATH)/cmlib/$(OBJDIR) -lcmlib
   NAG_REPLACEMENT:=libcm
   FEM:=-L$(VPATH)/fem/$(OBJDIR) -lfemlib
   CC:=cc
   CFLAGS:=$(CFLAGS) -DAdd_ ${SPRNG_INC} $(PETSC_FINC)
   FFLAGS:=$(FFLAGS) $(SPRNG_INC) $(NETCDF) $(PSPLINE_INC) $(PETSC_FINC) $(HDF5)

   #FFLAGS:=$(FFLAGS) -DUSERZ
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

LIB:=$(NAG) $(FEM) $(NETCDF) $(SPRNG_LIB) $(PSPLINE_LIB) $(HDF5) $(FFTLIB) $(PETSC_FLIB)

%.o : %.cpp
	$(CPP)  $(COPT) -c $< -o $@

%.o : %.c
	$(CC)  $(COPT) -I$(VPATH)/$(OBJDIR) $(GIT_FLAGS) -c $< -o $@

%.o: %.f
	$(F77) $(OPT) -c $< -o $@

%.o: %.F
	$(F77) $(OPT) -c $< -o $@

%.o: %.f90
	$(F90) $(OPT) -c $< -o $@

%.o: %.F90
	$(F90) $(OPT) -c $< -o $@