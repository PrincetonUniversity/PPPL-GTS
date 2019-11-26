
  CPP := mpicxx
  CC := mpicc
  F90 := mpif90
  F77 := mpif77
  LOADER := mpif90 -Wl,-z,muldefs -Wl,-z,noexecstack

  include  $(SRCDIR)/git_commands.mk

DFLAGS:= -D_MGY_TIMER_
ACCFLAG:=-acc -Minfo=accel

ifeq ($(F90),mpif90)
       ifeq ($(PE_ENV),PGI)   # Using PGI compiler
          FFLAGS:= -r8 -Kieee -Dflush_=flush -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR)
          OPT:=-fast
          COPT:=-fast
          CFLAGS:=
          OMPFLAG:=-mp=nonuma
          NOOMPFLAG:=
          SPRNG_DIR=/home/jinchen/LIB/PGI195
          SPRNG_LIB:=-L${SPRNG_DIR}/lib -lsprng
          SPRNG_INC:=-I${SPRNG_DIR}/include
          PSPLINE_DIR=/home/jinchen/LIB/PGI195
          PSPLINE_INC:=-I${PSPLINE_DIR}/mod
          PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline
          PETSC_DIR=/home/jinchen/LIB/petsc-3.12.0
          ifeq ($(CUDA),y)
				PETSC_DIR=/home/jinchen/LIB/petsc
				PETSC_ARCH=traverse-pgi-openmpi-gts-cuda-2019-2-32bit
				#PETSC_ARCH=traverse-pgi-openmpi-gts-cuda-2019
				CUDA_LIB:=-L/usr/local/cuda-10.1/lib64 -lcudart -lcusparse -lcusolver -lstdc++ -lcublas
				# CUDA_LIB:=-L/usr/local/cuda-10.1/lib64 -lcusparse
				PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lHYPRE -lhdf5 -lz -Wl,--start-group /opt/pgi/19.5/linuxpower/19.5/lib/liblapack.a /opt/pgi/19.5/linuxpower/19.5/lib/libblas.a -Wl,--end-group -lpthread -lm -ldl  -lparmetis -lmetis
			 else
				PETSC_ARCH=traverse-pgi-openmpi-gts-2019
			 	CUDA_LIB:=
				PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lHYPRE -lhdf5 -lz -Wl,--start-group /opt/pgi/19.5/linuxpower/19.5/lib/liblapack.a /opt/pgi/19.5/linuxpower/19.5/lib/libblas.a -Wl,--end-group -lpthread -lm -ldl
			 endif
          HDF5_LIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lhdf5_hl -lhdf5
          NETCDF_LIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lnetcdf
          PETSC_FINC:= -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
#-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
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
          HDF5_LIB:=-L${PETSC_DIR}/lib -lhdf5_hl -lhdf5 -lfftw3
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
   FFLAGS:=$(FFLAGS) $(SPRNG_INC) $(NETCDF_INC) $(PSPLINE_INC) $(PETSC_FINC) $(HDF5_INC)

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

   ifeq ($(ACC),y)
     OPT:=$(ACCFLAG) $(OPT)
     COPT:=$(ACCFLAG) $(COPT)
     GOPT:= $(ACCFLAG) $(GOPT)
	 endif

	OPT:=$(OPT) $(DFLAGS)
	COPT:=$(COPT) $(DFLAGS)
	GOPT:= -g  $(DFLAGS)


   ifeq ($(OPENMP),y)
     OPT:=$(OMPFLAG) $(OPT)
     COPT:=$(OMPFLAG) $(COPT)
     GOPT:= $(OMPFLAG) $(GOPT)
   else
     OPT:=$(NOOMPFLAG) $(OPT)
     COPT:=$(NOOMPFLAG) $(COPT)
     GOPT:= $(NOOMPFLAG) $(GOPT)
   endif

LIB:=-L$HDF5DIR/lib64 -lhdf5 -lhdf5_fortran -lhdf5hl_fortran $(CUDA_LIB) $(NAG) $(FEM) $(SPRNG_LIB) $(PSPLINE_LIB) $(HDF5_LIB) $(FFTLIB) $(PETSC_FLIB)

%.o : %.cpp
	$(CPP)  $(COPT) -c $< -o $@

%.o : %.c
	$(CC)  $(COPT) -c $< -o $@

%.o: %.f
	$(F77) $(OPT) -c $< -o $@

%.o: %.F
	$(F77) $(OPT) -c $< -o $@

%.o: %.f90
	$(F90) -cpp $(OPT) -c $< -o $@

%.o: %.F90
	$(F90) $(OPT) -c $< -o $@
