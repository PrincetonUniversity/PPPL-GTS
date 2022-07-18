
  CPP := CC
  CC := cc
  F90 := ftn
  F77 := ftn
  LOADER := ftn -Wl,-z,muldefs

  include  $(SRCDIR)/git_commands.mk


ifeq ($(F90),ftn)
       ifeq ($(PE_ENV),PGI)   # Using PGI compiler on Hopper
       else   # Using CRAY compiler
          ifeq ($(PE_ENV),NVIDIA)   # Using the INTEL compiler on Hopper or Edison
             FFLAGS:=-r8 -Dflush_=flush -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR)
             OPT:=-fast 
             CFLAGS:=
             COPT:=-fast 
             OMPFLAG:=-qopenmp -DUSE_OPENMP
             NOOMPFLAG:=
             SPRNG_DIR=/global/cfs/cdirs/mp288/jinchen/LIB/NVIDIA_perlmutter
             SPRNG_LIB:=-L${SPRNG_DIR}/lib -lsprng -lgmp          # for Edison
             SPRNG_INC:=-I${SPRNG_DIR}/include
             PSPLINE_DIR=/global/cfs/cdirs/mp288/jinchen/LIB/NVIDIA_perlmutter
             PSPLINE_INC:=-I${PSPLINE_DIR}/mod
             PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline #-L${NETCDF_DIR}/lib -lnetcdf
             PETSC_DIR=/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609
	     PETSC_ARCH=perlmuttercpu-nvidia-gts
             PETSC_FINC:= -I${PETSC_DIR}/include -I${PETSC_DIR}/$(PETSC_ARCH)/include
	     PETSC_FLIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/perlmuttercpu-nvidia-gts/lib -L/global/cfs/cdirs/mp288/jinchen/PETSC/petsc.20220609/perlmuttercpu-nvidia-gts/lib -L/opt/cray/pe/mpich/8.1.17/ofi/nvidia/20.7/lib -L/opt/cray/pe/libsci/21.08.1.2/NVIDIA/20.7/x86_64/lib -L/opt/cray/pe/dsmml/0.2.2/dsmml/lib -L/global/common/software/nersc/pm-2021q4/sw/darshan/3.3.1/lib -L/opt/cray/xpmem/2.3.2-2.2_7.5__g93dd7ee.shasta/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/compilers/lib -L/usr/lib64/gcc/x86_64-suse-linux/7 -Wl,-rpath,/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/compilers/lib -lpetsc -lHYPRE -lnetcdf -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lgmp -lstdc++ -ldl -ldarshan -lz -lsci_nvidia_mpi -lsci_nvidia -lmpifort_nvidia -lmpi_nvidia -ldsmml -lxpmem -lnvf -lnvomp -lnvhpcatm -latomic -lpthread -lnvcpumath -lnsnvc -lnvc -lrt -lgcc_s -lm -lquadmath -lstdc++ -ldl


          else  # Using the CRAY compiler (cce)
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
	$(CC)  $(COPT) -c $< -o $@

%.o: %.f
	$(F77) $(OPT) -c $< -o $@

%.o: %.F
	$(F77) $(OPT) -c $< -o $@

%.o: %.f90
	$(F90) $(OPT) -c $< -o $@

%.o: %.F90
	$(F90) $(OPT) -c $< -o $@
