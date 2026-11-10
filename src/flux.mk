ADIOS := n
CPP := mpiCC
CC := mpicc
F90 := mpif90
F77 := mpif77
LOADER := mpif90

include  $(SRCDIR)/git_commands.mk

# PE_ENV has to be set as an environment variable before running make.
# It does not get picked up since target.mk is called before eddy.mk
#PE_ENV := INTEL

ifeq ($(DEBUG),y)
else
   PSPLINE_LIB:=/u/jchen/GTS-GITREPO/LIB/INTEL18/lib
   PSPLINE_INC:=/u/jchen/GTS-GITREPO/LIB/INTEL18/include
   PSPLINE_MOD:=/u/jchen/GTS-GITREPO/LIB/INTEL18/mod
   SPRNG_LIB:=/u/jchen/GTS-GITREPO/LIB/INTEL18/lib
   SPRNG_INC:=/u/jchen/GTS-GITREPO/LIB/INTEL18/include
   SPRNG_MOD:=/u/jchen/GTS-GITREPO/LIB/INTEL18/mod
   NETCDFFDIR=${NETCDF_F_HOME}
   NETCDFCDIR=${NETCDF_C_HOME}
   PETSC_DIR=/p/swim/jchen/PETSC/petsc.20231214
   PETSC_ARCH=real-spark-gts
   GMP_DIR=/opt/pppl/spack-pkgs/linux-rocky9-zen4/gcc-11.3.1/gmp-6.2.1-h3i43rjfj72vm7pfdiq2krab65rfvokp
endif

FFLAGS:=-r8 -Dflush_=flush -assume no2underscore -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR)
OPT:=-O3
CFLAGS:=
COPT:=-O3
OMPFLAG:=-qopenmp
NOOMPFLAG:=
###LOCAL_LIB:=-L$(LIB_DIR) -lsprng -lpspline -lezcdf -lnetcdf
LOCAL_LIB:=-L${SPRNG_LIB} -lsprng -L${PSPLINE_LIB} -lpspline -L$(NETCDF_C_HOME)/lib64/ -lnetcdf -L$(NETCDF_F_HOME)/lib/ -lnetcdff

ifeq ($(DEBUG),y)
else
	PETSC_WITH_EXTERNAL_LIB = -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/p/swim/jchen/PETSC/petsc.20231214/real-spark-gts/lib -L/p/swim/jchen/PETSC/petsc.20231214/real-spark-gts/lib -Wl,-rpath,/opt/pppl/intel/2023.2.0/mkl/2023.2.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/mkl/2023.2.0/lib/intel64 -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/openmpi-4.1.5-u5n2mrzof35eduxvqlsilu2akf43ewup/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/hwloc-2.9.1-tfqm6yyezgbdjtpzlcsmpkcew2omhyok/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/libevent-2.1.12-mbhu3djdevlnislas2gssh5uy4obbuki/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/pmix-5.0.1-5pfyaosw4gg7t3avbartnn55yi2yhwru/lib -L/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-13.2.0/slurm-23-02-5-1-rjtbt7npfgbrtdjr4b2qnd3gfrospvcn/lib -L/opt/pppl/intel/2023.2.0/ipp/2021.9.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/ippcp/2021.8.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/dnnl/2023.2.0/cpu_dpcpp_gpu_dpcpp/lib -L/opt/pppl/intel/2023.2.0/dal/2023.2.0/lib/intel64 -L/opt/pppl/intel/2023.2.0/ccl/2021.10.0/lib/cpu_gpu_dpcpp -L/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/lib -L/opt/pppl/intel/2023.2.0/tbb/2021.10.0/lib/intel64/gcc4.8 -L/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin -L/usr/lib/gcc/x86_64-redhat-linux/11 -Wl,-rpath,/opt/pppl/intel/2023.2.0/compiler/2023.2.0/linux/compiler/lib/intel64_lin -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/openmpi-4.1.5-u5n2mrzof35eduxvqlsilu2akf43ewup/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/hwloc-2.9.1-tfqm6yyezgbdjtpzlcsmpkcew2omhyok/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-11.3.1/libevent-2.1.12-mbhu3djdevlnislas2gssh5uy4obbuki/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/intel-2023.2.0/pmix-5.0.1-5pfyaosw4gg7t3avbartnn55yi2yhwru/lib -Wl,-rpath,/opt/pppl/pkgs-spack/linux-rocky9-zen3/gcc-13.2.0/slurm-23-02-5-1-rjtbt7npfgbrtdjr4b2qnd3gfrospvcn/lib -lpetsc -lHYPRE -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lhdf5hl_fortran -lhdf5_fortran -lhdf5_hl -lhdf5 -lgmp -lX11 -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lifport -lifcoremt_pic -limf -lsvml -lm -lipgo -lirc -lpthread -lgcc_s -lirc_s -ldl -lmpi_cxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl -lquadmath -lmpi_cxx -lmpi -limf -lsvml -lirng -lstdc++ -lm -lipgo -ldecimal -lcilkrts -lgcc_s -lirc -lirc_s -ldl
	PETSC_FC_INCLUDES = -I/p/swim/jchen/PETSC/petsc.20231214/include -I/p/swim/jchen/PETSC/petsc.20231214/real-spark-gts/include
endif

   FFTLIB:=
   ifeq ($(NAG),y)   # Use nag libraries
	NAG_FILES:=fake_nag.o nag.o
	override NAG:=-L$(VPATH)/cmlib/$(OBJDIR) -lcmlib
   else   # Use nag replacement
	NAG_FILES:=fake_nag.o
	override NAG:=-L$(VPATH)/cmlib/$(OBJDIR) -lcmlib -L$(VPATH)/quadpack -lquadpack
   endif
   FEM:=-L$(VPATH)/fem/$(OBJDIR) -lfemlib

CFLAGS:=$(CFLAGS) -DAdd_ -I${PSPLINE_INC} $(SPRNG_INC) $(PETSC_FC_INCLUDES)
FFLAGS:=$(FFLAGS) -I${PSPLINE_INC} -I${SPRNG_INC} -I${PSPLINE_MOD} -I${SPRNG_MOD} $(PETSC_FC_INCLUDES)

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

LIB:=$(NAG) $(FEM) $(LOCAL_LIB) $(FFTLIB) $(PETSC_WITH_EXTERNAL_LIB) $(ADIOS_FLIB)

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

