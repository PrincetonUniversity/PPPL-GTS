
  CPP := mpicxx
  CC := mpicc
  F90 := mpif90
  F77 := mpif77
  LOADER := mpif90 -Wl,-z,muldefs -Wl,-z,noexecstack

#GIT information
GIT_CLEAN 	:= $(shell rm $(VPATH)/$(OBJDIR)/{git_info.o,gts_git_diff.out,gts_version_info.out})
GIT_VERSION := $(shell git --no-pager describe --tags --always)
GIT_BRANCH	:= $(shell git branch | grep \* | cut -c 3-)
GIT_COM_MSG := $(shell git show -s --format=%B | head -n 1 | sed "s/'//g" | sed 's/"//g')
GIT_COM_SHA := $(shell git rev-parse --verify HEAD)
GIT_WTREE   := $()
GIT_DATE    := $(firstword $(shell git --no-pager show --date=iso-strict --format="%ad" --name-only))
BUILD_DATE  := $(shell date --iso=seconds)

GIT_DIFF_RESULT:= $(shell git add -N -A)
GIT_DIFF_RESULT:= $(shell git diff > gts_git_diff.out)
GIT_DIFF_RESULT:= $(shell git reset HEAD)
GIT_DIFF_XXD 	:= $(shell xxd -i gts_git_diff.out > gts_git_diff.h)

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
 					-DBUILD_DATE="\"$(BUILD_DATE)\""


ifeq ($(F90),mpif90)
       ifeq ($(PE_ENV),PGI)   # Using PGI compiler
          FFLAGS:= -r8 -Kieee -Dflush_=flush -I$(VPATH)/$(OBJDIR) -I$(VPATH)/fem/$(OBJDIR)
          OPT:=-O2
          COPT:=-O2
          CFLAGS:=
          OMPFLAG:=-mp=nonuma
          NOOMPFLAG:=
          SPRNG_DIR=/projects/M3DC1/LIB/PGI195
          SPRNG_LIB:=-L${SPRNG_DIR}/lib -lsprng
          SPRNG_INC:=-I${SPRNG_DIR}/include
          PSPLINE_DIR=/projects/M3DC1/LIB/PGI195
          PSPLINE_INC:=-I${PSPLINE_DIR}/mod
          PSPLINE_LIB:=-L${PSPLINE_DIR}/lib -lpspline
          PETSC_DIR=/projects/M3DC1/PETSC/petsc-3.12.0
          PETSC_ARCH=traverse-pgi-openmpi-gts-2019
          HDF5_LIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lhdf5_hl -lhdf5
          NETCDF_LIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lnetcdf
          PETSC_FINC:= -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
          PETSC_FLIB:=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lHYPRE -lhdf5 -lz -Wl,--start-group /opt/pgi/19.5/linuxpower/19.5/lib/liblapack.a /opt/pgi/19.5/linuxpower/19.5/lib/libblas.a -Wl,--end-group -lpthread -lm -ldl
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

   ifeq ($(OPENMP),y)
     OPT:=$(OMPFLAG) $(OPT)
     COPT:=$(OMPFLAG) $(COPT)
     GOPT:= $(OMPFLAG) $(GOPT)
   else
     OPT:=$(NOOMPFLAG) $(OPT)
     COPT:=$(NOOMPFLAG) $(COPT)
     GOPT:= $(NOOMPFLAG) $(GOPT)
   endif

LIB:=$(NAG) $(FEM) $(SPRNG_LIB) $(PSPLINE_LIB) $(HDF5_LIB) $(FFTLIB) $(PETSC_FLIB)

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


all:
	@echo -e "\n\n"
	@echo "============ GTS VERSION INFORMATION ============" | tee gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@echo "       VERSION: $(GIT_VERSION)" | tee -a gts_version_info.out
	@echo "        BRANCH: $(GIT_BRANCH)" 	| tee -a gts_version_info.out
	@echo "    COMMIT MSG: $(GIT_COM_MSG)" | tee -a gts_version_info.out
	@echo "    COMMIT SHA: $(GIT_COM_SHA)" | tee -a gts_version_info.out
	@echo "      GIT DATE: $(GIT_DATE)" 	| tee -a gts_version_info.out
	@echo "    BUILD DATE: $(BUILD_DATE)" 	| tee -a gts_version_info.out
	@echo "  WORKING TREE: $(GIT_WTREE)" 	| tee -a gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@echo "=================================================" 	| tee -a gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@echo "If \"WORKING TREE\" is not clean, check \"git_diff.out\" for details" 	| tee -a gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
