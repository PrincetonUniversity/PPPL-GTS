#GIT information
MAKE_SRC_FOLDER :=  $(shell mkdir $(VPATH)/$(OBJDIR)/src_files)
COPY_SRC_FILES :=  $(shell cp $(VPATH)/* $(VPATH)/$(OBJDIR)/src_files)
GIT_CLEAN 	:= $(shell rm $(VPATH)/$(OBJDIR)/{git_info.o,gts_git_diff.out,gts_version_info.out})
GIT_VERSION := $(shell git --no-pager describe --tags --always)
GIT_BRANCH	:= $(shell git branch | grep \* | cut -c 3-)
GIT_COM_MSG := $(shell git show -s --format=%B | head -n 1 | sed "s/'//g" | sed 's/"//g')
GIT_COM_SHA := $(shell git rev-parse --verify HEAD)
GIT_WTREE   := $()
GIT_DATE    := $(firstword $(shell git --no-pager show --date=iso-strict --format="%ad" --name-only))
BUILD_DATE  := $(shell date --iso=seconds)
BUILD_DIR   := $(shell pwd)

GIT_STATUS_RESULT:= $(shell git status > gts_git_status.out)
GIT_DIFF_RESULT:= $(shell git add -N -u)
GIT_DIFF_RESULT:= $(shell git diff > gts_git_diff.out)
GIT_DIFF_RESULT:= $(shell git reset HEAD)
GIT_DIFF_XXD 	:= $(shell xxd -i gts_git_diff.out > gts_git_diff.h)

#Module Information
MODULE_FILES := $(shell module list &> gts_module_files.out)
MODULE_FILES := $(shell xxd -i gts_module_files.out > gts_module_files.h)

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


all:
	@echo -e "\n\n"
	@echo "============ GTS VERSION INFORMATION ============" | tee gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@echo "       VERSION: $(GIT_VERSION)" | tee -a gts_version_info.out
	@echo "        BRANCH: $(GIT_BRANCH)" 	| tee -a gts_version_info.out
	@echo "   COMMIT  MSG: $(GIT_COM_MSG)" | tee -a gts_version_info.out
	@echo "   COMMIT  SHA: $(GIT_COM_SHA)" | tee -a gts_version_info.out
	@echo "   COMMIT DATE: $(GIT_DATE)" 	| tee -a gts_version_info.out
	@echo "  WORKING TREE: $(GIT_WTREE)" 	| tee -a gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@echo "    BUILD DATE: $(BUILD_DATE)" 	| tee -a gts_version_info.out
	@echo "    BUILD  DIR: $(BUILD_DIR)" 	| tee -a gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@echo "=================================================" 	| tee -a gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@echo "If \"WORKING TREE\" is not clean, check \"git_diff.out\" for details" 	| tee -a gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@echo " 										" 	| tee -a gts_version_info.out
	@cat gts_module_files.out | tee -a gts_version_info.out

git_info.o: git_info.c
		$(CC) -c $(COPT) -I$(VPATH)/$(OBJDIR) $(GIT_FLAGS) $< -o $@