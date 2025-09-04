.SUFFIXES:

ifndef _ARCH
_ARCH := $(shell echo $(HOST) | awk '{ sub(/[0-9]+/,""); print }')
export _ARCH
endif

OBJDIR := _$(_ARCH)-$(PE_ENV)$(BIN_POSTFIX)
export OBJDIR 

MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/Makefile \
                 SRCDIR=$(CURDIR) $(MAKECMDGOALS)

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

Makefile : ;
%.mk :: ;

% :: $(OBJDIR) ; :

# The following tag is meant to "clean" the directory by removing the
# executable along with all the object files created by the compilation 
# One only has to run:  gmake clean
.PHONY: clean
clean:
	rm -rf $(OBJDIR)
	(cd cmlib; $(MAKE) clean ; cd ..)
	(cd fem; $(MAKE) clean ; cd ..)
	(cd quadpack; $(MAKE) clean ; cd ..)

# gmake cleangts to clean GTS code
.PHONY: cleangts
cleangts:
	rm -rf $(OBJDIR)

# gmake cleancmlib to clean cmlib directory
.PHONY: cleancmlib
cleancmlib:
	(cd cmlib; $(MAKE) clean ; cd ..)

.PHONY: cleanfem
cleanfem:
	(cd fem; $(MAKE) clean ; cd ..)

.PHONY: cleanquadpack
cleanquadpack:
	(cd quadpack; $(MAKE) clean ; cd ..)
