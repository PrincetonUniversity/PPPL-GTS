.SUFFIXES:

#we don't define here, we inherited from GTS/target.mk
#OBJDIR := object

#note: SRCDIR=GTS/fem, CURDIR=GTS/fem, option "-C $@" jump to $(OBJDIR)
#and VPATH in GTS/fem/Makefile will tell "make" to find source code from GTS/fem
#note: SRCDIR=GTS/cmlib, CURDIR=GTS/cmlib, option "-C $@" jump to $(OBJDIR)
#and VPATH in GTS/cmlib/Makefile will tell "make" to find source code from GTS/cmlib
MAKETARGET = $(MAKE) --no-print-directory -C $@ -f $(CURDIR)/Makefile \
	SRCDIR=$(CURDIR) \
	$(MAKECMDGOALS)

.PHONY: $(OBJDIR)
$(OBJDIR):
	+@[ -d $@ ] || mkdir -p $@
	+@$(MAKETARGET)

Makefile : ;
%.mk :: ;

% :: $(OBJDIR) ; :

.PHONY: clean
clean :
	rm -rf $(OBJDIR)
