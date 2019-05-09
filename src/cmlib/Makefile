####################################################################
ifeq (,$(filter $(OBJDIR),$(notdir $(CURDIR))))

# The following section determines the machine and code version
# and prepares a directory for object files and modules
#--------------------------------------------------------------

include ../subtarget.mk

else

# The machine-independent parts of the makefile are specified in this section.
# Machine-dependent parts should be put in makelib_f90
# ----------------------------------------------------------------------------

VPATH= $(SRCDIR)
vpath %.f .:$(SRCDIR)


####################################################################
#
#  The variables FC and FFLAGS needs to be set before executing make
#  or they have to be set on the make command line:
#         e.g.  make FC=lf95 FFLAGS="-O --dbl"
#
####################################################################

#
#ifeq ($(FC),ftn)
#   FFLAGS:= -r8 -Kieee -fastsse
#endif
#
#ifeq ($(FC),lf95)
#   FFLAGS:= -O --dbl --ap --pca --trace -C
#endif
#
#ifeq ($(FC),pgf90)
#   FFLAGS:= -r8 -Kieee -fastsse -Mbounds
#endif
#
#ifeq ($(DEBUG),y)
#   FFLAGS:=-g $(FFLAGS)
#endif

OBJ:= adapt.o bfqad.o bint4.o bintk.o bnfac.o bnslv.o bsgq8.o bspdr.o bspev.o bsppp.o bspvd.o bspvn.o bsqad.o bsrl.o bvalu.o fdump.o flsm.o gl15t.o i1mach.o intrv.o isamax.o j4save.o nxprt.o pfqad.o ppgq8.o ppqad.o ppval.o q1da.o q1dax.o q1db.o r1mach.o rank.o sort.o symrl.o uni.o wht.o xerabt.o xerctl.o xerprt.o xerror.o xerrwv.o xersav.o xgetua.o zeroin.o

libcmlib.a : $(OBJ)
	ar crv libcmlib.a *.o
	ranlib libcmlib.a

$(OBJ):


.f90.o :
	$(FC) $(FFLAGS) -c $<

.F90.o :
	$(FC) $(FFLAGS) $(PETSCFINCLUDE) -c $<

.F.o :
	$(FC) $(FFLAGS) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<

.c.o   :
	$(CC) $(CFLAGS) -c $<

#moved to ../subtarget.mk
#clean:
#	-rm *.o *.a *.lst
endif
