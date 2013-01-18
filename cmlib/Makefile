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

libcmlib.a :
	$(FC) -c $(FFLAGS) *.f
	ar crv libcmlib.a *.o
	ranlib libcmlib.a

clean:
	-rm *.o *.a
