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

FEMOBJ:=fem.o therm2d.o lagrange.o utility.o
MESHOBJ:=mesh_gen.o
TOBJ:=triangle.o

libfemlib.a: $(FEMOBJ) $(MESHOBJ) $(TOBJ)
	$(FC) $(FFLAGS) $(PETSCFINCLUDE) -I../ -c fem.F90
	$(FC) $(FFLAGS)  -c therm2d.F90
	ar crv libfemlib.a $(FEMOBJ) $(MESHOBJ) $(TOBJ)
	ranlib libfemlib.a

$(FEMOBJ):
$(MESHOBJ):
$(TOBJ): triangle.c
	$(CC) $(DEBUG) -DREDUCED -DCDT_ONLY -DTRILIBRARY -DANSI_DECLARATOR -c $< -o $@


.f90.o :
	$(FC) $(FFLAGS) -c $<

.F90.o :
	$(FC) $(FFLAGS) $(PETSCFINCLUDE) -c $<

.F.o :
	$(FC) $(FFLAGS) -c $<

.c.o   :
	$(CC) $(CFLAGS) -c $<

clean:
	-rm *.o *.a *.mod *.lst
