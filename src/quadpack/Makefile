FC = ifort
FFLAGS = -O2 -fPIC

#QUADPACK_URL = https://www.netlib.org/quadpack/quadpack.tar.gz
#QUADPACK_TAR = quadpack.tar.gz
#QUADPACK_DIR = quadpack

LIB = libquadpack.a

# Default target
all: $(LIB)

# Download QUADPACK
#$(QUADPACK_TAR):
#	curl -L -o $(QUADPACK_TAR) $(QUADPACK_URL)

# Extract QUADPACK
#$(QUADPACK_DIR): $(QUADPACK_TAR)
#	mkdir -p $(QUADPACK_DIR)
#	tar -xvzf $(QUADPACK_TAR) -C $(QUADPACK_DIR)

# List all .f files
F_SRC := $(wildcard *.f)
F90_SRC := $(wildcard *.f90)
OBJ   := $(F_SRC:.f=.o) $(F90_SRC:.f90=.o)

# Build library
$(LIB): $(OBJ)
	ar rcs $(LIB) $(OBJ)

# Pattern rule: compile .f → .o
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean
clean:
	rm -f *.o
	rm -f $(LIB)
#	rm -rf $(QUADPACK_DIR)
#	rm -f $(QUADPACK_TAR)

