
include ../../make.inc

.suffixes : .mod .o .f90

%.o:	%.f90
	$(FC) $(FLAGS) $(FFTWINCLUDE) -c $<

SRC = global_parameter.f90 MPI_mod.f90 MATH_mod.f90  
OBJ = $(SRC:.f90=.o)

compile: $(OBJ)
	$(ARCH) $(ARFLAG) $(LLIB) $(OBJ)
	$(RANLIB) $(LLIB)

clean:
	rm -f *.o *.mod *~ 


