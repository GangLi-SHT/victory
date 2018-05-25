include ../../make.inc

.suffixes : .mod .o .f90

%.o:    %.f90
	$(FC) $(FLAGS) $(MODPATH) -c $<

SRC = parquet_util.f90 parquet_kernel.f90 parquet_equation.f90 parquet_selfenergy.f90 parquet_phys.f90 main.f90
OBJ = $(SRC:.f90=.o)

app:  $(OBJ)
	$(FC) -O3 -o $(BIN)/victory $(OBJ) $(LIB) $(FLAGS) $(MODPATH)  

clean:
	rm -f *.o *.mod *~ 


