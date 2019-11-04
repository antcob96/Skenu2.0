COMPILER = gfortran

FLAGS = 

EXEC = mechanics_computational_set

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

parametric_curves.o: parametric_curves.f90 types.o 
	$(COMPILER) $(FLAGS) -c $<

satellite_orbit.o: satellite_orbit.f90 types.o  
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
